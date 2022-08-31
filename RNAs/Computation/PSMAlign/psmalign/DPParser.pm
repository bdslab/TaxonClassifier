package DPParser;

use BracketPairs;
use strict;

use constant DOT => '.';
use constant SPACE => ' ';

sub parse{
    my (undef, $dp_file_path) = @_;

    my ($base_seq, $secondary_structure) = ('', '');

    open (DP, "<$dp_file_path") or die "Cannot open file at $dp_file_path";
    while (<DP>) {
 	if ($_ =~ /^([A-Za-z]+)[\r]*\n$/) {
	    $base_seq = $base_seq . $1;
	}
	elsif ($_ =~ /^([\.\(\)\[\]\{\}<>A-Za-z]+)[\r]*\n$/) {
	    $secondary_structure = $secondary_structure . $1;
	}
	elsif ($_ !~ /^[#>].*/) {
	    die "Unknown input: $_";
	}
    }

    close DP or die "Cannot close file at $dp_file_path";

    if ($base_seq eq '') {
	die 'Base sequence is missing';
    }

    if ($secondary_structure eq '') {
	die 'Secondary structure is missing';
    }

    if (length($base_seq) != length($secondary_structure)) {
	die 'Base sequence length not equal to secondary structure length';
    }

    return $base_seq, $secondary_structure, _group_to_stems($secondary_structure);
}

sub _group_to_stems {
    my $secondary_structure = shift;

    my ($stems, $stem, $stem_outermost_base_pair) = ([], [], []);
    my $unsettled_bracket_upstream_pos = {};
    my $next_paired_pos = {};
    my $last_paired_pos = 0;

    my @structure_symbols = split(//, $secondary_structure);
    my $structure_length = scalar @structure_symbols;

    for (my $i = 0; $i < $structure_length; $i++) {
	my $symbol = $structure_symbols[$i];
	if ($symbol eq DOT) {
	    next;
	}
	elsif (BracketPairs->is_open_bracket($symbol)) {
	    my $unsettled_upstream_pos = $unsettled_bracket_upstream_pos->{$symbol};
	    if (!defined($unsettled_upstream_pos)) {
		$unsettled_upstream_pos = [];
		$unsettled_bracket_upstream_pos->{$symbol} = $unsettled_upstream_pos;
	    }

	    push @{$unsettled_upstream_pos}, $i;

	    if (defined($stem_outermost_base_pair->[0])) {
		push @{$stems}, $stem;
		($stem, $stem_outermost_base_pair) = ([], []);
	    }

	    $next_paired_pos->{$last_paired_pos} = $i;
	    $last_paired_pos = $i;
	}
	else {
	    my $pair_open_bracket = BracketPairs->get_open_bracket($symbol);
	    my $unsettled_upstream_pos = $unsettled_bracket_upstream_pos->{$pair_open_bracket};
	    if (defined($unsettled_upstream_pos) && defined($unsettled_upstream_pos->[0])) {
		my $paired_upstream_pos = pop @{$unsettled_upstream_pos};

		if (defined($stem_outermost_base_pair->[0])) {
		    if ($next_paired_pos->{$paired_upstream_pos} != $stem_outermost_base_pair->[0]) {
			push @{$stems}, $stem;
			($stem, $stem_outermost_base_pair) = ([], []);
		    }

		    $stem_outermost_base_pair = [$paired_upstream_pos, $i];
		    unshift @{$stem}, $stem_outermost_base_pair;
		}
		else {
		    $stem_outermost_base_pair = [$paired_upstream_pos, $i];
		    $stem = [$stem_outermost_base_pair];
		}

		$next_paired_pos->{$last_paired_pos} = $i;
		$last_paired_pos = $i;
	    }
	    else {
		die "Closing bracket $symbol not paired\n";
	    }
	}
    }

    if (!_is_all_open_bracket_settled($unsettled_bracket_upstream_pos)) {
	die "Unpaired open bracket remains\n";
    }

    if (defined($stem_outermost_base_pair->[0])) {
	push @{$stems}, $stem;
    }

    my @sorted_stems = sort {$a->[0][0] <=> $b->[0][0]} @{$stems};

    return \@sorted_stems;
}

sub _is_all_open_bracket_settled {
    my $unsettled_open_bracket_pos = shift;

    foreach (values %{$unsettled_open_bracket_pos}) {
	if (defined($_->[0])) {
	    return 0;
	}
    }

    return 1;
}

1;
