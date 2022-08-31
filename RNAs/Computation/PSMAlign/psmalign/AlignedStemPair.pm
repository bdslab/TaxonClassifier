package AlignedStemPair;

use strict;

=comment
use constant DEFAULT_OPEN_BRACKET => '(';
use constant DEFAULT_CLOSE_BRACKET => ')';

my $default_bracket_regex = '[' . DEFAULT_OPEN_BRACKET . DEFAULT_CLOSE_BRACKET . ']';
=cut

sub new {
    my (undef, $stem_alignment) = @_;

    my $self = {};

    $self->{cost} = ($stem_alignment->[0])->get(0);
    $self->{stem1_aligned_upstream} = $stem_alignment->[1];
    $self->{stem1_aligned_downstream} = $stem_alignment->[2];
    $self->{stem2_aligned_upstream} = $stem_alignment->[3];
    $self->{stem2_aligned_downstream} = $stem_alignment->[4];
=comment
    $self->{stem1_open_bracket} = DEFAULT_OPEN_BRACKET;
    $self->{stem1_close_bracket} = DEFAULT_CLOSE_BRACKET;
    $self->{stem2_open_bracket} = DEFAULT_OPEN_BRACKET;
    $self->{stem2_close_bracket} = DEFAULT_CLOSE_BRACKET;
=cut

    bless $self;

    return $self;
}

=comment
sub set_brackets {
    my ($self, $stem_open_bracket, $stem_close_bracket, $stem_num) = @_;

    $self->{'stem' . $stem_num . '_open_bracket'} = $stem_open_bracket;
    $self->{'stem' . $stem_num . '_close_bracket'} = $stem_close_bracket;
}
=cut

sub self_swap {
    my $self = shift;

    my $str_buffer = $self->{stem1_aligned_upstream};
    $self->{stem1_aligned_upstream} = $self->{stem2_aligned_upstream};
    $self->{stem2_aligned_upstream} = $str_buffer;

    $str_buffer = $self->{stem1_aligned_downstream};
    $self->{stem1_aligned_downstream} = $self->{stem2_aligned_downstream};
    $self->{stem2_aligned_downstream} = $str_buffer;
}

sub set_and_swap {
    my ($self, $input_stem_pair) = @_;

    $self->{cost} = $input_stem_pair->get_cost();
    $self->{stem1_aligned_upstream} = $input_stem_pair->{stem2_aligned_upstream};
    $self->{stem1_aligned_downstream} = $input_stem_pair->{stem2_aligned_downstream};
    $self->{stem2_aligned_upstream} = $input_stem_pair->{stem1_aligned_upstream};
    $self->{stem2_aligned_downstream} = $input_stem_pair->{stem1_aligned_downstream};
=comment
    $self->{stem1_open_bracket} = $input_stem_pair->{stem2_open_bracket};
    $self->{stem1_close_bracket} = $input_stem_pair->{stem2_close_bracket};
    $self->{stem2_open_bracket} = $input_stem_pair->{stem1_open_bracket};
    $self->{stem2_close_bracket} = $input_stem_pair->{stem1_close_bracket};
=cut
}

sub restore_seq_bases {
    my ($self, $stem1_upstream_org_seq, $stem1_downstream_org_seq, $stem2_upstream_org_seq, $stem2_downstream_org_seq) = @_;

    my $stem1_upstream_aligned_seq = ($self->{stem1_aligned_upstream})->get(1);
    my $stem1_downstream_aligned_seq = ($self->{stem1_aligned_downstream})->get(1);
    my $stem2_upstream_aligned_seq = ($self->{stem2_aligned_upstream})->get(1);
    my $stem2_downstream_aligned_seq = ($self->{stem2_aligned_downstream})->get(1);

    if ($stem1_upstream_aligned_seq =~ /\?/) {
	my $restored_stem1_upstream_aligned_seq = _restore_stream_seq_bases($stem1_upstream_aligned_seq, $stem1_upstream_org_seq);
	$self->{stem1_aligned_upstream}->set(1, $restored_stem1_upstream_aligned_seq);
    }

    if ($stem1_downstream_aligned_seq =~ /\?/) {
	my $restored_stem1_downstream_aligned_seq = _restore_stream_seq_bases($stem1_downstream_aligned_seq, $stem1_downstream_org_seq);
	$self->{stem1_aligned_downstream}->set(1, $restored_stem1_downstream_aligned_seq);
    }

    if ($stem2_upstream_aligned_seq =~ /\?/) {
	my $restored_stem2_upstream_aligned_seq = _restore_stream_seq_bases($stem2_upstream_aligned_seq, $stem2_upstream_org_seq);
	$self->{stem2_aligned_upstream}->set(1, $restored_stem2_upstream_aligned_seq);
    }

    if ($stem2_downstream_aligned_seq =~ /\?/) {
	my $restored_stem2_downstream_aligned_seq = _restore_stream_seq_bases($stem2_downstream_aligned_seq, $stem2_downstream_org_seq);
	$self->{stem2_aligned_downstream}->set(1, $restored_stem2_downstream_aligned_seq);
    }
}

sub _restore_stream_seq_bases {
    my ($stream_aligned_seq, $stream_org_seq) = @_;

    my $restored_stream_aligned_seq_bases = [];
    my @stream_aligned_seq_bases = split(//, $stream_aligned_seq);
    my $org_seq_base_ptr = 0;

    for (my $i = 0; $i < @stream_aligned_seq_bases; $i++) {
	if ($stream_aligned_seq_bases[$i] eq DefaultAlignParams->BASE_SEQ_GAP) {
	    push @{$restored_stream_aligned_seq_bases}, DefaultAlignParams->BASE_SEQ_GAP;
	}
	else {
	    push @{$restored_stream_aligned_seq_bases}, substr($stream_org_seq, $org_seq_base_ptr, 1);
	    $org_seq_base_ptr += 1;
	}
    }

    return join('', @{$restored_stream_aligned_seq_bases});
}

sub get_cost {
    my $self = shift;

    return $self->{cost};
}

=comment
sub get_struct_dp {
    my ($self, $seq_num, $stream_opt) = @_;

    my $selected_stem = 'stem' . $seq_num;

    if ($stream_opt eq 'upstream') {
	return _restore_original_brackets(($self->{$selected_stem . '_aligned_' . $stream_opt})->get(0), $self->{$selected_stem . '_open_bracket'});
    }
    elsif ($stream_opt eq 'downstream') {
	return _restore_original_brackets(($self->{$selected_stem . '_aligned_' . $stream_opt})->get(0), $self->{$selected_stem . '_close_bracket'});
    }

    return undef;
}
=cut

sub get_base_seq {
    my ($self, $seq_num, $stream_opt) = @_;

    my $selected_stem = 'stem' . $seq_num;

    return ($self->{$selected_stem . '_aligned_' . $stream_opt})->get(1);
}

sub get_del_edge_len {
    my ($self, $seq_num, $stream_opt, $end_opt) = @_;

    my $selected_stem = 'stem' . $seq_num;
    if ($end_opt eq 'left') {
	return ($self->{$selected_stem . '_aligned_' . $stream_opt})->get(2);
    }
    elsif ($end_opt eq 'right') {
	return ($self->{$selected_stem . '_aligned_' . $stream_opt})->get(3);
    }

    return undef;
}

sub get_del_edge_seq {
    my ($self, $seq_num, $stream_opt, $end_opt) = @_;

    my $del_edge_len = get_del_edge_len($self, $seq_num, $stream_opt, $end_opt);
    if (!defined($del_edge_len) || $del_edge_len == 0) {
	return '';
    }

    my $base_seq = get_base_seq($self, $seq_num, $stream_opt);
    if ($end_opt eq 'left') {
	return substr($base_seq, 0, $del_edge_len);
    }
    elsif ($end_opt eq 'right') {
	return substr($base_seq, (-1 * $del_edge_len));
    }

    return '';
}

sub get_total_del_edge_len {
    my $self = shift;

    return (($self->{stem1_aligned_upstream})->get(2) + ($self->{stem1_aligned_upstream})->get(3) +
	    ($self->{stem1_aligned_downstream})->get(2) + ($self->{stem1_aligned_downstream})->get(3) +
	    ($self->{stem2_aligned_upstream})->get(2) + ($self->{stem2_aligned_upstream})->get(3) +
	    ($self->{stem2_aligned_downstream})->get(2) + ($self->{stem2_aligned_downstream})->get(3));
}

=comment
sub _restore_original_brackets {
    my ($region_dp, $org_brac) = @_;

    if ($org_brac ne DEFAULT_OPEN_BRACKET && $org_brac ne DEFAULT_CLOSE_BRACKET) {
	$region_dp =~ s/$default_bracket_regex/$org_brac/g;
    }

    return $region_dp;
}
=cut

1;
