package StemGraph;

use rna_stem_align;
use strict;

sub new {
    my (undef, $base_seq, $secondary_structure, $stems) = @_;

    my ($vertex_attrs, $stem_types, $nest_props) = ([], [], []);
    my $edge_labels = {};

    my $split_stems = _split_base_pair_stems($stems);

    my $vertex_count = @{$split_stems};
    for (my $i = 0; $i < $vertex_count; $i++) {
	$stem_types->[$i] = 'H';
	$nest_props->[$i] = 0;

	my ($stem_outermost_pair_start, $stem_outermost_pair_end) = @{$split_stems->[$i][0]};
	for (my $j = 0; $j < $i; $j++) {
	    my ($prev_stem_outermost_pair_start, $prev_stem_outermost_pair_end) = @{$split_stems->[$j][0]};
	    if ($stem_outermost_pair_start > $prev_stem_outermost_pair_end) {
		$edge_labels->{$j . '-' . $i} = 'P';
	    }
	    elsif ($stem_outermost_pair_end < $prev_stem_outermost_pair_end) {
		$edge_labels->{$j . '-' . $i} = 'N';
		if ($stem_types->[$j] ne 'T') {
		    $stem_types->[$j] = 'M';
		}

		$nest_props->[$j] = 1;
	    }
	    else {
		$edge_labels->{$j . '-' . $i} = 'K';
		$stem_types->[$i] = 'T';
		$stem_types->[$j] = 'T';
	    }
	}

	my $base_pairs = $split_stems->[$i];
	my $base_pair_count = @{$base_pairs};

	my $stem_upstream_start = $base_pairs->[0][0];
	my $stem_upstream_length = $base_pairs->[$base_pair_count - 1][0] - $stem_upstream_start + 1;
	my $stem_downstream_start = $base_pairs->[$base_pair_count - 1][1];
	my $stem_downstream_length = $base_pairs->[0][1] - $stem_downstream_start + 1;

	my $attrs = {};
	$attrs->{base_pair_count} = $base_pair_count;
	$attrs->{base_pairs} = $base_pairs;
	$attrs->{upstream_base_seq} = substr($base_seq, $stem_upstream_start ,$stem_upstream_length);
	$attrs->{downstream_base_seq} = substr($base_seq, $stem_downstream_start, $stem_downstream_length);
	$attrs->{upstream_start} = $stem_upstream_start;
	$attrs->{upstream_length} = $stem_upstream_length;
	$attrs->{downstream_start} = $stem_downstream_start;
	$attrs->{downstream_length} = $stem_downstream_length;
#	($attrs->{pseudo_stem_pairs}, $attrs->{pseudo_pair_pos_map}) = _convert_to_pseudo_base_pairs($base_pairs, $base_pair_count);
#	$attrs->{pseudo_base_pairs} = _convert_to_pseudo_base_pairs($base_pairs, $base_pair_count);
	$attrs->{pseudo_base_pairs} = GraphMatchTools->convert_to_pseudo_base_pairs($base_pairs);
#	$attrs->{open_bracket} = substr($secondary_structure, $stem_upstream_start, 1);
#	$attrs->{close_bracket} = substr($secondary_structure, $stem_downstream_start, 1);
	$vertex_attrs->[$i] = $attrs;
    }

#    for (my $i = 0; $i < $vertex_count; $i++) {
#	my $attrs = $vertex_attrs->[$i];
#	$attrs->{stem_type} = $stem_types->[$i];
#    }

    my ($nest_stem_vertices, $non_nest_stem_vertices) = ([], []);
    for (my $i = 0; $i < $vertex_count; $i++) {
	if ($nest_props->[$i]) {
	    push @{$nest_stem_vertices}, $i;
	}
	else {
	    push @{$non_nest_stem_vertices}, $i;
	}
    }

    my $self = {};
    $self->{vertex_count} = $vertex_count;
    $self->{edge_labels} = $edge_labels;
    $self->{vertex_attrs} = $vertex_attrs;
    $self->{stem_types} = $stem_types;
    $self->{nest_props} = $nest_props;
    $self->{nest_stem_vertices} = $nest_stem_vertices;
    $self->{non_nest_stem_vertices} = $non_nest_stem_vertices;

    bless $self;

    return $self;
}

sub _split_base_pair_stems {
    my $stems = shift;

=comment
    my $split_stems = [];

    my $stem_count = @{$stems};
    for (my $i = 0; $i < $stem_count; $i++) {
	my $is_non_nest_stem = 1;
	for (my $j = $i + 1; $j < $stem_count; $j++) {
	    if ($stems->[$j][0][1] < $stems->[$i][0][1]) {
		$is_non_nest_stem = 0;
		last;
	    }
	    
	    if ($stems->[$i][0][1] < $stems->[$j][0][0]) {
		last;
	    }
	}

	my $stem = $stems->[$i];
	my $stem_size = @{$stem};
	my $split_stem = [$stem->[0]];
#	my $acc_stem_size = 1;

	for (my $j = 1; $j < $stem_size; $j++) {
	    my $upstream_stem_gap_size = $stem->[$j][0] - $stem->[$j - 1][0] - 1;
	    my $downstream_stem_gap_size = $stem->[$j - 1][1] - $stem->[$j][1] - 1;
#	    if (($is_non_nest_stem && ($acc_stem_size > 1) && ($upstream_stem_gap_size >= $acc_stem_size) && ($downstream_stem_gap_size >= $acc_stem_size)) ||
#		($upstream_stem_gap_size > $stem_size || $downstream_stem_gap_size > $stem_size)) {
	    if (($is_non_nest_stem && ($j > 1) && ($upstream_stem_gap_size >= $j) && ($downstream_stem_gap_size >= $j)) ||
		($upstream_stem_gap_size > $stem_size || $downstream_stem_gap_size > $stem_size)) {
		push @{$split_stems}, $split_stem;
		$split_stem = [];
	    }
#	    else {
#		$acc_stem_size++;
#	    }

	    push @{$split_stem}, $stem->[$j];
	}

	if (defined($split_stem->[0])) {
	    push @{$split_stems}, $split_stem;
	}
    }

    return $split_stems;
=cut
    my $split_stems = [];

    my $stem_count = @{$stems};
    for (my $i = 0; $i < $stem_count; $i++) {
	my $is_non_nest_stem = 1;
	for (my $j = $i + 1; $j < $stem_count; $j++) {
	    if ($stems->[$j][0][1] < $stems->[$i][0][1]) {
		$is_non_nest_stem = 0;
		last;
	    }
	    
	    if ($stems->[$i][0][1] < $stems->[$j][0][0]) {
		last;
	    }
	}

	my $stacks = _get_stacks_from_stem($stems->[$i]);
	my $merged_stacks = _merge_stacks($stacks, $is_non_nest_stem);
	my $num_of_stacks = @{$stacks};
	my $num_of_merged_stacks = @{$merged_stacks};

	while ($num_of_merged_stacks < $num_of_stacks) {
	    $stacks = $merged_stacks;
	    $merged_stacks = _merge_stacks($stacks, $is_non_nest_stem);
	    $num_of_stacks = @{$stacks};
	    $num_of_merged_stacks = @{$merged_stacks};
	}

	foreach my $merged_stack (@{$merged_stacks}) {
	    foreach (@{$merged_stack}) {
#		print '[' . $_->[0] . ', ' . $_->[1] . '] ';
	    }
#	    print "\n";
	}
#	print "-------------\n";

	push @{$split_stems}, @{$merged_stacks};
    }

    return $split_stems;
}

sub _get_stacks_from_stem {
    my $stem = shift;

    my $stacks = [];
    my $stack = [$stem->[0]];

    for (my $i = 1; $i < @{$stem}; $i++) {
	if (($stem->[$i][0] != ($stem->[$i - 1][0] + 1)) || ($stem->[$i][1] != ($stem->[$i - 1][1] - 1))) {
	    push @{$stacks}, $stack;
	    $stack = [];
	}

	push @{$stack}, $stem->[$i];
    }

    push @{$stacks}, $stack;

=comment
    foreach $stack (@{$stacks}) {
	foreach (@{$stack}) {
	    print '[' . $_->[0] . ', ' . $_->[1] . '] ';
	}
	print "\n";
    }
    print "-------------\n";
=cut

    return $stacks;
}

sub _merge_stacks {
    my ($stacks, $is_non_nest_stem) = @_;

    my $merged_stacks = [];

    my $merged_stack = $stacks->[0];

    for (my $i = 1; $i < @{$stacks}; $i++) {
	my $prev_stack_upstream_size = $stacks->[$i - 1][-1][0] - $stacks->[$i - 1][0][0] + 1;
	my $prev_stack_downstream_size = $stacks->[$i - 1][0][1] - $stacks->[$i - 1][-1][1] + 1;
	my $curr_stack_upstream_size = $stacks->[$i][-1][0] - $stacks->[$i][0][0] + 1;
	my $curr_stack_downstream_size = $stacks->[$i][0][1] - $stacks->[$i][-1][1] + 1;

#	my $prev_stack_size = @{$stacks->[$i - 1]};
#	my $curr_stack_size = @{$stacks->[$i]};
	my $upstream_gap_size = $stacks->[$i][0][0] - $stacks->[$i - 1][-1][0] - 1;
	my $downstream_gap_size = $stacks->[$i - 1][-1][1] - $stacks->[$i][0][1] - 1;

	if ($is_non_nest_stem && ($prev_stack_upstream_size > 1 || $prev_stack_downstream_size > 1) && ($curr_stack_upstream_size > 1 || $curr_stack_downstream_size > 1)) {
	    if ($upstream_gap_size > $prev_stack_upstream_size && $upstream_gap_size > $curr_stack_upstream_size &&
		$downstream_gap_size > $prev_stack_downstream_size && $downstream_gap_size > $curr_stack_downstream_size) {
		push @{$merged_stacks}, $merged_stack;
		$merged_stack = [];
	    }
	}
	else {
	    if ($upstream_gap_size > $prev_stack_upstream_size || $upstream_gap_size > $curr_stack_upstream_size ||
		$downstream_gap_size > $prev_stack_downstream_size || $downstream_gap_size > $curr_stack_downstream_size) {
		push @{$merged_stacks}, $merged_stack;
		$merged_stack = [];
	    }
	}

	push @{$merged_stack}, @{$stacks->[$i]};
    }

    push @{$merged_stacks}, $merged_stack;

    return $merged_stacks;
}

=comment
sub _convert_to_pseudo_base_pairs {
    my ($base_pairs, $base_pair_count) = @_;

    my $pseudo_base_pairs = new rna_stem_align::vector_i2();
#    my $pseudo_pair_pos_map = {};

    my ($stem_outermost_pair, $stem_innermost_pair) = ($base_pairs->[0], $base_pairs->[$base_pair_count - 1]);
    my $stem_region_length = $stem_innermost_pair->[0] - $stem_outermost_pair->[0] + $stem_outermost_pair->[1] - $stem_innermost_pair->[1] + 2;
    my $last_pseudo_base_pair = new rna_stem_align::vector_i();
    $last_pseudo_base_pair->push(0);
    $last_pseudo_base_pair->push($stem_region_length - 1);
    $pseudo_base_pairs->push($last_pseudo_base_pair);
#    $pseudo_pair_pos_map->{1} = $stem_outermost_pair->[0];
#    $pseudo_pair_pos_map->{$stem_region_length} = $stem_outermost_pair->[1];

    for (my $i = 1; $i < $base_pair_count; $i++) {
	my $curr_pseudo_base_pair_upstream_pos = $last_pseudo_base_pair->get(0) + $base_pairs->[$i][0] - $base_pairs->[$i - 1][0];
	my $curr_pseudo_base_pair_downstream_pos = $last_pseudo_base_pair->get(1) + $base_pairs->[$i][1] - $base_pairs->[$i - 1][1];
	$last_pseudo_base_pair = new rna_stem_align::vector_i();
	$last_pseudo_base_pair->push($curr_pseudo_base_pair_upstream_pos);
	$last_pseudo_base_pair->push($curr_pseudo_base_pair_downstream_pos);
	$pseudo_base_pairs->push($last_pseudo_base_pair);
#	$pseudo_pair_pos_map->{$curr_pseudo_base_pair_upstream_pos} = $base_pairs->[$i][0];
#	$pseudo_pair_pos_map->{$curr_pseudo_base_pair_downstream_pos} = $base_pairs->[$i][1];
    }

#    return $pseudo_base_pairs, $pseudo_pair_pos_map;
    return $pseudo_base_pairs;
}
=cut

sub get_vertex_count {
    my $self = shift;

    return $self->{vertex_count};
}

sub get_stem_type_at {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return '';
    }

    my $stem_types = $self->{stem_types};

    return $stem_types->[$vertex_num];
}

sub is_nest_stem_vertex {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return 0;
    }

    my $nest_props = $self->{nest_props};

    return $nest_props->[$vertex_num];
}

sub get_nest_stem_vertices {
    my $self = shift;

    return $self->{nest_stem_vertices};
}

sub get_non_nest_stem_vertices {
    my $self = shift;

    return $self->{non_nest_stem_vertices};
}

sub get_edge_label {
    my ($self, $source_vertex, $dest_vertex) = @_;

    my $edge_labels = $self->{edge_labels};
    my $edge_label = $edge_labels->{$source_vertex . '-' . $dest_vertex};

    if (defined($edge_label)) {
	return $edge_label;
    }

    return '';
}

sub get_vertex_attrs_at {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return [];
    }

    my $vertex_attrs = $self->{vertex_attrs};

    return $vertex_attrs->[$vertex_num];
}

1;
