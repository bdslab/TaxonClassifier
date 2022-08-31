package GraphMatchTools;

use strict;

sub convert_to_pseudo_base_pairs {
    my (undef, $base_pairs) = @_;

    my $pseudo_base_pairs = new rna_stem_align::vector_i2();

    my ($stem_outermost_pair, $stem_innermost_pair) = ($base_pairs->[0], $base_pairs->[-1]);
    my $stem_region_length = $stem_innermost_pair->[0] - $stem_outermost_pair->[0] + $stem_outermost_pair->[1] - $stem_innermost_pair->[1] + 2;
    my $last_pseudo_base_pair = new rna_stem_align::vector_i();
    $last_pseudo_base_pair->push(0);
    $last_pseudo_base_pair->push($stem_region_length - 1);
    $pseudo_base_pairs->push($last_pseudo_base_pair);

    for (my $i = 1; $i < @{$base_pairs}; $i++) {
	my $curr_pseudo_base_pair_upstream_pos = $last_pseudo_base_pair->get(0) + $base_pairs->[$i][0] - $base_pairs->[$i - 1][0];
	my $curr_pseudo_base_pair_downstream_pos = $last_pseudo_base_pair->get(1) + $base_pairs->[$i][1] - $base_pairs->[$i - 1][1];
	$last_pseudo_base_pair = new rna_stem_align::vector_i();
	$last_pseudo_base_pair->push($curr_pseudo_base_pair_upstream_pos);
	$last_pseudo_base_pair->push($curr_pseudo_base_pair_downstream_pos);
	$pseudo_base_pairs->push($last_pseudo_base_pair);
    }

    return $pseudo_base_pairs;
}

=comment
sub get_unmapped_stem_vertices {
    my (undef, $stem_vertices1, $stem_vertices2, $fixed_vertex_mapping) = @_;

    my ($mapped_stem_vertices1, $mapped_stem_vertices2) = ({}, {});
    foreach (@{$fixed_vertex_mapping}) {
	$mapped_stem_vertices1->{$_->[0]} = 1;
	$mapped_stem_vertices2->{$_->[1]} = 1;
    }

    my ($unmapped_stem_vertices1, $unmapped_stem_vertices2) = ([], []);
    foreach (@{$stem_vertices1}) {
	if (!exists($mapped_stem_vertices1->{$_})) {
	    push @{$unmapped_stem_vertices1}, $_;
	}
    }

    foreach (@{$stem_vertices2}) {
	if (!exists($mapped_stem_vertices2->{$_})) {
	    push @{$unmapped_stem_vertices2}, $_;
	}
    }

    return $unmapped_stem_vertices1, $unmapped_stem_vertices2;
}

sub partition_unmapped_stem_vertices {
    my (undef, $stem_vertices1, $stem_vertices2, $fixed_vertex_mapping) = @_;

    my $mapped_vertex_pair_count = @{$fixed_vertex_mapping};
    my $ptr = 0;

    my ($unmapped_stem_vertex_partitions1, $unmapped_stem_vertex_partition) = ([], []);

    for (my $i = 0; $i < @{$stem_vertices1}; $i++) {
	if (($ptr == $mapped_vertex_pair_count) || ($stem_vertices1->[$i] < $fixed_vertex_mapping->[$ptr][0])) {
	    push @{$unmapped_stem_vertex_partition}, $stem_vertices1->[$i];
	}
	elsif ($stem_vertices1->[$i] == $fixed_vertex_mapping->[$ptr][0]) {
	    push @{$unmapped_stem_vertex_partitions1}, $unmapped_stem_vertex_partition;
	    $unmapped_stem_vertex_partition = [];
	    $ptr++;
	}
	else {
	    push @{$unmapped_stem_vertex_partitions1}, $unmapped_stem_vertex_partition;
	    $unmapped_stem_vertex_partition = [];
	    $ptr++;
	    $i--;
	}
    }

    push @{$unmapped_stem_vertex_partitions1}, $unmapped_stem_vertex_partition;

    my $unmapped_stem_vertex_partitions2 = [];
    $unmapped_stem_vertex_partition = [];
    $ptr = 0;

    for (my $i = 0; $i < @{$stem_vertices2}; $i++) {
	if (($ptr == $mapped_vertex_pair_count) || ($stem_vertices2->[$i] < $fixed_vertex_mapping->[$ptr][1])) {
	    push @{$unmapped_stem_vertex_partition}, $stem_vertices2->[$i];
	}
	elsif ($stem_vertices2->[$i] == $fixed_vertex_mapping->[$ptr][1]) {
	    push @{$unmapped_stem_vertex_partitions2}, $unmapped_stem_vertex_partition;
	    $unmapped_stem_vertex_partition = [];
	    $ptr++;
	}
	else {
	    push @{$unmapped_stem_vertex_partitions2}, $unmapped_stem_vertex_partition;
	    $unmapped_stem_vertex_partition = [];
	    $ptr++;
	    $i--;
	}
    }

    push @{$unmapped_stem_vertex_partitions2}, $unmapped_stem_vertex_partition;

    return $unmapped_stem_vertex_partitions1, $unmapped_stem_vertex_partitions2;
}
=cut

sub generate_unmapped_stem_vertex_pairs {
    my (undef, $unmapped_stem_vertices1, $unmapped_stem_vertices2, $predefined_vertex_mapping, $stem_graph1, $stem_graph2, $is_check_stem_type) = @_;

    my ($external_edge_label_tokens1, $external_edge_label_tokens2) = ({}, {});

    foreach (@{$unmapped_stem_vertices1}) {
	$external_edge_label_tokens1->{$_} = '';
    }

    foreach (@{$unmapped_stem_vertices2}) {
	$external_edge_label_tokens2->{$_} = '';
    }

    foreach my $mapped_vertex_pair (@{$predefined_vertex_mapping}) {
	foreach (@{$unmapped_stem_vertices1}) {
	    if ($mapped_vertex_pair->[0] < $_) {
		$external_edge_label_tokens1->{$_} = $external_edge_label_tokens1->{$_} . $stem_graph1->get_edge_label($mapped_vertex_pair->[0], $_);
	    }
	    else {
		$external_edge_label_tokens1->{$_} = $external_edge_label_tokens1->{$_} . lc($stem_graph1->get_edge_label($_, $mapped_vertex_pair->[0]));
	    }
	}

	foreach (@{$unmapped_stem_vertices2}) {
	    if ($mapped_vertex_pair->[1] < $_) {
		$external_edge_label_tokens2->{$_} = $external_edge_label_tokens2->{$_} . $stem_graph2->get_edge_label($mapped_vertex_pair->[1], $_);
	    }
	    else {
		$external_edge_label_tokens2->{$_} = $external_edge_label_tokens2->{$_} . lc($stem_graph2->get_edge_label($_, $mapped_vertex_pair->[1]));
	    }
	}
    }

    my $unmapped_stem_vertex_pairs = [];
    if (defined($is_check_stem_type) && $is_check_stem_type) {
	foreach my $unmapped_stem_vertex1 (@{$unmapped_stem_vertices1}) {
	my $is_vertex1_nest = $stem_graph1->is_nest_stem_vertex($unmapped_stem_vertex1);
	    foreach my $unmapped_stem_vertex2 (@{$unmapped_stem_vertices2}) {
		if ($is_vertex1_nest) {
		    if (!$stem_graph2->is_nest_stem_vertex($unmapped_stem_vertex2)) {
			next;
		    }
		}
		else {
		    if ($stem_graph2->is_nest_stem_vertex($unmapped_stem_vertex2)) {
			next;
		    }
		}

		if ($external_edge_label_tokens1->{$unmapped_stem_vertex1} eq $external_edge_label_tokens2->{$unmapped_stem_vertex2}) {
		    push @{$unmapped_stem_vertex_pairs}, [$unmapped_stem_vertex1, $unmapped_stem_vertex2];
		}
	    }
	}
    }
    else {
	foreach my $unmapped_stem_vertex1 (@{$unmapped_stem_vertices1}) {
	    foreach my $unmapped_stem_vertex2 (@{$unmapped_stem_vertices2}) {
		if ($external_edge_label_tokens1->{$unmapped_stem_vertex1} eq $external_edge_label_tokens2->{$unmapped_stem_vertex2}) {
		    push @{$unmapped_stem_vertex_pairs}, [$unmapped_stem_vertex1, $unmapped_stem_vertex2];
		}
	    }
	}
    }

    return $unmapped_stem_vertex_pairs;
}

=comment
sub expand_mscs_mapping {
    my (undef, $candidate_stem_vertices1, $candidate_stem_vertices2, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2,
	$stem_removal_costs1, $stem_removal_costs2) = @_;

#    my $candidate_stem_vertex_pairs = generate_unmapped_stem_vertex_pairs(undef, $candidate_stem_vertices1, $candidate_stem_vertices2,
#									  $mscs_mapping, $stem_graph1, $stem_graph2);
    my $independent_partition_mscs_candidates = _generate_unmapped_non_nest_stem_vertex_pairs($mscs_mapping, $stem_graph1, $stem_graph2);

    my $curr_mscs_mappings = [$mscs_mapping];
    my $total_reduced_cost = 0;

    my $all_partition_mscs_mappings = [];
    foreach my $partition_mscs_candidates (@{$independent_partition_mscs_candidates}) {
	(my ($partition_mscs_mappings, $mscs_reduced_cost), $aligned_stem_pairs) =
	    expand_mscs_mapping_with_candidate_pairs(undef, $partition_mscs_candidates, [], $aligned_stem_pairs, $stem_graph1, $stem_graph2,
						     $stem_removal_costs1, $stem_removal_costs2);
	push @{$all_partition_mscs_mappings}, $partition_mscs_mappings;
	$total_reduced_cost += $mscs_reduced_cost;

#=comment
	foreach (@{$_}) {
	    print '[' . $_->[0] . ', ' . $_->[1] . '] ';
	}
	print "\n---------\n";
#=cut
    }

    my $curr_mscs_mappings = [$mscs_mapping];
    foreach my $partition_mscs_mappings (@{$all_partition_mscs_mappings}) {
	my $expanded_mscs_mappings = [];
	foreach my $partition_mscs_mapping (@{$partition_mscs_mappings}) {
	    foreach my $curr_mscs_mapping (@{$curr_mscs_mappings}) {
		my @expanded_mscs_mapping = (@{$partition_mscs_mapping}, @{$curr_mscs_mapping});
		push @{$expanded_mscs_mappings}, \@expanded_mscs_mapping;
	    }
	}

	$curr_mscs_mappings = $expanded_mscs_mappings;
    }

    return $curr_mscs_mappings, $total_reduced_cost, $aligned_stem_pairs;

#    expand_mscs_mapping_with_candidate_pairs(undef, $candidate_stem_vertex_pairs, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2);
}
=cut

sub expand_mscs_mapping {
    my (undef, $candidate_vertex_pairs, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2,
	$base_seq1, $base_seq2) = @_;

    $candidate_vertex_pairs = _filter_candidate_vertex_pairs($candidate_vertex_pairs, $mscs_mapping, $stem_graph1, $stem_graph2);
    $aligned_stem_pairs = align_stem_pairs(undef, $candidate_vertex_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2);

    my ($expanded_mscs_mappings, $mapping_cost) = _get_mscs($candidate_vertex_pairs, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2,
							    $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);

    return $expanded_mscs_mappings, $aligned_stem_pairs, $mapping_cost;
}

sub _filter_candidate_vertex_pairs {
    my ($candidate_vertex_pairs, $mscs_mapping, $stem_graph1, $stem_graph2) = @_;

    my $filtered_candidate_vertex_pairs = [];
    foreach my $candidate_vertex_pair (@{$candidate_vertex_pairs}) {
	my $is_mscs_mapping_violated = 0;

	foreach (@{$mscs_mapping}) {
	    if ($_->[0] < $candidate_vertex_pair->[0] && $_->[1] < $candidate_vertex_pair->[1]) {
		if ($stem_graph1->get_edge_label($_->[0], $candidate_vertex_pair->[0]) ne $stem_graph2->get_edge_label($_->[1], $candidate_vertex_pair->[1])) {
		    $is_mscs_mapping_violated = 1;
		    last;
		}
	    }
	    elsif ($_->[0] > $candidate_vertex_pair->[0] && $_->[1] > $candidate_vertex_pair->[1]) {
		if ($stem_graph1->get_edge_label($candidate_vertex_pair->[0], $_->[0]) ne $stem_graph2->get_edge_label($candidate_vertex_pair->[1], $_->[1])) {
		    $is_mscs_mapping_violated = 1;
		    last;
		}
	    }
	    else {
		$is_mscs_mapping_violated = 1;
		last;
	    }
	}

	if (!$is_mscs_mapping_violated) {
	    push @{$filtered_candidate_vertex_pairs}, $candidate_vertex_pair;
	}
    }

    return $filtered_candidate_vertex_pairs;
}

sub expand_non_nest_motif_mscs_mapping {
    my (undef, $candidate_vertex_pairs, $non_nest_motif_mscs_mapping, $aligned_motif_pairs, $stem_graph1, $stem_graph2,
	$motif_removal_costs1, $motif_removal_costs2) = @_;

    $candidate_vertex_pairs = _filter_with_motif_mscs_mapping($candidate_vertex_pairs, $non_nest_motif_mscs_mapping);

    (my $expanded_non_nest_motif_mscs_mapping, undef) = _get_mscs($candidate_vertex_pairs, $non_nest_motif_mscs_mapping,
								  $aligned_motif_pairs, $stem_graph1, $stem_graph2,
								  $motif_removal_costs1, $motif_removal_costs2);

    return $expanded_non_nest_motif_mscs_mapping;

=comment
    my @expanded_non_nest_motif_mscs_mappings;
    foreach (@{$new_non_nest_motif_mscs_mappings}) {
	if (defined($non_nest_motif_mscs_mapping->[0])) {
	   my @merged_non_nest_motif_mscs_mapping = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} (@{$non_nest_motif_mscs_mapping}, @{$_});
	   push @expanded_non_nest_motif_mscs_mappings, \@merged_non_nest_motif_mscs_mapping;
	}
	else {
	    push @expanded_non_nest_motif_mscs_mappings, $_;
	}
    }

    if (@expanded_non_nest_motif_mscs_mappings == 0) {
	push @expanded_non_nest_motif_mscs_mappings, $non_nest_motif_mscs_mapping;
    }

    return \@expanded_non_nest_motif_mscs_mappings;
=cut
}

sub _filter_with_motif_mscs_mapping {
    my ($vertex_pairs, $mscs_mapping) = @_;

    if (!defined($mscs_mapping->[0])) {
	return $vertex_pairs;
    }

    my $vertex_mapping_ranges = {};
    for (my $i = 1; $i < @{$mscs_mapping}; $i++) {
	my $mapping_range = [$mscs_mapping->[$i - 1][1] + 1, $mscs_mapping->[$i][1] - 1];
	for (my $j = $mscs_mapping->[$i - 1][0] + 1; $j < $mscs_mapping->[$i][0]; $j++) {
	    $vertex_mapping_ranges->{$j} = $mapping_range;
	}

	$vertex_mapping_ranges->{$mscs_mapping->[$i - 1][0]} = [-1, -1];
	$vertex_mapping_ranges->{$mscs_mapping->[$i][0]} = [-1, -1];
    }

    my $filtered_vertex_pairs = [];

    foreach (@{$vertex_pairs}) {
	if (exists($vertex_mapping_ranges->{$_->[0]})) {
	    my $mapping_range = $vertex_mapping_ranges->{$_->[0]};
	    if ($_->[1] >= $mapping_range->[0] && $_->[1] <= $mapping_range->[1]) {
		push @{$filtered_vertex_pairs}, $_;
	    }
	}
	else {
	    if (($_->[0] < $mscs_mapping->[0][0] && $_->[1] < $mscs_mapping->[0][1]) ||
		($_->[0] > $mscs_mapping->[-1][0] && $_->[1] > $mscs_mapping->[-1][1])) {
		push @{$filtered_vertex_pairs}, $_;
	    }
	}
    }

    return $filtered_vertex_pairs;
}

=comment
sub _filter_non_nest_mapping_candidates {
    my ($candidate_stem_vertex_pairs, $mscs_mapping, $stem_graph1, $stem_graph2) = @_;

    my ($non_nest_stem_partition_map1, $last_partition_end1) = _get_non_nest_stem_partition_map($mscs_mapping, $stem_graph1, 0);
    my ($non_nest_stem_partition_map2, $last_partition_end2) = _get_non_nest_stem_partition_map($mscs_mapping, $stem_graph2, 1);

    my ($in_partition_candidates, $cross_partition_candidates) = ({}, {});
    foreach (@{$candidate_stem_vertex_pairs}) {
	my $vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($_->[0]);
	my $vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($_->[1]);

	my $vertex_upstream_part1 = _get_non_nest_stem_partition($non_nest_stem_partition_map1, $last_partition_end1, $vertex_attrs1->{upstream_start});
	my $vertex_upstream_part2 = _get_non_nest_stem_partition($non_nest_stem_partition_map2, $last_partition_end2, $vertex_attrs2->{upstream_start});
	if ($vertex_upstream_part1 != $vertex_upstream_part2) {
	    next;
	}

	my $vertex_downstream_part1 = _get_non_nest_stem_partition($non_nest_stem_partition_map1, $last_partition_end1, $vertex_attrs1->{downstream_start});
	my $vertex_downstream_part2 = _get_non_nest_stem_partition($non_nest_stem_partition_map2, $last_partition_end2, $vertex_attrs2->{downstream_start});
	if ($vertex_downstream_part1 != $vertex_downstream_part2) {
	    next;
	}

	if ($vertex_upstream_part1 == $vertex_downstream_part1) {
	    my $same_partition_candidates = $in_partition_candidates->{$vertex_upstream_part1};
	    if (defined($same_partition_candidates)) {
		push @{$same_partition_candidates}, $_;
	    }
	    else {
		$in_partition_candidates->{$vertex_upstream_part1} = [$_];
	    }
	}
	else {
	    my $same_cross_candidates = $cross_partition_candidates->{$vertex_upstream_part1 . '-' . $vertex_downstream_part1};
	    if (defined($same_cross_candidates)) {
		push @{$same_cross_candidates}, $_;
	    }
	    else {
		$cross_partition_candidates->{$vertex_upstream_part1 . '-' . $vertex_downstream_part1} = [$_];
	    }
	}
    }

    my $linked_partitions = [];
    while (my ($cross_partition_id, $same_cross_candidates) = each %{$crossed_partition_candidates}) {
	my ($upstream_partition, $downstream_partition) = split(/-/, $cross_partition_id);
	my $upstream_partition_candidates = $in_partition_candidates->{$upstream_partition};
	for (my $i = 0; $i < @{$upstream_partition_candidates}; $i++) {

	}
    }

    return 
}
=cut

=comment
sub _is_in_same_non_nest_stem_partition {
    my ($non_nest_stem_partition_map1, $non_nest_stem_partition_map2, $last_partition_end1, $last_partition_end2, $non_nest_stem_pos1, $non_nest_stem_pos2) = @_;

    my $is_in_last_partition1 = ($non_nest_stem_pos1 > $last_partition_end1);
    my $is_in_last_partition2 = ($non_nest_stem_pos2 > $last_partition_end2);

    if ($is_in_last_partition1 != $is_in_last_partition2) {
	return 0;
    }

    if ($is_in_last_partition1) {
	return 1;
    }

    return ($non_nest_stem_partition_map1->{$non_nest_stem_pos1} == $non_nest_stem_partition_map2->{$non_nest_stem_pos2});
}
=cut

sub _get_mscs {
    my ($candidate_stem_vertex_pairs, $predefined_vertex_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2,
	$base_seq1, $base_seq2) = @_;

    my @sorted_vertex_pairs = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} (@{$predefined_vertex_mapping}, @{$candidate_stem_vertex_pairs});
    my $diff_edge_label_vertex_pairs = _create_diff_edge_label_vertex_pairs(\@sorted_vertex_pairs, $stem_graph1, $stem_graph2);
    my $exclusive_vertex_pair_groups = _group_exclusive_vertex_pairs(\@sorted_vertex_pairs);

    my $mscs_mappings = [];
    my $mscs_candidates = {};
    my ($new_candidates, $excess_candidate_signatures) = ({}, {});

    my $init_mapping_cost = _get_init_mapping_cost(\@sorted_vertex_pairs, $stem_removal_costs1, $stem_removal_costs2);
    my $min_mapping_cost = $init_mapping_cost;

    my $is_merge_stem = 0;
    if (defined($base_seq1) && defined($base_seq2)) {
	$is_merge_stem = 1;
    }

    foreach (@{$exclusive_vertex_pair_groups}) {
	my ($new_candidates, $excess_candidate_signatures) = ({}, {});

	my ($vertex_pair_graph1_vertex, $vertex_pair_graph2_vertices) = ($_->[0], $_->[1]);
	foreach (@{$vertex_pair_graph2_vertices}) {
	    my $vertex_pair = [$vertex_pair_graph1_vertex, $_];

	    if (scalar(keys %{$mscs_candidates}) == 0) {
		$new_candidates->{$vertex_pair->[0] . '-' . $vertex_pair->[1]} = [$vertex_pair];
	    }
	    else {
		while (my ($candidate_signature, $mscs_candidate) = each %{$mscs_candidates}) {
		    my ($tentative_candidate, $tentative_candidate_signature) = _generate_tentative_candidate($vertex_pair, $mscs_candidate, $diff_edge_label_vertex_pairs);
		    $new_candidates->{$tentative_candidate_signature} = $tentative_candidate;
		    if (index($tentative_candidate_signature, $candidate_signature) > -1) {
			$excess_candidate_signatures->{$candidate_signature} = 1;
		    }
		}
	    }
	}

	foreach (keys %{$excess_candidate_signatures}) {
	    delete $mscs_candidates->{$_};
	}

	%{$mscs_candidates} = (%{$mscs_candidates}, %{$new_candidates});
    }

    while (my ($mscs_candidate_signature, $mscs_candidate) = each %{$mscs_candidates}) {
	my $candidate_mapping_cost = _get_candidate_mapping_cost($mscs_candidate, $aligned_stem_pairs, $stem_removal_costs1, $stem_removal_costs2, $init_mapping_cost);

	if ($is_merge_stem) {
	    (my $total_reduced_cost) = StemMatchExtension->extend_stem_match($mscs_candidate, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
	    $candidate_mapping_cost += $total_reduced_cost;
	}

	if ($candidate_mapping_cost < $min_mapping_cost) {
	    $mscs_mappings = [$mscs_candidate];
	    $min_mapping_cost = $candidate_mapping_cost;
	}
	elsif ($candidate_mapping_cost == $min_mapping_cost) {
	    push @{$mscs_mappings}, $mscs_candidate;
	}
    }

    return $mscs_mappings, $min_mapping_cost;
}

=comment
sub _get_mapping_subset {
    my ($mscs_mapping_signature1, $mscs_mapping1, $mscs_mapping_signature2, $mscs_mapping2) = @_;

    my $mscs_vertex_pair_count1 = @{$mscs_mapping1};
    my $mscs_vertex_pair_count2 = @{$mscs_mapping2};

    if ($mscs_vertex_pair_count1 == $mscs_vertex_pair_count2) {
	return;
    }

    if ($mscs_vertex_pair_count1 > $mscs_vertex_pair_count2) {
	foreach (@{$mscs_mapping2}) {
	    if (index($mscs_mapping_signature1, $_->[0] . '-' . $_->[1]) == -1) {
		return;
	    }
	}

	return $mscs_mapping_signature2;
    }
    else {
	foreach (@{$mscs_mapping1}) {
	    if (index($mscs_mapping_signature2, $_->[0] . '-' . $_->[1]) == -1) {
		return;
	    }
	}

	return $mscs_mapping_signature1;
    }
}
=cut

sub _group_exclusive_vertex_pairs {
    my $sorted_vertex_pairs = shift;

    my ($exclusive_vertex_pair_groups, $same_graph1_partner_vertices) = ([], []);
    my $last_graph1_vertex = -1;

    foreach (@{$sorted_vertex_pairs}) {
	if ($_->[0] != $last_graph1_vertex && $last_graph1_vertex >= 0) {
	    push @{$exclusive_vertex_pair_groups}, [$last_graph1_vertex, $same_graph1_partner_vertices];
	    $same_graph1_partner_vertices = [];
	    $last_graph1_vertex = $_->[0];
	}

	push @{$same_graph1_partner_vertices}, $_->[1];
	if ($last_graph1_vertex < 0) {
	    $last_graph1_vertex = $_->[0];
	}
    }

    if ($last_graph1_vertex >= 0) {
	push @{$exclusive_vertex_pair_groups}, [$last_graph1_vertex, $same_graph1_partner_vertices];
    }

    return $exclusive_vertex_pair_groups;
}

sub _create_diff_edge_label_vertex_pairs {
    my ($vertex_pairs, $stem_graph1, $stem_graph2) = @_;

    my $diff_edge_label_vertex_pairs = {};

    for (my $i = 0; $i < @{$vertex_pairs}; $i++) {
	for (my $j = $i + 1; $j < @{$vertex_pairs}; $j++) {
	    if ($vertex_pairs->[$i][0] == $vertex_pairs->[$j][0] || $vertex_pairs->[$i][1] == $vertex_pairs->[$j][1]) {
		next;
	    }

	    if ($stem_graph1->get_edge_label($vertex_pairs->[$i][0], $vertex_pairs->[$j][0]) ne
		$stem_graph2->get_edge_label($vertex_pairs->[$i][1], $vertex_pairs->[$j][1])) {
		$diff_edge_label_vertex_pairs->{$vertex_pairs->[$i][0] . '-' . $vertex_pairs->[$j][0] . '_' .
						    $vertex_pairs->[$i][1] . '-' . $vertex_pairs->[$j][1]} = 1;
	    }
	}
    }

    return $diff_edge_label_vertex_pairs;
}

sub _get_init_mapping_cost {
    my ($selected_vertex_pairs, $stem_removal_costs1, $stem_removal_costs2) = @_;

    my ($included_vertices1, $included_vertices2) = ({}, {});
    my $init_mapping_cost = 0;

    foreach (@{$selected_vertex_pairs}) {
	if (!exists($included_vertices1->{$_->[0]})) {
	    $init_mapping_cost += $stem_removal_costs1->[$_->[0]];
	    $included_vertices1->{$_->[0]} = 1;
	}

	if (!exists($included_vertices2->{$_->[1]})) {
	    $init_mapping_cost += $stem_removal_costs2->[$_->[1]];
	    $included_vertices2->{$_->[1]} = 1;
	}
    }

    return $init_mapping_cost;
}

sub align_stem_pairs {
    my (undef, $stem_vertex_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2) = @_;

    foreach (@{$stem_vertex_pairs}) {
	if (exists($aligned_stem_pairs->{$_->[0] . '-' . $_->[1]})) {
	    next;
	}

	my $vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($_->[0]);
	my $vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($_->[1]);
	my ($aligned_stem_pair, $stem_align_cost1, $stem_align_cost2) = align_stem_pair(undef, $vertex_attrs1, $vertex_attrs2);
	if (defined($stem_align_cost1) && defined($stem_align_cost2)) {
	    print "Inconsistent stem alignment score: $stem_align_cost1 vs $stem_align_cost2 (" . $_->[0] . ', ' . $_->[1] . ")\n";
	}

	$aligned_stem_pairs->{$_->[0] . '-' . $_->[1]} = $aligned_stem_pair;
    }

    return $aligned_stem_pairs;
}

sub align_stem_pair {
    my (undef, $vertex_attrs1, $vertex_attrs2, $appended_seqs, $is_symm_align) = @_;

    my $matched_stem_pair = MatchedStemPair->new($vertex_attrs1, $vertex_attrs2);

    if (defined($appended_seqs)) {
	if (!defined($is_symm_align)) {
	    $is_symm_align = 0;
	}

	if ($appended_seqs->[0][0] ne '') {
	    $matched_stem_pair->append_seq_to_left_end($appended_seqs->[0][0], 1, 'upstream');
	}

	if ($appended_seqs->[0][1] ne '') {
	    $matched_stem_pair->append_seq_to_right_end($appended_seqs->[0][1], 1, 'upstream');
	}

	if ($appended_seqs->[0][2] ne '') {
	    $matched_stem_pair->append_seq_to_left_end($appended_seqs->[0][2], 1, 'downstream');
	}

	if ($appended_seqs->[0][3] ne '') {
	    $matched_stem_pair->append_seq_to_right_end($appended_seqs->[0][3], 1, 'downstream');
	}

	if ($appended_seqs->[1][0] ne '') {
	    $matched_stem_pair->append_seq_to_left_end($appended_seqs->[1][0], 2, 'upstream');
	}

	if ($appended_seqs->[1][1] ne '') {
	    $matched_stem_pair->append_seq_to_right_end($appended_seqs->[1][1], 2, 'upstream');
	}

	if ($appended_seqs->[1][2] ne '') {
	    $matched_stem_pair->append_seq_to_left_end($appended_seqs->[1][2], 2, 'downstream');
	}

	if ($appended_seqs->[1][3] ne '') {
	    $matched_stem_pair->append_seq_to_right_end($appended_seqs->[1][3], 2, 'downstream');
	}
    }
    else {
	if (!defined($is_symm_align)) {
	    $is_symm_align = 1;
	}
    }

    my $stem1_base_seq = $matched_stem_pair->get_base_seq(1, 'upstream') . $matched_stem_pair->get_base_seq(1, 'downstream');
    my $stem2_base_seq = $matched_stem_pair->get_base_seq(2, 'upstream') . $matched_stem_pair->get_base_seq(2, 'downstream');
    my ($stem1_upstream_align_bound, $stem1_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(1);
    my ($stem2_upstream_align_bound, $stem2_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(2);
    my $stem_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $matched_stem_pair->get_pseudo_base_pairs(1),
						       $matched_stem_pair->get_pseudo_base_pairs(2), $stem1_upstream_align_bound,
						       $stem1_downstream_align_bound, $stem2_upstream_align_bound, $stem2_downstream_align_bound);

    my $aligned_stem_pair1 = AlignedStemPair->new($stem_align_result);
    $aligned_stem_pair1->restore_seq_bases($matched_stem_pair->get_base_seq(1, 'upstream'), $matched_stem_pair->get_base_seq(1, 'downstream'),
					   $matched_stem_pair->get_base_seq(2, 'upstream'), $matched_stem_pair->get_base_seq(2, 'downstream'));

    if ($is_symm_align) {
	$stem_align_result = rna_stem_align::align_stem($stem2_base_seq, $stem1_base_seq, $matched_stem_pair->get_pseudo_base_pairs(2),
							$matched_stem_pair->get_pseudo_base_pairs(1), $stem2_upstream_align_bound,
							$stem2_downstream_align_bound, $stem1_upstream_align_bound, $stem1_downstream_align_bound);
	my $aligned_stem_pair2 = AlignedStemPair->new($stem_align_result);
	$aligned_stem_pair2->restore_seq_bases($matched_stem_pair->get_base_seq(2, 'upstream'), $matched_stem_pair->get_base_seq(2, 'downstream'),
					       $matched_stem_pair->get_base_seq(1, 'upstream'), $matched_stem_pair->get_base_seq(1, 'downstream'));
	
	my $stem_align_cost1 = $aligned_stem_pair1->get_cost();
	my $stem_align_cost2 = $aligned_stem_pair2->get_cost();
	if ($stem_align_cost1 == $stem_align_cost2) {
	    if ($aligned_stem_pair2->get_total_del_edge_len() > $aligned_stem_pair1->get_total_del_edge_len()) {
		$aligned_stem_pair1->set_and_swap($aligned_stem_pair2);
	    }
	}
	elsif ($stem_align_cost1 > $stem_align_cost2) {
	    return $aligned_stem_pair1, $stem_align_cost1, $stem_align_cost2;
	}
	else {
	    return $aligned_stem_pair2, $stem_align_cost1, $stem_align_cost2;
	}
    }

    return $aligned_stem_pair1;
}

=comment
sub _align_stem_pair {
    my ($vertex_attrs1, $vertex_attrs2) = @_;

    my $stem1_base_seq = $vertex_attrs1->{upstream_base_seq} . $vertex_attrs1->{downstream_base_seq};
    my $stem2_base_seq = $vertex_attrs2->{upstream_base_seq} . $vertex_attrs2->{downstream_base_seq};
    my $stem1_pseudo_base_pairs = $vertex_attrs1->{pseudo_base_pairs};
    my $stem2_pseudo_base_pairs = $vertex_attrs2->{pseudo_base_pairs};
    my $stem1_downstream_start = $vertex_attrs1->{upstream_length};
    my $stem1_upstream_end = $stem1_downstream_start - 1;
    my $stem2_downstream_start = $vertex_attrs2->{upstream_length};
    my $stem2_upstream_end = $stem2_downstream_start - 1;

    my $stem_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $stem1_pseudo_base_pairs, $stem2_pseudo_base_pairs,
						       $stem1_upstream_end, $stem1_downstream_start, $stem2_upstream_end, $stem2_downstream_start);
    my $aligned_stem_pair1 = AlignedStemPair->new($stem_align_result);

    $stem_align_result = rna_stem_align::align_stem($stem2_base_seq, $stem1_base_seq, $stem2_pseudo_base_pairs, $stem1_pseudo_base_pairs,
						    $stem2_upstream_end, $stem2_downstream_start, $stem1_upstream_end, $stem1_downstream_start);
    my $aligned_stem_pair2 = AlignedStemPair->new($stem_align_result);

    my $stem_align_cost1 = $aligned_stem_pair1->get_cost();
    my $stem_align_cost2 = $aligned_stem_pair2->get_cost();
    if ($stem_align_cost1 == $stem_align_cost2) {
	if ($aligned_stem_pair2->get_total_del_edge_len() > $aligned_stem_pair1->get_total_del_edge_len()) {
	    $aligned_stem_pair1->set_and_swap($aligned_stem_pair2);
	}

	return $aligned_stem_pair1;
    }
    elsif ($stem_align_cost1 > $stem_align_cost2) {
	return $aligned_stem_pair1, $stem_align_cost1, $stem_align_cost2;
    }
    else {
	return $aligned_stem_pair2, $stem_align_cost1, $stem_align_cost2;
    }
}

sub align_stem_with_merged_stem {
    my (undef, $stem_vertex, $stem_graph, $merged_stem_vertex_attrs, $merged_stem_vertex, $merge_graph_index) = @_;

    my $vertex_attrs = $stem_graph->get_vertex_attrs_at($stem_vertex);

    my $stem_base_seq = $vertex_attrs->{upstream_base_seq} . $vertex_attrs->{downstream_base_seq};
    my $merged_stem_base_seq = $merged_stem_vertex_attrs->{upstream_base_seq} . $merged_stem_vertex_attrs->{downstream_base_seq};
    my $stem_pseudo_base_pairs = $vertex_attrs->{pseudo_base_pairs};
    my $merged_stem_pseudo_base_pairs = $merged_stem_vertex_attrs->{pseudo_base_pairs};
    my $stem_downstream_start = $vertex_attrs->{upstream_length};
    my $stem_upstream_end = $stem_downstream_start - 1;
    my $merged_stem_downstream_start = $merged_stem_vertex_attrs->{upstream_length};
    my $merged_stem_upstream_end = $merged_stem_downstream_start - 1;

    my $stem_align_result = rna_stem_align::align_stem($stem_base_seq, $merged_stem_base_seq, $stem_pseudo_base_pairs, $merged_stem_pseudo_base_pairs,
						       $stem_upstream_end, $stem_downstream_start, $merged_stem_upstream_end, $merged_stem_downstream_start);
    my $aligned_stem_pair1 = AlignedStemPair->new($stem_align_result);

    $stem_align_result = rna_stem_align::align_stem($merged_stem_base_seq, $stem_base_seq, $merged_stem_pseudo_base_pairs, $stem_pseudo_base_pairs,
						    $merged_stem_upstream_end, $merged_stem_downstream_start, $stem_upstream_end, $stem_downstream_start);
    my $aligned_stem_pair2 = AlignedStemPair->new($stem_align_result);
    my $stem_align_cost1 = $aligned_stem_pair1->get_cost();
    my $stem_align_cost2 = $aligned_stem_pair2->get_cost();
    if ($stem_align_cost1 != $stem_align_cost2) {
	if ($merge_graph_index = 0) {
	    print "Inconsistent stem alignment score: $stem_align_cost1 vs $stem_align_cost2 (" . $merged_stem_vertex . ' (merged), ' . $stem_vertex . ")\n";
	}
	else {
	    print "Inconsistent stem alignment score: $stem_align_cost1 vs $stem_align_cost2 (" . $stem_vertex . ', ' . $merged_stem_vertex . " (merged))\n";
	}
    }

    if ($aligned_stem_pair2->get_total_del_edge_len() > $aligned_stem_pair1->get_total_del_edge_len()) {
	$aligned_stem_pair1->set_and_swap($aligned_stem_pair2);
    }

    return $aligned_stem_pair1;
}
=cut

sub align_motif_pair {
    my (undef, $non_nest_stem_vertex1, $non_nest_stem_vertex2, $aligned_motif_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $ext_appended_seqs) = @_;

    if (exists($aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2})) {
	return $aligned_motif_pairs;
    }

    my $non_nest_stem_vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($non_nest_stem_vertex1);
    my $motif1_loop_seq = get_bounded_loop_seq(undef, $non_nest_stem_vertex_attrs1, $base_seq1);

    my $non_nest_stem_vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($non_nest_stem_vertex2);
    my $motif2_loop_seq = get_bounded_loop_seq(undef, $non_nest_stem_vertex_attrs2, $base_seq2);

    my $motif1_size = $non_nest_stem_vertex_attrs1->{downstream_start} + $non_nest_stem_vertex_attrs1->{downstream_length} - $non_nest_stem_vertex_attrs1->{upstream_start};
    my $motif2_size = $non_nest_stem_vertex_attrs2->{downstream_start} + $non_nest_stem_vertex_attrs2->{downstream_length} - $non_nest_stem_vertex_attrs2->{upstream_start};

    if ($motif1_size >= DefaultAlignParams->MAX_NON_NEST_MOTIF_SIZE_FOR_EXACT_ALIGN || $motif2_size >= DefaultAlignParams->MAX_NON_NEST_MOTIF_SIZE_FOR_EXACT_ALIGN) {
	my ($aligned_stem_pair) = align_stem_pair(undef, $non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2, $ext_appended_seqs);
	my $expected_start_edge_gap_size = get_expected_edge_gap_size_from_stem(undef, $aligned_stem_pair, 'upstream', 'right');
	my $expected_end_edge_gap_size = get_expected_edge_gap_size_from_stem(undef, $aligned_stem_pair, 'downstream', 'left');

	my $seq_align_result = SequenceDP->align($motif1_loop_seq, $motif2_loop_seq, $expected_start_edge_gap_size, $expected_end_edge_gap_size);
	my $aligned_seq_pair = AlignedSequencePair->new($seq_align_result);

#	my $matched_stem_pair = MatchedStemPair->new($non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2);
	my $int_appended_seqs = [['', '', '', ''], ['', '', '', '']];
	my $is_bad_align = 0;

	if ($expected_start_edge_gap_size > 0) {
	    if ($aligned_seq_pair->get_del_edge_len(2, 'left') > 0) {
=comment
		my $aligned_seq_edge_gap_cutoff_len = $aligned_seq_pair->get_del_edge_len(2, 'left');
		my $eff_aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(1, 'upstream', 'right'),
												  $motif2_loop_seq, 'left');
		if ($eff_aligned_seq_edge_gap_cutoff_len < $aligned_seq_edge_gap_cutoff_len) {
		    $aligned_seq_edge_gap_cutoff_len = $eff_aligned_seq_edge_gap_cutoff_len;
		}
=cut
#		my $del_seq = $aligned_seq_pair->trim_left_edge_gap(DefaultAlignParams->BASE_REMOVAL_COST);
#		$matched_stem_pair->append_seq_to_right_end($del_seq, 2, 'upstream');

#		$int_appended_seqs->[1][1] = $aligned_seq_pair->trim_left_edge_gap();

		my $aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(1, 'upstream', 'right'),
											      $aligned_seq_pair->get_del_edge_seq(2, 'left'), 'left');
		$int_appended_seqs->[1][1] = $aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);

		$is_bad_align = 1;
	    }
	}
	elsif ($expected_start_edge_gap_size < 0) {
	    if ($aligned_seq_pair->get_del_edge_len(1, 'left') > 0) {
=comment
		my $aligned_seq_edge_gap_cutoff_len = $aligned_seq_pair->get_del_edge_len(1, 'left');
		my $eff_aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(2, 'upstream', 'right'),
												  $motif1_loop_seq, 'left');
		if ($eff_aligned_seq_edge_gap_cutoff_len < $aligned_seq_edge_gap_cutoff_len) {
		    $aligned_seq_edge_gap_cutoff_len = $eff_aligned_seq_edge_gap_cutoff_len;
		}
=cut
#		my $del_seq = $aligned_seq_pair->trim_left_edge_gap(DefaultAlignParams->BASE_REMOVAL_COST);
#		$matched_stem_pair->append_seq_to_right_end($del_seq, 1, 'upstream');

#		$int_appended_seqs->[0][1] = $aligned_seq_pair->trim_left_edge_gap();

		my $aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(2, 'upstream', 'right'),
											      $aligned_seq_pair->get_del_edge_seq(1, 'left'), 'left');
		$int_appended_seqs->[0][1] = $aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);

		$is_bad_align = 1;
	    }
	}

	if ($expected_end_edge_gap_size > 0) {
	    if ($aligned_seq_pair->get_del_edge_len(2, 'right') > 0) {
=comment
		my $aligned_seq_edge_gap_cutoff_len = $aligned_seq_pair->get_del_edge_len(2, 'right');
		my $eff_aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(1, 'downstream', 'left'),
												  $motif2_loop_seq, 'right');
		if ($eff_aligned_seq_edge_gap_cutoff_len < $aligned_seq_edge_gap_cutoff_len) {
		    $aligned_seq_edge_gap_cutoff_len = $eff_aligned_seq_edge_gap_cutoff_len;
		}
=cut
#		my $del_seq = $aligned_seq_pair->trim_right_edge_gap(DefaultAlignParams->BASE_REMOVAL_COST);
#		$matched_stem_pair->append_seq_to_left_end($del_seq, 2, 'downstream');

#		$int_appended_seqs->[1][2] = $aligned_seq_pair->trim_right_edge_gap();

		my $aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(1, 'downstream', 'left'),
											      $aligned_seq_pair->get_del_edge_seq(2, 'right'), 'right');
		$int_appended_seqs->[1][2] = $aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		$is_bad_align = 1;
	    }
	}
	elsif ($expected_end_edge_gap_size < 0) {
	    if ($aligned_seq_pair->get_del_edge_len(1, 'right') > 0) {
=comment
		my $aligned_seq_edge_gap_cutoff_len = $aligned_seq_pair->get_del_edge_len(1, 'right');
		my $eff_aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(2, 'downstream', 'left'),
												  $motif1_loop_seq, 'right');
		if ($eff_aligned_seq_edge_gap_cutoff_len < $aligned_seq_edge_gap_cutoff_len) {
		    $aligned_seq_edge_gap_cutoff_len = $eff_aligned_seq_edge_gap_cutoff_len;
		}
=cut
#		my $del_seq = $aligned_seq_pair->trim_right_edge_gap(DefaultAlignParams->BASE_REMOVAL_COST);
#		$matched_stem_pair->append_seq_to_left_end($del_seq, 1, 'downstream');

#		$int_appended_seqs->[0][2] = $aligned_seq_pair->trim_right_edge_gap();

		my $aligned_seq_edge_gap_cutoff_len = get_eff_aligned_seq_edge_gap_cutoff_len(undef, $aligned_stem_pair->get_del_edge_seq(2, 'downstream', 'left'),
											      $aligned_seq_pair->get_del_edge_seq(1, 'right'), 'right');
		$int_appended_seqs->[0][2] = $aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);

		$is_bad_align = 1;
	    }
	}

	if ($is_bad_align) {
	    my $appended_seqs;
	    if (defined($ext_appended_seqs)) {
		$appended_seqs = [[$ext_appended_seqs->[0][0], $int_appended_seqs->[0][1], $int_appended_seqs->[0][2], $ext_appended_seqs->[0][3]],
				 [$ext_appended_seqs->[1][0], $int_appended_seqs->[1][1], $int_appended_seqs->[1][2], $ext_appended_seqs->[1][3]]];
	    }
	    else {
		$appended_seqs = [['', $int_appended_seqs->[0][1], $int_appended_seqs->[0][2], ''],
				 ['', $int_appended_seqs->[1][1], $int_appended_seqs->[1][2], '']];
	    }

	    my ($realigned_stem_pair) = align_stem_pair(undef, $non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2, $appended_seqs);

	    $aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2} = AlignedMotifPair->new($realigned_stem_pair, $aligned_seq_pair);
=comment
	    my $stem1_base_seq = $matched_stem_pair->get_base_seq(1, 'upstream') . $matched_stem_pair->get_base_seq(1, 'downstream');
	    my $stem2_base_seq = $matched_stem_pair->get_base_seq(2, 'upstream') . $matched_stem_pair->get_base_seq(2, 'downstream');
	    my ($stem1_upstream_align_bound, $stem1_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(1);
	    my ($stem2_upstream_align_bound, $stem2_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(2);

	    my $stem_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $matched_stem_pair->get_pseudo_base_pairs(1),
							       $matched_stem_pair->get_pseudo_base_pairs(2), $stem1_upstream_align_bound,
							       $stem1_downstream_align_bound, $stem2_upstream_align_bound, $stem2_downstream_align_bound);

	    my $realigned_stem_pair = AlignedStemPair->new($stem_align_result);

	    $aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2} = AlignedMotifPair->new_large_motif($realigned_stem_pair, $aligned_seq_pair);
=cut
	}
	else {
	    $aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2} = AlignedMotifPair->new($aligned_stem_pair, $aligned_seq_pair);
	}
    }
    else {
	my ($appended_seqs, $is_symm_align);
	if (defined($ext_appended_seqs)) {
	    $appended_seqs = [[$ext_appended_seqs->[0][0], $motif1_loop_seq, '', $ext_appended_seqs->[0][3]],
			      [$ext_appended_seqs->[1][0], $motif2_loop_seq, '', $ext_appended_seqs->[1][3]]];
	    $is_symm_align = 0;
	}
	else {
	    $appended_seqs = [['', $motif1_loop_seq, '', ''],
			      ['', $motif2_loop_seq, '', '']];
	    $is_symm_align = 1;
	}

	my ($aligned_stem_pair) = align_stem_pair(undef, $non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2, $appended_seqs, $is_symm_align);
	$aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2} = AlignedMotifPair->new($aligned_stem_pair);
=comment
	my $matched_stem_pair = MatchedStemPair->new($non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2);
	$matched_stem_pair->append_seq_to_right_end($motif1_loop_seq, 1, 'upstream');
	$matched_stem_pair->append_seq_to_right_end($motif2_loop_seq, 2, 'upstream');

	my $stem1_base_seq = $matched_stem_pair->get_base_seq(1, 'upstream') . $matched_stem_pair->get_base_seq(1, 'downstream');
	my $stem2_base_seq = $matched_stem_pair->get_base_seq(2, 'upstream') . $matched_stem_pair->get_base_seq(2, 'downstream');
	my ($stem1_upstream_align_bound, $stem1_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(1);
	my ($stem2_upstream_align_bound, $stem2_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(2);

	my $motif_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $matched_stem_pair->get_pseudo_base_pairs(1), $matched_stem_pair->get_pseudo_base_pairs(2),
							   $stem1_upstream_align_bound, $stem1_downstream_align_bound, $stem2_upstream_align_bound, $stem2_downstream_align_bound);

	$aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2} = AlignedMotifPair->new($motif_align_result);
=cut
    }

    return $aligned_motif_pairs;
}

sub get_bounded_loop_seq {
    my (undef, $vertex_attrs, $base_seq) = @_;

    my $stem_base_pairs = $vertex_attrs->{base_pairs};

    return substr($base_seq, $stem_base_pairs->[-1][0] + 1, ($stem_base_pairs->[-1][1] - $stem_base_pairs->[-1][0] - 1));
}

sub get_expected_edge_gap_size_from_stem {
    my (undef, $aligned_stem_pair, $stream_opt, $end_opt) = @_;

    my $expected_align_seq_edge_gap_size = 0;

    my $stem1_del_edge_len = $aligned_stem_pair->get_del_edge_len(1, $stream_opt, $end_opt);
    my $stem2_del_edge_len = $aligned_stem_pair->get_del_edge_len(2, $stream_opt, $end_opt);
    if ($stem1_del_edge_len >= $stem2_del_edge_len) {
	$expected_align_seq_edge_gap_size = $stem1_del_edge_len;
    }
    else {
	$expected_align_seq_edge_gap_size = $stem2_del_edge_len * -1;
    }

    return $expected_align_seq_edge_gap_size;
}

sub get_expected_edge_gap_size_from_motif {
    my (undef, $aligned_motif_pair, $end_opt) = @_;

    my $expected_align_seq_edge_gap_size = 0;

    my $stem1_del_edge_len = $aligned_motif_pair->get_del_edge_len(1, $end_opt);
    my $stem2_del_edge_len = $aligned_motif_pair->get_del_edge_len(2, $end_opt);
    if ($stem1_del_edge_len >= $stem2_del_edge_len) {
	$expected_align_seq_edge_gap_size = $stem1_del_edge_len;
    }
    else {
	$expected_align_seq_edge_gap_size = $stem2_del_edge_len * -1;
    }

    return $expected_align_seq_edge_gap_size;
}

sub get_eff_aligned_seq_edge_gap_cutoff_len {
    my (undef, $stem_del_edge_seq, $base_seq, $end_opt) = @_;

    if ($end_opt eq 'left') {
	my $seq_align_result = SequenceDP->align($stem_del_edge_seq, $base_seq, 0, 1);
	my $aligned_seq_pair = AlignedSequencePair->new($seq_align_result);
	my $seq_del_edge_len = $aligned_seq_pair->get_del_edge_len(2, 'right');

	return length($base_seq) - $seq_del_edge_len;
    }
    elsif ($end_opt eq 'right') {
	my $seq_align_result = SequenceDP->align($stem_del_edge_seq, $base_seq, 1, 0);
	my $aligned_seq_pair = AlignedSequencePair->new($seq_align_result);
	my $seq_del_edge_len = $aligned_seq_pair->get_del_edge_len(2, 'left');

	return length($base_seq) - $seq_del_edge_len;
    }

    return 0;
}

sub _generate_tentative_candidate {
    my ($vertex_pair, $solution_candidate, $diff_edge_label_vertex_pairs) = @_;

    my $tentative_candidate = [];
    my $tentative_candidate_signature = '';

    foreach (@{$solution_candidate}) {
	if ($_->[1] >= $vertex_pair->[1]) {
	    last;
	}

	if (exists($diff_edge_label_vertex_pairs->{$_->[0] . '-' . $vertex_pair->[0] . '_' . $_->[1] . '-' . $vertex_pair->[1]})) {
	    next;
	}

	push @{$tentative_candidate}, $_;
	$tentative_candidate_signature = $tentative_candidate_signature . $_->[0] . '-' . $_->[1] . '_';
    }

    push @{$tentative_candidate}, $vertex_pair;
    $tentative_candidate_signature = $tentative_candidate_signature . $vertex_pair->[0] . '-' . $vertex_pair->[1];

    return $tentative_candidate, $tentative_candidate_signature;
}

sub _get_candidate_mapping_cost {
    my ($candidate, $aligned_stem_pairs, $stem_removal_costs1, $stem_removal_costs2, $init_mapping_cost) = @_;

    my $mapping_cost = $init_mapping_cost;

    foreach (@{$candidate}) {
	$mapping_cost += (($aligned_stem_pairs->{$_->[0] . '-' . $_->[1]})->get_cost() - $stem_removal_costs1->[$_->[0]] -
			  $stem_removal_costs2->[$_->[1]]);
    }

    return $mapping_cost;
}

=comment
sub expand_non_nest_stem_vertex_align_trend {
    my (undef, $candidate_stem_vertex_pairs, $align_trend, $motif_mapping_cost_ratios, $stem_graph1, $stem_graph2) = @_;

    $candidate_stem_vertex_pairs = _filter_mapping_candidates($candidate_stem_vertex_pairs, $align_trend, $stem_graph1, $stem_graph2);
    my $init_total_mapping_cost_ratio = _get_init_total_mapping_cost_ratio($candidate_stem_vertex_pairs);

    my ($align_trend_extensions, $extension_mapping_cost_ratio) =
	_get_non_nest_stem_vertex_align_trends($candidate_stem_vertex_pairs, $motif_mapping_cost_ratios, $stem_graph1, $stem_graph2, $init_total_mapping_cost_ratio);

    my @expanded_align_trends;
    foreach (@{$align_trend_extensions}) {
	my @expanded_align_trend = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} (@{$align_trend}, @{$_});
	push @expanded_align_trends, \@expanded_align_trend;
    }

    if (@expanded_align_trends == 0) {
	push @expanded_align_trends, $align_trend;
    }

    return \@expanded_align_trends, ($extension_mapping_cost_ratio - $init_total_mapping_cost_ratio), $motif_mapping_cost_ratios;
}

sub _get_non_nest_stem_vertex_align_trends {
    my ($candidate_stem_vertex_pairs, $motif_mapping_cost_ratios, $stem_graph1, $stem_graph2, $init_total_mapping_cost_ratio) = @_;

    my @sorted_vertex_pairs = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$candidate_stem_vertex_pairs};
    my $diff_edge_label_vertex_pairs = _create_diff_edge_label_vertex_pairs(\@sorted_vertex_pairs, $stem_graph1, $stem_graph2);
    my $exclusive_vertex_pair_groups = _group_exclusive_vertex_pairs(\@sorted_vertex_pairs);

    my $align_trends = [];
    my $align_trend_candidates = {};
    my ($new_candidates, $excess_candidate_signatures) = ({}, {});
    my $min_mapping_cost_ratio = $init_total_mapping_cost_ratio;
    my $is_align_trend_candidates_empty = 1;

    foreach (@{$exclusive_vertex_pair_groups}) {
	my ($new_candidates, $excess_candidate_signatures) = ({}, {});

	my ($stem_vertex1, $stem_vertices2) = ($_->[0], $_->[1]);
	foreach (@{$stem_vertices2}) {
	    my $vertex_pair = [$stem_vertex1, $_];

	    if ($is_align_trend_candidates_empty) {
		$new_candidates->{$vertex_pair->[0] . '-' . $vertex_pair->[1]} = [$vertex_pair];
	    }
	    else {
		while (my ($candidate_signature, $align_trend_candidate) = each %{$align_trend_candidates}) {
		    my ($tentative_candidate, $tentative_candidate_signature, $is_direct_expansion) =
			_generate_tentative_candidate($vertex_pair, $align_trend_candidate, $diff_edge_label_vertex_pairs);
		    $new_candidates->{$tentative_candidate_signature} = $tentative_candidate;
		    if ($is_direct_expansion) {
			$excess_candidate_signatures->{$candidate_signature} = 1;
		    }
		}
	    }
	}

	foreach (keys %{$excess_candidate_signatures}) {
	    delete $align_trend_candidates->{$_};
	}

	while (my ($new_candidate_signature, $new_candidate) = each %{$new_candidates}) {
	    my $candidate_total_mapping_cost_ratio = _get_candidate_total_mapping_cost_ratio($new_candidate, $motif_mapping_cost_ratios, $init_total_mapping_cost_ratio);
	    if ($candidate_total_mapping_cost_ratio < $min_mapping_cost_ratio) {
		$align_trends = [$new_candidate];
		$min_mapping_cost_ratio = $candidate_total_mapping_cost_ratio;
	    }
	    elsif ($candidate_total_mapping_cost_ratio == $min_mapping_cost_ratio) {
		push @{$align_trends}, $new_candidate;
	    }

	    $align_trend_candidates->{$new_candidate_signature} = $new_candidate;
	}

	$is_align_trend_candidates_empty = 0;
    }

    return $align_trends, $min_mapping_cost_ratio;
}

sub _get_candidate_total_mapping_cost_ratio {
    my ($candidate, $motif_mapping_cost_ratios, $init_total_mapping_cost_ratio) = @_;

    my $total_mapping_cost_ratio = $init_total_mapping_cost_ratio;

    foreach (@{$candidate}) {
	$total_mapping_cost_ratio += ($motif_mapping_cost_ratios->{$_->[0] . '-' . $_->[1]} - 1);
    }

    return $total_mapping_cost_ratio;
}

sub _get_init_total_mapping_cost_ratio {
    my $selected_vertex_pairs = shift;

    my ($included_vertices1, $included_vertices2) = ({}, {});
    my $init_total_mapping_cost_ratio = 0;

    foreach (@{$selected_vertex_pairs}) {
	if (!exists($included_vertices1->{$_->[0]})) {
	    $init_total_mapping_cost_ratio += 1;
	    $included_vertices1->{$_->[0]} = 1;
	}

	if (!exists($included_vertices2->{$_->[1]})) {
	    $init_total_mapping_cost_ratio += 1;
	    $included_vertices2->{$_->[1]} = 1;
	}
    }

    return $init_total_mapping_cost_ratio;
}
=cut

1;
