package StemAlignment;

use strict;

sub align {
    my (undef, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $align_params) = @_;

    if ($align_params->is_no_progressive_stem_align()) {
	return _align_all_stems($stem_graph1, $stem_graph2, {}, $base_seq1, $base_seq2);
    }
    else {
	my ($mscs_mappings, $aligned_stem_pairs, $aligned_motif_pairs) = _align_by_progressive_stem_match($stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $align_params);
	if (defined($mscs_mappings->[0])) {
	    return $mscs_mappings, $aligned_stem_pairs, $aligned_motif_pairs;
	}
	elsif (!$align_params->is_progressive_stem_align_only()) {
	    return _align_all_stems($stem_graph1, $stem_graph2, $aligned_stem_pairs, $base_seq1, $base_seq2);
	}
	else {
	    return [], {};
	}
    }
}

sub _align_by_progressive_stem_match {
    my ($stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $align_params) = @_;

    my $stem_removal_costs1 = _precalculate_stem_removal_costs($stem_graph1);
    my $stem_removal_costs2 = _precalculate_stem_removal_costs($stem_graph2);

    my $nest_mscs_mappings = [];
    my $aligned_stem_pairs = {};

    my ($best_non_nest_stem_vertex_pairs, $aligned_motif_pairs, $motif_removal_costs1, $motif_removal_costs2) =
	NonNestMotifProbeAlign->match_non_nest_motifs($stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2, $align_params);

    if (defined($best_non_nest_stem_vertex_pairs->[0])) {
	($nest_mscs_mappings, $aligned_stem_pairs) = SelectiveGraphMatch->match_nest_stem_vertices($best_non_nest_stem_vertex_pairs, $stem_graph1, $stem_graph2,
												   $stem_removal_costs1, $stem_removal_costs2, $align_params);
    }

    my $mscs_mappings = [];
    my $final_mscs_mapping_cost = -1;
#    my $max_reduced_mscs_cost = 0;

    foreach my $nest_mscs_mapping (@{$nest_mscs_mappings}) {
#	foreach my $t1 (@{$nest_mscs_mapping}) {
#	    print '[' . $t1->[0] . ', ' . $t1->[1] . '] ';
#	}
#	print "\n";

	my $expanded_mscs_mappings;
	if (defined($nest_mscs_mapping->[0])) {
	    ($expanded_mscs_mappings, $aligned_stem_pairs) =
		_expand_mscs_mapping_with_non_nest_vertices($nest_mscs_mapping, $aligned_motif_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1,
							    $motif_removal_costs2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);
	}
	else {
	    $expanded_mscs_mappings = [$best_non_nest_stem_vertex_pairs];
	}

	foreach (@{$expanded_mscs_mappings}) {
	    $aligned_stem_pairs = GraphMatchTools->align_stem_pairs($_, $aligned_stem_pairs, $stem_graph1, $stem_graph2);
	    (my $final_mscs_mappings, $aligned_stem_pairs, my $restored_mapping_cost) =
#		_restore_missing_nest_stem_vertex_mscs_mapping($_, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2,
#							       $base_seq1, $base_seq2);
		_restore_missing_mscs_pairs($_, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);

=comment
	    if ($reduced_mscs_cost == 0) {
		($reduced_mscs_cost) = StemMatchExtension->extend_stem_match($_, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
	    }

	    if ($reduced_mscs_cost < $max_reduced_mscs_cost) {
		$mscs_mappings = $final_mscs_mappings;
		$max_reduced_mscs_cost = $reduced_mscs_cost;
	    }
	    elsif ($reduced_mscs_cost == $max_reduced_mscs_cost) {
		push @{$mscs_mappings}, @{$final_mscs_mappings};
	    }
=cut

	    if ($restored_mapping_cost < $final_mscs_mapping_cost) {
		$mscs_mappings = [];
	    }
	    elsif ($final_mscs_mapping_cost > -1 && $restored_mapping_cost > $final_mscs_mapping_cost) {
		next;
	    }

	    push @{$mscs_mappings}, @{$final_mscs_mappings};
	    $final_mscs_mapping_cost = $restored_mapping_cost;
	}
    }

    my $sorted_mscs_mappings = [];
    foreach (@{$mscs_mappings}) {
	my @sorted_mscs_mapping = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$_};
	push @{$sorted_mscs_mappings}, \@sorted_mscs_mapping;
    }

    return $sorted_mscs_mappings, $aligned_stem_pairs, $aligned_motif_pairs;
}

#sub _restore_missing_nest_stem_vertex_mscs_mapping {
sub _restore_missing_mscs_pairs {
    my ($mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2) = @_;

=comment
    my %unmapped_nest_stem_vertices1 = map {$_ => 1} @{$stem_graph1->get_nest_stem_vertices()};
    my %unmapped_nest_stem_vertices2 = map {$_ => 1} @{$stem_graph2->get_nest_stem_vertices()};

    foreach (@{$mscs_mapping}) {
	delete $unmapped_nest_stem_vertices1{$_->[0]};
	delete $unmapped_nest_stem_vertices2{$_->[1]};
    }

    my @unmapped_vertex_ids1 = keys %unmapped_nest_stem_vertices1;
    my @unmapped_vertex_ids2 = keys %unmapped_nest_stem_vertices2;
=cut

    my %unmapped_stem_vertices1 = (map {$_ => 1} @{$stem_graph1->get_nest_stem_vertices()}, map {$_ => 1} @{$stem_graph1->get_non_nest_stem_vertices()});
    my %unmapped_stem_vertices2 = (map {$_ => 1} @{$stem_graph2->get_nest_stem_vertices()}, map {$_ => 1} @{$stem_graph2->get_non_nest_stem_vertices()});

    foreach (@{$mscs_mapping}) {
	delete $unmapped_stem_vertices1{$_->[0]};
	delete $unmapped_stem_vertices2{$_->[1]};
    }

    my @unmapped_vertex_ids1 = sort {$a <=> $b} keys %unmapped_stem_vertices1;
    my @unmapped_vertex_ids2 = sort {$a <=> $b} keys %unmapped_stem_vertices2;

    my $candidate_vertex_pairs = GraphMatchTools->generate_unmapped_stem_vertex_pairs(\@unmapped_vertex_ids1, \@unmapped_vertex_ids2, $mscs_mapping,
										      $stem_graph1, $stem_graph2, 0);
    my $missing_vertex_pair_candidate_partitions = _partition_missing_vertex_pair_candidates($candidate_vertex_pairs, \@unmapped_vertex_ids1, \@unmapped_vertex_ids2,
											     $stem_graph1, $stem_graph2);
    my $filtered_candidate_vertex_pairs = [];

    foreach (@{$missing_vertex_pair_candidate_partitions}) {
	if (@{$_} == 1) {
	    push @{$filtered_candidate_vertex_pairs}, $_->[0];
	}
	else {
	    (my $local_recovered_ecgm_mappings, $aligned_stem_pairs) =
		GraphMatchTools->expand_mscs_mapping($_, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2,
						     $base_seq1, $base_seq2);

	    my $local_recovered_ecgm_vertex_pairs = {};
	    foreach my $local_recovered_ecgm_mapping (@{$local_recovered_ecgm_mappings}) {
		foreach my $recovered_ecgm_vertex_pair (@{$local_recovered_ecgm_mapping}) {
		    if (exists($unmapped_stem_vertices1{$recovered_ecgm_vertex_pair->[0]}) || exists($unmapped_stem_vertices2{$recovered_ecgm_vertex_pair->[1]})) {
			$local_recovered_ecgm_vertex_pairs->{$recovered_ecgm_vertex_pair->[0] . '-' . $recovered_ecgm_vertex_pair->[1]} = $recovered_ecgm_vertex_pair;
		    }
		}
	    }

	    @{$filtered_candidate_vertex_pairs} = (@{$filtered_candidate_vertex_pairs}, values %{$local_recovered_ecgm_vertex_pairs});
	}
    }

    return GraphMatchTools->expand_mscs_mapping($filtered_candidate_vertex_pairs, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2,
						$stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);
}

sub _partition_missing_vertex_pair_candidates {
    my ($candidate_vertex_pairs, $unmapped_vertex_ids1, $unmapped_vertex_ids2, $stem_graph1, $stem_graph2) = @_;

    my ($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map) = ({}, {}, {});
    my $partition_count = 0;

    my @sorted_candidate_vertex_pairs = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$candidate_vertex_pairs};
    foreach (@sorted_candidate_vertex_pairs) {
	my $graph1_vertex_partition_id = $graph1_vertex_to_partition_map->{$_->[0]};
	my $graph2_vertex_partition_id = $graph2_vertex_to_partition_map->{$_->[1]};

	if (defined($graph1_vertex_partition_id)) {
	    if (defined($graph2_vertex_partition_id)) {
		if ($graph1_vertex_partition_id == $graph2_vertex_partition_id) {
		    push @{$candidate_partitions->{$graph1_vertex_partition_id}}, $_;
		}
		else {
		    ($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map) =
			_merge_missing_vertex_pair_candidate_partitions($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map,
									$graph1_vertex_partition_id, $graph2_vertex_partition_id);
		}
	    }
	    else {
		push @{$candidate_partitions->{$graph1_vertex_partition_id}}, $_;
		$graph2_vertex_to_partition_map->{$_->[1]} = $graph1_vertex_partition_id;
	    }
	}
	else {
	    if (defined($graph2_vertex_partition_id)) {
		push @{$candidate_partitions->{$graph2_vertex_partition_id}}, $_;
		$graph1_vertex_to_partition_map->{$_->[0]} = $graph2_vertex_partition_id;
	    }
	    else {
		$candidate_partitions->{$partition_count} = [$_];
		$graph1_vertex_to_partition_map->{$_->[0]} = $partition_count;
		$graph2_vertex_to_partition_map->{$_->[1]} = $partition_count;
		$partition_count++;
	    }
	}
    }

    my ($unmapped_graph1_crossing_vertices, $unmapped_graph2_crossing_vertices) = ([], []);
    foreach (@{$unmapped_vertex_ids1}) {
	if (($stem_graph1->get_stem_type_at($_) eq 'T') && exists($graph1_vertex_to_partition_map->{$_})) {
	    push @{$unmapped_graph1_crossing_vertices}, $_;
	}
    }

    foreach (@{$unmapped_vertex_ids2}) {
	if (($stem_graph2->get_stem_type_at($_) eq 'T') && exists($graph2_vertex_to_partition_map->{$_})) {
	    push @{$unmapped_graph2_crossing_vertices}, $_;
	}
    }

    for (my $i = 0; $i < @{$unmapped_graph1_crossing_vertices}; $i++) {
	my $merge_partition_id = $graph1_vertex_to_partition_map->{$unmapped_graph1_crossing_vertices->[$i]};
	for (my $j = $i + 1; $j < @{$unmapped_graph1_crossing_vertices}; $j++) {
	    if ($stem_graph1->get_edge_label($unmapped_graph1_crossing_vertices->[$i], $unmapped_graph1_crossing_vertices->[$j]) eq 'K') {
		my $delete_partition_id = $graph1_vertex_to_partition_map->{$unmapped_graph1_crossing_vertices->[$j]};
		if ($merge_partition_id == $delete_partition_id) {
		    next;
		}

		($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map) =
		    _merge_missing_vertex_pair_candidate_partitions($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map,
								    $merge_partition_id, $delete_partition_id);
	    }
	}
    }

    for (my $i = 0; $i < @{$unmapped_graph2_crossing_vertices}; $i++) {
	my $merge_partition_id = $graph2_vertex_to_partition_map->{$unmapped_graph2_crossing_vertices->[$i]};
	for (my $j = $i + 1; $j < @{$unmapped_graph2_crossing_vertices}; $j++) {
	    if ($stem_graph2->get_edge_label($unmapped_graph2_crossing_vertices->[$i], $unmapped_graph2_crossing_vertices->[$j]) eq 'K') {
		my $delete_partition_id = $graph2_vertex_to_partition_map->{$unmapped_graph2_crossing_vertices->[$j]};
		if ($merge_partition_id == $delete_partition_id) {
		    next;
		}

		($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map) =
		    _merge_missing_vertex_pair_candidate_partitions($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map,
								    $merge_partition_id, $delete_partition_id);
	    }
	}
    }

    my @missing_vertex_pair_candidate_partitions = values %{$candidate_partitions};

    return \@missing_vertex_pair_candidate_partitions;
}

sub _merge_missing_vertex_pair_candidate_partitions {
    my ($candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map, $merge_partition_id, $delete_partition_id) = @_;

    my $partition_to_merge = $candidate_partitions->{$merge_partition_id};
    my $partition_to_delete = $candidate_partitions->{$delete_partition_id};

    my @merged_partitions = (@{$partition_to_merge}, @{$partition_to_delete});
    $candidate_partitions->{$merge_partition_id} = \@merged_partitions;

    foreach (@{$partition_to_delete}) {
	$graph1_vertex_to_partition_map->{$_->[0]} = $merge_partition_id;
	$graph2_vertex_to_partition_map->{$_->[1]} = $merge_partition_id;
    }

    delete $candidate_partitions->{$delete_partition_id};

    return $candidate_partitions, $graph1_vertex_to_partition_map, $graph2_vertex_to_partition_map;
}

sub _expand_mscs_mapping_with_non_nest_vertices {
    my ($nest_mscs_mapping, $aligned_motif_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2, $stem_removal_costs1,
	$stem_removal_costs2, $base_seq1, $base_seq2) = @_;

    my $independent_partition_candidate_vertex_pairs =
	_generate_independent_partition_candidate_vertex_pairs($nest_mscs_mapping, $aligned_motif_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2,
							       $motif_removal_costs1, $motif_removal_costs2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);

    my $updated_mscs_mappings = [$nest_mscs_mapping];

    my $all_partition_mscs_mappings = [];
    foreach my $partition_candidate_vertex_pairs (@{$independent_partition_candidate_vertex_pairs}) {
	my $non_redundant_part_mscs_mappings = {};

	foreach (@{$partition_candidate_vertex_pairs}) {
	    $aligned_motif_pairs = GraphMatchTools->align_motif_pair($_->[0], $_->[1], $aligned_motif_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
	}

	my $partition_mscs_mappings = GraphMatchTools->expand_non_nest_motif_mscs_mapping($partition_candidate_vertex_pairs, [], $aligned_motif_pairs, $stem_graph1,
											  $stem_graph2, $motif_removal_costs1, $motif_removal_costs2);
	$non_redundant_part_mscs_mappings = _add_partition_mscs_mappings($non_redundant_part_mscs_mappings, $partition_mscs_mappings);

	($partition_mscs_mappings, $aligned_stem_pairs) =
	    GraphMatchTools->expand_mscs_mapping($partition_candidate_vertex_pairs, [], $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2);
	$non_redundant_part_mscs_mappings = _add_partition_mscs_mappings($non_redundant_part_mscs_mappings, $partition_mscs_mappings);

	@{$partition_mscs_mappings} = values %{$non_redundant_part_mscs_mappings};

	push @{$all_partition_mscs_mappings}, $partition_mscs_mappings;
    }

    my $updated_mscs_mappings = [$nest_mscs_mapping];
    foreach my $partition_mscs_mappings (@{$all_partition_mscs_mappings}) {
	my $expanded_mscs_mappings = [];
	foreach my $partition_mscs_mapping (@{$partition_mscs_mappings}) {
	    foreach my $updated_mscs_mapping (@{$updated_mscs_mappings}) {
		my @expanded_mscs_mapping = (@{$partition_mscs_mapping}, @{$updated_mscs_mapping});
		push @{$expanded_mscs_mappings}, \@expanded_mscs_mapping;
	    }
	}

	$updated_mscs_mappings = $expanded_mscs_mappings;
    }

    for (my $i = 0; $i < @{$updated_mscs_mappings}; $i++) {
	my @sorted_mscs_mapping = sort {$a->[0] <=> $b->[0]} @{$updated_mscs_mappings->[$i]};
	$updated_mscs_mappings->[$i] = \@sorted_mscs_mapping;
    }

    return $updated_mscs_mappings, $aligned_stem_pairs;
}

sub _generate_independent_partition_candidate_vertex_pairs {
#    my ($mscs_mapping, $aligned_motif_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2, $base_seq1, $base_seq2) = @_;
    my ($mscs_mapping, $aligned_motif_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2,
	$stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2) = @_;

    my ($in_partition_stem_vertices1, $cross_partition_stem_vertices1) = _partition_non_nest_stem_vertices($mscs_mapping, $stem_graph1, 0);
    my ($in_partition_stem_vertices2, $cross_partition_stem_vertices2) = _partition_non_nest_stem_vertices($mscs_mapping, $stem_graph2, 1);

    my $in_partition_mscs_candidates = {};
    for (my $i = 0; $i < @{$in_partition_stem_vertices1}; $i++) {
	if (defined($in_partition_stem_vertices1->[$i][0]) && defined($in_partition_stem_vertices2->[$i][0])) {
	    $in_partition_mscs_candidates->{$i} = _generate_mscs_candidates($in_partition_stem_vertices1->[$i], $in_partition_stem_vertices2->[$i], $aligned_motif_pairs,
									    $aligned_stem_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2,
									    $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);
#	    $in_partition_mscs_candidates->{$i} = _generate_mscs_candidates($in_partition_stem_vertices1->[$i], $in_partition_stem_vertices2->[$i]);
	}
    }

    my $cross_partition_mscs_candidates = {};
    my $potential_mapped_cross_partition_ids = [];
    foreach (keys %{$cross_partition_stem_vertices1}) {
	if (exists($cross_partition_stem_vertices2->{$_})) {
	    $cross_partition_mscs_candidates->{$_} = _generate_mscs_candidates($cross_partition_stem_vertices1->{$_}, $cross_partition_stem_vertices2->{$_}, $aligned_motif_pairs,
									       $aligned_stem_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2,
									       $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);
#	    $cross_partition_mscs_candidates->{$_} = _generate_mscs_candidates($cross_partition_stem_vertices1->{$_}, $cross_partition_stem_vertices2->{$_});
	    push @{$potential_mapped_cross_partition_ids}, $_;
	}
    }

    my $linked_partition_map = {};
    my $max_linked_partition_id = -1;
    for (my $i = 0; $i < @{$potential_mapped_cross_partition_ids}; $i++) {
	my ($upstream_partition, $downstream_partition) = split(/-/, $potential_mapped_cross_partition_ids->[$i]);
	for (my $j = $i + 1; $j < @{$potential_mapped_cross_partition_ids}; $j++) {
	    my ($opponent_upstream_partition, $opponent_downstream_partition) = split(/-/, $potential_mapped_cross_partition_ids->[$j]);
	    if ($upstream_partition == $opponent_upstream_partition || $downstream_partition == $opponent_downstream_partition) {
		($linked_partition_map, $max_linked_partition_id) = _update_linked_partition_map($linked_partition_map, $max_linked_partition_id,
												 $potential_mapped_cross_partition_ids->[$i], $potential_mapped_cross_partition_ids->[$j]);
	    }
	    elsif (($upstream_partition < $opponent_upstream_partition && $downstream_partition > $opponent_upstream_partition) ||
		   ($opponent_upstream_partition < $upstream_partition && $opponent_downstream_partition > $upstream_partition)) {
		($linked_partition_map, $max_linked_partition_id) = _update_linked_partition_map($linked_partition_map, $max_linked_partition_id,
												 $potential_mapped_cross_partition_ids->[$i], $potential_mapped_cross_partition_ids->[$j]);
	    }
	    elsif ($downstream_partition == $opponent_upstream_partition) {
		if (_is_partitions_vertices_crossed($cross_partition_stem_vertices1->{$potential_mapped_cross_partition_ids->[$i]},
					  $cross_partition_stem_vertices1->{$potential_mapped_cross_partition_ids->[$j]}, $stem_graph1) ||
		    _is_partitions_vertices_crossed($cross_partition_stem_vertices2->{$potential_mapped_cross_partition_ids->[$i]},
					  $cross_partition_stem_vertices2->{$potential_mapped_cross_partition_ids->[$j]}, $stem_graph2)) {
		    ($linked_partition_map, $max_linked_partition_id) = _update_linked_partition_map($linked_partition_map, $max_linked_partition_id,
												     $potential_mapped_cross_partition_ids->[$i], $potential_mapped_cross_partition_ids->[$j]);
		}
	    }
	    elsif ($upstream_partition == $opponent_downstream_partition) {
		if (_is_partitions_vertices_crossed($cross_partition_stem_vertices1->{$potential_mapped_cross_partition_ids->[$j]},
					  $cross_partition_stem_vertices1->{$potential_mapped_cross_partition_ids->[$i]}, $stem_graph1) ||
		    _is_partitions_vertices_crossed($cross_partition_stem_vertices2->{$potential_mapped_cross_partition_ids->[$j]},
					  $cross_partition_stem_vertices2->{$potential_mapped_cross_partition_ids->[$i]}, $stem_graph2)) {
		    ($linked_partition_map, $max_linked_partition_id) = _update_linked_partition_map($linked_partition_map, $max_linked_partition_id,
												     $potential_mapped_cross_partition_ids->[$i], $potential_mapped_cross_partition_ids->[$j]);
		}
	    }
	}
    }

    for (my $i = 0; $i < @{$potential_mapped_cross_partition_ids}; $i++) {
	my ($upstream_partition, $downstream_partition) = split(/-/, $potential_mapped_cross_partition_ids->[$i]);

#	print "$upstream_partition, $downstream_partition\n";

	if (_is_partitions_vertices_crossed($in_partition_stem_vertices1->[$upstream_partition], $cross_partition_stem_vertices1->{$potential_mapped_cross_partition_ids->[$i]}, $stem_graph1) ||
	    _is_partitions_vertices_crossed($in_partition_stem_vertices2->[$upstream_partition], $cross_partition_stem_vertices2->{$potential_mapped_cross_partition_ids->[$i]}, $stem_graph2)) {
	    ($linked_partition_map, $max_linked_partition_id) = _update_linked_partition_map($linked_partition_map, $max_linked_partition_id,
											     $upstream_partition, $potential_mapped_cross_partition_ids->[$i]);
	}

	if (_is_partitions_vertices_crossed($cross_partition_stem_vertices1->{$potential_mapped_cross_partition_ids->[$i]}, $in_partition_stem_vertices1->[$downstream_partition], $stem_graph1) ||
	    _is_partitions_vertices_crossed($cross_partition_stem_vertices2->{$potential_mapped_cross_partition_ids->[$i]}, $in_partition_stem_vertices2->[$downstream_partition], $stem_graph2)) {
	    ($linked_partition_map, $max_linked_partition_id) = _update_linked_partition_map($linked_partition_map, $max_linked_partition_id,
											     $downstream_partition, $potential_mapped_cross_partition_ids->[$i]);
	}
    }

    my $independent_partition_mscs_candidates = [];
    my $all_linked_partition_mscs_candidates = {};
    while (my ($partition_id, $partition_mscs_candidates) = each %{$cross_partition_mscs_candidates}) {
	my $linked_partition_id = $linked_partition_map->{$partition_id};
	if (defined($linked_partition_id)) {
	    my $linked_partition_mscs_candidates = $all_linked_partition_mscs_candidates->{$linked_partition_id};
	    if (defined($linked_partition_mscs_candidates)) {
		push @{$linked_partition_mscs_candidates}, @{$partition_mscs_candidates};
	    }
	    else {
		$all_linked_partition_mscs_candidates->{$linked_partition_id} = $partition_mscs_candidates;
	    }
	}
	else {
	    push @{$independent_partition_mscs_candidates}, $partition_mscs_candidates;
	}
    }

    while (my ($partition_id, $partition_mscs_candidates) = each %{$in_partition_mscs_candidates}) {
	my $linked_partition_id = $linked_partition_map->{$partition_id};
	if (defined($linked_partition_id)) {
	    my $linked_partition_mscs_candidates = $all_linked_partition_mscs_candidates->{$linked_partition_id};
	    if (defined($linked_partition_mscs_candidates)) {
		push @{$linked_partition_mscs_candidates}, @{$partition_mscs_candidates};
	    }
	    else {
		$all_linked_partition_mscs_candidates->{$linked_partition_id} = $partition_mscs_candidates;
	    }
	}
	else {
	    push @{$independent_partition_mscs_candidates}, $partition_mscs_candidates;
	}
    }

    foreach (values %{$all_linked_partition_mscs_candidates}) {
	push @{$independent_partition_mscs_candidates}, $_;
    }

    return $independent_partition_mscs_candidates;
}

sub _partition_non_nest_stem_vertices {
    my ($mscs_mapping, $stem_graph, $graph_index) = @_;

    my ($non_nest_stem_partition_map, $non_nest_stem_partition_count, $last_partition_end) = _generate_non_nest_stem_partition_map($mscs_mapping, $stem_graph, $graph_index);

    my $in_partition_stem_vertices = [];
    my $cross_partition_stem_vertices = {};

    foreach (@{$stem_graph->get_non_nest_stem_vertices()}) {
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($_);
	my $vertex_upstream_partition = _get_non_nest_stem_partition($non_nest_stem_partition_map, $non_nest_stem_partition_count, $last_partition_end, $vertex_attrs->{upstream_start});
	my $vertex_downstream_partition = _get_non_nest_stem_partition($non_nest_stem_partition_map, $non_nest_stem_partition_count, $last_partition_end, $vertex_attrs->{downstream_start});
	if ($vertex_upstream_partition == $vertex_downstream_partition) {
	    my $same_partition_vertices = $in_partition_stem_vertices->[$vertex_upstream_partition];
	    if (defined($same_partition_vertices)) {
		push @{$same_partition_vertices}, $_;
	    }
	    else {
		$in_partition_stem_vertices->[$vertex_upstream_partition] = [$_];
	    }
	}
	else {
	    my $same_cross_vertices = $cross_partition_stem_vertices->{$vertex_upstream_partition . '-' . $vertex_downstream_partition};
	    if (defined($same_cross_vertices)) {
		push @{$same_cross_vertices}, $_;
	    }
	    else {
		$cross_partition_stem_vertices->{$vertex_upstream_partition . '-' . $vertex_downstream_partition} = [$_];
	    }
	}
    }

    return $in_partition_stem_vertices, $cross_partition_stem_vertices;
}

sub _generate_mscs_candidates {
    my ($stem_vertices1, $stem_vertices2, $aligned_motif_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2,
	$stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2) = @_;

    if (@{$stem_vertices1} <= DefaultAlignParams->PARTITION_STEM_COUNT_THRESHOLD &&
	@{$stem_vertices2} <= DefaultAlignParams->PARTITION_STEM_COUNT_THRESHOLD) {
	my $stem_vertex_pairs = [];

	foreach my $stem_vertex1 (@{$stem_vertices1}) {
	    foreach my $stem_vertex2 (@{$stem_vertices2}) {
		push @{$stem_vertex_pairs}, [$stem_vertex1, $stem_vertex2];
	    }
	}

	return $stem_vertex_pairs;
    }

    my ($crossing_stem_vertices1, $crossing_stem_vertices2) = ([], []);
    my ($non_crossing_stem_vertices1, $non_crossing_stem_vertices2) = ([], []);
    foreach (sort {$a <=> $b} @{$stem_vertices1}) {
	if ($stem_graph1->get_stem_type_at($_) eq 'T') {
	    push @{$crossing_stem_vertices1}, $_;
	}
	else {
	    push @{$non_crossing_stem_vertices1}, $_;
	}
    }

    foreach (sort {$a <=> $b} @{$stem_vertices2}) {
	if ($stem_graph2->get_stem_type_at($_) eq 'T') {
	    push @{$crossing_stem_vertices2}, $_;
	}
	else {
	    push @{$non_crossing_stem_vertices2}, $_;
	}
    }

    my $non_crossing_candidate_vertex_pairs = {};
    my $probe_dp_table = _create_probe_dp_table($non_crossing_stem_vertices1, $non_crossing_stem_vertices2, $aligned_motif_pairs, $stem_graph1, $stem_graph2,
						$motif_removal_costs1, $motif_removal_costs2, $base_seq1, $base_seq2);
    my $dp_motif_aligned_vertex_pairs = _backtrack_probe_dp_table($probe_dp_table, $non_crossing_stem_vertices1, $non_crossing_stem_vertices2);
    foreach (@{$dp_motif_aligned_vertex_pairs}) {
	$non_crossing_candidate_vertex_pairs->{$_->[0] . '-' . $_->[1]} = $_;
    }

    $probe_dp_table = _create_probe_stem_align_dp_table($non_crossing_stem_vertices1, $non_crossing_stem_vertices2, $aligned_stem_pairs, $stem_graph1, $stem_graph2,
							$stem_removal_costs1, $stem_removal_costs2);
    my $dp_stem_aligned_vertex_pairs = _backtrack_probe_dp_table($probe_dp_table, $non_crossing_stem_vertices1, $non_crossing_stem_vertices2);
    foreach (@{$dp_stem_aligned_vertex_pairs}) {
	$non_crossing_candidate_vertex_pairs->{$_->[0] . '-' . $_->[1]} = $_;
    }

    my @candidate_vertex_pairs = values %{$non_crossing_candidate_vertex_pairs};

    my $motif_aligned_vertex_pair_count = @{$dp_motif_aligned_vertex_pairs};
    my $stem_aligned_vertex_pair_count = @{$dp_stem_aligned_vertex_pairs};
    my $candidate_vertex_pair_count = @candidate_vertex_pairs;
    if (($candidate_vertex_pair_count - $motif_aligned_vertex_pair_count > DefaultAlignParams->PARTITION_STEM_COUNT_THRESHOLD) &&
	($candidate_vertex_pair_count - $stem_aligned_vertex_pair_count > DefaultAlignParams->PARTITION_STEM_COUNT_THRESHOLD)) {
	if ($motif_aligned_vertex_pair_count >= $stem_aligned_vertex_pair_count) {
	    @candidate_vertex_pairs = @{$dp_motif_aligned_vertex_pairs};
	}
	else {
	    @candidate_vertex_pairs = @{$dp_stem_aligned_vertex_pairs};
	}
    }

=comment
    foreach my $dp_row (@{$probe_dp_table}) {
	foreach (@{$dp_row}) {
	    print $_ . "\t";
	}
	print "\n";
    }
=cut

    foreach my $crossing_stem_vertex1 (@{$crossing_stem_vertices1}) {
	foreach my $crossing_stem_vertex2 (@{$crossing_stem_vertices2}) {
	    push @candidate_vertex_pairs, [$crossing_stem_vertex1, $crossing_stem_vertex2];
	}
    }

    my $mixed_type_candidate_vertex_pairs = _get_mixed_type_candidate_vertex_pairs($crossing_stem_vertices1, $non_crossing_stem_vertices2, $aligned_motif_pairs,
										   $aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2, 0);
    push @candidate_vertex_pairs, @{$mixed_type_candidate_vertex_pairs};

    $mixed_type_candidate_vertex_pairs = _get_mixed_type_candidate_vertex_pairs($crossing_stem_vertices2, $non_crossing_stem_vertices1, $aligned_motif_pairs,
										$aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2, 1);
    push @candidate_vertex_pairs, @{$mixed_type_candidate_vertex_pairs};

#    return $candidate_vertex_pairs;
    return \@candidate_vertex_pairs;
}

sub _create_probe_stem_align_dp_table {
    my ($non_crossing_stem_vertices1, $non_crossing_stem_vertices2, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2) = @_;

    my $first_row = [0];
    foreach (@{$non_crossing_stem_vertices1}) {
	push @{$first_row}, ($first_row->[-1] + $stem_removal_costs1->[$_]);
    }

    my $dp_table = [$first_row];

    foreach (@{$non_crossing_stem_vertices2}) {
	my $init_value = $dp_table->[-1][0] + $stem_removal_costs2->[$_];
	push @{$dp_table}, [$init_value];
    }

    my $dp_cell_vertex_pairs = [];
    for (my $i = 0; $i < @{$non_crossing_stem_vertices2}; $i++) {
	for (my $j = 0; $j < @{$non_crossing_stem_vertices1}; $j++) {
	    push @{$dp_cell_vertex_pairs}, [$non_crossing_stem_vertices1->[$j], $non_crossing_stem_vertices2->[$i]];
	}
    }

    $aligned_stem_pairs = GraphMatchTools->align_stem_pairs($dp_cell_vertex_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2);

    for (my $i = 0; $i < @{$non_crossing_stem_vertices2}; $i++) {
	for (my $j = 0; $j < @{$non_crossing_stem_vertices1}; $j++) {
	    my $stem_align_cost = $aligned_stem_pairs->{$non_crossing_stem_vertices1->[$j] . '-' . $non_crossing_stem_vertices2->[$i]}->get_cost();

	    my $dp_cell_cost = $dp_table->[$i][$j] + $stem_align_cost;

	    my $dp_cell_candidate_cost1 = $dp_table->[$i + 1][$j] + $stem_removal_costs1->[$non_crossing_stem_vertices1->[$j]];
	    if ($dp_cell_candidate_cost1 < $dp_cell_cost) {
		$dp_cell_cost = $dp_cell_candidate_cost1;
	    }

	    my $dp_cell_candidate_cost2 = $dp_table->[$i][$j + 1] + $stem_removal_costs2->[$non_crossing_stem_vertices2->[$i]];
	    if ($dp_cell_candidate_cost2 < $dp_cell_cost) {
		$dp_cell_cost = $dp_cell_candidate_cost2;
	    }

	    $dp_table->[$i + 1][$j + 1] = $dp_cell_cost;
	}
    }

    return $dp_table;
}

sub _create_probe_dp_table {
    my ($non_crossing_stem_vertices1, $non_crossing_stem_vertices2, $aligned_motif_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2,
	$base_seq1, $base_seq2) = @_;

    my $first_row = [0];
    foreach (@{$non_crossing_stem_vertices1}) {
	push @{$first_row}, ($first_row->[-1] + $motif_removal_costs1->[$_]);
    }

    my $dp_table = [$first_row];

    foreach (@{$non_crossing_stem_vertices2}) {
	my $init_value = $dp_table->[-1][0] + $motif_removal_costs2->[$_];
	push @{$dp_table}, [$init_value];
    }

    for (my $i = 0; $i < @{$non_crossing_stem_vertices2}; $i++) {
	for (my $j = 0; $j < @{$non_crossing_stem_vertices1}; $j++) {
	    $aligned_motif_pairs = GraphMatchTools->align_motif_pair($non_crossing_stem_vertices1->[$j], $non_crossing_stem_vertices2->[$i], $aligned_motif_pairs,
								     $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
	    my $motif_align_cost = $aligned_motif_pairs->{$non_crossing_stem_vertices1->[$j] . '-' . $non_crossing_stem_vertices2->[$i]}->get_cost();

	    my $dp_cell_cost = $dp_table->[$i][$j] + $motif_align_cost;

	    my $dp_cell_candidate_cost1 = $dp_table->[$i + 1][$j] + $motif_removal_costs1->[$non_crossing_stem_vertices1->[$j]];
	    if ($dp_cell_candidate_cost1 < $dp_cell_cost) {
		$dp_cell_cost = $dp_cell_candidate_cost1;
	    }

	    my $dp_cell_candidate_cost2 = $dp_table->[$i][$j + 1] + $motif_removal_costs2->[$non_crossing_stem_vertices2->[$i]];
	    if ($dp_cell_candidate_cost2 < $dp_cell_cost) {
		$dp_cell_cost = $dp_cell_candidate_cost2;
	    }

	    $dp_table->[$i + 1][$j + 1] = $dp_cell_cost;
	}
    }

    return $dp_table;
}

sub _backtrack_probe_dp_table {
    my ($dp_table, $non_crossing_stem_vertices1, $non_crossing_stem_vertices2) = @_;

    my $candidate_vertex_pairs = [];
    my ($last_tracked_cell_row_index, $last_tracked_cell_col_index) = (scalar(@{$non_crossing_stem_vertices2}), scalar(@{$non_crossing_stem_vertices1}));

    while ($last_tracked_cell_row_index > 0 || $last_tracked_cell_col_index > 0) {
	my ($tracked_cell_row_index, $tracked_cell_col_index);
	my $lowest_cost;
	if ($last_tracked_cell_row_index > 0 && $last_tracked_cell_col_index > 0) {
	    $tracked_cell_row_index = $last_tracked_cell_row_index - 1;
	    $tracked_cell_col_index = $last_tracked_cell_col_index - 1;
	    $lowest_cost = $dp_table->[$tracked_cell_row_index][$tracked_cell_col_index];
	}

	if ($last_tracked_cell_row_index > 0) {
	    if (!defined($lowest_cost) || $dp_table->[$last_tracked_cell_row_index - 1][$last_tracked_cell_col_index] < $lowest_cost) {
		$tracked_cell_row_index = $last_tracked_cell_row_index - 1;
		$tracked_cell_col_index = $last_tracked_cell_col_index;
		$lowest_cost = $dp_table->[$tracked_cell_row_index][$tracked_cell_col_index];
	    }
	}

	if ($last_tracked_cell_col_index > 0) {
	    if (!defined($lowest_cost) || $dp_table->[$last_tracked_cell_row_index][$last_tracked_cell_col_index - 1] < $lowest_cost) {
		$tracked_cell_row_index = $last_tracked_cell_row_index;
		$tracked_cell_col_index = $last_tracked_cell_col_index - 1;
		$lowest_cost = $dp_table->[$tracked_cell_row_index][$tracked_cell_col_index];
	    }
	}

	if ($tracked_cell_row_index < $last_tracked_cell_row_index && $tracked_cell_col_index < $last_tracked_cell_col_index) {
	    my $candidate_vertex_pair = [$non_crossing_stem_vertices1->[$tracked_cell_col_index], $non_crossing_stem_vertices2->[$tracked_cell_row_index]];
	    unshift @{$candidate_vertex_pairs}, $candidate_vertex_pair;
	}

	$last_tracked_cell_row_index = $tracked_cell_row_index;
	$last_tracked_cell_col_index = $tracked_cell_col_index;
    }

    return $candidate_vertex_pairs;
}

sub _get_mixed_type_candidate_vertex_pairs {
    my ($crossing_stem_vertices, $non_crossing_stem_vertices, $aligned_motif_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $is_reverse) = @_;

    my $candidate_vertex_pairs = [];
    my ($stem_vertex1, $stem_vertex2);

    foreach my $crossing_stem_vertex (@{$crossing_stem_vertices}) {
	if ($is_reverse) {
	    $stem_vertex2 = $crossing_stem_vertex;
	}
	else {
	    $stem_vertex1 = $crossing_stem_vertex;
	}

	my ($aligned_motif_cost_rank, $aligned_stem_cost_rank) = ({}, {});
	my $tentative_candidate_vertex_pairs = [];

	foreach my $non_crossing_stem_vertex (@{$non_crossing_stem_vertices}) {
	    if ($is_reverse) {
		$stem_vertex1 = $non_crossing_stem_vertex;
	    }
	    else {
		$stem_vertex2 = $non_crossing_stem_vertex;
	    }

	    $aligned_motif_pairs = GraphMatchTools->align_motif_pair($stem_vertex1, $stem_vertex2, $aligned_motif_pairs, $stem_graph1, $stem_graph2,
								     $base_seq1, $base_seq2);
	    my $motif_align_cost = $aligned_motif_pairs->{$stem_vertex1 . '-' . $stem_vertex2}->get_cost();

	    my $same_rank_motif_pairs = $aligned_motif_cost_rank->{$motif_align_cost};
	    if (!defined($same_rank_motif_pairs)) {
		$same_rank_motif_pairs = [];
		$aligned_motif_cost_rank->{$motif_align_cost} = $same_rank_motif_pairs;
	    }

	    push @{$same_rank_motif_pairs}, [$stem_vertex1, $stem_vertex2];
	    push @{$tentative_candidate_vertex_pairs}, [$stem_vertex1, $stem_vertex2];
	}

	$aligned_stem_pairs = GraphMatchTools->align_stem_pairs($tentative_candidate_vertex_pairs, $aligned_stem_pairs, $stem_graph1, $stem_graph2);
	foreach (@{$tentative_candidate_vertex_pairs}) {
	    my $stem_align_cost = $aligned_stem_pairs->{$_->[0] . '-' . $_->[1]}->get_cost();

	    my $same_rank_stem_pairs = $aligned_stem_cost_rank->{$stem_align_cost};
	    if (!defined($same_rank_stem_pairs)) {
		$same_rank_stem_pairs = [];
		$aligned_stem_cost_rank->{$stem_align_cost} = $same_rank_stem_pairs;
	    }

	    push @{$same_rank_stem_pairs}, [$_->[0], $_->[1]];
	}

	my $selected_vertex_pairs = {};
	my $selected_vertex_pair_count = 0;

	foreach my $motif_align_cost (sort {$a <=> $b} keys %{$aligned_motif_cost_rank}) {
	    my $same_rank_motif_pairs = $aligned_motif_cost_rank->{$motif_align_cost};
	    foreach (@{$same_rank_motif_pairs}) {
		$selected_vertex_pairs->{$_->[0] . '-' . $_->[1]} = [$_->[0], $_->[1]];
	    }

	    $selected_vertex_pair_count += @{$same_rank_motif_pairs};
	    if ($selected_vertex_pair_count >= DefaultAlignParams->MAX_MIXED_STEM_TYPE_CANDIDATES) {
		last;
	    }
	}

	$selected_vertex_pair_count = 0;
	foreach my $stem_align_cost (sort {$a <=> $b} keys %{$aligned_stem_cost_rank}) {
	    my $same_rank_stem_pairs = $aligned_stem_cost_rank->{$stem_align_cost};
	    foreach (@{$same_rank_stem_pairs}) {
		$selected_vertex_pairs->{$_->[0] . '-' . $_->[1]} = [$_->[0], $_->[1]];
	    }

	    $selected_vertex_pair_count += @{$same_rank_stem_pairs};
	    if ($selected_vertex_pair_count >= DefaultAlignParams->MAX_MIXED_STEM_TYPE_CANDIDATES) {
		last;
	    }
	}

	push @{$candidate_vertex_pairs}, values %{$selected_vertex_pairs};
    }

    return $candidate_vertex_pairs;
}

=comment
sub _generate_mscs_candidates {
    my ($stem_vertices1, $stem_vertices2) = @_;

    my $stem_vertex_pairs = [];

    foreach my $stem_vertex1 (@{$stem_vertices1}) {
	foreach my $stem_vertex2 (@{$stem_vertices2}) {
	    push @{$stem_vertex_pairs}, [$stem_vertex1, $stem_vertex2];
	}
    }

    return $stem_vertex_pairs;
}
=cut

sub _update_linked_partition_map {
    my ($linked_partition_map, $max_linked_partition_id, $first_partition_id, $second_partition_id) = @_;

    my $first_linked_partition_id = $linked_partition_map->{$first_partition_id};
    my $second_linked_partition_id = $linked_partition_map->{$second_partition_id};
    if (defined($first_linked_partition_id)) {
	if (defined($second_linked_partition_id)) {
	    if ($first_linked_partition_id != $second_linked_partition_id) {
		while (my ($partition_id, $linked_partition_id) = each %{$linked_partition_map}) {
		    if ($linked_partition_id == $second_linked_partition_id)  {
			$linked_partition_map->{$partition_id} = $first_linked_partition_id;
		    }
		}
	    }
	}
	else {
	    $linked_partition_map->{$second_partition_id} = $first_linked_partition_id;
	}
    }
    else {
	if (defined($second_linked_partition_id)) {
	    $linked_partition_map->{$first_partition_id} = $second_linked_partition_id;
	}
	else {
	    $linked_partition_map->{$first_partition_id} = ++$max_linked_partition_id;
	    $linked_partition_map->{$second_partition_id} = $max_linked_partition_id;
	}
    }

    return $linked_partition_map, $max_linked_partition_id;
}

sub _is_partitions_vertices_crossed {
    my ($precede_partition_stem_vertices, $subseq_partition_stem_vertices, $stem_graph) = @_;

    if (!defined($precede_partition_stem_vertices->[0]) || !defined($subseq_partition_stem_vertices->[0])) {
	return 0;
    }

    return ($stem_graph->get_edge_label($precede_partition_stem_vertices->[-1], $subseq_partition_stem_vertices->[0]) eq 'K');
}

sub _generate_non_nest_stem_partition_map {
    my ($nest_mscs_mapping, $stem_graph, $graph_index) = @_;

    my ($non_nest_stem_partition_map, $nest_stem_pos) = ({}, {});

    foreach (@{$nest_mscs_mapping}) {
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($_->[$graph_index]);
	$nest_stem_pos->{$vertex_attrs->{upstream_start}} = 1;
	$nest_stem_pos->{$vertex_attrs->{downstream_start}} = 1;
    }

    my @sorted_nest_stem_pos = sort {$a <=> $b} keys %{$nest_stem_pos};
    
    for (my $i = 0; $i < @sorted_nest_stem_pos; $i++) {
	if ($i == 0) {
	    for (my $j = 0; $j < $sorted_nest_stem_pos[$i]; $j++) {
		$non_nest_stem_partition_map->{$j} = $i;
	    }
	}
	else {
	    for (my $j = $sorted_nest_stem_pos[$i - 1] + 1; $j < $sorted_nest_stem_pos[$i]; $j++) {
		$non_nest_stem_partition_map->{$j} = $i;
	    }
	}
    }

    return $non_nest_stem_partition_map, scalar(@sorted_nest_stem_pos), $sorted_nest_stem_pos[-1];
}

sub _get_non_nest_stem_partition {
    my ($non_nest_stem_partition_map, $non_nest_stem_partition_count, $last_partition_end, $non_nest_stem_pos) = @_;

    if ($non_nest_stem_pos > $last_partition_end) {
	return $non_nest_stem_partition_count;
    }

    return $non_nest_stem_partition_map->{$non_nest_stem_pos};
}

sub _add_partition_mscs_mappings {
    my ($non_redundant_part_mscs_mappings, $partition_mscs_mappings) = @_;

    foreach my $partition_mscs_mapping (@{$partition_mscs_mappings}) {
	my $part_mscs_mapping_id = '';
	foreach (@{$partition_mscs_mapping}) {
	    $part_mscs_mapping_id = $part_mscs_mapping_id . $_->[0] . '-' . $_->[1] . '_';
	}

	$non_redundant_part_mscs_mappings->{$part_mscs_mapping_id} = $partition_mscs_mapping;
    }

    return $non_redundant_part_mscs_mappings;
}

=comment
sub _get_mapping_crossed_stem_vertices {
    my ($stem_graph1, $stem_graph2, $mscs_mapping, $is_nest_stem_vertex) = @_;

    my $selected_cross_stem_vertices1 = _get_crossing_stem_vertices($stem_graph1, $is_nest_stem_vertex);
    my $selected_cross_stem_vertices2 = _get_crossing_stem_vertices($stem_graph2, $is_nest_stem_vertex);

    my ($mapping_crossed_stem_vertices1, $mapping_crossed_stem_vertices2) = ([], []);
    foreach my $stem_vertex (@{$selected_cross_stem_vertices1}) {
	foreach (@{$mscs_mapping}) {
	    if (($_->[0] > $stem_vertex && $stem_graph1->get_edge_label($stem_vertex, $_->[0]) eq 'K') ||
		($_->[0] < $stem_vertex && $stem_graph1->get_edge_label($_->[0], $stem_vertex) eq 'K')) {
		push @{$mapping_crossed_stem_vertices1}, $stem_vertex;
		last;
	    }
	}
    }

    foreach my $stem_vertex (@{$selected_cross_stem_vertices2}) {
	foreach (@{$mscs_mapping}) {
	    if (($_->[1] > $stem_vertex && $stem_graph2->get_edge_label($stem_vertex, $_->[1]) eq 'K') ||
		($_->[1] < $stem_vertex && $stem_graph2->get_edge_label($_->[1], $stem_vertex) eq 'K')) {
		push @{$mapping_crossed_stem_vertices2}, $stem_vertex;
		last;
	    }
	}
    }

    return $mapping_crossed_stem_vertices1, $mapping_crossed_stem_vertices2;
}

sub _get_crossing_stem_vertices {
    my ($stem_graph, $is_nest_stem_vertex) = @_;

    my $stem_vertices;
    if ($is_nest_stem_vertex) {
	$stem_vertices = $stem_graph->get_nest_stem_vertices();
    }
    else {
	$stem_vertices = $stem_graph->get_non_nest_stem_vertices();
    }

    my $crossing_stem_vertices = [];

    foreach (@{$stem_vertices}) {
	if ($stem_graph->get_stem_type_at($_) eq 'T') {
	    push @{$crossing_stem_vertices}, $_;
	}
    }

    return $crossing_stem_vertices;
}
=cut

sub _precalculate_stem_removal_costs {
    my $stem_graph = shift;

    my $stem_removal_costs = [];

    for (my $i = 0; $i < $stem_graph->get_vertex_count(); $i++) {
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($i);
	$stem_removal_costs->[$i] = $vertex_attrs->{base_pair_count} * DefaultAlignParams->PAIR_REMOVAL_COST + ($vertex_attrs->{upstream_length} +
	   $vertex_attrs->{downstream_length} - $vertex_attrs->{base_pair_count} * 2) * DefaultAlignParams->BASE_REMOVAL_COST;
    }

    return $stem_removal_costs;
}

sub _merge_mscs_mappings {
    my ($nest_crossed_non_nest_mscs_mapping, $all_part_non_nest_mscs_mappings) = @_;

    my $mscs_mappings = [$nest_crossed_non_nest_mscs_mapping];

    foreach my $part_non_nest_mscs_mappings (@{$all_part_non_nest_mscs_mappings}) {
	my $expanded_mscs_mappings = [];
	foreach my $part_non_nest_mscs_mapping (@{$part_non_nest_mscs_mappings}) {
	    foreach (@{$mscs_mappings}) {
		my @merged_mscs_mapping = (@{$_}, @{$part_non_nest_mscs_mapping});
		push @{$expanded_mscs_mappings}, \@merged_mscs_mapping;
	    }
	}

	$mscs_mappings = $expanded_mscs_mappings;
    }

    return $mscs_mappings;
}

sub _align_all_stems {
    my ($stem_graph1, $stem_graph2, $aligned_stem_pairs, $base_seq1, $base_seq2) = @_;

    my $nest_stem_vertices = $stem_graph1->get_nest_stem_vertices();
    my $non_nest_stem_vertices = $stem_graph1->get_non_nest_stem_vertices();
    my @stem_vertices1 = sort {$a <=> $b} (@{$nest_stem_vertices}, @{$non_nest_stem_vertices});

    $nest_stem_vertices = $stem_graph2->get_nest_stem_vertices();
    $non_nest_stem_vertices = $stem_graph2->get_non_nest_stem_vertices();
    my @stem_vertices2 = sort {$a <=> $b} (@{$nest_stem_vertices}, @{$non_nest_stem_vertices});

    my $mscs_candidate_vertex_pairs = GraphMatchTools->generate_unmapped_stem_vertex_pairs(\@stem_vertices1, \@stem_vertices2, [],
											   $stem_graph1, $stem_graph2, 1);
    my $stem_removal_costs1 = _precalculate_stem_removal_costs($stem_graph1);
    my $stem_removal_costs2 = _precalculate_stem_removal_costs($stem_graph2);

    (my ($mscs_mapping, $aligned_stem_pairs), undef) =
	GraphMatchTools->expand_mscs_mapping($mscs_candidate_vertex_pairs, [], $aligned_stem_pairs, $stem_graph1, $stem_graph2,
					     $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2);

    return $mscs_mapping, $aligned_stem_pairs;
}

1;
