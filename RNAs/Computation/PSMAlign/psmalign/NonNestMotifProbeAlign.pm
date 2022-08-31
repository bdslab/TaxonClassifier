package NonNestMotifProbeAlign;

use strict;

sub match_non_nest_motifs {
    my (undef, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2, $align_params) = @_;

    my $non_nest_stem_vertices1 = $stem_graph1->get_non_nest_stem_vertices();
    my $non_nest_stem_vertices2 = $stem_graph2->get_non_nest_stem_vertices();

    my $non_nest_stem_count1 = @{$non_nest_stem_vertices1};
    my $non_nest_stem_count2 = @{$non_nest_stem_vertices2};
    my $max_non_nest_motif_match_count = ($non_nest_stem_count1, $non_nest_stem_count2)[$non_nest_stem_count1 > $non_nest_stem_count2];
    if ($max_non_nest_motif_match_count < $align_params->get_non_nest_motif_match_threshold()) {
	return [];
    }

    my $motif_removal_costs1 = _get_motif_removal_costs($non_nest_stem_vertices1, $stem_graph1, $stem_removal_costs1, $base_seq1);
    my $motif_removal_costs2 = _get_motif_removal_costs($non_nest_stem_vertices2, $stem_graph2, $stem_removal_costs2, $base_seq2);

    my $cost_ratio_selection_factor = $align_params->get_cost_ratio_selection_factor();

    my ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs) =
	_align_similar_len_motifs($non_nest_stem_vertices1, $non_nest_stem_vertices2, $motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);

    foreach my $non_nest_stem_vertex1 (@{$non_nest_stem_vertices1}) {
	my $vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($non_nest_stem_vertex1);
	my $base_pair_count1 = $vertex_attrs1->{base_pair_count};
	my $unpaired_base_count1 = $vertex_attrs1->{upstream_length} + $vertex_attrs1->{downstream_length} - 2 * $base_pair_count1;

	foreach my $non_nest_stem_vertex2 (@{$non_nest_stem_vertices2}) {
	    my $candidate_id = $non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2;
	    if (exists($motif_cost_ratios->{$candidate_id})) {
		next;
	    }

	    my $min_motif_removal_cost = $motif_removal_costs1->[$non_nest_stem_vertex1];
	    if ($motif_removal_costs2->[$non_nest_stem_vertex2] < $min_motif_removal_cost) {
		$min_motif_removal_cost = $motif_removal_costs2->[$non_nest_stem_vertex2];
	    }

	    my $vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($non_nest_stem_vertex2);
	    my $base_pair_count2 = $vertex_attrs2->{base_pair_count};
	    my $unpaired_base_count2 = $vertex_attrs2->{upstream_length} + $vertex_attrs2->{downstream_length} - 2 * $base_pair_count2;

	    my $best_motif_align_cost = _calculate_best_motif_align_cost($base_pair_count1, $base_pair_count2, $unpaired_base_count1, $unpaired_base_count2);
	    my $best_motif_cost_ratio = $best_motif_align_cost / $min_motif_removal_cost;
	    if ($best_motif_cost_ratio > $cost_ratio_selection_factor * $min_vertex_cost_ratios1->{$non_nest_stem_vertex1} &&
		$best_motif_cost_ratio > $cost_ratio_selection_factor * $min_vertex_cost_ratios2->{$non_nest_stem_vertex2}) {
		next;
	    }

	    $aligned_motif_pairs = GraphMatchTools->align_motif_pair($non_nest_stem_vertex1, $non_nest_stem_vertex2, $aligned_motif_pairs, $stem_graph1,
								     $stem_graph2, $base_seq1, $base_seq2);
	    my $motif_mapping_cost_ratio = $aligned_motif_pairs->{$candidate_id}->get_cost() / $min_motif_removal_cost;
	    $motif_cost_ratios->{$candidate_id} = $motif_mapping_cost_ratio;

	    if ($motif_mapping_cost_ratio < $min_vertex_cost_ratios1->{$non_nest_stem_vertex1}) {
		$min_vertex_cost_ratios1->{$non_nest_stem_vertex1} = $motif_mapping_cost_ratio;
	    }

	    if ($motif_mapping_cost_ratio < $min_vertex_cost_ratios2->{$non_nest_stem_vertex2}) {
		$min_vertex_cost_ratios2->{$non_nest_stem_vertex2} = $motif_mapping_cost_ratio;
	    }

#	    my $same_rank_vertex_pairs = $vertex_pair_ranks->{$motif_mapping_cost_ratio};
#	    if (!defined($same_rank_vertex_pairs)) {
#		$same_rank_vertex_pairs = [];
#		$vertex_pair_ranks->{$motif_mapping_cost_ratio} = $same_rank_vertex_pairs;
#	    }

#	    push @{$same_rank_vertex_pairs}, [$non_nest_stem_vertex1, $non_nest_stem_vertex2];
	}
    }

    my $vertex_pair_ranks = {};

    foreach my $non_nest_stem_vertex1 (@{$non_nest_stem_vertices1}) {
	foreach my $non_nest_stem_vertex2 (@{$non_nest_stem_vertices2}) {
	    my $candidate_id = $non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2;
	    if (!exists($motif_cost_ratios->{$candidate_id})) {
		next;
	    }

	    my $motif_cost_ratio = $motif_cost_ratios->{$candidate_id};
	    if ($motif_cost_ratio <= $cost_ratio_selection_factor * $min_vertex_cost_ratios1->{$non_nest_stem_vertex1} ||
		$motif_cost_ratio <= $cost_ratio_selection_factor * $min_vertex_cost_ratios2->{$non_nest_stem_vertex2}) {
		my $same_rank_vertex_pairs = $vertex_pair_ranks->{$motif_cost_ratio};
		if (!defined($same_rank_vertex_pairs)) {
		    $same_rank_vertex_pairs = [];
		    $vertex_pair_ranks->{$motif_cost_ratio} = $same_rank_vertex_pairs;
		}

		push @{$same_rank_vertex_pairs}, [$non_nest_stem_vertex1, $non_nest_stem_vertex2];
	    }
	}
    }

    my @sorted_cost_ratios = sort {$a <=> $b} keys %{$vertex_pair_ranks};

    my $non_nest_motif_mscs_mappings = [[]];
    my $max_non_nest_motif_mscs_mapping = [];

    my $select_range = [-1, $align_params->get_init_cost_ratio_select_range_end()];
#    while ($select_range->[1] <= $align_params->get_max_select_range_cost_ratio()) {
    while (1) {
	my $candidate_vertex_pairs = [];
	my $candidate_vertex_pair_cost_ratios = {};
	foreach my $mapping_cost_ratio (@sorted_cost_ratios) {
	    if ($mapping_cost_ratio <= $select_range->[0]) {
		next;
	    }

	    if ($mapping_cost_ratio > $select_range->[1]) {
		last;
	    }

	    push @{$candidate_vertex_pairs}, @{$vertex_pair_ranks->{$mapping_cost_ratio}};

	    foreach (@{$vertex_pair_ranks->{$mapping_cost_ratio}}) {
		$candidate_vertex_pair_cost_ratios->{$_->[0] . '-' . $_->[1]} = $mapping_cost_ratio;
	    }
=comment
	    print $mapping_cost_ratio . ': ';
	    foreach (@{$vertex_pair_ranks->{$mapping_cost_ratio}}) {
		print '[' . $_->[0] . ', ' . $_->[1] . '] ';
	    }
	    print "\n";
=cut
	}

=comment
	print "Round new selected pairs\n";
	foreach (sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$candidate_vertex_pairs}) {
	    print '[' . $_->[0] . ', ' . $_->[1] . ' (' . ($_->[0] - $_->[1]) . ')] ';
	}
	print "\nEnd\n";
=cut
	my $updated_non_nest_motif_mscs_mappings = [];

	foreach (@{$non_nest_motif_mscs_mappings}) {
=comment
	    print "After removal\n";
	    foreach (sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$candidate_vertex_pairs}) {
		print '[' . $_->[0] . ', ' . $_->[1] . ' (' . ($_->[0] - $_->[1]) . ')] ';
	    }
	    print "\nEnd\n";
=cut

#	    my ($expanded_mscs_mappings, undef) = GraphMatchTools->expand_non_nest_stem_vertex_align_trend($candidate_vertex_pairs, $_, $selected_vertex_pair_cost_ratios, $stem_graph1, $stem_graph2);
	    my $expanded_motif_mscs_mappings = GraphMatchTools->expand_non_nest_motif_mscs_mapping($candidate_vertex_pairs, $_, $aligned_motif_pairs, $stem_graph1,
												   $stem_graph2, $motif_removal_costs1, $motif_removal_costs2);
	    push @{$updated_non_nest_motif_mscs_mappings}, @{$expanded_motif_mscs_mappings};
	}

	$non_nest_motif_mscs_mappings = [];
	foreach (@{$updated_non_nest_motif_mscs_mappings}) {
	    if (@{$_} >= $align_params->get_non_nest_motif_match_threshold()) {
		push @{$non_nest_motif_mscs_mappings}, $_;
	    }
	}

	if (!defined($non_nest_motif_mscs_mappings->[0])) {
	    $non_nest_motif_mscs_mappings = [[]];
	}

#	$non_nest_motif_mscs_mappings = $updated_non_nest_mscs_mappings;
#	print $select_range->[0] . ' - ' . $select_range->[1] . "\n";
=comment
	print "======\n";
	print $select_range->[0] . ' - ' . $select_range->[1] . "\n";
	foreach my $best_mapping (@{$non_nest_motif_mscs_mappings}) {
	    foreach (@{$best_mapping}) {
		print '[' . $_->[0] . ', ' . $_->[1] . ' (' . ($_->[0] - $_->[1]) . ')] ';
	    }
	    print "\n";
	}
	print "======\n";
=cut

	foreach my $non_nest_motif_mscs_mapping (@{$non_nest_motif_mscs_mappings}) {
	    if (scalar(@{$non_nest_motif_mscs_mapping}) > scalar(@{$max_non_nest_motif_mscs_mapping})) {
		$max_non_nest_motif_mscs_mapping = $non_nest_motif_mscs_mapping;
	    }
	}

	if ($select_range->[1] >= $sorted_cost_ratios[-1]) {
	    last;
	}

	if (defined($non_nest_motif_mscs_mappings->[0][0])) {
	    $select_range->[0] = $select_range->[1];
	}

	$select_range->[1] += $align_params->get_select_range_size();
=comment
	print "======\n";
	print $select_range->[0] . ' - ' . $select_range->[1] . "\n";
	foreach (@{$max_non_nest_motif_mscs_mapping}) {
	    print '[' . $_->[0] . ', ' . $_->[1] . ' (' . ($_->[0] - $_->[1]) . ')] ';
	}
	print "\n";
	print "======\n";
=cut
    }

    return $max_non_nest_motif_mscs_mapping, $aligned_motif_pairs, $motif_removal_costs1, $motif_removal_costs2;
}

sub _get_motif_removal_costs {
    my ($non_nest_stem_vertices, $stem_graph, $stem_removal_costs, $base_seq) = @_;

    my $motif_removal_costs = [];

    foreach (@{$non_nest_stem_vertices}) {
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($_);
	my $stem_base_pairs = $vertex_attrs->{base_pairs};
	my $motif_loop_seq = substr($base_seq, $stem_base_pairs->[-1][0] + 1, ($stem_base_pairs->[-1][1] - $stem_base_pairs->[-1][0] - 1));
	$motif_removal_costs->[$_] = $stem_removal_costs->[$_] + length($motif_loop_seq) * DefaultAlignParams->BASE_REMOVAL_COST;
    }

    return $motif_removal_costs;
}

sub _align_similar_len_motifs {
    my ($non_nest_stem_vertices1, $non_nest_stem_vertices2, $motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2) = @_;

    my ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs) = ({}, {}, {}, {});

    my $non_nest_stem_vertex_parts1 = _partition_by_motif_len($non_nest_stem_vertices1, $stem_graph1);
    my $non_nest_stem_vertex_parts2 = _partition_by_motif_len($non_nest_stem_vertices2, $stem_graph2);

    my @motif_lens2 = sort {$a <=> $b} keys %{$non_nest_stem_vertex_parts2};
    my $min_motif_len_part2 = $motif_lens2[0];
    my $max_motif_len_part2 = $motif_lens2[-1];

    foreach my $motif_len (sort {$a <=> $b} keys %{$non_nest_stem_vertex_parts1}) {
	foreach my $non_nest_stem_vertex1 (@{$non_nest_stem_vertex_parts1->{$motif_len}}) {
	    for (my $i = $motif_len; $i >= $min_motif_len_part2; $i--) {
		if (!exists($non_nest_stem_vertex_parts2->{$i})) {
		    next;
		}

		foreach my $non_nest_stem_vertex2 (@{$non_nest_stem_vertex_parts2->{$i}}) {
		    ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs) =
			_set_motif_align_cost_ratios($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs, $non_nest_stem_vertex1,
						     $non_nest_stem_vertex2, $motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
		}

		if ($i == $motif_len) {
		    next;
		}

		last;
	    }

	    for (my $i = $motif_len + 1; $i <= $max_motif_len_part2; $i++) {
		if (!exists($non_nest_stem_vertex_parts2->{$i})) {
		    next;
		}

		foreach my $non_nest_stem_vertex2 (@{$non_nest_stem_vertex_parts2->{$i}}) {
		    ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs) =
			_set_motif_align_cost_ratios($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs, $non_nest_stem_vertex1,
						     $non_nest_stem_vertex2, $motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
		}

		last;
	    }
	}
    }

    my @motif_lens1 = sort {$a <=> $b} keys %{$non_nest_stem_vertex_parts1};
    my $min_motif_len_part1 = $motif_lens1[0];
    my $max_motif_len_part1 = $motif_lens1[-1];

    foreach my $motif_len (sort {$a <=> $b} keys %{$non_nest_stem_vertex_parts2}) {
	foreach my $non_nest_stem_vertex2 (@{$non_nest_stem_vertex_parts2->{$motif_len}}) {
	    for (my $i = $motif_len; $i >= $min_motif_len_part1; $i--) {
		if (!exists($non_nest_stem_vertex_parts1->{$i})) {
		    next;
		}

		foreach my $non_nest_stem_vertex1 (@{$non_nest_stem_vertex_parts1->{$i}}) {
		    ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs) =
			_set_motif_align_cost_ratios($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs, $non_nest_stem_vertex1,
						     $non_nest_stem_vertex2, $motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
		}

		if ($i == $motif_len) {
		    next;
		}

		last;
	    }

	    for (my $i = $motif_len + 1; $i <= $max_motif_len_part1; $i++) {
		if (!exists($non_nest_stem_vertex_parts1->{$i})) {
		    next;
		}

		foreach my $non_nest_stem_vertex1 (@{$non_nest_stem_vertex_parts1->{$i}}) {
		    ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs) =
			_set_motif_align_cost_ratios($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs, $non_nest_stem_vertex1,
						     $non_nest_stem_vertex2, $motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);
		}

		last;
	    }
	}
    }

    return $motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs;
}

sub _set_motif_align_cost_ratios {
    my ($motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs, $non_nest_stem_vertex1, $non_nest_stem_vertex2,
	$motif_removal_costs1, $motif_removal_costs2, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2) = @_;

    my $candidate_id = $non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2;
    if (exists($motif_cost_ratios->{$candidate_id})) {
	return $motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs;
    }

    my $min_motif_removal_cost = $motif_removal_costs1->[$non_nest_stem_vertex1];
    if ($motif_removal_costs2->[$non_nest_stem_vertex2] < $min_motif_removal_cost) {
	$min_motif_removal_cost = $motif_removal_costs2->[$non_nest_stem_vertex2];
    }

    $aligned_motif_pairs = GraphMatchTools->align_motif_pair($non_nest_stem_vertex1, $non_nest_stem_vertex2, $aligned_motif_pairs, $stem_graph1,
							     $stem_graph2, $base_seq1, $base_seq2);
    my $motif_mapping_cost_ratio = $aligned_motif_pairs->{$candidate_id}->get_cost() / $min_motif_removal_cost;
    $motif_cost_ratios->{$candidate_id} = $motif_mapping_cost_ratio;

    if (!defined($min_vertex_cost_ratios1->{$non_nest_stem_vertex1}) || $motif_mapping_cost_ratio < $min_vertex_cost_ratios1->{$non_nest_stem_vertex1}) {
	$min_vertex_cost_ratios1->{$non_nest_stem_vertex1} = $motif_mapping_cost_ratio;
    }

    if (!defined($min_vertex_cost_ratios2->{$non_nest_stem_vertex2}) || $motif_mapping_cost_ratio < $min_vertex_cost_ratios2->{$non_nest_stem_vertex2}) {
	$min_vertex_cost_ratios2->{$non_nest_stem_vertex2} = $motif_mapping_cost_ratio;
    }

    return $motif_cost_ratios, $min_vertex_cost_ratios1, $min_vertex_cost_ratios2, $aligned_motif_pairs;
}

sub _partition_by_motif_len {
    my ($non_nest_stem_vertices, $stem_graph) = @_;

    my $non_nest_stem_vertex_partitions = {};

    foreach (@{$non_nest_stem_vertices}) {
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($_);
	my $stem_base_pairs = $vertex_attrs->{base_pairs};
	my $motif_len = $stem_base_pairs->[0][1] - $stem_base_pairs->[0][0];

	my $same_motif_len_vertices = $non_nest_stem_vertex_partitions->{$motif_len};
	if (defined($same_motif_len_vertices)) {
	    push @{$same_motif_len_vertices}, $_;
	}
	else {
	    $non_nest_stem_vertex_partitions->{$motif_len} = [$_];
	}
    }

    return $non_nest_stem_vertex_partitions;
}

sub _calculate_best_motif_align_cost {
    my ($base_pair_count1, $base_pair_count2, $unpaired_base_count1, $unpaired_base_count2) = @_;

    my $best_motif_align_cost;

    if ($unpaired_base_count1 >= $unpaired_base_count2) {
	$best_motif_align_cost = ($base_pair_count1 - $base_pair_count2) * DefaultAlignParams->PAIR_REMOVAL_COST +
	    ($unpaired_base_count1 - $unpaired_base_count2) * DefaultAlignParams->BASE_REMOVAL_COST;
    }
    else {
	my $unmatched_base_pair_count = $base_pair_count1 - $base_pair_count2;
	my $surplus_unpaired_base_count = $unpaired_base_count2 - $unpaired_base_count1;
	if ($surplus_unpaired_base_count >= $unmatched_base_pair_count * 2) {
	    $best_motif_align_cost = $unmatched_base_pair_count * DefaultAlignParams->BOND_BREAKING_COST +
		($surplus_unpaired_base_count - $unmatched_base_pair_count * 2) * DefaultAlignParams->BASE_REMOVAL_COST;
	}
	else {
	    my $bond_broken_base_pair_count = int($surplus_unpaired_base_count / 2);
	    $best_motif_align_cost = $bond_broken_base_pair_count * DefaultAlignParams->BOND_BREAKING_COST;

	    if ($surplus_unpaired_base_count % 2) {
		$bond_broken_base_pair_count++;
		$best_motif_align_cost += DefaultAlignParams->BOND_BREAKING_AND_BASE_REMOVAL_COST;
	    }

	    $best_motif_align_cost += ($unmatched_base_pair_count - $bond_broken_base_pair_count) * DefaultAlignParams->PAIR_REMOVAL_COST;
	}
    }

    return $best_motif_align_cost;
}

=comment
sub _align_motif_pair {
    my ($non_nest_stem_vertex1, $non_nest_stem_vertex2, $aligned_motif_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2) = @_;

    my $non_nest_stem_vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($non_nest_stem_vertex1);
    my $stem1_base_pairs = $non_nest_stem_vertex_attrs1->{base_pairs};
    my $motif1_loop_seq = substr($base_seq1, $stem1_base_pairs->[-1][0] + 1, ($stem1_base_pairs->[-1][1] - $stem1_base_pairs->[-1][0] - 1));

    my $non_nest_stem_vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($non_nest_stem_vertex2);
    my $stem2_base_pairs = $non_nest_stem_vertex_attrs2->{base_pairs};
    my $motif2_loop_seq = substr($base_seq2, $stem2_base_pairs->[-1][0] + 1, ($stem2_base_pairs->[-1][1] - $stem2_base_pairs->[-1][0] - 1));

    my $matched_stem_pair = MatchedStemPair->new($non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2);
    $matched_stem_pair->append_seq_to_right_end($motif1_loop_seq, 1, 'upstream');
    $matched_stem_pair->append_seq_to_right_end($motif2_loop_seq, 2, 'upstream');

    my $stem1_base_seq = $matched_stem_pair->get_base_seq(1, 'upstream') . $matched_stem_pair->get_base_seq(1, 'downstream');
    my $stem2_base_seq = $matched_stem_pair->get_base_seq(2, 'upstream') . $matched_stem_pair->get_base_seq(2, 'downstream');
    my ($stem1_upstream_align_bound, $stem1_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(1);
    my ($stem2_upstream_align_bound, $stem2_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(2);

    my $stem_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $matched_stem_pair->get_pseudo_base_pairs(1), $matched_stem_pair->get_pseudo_base_pairs(2),
						       $stem1_upstream_align_bound, $stem1_downstream_align_bound, $stem2_upstream_align_bound, $stem2_downstream_align_bound);
    $aligned_motif_pairs->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2} = AlignedStemPair->new($stem_align_result);

    return $aligned_motif_pairs;
}
=cut

=comment
sub _update_cost_ratio_list {
    my ($vertex, $cost_ratio_lists, $cost_ratio) = @_;

    my $max_list_size = 3;

    my $cost_ratio_list = $cost_ratio_lists->{$vertex};
    if (!defined($cost_ratio_list)) {
	$cost_ratio_lists->{$vertex} = [$cost_ratio];
	return $cost_ratio_lists;
    }

    if (@{$cost_ratio_list} == $max_list_size && $cost_ratio >= $cost_ratio_list->[-1]) {
	return $cost_ratio_lists;
    }

    push @{$cost_ratio_list}, $cost_ratio;
    my @sorted_list = sort {$a <=> $b} (@{$cost_ratio_list});
    my $sorted_list_size = @sorted_list;

    if ($sorted_list_size > $max_list_size) {
	my @new_list = @sorted_list[0..($max_list_size - 1)];
	$cost_ratio_lists->{$vertex} = \@new_list;
    }
    else {
	$cost_ratio_lists->{$vertex} = \@sorted_list;
    }

    return $cost_ratio_lists;
}

sub _get_motif_align_cost_ratio {
    my ($non_nest_stem_vertex1, $non_nest_stem_vertex2, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2, $base_seq1, $base_seq2) = @_;

    my $non_nest_stem_vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($non_nest_stem_vertex1);
    my $stem1_base_pairs = $non_nest_stem_vertex_attrs1->{base_pairs};
    my $motif1_loop_seq = substr($base_seq1, $stem1_base_pairs->[-1][0] + 1, ($stem1_base_pairs->[-1][1] - $stem1_base_pairs->[-1][0] - 1));

    my $non_nest_stem_vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($non_nest_stem_vertex2);
    my $stem2_base_pairs = $non_nest_stem_vertex_attrs2->{base_pairs};
    my $motif2_loop_seq = substr($base_seq2, $stem2_base_pairs->[-1][0] + 1, ($stem2_base_pairs->[-1][1] - $stem2_base_pairs->[-1][0] - 1));

    my $matched_stem_pair = MatchedStemPair->new($non_nest_stem_vertex_attrs1, $non_nest_stem_vertex_attrs2);
    $matched_stem_pair->append_seq_to_right_end($motif1_loop_seq, 1, 'upstream');
    $matched_stem_pair->append_seq_to_right_end($motif2_loop_seq, 2, 'upstream');

    my $stem1_base_seq = $matched_stem_pair->get_base_seq(1, 'upstream') . $matched_stem_pair->get_base_seq(1, 'downstream');
    my $stem2_base_seq = $matched_stem_pair->get_base_seq(2, 'upstream') . $matched_stem_pair->get_base_seq(2, 'downstream');
    my ($stem1_upstream_align_bound, $stem1_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(1);
    my ($stem2_upstream_align_bound, $stem2_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(2);

    my $stem_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $matched_stem_pair->get_pseudo_base_pairs(1), $matched_stem_pair->get_pseudo_base_pairs(2),
						       $stem1_upstream_align_bound, $stem1_downstream_align_bound, $stem2_upstream_align_bound, $stem2_downstream_align_bound);
    my $aligned_stem_pair = AlignedStemPair->new($stem_align_result);
    my $mapping_cost = $aligned_stem_pair->get_cost();

    my $motif1_removal_cost = $stem_removal_costs1->[$non_nest_stem_vertex1] + length($motif1_loop_seq) * AlignParams->BASE_REMOVAL_COST;
    my $motif2_removal_cost = $stem_removal_costs2->[$non_nest_stem_vertex2] + length($motif2_loop_seq) * AlignParams->BASE_REMOVAL_COST;

    return ($mapping_cost / ($motif1_removal_cost + $motif2_removal_cost));
}
=cut

=comment
sub _expand_non_nest_motif_mscs_mapping {
    my ($candidate_vertex_pairs, $non_nest_motif_mscs_mapping, $aligned_motif_pairs, $stem_graph1, $stem_graph2, $motif_removal_costs1, $motif_removal_costs2) = @_;

    $candidate_vertex_pairs = _filter_with_motif_mscs_mapping($candidate_vertex_pairs, $non_nest_motif_mscs_mapping);
    my $init_mapping_cost = GraphMatchTools->get_init_mapping_cost($candidate_vertex_pairs, $motif_removal_costs1, $motif_removal_costs2);
    my ($new_non_nest_motif_mscs_mappings, undef) = GraphMatchTools->get_mscs($candidate_vertex_pairs, $aligned_motif_pairs, $stem_graph1, $stem_graph2,
									      $init_mapping_cost, $motif_removal_costs1, $motif_removal_costs2);

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
=cut

=comment
sub _get_motif_align_trends {
    my ($candidate_vertex_pairs, $curr_motif_align_trend, $candidate_mapping_cost_ratios) = @_;

    my @sorted_vertex_pairs = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$candidate_vertex_pairs};
    my $exclusive_vertex_pair_groups = _group_exclusive_vertex_pairs(\@sorted_vertex_pairs);

    my $expanded_motif_align_trends = [$curr_motif_align_trend];

    foreach (@{$exclusive_vertex_pair_groups}) {
	my ($stem_vertex1, $stem_vertices2) = ($_->[0], $_->[1]);
	my $new_motif_align_trends = [];
	foreach my $stem_vertex2 (@{$stem_vertices2}) {
	    foreach my $expanded_motif_align_trend (@{$expanded_motif_align_trends}) {
		my @new_motif_align_trend = (@{$expanded_motif_align_trend}, ([$stem_vertex1, $stem_vertex2]));
		push @{$new_motif_align_trends}, \@new_motif_align_trend;
	    }
	}

	$expanded_motif_align_trends = $new_motif_align_trends;
    }
}

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

sub _remove_outliers {
    my $probable_vertex_pairs = shift;

    my $stem_align_shifts = [];
    foreach (@{$probable_vertex_pairs}) {
	push @{$stem_align_shifts}, ($_->[0] - $_->[1]);
    }

    my $histogram = {};
    foreach (@{$stem_align_shifts}) {
	my $count = $histogram->{$_};
	if (defined($count)) {
	    $count++;
	}
	else {
	    $count = 1;
	}

	$histogram->{$_} = $count;
    }

    my @sorted_shift_values = sort {$a <=> $b} keys %{$histogram};
    my $shift_value_count = @sorted_shift_values;
    my ($left_highest_point, $right_highest_point, $highest_point_count);

    for (my $i = 0; $i < @sorted_shift_values; $i++) {
	my $shift_value = $sorted_shift_values[$i];
	my $point_count = $histogram->{$shift_value};

	if ($point_count > $highest_point_count) {
	    $left_highest_point = $i;
	    $right_highest_point = $i;
	    $highest_point_count = $point_count;
	}
	elsif ($point_count == $highest_point_count) {
	    $right_highest_point = $i;
	}
    }

    my ($max_left_cutoff_eff, $max_right_cutoff_eff);
    my $cutoff_window = [-1, -1];
    my $cumulative_count = $histogram->{$sorted_shift_values[0]};
    for (my $i = 1; $i <= $left_highest_point; $i++) {
	my $left_cutoff_eff = ($sorted_shift_values[$i] - $sorted_shift_values[$i - 1]) * $cumulative_count;
	if (!defined($max_left_cutoff_eff) || $left_cutoff_eff > $max_left_cutoff_eff) {
	    $cutoff_window->[0] = $i;
	    $max_left_cutoff_eff = $left_cutoff_eff;
	}

	$cumulative_count += $histogram->{$sorted_shift_values[$i]};
    }

    $cumulative_count = $histogram->{$sorted_shift_values[-1]};
    for (my $i = @sorted_shift_values - 2; $i >= $right_highest_point; $i--) {
	my $right_cutoff_eff = ($sorted_shift_values[$i + 1] - $sorted_shift_values[$i]) * $cumulative_count;
	if (!defined($max_right_cutoff_eff) || $right_cutoff_eff > $max_right_cutoff_eff) {
	    $cutoff_window->[1] = $i;
	    $max_right_cutoff_eff = $right_cutoff_eff;
	}

	$cumulative_count += $histogram->{$sorted_shift_values[$i]};
    }

    my ($min_value, $max_value) = ($sorted_shift_values[$cutoff_window->[0]], $sorted_shift_values[$cutoff_window->[1]]);

    my $cleansed_results = [];
    for (my $i = 0; $i < @{$stem_align_shifts}; $i++) {
	if ($stem_align_shifts->[$i] >= $min_value && $stem_align_shifts->[$i] <= $max_value) {
	    push @{$cleansed_results}, $probable_vertex_pairs->[$i];
	}
    }

    return $cleansed_results;
}
=cut

1;
