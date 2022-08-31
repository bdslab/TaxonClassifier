package StemMatchExtension;

use strict;

my ($merge_candidate_align_cache, $merge_candidate_vertex_attrs_cache1, $merge_candidate_vertex_attrs_cache2) = ({}, {}, {});

sub extend_stem_match {
    my (undef, $mscs_mapping, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2) = @_;

#    $mscs_mapping = _reduce_mscs_mapping($mscs_mapping, $stem_graph1);

    my ($merged_stem_vertices1, $merged_stem_vertices2, $aligned_merge_stem_pairs, $merge_vertex_attrs_sets1, $merge_vertex_attrs_sets2) = ({}, {}, {}, {}, {});
    my ($merged_unmapped_vertices1, $merged_unmapped_vertices2) = ({}, {});
    my $total_reduced_cost = 0;

    my $merge_candidate_sets = _generate_merge_candidates($mscs_mapping, $stem_graph1, $stem_graph2);
    foreach (@{$merge_candidate_sets}) {
	my $merge_candidates1 = $_->[1];
	my $merge_candidates2 = $_->[2];
	if (defined($merge_candidates1->[0])) {
	    my $vertex_attrs = $stem_graph2->get_vertex_attrs_at($_->[0][1]);

	    foreach my $merge_candidate (@{$merge_candidates1}) {
		if (!_is_candidate_valid($_->[0][0], $merge_candidate, $merged_unmapped_vertices1)) {
		    next;
		}

		my $aligned_stem_pair = $aligned_stem_pairs->{$_->[0][0] . '-' . $_->[0][1]};
		my $min_align_cost = $aligned_stem_pair->get_cost();

		my $merge_candidate_id = join('_', @{$merge_candidate});

		my $merge_candidate_vertex_attrs;
		if (exists($merge_candidate_vertex_attrs_cache1->{$merge_candidate_id})) {
		    $merge_candidate_vertex_attrs = $merge_candidate_vertex_attrs_cache1->{$merge_candidate_id};
		}
		else {
		    $merge_candidate_vertex_attrs = _create_merged_stem_vertex_attrs($merge_candidate, $stem_graph1, $base_seq1);
		    $merge_candidate_vertex_attrs_cache1->{$merge_candidate_id} = $merge_candidate_vertex_attrs;
		}

		my $merge_candidate_align_result;
		if (exists($merge_candidate_align_cache->{$merge_candidate_id . '-' . $_->[0][1]})) {
		    $merge_candidate_align_result = $merge_candidate_align_cache->{$merge_candidate_id . '-' . $_->[0][1]};
		}
		else {
		    ($merge_candidate_align_result, my ($stem_align_cost1, $stem_align_cost2)) = GraphMatchTools->align_stem_pair($vertex_attrs, $merge_candidate_vertex_attrs);
		    if (defined($stem_align_cost1) && defined($stem_align_cost2)) {
			print "Inconsistent stem alignment score: $stem_align_cost1 vs $stem_align_cost2 (" . $_->[0][0] . ' (merged), ' . $_->[0][1] . ")\n";
		    }

		    $merge_candidate_align_result->self_swap();
		    $merge_candidate_align_cache->{$merge_candidate_id . '-' . $_->[0][1]} = $merge_candidate_align_result;
		}

		my $candidate_align_cost = $merge_candidate_align_result->get_cost();
		if ($candidate_align_cost < $min_align_cost) {
		    $merged_unmapped_vertices1 = _update_merged_unmapped_vertices($_->[0][0], $merged_stem_vertices1->{$_->[0][0]}, $merge_candidate,
										  $merged_unmapped_vertices1);
		    $merged_stem_vertices1->{$_->[0][0]} = $merge_candidate;
		    $aligned_merge_stem_pairs->{$_->[0][0] . '-' . $_->[0][1]} = $merge_candidate_align_result;
		    $merge_vertex_attrs_sets1->{$_->[0][0]} = $merge_candidate_vertex_attrs_cache1->{$merge_candidate_id};

		    $total_reduced_cost += ($candidate_align_cost - $min_align_cost);

		    $min_align_cost = $candidate_align_cost;
		}
	    }
	}
	elsif (defined($merge_candidates2->[0])) {
	    my $vertex_attrs = $stem_graph1->get_vertex_attrs_at($_->[0][0]);

	    foreach my $merge_candidate (@{$merge_candidates2}) {
		if (!_is_candidate_valid($_->[0][1], $merge_candidate, $merged_unmapped_vertices2)) {
		    next;
		}

		my $aligned_stem_pair = $aligned_stem_pairs->{$_->[0][0] . '-' . $_->[0][1]};
		my $min_align_cost = $aligned_stem_pair->get_cost();

		my $merge_candidate_id = join('_', @{$merge_candidate});

		my $merge_candidate_vertex_attrs;
		if (exists($merge_candidate_vertex_attrs_cache2->{$merge_candidate_id})) {
		    $merge_candidate_vertex_attrs = $merge_candidate_vertex_attrs_cache2->{$merge_candidate_id};
		}
		else {
		    $merge_candidate_vertex_attrs = _create_merged_stem_vertex_attrs($merge_candidate, $stem_graph2, $base_seq2);
		    $merge_candidate_vertex_attrs_cache2->{$merge_candidate_id} = $merge_candidate_vertex_attrs;
		}

		my $merge_candidate_align_result;
		if (exists($merge_candidate_align_cache->{$_->[0][0] . '-' . $merge_candidate_id})) {
		    $merge_candidate_align_result = $merge_candidate_align_cache->{$_->[0][0] . '-' . $merge_candidate_id};
		}
		else {
		    ($merge_candidate_align_result, my ($stem_align_cost1, $stem_align_cost2)) = GraphMatchTools->align_stem_pair($vertex_attrs, $merge_candidate_vertex_attrs);
		    if (defined($stem_align_cost1) && defined($stem_align_cost2)) {
			print "Inconsistent stem alignment score: $stem_align_cost1 vs $stem_align_cost2 (" . $_->[0][0] . ', ' . $_->[0][1] . " (merged))\n";
		    }

		    $merge_candidate_align_cache->{$_->[0][0] . '-' . $merge_candidate_id} = $merge_candidate_align_result;
		}

		my $candidate_align_cost = $merge_candidate_align_result->get_cost();
		if ($candidate_align_cost < $min_align_cost) {
		    $merged_unmapped_vertices2 = _update_merged_unmapped_vertices($_->[0][1], $merged_stem_vertices2->{$_->[0][1]}, $merge_candidate,
										  $merged_unmapped_vertices2);
		    $merged_stem_vertices2->{$_->[0][1]} = $merge_candidate;
		    $aligned_merge_stem_pairs->{$_->[0][0] . '-' . $_->[0][1]} = $merge_candidate_align_result;
		    $merge_vertex_attrs_sets2->{$_->[0][1]} = $merge_candidate_vertex_attrs_cache2->{$merge_candidate_id};

		    $total_reduced_cost += ($candidate_align_cost - $min_align_cost);

		    $min_align_cost = $candidate_align_cost;
		}
	    }
	}
    }

    return $total_reduced_cost, $merged_stem_vertices1, $merged_stem_vertices2, $aligned_merge_stem_pairs, $merge_vertex_attrs_sets1, $merge_vertex_attrs_sets2;
}

sub _is_candidate_valid {
    my ($org_mapped_vertex, $merge_candidate, $merged_unmapped_vertices) = @_;

    foreach my $candidate_vertex (@{$merge_candidate}) {
	if (exists($merged_unmapped_vertices->{$candidate_vertex}) && $merged_unmapped_vertices->{$candidate_vertex} != $org_mapped_vertex) {
	    return 0;
	}
    }

    return 1;
}

sub _update_merged_unmapped_vertices {
    my ($org_mapped_vertex, $org_merge_candidate, $new_merge_candidate, $merged_unmapped_vertices) = @_;

    foreach my $org_candidate_vertex (@{$org_merge_candidate}) {
	delete $merged_unmapped_vertices->{$org_candidate_vertex};
    }

    foreach my $new_candidate_vertex (@{$new_merge_candidate}) {
	$merged_unmapped_vertices->{$new_candidate_vertex} = $org_mapped_vertex;
    }

    return $merged_unmapped_vertices;
}

sub _generate_merge_candidates {
    my ($mscs_mapping, $stem_graph1, $stem_graph2) = @_;

    my $merge_candidate_sets = [];

    my $vertex_count1 = $stem_graph1->get_vertex_count();
    my $vertex_count2 = $stem_graph2->get_vertex_count();
    my $mapping_size = @{$mscs_mapping};

    my $vertex_attr = $stem_graph1->get_vertex_attrs_at($mscs_mapping->[0][0]);
    my $base_pair_count1 = $vertex_attr->{base_pair_count};

    $vertex_attr = $stem_graph2->get_vertex_attrs_at($mscs_mapping->[0][1]);
    my $base_pair_count2 = $vertex_attr->{base_pair_count};

    if ($mapping_size > 1) {
	my $merge_candidates1 = _generate_candidates_for_mapped_pair($mscs_mapping->[0][0], 0, $mscs_mapping->[1][0], $base_pair_count1, $base_pair_count2, $mscs_mapping, $stem_graph1, 0);
	my $merge_candidates2 = _generate_candidates_for_mapped_pair($mscs_mapping->[0][1], 0, $mscs_mapping->[1][1], $base_pair_count2, $base_pair_count1, $mscs_mapping, $stem_graph2, 1);
	if (defined($merge_candidates1->[0]) || defined($merge_candidates2->[0])) {
	    push @{$merge_candidate_sets}, [$mscs_mapping->[0], $merge_candidates1, $merge_candidates2];
	}
    }
    else {
	my $merge_candidates1 = _generate_candidates_for_mapped_pair($mscs_mapping->[0][0], 0, $vertex_count1, $base_pair_count1, $base_pair_count2, $mscs_mapping, $stem_graph1, 0);
	my $merge_candidates2 = _generate_candidates_for_mapped_pair($mscs_mapping->[0][1], 0, $vertex_count2, $base_pair_count2, $base_pair_count1, $mscs_mapping, $stem_graph2, 1);
	if (defined($merge_candidates1->[0]) || defined($merge_candidates2->[0])) {
	    push @{$merge_candidate_sets}, [$mscs_mapping->[0], $merge_candidates1, $merge_candidates2];
	}
    }

    for (my $i = 1; $i < ($mapping_size - 1); $i++) {
	$vertex_attr = $stem_graph1->get_vertex_attrs_at($mscs_mapping->[$i][0]);
	my $base_pair_count1 = $vertex_attr->{base_pair_count};

	$vertex_attr = $stem_graph2->get_vertex_attrs_at($mscs_mapping->[$i][1]);
	my $base_pair_count2 = $vertex_attr->{base_pair_count};

	my $merge_candidates1 = _generate_candidates_for_mapped_pair($mscs_mapping->[$i][0], $mscs_mapping->[$i - 1][0] + 1, $mscs_mapping->[$i + 1][0],
								     $base_pair_count1, $base_pair_count2, $mscs_mapping, $stem_graph1, 0);
	my $merge_candidates2 = _generate_candidates_for_mapped_pair($mscs_mapping->[$i][1], $mscs_mapping->[$i - 1][1] + 1, $mscs_mapping->[$i + 1][1],
								     $base_pair_count2, $base_pair_count1, $mscs_mapping, $stem_graph2, 1);
	if (defined($merge_candidates1->[0]) || defined($merge_candidates2->[0])) {
	    push @{$merge_candidate_sets}, [$mscs_mapping->[$i], $merge_candidates1, $merge_candidates2];
	}
    }

    if ($mapping_size > 1) {
	$vertex_attr = $stem_graph1->get_vertex_attrs_at($mscs_mapping->[-1][0]);
	my $base_pair_count1 = $vertex_attr->{base_pair_count};

	$vertex_attr = $stem_graph2->get_vertex_attrs_at($mscs_mapping->[-1][1]);
	my $base_pair_count2 = $vertex_attr->{base_pair_count};
	
	my $merge_candidates1 = _generate_candidates_for_mapped_pair($mscs_mapping->[-1][0], $mscs_mapping->[-2][0] + 1, $vertex_count1, $base_pair_count1, $base_pair_count2,
								     $mscs_mapping, $stem_graph1, 0);
	my $merge_candidates2 = _generate_candidates_for_mapped_pair($mscs_mapping->[-1][1], $mscs_mapping->[-2][1] + 1, $vertex_count2, $base_pair_count2, $base_pair_count1,
								     $mscs_mapping, $stem_graph2, 1);
	if (defined($merge_candidates1->[0]) || defined($merge_candidates2->[0])) {
	    push @{$merge_candidate_sets}, [$mscs_mapping->[-1], $merge_candidates1, $merge_candidates2];
	}
    }

    return $merge_candidate_sets;
}

sub _generate_candidates_for_mapped_pair {
    my ($mapped_stem_vertex, $search_range_start, $search_range_end, $mapped_vertex_base_pair_count, $mapped_partner_base_pair_count, $mscs_mapping, $stem_graph, $mapping_index) = @_;

    my $org_base_pair_count_diff = $mapped_vertex_base_pair_count - $mapped_partner_base_pair_count;
    if ($org_base_pair_count_diff >= 0) {
	return [];
    }

    if ($mapped_stem_vertex == $search_range_start && $mapped_stem_vertex == $search_range_end - 1) {
	return [];
    }

    my $candidates = [];

    my $limit = abs($org_base_pair_count_diff);
    my $base_pair_count_cache = {};

    my $unmapped_nest_stem_vertex_pool = _get_unmapped_nest_stem_vertex_pool($mapped_stem_vertex, $search_range_start, $search_range_end, $mscs_mapping, $stem_graph, $mapping_index);
    my $unmapped_nest_stem_vertex_sets = _group_unmapped_nest_stem_vertices($unmapped_nest_stem_vertex_pool, $mapped_stem_vertex, $stem_graph);

    foreach (@{$unmapped_nest_stem_vertex_sets}) {
	my ($prev_unmapped_vertices, $succ_unmapped_vertices) = @{$_};

	my $dp_matrix = [[$org_base_pair_count_diff]];

	for (my $i = 0; $i < @{$prev_unmapped_vertices}; $i++) {
	    my $base_pair_count = $base_pair_count_cache->{$prev_unmapped_vertices->[-1 - $i]};
	    if (!defined($base_pair_count)) {
		my $vertex_attrs = $stem_graph->get_vertex_attrs_at($prev_unmapped_vertices->[-1 - $i]);
		$base_pair_count = $vertex_attrs->{base_pair_count};
		$base_pair_count_cache->{$prev_unmapped_vertices->[-1 - $i]} = $base_pair_count;
	    }

	    $dp_matrix->[$i + 1][0] = $dp_matrix->[$i][0] + $base_pair_count;
	    if ($dp_matrix->[$i + 1][0] > 0) {
		if ($dp_matrix->[$i + 1][0] > $limit) {
		    undef $dp_matrix->[$i + 1][0];
		}

		last;
	    }
	}

	for (my $i = 0; $i < @{$succ_unmapped_vertices}; $i++) {
	    my $base_pair_count = $base_pair_count_cache->{$succ_unmapped_vertices->[$i]};
	    if (!defined($base_pair_count)) {
		my $vertex_attrs = $stem_graph->get_vertex_attrs_at($succ_unmapped_vertices->[$i]);
		$base_pair_count = $vertex_attrs->{base_pair_count};
		$base_pair_count_cache->{$succ_unmapped_vertices->[$i]} = $base_pair_count;
	    }

	    $dp_matrix->[0][$i + 1] = $dp_matrix->[0][$i] + $base_pair_count;
	    if ($dp_matrix->[0][$i + 1] > 0) {
		if ($dp_matrix->[0][$i + 1] > $limit) {
		    undef $dp_matrix->[0][$i + 1];
		}

		last;
	    }
	}

	for (my $i = 0; $i < @{$prev_unmapped_vertices}; $i++) {
	    for (my $j = 0; $j < @{$succ_unmapped_vertices}; $j++) {
		if (!defined($dp_matrix->[$i][$j + 1]) || !defined($dp_matrix->[$i + 1][$j])) {
		    next;
		}

		if ($dp_matrix->[$i][$j + 1] > 0 || $dp_matrix->[$i + 1][$j] > 0) {
		    next;
		}

		$dp_matrix->[$i + 1][$j + 1] = $dp_matrix->[$i + 1][$j] + $base_pair_count_cache->{$succ_unmapped_vertices->[$j]};
		if ($dp_matrix->[$i + 1][$j + 1] > 0) {
		    if ($dp_matrix->[$i + 1][$j + 1] > $limit) {
			undef $dp_matrix->[$i + 1][$j + 1];
		    }

		    last;
		}
	    }
	}

	for (my $i = 0; $i < @{$dp_matrix}; $i++) {
	    for (my $j = 0; $j < @{$dp_matrix->[$i]}; $j++) {
		if ($i == 0 && $j == 0) {
		    next;
		}

		if (defined($dp_matrix->[$i][$j])) {
		    my $candidate = [];
		    if ($i > 0) {
			for (my $k = $i - 1; $k >= 0; $k--) {
			    push @{$candidate}, $prev_unmapped_vertices->[-1 - $k];
			}
		    }

		    push @{$candidate}, $mapped_stem_vertex;

		    if ($j > 0) {
			for (my $k = 0; $k <= $j - 1; $k++) {
			    push @{$candidate}, $succ_unmapped_vertices->[$k];
			}
		    }

		    push @{$candidates}, $candidate;
		}
	    }
	}
    }

    return $candidates;
}

sub _get_unmapped_nest_stem_vertex_pool {
    my ($mapped_stem_vertex, $search_range_start, $search_range_end, $mscs_mapping, $stem_graph, $mapping_index) = @_;

    my $mapped_vertex_ext_edge_label_token = _get_external_edge_label_token($mapped_stem_vertex, $mapped_stem_vertex, $mscs_mapping, $stem_graph, $mapping_index);

    my $unmapped_nest_stem_vertex_pool = [];

    for (my $i = $search_range_start; $i < $mapped_stem_vertex; $i++) {
	if ($stem_graph->is_nest_stem_vertex($i) && $stem_graph->get_edge_label($i, $mapped_stem_vertex) eq 'N') {
	    if (_get_external_edge_label_token($i, $mapped_stem_vertex, $mscs_mapping, $stem_graph, $mapping_index) eq $mapped_vertex_ext_edge_label_token) {
		push @{$unmapped_nest_stem_vertex_pool}, $i;
	    }
	}
    }

    if ($stem_graph->is_nest_stem_vertex($mapped_stem_vertex)) {
	for (my $i = $mapped_stem_vertex + 1; $i < $search_range_end; $i++) {
	    if ($stem_graph->is_nest_stem_vertex($i) && $stem_graph->get_edge_label($mapped_stem_vertex, $i) eq 'N') {
		if (_get_external_edge_label_token($i, $mapped_stem_vertex, $mscs_mapping, $stem_graph, $mapping_index) eq $mapped_vertex_ext_edge_label_token) {
		    push @{$unmapped_nest_stem_vertex_pool}, $i;
		}
	    }
	}
    }

    return $unmapped_nest_stem_vertex_pool;
}

sub _get_external_edge_label_token {
    my ($stem_vertex, $skip_vertex, $mscs_mapping, $stem_graph, $mapping_index) = @_;

    my $external_edge_label_token = '';

    foreach (@{$mscs_mapping}) {
	my $mapped_stem_vertex = $_->[$mapping_index];
	if ($mapped_stem_vertex == $skip_vertex) {
	    next;
	}

	if ($mapped_stem_vertex < $stem_vertex) {
	    $external_edge_label_token = $external_edge_label_token . $stem_graph->get_edge_label($mapped_stem_vertex, $stem_vertex);
	}
	elsif ($mapped_stem_vertex > $stem_vertex) {
	    $external_edge_label_token = $external_edge_label_token . lc($stem_graph->get_edge_label($stem_vertex, $mapped_stem_vertex));
	}
    }

    return $external_edge_label_token;
}

sub _group_unmapped_nest_stem_vertices {
    my ($unmapped_nest_stem_vertex_pool, $mapped_stem_vertex, $stem_graph) = @_;

    my $non_crossing_unmapped_vertex_sets = [];

    if (defined($unmapped_nest_stem_vertex_pool->[0])) {
	$non_crossing_unmapped_vertex_sets = [[$unmapped_nest_stem_vertex_pool->[0]]];

	for (my $i = 1; $i < @{$unmapped_nest_stem_vertex_pool}; $i++) {
	    my $new_vertex_sets = [];
	    for (my $j = 0; $j < @{$non_crossing_unmapped_vertex_sets}; $j++) {
		my $new_vertex_set = [];
		my $is_vertex_removed = 0;
		foreach (@{$non_crossing_unmapped_vertex_sets->[$j]}) {
		    if ($stem_graph->get_edge_label($_, $unmapped_nest_stem_vertex_pool->[$i]) eq 'N') {
			push @{$new_vertex_set}, $_;
		    }
		    else {
			$is_vertex_removed = 1;
		    }
		}

		push @{$new_vertex_set}, $unmapped_nest_stem_vertex_pool->[$i];

		if ($is_vertex_removed) {
		    push @{$new_vertex_sets}, $new_vertex_set;
		}
		else {
		    $non_crossing_unmapped_vertex_sets->[$j] = $new_vertex_set;
		}
	    }

	    push @{$non_crossing_unmapped_vertex_sets}, @{$new_vertex_sets};
	}
    }

    my $unmapped_nest_stem_vertex_sets = [];
    foreach my $non_crossing_unmapped_vertex_set (@{$non_crossing_unmapped_vertex_sets}) {
	my $nest_stem_vertex_set = [[], []];
	foreach (@{$non_crossing_unmapped_vertex_set}) {
	    if ($_ < $mapped_stem_vertex) {
		push @{$nest_stem_vertex_set->[0]}, $_;
	    }
	    else {
		push @{$nest_stem_vertex_set->[1]}, $_;
	    }
	}

	push @{$unmapped_nest_stem_vertex_sets}, $nest_stem_vertex_set;
    }

    return $unmapped_nest_stem_vertex_sets;
}

sub _create_merged_stem_vertex_attrs {
    my ($merge_candidate, $stem_graph, $base_seq) = @_;

    my $merged_stem_vertex_attrs = {};
    my $merged_stem_base_pairs = _merge_stem_base_pairs($merge_candidate, $stem_graph);
    my $merged_stem_upstream_start = $merged_stem_base_pairs->[0][0];
    my $merged_stem_downstream_start = $merged_stem_base_pairs->[-1][1];
    my $merged_stem_upstream_len = $merged_stem_base_pairs->[-1][0] - $merged_stem_upstream_start + 1;
    my $merged_stem_downstream_len = $merged_stem_base_pairs->[0][1] - $merged_stem_downstream_start + 1;

    my $vertex_attrs = {};
    $vertex_attrs->{base_pair_count} = @{$merged_stem_base_pairs};
    $vertex_attrs->{base_pairs} = $merged_stem_base_pairs;
    $vertex_attrs->{upstream_base_seq} = substr($base_seq, $merged_stem_upstream_start, $merged_stem_upstream_len);
    $vertex_attrs->{downstream_base_seq} = substr($base_seq, $merged_stem_downstream_start, $merged_stem_downstream_len);
    $vertex_attrs->{upstream_start} = $merged_stem_upstream_start;
    $vertex_attrs->{upstream_length} = $merged_stem_upstream_len;
    $vertex_attrs->{downstream_start} = $merged_stem_downstream_start;
    $vertex_attrs->{downstream_length} = $merged_stem_downstream_len;
    $vertex_attrs->{pseudo_base_pairs} = GraphMatchTools->convert_to_pseudo_base_pairs($merged_stem_base_pairs);

    return $vertex_attrs;
}

sub _merge_stem_base_pairs {
    my ($merge_candidate, $stem_graph) = @_;

    my $merged_stem_base_pairs = [];
    foreach (@{$merge_candidate}) {
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($_);
	push @{$merged_stem_base_pairs}, @{$vertex_attrs->{base_pairs}};
    }

    return $merged_stem_base_pairs;
}

=comment
sub merge_mapped_stems {
    my (undef, $mscs_mappings, $nest_stem_vertex_neighbor_profiles1, $nest_stem_vertex_neighbor_profiles2, $stem_graph1, $stem_graph2) = @_;

    foreach my $mscs_mapping (@{$mscs_mappings}) {
	my ($collapsible_stem_vertex_groups1, $collapsible_stem_group2) = _collapse_mapped_stem_vertices($mscs_mapping, $nest_stem_vertex_neighbor_profiles1, $stem_graph1);
    }
}

sub _collapse_mapped_stem_vertices {
    my ($mscs_mapping, $nest_stem_vertex_neighbor_profiles1, $stem_graph1) = @_;

    my $mapped_stem_vertices = {};
    foreach (@{$mscs_mapping}) {
	$mapped_stem_vertices->{$_->[0]} = 1;
    }

    my $eff_neighbor_profiles = {};
    foreach (@{$mscs_mapping}) {
	my ($eff_nested_neighbors, $eff_prev_crossed_neighbor_str, $eff_succ_crossed_neighbor_str) = _get_effective_neighbor_profile($nest_stem_vertex_neighbor_profiles1->{$_->[0]}, $mapped_stem_vertices);
	$eff_neighbor_profiles->{$_->[0]} = [$eff_nested_neighbors, $eff_prev_crossed_neighbor_str, $eff_succ_crossed_neighbor_str];
    }

    my ($collapsible_stem_vertex_groups1, $collapsible_stem_vertex_groups2) = ([], []);
    for (my $i = 0; $i < @{$mscs_mapping} - 1; $i++) {
	my $collapsible_stem_vertex_group1 = [$mscs_mapping->[$i][0]];
	push @{$collapsible_stem_vertex_groups1}, $collapsible_stem_vertex_group1;

	my $collapsible_stem_vertex_group2 = [$mscs_mapping->[$i][1]];
	push @{$collapsible_stem_vertex_groups2}, $collapsible_stem_vertex_group2;

	if (!$stem_graph1->is_nest_stem_vertex($mscs_mapping->[$i][0])) {
	    next;
	}

	my $mapped_parent_stem_vertex = $mscs_mapping->[$i][0];
	my $parent_eff_neighbor_profile = $eff_neighbor_profiles->{$mapped_parent_stem_vertex};

	for (my $j = $i + 1; $j < @{$mscs_mapping}; $j++) {
	    if ($stem_graph1->get_edge_label($mapped_parent_stem_vertex, $mscs_mapping->[$j][0]) ne 'N') {
		last;
	    }

	    my $mapped_child_stem_vertex = $mscs_mapping->[$j][0];
	    my $child_eff_neighbor_profile = $eff_neighbor_profiles->{$mapped_child_stem_vertex};

	    if ($parent_eff_neighbor_profile->[1] ne $child_eff_neighbor_profile->[1] || $parent_eff_neighbor_profile->[2] ne $child_eff_neighbor_profile->[2]) {
		last;
	    }

	    my $is_nested_neighbor_matched = 1;
	    foreach (@{$parent_eff_neighbor_profile->[0]}) {
		my $edge_label;
		if ($_ < $mapped_child_stem_vertex) {
		    $edge_label = $stem_graph1->get_edge_label($_, $mapped_child_stem_vertex);
		}
		elsif ($_ > $mapped_child_stem_vertex) {
		    $edge_label = $stem_graph1->get_edge_label($mapped_child_stem_vertex, $_);
		}
		else {
		    next;
		}

		if ($edge_label ne 'N') {
		    $is_nested_neighbor_matched = 0;
		    last;
		}
	    }

	    if (!$is_nested_neighbor_matched) {
		last;
	    }

	    push @{$collapsible_stem_vertex_group1}, $mapped_child_stem_vertex;
	    push @{$collapsible_stem_vertex_group2}, $mscs_mapping->[$j][1];
	    $i++;
	}
    }

    return $collapsible_stem_vertex_groups1, $collapsible_stem_vertex_groups2;
}

sub _get_effective_neighbor_profile {
    my ($stem_vertex_neighbor_profile, $mapped_stem_vertices) = @_;

    my ($eff_nested_neighbors, $eff_prev_crossed_neighbors, $eff_succ_crossed_neighbors) = ([], [], []);
    my $eff_neighbor_profile = [$eff_nested_neighbors, $eff_prev_crossed_neighbors, $eff_succ_crossed_neighbors];

    for (my $i = 0; $i < 2; $i++) {
	for (my $j = 0; $j < @{$stem_vertex_neighbor_profile->[$i]}; $j++) {
	    foreach (@{$stem_vertex_neighbor_profile->[$i][$j]}) {
		if (exists($mapped_stem_vertices->{$_})) {
		    push @{$eff_neighbor_profile->[$j]}, $_;
		}
	    }
	}
    }

    for (my $i = 0; $i < @{$eff_neighbor_profile}; $i++) {
	my @sorted_eff_neighbors = sort {$a <=> $b} @{$eff_neighbor_profile->[$i]};
	$eff_neighbor_profile->[$i] = \@sorted_eff_neighbors;
    }

    return $eff_nested_neighbors, join('-', @{$eff_prev_crossed_neighbors}), join('-', @{$eff_succ_crossed_neighbors});
}

sub _interpolate_collapsible_stem_vertices {
    my ($collapsible_stem_vertex_groups, $stem_graph) = @_;

    my $extended_collapsible_stem_vertex_groups = [];
    foreach my $collapsible_stem_vertex_group (@{$collapsible_stem_vertex_groups}) {
	my $extended_collapsible_stem_vertex_group = [];
	push @{$extended_collapsible_stem_vertex_groups}, $extended_collapsible_stem_vertex_group;

	for (my $i = 1; $i < @{$collapsible_stem_vertex_group}; $i++) {
	    push @{$extended_collapsible_stem_vertex_group}, $collapsible_stem_vertex_group->[$i - 1];

	    for my ($j = $collapsible_stem_vertex_group->[$i - 1] + 1; $j < $collapsible_stem_vertex_group->[$i]; $j++) {
		if ($stem_graph->is_nest_stem_vertex($j)) {
		    if ($stem_graph->get_edge_label($collapsible_stem_vertex_group->[$i - 1], $j) eq 'N' &&
			$stem_graph->get_edge_label($j, $collapsible_stem_vertex_group->[$i]) eq 'N') {
			push @{$extended_collapsible_stem_vertex_group}, $j;
		    }
		}
	    }

	    push @{$extended_collapsible_stem_vertex_group}, $collapsible_stem_vertex_group->[$i];
	}
    }

    return $extended_collapsible_stem_vertex_groups;
}

sub _exterpolate_collapsible_stem_vertices {
    my ($collapsible_stem_vertex_groups, $stem_graph) = @_;


}
=cut

1;
