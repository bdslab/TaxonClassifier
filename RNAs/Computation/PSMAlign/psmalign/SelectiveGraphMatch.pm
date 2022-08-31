package SelectiveGraphMatch;

use strict;

sub match_nest_stem_vertices {
    my (undef, $best_non_nest_stem_vertex_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2, $align_params) = @_;

    my ($matched_non_nest_stem_vertices1, $matched_non_nest_stem_vertices2, $matched_vertex_pair_hash) = ({}, {}, {});
    foreach (@{$best_non_nest_stem_vertex_pairs}) {
	$matched_non_nest_stem_vertices1->{$_->[0]} = 1;
	$matched_non_nest_stem_vertices2->{$_->[1]} = 1;
	$matched_vertex_pair_hash->{$_->[0] . '-' . $_->[1]} = 1;
    }

    my $nest_stem_vertex_neighbor_profiles1 = _compute_nest_stem_vertex_neighbor_profiles($stem_graph1);
    my $nest_stem_vertex_neighbor_profiles2 = _compute_nest_stem_vertex_neighbor_profiles($stem_graph2);

    my ($nest_stem_match_candidates1, $nest_stem_match_candidates2) = ({}, {});
    my $same_neighbor_match_vertex_pairs = {};
    my (@unique_matched_neighbors1, @unique_matched_neighbors2);
#    my $threshold = int(@{$best_non_nest_stem_vertex_pairs} / 2 + 0.5);

    while (my ($nest_stem_vertex1, $neighbor_profile1) = each %{$nest_stem_vertex_neighbor_profiles1}) {
	if (exists($matched_non_nest_stem_vertices1->{$nest_stem_vertex1})) {
	    next;
	}

	while (my ($nest_stem_vertex2, $neighbor_profile2) = each %{$nest_stem_vertex_neighbor_profiles2}) {
	    if (exists($matched_non_nest_stem_vertices2->{$nest_stem_vertex2})) {
		next;
	    }

	    my $matched_neighbor_pairs_str = '';
	    my $matched_neighbor_pair_count = 0;

#	    my $is_contain_best_non_nest_vertex_pair = 0;

	    for (my $i = 0; $i < @{$neighbor_profile1}; $i++) {
		my ($matched_neighbors1, $matched_neighbors2) = ({}, {});

		foreach my $non_nest_stem_vertex1 (@{$neighbor_profile1->[$i]}) {
		    foreach my $non_nest_stem_vertex2 (@{$neighbor_profile2->[$i]}) {
			if (exists($matched_vertex_pair_hash->{$non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2})) {
			    if (!exists($matched_neighbors1->{$non_nest_stem_vertex1}) && !exists($matched_neighbors2->{$non_nest_stem_vertex2})) {
				$matched_neighbors1->{$non_nest_stem_vertex1} = 1;
				$matched_neighbors2->{$non_nest_stem_vertex2} = 1;
				$matched_neighbor_pair_count++;
			    }
#			    $matched_neighbor_pairs_str = $matched_neighbor_pairs_str . '_' . $non_nest_stem_vertex1 . '-' . $non_nest_stem_vertex2;
#			    $matched_neighbor_pair_count++;
			}
		    }
		}

		foreach my $non_nest_stem_vertex1 (@{$neighbor_profile1->[$i]}) {
		    if (!exists($matched_neighbors1->{$non_nest_stem_vertex1}) && exists($matched_non_nest_stem_vertices1->{$non_nest_stem_vertex1})) {
			$matched_neighbor_pair_count -= 0.5;
		    }
		}

		foreach my $non_nest_stem_vertex2 (@{$neighbor_profile2->[$i]}) {
		    if (!exists($matched_neighbors2->{$non_nest_stem_vertex2}) && exists($matched_non_nest_stem_vertices2->{$non_nest_stem_vertex2})) {
			$matched_neighbor_pair_count -= 0.5;
		    }
		}

#		$matched_neighbor_pairs_str = $matched_neighbor_pairs_str . '|';

		@unique_matched_neighbors1 = sort {$a <=> $b} keys %{$matched_neighbors1};
		@unique_matched_neighbors2 = sort {$a <=> $b} keys %{$matched_neighbors2};
#		my $unique_matched_neighbor_count1 = scalar(@unique_matched_neighbors1);
#		my $unique_matched_neighbor_count2 = scalar(@unique_matched_neighbors2);
#		$matched_neighbor_count += ($unique_matched_neighbor_count1, $unique_matched_neighbor_count2)[$unique_matched_neighbor_count1 > $unique_matched_neighbor_count2];
		$matched_neighbor_pairs_str = $matched_neighbor_pairs_str . join('_', @unique_matched_neighbors1) . '-' . join('_', @unique_matched_neighbors2) .  '|';

#		if ($i == 0 && $matched_neighbor_pair_count > 0) {
#		    $is_contain_best_non_nest_vertex_pair = 1;
#		}
	    }

=comment
	    if ($matched_neighbor_pairs_str ne '-|-|-|') {
		print "$nest_stem_vertex1, $nest_stem_vertex2: $matched_neighbor_pairs_str\n";
	    }
=cut

#	    if (!$is_contain_best_non_nest_vertex_pair || ($matched_neighbor_pair_count < $threshold)) {
#		next;
#	    }

	    if ($matched_neighbor_pair_count <= 0) {
#		print "$nest_stem_vertex1, $nest_stem_vertex2\n";
		next;
	    }

	    if (exists($same_neighbor_match_vertex_pairs->{$matched_neighbor_pairs_str})) {
		my $vertex_pairs = $same_neighbor_match_vertex_pairs->{$matched_neighbor_pairs_str};
		push @{$vertex_pairs}, [$nest_stem_vertex1, $nest_stem_vertex2];
	    }
	    else {
		$same_neighbor_match_vertex_pairs->{$matched_neighbor_pairs_str} = [[$nest_stem_vertex1, $nest_stem_vertex2]];
	    }

	    if (!exists($nest_stem_match_candidates1->{$nest_stem_vertex1})) {
		$nest_stem_match_candidates1->{$nest_stem_vertex1} = [[$matched_neighbor_pair_count, $nest_stem_vertex2]];
	    }
	    else {
		push @{$nest_stem_match_candidates1->{$nest_stem_vertex1}}, [$matched_neighbor_pair_count, $nest_stem_vertex2];
	    }

	    if (!exists($nest_stem_match_candidates2->{$nest_stem_vertex2})) {
		$nest_stem_match_candidates2->{$nest_stem_vertex2} = [[$matched_neighbor_pair_count, $nest_stem_vertex1]];
	    }
	    else {
		push @{$nest_stem_match_candidates2->{$nest_stem_vertex2}}, [$matched_neighbor_pair_count, $nest_stem_vertex1];
	    }
	}
    }

    my $top_k_matches1 = _get_top_k_matches($nest_stem_match_candidates1, $align_params->get_k());
    my $top_k_matches2 = _get_top_k_matches($nest_stem_match_candidates2, $align_params->get_k());

#    my $matched_nest_vertex_pairs = [];
    my $candidate_nest_vertex_pairs_groups = [];
    my ($graph1_candidate_nest_vertices, $graph2_candidate_nest_vertices) = ({}, {});
    my $aligned_stem_pairs = {};

    foreach (values %{$same_neighbor_match_vertex_pairs}) {
	my $grouped_candidate_nest_vertex_pairs = [];
	foreach (@{$_}) {
	    if (exists($top_k_matches1->{$_->[0] . '-' . $_->[1]}) && exists($top_k_matches2->{$_->[1] . '-' . $_->[0]})) {
		push @{$grouped_candidate_nest_vertex_pairs}, $_;
		$graph1_candidate_nest_vertices->{$_->[0]} = 1;
		$graph2_candidate_nest_vertices->{$_->[1]} = 1;
	    }
	}

	push @{$candidate_nest_vertex_pairs_groups}, $grouped_candidate_nest_vertex_pairs;
    }

    my $candidate_nest_vertex_profiles1 = _compute_candidate_nest_vertex_profiles($graph1_candidate_nest_vertices, $stem_graph1);
    my $candidate_nest_vertex_profiles2 = _compute_candidate_nest_vertex_profiles($graph2_candidate_nest_vertices, $stem_graph2);
    my ($filtered_candidate_nest_vertex_pairs, $aligned_stem_pairs) =
	_filter_candidate_nest_vertex_pairs($candidate_nest_vertex_pairs_groups, $candidate_nest_vertex_profiles1,$candidate_nest_vertex_profiles2,
					    $aligned_stem_pairs, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2);

#    my ($mscs_mappings, undef, $aligned_stem_pairs) = GraphMatchTools->expand_mscs_mapping($matched_nest_vertex_pairs, [], {}, $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2);
    (my $mscs_mappings, $aligned_stem_pairs) = GraphMatchTools->expand_mscs_mapping($filtered_candidate_nest_vertex_pairs, [], $aligned_stem_pairs,
										    $stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2);

    if (!defined($mscs_mappings->[0])) {
	$mscs_mappings = [[]];
    }

=comment
    print "After nest stem match:\n";
    foreach my $mscs_mapping (@{$mscs_mappings}) {
	foreach (@{$mscs_mapping}) {
	    print '[' . $_->[0] . ', ' . $_->[1] . '] ';
	}
	print "\n";
    }
=cut
#    my @probable_stem_vertex_pairs = (@{$local_aligned_vertex_pairs}, @{$matched_nest_stem_vertex_pairs});

#    return \@probable_stem_vertex_pairs, $aligned_stem_pairs;
    return $mscs_mappings, $aligned_stem_pairs;
}

sub _compute_nest_stem_vertex_neighbor_profiles {
    my $stem_graph = shift;

    my $vertex_neighbor_profiles = {};

    my $vertex_count = $stem_graph->get_vertex_count();
    my $nest_stem_vertices = $stem_graph->get_nest_stem_vertices();
    foreach (@{$nest_stem_vertices}) {
#	my ($nested_non_nest_stem_neighbors, $prev_crossed_non_nest_stem_neighbors, $succ_crossed_non_nest_stem_neighbors, $prev_parallel_non_nest_stem_neighbors, 
#	    $succ_parallel_non_nest_stem_neighbors) = ([], [], [], [], []);
	my ($nested_non_nest_stem_neighbors, $prev_crossed_non_nest_stem_neighbors, $succ_crossed_non_nest_stem_neighbors) = ([], [], []);
	for (my $i = $_ - 1; $i >= 0; $i--) {
	    my $edge_label = $stem_graph->get_edge_label($i, $_);
	    if ($edge_label eq 'P') {
		last;
#		push @{$prev_parallel_non_nest_stem_neighbors}, $i;
	    }

	    if ($edge_label eq 'K') {
		if (!$stem_graph->is_nest_stem_vertex($i)) {
		    push @{$prev_crossed_non_nest_stem_neighbors}, $i;
		}
	    }
	}

	for (my $i = $_ + 1; $i < $vertex_count; $i++) {
	    my $edge_label = $stem_graph->get_edge_label($_, $i);
	    if ($edge_label eq 'P') {
		last;
#		push @{$succ_parallel_non_nest_stem_neighbors}, $i;
	    }

	    if ($edge_label eq 'N') {
		if (!$stem_graph->is_nest_stem_vertex($i)) {
		    push @{$nested_non_nest_stem_neighbors}, $i;
		}
	    }
	    elsif ($edge_label eq 'K') {
		if (!$stem_graph->is_nest_stem_vertex($i)) {
		    push @{$succ_crossed_non_nest_stem_neighbors}, $i;
		}
	    }
	}

#	my $non_parallel_neighbor_count = @{$nested_non_nest_stem_neighbors} + @{$prev_crossed_non_nest_stem_neighbors} + @{$succ_crossed_non_nest_stem_neighbors} +
#	    @{$prev_parallel_non_nest_stem_neighbors} + @{$succ_parallel_non_nest_stem_neighbors};

#	$vertex_neighbor_profiles->{$_} = [$non_parallel_neighbor_count, $nested_non_nest_stem_neighbors, $prev_crossed_non_nest_stem_neighbors, $succ_crossed_non_nest_stem_neighbors];

	$vertex_neighbor_profiles->{$_} = [$nested_non_nest_stem_neighbors, $prev_crossed_non_nest_stem_neighbors, $succ_crossed_non_nest_stem_neighbors];

#	$vertex_neighbor_profiles->{$_} = [$nested_non_nest_stem_neighbors, $prev_crossed_non_nest_stem_neighbors, $succ_crossed_non_nest_stem_neighbors,
#	    $prev_parallel_non_nest_stem_neighbors, $succ_parallel_non_nest_stem_neighbors];

    }

    return $vertex_neighbor_profiles;
}

sub _get_top_k_matches {
    my ($vertex_matches, $k) = @_;

    my $top_k_matches = {};

    while (my ($stem_vertex, $matched_vertices) = each %{$vertex_matches}) {
	my $selected_vertex_count = 0;
	my $last_match_score;

	foreach (sort {$b->[0] <=> $a->[0] || $a->[1] <=> $b->[1]} @{$matched_vertices}) {
	    if (!defined($last_match_score) || $_->[0] == $last_match_score || $selected_vertex_count < $k) {
		$top_k_matches->{$stem_vertex . '-' . $_->[1]} = [$stem_vertex, $_->[1]];
		$selected_vertex_count++;
		$last_match_score = $_->[0];
	    }
	}
    }

    return $top_k_matches;
}

sub _compute_candidate_nest_vertex_profiles {
    my ($candidate_nest_vertices, $stem_graph) = @_;

    my @sorted_nest_vertices = sort {$a <=> $b} keys %{$candidate_nest_vertices};
    my $candidate_nest_vertex_profiles = {};

    for (my $i = 0; $i < @sorted_nest_vertices; $i++) {
	my $curr_candidate_nest_vertex = $sorted_nest_vertices[$i];

	if ($stem_graph->get_stem_type_at($curr_candidate_nest_vertex) eq 'M') {
	    $candidate_nest_vertex_profiles->{$curr_candidate_nest_vertex} = [[], []];
	    next;
	}

	my ($prev_crossed_nest_stem_neighbors, $succ_crossed_nest_stem_neighbors) = ([], []);

	for (my $j = $i - 1; $j >= 0; $j--) {
	    if ($stem_graph->get_edge_label($sorted_nest_vertices[$j], $curr_candidate_nest_vertex) eq 'K') {
		push @{$prev_crossed_nest_stem_neighbors}, $sorted_nest_vertices[$j];
	    }
	}

	for (my $j = $i + 1; $j < @sorted_nest_vertices; $j++) {
	    my $edge_label = $stem_graph->get_edge_label($curr_candidate_nest_vertex, $sorted_nest_vertices[$j]);
	    if ($edge_label eq 'P') {
		last;
	    }

	    if ($edge_label eq 'K') {
		push @{$succ_crossed_nest_stem_neighbors}, $sorted_nest_vertices[$j];
	    }
	}

	$candidate_nest_vertex_profiles->{$curr_candidate_nest_vertex} = [$prev_crossed_nest_stem_neighbors, $succ_crossed_nest_stem_neighbors];
    }

    return $candidate_nest_vertex_profiles;
}

sub _filter_candidate_nest_vertex_pairs {
    my ($candidate_nest_vertex_pairs_groups, $candidate_nest_vertex_profiles1, $candidate_nest_vertex_profiles2, $aligned_stem_pairs,
	$stem_graph1, $stem_graph2, $stem_removal_costs1, $stem_removal_costs2) = @_;

    my $matched_nest_vertex_candidate_pairs = [];

    foreach my $candidate_nest_vertex_pairs_group (@{$candidate_nest_vertex_pairs_groups}) {
	if (@{$candidate_nest_vertex_pairs_group} == 1) {
	    push @{$matched_nest_vertex_candidate_pairs}, $candidate_nest_vertex_pairs_group->[0];
	}
	else {
	    my $same_profile_token_vertex_pairs_groups = {};

	    foreach my $candidate_nest_vertex_pair (@{$candidate_nest_vertex_pairs_group}) {
		my $nest_vertex_profile1 = $candidate_nest_vertex_profiles1->{$candidate_nest_vertex_pair->[0]};
		my $nest_vertex_profile2 = $candidate_nest_vertex_profiles2->{$candidate_nest_vertex_pair->[1]};
		my $pair_profile_token = join(@{$nest_vertex_profile1->[0]}, '_') . '-' . join(@{$nest_vertex_profile1->[1]}, '_') . '#' .
		    join(@{$nest_vertex_profile2->[0]}, '_') . '-' . join(@{$nest_vertex_profile2->[1]}, '_');
		my $same_profile_token_vertex_pairs = $same_profile_token_vertex_pairs_groups->{$pair_profile_token};
		if (!defined($same_profile_token_vertex_pairs)) {
		    $same_profile_token_vertex_pairs = [];
		    $same_profile_token_vertex_pairs_groups->{$pair_profile_token} = $same_profile_token_vertex_pairs;
		}

		push @{$same_profile_token_vertex_pairs}, $candidate_nest_vertex_pair;
	    }

	    foreach my $same_profile_token_vertex_pairs (values %{$same_profile_token_vertex_pairs_groups}) {
		my ($local_ecgm_mappings, $aligned_stem_pairs) =
		    GraphMatchTools->expand_mscs_mapping($same_profile_token_vertex_pairs, [], $aligned_stem_pairs, $stem_graph1, $stem_graph2,
							 $stem_removal_costs1, $stem_removal_costs2);
		my $local_ecgm_mapped_vertex_pairs = {};
		foreach my $local_ecgm_mapping (@{$local_ecgm_mappings}) {
		    foreach (@{$local_ecgm_mapping}) {
			$local_ecgm_mapped_vertex_pairs->{$_->[0] . '-' . $_->[1]} = $_;
		    }
		}

		foreach (values %{$local_ecgm_mapped_vertex_pairs}) {
		    push @{$matched_nest_vertex_candidate_pairs}, $_;
		}
	    }
	}
    }

    return $matched_nest_vertex_candidate_pairs, $aligned_stem_pairs;
}

1;
