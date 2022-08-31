package MappedVertexPairCluster;

use strict;

sub cluster {
    my (undef, $mapped_vertex_pairs) = @_;

    my ($mapped_vertex_pair_vectors, $mask) = ([], []);
    my $vector_sum = 0;

    foreach (@{$mapped_vertex_pairs}) {
	push @{$mapped_vertex_pair_vectors}, [$_->[0] - $_->[1]];
	push @{$mask}, [1];
	$vector_sum += $mapped_vertex_pair_vectors->[-1][0];
    }

    my $global_centroid = $vector_sum / scalar(@{$mapped_vertex_pairs});

#    print "Global centroid: $global_centroid\n";

    my %params = (
	transpose => 0,
	method => 'm',
	dist => 'b',
	data => $mapped_vertex_pair_vectors,
	mask => $mask,
	weight => [1],
    );

    my $tree = Algorithm::Cluster::treecluster(%params);

    my ($clusters, $cluster_centroids) = ({}, {});
    my %best_clusters;
    my $max_cluster_gain = 0;

    for (my $i = 0; $i < $tree->length; $i++) {
	($clusters, $cluster_centroids) = _merge_clusters($clusters, $cluster_centroids, $tree->get($i), (-1 - $i), $mapped_vertex_pair_vectors);

	my $cluster_gain = 0;
	my @cluster_ids = keys %{$clusters};

	for (my $j = 0; $j < @cluster_ids; $j++) {
	    my $cluster_id = $cluster_ids[$j];

#	    print "Cluster centroid: " . $cluster_centroids->{$cluster_id} . ' cluster size: ' . scalar(@{$clusters->{$cluster_id}}) . "\n";

	    $cluster_gain += (scalar(@{$clusters->{$cluster_id}}) - 1) * abs($global_centroid - $cluster_centroids->{$cluster_id});
	}

#	print "$cluster_gain, $max_cluster_gain\n";

	if (($cluster_gain > $max_cluster_gain) || ($cluster_gain == 0 && $max_cluster_gain == 0)) {
	    %best_clusters = %{$clusters};
	    $max_cluster_gain = $cluster_gain;
	}
    }

    my $max_cluster_size = _get_max_cluster_size(\%best_clusters);

    my $mapped_vertex_pair_clusters = [];

    foreach (values %best_clusters) {
	my $cluster_size = @{$_};
	if ($cluster_size < $max_cluster_size / 2) {
	    next;
	}

	my $mapped_vertex_pair_cluster = [];
	foreach (@{$_}) {
	    push @{$mapped_vertex_pair_cluster}, $mapped_vertex_pairs->[$_];
	}

	push @{$mapped_vertex_pair_clusters}, $mapped_vertex_pair_cluster;
    }

    return $mapped_vertex_pair_clusters;
}

sub _merge_clusters {
    my ($clusters, $cluster_centroids, $node, $new_cluster_id, $mapped_vertex_pair_vectors) = @_;

    if ($node->left < 0 && $node->right < 0) {
	my @new_cluster = @{$clusters->{$node->left}};
	push @new_cluster, @{$clusters->{$node->right}};
	$clusters->{$new_cluster_id} = \@new_cluster;
	delete $clusters->{$node->left};
	delete $clusters->{$node->right};

	$cluster_centroids->{$new_cluster_id} = _calculate_cluster_centroid(\@new_cluster, $mapped_vertex_pair_vectors);
	delete $cluster_centroids->{$node->left};
	delete $cluster_centroids->{$node->right};
    }
    elsif ($node->left < 0) {
	my @new_cluster = @{$clusters->{$node->left}};
	push @new_cluster, $node->right;
	$clusters->{$new_cluster_id} = \@new_cluster;
	delete $clusters->{$node->left};

	$cluster_centroids->{$new_cluster_id} = _calculate_cluster_centroid(\@new_cluster, $mapped_vertex_pair_vectors);
	delete $cluster_centroids->{$node->left};
    }
    elsif ($node->right < 0) {
	my @new_cluster = @{$clusters->{$node->right}};
	push @new_cluster, $node->left;
	$clusters->{$new_cluster_id} = \@new_cluster;
	delete $clusters->{$node->right};

	$cluster_centroids->{$new_cluster_id} = _calculate_cluster_centroid(\@new_cluster, $mapped_vertex_pair_vectors);
	delete $cluster_centroids->{$node->right};
    }
    else {
	$clusters->{$new_cluster_id} = [$node->left, $node->right];

	$cluster_centroids->{$new_cluster_id} = _calculate_cluster_centroid($clusters->{$new_cluster_id}, $mapped_vertex_pair_vectors);
    }

    return $clusters, $cluster_centroids;
}

sub _calculate_cluster_centroid {
    my ($cluster, $mapped_vertex_pair_vectors) = @_;

    my $cluster_vector_sum = 0;
    foreach (@{$cluster}) {
	$cluster_vector_sum += $mapped_vertex_pair_vectors->[$_][0];
    }

    return ($cluster_vector_sum / scalar(@{$cluster}));
}

sub _get_max_cluster_size {
    my $clusters = shift;

    my $max_cluster_size = 0;

    foreach (values %{$clusters}) {
	my $cluster_size = @{$_};
	if ($cluster_size > $max_cluster_size) {
	    $max_cluster_size = $cluster_size;
	}
    }

    return $max_cluster_size;
}

=comment
sub _get_outliers {
    my $data_vectors = shift;

    my $data_vector_histogram = _get_data_vector_histogram($data_vectors);
    my @sorted_data_vector_index = sort {$a <=> $b} keys %{$data_vector_histogram};
    my $avg_vector_separation = _calculate_avg_vector_separation(\@sorted_data_vector_index);

    print "Avg: $avg_vector_separation\n";

    my $isolated_vectors = {};
    my $index_size = @sorted_data_vector_index;
    if ($index_size > 1) {
	my $curr_data_vector = $sorted_data_vector_index[0];
	my $next_data_vector = $sorted_data_vector_index[1];
	if (($next_data_vector - $curr_data_vector > $avg_vector_separation) &&
	    ($data_vector_histogram->{$curr_data_vector} <= $data_vector_histogram->{$next_data_vector})) {
	    $isolated_vectors->{$curr_data_vector} = 1;
	}

	$curr_data_vector = $sorted_data_vector_index[$index_size - 1];
	my $prev_data_vector = $sorted_data_vector_index[$index_size - 2];
	if (($curr_data_vector - $prev_data_vector > $avg_vector_separation) &&
	    ($data_vector_histogram->{$curr_data_vector} <= $data_vector_histogram->{$prev_data_vector})) {
	    $isolated_vectors->{$curr_data_vector} = 1;
	}

	for (my $i = 1; $i < $index_size - 1; $i++) {
	    $curr_data_vector = $sorted_data_vector_index[$i];
	    $prev_data_vector = $sorted_data_vector_index[$i - 1];
	    $next_data_vector = $sorted_data_vector_index[$i + 1];

	    print 'YYY: ' . ($curr_data_vector - $prev_data_vector) . ' ' . ($next_data_vector - $curr_data_vector) . "\n";

	    if (($curr_data_vector - $prev_data_vector > $avg_vector_separation) &&
		($next_data_vector - $curr_data_vector > $avg_vector_separation)) {
		print 'HHH ' . $data_vector_histogram->{$prev_data_vector} . ' ' . $data_vector_histogram->{$curr_data_vector} . ' ' .
		    $data_vector_histogram->{$next_data_vector} . "\n";
	    }

	    if (($curr_data_vector - $prev_data_vector > $avg_vector_separation) &&
		($next_data_vector - $curr_data_vector > $avg_vector_separation) &&
		($data_vector_histogram->{$curr_data_vector} <= $data_vector_histogram->{$prev_data_vector}) &&
		($data_vector_histogram->{$curr_data_vector} <= $data_vector_histogram->{$next_data_vector})) {
		$isolated_vectors->{$curr_data_vector} = 1;
	    }
	}
    }

    return $isolated_vectors;
}

sub _get_data_vector_histogram {
    my $data_vectors = shift;

    my $data_vector_histogram = {};
    foreach (@{$data_vectors}) {
	my $vector_amplitude = $_->[0];
	if (exists($data_vector_histogram->{$vector_amplitude})) {
	    $data_vector_histogram->{$vector_amplitude}++;
	}
	else {
	    $data_vector_histogram->{$vector_amplitude} = 1;
	}
    }

    return $data_vector_histogram;
}

sub _calculate_avg_vector_separation {
    my $sorted_data_vector_index = shift;

    my $total_separation = 0;
    for (my $i = 1; $i < @{$sorted_data_vector_index}; $i++) {
	$total_separation += $sorted_data_vector_index->[$i] - $sorted_data_vector_index->[$i - 1];
    }

    print "Total separation: $total_separation\n";

    return $total_separation / scalar(@{$sorted_data_vector_index});
}

sub _get_histogram_peaks {
    my ($sorted_data_vector_index, $data_vector_histogram) = @_;

    my $histogram_peak_index = {};
    my $index_size = @{$sorted_data_vector_index};

    if ($index_size > 1) {
	if ($data_vector_histogram->{$sorted_data_vector_index->[0]} > $data_vector_histogram->{$sorted_data_vector_index->[1]}) {
	    $histogram_peak_index->{0} = 1;
	}

	if ($data_vector_histogram->{$sorted_data_vector_index->[$index_size - 1]} > $data_vector_histogram->{$sorted_data_vector_index->[$index_size - 2]}) {
	    $histogram_peak_index->{$index_size - 1} = 1;
	}
    }

    foreach (my $i = 1; $i < $index_size - 1; $i++) {
	if (($data_vector_histogram->{$sorted_data_vector_index->[$i]} > $data_vector_histogram->{$sorted_data_vector_index->[$i + 1]}) &&
	    ($data_vector_histogram->{$sorted_data_vector_index->[$i]} > $data_vector_histogram->{$sorted_data_vector_index->[$i - 1]})) {
	    $histogram_peak_index->{$i} = 1;
	}
    }

    return $histogram_peak_index;
}
=cut

1;
