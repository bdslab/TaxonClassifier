#!/usr/bin/perl

use lib 'rna_align/swig_version';
use rna_stem_align;
use AlignParams;
use AlignedMotifPair;
use AlignedSequencePair;
use AlignedStemPair;
use DefaultAlignParams;
use DPParser;
use GraphMatchTools;
use MappedVertexPairCluster;
use MatchedStemPair;
use NonNestMotifProbeAlign;
use SelectiveGraphMatch;
use SequenceDP;
use StemAlignment;
use StemGraph;
use StemMatchExtension;

use strict;

if (@ARGV < 2) {
    print "Usage: perl $0 <sequence 1 DP file> <sequence 2 DP file> [-o <output file>] [-t <selection threshold>] [-k <top K value>] [-d|-p]\n";
    exit;
}

my $dp_file_path1 = $ARGV[0];
my $dp_file_path2 = $ARGV[1];
my $output_file_path;

my $align_params = AlignParams->new();
($align_params, $output_file_path, my $is_param_valid) = _parse_input_params(\@ARGV, $align_params);
if (!$is_param_valid) {
    exit;
}

my ($base_seq1, $secondary_structure1, $stem_outermost_base_pairs, $stems) = DPParser->parse($dp_file_path1);
my $stem_graph1 = StemGraph->new($base_seq1, $secondary_structure1, $stem_outermost_base_pairs, $stems);
my $struct1_length = length($base_seq1);

my ($base_seq2, $secondary_structure2, $stem_outermost_base_pairs, $stems) = DPParser->parse($dp_file_path2);
my $stem_graph2 = StemGraph->new($base_seq2, $secondary_structure2, $stem_outermost_base_pairs, $stems);
my $struct2_length = length($base_seq2);

undef $stem_outermost_base_pairs;
undef $stems;

my ($graph_vertex_mappings, $aligned_stem_pairs, $aligned_motif_pairs) = StemAlignment->align($stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $align_params);
if (!defined($aligned_motif_pairs)) {
    $aligned_motif_pairs = {};
}

my $adj_stem_alignments = {};
my $best_cost_alignments = [];
my $min_cost = 999999;

foreach (@{$graph_vertex_mappings}) {
    my (undef, $merged_stem_vertices1, $merged_stem_vertices2, $aligned_merge_stem_pairs, $merge_vertex_attrs_sets1, $merge_vertex_attrs_sets2) =
	StemMatchExtension->extend_stem_match($_, $aligned_stem_pairs, $stem_graph1, $stem_graph2, $base_seq1, $base_seq2);

    my ($struct1_regions, $struct2_regions, $stem_region_map, $struct1_unmatched_stem_bp_pos, $struct2_unmatched_stem_bp_pos) =
	generate_structure_regions_and_unmatched_stem_pair_pos($_, $merged_stem_vertices1, $merged_stem_vertices2, $merge_vertex_attrs_sets1, $merge_vertex_attrs_sets2,
							       $stem_graph1, $stem_graph2, length($base_seq1), length($base_seq2));

    my $unextended_hairpin_regions = {};
    if (defined($aligned_motif_pairs)) {
	$unextended_hairpin_regions = get_unextended_hairpin_regions($struct1_regions, $struct2_regions, $merged_stem_vertices1, $merged_stem_vertices2,
								     $aligned_motif_pairs, $stem_graph1, $stem_graph2);
    }

    my $region_count = @{$struct1_regions};
    my $aligned_seq_pairs = {};

    for (my $i = 0; $i < $region_count; $i += 2) {
	if (exists($unextended_hairpin_regions->{$i})) {
	    next;
	}

	my ($struct1_region_start, $struct1_region_len) = ($struct1_regions->[$i][0], $struct1_regions->[$i][1]);
	my ($struct2_region_start, $struct2_region_len) = ($struct2_regions->[$i][0], $struct2_regions->[$i][1]);
	my $struct1_region_base_seq = substr($base_seq1, $struct1_region_start, $struct1_region_len);
	my $struct2_region_base_seq = substr($base_seq2, $struct2_region_start, $struct2_region_len);

	my $seq_align_result;
	if ($struct1_region_base_seq eq '') {
	    my $base_seq_gap = create_gap_alignment($struct2_region_len);
	    $seq_align_result = [$struct2_region_len * DefaultAlignParams->BASE_REMOVAL_COST, $base_seq_gap, $struct2_region_base_seq];
	}
	elsif ($struct2_region_base_seq eq '') {
	    my $base_seq_gap = create_gap_alignment($struct1_region_len);
	    $seq_align_result = [$struct1_region_len * DefaultAlignParams->BASE_REMOVAL_COST, $struct1_region_base_seq, $base_seq_gap];
	}
	else {
	    my ($expected_start_edge_gap_size, $expected_end_edge_gap_size) = (0, 0);

	    if ($i > 0) {
		my $prev_struct1_stem_region = $struct1_regions->[$i - 1];
		my $prev_struct2_stem_region = $struct2_regions->[$i - 1];

		if (exists($unextended_hairpin_regions->{($i - 1)})) {
		    my $prev_aligned_motif_pair = $aligned_motif_pairs->{$prev_struct1_stem_region->[3] . '-' . $prev_struct2_stem_region->[3]};
		    if ($prev_struct1_stem_region->[0] < $prev_struct1_stem_region->[2]) {
			$expected_start_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_motif($prev_aligned_motif_pair, 'right');

			if ($expected_start_edge_gap_size == 0) {
			    $expected_start_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_motif($prev_aligned_motif_pair, 'left');
			}
		    }
		    else {
			$expected_start_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_motif($prev_aligned_motif_pair, 'right');
		    }
		}
		else {
		    my $prev_aligned_stem_pair;
		    if (exists($aligned_merge_stem_pairs->{$prev_struct1_stem_region->[3] . '-' . $prev_struct2_stem_region->[3]})) {
			$prev_aligned_stem_pair = $aligned_merge_stem_pairs->{$prev_struct1_stem_region->[3] . '-' . $prev_struct2_stem_region->[3]};
		    }
		    else {
			$prev_aligned_stem_pair = $aligned_stem_pairs->{$prev_struct1_stem_region->[3] . '-' . $prev_struct2_stem_region->[3]};
		    }

		    if ($prev_struct1_stem_region->[0] < $prev_struct1_stem_region->[2]) {
			$expected_start_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_stem($prev_aligned_stem_pair, 'upstream', 'right');

			if ($expected_start_edge_gap_size == 0) {
			    $expected_start_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_stem($prev_aligned_stem_pair, 'upstream', 'left');
			}
		    }
		    else {
			$expected_start_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_stem($prev_aligned_stem_pair, 'downstream', 'right');
		    }
		}
	    }

	    if ($i < ($region_count - 1)) {
		my $next_struct1_stem_region = $struct1_regions->[$i + 1];
		my $next_struct2_stem_region = $struct2_regions->[$i + 1];

		if (exists($unextended_hairpin_regions->{($i + 1)})) {
		    my $next_aligned_motif_pair = $aligned_motif_pairs->{$next_struct1_stem_region->[3] . '-' . $next_struct2_stem_region->[3]};

		    if ($next_struct1_stem_region->[0] > $next_struct1_stem_region->[2]) {
			$expected_end_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_motif($next_aligned_motif_pair, 'left');

			if ($expected_end_edge_gap_size == 0) {
			    $expected_end_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_motif($next_aligned_motif_pair, 'right');
			}
		    }
		    else {
			$expected_end_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_motif($next_aligned_motif_pair, 'left');
		    }
		}
		else {
		    my $next_aligned_stem_pair;
		    if (exists($aligned_merge_stem_pairs->{$next_struct1_stem_region->[3] . '-' . $next_struct2_stem_region->[3]})) {
			$next_aligned_stem_pair = $aligned_merge_stem_pairs->{$next_struct1_stem_region->[3] . '-' . $next_struct2_stem_region->[3]};
		    }
		    else {
			$next_aligned_stem_pair = $aligned_stem_pairs->{$next_struct1_stem_region->[3] . '-' . $next_struct2_stem_region->[3]};
		    }

		    if ($next_struct1_stem_region->[0] > $next_struct1_stem_region->[2]) {
			$expected_end_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_stem($next_aligned_stem_pair, 'downstream', 'left');

			if ($expected_end_edge_gap_size == 0) {
			    $expected_end_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_stem($next_aligned_stem_pair, 'downstream', 'right');
			}
		    }
		    else {
			$expected_end_edge_gap_size = GraphMatchTools->get_expected_edge_gap_size_from_stem($next_aligned_stem_pair, 'upstream', 'left');
		    }
		}
	    }

	    $seq_align_result = SequenceDP->align($struct1_region_base_seq, $struct2_region_base_seq, $expected_start_edge_gap_size, $expected_end_edge_gap_size);
	}

	$aligned_seq_pairs->{$i} = AlignedSequencePair->new($seq_align_result);
    }

    my ($realigned_motif_pairs, $realigned_stem_pairs) = ({}, {});

    for (my $i = 1; $i < @{$struct1_regions}; $i += 2) {
	if ($struct1_regions->[$i][0] > $struct1_regions->[$i][2]) {
	    next;
	}

	if (exists($unextended_hairpin_regions->{$i})) {
	    my $appended_seqs = [['', '', '', ''], ['', '', '', '']];
	    my $is_bad_align = 0;

	    my $aligned_motif_pair = $aligned_motif_pairs->{$struct1_regions->[$i][3] . '-' . $struct2_regions->[$i][3]};
	    my ($prev_aligned_seq_pair, $next_aligned_seq_pair) = ($aligned_seq_pairs->{$i - 1}, $aligned_seq_pairs->{$i + 3});

	    if (($aligned_motif_pair->get_del_edge_len(1, 'left') > 0) && ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > 0)) {
		if ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_motif_pair->get_del_edge_seq(1, 'left'),
														   $prev_aligned_seq_pair->get_del_edge_seq(2, 'right'), 'right');
		    $appended_seqs->[1][0] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[1][0] = $prev_aligned_seq_pair->trim_right_edge_gap();
		}

		$is_bad_align = 1;
	    }
	    elsif (($aligned_motif_pair->get_del_edge_len(2, 'left') > 0) && ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > 0)) {
		if ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_motif_pair->get_del_edge_seq(2, 'left'),
														   $prev_aligned_seq_pair->get_del_edge_seq(1, 'right'), 'right');
		    $appended_seqs->[0][0] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[0][0] = $prev_aligned_seq_pair->trim_right_edge_gap();
		}

		$is_bad_align = 1;
	    }

	    if (($aligned_motif_pair->get_del_edge_len(1, 'right') > 0) && ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > 0)) {
		if ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_motif_pair->get_del_edge_seq(1, 'right'),
														   $next_aligned_seq_pair->get_del_edge_seq(2, 'left'), 'left');
		    $appended_seqs->[1][3] = $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[1][3] = $next_aligned_seq_pair->trim_left_edge_gap();
		}

		$is_bad_align = 1;
	    }
	    elsif (($aligned_motif_pair->get_del_edge_len(2, 'right') > 0) && ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > 0)) {
		if ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_motif_pair->get_del_edge_seq(2, 'right'),
														   $next_aligned_seq_pair->get_del_edge_seq(1, 'left'), 'left');
		    $appended_seqs->[0][3] = $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[0][3] = $next_aligned_seq_pair->trim_left_edge_gap();
		}

		$is_bad_align = 1;
	    }

	    if ($is_bad_align) {
		$realigned_motif_pairs = GraphMatchTools->align_motif_pair($struct1_regions->[$i][3], $struct2_regions->[$i][3], $realigned_motif_pairs,
									   $stem_graph1, $stem_graph2, $base_seq1, $base_seq2, $appended_seqs);
	    }

	    $i += 2;
	}
	else {
	    my ($vertex_attrs1, $vertex_attrs2);
	    if (exists($merge_vertex_attrs_sets1->{$struct1_regions->[$i][3]})) {
		$vertex_attrs1 = $merge_vertex_attrs_sets1->{$struct1_regions->[$i][3]};
	    }
	    else {
		$vertex_attrs1 = $stem_graph1->get_vertex_attrs_at($struct1_regions->[$i][3]);
	    }

	    if (exists($merge_vertex_attrs_sets2->{$struct2_regions->[$i][3]})) {
		$vertex_attrs2 = $merge_vertex_attrs_sets2->{$struct2_regions->[$i][3]};
	    }
	    else {
		$vertex_attrs2 = $stem_graph2->get_vertex_attrs_at($struct2_regions->[$i][3]);
	    }

#	my $matched_stem_pair = MatchedStemPair->new($vertex_attrs1, $vertex_attrs2);
	    my $appended_seqs = [['', '', '', ''], ['', '', '', '']];
	    my $is_bad_align = 0;

	    my $aligned_stem_pair;
	    if (exists($aligned_merge_stem_pairs->{$struct1_regions->[$i][3] . '-' . $struct2_regions->[$i][3]})) {
		$aligned_stem_pair = $aligned_merge_stem_pairs->{$struct1_regions->[$i][3] . '-' . $struct2_regions->[$i][3]};
	    }
	    else {
		$aligned_stem_pair = $aligned_stem_pairs->{$struct1_regions->[$i][3] . '-' . $struct2_regions->[$i][3]};
	    }

	    my ($prev_aligned_seq_pair, $next_aligned_seq_pair) = ($aligned_seq_pairs->{$i - 1}, $aligned_seq_pairs->{$i + 1});

=comment
	my $aligned_stem_base_seq1 = $aligned_stem_pair->get_base_seq(1, 'upstream');
	my $aligned_stem_base_seq2 = $aligned_stem_pair->get_base_seq(2, 'upstream');
	if (is_gap_dominant_alignment($aligned_stem_base_seq1) || is_gap_dominant_alignment($aligned_stem_base_seq2)) {
	    if ($prev_aligned_seq_pair->get_cost() > 0) {
		$matched_stem_pair->append_seq_to_left_end(substr($base_seq1, $struct1_regions->[$i - 1][0], $struct1_regions->[$i - 1][1]), 1, 'upstream');
		$matched_stem_pair->append_seq_to_left_end(substr($base_seq2, $struct2_regions->[$i - 1][0], $struct2_regions->[$i - 1][1]), 2, 'upstream');
		$aligned_seq_pairs->{$i - 1} = AlignedSequencePair->new([0, '', '']);
		$is_bad_align = 1;
	    }

	    if ($next_aligned_seq_pair->get_cost() > 0) {
		$matched_stem_pair->append_seq_to_right_end(substr($base_seq1, $struct1_regions->[$i + 1][0], $struct1_regions->[$i + 1][1]), 1, 'upstream');
		$matched_stem_pair->append_seq_to_right_end(substr($base_seq2, $struct2_regions->[$i + 1][0], $struct2_regions->[$i + 1][1]), 2, 'upstream');
		$aligned_seq_pairs->{$i + 1} = AlignedSequencePair->new([0, '', '']);
		$is_bad_align = 1;
	    }
	}
	else {
=cut
	    if ($aligned_stem_pair->get_del_edge_len(1, 'upstream', 'left') > 0) {
		if ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > 0) {
		    if ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(1, 'upstream', 'left'),
														       $prev_aligned_seq_pair->get_del_edge_seq(2, 'right'), 'right');
			$appended_seqs->[1][0] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		    }
		    else {
			$appended_seqs->[1][0] = $prev_aligned_seq_pair->trim_right_edge_gap();
		    }

		    $is_bad_align = 1;
		}

		if ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > 0) {
		    if ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			$appended_seqs->[1][1] = $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_stem_pair->get_del_edge_len(1, 'upstream', 'left') * 2);
		    }
		    else {
			$appended_seqs->[1][1] = $next_aligned_seq_pair->trim_left_edge_gap();
		    }

		    $is_bad_align = 1;
		}
	    }
	    elsif ($aligned_stem_pair->get_del_edge_len(2, 'upstream', 'left') > 0) {
		if ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > 0) {
		    if ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(2, 'upstream', 'left'),
														       $prev_aligned_seq_pair->get_del_edge_seq(1, 'right'), 'right');
			$appended_seqs->[0][0] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		    }
		    else {
			$appended_seqs->[0][0] = $prev_aligned_seq_pair->trim_right_edge_gap();
		    }

		    $is_bad_align = 1;
		}

		if ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > 0) {
		    if ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			$appended_seqs->[0][1] = $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_stem_pair->get_del_edge_len(2, 'upstream', 'left') * 2);
		    }
		    else {
			$appended_seqs->[0][1] = $next_aligned_seq_pair->trim_left_edge_gap();
		    }

		    $is_bad_align = 1;
		}
	    }

	    if (($aligned_stem_pair->get_del_edge_len(1, 'upstream', 'right') > 0) && ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > 0)) {
		if ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(1, 'upstream', 'right'),
														   $next_aligned_seq_pair->get_del_edge_seq(2, 'left'), 'left');
		    $appended_seqs->[1][1] = $appended_seqs->[1][1] . $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[1][1] = $appended_seqs->[1][1] . $next_aligned_seq_pair->trim_left_edge_gap();
		}

		$is_bad_align = 1;
	    }
	    elsif (($aligned_stem_pair->get_del_edge_len(2, 'upstream', 'right') > 0) && ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > 0)) {
		if ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(2, 'upstream', 'right'),
														   $next_aligned_seq_pair->get_del_edge_seq(1, 'left'), 'left');
		    $appended_seqs->[0][1] = $appended_seqs->[0][1] . $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[0][1] = $appended_seqs->[0][1] . $next_aligned_seq_pair->trim_left_edge_gap();
		}

		$is_bad_align = 1;
	    }
#	}

	    my $stem_downstream_region_index = $stem_region_map->{$i};
	    ($prev_aligned_seq_pair, $next_aligned_seq_pair) = ($aligned_seq_pairs->{$stem_downstream_region_index - 1}, $aligned_seq_pairs->{$stem_downstream_region_index + 1});

=comment
	$aligned_stem_base_seq1 = $aligned_stem_pair->get_base_seq(1, 'downstream');
	$aligned_stem_base_seq2 = $aligned_stem_pair->get_base_seq(2, 'downstream');
	if (is_gap_dominant_alignment($aligned_stem_base_seq1) || is_gap_dominant_alignment($aligned_stem_base_seq2)) {
	    if ($prev_aligned_seq_pair->get_cost() > 0) {
		$matched_stem_pair->append_seq_to_left_end(substr($base_seq1, $struct1_regions->[$stem_downstream_region_index - 1][0], $struct1_regions->[$stem_downstream_region_index - 1][1]),
							   1, 'downstream');
		$matched_stem_pair->append_seq_to_left_end(substr($base_seq2, $struct2_regions->[$stem_downstream_region_index - 1][0], $struct2_regions->[$stem_downstream_region_index - 1][1]),
							   2, 'downstream');
		$aligned_seq_pairs->{$stem_downstream_region_index - 1} = AlignedSequencePair->new([0, '', '']);
		$is_bad_align = 1;
	    }

	    if ($next_aligned_seq_pair->get_cost() > 0) {
		$matched_stem_pair->append_seq_to_right_end(substr($base_seq1, $struct1_regions->[$stem_downstream_region_index + 1][0], $struct1_regions->[$stem_downstream_region_index + 1][1]),
							    1, 'downstream');
		$matched_stem_pair->append_seq_to_right_end(substr($base_seq2, $struct2_regions->[$stem_downstream_region_index + 1][0], $struct2_regions->[$stem_downstream_region_index - 1][1]),
							    2, 'downstream');
		$aligned_seq_pairs->{$stem_downstream_region_index + 1} = AlignedSequencePair->new([0, '', '']);
		$is_bad_align = 1;
	    }
	}
	else {
=cut
	    if (($aligned_stem_pair->get_del_edge_len(1, 'downstream', 'left') > 0) && ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > 0)) {
		if ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(1, 'downstream', 'left'),
														   $prev_aligned_seq_pair->get_del_edge_seq(2, 'right'), 'right');
		    $appended_seqs->[1][2] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[1][2] = $prev_aligned_seq_pair->trim_right_edge_gap();
		}

		$is_bad_align = 1;
	    }
	    elsif (($aligned_stem_pair->get_del_edge_len(2, 'downstream', 'left') > 0) && ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > 0)) {
		if ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
		    my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(2, 'downstream', 'left'),
														   $prev_aligned_seq_pair->get_del_edge_seq(1, 'right'), 'right');
		    $appended_seqs->[0][2] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		}
		else {
		    $appended_seqs->[0][2] = $prev_aligned_seq_pair->trim_right_edge_gap();
		}

		$is_bad_align = 1;
	    }

	    if ($aligned_stem_pair->get_del_edge_len(1, 'downstream', 'right') > 0) {
		if ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > 0) {
		    if ($prev_aligned_seq_pair->get_del_edge_len(2, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			$appended_seqs->[1][2] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_stem_pair->get_del_edge_len(1, 'downstream', 'right') * 2) .
			    $appended_seqs->[1][2];
		    }
		    else {
			$appended_seqs->[1][2] = $prev_aligned_seq_pair->trim_right_edge_gap() . $appended_seqs->[1][2];
		    }

		    $is_bad_align = 1;
		}

		if ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > 0) {
		    if ($next_aligned_seq_pair->get_del_edge_len(2, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(1, 'downstream', 'right'),
														       $next_aligned_seq_pair->get_del_edge_seq(2, 'left'), 'left');
			$appended_seqs->[1][3] = $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		    }
		    else {
			$appended_seqs->[1][3] = $next_aligned_seq_pair->trim_left_edge_gap();
		    }

		    $is_bad_align = 1;
		}
	    }
	    elsif ($aligned_stem_pair->get_del_edge_len(2, 'downstream', 'right') > 0) {
		if ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > 0) {
		    if ($prev_aligned_seq_pair->get_del_edge_len(1, 'right') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			$appended_seqs->[0][2] = $prev_aligned_seq_pair->trim_right_edge_gap_by_len($aligned_stem_pair->get_del_edge_len(2, 'downstream', 'right') * 2) . 
			    $appended_seqs->[0][2];
		    }
		    else {
			$appended_seqs->[0][2] = $prev_aligned_seq_pair->trim_right_edge_gap() . $appended_seqs->[0][2];
		    }

		    $is_bad_align = 1;
		}

		if ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > 0) {
		    if ($next_aligned_seq_pair->get_del_edge_len(1, 'left') > DefaultAlignParams->MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN) {
			my $aligned_seq_edge_gap_cutoff_len = GraphMatchTools->get_eff_aligned_seq_edge_gap_cutoff_len($aligned_stem_pair->get_del_edge_seq(2, 'downstream', 'right'),
														       $next_aligned_seq_pair->get_del_edge_seq(1, 'left'), 'left');
			$appended_seqs->[0][3] = $next_aligned_seq_pair->trim_left_edge_gap_by_len($aligned_seq_edge_gap_cutoff_len);
		    }
		    else {
			$appended_seqs->[0][3] = $next_aligned_seq_pair->trim_left_edge_gap();
		    }

		    $is_bad_align = 1;
		}
	    }
#	}

	    if ($is_bad_align) {
=comment
		my $stem1_base_seq = $matched_stem_pair->get_base_seq(1, 'upstream') . $matched_stem_pair->get_base_seq(1, 'downstream');
		my $stem2_base_seq = $matched_stem_pair->get_base_seq(2, 'upstream') . $matched_stem_pair->get_base_seq(2, 'downstream');
		my ($stem1_upstream_align_bound, $stem1_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(1);
		my ($stem2_upstream_align_bound, $stem2_downstream_align_bound) = $matched_stem_pair->get_region_align_bound(2);
		my $stem_align_result = rna_stem_align::align_stem($stem1_base_seq, $stem2_base_seq, $matched_stem_pair->get_pseudo_base_pairs(1),
								   $matched_stem_pair->get_pseudo_base_pairs(2), $stem1_upstream_align_bound,
								   $stem1_downstream_align_bound, $stem2_upstream_align_bound, $stem2_downstream_align_bound);

		$realigned_stem_pairs->{$struct1_regions->[$i][3] . '-' . $struct2_regions->[$i][3]} = AlignedStemPair->new($stem_align_result);
=cut
                $realigned_stem_pairs->{$struct1_regions->[$i][3] . '-' . $struct2_regions->[$i][3]} =
		    GraphMatchTools->align_stem_pair($vertex_attrs1, $vertex_attrs2, $appended_seqs);
	    }
	}
    }

    my $overall_alignment_cost = 0;
    my ($aligned_seq1, $aligned_seq2) = ('', '');
    for (my $i = 0; $i < @{$struct1_regions}; $i += 2) {
	if ($i > 0) {
	    if (exists($unextended_hairpin_regions->{$i})) {
		my $aligned_motif_pair;
		if (exists($realigned_motif_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]})) {
		    $aligned_motif_pair = $realigned_motif_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]};
		}
		else {
		    $aligned_motif_pair = $aligned_motif_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]};
		}

		$overall_alignment_cost += $aligned_motif_pair->get_cost();
		$aligned_seq1 = $aligned_seq1 . $aligned_motif_pair->get_base_seq(1);
		$aligned_seq2 = $aligned_seq2 . $aligned_motif_pair->get_base_seq(2);

		$i += 2;
	    }
	    else {
		my $aligned_stem_pair;
		if (exists($realigned_stem_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]})) {
		    $aligned_stem_pair = $realigned_stem_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]};
		}
		elsif (exists($aligned_merge_stem_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]})) {
		    $aligned_stem_pair = $aligned_merge_stem_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]};
		}
		else {
		    $aligned_stem_pair = $aligned_stem_pairs->{$struct1_regions->[$i - 1][3] . '-' . $struct2_regions->[$i - 1][3]};
		}

		if ($struct1_regions->[$i - 1][0] < $struct1_regions->[$i - 1][2]) {
		    $overall_alignment_cost += $aligned_stem_pair->get_cost();
		    $aligned_seq1 = $aligned_seq1 . $aligned_stem_pair->get_base_seq(1, 'upstream');
		    $aligned_seq2 = $aligned_seq2 . $aligned_stem_pair->get_base_seq(2, 'upstream');
		}
		else {
		    $aligned_seq1 = $aligned_seq1 . $aligned_stem_pair->get_base_seq(1, 'downstream');
		    $aligned_seq2 = $aligned_seq2 . $aligned_stem_pair->get_base_seq(2, 'downstream');
		}
	    }
	}

	my $aligned_seq_pair = $aligned_seq_pairs->{$i};
	$overall_alignment_cost += $aligned_seq_pair->get_cost();
	$aligned_seq1 = $aligned_seq1 . $aligned_seq_pair->get_aligned_seq(1);
	$aligned_seq2 = $aligned_seq2 . $aligned_seq_pair->get_aligned_seq(2);
    }

    my $aligned_dp1 = restore_aligned_structure($aligned_seq1, $secondary_structure1);
    my $aligned_dp2 = restore_aligned_structure($aligned_seq2, $secondary_structure2);

    my $struct1_unmatched_stem_bp_count = scalar(keys %{$struct1_unmatched_stem_bp_pos}) / 2;
    my $struct2_unmatched_stem_bp_count = scalar(keys %{$struct2_unmatched_stem_bp_pos}) / 2;

    $overall_alignment_cost += ($struct1_unmatched_stem_bp_count + $struct2_unmatched_stem_bp_count) * DefaultAlignParams->BOND_BREAKING_COST;
    if ($overall_alignment_cost < $min_cost) {
	$best_cost_alignments = [[$aligned_seq1, $aligned_seq2, $aligned_dp1, $aligned_dp2, $_, $merged_stem_vertices1, $merged_stem_vertices2]];
	$min_cost = $overall_alignment_cost;
    }
    elsif ($overall_alignment_cost == $min_cost) {
	push @{$best_cost_alignments}, [$aligned_seq1, $aligned_seq2, $aligned_dp1, $aligned_dp2, $_, $merged_stem_vertices1, $merged_stem_vertices2];
    }
}

for (my $i = 0; $i < @{$best_cost_alignments}; $i++) {
    print 'Stem mapping for alignment ' . ($i + 1) . "\n";
    my ($graph_vertex_mapping, $merged_stem_vertices1, $merged_stem_vertices2) = ($best_cost_alignments->[$i][4], $best_cost_alignments->[$i][5], $best_cost_alignments->[$i][6]);

    foreach (@{$graph_vertex_mapping}) {
	my $output_str = format_stem_match_output_str($_->[0], $merged_stem_vertices1, $stem_graph1);
	$output_str = $output_str . '<=> ';
	$output_str = $output_str . format_stem_match_output_str($_->[1], $merged_stem_vertices2, $stem_graph2);
	chop $output_str;
	print $output_str . "\n";
    }

    print "\n";
}

if (defined($output_file_path)) {
    open (OUTFILE, ">$output_file_path") or die "Cannot open file at $output_file_path";

    print OUTFILE "Alignment cost: $min_cost\n";
    for (my $i = 0; $i < @{$best_cost_alignments}; $i++) {
	print 'Structural alignment ' . ($i + 1) . "\n";
	my $output_lines = format_alignment_output($best_cost_alignments->[$i]);
	foreach (@{$output_lines}) {
	    print OUTFILE $_;
	}
    }

    close OUTFILE or die "Cannot close file at $output_file_path";
}
else {
    print "Alignment cost: $min_cost\n";
    for (my $i = 0; $i < @{$best_cost_alignments}; $i++) {
	print 'Structural alignment ' . ($i + 1) . "\n";
	my $output_lines = format_alignment_output($best_cost_alignments->[$i]);
	foreach (@{$output_lines}) {
	    print $_;
	}
    }
}

sub format_alignment_output {
    my $best_cost_alignment = shift;

    my $output_lines = [];

    my $total_len = length($best_cost_alignment->[0]);
    my $curr_line_start = 0;

    while ($curr_line_start < $total_len) {
	my $curr_line_width;
	if ($curr_line_start + DefaultAlignParams->LINE_WIDTH > $total_len) {
	    $curr_line_width = $total_len - $curr_line_start;
	}
	else {
	    $curr_line_width = DefaultAlignParams->LINE_WIDTH;
	}

	for (my $i = 0; $i < 4; $i++) {
	    my $output_line;
	    if ($i == 0) {
		$output_line = 'Sequence  1: ';
	    }
	    elsif ($i == 1) {
		$output_line = 'Sequence  2: ';
	    }
	    elsif ($i == 2) {
		$output_line = 'Structure 1: ';
	    }
	    elsif ($i == 3) {
		$output_line = 'Structure 2: ';
	    }

	    $output_line = $output_line . substr($best_cost_alignment->[$i], $curr_line_start, $curr_line_width) . "\n";

	    push @{$output_lines}, $output_line;
	}

	push @{$output_lines}, "\n";

	$curr_line_start += $curr_line_width;
    }

    return $output_lines;
}

sub restore_aligned_structure {
    my ($aligned_seq, $org_secondary_structure) = @_;

    my @aligned_seq_arr = split(//, $aligned_seq);
    my @secondary_structure_arr = split(//, $org_secondary_structure);

    my $aligned_structure = [];

    my $index = 0;
    foreach (@aligned_seq_arr) {
	if ($_ eq DefaultAlignParams->BASE_SEQ_GAP) {
	    push @{$aligned_structure}, DefaultAlignParams->STRUCTURE_DOT;
	}
	else {
	    push @{$aligned_structure}, $secondary_structure_arr[$index++];
	}
    }

    return join('', @{$aligned_structure});
}

sub generate_structure_regions_and_unmatched_stem_pair_pos {
    my ($graph_vertex_mapping, $merged_stem_vertices1, $merged_stem_vertices2, $merge_vertex_attrs_sets1, $merge_vertex_attrs_sets2, $stem_graph1, $stem_graph2, $seq1_len, $seq2_len) = @_;

    my ($struct1_stem_regions, $struct1_unmatched_stem_bp_pos) = ({}, {});
    my ($struct2_stem_regions, $struct2_unmatched_stem_bp_pos) = ({}, {});
    my ($last_read_vertex1, $last_read_vertex2) = (-1, -1);

    foreach (@{$graph_vertex_mapping}) {
	if (exists($merged_stem_vertices1->{$_->[0]})) {
	    foreach my $merged_stem_vertex (@{$merged_stem_vertices1->{$_->[0]}}) {
		for (my $i = $last_read_vertex1 + 1; $i < $merged_stem_vertex; $i++) {
		    my $vertex_attrs = $stem_graph1->get_vertex_attrs_at($i);
		    $struct1_unmatched_stem_bp_pos = update_unmatched_stem_bp_pos($struct1_unmatched_stem_bp_pos, $vertex_attrs);
		}

		$last_read_vertex1 = $merged_stem_vertex;
	    }
	}
	else {
	    for (my $i = $last_read_vertex1 + 1; $i < $_->[0]; $i++) {
		my $vertex_attrs = $stem_graph1->get_vertex_attrs_at($i);
		$struct1_unmatched_stem_bp_pos = update_unmatched_stem_bp_pos($struct1_unmatched_stem_bp_pos, $vertex_attrs);
	    }

	    $last_read_vertex1 = $_->[0];
	}

	if (exists($merged_stem_vertices2->{$_->[1]})) {
	    foreach my $merged_stem_vertex (@{$merged_stem_vertices2->{$_->[1]}}) {
		for (my $i = $last_read_vertex2 + 1; $i < $merged_stem_vertex; $i++) {
		    my $vertex_attrs = $stem_graph2->get_vertex_attrs_at($i);
		    $struct2_unmatched_stem_bp_pos = update_unmatched_stem_bp_pos($struct2_unmatched_stem_bp_pos, $vertex_attrs);
		}

		$last_read_vertex2 = $merged_stem_vertex;
	    }
	}
	else {
	    for (my $i = $last_read_vertex2 + 1; $i < $_->[1]; $i++) {
		my $vertex_attrs = $stem_graph2->get_vertex_attrs_at($i);
		$struct2_unmatched_stem_bp_pos = update_unmatched_stem_bp_pos($struct2_unmatched_stem_bp_pos, $vertex_attrs);
	    }

	    $last_read_vertex2 = $_->[1];
	}

	my $vertex_attrs;
	if (exists($merge_vertex_attrs_sets1->{$_->[0]})) {
	    $vertex_attrs = $merge_vertex_attrs_sets1->{$_->[0]};
	}
	else {
	    $vertex_attrs = $stem_graph1->get_vertex_attrs_at($_->[0]);
	}

	$struct1_stem_regions->{$vertex_attrs->{upstream_start}} = [$vertex_attrs->{upstream_start}, $vertex_attrs->{upstream_length},
								    $vertex_attrs->{downstream_start}, $_->[0]];
	$struct1_stem_regions->{$vertex_attrs->{downstream_start}} = [$vertex_attrs->{downstream_start}, $vertex_attrs->{downstream_length},
								      $vertex_attrs->{upstream_start}, $_->[0]];

	if (exists($merge_vertex_attrs_sets2->{$_->[1]})) {
	    $vertex_attrs = $merge_vertex_attrs_sets2->{$_->[1]};
	}
	else {
	    $vertex_attrs = $stem_graph2->get_vertex_attrs_at($_->[1]);
	}

	$struct2_stem_regions->{$vertex_attrs->{upstream_start}} = [$vertex_attrs->{upstream_start}, $vertex_attrs->{upstream_length},
								    $vertex_attrs->{downstream_start}, $_->[1]];
	$struct2_stem_regions->{$vertex_attrs->{downstream_start}} = [$vertex_attrs->{downstream_start}, $vertex_attrs->{downstream_length},
								      $vertex_attrs->{upstream_start}, $_->[1]];
    }

    for (my $i = $last_read_vertex1 + 1; $i < $stem_graph1->get_vertex_count(); $i++) {
	my $vertex_attrs = $stem_graph1->get_vertex_attrs_at($i);
	$struct1_unmatched_stem_bp_pos = update_unmatched_stem_bp_pos($struct1_unmatched_stem_bp_pos, $vertex_attrs);
    }

    for (my $i = $last_read_vertex2 + 1; $i < $stem_graph2->get_vertex_count(); $i++) {
	my $vertex_attrs = $stem_graph2->get_vertex_attrs_at($i);
	$struct2_unmatched_stem_bp_pos = update_unmatched_stem_bp_pos($struct2_unmatched_stem_bp_pos, $vertex_attrs);
    }

    my $struct1_regions = [];
    my $stem_region_start_to_index_map = {};
    my ($last_stem_region_end, $last_stem_region_index) = (0, 1);

    foreach (sort {$a <=> $b} keys %{$struct1_stem_regions}) {
	push @{$struct1_regions}, [$last_stem_region_end, ($_ - $last_stem_region_end)];

	my $stem_region = $struct1_stem_regions->{$_};
	push @{$struct1_regions}, $stem_region;
	$last_stem_region_end = $stem_region->[0] + $stem_region->[1];

	$stem_region_start_to_index_map->{$_} = $last_stem_region_index;
	$last_stem_region_index += 2;
    }

    push @{$struct1_regions}, [$last_stem_region_end, ($seq1_len - $last_stem_region_end)];

    my $stem_region_map = {};
    foreach (keys %{$struct1_stem_regions}) {
	my $stem_region = $struct1_stem_regions->{$_};
	if ($_ > $stem_region->[2]) {
	    next;
	}

	my $stem_upstream_region_index = $stem_region_start_to_index_map->{$_};
	my $stem_downstream_region_index = $stem_region_start_to_index_map->{$stem_region->[2]};
	$stem_region_map->{$stem_upstream_region_index} = $stem_downstream_region_index;
	$stem_region_map->{$stem_downstream_region_index} = $stem_upstream_region_index;
    }

    my $struct2_regions = [];
    $last_stem_region_end = 0;

    foreach (sort {$a <=> $b} keys %{$struct2_stem_regions}) {
	push @{$struct2_regions}, [$last_stem_region_end, ($_ - $last_stem_region_end)];

	my $stem_region = $struct2_stem_regions->{$_};
	push @{$struct2_regions}, $stem_region;
	$last_stem_region_end = $stem_region->[0] + $stem_region->[1];
    }

    push @{$struct2_regions}, [$last_stem_region_end, ($seq2_len - $last_stem_region_end)];

    return $struct1_regions, $struct2_regions, $stem_region_map, $struct1_unmatched_stem_bp_pos, $struct2_unmatched_stem_bp_pos;
}

sub get_unextended_hairpin_regions {
    my ($struct1_regions, $struct2_regions, $merged_stem_vertices1, $merged_stem_vertices2, $aligned_motif_pairs, $stem_graph1, $stem_graph2) = @_;

    my $unextended_hairpin_regions = {};

    for (my $i = 1; $i < @{$struct1_regions}; $i += 2) {
	if ($struct1_regions->[$i][0] > $struct1_regions->[$i][2]) {
	    next;
	}

	my $stem_vertex1 = $struct1_regions->[$i][3];
	my $stem_vertex2 = $struct2_regions->[$i][3];
	if (($stem_graph1->get_stem_type_at($stem_vertex1) eq 'T') || ($stem_graph2->get_stem_type_at($stem_vertex2) eq 'T')) {
	    next;
	}

	if (exists($merged_stem_vertices1->{$stem_vertex1}) || exists($merged_stem_vertices2->{$stem_vertex2})) {
	    next;
	}

	if (exists($aligned_motif_pairs->{$stem_vertex1 . '-' . $stem_vertex2})) {
	    $unextended_hairpin_regions->{$i} = 1;
	    $unextended_hairpin_regions->{$i + 1} = 1;
	    $unextended_hairpin_regions->{$i + 2} = 1;
	}
    }

    return $unextended_hairpin_regions;
}

sub update_unmatched_stem_bp_pos {
    my ($unmatched_stem_bp_pos, $vertex_attrs) = @_;

    foreach (@{$vertex_attrs->{base_pairs}}) {
	$unmatched_stem_bp_pos->{$_->[0]} = $_->[1];
	$unmatched_stem_bp_pos->{$_->[1]} = $_->[0];
    }

    return $unmatched_stem_bp_pos;
}

sub create_gap_alignment {
    my $gap_length = shift;

    my $base_seq_gap = '';
    for (my $i = 0; $i < $gap_length; $i++) {
	$base_seq_gap = $base_seq_gap . DefaultAlignParams->BASE_SEQ_GAP;
    }

    return $base_seq_gap;
}

sub _parse_input_params {
    my ($input_args, $align_params) = @_;

    my $output_file_path;
    my $is_valid_param = 1;

    for (my $i = 2; $i < @{$input_args}; $i++) {
	my $param_opt = $input_args->[$i];

	if ($param_opt eq '-o') {
	    if (defined($input_args->[$i + 1])) {
		$output_file_path = $input_args->[++$i];
	    }
	    else {
		print "Output file path missing\n";
		$is_valid_param = 0;
		last;
	    }
	}
	elsif ($param_opt eq '-t') {
	    if (defined($input_args->[$i + 1])) {
		my $param_value = $input_args->[++$i];
		if ($param_value =~ /^\d+$/ && $param_value > 0) {
		    $align_params->set_non_nest_motif_match_threshold($param_value);
		}
		else {
		    print "Non-nest motif matching threshold must be positive integer\n";
		    $is_valid_param = 0;
		    last;
		}
	    }
	    else {
		print "Non-nest motif matching threshold missing\n";
		$is_valid_param = 0;
		last;
	    }
	}
=comment
	elsif ($param_opt eq '-r') {
	    if (defined($input_args->[$i + 1])) {
		my $param_value = $input_args->[++$i];
		if ($param_value =~ /^\d+(\.\d+)?$/ && $param_value >= 0) {
		    $align_params->set_max_select_range_cost_ratio($param_value);
		}
		else {
		    print "Max. selection range cost ratio must be non-negative value\n";
		    $is_valid_param = 0;
		    last;
		}
	    }
	    else {
		print "Max. selection range cost ratio missing\n";
		$is_valid_param = 0;
		last;
	    }
	}
=cut
	elsif ($param_opt eq '-k') {
	    if (defined($input_args->[$i + 1])) {
		my $param_value = $input_args->[++$i];
		if ($param_value =~ /^\d+$/ && $param_value > 0) {
		    $align_params->set_k($param_value);
		}
		else {
		    print "Top K for nesting stem match selection must be positive integer\n";
		    $is_valid_param = 0;
		    last;
		}
	    }
	    else {
		print "Top K for nesting stem match selection missing\n";
		$is_valid_param = 0;
		last;
	    }
	}
	elsif ($param_opt eq '-d') {
	    if ($align_params->is_progressive_stem_align_only()) {
		print "Progressive stem alignment cannot be disabled as it has been selected to be the only stem alignment method\n";
		$is_valid_param = 0;
	    }
	    else {
		$align_params->set_no_progressive_stem_align();
	    }
	}
	elsif ($param_opt eq '-p') {
	    if ($align_params->is_no_progressive_stem_align()) {
		print "Progressive stem alignment is not available as it has been disabled\n";
		$is_valid_param = 0;
	    }
	    else {
		$align_params->set_progressive_stem_align_only();
	    }
	}
	else {
	    print "Unknow input parameter $param_opt\n";
	    $is_valid_param = 0;
	    last;
	}
    }

    return $align_params, $output_file_path, $is_valid_param;
}

sub format_stem_match_output_str {
    my ($stem_vertex, $merged_stem_vertices, $stem_graph) = @_;

    my $output_str = '';

    if (exists($merged_stem_vertices->{$stem_vertex})) {
	foreach (@{$merged_stem_vertices->{$stem_vertex}}) {
	    $output_str = $output_str . $_ . ' ';
	    my $vertex_attrs = $stem_graph->get_vertex_attrs_at($_);
	    my $stem_base_pairs = $vertex_attrs->{base_pairs};
	    $output_str = $output_str . '(' . ($stem_base_pairs->[0][0] + 1) . ' - ' . ($stem_base_pairs->[-1][0] + 1) . ' : ' .
		($stem_base_pairs->[-1][1] + 1) . ' - ' . ($stem_base_pairs->[0][1] + 1) . ') ';
	}
    }
    else {
	$output_str = $output_str . $stem_vertex . ' ';
	my $vertex_attrs = $stem_graph->get_vertex_attrs_at($stem_vertex);
	my $stem_base_pairs = $vertex_attrs->{base_pairs};
	$output_str = $output_str . '(' . ($stem_base_pairs->[0][0] + 1) . ' - ' . ($stem_base_pairs->[-1][0] + 1) . ' : ' .
	    ($stem_base_pairs->[-1][1] + 1) . ' - ' . ($stem_base_pairs->[0][1] + 1) . ') ';
    }

    return $output_str;
}
