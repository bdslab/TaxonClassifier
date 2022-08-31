package SequenceDP;

use strict;

#+ve edge gap size => expect sequence 1 to have gap at the corresponding edge
#-ve edge gap size => expect sequence 2 to have gap at the corresponding edge

use constant GAP_OPEN_COST => 0.00001;
use constant TOLERANCE => 1e-8;

my $edge_gap_align_pattern = DefaultAlignParams->EDGE_GAP_ALIGN_PATTERN;

sub align {
    my (undef, $seq1, $seq2, $expected_start_edge_gap_size, $expected_end_edge_gap_size) = @_;

    $seq1 = uc($seq1);
    $seq2 = uc($seq2);

    my $num_of_cols = length($seq1);
    my $num_of_rows = length($seq2);

    my @seq1_bases = split(//, $seq1);
    my @seq2_bases = split(//, $seq2);

    if ($expected_start_edge_gap_size == 0 && $expected_end_edge_gap_size == 0) {
	my ($match_table, $del_table1, $del_table2) = _initialize_dp_table_with_gap_open_cost($num_of_rows, $num_of_cols);

#	_display_dp_table($match_table, \@seq1_bases, \@seq2_bases);
#	_display_dp_table($del_table1, \@seq1_bases, \@seq2_bases);
#	_display_dp_table($del_table2, \@seq1_bases, \@seq2_bases);

	($match_table, $del_table1, $del_table2) = _fill_dp_table_with_gap_open_cost($match_table, $del_table1, $del_table2,
										     \@seq1_bases, \@seq2_bases);

#	_display_dp_table($match_table, \@seq1_bases, \@seq2_bases);
#	_display_dp_table($del_table1, \@seq1_bases, \@seq2_bases);
#	_display_dp_table($del_table2, \@seq1_bases, \@seq2_bases);

	my $min_cost = _get_min_cost_among_tables($match_table, $del_table1, $del_table2, -1, -1);
	my $table_ptr;
	if ($match_table->[-1][-1] == $min_cost) {
	    $table_ptr = 'MATCH';
	}
	elsif ($del_table1->[-1][-1] == $min_cost) {
	    $table_ptr = 'DEL1';
	}
	elsif ($del_table2->[-1][-1] == $min_cost) {
	    $table_ptr = 'DEL2';
	}

	my ($selected_align_seq_pair, $total_gap_open_cost) =
	    _back_track_with_gap_open_cost($match_table, $del_table1, $del_table2, $table_ptr, $num_of_cols, $num_of_rows,
					   \@seq1_bases, \@seq2_bases);

	my $adj_align_cost = $min_cost - $total_gap_open_cost;

	my ($aligned_seq1_start_gap_size, $aligned_seq1_end_gap_size) = _get_aligned_seq_edge_gap_sizes($selected_align_seq_pair->[0]);
	my ($aligned_seq2_start_gap_size, $aligned_seq2_end_gap_size) = _get_aligned_seq_edge_gap_sizes($selected_align_seq_pair->[1]);

	return [$adj_align_cost, $selected_align_seq_pair->[0], $selected_align_seq_pair->[1], $aligned_seq1_start_gap_size,
		$aligned_seq1_end_gap_size, $aligned_seq2_start_gap_size, $aligned_seq2_end_gap_size];
    }
    else {
	my $dp_table = _initialize_dp_table($num_of_rows, $num_of_cols);
	$dp_table = _fill_dp_table ($dp_table, \@seq1_bases, \@seq2_bases);

	my $selected_align_seq_pair;
	if ($expected_start_edge_gap_size == 0) {
	    $selected_align_seq_pair = _back_track_one_side_edge_gap($dp_table, $num_of_cols, $num_of_rows, \@seq1_bases, \@seq2_bases, 0);
	}
	elsif ($expected_end_edge_gap_size == 0) {
	    $selected_align_seq_pair = _back_track_one_side_edge_gap($dp_table, $num_of_cols, $num_of_rows, \@seq1_bases, \@seq2_bases, 1);
	}
	elsif ($num_of_cols > DefaultAlignParams->MAX_SEQ_SIZE_FOR_ENUM || $num_of_rows > DefaultAlignParams->MAX_SEQ_SIZE_FOR_ENUM) {
	    $selected_align_seq_pair = _back_track_both_side_edge_gap($dp_table, $num_of_cols, $num_of_rows, \@seq1_bases, \@seq2_bases, 0);
	}
	else {
	    my $aligned_seq_pairs = _back_track($dp_table, $num_of_cols, $num_of_rows, \@seq1_bases, \@seq2_bases);

	    my $max_seq_len = ($num_of_cols, $num_of_rows)[$num_of_cols < $num_of_rows];

	    $selected_align_seq_pair =
		_select_aligned_seq_pair($aligned_seq_pairs, $expected_start_edge_gap_size, $expected_end_edge_gap_size, $max_seq_len);
	}

	my ($aligned_seq1_start_gap_size, $aligned_seq1_end_gap_size) = _get_aligned_seq_edge_gap_sizes($selected_align_seq_pair->[0]);
	my ($aligned_seq2_start_gap_size, $aligned_seq2_end_gap_size) = _get_aligned_seq_edge_gap_sizes($selected_align_seq_pair->[1]);

	return [$dp_table->[$num_of_rows][$num_of_cols], $selected_align_seq_pair->[0], $selected_align_seq_pair->[1], $aligned_seq1_start_gap_size,
		$aligned_seq1_end_gap_size, $aligned_seq2_start_gap_size, $aligned_seq2_end_gap_size];
    }
}

sub _initialize_dp_table {
    my ($num_of_rows, $num_of_cols) = @_;

    my $dp_table = [[0]];

    for (my $i = 0; $i < $num_of_cols; $i++) {
	$dp_table->[0][$i + 1] = $dp_table->[0][$i] + DefaultAlignParams->BASE_REMOVAL_COST;
    }

    for (my $i = 0; $i < $num_of_rows; $i++) {
	push @{$dp_table}, [];
	$dp_table->[$i + 1][0] = $dp_table->[$i][0] + DefaultAlignParams->BASE_REMOVAL_COST;
    }

    return $dp_table;
}

sub _initialize_dp_table_with_gap_open_cost {
   my ($num_of_rows, $num_of_cols) = @_;

   my $max_cost = ($num_of_rows + $num_of_cols) * (DefaultAlignParams->BASE_REMOVAL_COST + GAP_OPEN_COST);

   my $match_table = [[0]];
   my $del_table1 = [[GAP_OPEN_COST]];
   my $del_table2 = [[GAP_OPEN_COST]];

   for (my $i = 1; $i <= $num_of_cols; $i++) {
       $del_table1->[0][$i] = $del_table1->[0][$i - 1] + DefaultAlignParams->BASE_REMOVAL_COST;
       $match_table->[0][$i] = $max_cost;
       $del_table2->[0][$i] = $max_cost;
   }

   for (my $i = 1; $i <= $num_of_rows; $i++) {
       push @{$del_table2}, [];
       $del_table2->[$i][0] = $del_table2->[$i - 1][0] + DefaultAlignParams->BASE_REMOVAL_COST;

       push @{$match_table}, [];
       $match_table->[$i][0] = $max_cost;

       push @{$del_table1}, [];
       $del_table1->[$i][0] = $max_cost;
   }

   return $match_table, $del_table1, $del_table2;
}

sub _fill_dp_table {
    my ($dp_table, $seq1_bases, $seq2_bases) = @_;

    for (my $i = 1; $i <= @{$seq2_bases}; $i++) {
	for (my $j = 1; $j <= @{$seq1_bases}; $j++) {
	    $dp_table->[$i][$j] = $dp_table->[$i][$j - 1] + DefaultAlignParams->BASE_REMOVAL_COST;
	    my $seq2_base_removal_cost = $dp_table->[$i - 1][$j] + DefaultAlignParams->BASE_REMOVAL_COST;
	    if ($seq2_base_removal_cost < $dp_table->[$i][$j]) {
		$dp_table->[$i][$j] = $seq2_base_removal_cost;
	    }

	    my $base_match_cost = $dp_table->[$i - 1][$j - 1] + ($seq1_bases->[$j - 1] ne $seq2_bases->[$i - 1]) * DefaultAlignParams->BASE_MISMATCH_COST;
	    if ($base_match_cost < $dp_table->[$i][$j]) {
		$dp_table->[$i][$j] = $base_match_cost;
	    }
	}
    }

    return $dp_table;
}

sub _fill_dp_table_with_gap_open_cost {
    my ($match_table, $del_table1, $del_table2, $seq1_bases, $seq2_bases) = @_;

    for (my $i = 1; $i <= @{$seq2_bases}; $i++) {
	for (my $j = 1; $j <= @{$seq1_bases}; $j++) {
#	    print "$i, $j:\n";

	    $del_table1->[$i][$j] = _get_min_base_removal_cost($del_table1->[$i][$j - 1], $match_table->[$i][$j - 1]);
	    $del_table2->[$i][$j] = _get_min_base_removal_cost($del_table2->[$i - 1][$j], $match_table->[$i - 1][$j]);

#	    print $del_table1->[$i][$j] . "\n";
#	    print $del_table2->[$i][$j] . "\n";

	    my $match_cost = ($seq1_bases->[$j - 1] ne $seq2_bases->[$i - 1]) * DefaultAlignParams->BASE_MISMATCH_COST;
	    $match_table->[$i][$j] = $match_table->[$i - 1][$j - 1] + $match_cost;

	    my $del_match_cost1 = $del_table1->[$i - 1][$j - 1] + $match_cost;
	    my $del_match_cost2 = $del_table2->[$i - 1][$j - 1] + $match_cost;
	    my $min_del_match_cost = ($del_match_cost1, $del_match_cost2)[$del_match_cost1 > $del_match_cost2];

	    if ($min_del_match_cost < $match_table->[$i][$j]) {
		$match_table->[$i][$j] = $min_del_match_cost;
	    }

# 	    print $match_table->[$i][$j] . "\n\n";

#	    _display_dp_table($match_table, $seq1_bases, $seq2_bases);
#	    _display_dp_table($del_table1, $seq1_bases, $seq2_bases);
#	    _display_dp_table($del_table2, $seq1_bases, $seq2_bases);
	}
    }

    return $match_table, $del_table1, $del_table2;
}

sub _get_min_base_removal_cost {
    my ($dp_table_base_del_cell, $dp_table_match_cell) = @_;

    my $base_removal_cost1 = $dp_table_base_del_cell + DefaultAlignParams->BASE_REMOVAL_COST;
    my $base_removal_cost2 = $dp_table_match_cell + GAP_OPEN_COST + DefaultAlignParams->BASE_REMOVAL_COST;

    return ($base_removal_cost1, $base_removal_cost2)[$base_removal_cost1 > $base_removal_cost2];
}

sub _back_track_with_gap_open_cost {
    my ($match_table, $del_table1, $del_table2, $table_ptr, $col_ptr, $row_ptr, $seq1_bases, $seq2_bases) = @_;

    my $aligned_seq_pair = ['', ''];
    my $total_gap_open_cost = 0;

#    print "$table_ptr, $row_ptr, $col_ptr\n";

    if ($col_ptr > 0 || $row_ptr > 0) {
	if ($table_ptr eq 'MATCH') {
	    my $recur_table_ptr;
	    my $base_match_cost = ($seq1_bases->[$col_ptr - 1] ne $seq2_bases->[$row_ptr - 1]) *
		DefaultAlignParams->BASE_MISMATCH_COST;
	    if (abs($match_table->[$row_ptr - 1][$col_ptr - 1] + $base_match_cost - $match_table->[$row_ptr][$col_ptr]) < TOLERANCE) {
		$recur_table_ptr = 'MATCH';
	    }
	    elsif (abs($del_table1->[$row_ptr - 1][$col_ptr - 1] + $base_match_cost - $match_table->[$row_ptr][$col_ptr]) < TOLERANCE) {
		$recur_table_ptr = 'DEL1';
	    }
	    elsif (abs($del_table2->[$row_ptr - 1][$col_ptr - 1] + $base_match_cost - $match_table->[$row_ptr][$col_ptr]) < TOLERANCE) {
		$recur_table_ptr = 'DEL2';
	    }

	    my ($recur_aligned_seq_pair, $recur_gap_open_cost) =
		_back_track_with_gap_open_cost($match_table, $del_table1, $del_table2, $recur_table_ptr, $col_ptr - 1, $row_ptr - 1,
					       $seq1_bases, $seq2_bases);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	    $total_gap_open_cost += $recur_gap_open_cost;
	}
	elsif ($table_ptr eq 'DEL1') {
	    my $recur_table_ptr;
	    if (abs($match_table->[$row_ptr][$col_ptr - 1] + DefaultAlignParams->BASE_REMOVAL_COST + GAP_OPEN_COST - $del_table1->[$row_ptr][$col_ptr])
		< TOLERANCE) {
		$recur_table_ptr = 'MATCH';
		$total_gap_open_cost += GAP_OPEN_COST;
	    }
	    elsif (abs($del_table1->[$row_ptr][$col_ptr - 1] + DefaultAlignParams->BASE_REMOVAL_COST - $del_table1->[$row_ptr][$col_ptr]) < TOLERANCE) {
		$recur_table_ptr = 'DEL1';
	    }

	    my ($recur_aligned_seq_pair, $recur_gap_open_cost) =
		_back_track_with_gap_open_cost($match_table, $del_table1, $del_table2, $recur_table_ptr, $col_ptr - 1, $row_ptr,
					       $seq1_bases, $seq2_bases);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . DefaultAlignParams->BASE_SEQ_GAP;
	    $total_gap_open_cost += $recur_gap_open_cost;
	}
	elsif ($table_ptr eq 'DEL2') {
	    my $recur_table_ptr;
	    if (abs($match_table->[$row_ptr - 1][$col_ptr] + DefaultAlignParams->BASE_REMOVAL_COST + GAP_OPEN_COST - $del_table2->[$row_ptr][$col_ptr])
		< TOLERANCE) {
		$recur_table_ptr = 'MATCH';
		$total_gap_open_cost += GAP_OPEN_COST;
	    }
	    elsif (abs($del_table2->[$row_ptr - 1][$col_ptr] + DefaultAlignParams->BASE_REMOVAL_COST - $del_table2->[$row_ptr][$col_ptr]) < TOLERANCE) {
		$recur_table_ptr = 'DEL2';
	    }

	    my ($recur_aligned_seq_pair, $recur_gap_open_cost) =
		_back_track_with_gap_open_cost($match_table, $del_table1, $del_table2, $recur_table_ptr, $col_ptr, $row_ptr - 1,
					       $seq1_bases, $seq2_bases);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . DefaultAlignParams->BASE_SEQ_GAP;
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	    $total_gap_open_cost += $recur_gap_open_cost;
	}
    }
    elsif ($table_ptr eq 'DEL1' || $table_ptr eq 'DEL2') {
	$total_gap_open_cost += GAP_OPEN_COST;
    }

    return $aligned_seq_pair, $total_gap_open_cost;
}

sub _get_min_cost_among_tables {
    my ($match_table, $del_table1, $del_table2, $col_ptr, $row_ptr) = @_;

    my $min_cost = $match_table->[$row_ptr][$col_ptr];

    if ($del_table1->[$row_ptr][$col_ptr] < $min_cost) {
	$min_cost = $del_table1->[$row_ptr][$col_ptr];
    }

    if ($del_table2->[$row_ptr][$col_ptr] < $min_cost) {
	$min_cost = $del_table2->[$row_ptr][$col_ptr];
    }

    return $min_cost;
}

sub _back_track_one_side_edge_gap {
    my ($dp_table, $col_ptr, $row_ptr, $seq1_bases, $seq2_bases, $is_prefer_gap_at_left) = @_;

    my $seq1_base_removal_cost = $dp_table->[-1][-1] + 1;
    my $seq2_base_removal_cost = $seq1_base_removal_cost;
    my $base_match_cost = $seq1_base_removal_cost;

    if ($col_ptr > 0) {
	$seq1_base_removal_cost = $dp_table->[$row_ptr][$col_ptr - 1] + DefaultAlignParams->BASE_REMOVAL_COST;

	if ($row_ptr > 0) {
	    $base_match_cost = $dp_table->[$row_ptr - 1][$col_ptr - 1] + ($seq1_bases->[$col_ptr - 1] ne $seq2_bases->[$row_ptr - 1]) * DefaultAlignParams->BASE_MISMATCH_COST;
	}
    }

    if ($row_ptr > 0) {
	$seq2_base_removal_cost = $dp_table->[$row_ptr - 1][$col_ptr] + DefaultAlignParams->BASE_REMOVAL_COST;
    }

    my $aligned_seq_pair = ['', ''];

    if ($is_prefer_gap_at_left) {
	if ($base_match_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_one_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr - 1, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
	elsif ($seq1_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_one_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . DefaultAlignParams->BASE_SEQ_GAP;
	}
	elsif ($seq2_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_one_side_edge_gap($dp_table, $col_ptr, $row_ptr - 1, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . DefaultAlignParams->BASE_SEQ_GAP;
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
    }
    else {
	if ($seq1_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_one_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr, $seq1_bases, $seq2_bases, 0);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . DefaultAlignParams->BASE_SEQ_GAP;
	}
	elsif ($seq2_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_one_side_edge_gap($dp_table, $col_ptr, $row_ptr - 1, $seq1_bases, $seq2_bases, 0);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . DefaultAlignParams->BASE_SEQ_GAP;
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
	elsif ($base_match_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_one_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr - 1, $seq1_bases, $seq2_bases, 0);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
    }

    return $aligned_seq_pair;
}

sub _back_track_both_side_edge_gap {
    my ($dp_table, $col_ptr, $row_ptr, $seq1_bases, $seq2_bases, $is_prefer_gap_at_left) = @_;

    my $seq1_base_removal_cost = $dp_table->[-1][-1] + 1;
    my $seq2_base_removal_cost = $seq1_base_removal_cost;
    my $base_match_cost = $seq1_base_removal_cost;

    if ($col_ptr > 0) {
	$seq1_base_removal_cost = $dp_table->[$row_ptr][$col_ptr - 1] + DefaultAlignParams->BASE_REMOVAL_COST;

	if ($row_ptr > 0) {
	    $base_match_cost = $dp_table->[$row_ptr - 1][$col_ptr - 1] + ($seq1_bases->[$col_ptr - 1] ne $seq2_bases->[$row_ptr - 1]) * DefaultAlignParams->BASE_MISMATCH_COST;
	}
    }

    if ($row_ptr > 0) {
	$seq2_base_removal_cost = $dp_table->[$row_ptr - 1][$col_ptr] + DefaultAlignParams->BASE_REMOVAL_COST;
    }

    my $aligned_seq_pair = ['', ''];

    if ($is_prefer_gap_at_left) {
	if ($base_match_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_both_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr - 1, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
	elsif ($seq1_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_both_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . DefaultAlignParams->BASE_SEQ_GAP;
	}
	elsif ($seq2_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_both_side_edge_gap($dp_table, $col_ptr, $row_ptr - 1, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . DefaultAlignParams->BASE_SEQ_GAP;
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
    }
    else {
	if ($seq1_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_both_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr, $seq1_bases, $seq2_bases, 0);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . DefaultAlignParams->BASE_SEQ_GAP;
	}
	elsif ($seq2_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_both_side_edge_gap($dp_table, $col_ptr, $row_ptr - 1, $seq1_bases, $seq2_bases, 0);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . DefaultAlignParams->BASE_SEQ_GAP;
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
	elsif ($base_match_cost == $dp_table->[$row_ptr][$col_ptr]) {
	    my $recur_aligned_seq_pair = _back_track_both_side_edge_gap($dp_table, $col_ptr - 1, $row_ptr - 1, $seq1_bases, $seq2_bases, 1);
	    $aligned_seq_pair->[0] = $recur_aligned_seq_pair->[0] . $seq1_bases->[$col_ptr - 1];
	    $aligned_seq_pair->[1] = $recur_aligned_seq_pair->[1] . $seq2_bases->[$row_ptr - 1];
	}
    }

    return $aligned_seq_pair;
}

sub _back_track {
    my ($dp_table, $col_ptr, $row_ptr, $seq1_bases, $seq2_bases) = @_;

    my $aligned_seq_pairs = [];

    my $seq1_base_removal_cost = $dp_table->[-1][-1] + 1;
    my $seq2_base_removal_cost = $seq1_base_removal_cost;
    my $base_match_cost = $seq1_base_removal_cost;

    if ($col_ptr > 0) {
	$seq1_base_removal_cost = $dp_table->[$row_ptr][$col_ptr - 1] + DefaultAlignParams->BASE_REMOVAL_COST;

	if ($row_ptr > 0) {
	    $base_match_cost = $dp_table->[$row_ptr - 1][$col_ptr - 1] + ($seq1_bases->[$col_ptr - 1] ne $seq2_bases->[$row_ptr - 1]) * DefaultAlignParams->BASE_MISMATCH_COST;
	}
    }

    if ($row_ptr > 0) {
	$seq2_base_removal_cost = $dp_table->[$row_ptr - 1][$col_ptr] + DefaultAlignParams->BASE_REMOVAL_COST;
    }

    if ($seq1_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	my $recur_aligned_seq_pairs = _back_track($dp_table, $col_ptr - 1, $row_ptr, $seq1_bases, $seq2_bases);
	$aligned_seq_pairs = _extend_aligned_seq_pairs($aligned_seq_pairs, $recur_aligned_seq_pairs, $seq1_bases->[$col_ptr - 1], DefaultAlignParams->BASE_SEQ_GAP);
    }

    if ($seq2_base_removal_cost == $dp_table->[$row_ptr][$col_ptr]) {
	my $recur_aligned_seq_pairs = _back_track($dp_table, $col_ptr, $row_ptr - 1, $seq1_bases, $seq2_bases);
	$aligned_seq_pairs = _extend_aligned_seq_pairs($aligned_seq_pairs, $recur_aligned_seq_pairs, DefaultAlignParams->BASE_SEQ_GAP, $seq2_bases->[$row_ptr - 1]);
    }

    if ($base_match_cost == $dp_table->[$row_ptr][$col_ptr]) {
	my $recur_aligned_seq_pairs = _back_track($dp_table, $col_ptr - 1, $row_ptr - 1, $seq1_bases, $seq2_bases);
	$aligned_seq_pairs = _extend_aligned_seq_pairs($aligned_seq_pairs, $recur_aligned_seq_pairs, $seq1_bases->[$col_ptr - 1], $seq2_bases->[$row_ptr - 1]);
    }

    return $aligned_seq_pairs;
}

sub _extend_aligned_seq_pairs {
    my ($aligned_seq_pairs, $recur_aligned_seq_pairs, $aligned_seq1_base, $aligned_seq2_base) = @_;

    if (defined($recur_aligned_seq_pairs->[0])) {
	foreach (@{$recur_aligned_seq_pairs}) {
	    my @aligned_seq_pair = @{$_};
	    $aligned_seq_pair[0] = $aligned_seq_pair[0] . $aligned_seq1_base;
	    $aligned_seq_pair[1] = $aligned_seq_pair[1] . $aligned_seq2_base;
	    push @{$aligned_seq_pairs}, \@aligned_seq_pair;
	}
    }
    else {
	push @{$aligned_seq_pairs}, [$aligned_seq1_base, $aligned_seq2_base];
    }

    return $aligned_seq_pairs;
}

sub _select_aligned_seq_pair {
    my ($aligned_seq_pairs, $expected_start_edge_gap_size, $expected_end_edge_gap_size, $max_seq_len) = @_;

    my $selected_seq_pair;

    if (!defined($expected_start_edge_gap_size)) {
	$expected_start_edge_gap_size = 0;
    }

    if (!defined($expected_end_edge_gap_size)) {
	$expected_end_edge_gap_size = 0;
    }

    my ($eff_start_edge_gap_size, $eff_end_edge_gap_size) = (0, 0);
    if ($expected_start_edge_gap_size > 0) {
	if ($expected_end_edge_gap_size > 0) {
	    $eff_start_edge_gap_size = int($expected_start_edge_gap_size / ($expected_start_edge_gap_size + $expected_end_edge_gap_size) * $max_seq_len + 0.5);
	    $eff_end_edge_gap_size = $max_seq_len - $eff_start_edge_gap_size;
	}
	elsif ($expected_end_edge_gap_size < 0) {
	    $eff_start_edge_gap_size = $max_seq_len;
	    $eff_end_edge_gap_size = $max_seq_len * -1;
	}
	else {
	    $eff_start_edge_gap_size = $max_seq_len;
	}
    }
    elsif ($expected_start_edge_gap_size < 0) {
	if ($expected_end_edge_gap_size < 0) {
	    $eff_start_edge_gap_size = int($expected_start_edge_gap_size / ($expected_start_edge_gap_size + $expected_end_edge_gap_size) * $max_seq_len + 0.5) * -1;
	    $eff_end_edge_gap_size = ($max_seq_len + $eff_start_edge_gap_size) * -1;
	}
	elsif ($expected_end_edge_gap_size > 0) {
	    $eff_start_edge_gap_size = $max_seq_len * -1;
	    $eff_end_edge_gap_size = $max_seq_len;
	}
	else {
	    $eff_start_edge_gap_size = $max_seq_len * -1;
	}
    }
    else {
	if ($expected_end_edge_gap_size > 0) {
	    $eff_end_edge_gap_size = $max_seq_len;
	}
	elsif ($expected_end_edge_gap_size < 0) {
	    $eff_end_edge_gap_size = $max_seq_len * -1;
	}
    }

    my $best_total_cost = _calculate_euclidean_sqr_dist($eff_start_edge_gap_size, $eff_end_edge_gap_size);

    foreach (@{$aligned_seq_pairs}) {
	my $candidate_total_cost = _get_total_edge_gap_cost($_->[0], $_->[1], $eff_start_edge_gap_size, $eff_end_edge_gap_size);
	if ($candidate_total_cost < $best_total_cost || !defined($selected_seq_pair)) {
	    $selected_seq_pair = $_;
	    $best_total_cost = $candidate_total_cost;
	}
    }

    return $selected_seq_pair;
}

sub _get_total_edge_gap_cost {
    my ($aligned_seq1, $aligned_seq2, $eff_start_edge_gap_size, $eff_end_edge_gap_size) = @_;

    my ($seq1_start_edge_gap_size, $seq1_end_edge_gap_size) = _get_aligned_seq_edge_gap_sizes($aligned_seq1);
    my ($seq2_start_edge_gap_size, $seq2_end_edge_gap_size) = _get_aligned_seq_edge_gap_sizes($aligned_seq2);

    my $start_edge_gap_cost = 0;
    if ($eff_start_edge_gap_size > 0) {
	$start_edge_gap_cost = ($eff_start_edge_gap_size - $seq1_start_edge_gap_size);
	if ($start_edge_gap_cost < 0) {
	    $start_edge_gap_cost = 0;
	}
    }
    elsif ($eff_start_edge_gap_size < 0) {
	$start_edge_gap_cost = ($eff_start_edge_gap_size + $seq2_start_edge_gap_size); 
	if ($start_edge_gap_cost > 0) {
	    $start_edge_gap_cost = 0;
	}
    }

    my $end_edge_gap_cost = 0;
    if ($eff_end_edge_gap_size > 0) {
	$end_edge_gap_cost = ($eff_end_edge_gap_size - $seq1_end_edge_gap_size);
	if ($end_edge_gap_cost < 0) {
	    $end_edge_gap_cost = 0;
	}
    }
    elsif ($eff_end_edge_gap_size < 0) {
	$end_edge_gap_cost = ($eff_end_edge_gap_size + $seq2_end_edge_gap_size);
	if ($end_edge_gap_cost > 0) {
	    $end_edge_gap_cost = 0;
	}
    }

    return _calculate_euclidean_sqr_dist($start_edge_gap_cost, $end_edge_gap_cost);
}

sub _calculate_euclidean_sqr_dist {
    my ($x_dist, $y_dist) = @_;

    return $x_dist ** 2 + $y_dist ** 2;
}

sub _get_aligned_seq_edge_gap_sizes {
    my $aligned_seq = shift;

    $aligned_seq =~ /$edge_gap_align_pattern/;

    return length($1), length($2);
}

sub _display_dp_table {
    my ($dp_table, $seq1_bases, $seq2_bases) = @_;

    print "\t_\t";
    foreach (@{$seq1_bases}) {
	print "$_\t";
    }
    print "\n";

    for (my $i = 0; $i <= @{$seq2_bases}; $i++) {
	if ($i > 0) {
	    print $seq2_bases->[$i - 1] . "\t";
	}
	else {
	    print "_\t";
	}

	foreach (@{$dp_table->[$i]}) {
	    print "$_\t";
	}

	print "\n";
    }
    print "\n";
}

1;
