package AlignedSequencePair;

use strict;

use constant BASE_SEQ_GAP => '-';
#use constant STRUCTURE_DOT => '.';

sub new {
    my (undef, $seq_alignment) = @_;

    my $self = {};

    if (ref($seq_alignment->[0]) eq '') {
	$self->{cost} = $seq_alignment->[0];
	$self->{seq1} = $seq_alignment->[1];
	$self->{seq2} = $seq_alignment->[2];

	my $seq1_len = length($seq_alignment->[1]);
	my $seq2_len = length($seq_alignment->[2]);
=comment
	$self->{seq1_dp} = _create_unpaired_dp($seq1_len);
	$self->{seq2_dp} = _create_unpaired_dp($seq2_len);
=cut

	if (defined($seq_alignment->[6])) {
	    $self->{seq1_left_del_edge_len} = $seq_alignment->[5];
	    $self->{seq1_right_del_edge_len} = $seq_alignment->[6];
	    $self->{seq2_left_del_edge_len} = $seq_alignment->[3];
	    $self->{seq2_right_del_edge_len} = $seq_alignment->[4];
	}
	else {
	    if (substr($self->{seq1}, 0, 1) ne BASE_SEQ_GAP) {
		$self->{seq1_left_del_edge_len} = $seq1_len;
		$self->{seq1_right_del_edge_len} = $self->{seq1_left_del_edge_len};
		$self->{seq2_left_del_edge_len} = 0;
		$self->{seq2_right_del_edge_len} = 0;
	    }
	    elsif (substr($self->{seq2}, 0, 1) ne BASE_SEQ_GAP) {
		$self->{seq1_left_del_edge_len} = 0;
		$self->{seq1_right_del_edge_len} = 0;
		$self->{seq2_left_del_edge_len} = $seq2_len;
		$self->{seq2_right_del_edge_len} = $self->{seq2_left_del_edge_len};
	    }
	}
    }
    else {
	$self->{cost} = ($seq_alignment->[0])->get(0);

	my ($seq1_alignment, $seq2_alignment) = ($seq_alignment->[1], $seq_alignment->[3]);

	$self->{seq1} = $seq1_alignment->get(1);
	$self->{seq2} = $seq2_alignment->get(1);
=comment
	$self->{seq1_dp} = $seq1_alignment->get(0);
	$self->{seq2_dp} = $seq2_alignment->get(0);
=cut
	$self->{seq1_left_del_edge_len} = $seq1_alignment->get(2);
	$self->{seq1_right_del_edge_len} = $seq1_alignment->get(3);
	$self->{seq2_left_del_edge_len} = $seq2_alignment->get(2);
	$self->{seq2_right_del_edge_len} = $seq2_alignment->get(3);
    }

    bless $self;

    return $self;
}

=comment
sub _create_unpaired_dp {
    my $length = shift;

    my $struct_dp = '';
    for (my $i = 0; $i < $length; $i++) {
	$struct_dp = $struct_dp . STRUCTURE_DOT;
    }

    return $struct_dp;
}

sub set_seq_dp {
    my ($self, $seq_dp, $seq_num) = @_;

    my $aligned_seq_dp = '';
    my @seq_dp_arr = split(//, $seq_dp);
    my $next_append_seq_dp_pos = 0;

    foreach (split(//, $self->{'seq' . $seq_num})) {
	if ($_ eq BASE_SEQ_GAP) {
	    $aligned_seq_dp = $aligned_seq_dp . STRUCTURE_DOT;
	}
	else {
	    $aligned_seq_dp = $aligned_seq_dp . $seq_dp_arr[$next_append_seq_dp_pos++]; 
	}
    }

    $self->{'seq' . $seq_num . '_dp'} = $aligned_seq_dp;
}
=cut

sub clone {
    my (undef, $input_seq_pair) = @_;

    my $self = {};

    $self->{cost} = $input_seq_pair->get_cost();
    $self->{seq1} = $input_seq_pair->get_aligned_seq(1);
    $self->{seq2} = $input_seq_pair->get_aligned_seq(2);
=comment
    $self->{seq1_dp} = $input_seq_pair->get_seq_dp(1);
    $self->{seq2_dp} = $input_seq_pair->get_seq_dp(2);
=cut
    $self->{seq1_left_del_edge_len} = $input_seq_pair->get_del_edge_len(1, 'left');
    $self->{seq1_right_del_edge_len} = $input_seq_pair->get_del_edge_len(1, 'right');
    $self->{seq2_left_del_edge_len} = $input_seq_pair->get_del_edge_len(2, 'left');
    $self->{seq2_right_del_edge_len} = $input_seq_pair->get_del_edge_len(2, 'right');

    bless $self;

    return $self;
}

sub set_and_swap {
    my ($self, $input_seq_pair) = @_;

    $self->{cost} = $input_seq_pair->get_cost();
    $self->{seq1} = $input_seq_pair->get_aligned_seq(2);
    $self->{seq2} = $input_seq_pair->get_aligned_seq(1);
=comment
    $self->{seq1_dp} = $input_seq_pair->get_seq_dp(2);
    $self->{seq2_dp} = $input_seq_pair->get_seq_dp(1);
=cut
    $self->{seq1_left_del_edge_len} = $input_seq_pair->get_del_edge_len(2, 'left');
    $self->{seq1_right_del_edge_len} = $input_seq_pair->get_del_edge_len(2, 'right');
    $self->{seq2_left_del_edge_len} = $input_seq_pair->get_del_edge_len(1, 'left');
    $self->{seq2_right_del_edge_len} = $input_seq_pair->get_del_edge_len(1, 'right');
}

sub trim_left_edge_gap {
    my $self = shift;

    my $trimmed_seq = '';
    my $trim_len = 0;

    if ($self->{seq1_left_del_edge_len} > 0)  {
	$trim_len = $self->{seq1_left_del_edge_len};
	$trimmed_seq = substr($self->{seq1}, 0, $trim_len);
	$self->{seq1_left_del_edge_len} = 0;
    }
    elsif ($self->{seq2_left_del_edge_len} > 0) {
	$trim_len = $self->{seq2_left_del_edge_len};
	$trimmed_seq = substr($self->{seq2}, 0, $trim_len);
	$self->{seq2_left_del_edge_len} = 0;
    }

    $self->{seq1} = substr($self->{seq1}, $trim_len);
    $self->{seq2} = substr($self->{seq2}, $trim_len);
=comment
    $self->{seq1_dp} = substr($self->{seq1_dp}, $trim_len);
    $self->{seq2_dp} = substr($self->{seq2_dp}, $trim_len);
=cut
    $self->{cost} -= $trim_len * DefaultAlignParams->BASE_REMOVAL_COST;

    if (length($self->{seq1}) == 0) {
	$self->{seq1_right_del_edge_len} = 0;
    }

    if (length($self->{seq2}) == 0) {
	$self->{seq2_right_del_edge_len} = 0;
    }

    return $trimmed_seq;
}

sub trim_left_edge_gap_by_len {
    my ($self, $expected_trim_len) = @_;

    my $trimmed_seq = '';
    my $trim_len = 0;

    if ($self->{seq1_left_del_edge_len} > 0)  {
	$trim_len = $self->{seq1_left_del_edge_len};
	if ($trim_len > $expected_trim_len) {
	    $trim_len = $expected_trim_len;
	}

	$trimmed_seq = substr($self->{seq1}, 0, $trim_len);
	$self->{seq1_left_del_edge_len} -= $trim_len;
    }
    elsif ($self->{seq2_left_del_edge_len} > 0) {
	$trim_len = $self->{seq2_left_del_edge_len};
	if ($trim_len > $expected_trim_len) {
	    $trim_len = $expected_trim_len;
	}

	$trimmed_seq = substr($self->{seq2}, 0, $trim_len);
	$self->{seq2_left_del_edge_len} -= $trim_len;
    }

    $self->{seq1} = substr($self->{seq1}, $trim_len);
    $self->{seq2} = substr($self->{seq2}, $trim_len);
    $self->{cost} -= $trim_len * DefaultAlignParams->BASE_REMOVAL_COST;

    if (length($self->{seq1}) == 0) {
	$self->{seq1_right_del_edge_len} = 0;
    }

    if (length($self->{seq2}) == 0) {
	$self->{seq2_right_del_edge_len} = 0;
    }

    return $trimmed_seq;
}

sub trim_right_edge_gap {
    my $self = shift;

    my $trimmed_seq = '';
    my $trim_len = 0;

    if ($self->{seq1_right_del_edge_len} > 0)  {
	$trim_len = -1 * $self->{seq1_right_del_edge_len};
	$trimmed_seq = substr($self->{seq1}, $trim_len);
	$self->{seq1_right_del_edge_len} = 0;
    }
    else {
	$trim_len = -1 * $self->{seq2_right_del_edge_len};
	$trimmed_seq = substr($self->{seq2}, $trim_len);
	$self->{seq2_right_del_edge_len} = 0;
    }

    $self->{seq1} = substr($self->{seq1}, 0, $trim_len);
    $self->{seq2} = substr($self->{seq2}, 0, $trim_len);
=comment
    $self->{seq1_dp} = substr($self->{seq1_dp}, 0, $trim_len);
    $self->{seq2_dp} = substr($self->{seq2_dp}, 0, $trim_len);
=cut
    $self->{cost} += $trim_len * DefaultAlignParams->BASE_REMOVAL_COST;

    if (length($self->{seq1}) == 0) {
	$self->{seq1_left_del_edge_len} = 0;
    }

    if (length($self->{seq2}) == 0) {
	$self->{seq2_left_del_edge_len} = 0;
    }

    return $trimmed_seq;
}

sub trim_right_edge_gap_by_len {
    my ($self, $expected_trim_len) = @_;

    my $trimmed_seq = '';
    my $trim_len = 0;
    $expected_trim_len *= -1;

    if ($self->{seq1_right_del_edge_len} > 0)  {
	$trim_len = -1 * $self->{seq1_right_del_edge_len};
	if ($trim_len < $expected_trim_len) {
	    $trim_len = $expected_trim_len;
	}

	$trimmed_seq = substr($self->{seq1}, $trim_len);
	$self->{seq1_right_del_edge_len} += $trim_len;
    }
    else {
	$trim_len = -1 * $self->{seq2_right_del_edge_len};
	if ($trim_len < $expected_trim_len) {
	    $trim_len = $expected_trim_len;
	}

	$trimmed_seq = substr($self->{seq2}, $trim_len);
	$self->{seq2_right_del_edge_len} += $trim_len;
    }

    $self->{seq1} = substr($self->{seq1}, 0, $trim_len);
    $self->{seq2} = substr($self->{seq2}, 0, $trim_len);
    $self->{cost} += $trim_len * DefaultAlignParams->BASE_REMOVAL_COST;

    if (length($self->{seq1}) == 0) {
	$self->{seq1_left_del_edge_len} = 0;
    }

    if (length($self->{seq2}) == 0) {
	$self->{seq2_left_del_edge_len} = 0;
    }

    return $trimmed_seq;
}

sub get_cost {
    my $self = shift;

    return $self->{cost};
}

sub get_aligned_seq {
    my ($self, $seq_num) = @_;

    return $self->{'seq' . $seq_num};
}

=comment
sub get_seq_dp {
    my ($self, $seq_num) = @_;

    return $self->{'seq' . $seq_num . '_dp'};
}
=cut

sub get_del_edge_len {
    my ($self, $seq_num, $end_opt) = @_;

    return $self->{'seq' . $seq_num . '_' . $end_opt . '_del_edge_len'};
}

sub get_del_edge_seq {
   my ($self, $seq_num, $end_opt) = @_;

    my $del_edge_len = get_del_edge_len($self, $seq_num, $end_opt);
    if (!defined($del_edge_len) || $del_edge_len == 0) {
	return '';
    }

    my $aligned_seq = get_aligned_seq($self, $seq_num);
    if ($end_opt eq 'left') {
	return substr($aligned_seq, 0, $del_edge_len);
    }
    elsif ($end_opt eq 'right') {
	return substr($aligned_seq, (-1 * $del_edge_len));
    }

    return '';
}

sub get_total_del_edge_len {
    my $self = shift;

    return ($self->{seq1_left_del_edge_len} + $self->{seq1_right_del_edge_len} + $self->{seq2_left_del_edge_len} +
	    $self->{seq2_right_del_edge_len});
}

1;
