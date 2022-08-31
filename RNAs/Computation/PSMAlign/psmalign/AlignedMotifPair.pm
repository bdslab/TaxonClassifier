package AlignedMotifPair;

use strict;

sub new {
    my (undef, $aligned_stem_pair, $aligned_seq_pair) = @_;
#    my (undef, $motif_alignment) = @_;

    my $self = {};

#    $self->{cost} = ($motif_alignment->[0])->get(0);
#    $self->{seq1} = $motif_alignment->[1]->get(1) . ' ' . $motif_alignment->[2]->get(1);
#    $self->{seq2} = $motif_alignment->[3]->get(1) . ' ' . $motif_alignment->[4]->get(1);
    $self->{aligned_stem_pair} = $aligned_stem_pair;
    $self->{aligned_seq_pair} = $aligned_seq_pair;
#    $self->{cost} = $aligned_stem_pair->get_cost();
#    $self->{seq1} = $aligned_stem_pair->get_base_seq(1, 'upstream') . $aligned_stem_pair->get_base_seq(1, 'downstream');
#    $self->{seq2} = $aligned_stem_pair->get_base_seq(2, 'upstream') . $aligned_stem_pair->get_base_seq(2, 'downstream');

    bless $self;

    return $self;
}

=comment
sub new_large_motif {
    my (undef, $aligned_stem_pair, $aligned_seq_pair) = @_;

    my $self = {};

    $self->{aligned_stem_pair} = $aligned_stem_pair;
    $self->{aligned_seq_pair} = $aligned_seq_pair;
#    $self->{cost} = $aligned_stem_pair->get_cost() + $aligned_seq_pair->get_cost();
#    $self->{seq1} = $aligned_stem_pair->get_base_seq(1, 'upstream') . $aligned_seq_pair->get_aligned_seq(1) . $aligned_stem_pair->get_base_seq(1, 'downstream');
#    $self->{seq2} = $aligned_stem_pair->get_base_seq(2, 'upstream') . $aligned_seq_pair->get_aligned_seq(2) . $aligned_stem_pair->get_base_seq(2, 'downstream');

    bless $self;

    return $self;
}
=cut

=comment
sub new_by_cost {
    my (undef, $motif_align_cost) = @_;

    my $self = {};

    $self->{cost} = $motif_align_cost;

    bless $self;

    return $self;
}
=cut

sub get_cost {
    my $self = shift;

    if (defined($self->{aligned_seq_pair})) {
	return ($self->{aligned_stem_pair})->get_cost() + ($self->{aligned_seq_pair})->get_cost();
    }

    return ($self->{aligned_stem_pair})->get_cost();
#    return $self->{cost};
}

sub get_base_seq {
    my ($self, $seq_num) = @_;

#    my $selected_seq = 'seq' . $seq_num;

#    return $self->{$selected_seq};
    my $seq;
    if (defined($self->{aligned_seq_pair})) {
	$seq = ($self->{aligned_stem_pair})->get_base_seq($seq_num, 'upstream') . ($self->{aligned_seq_pair})->get_aligned_seq($seq_num) .
	    ($self->{aligned_stem_pair})->get_base_seq($seq_num, 'downstream');
    }
    else {
	$seq = ($self->{aligned_stem_pair})->get_base_seq($seq_num, 'upstream') . ($self->{aligned_stem_pair})->get_base_seq($seq_num, 'downstream');
    }

    return $seq;
}

sub get_del_edge_len {
    my ($self, $seq_num, $end_opt) = @_;

    if ($end_opt eq 'left') {
	return ($self->{aligned_stem_pair})->get_del_edge_len($seq_num, 'upstream', $end_opt);
    }
    elsif ($end_opt eq 'right') {
	return ($self->{aligned_stem_pair})->get_del_edge_len($seq_num, 'downstream', $end_opt);
    }

    return undef;
}

sub get_del_edge_seq {
    my ($self, $seq_num, $end_opt) = @_;

    if ($end_opt eq 'left') {
	return ($self->{aligned_stem_pair})->get_del_edge_seq($seq_num, 'upstream', $end_opt);
    }
    elsif ($end_opt eq 'right') {
	return ($self->{aligned_stem_pair})->get_del_edge_seq($seq_num, 'downstream', $end_opt);
    }

    return '';
}

1;
