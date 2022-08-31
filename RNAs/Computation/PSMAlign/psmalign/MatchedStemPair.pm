package MatchedStemPair;

use rna_stem_align;
use strict;

sub new {
    my (undef, $vertex_attrs1, $vertex_attrs2) = @_;

    my $self = {};
    $self->{stem1_upstream_base_seq} = $vertex_attrs1->{upstream_base_seq};
    $self->{stem1_downstream_base_seq} = $vertex_attrs1->{downstream_base_seq};
    $self->{stem2_upstream_base_seq} = $vertex_attrs2->{upstream_base_seq};
    $self->{stem2_downstream_base_seq} = $vertex_attrs2->{downstream_base_seq};

    $self->{stem1_pseudo_base_pairs} = $vertex_attrs1->{pseudo_base_pairs};
    $self->{stem2_pseudo_base_pairs} = $vertex_attrs2->{pseudo_base_pairs};

    my $stem1_upstream_len = $vertex_attrs1->{upstream_length};
    my $stem2_upstream_len = $vertex_attrs2->{upstream_length};

    $self->{stem1_upstream_last_pseudo_bp_pos} = $stem1_upstream_len - 1;
    $self->{stem1_downstream_first_pseudo_bp_pos} = $stem1_upstream_len;
    $self->{stem2_upstream_last_pseudo_bp_pos} = $stem2_upstream_len - 1;
    $self->{stem2_downstream_first_pseudo_bp_pos} = $stem2_upstream_len;

    $self->{stem1_upstream_len} = $stem1_upstream_len;
    $self->{stem2_upstream_len} = $stem2_upstream_len;

    $self->{stem1_upstream_end} = $vertex_attrs1->{upstream_start} + $stem1_upstream_len;
    $self->{stem1_downstream_start} = $vertex_attrs1->{downstream_start};
    $self->{stem2_upstream_end} = $vertex_attrs2->{upstream_start} + $stem2_upstream_len;
    $self->{stem2_downstream_start} = $vertex_attrs2->{downstream_start};
    _check_region_cont($self, 1);
    _check_region_cont($self, 2);

    bless $self;

    return $self;
}

sub _check_region_cont {
    my ($self, $seq_num) = @_;

    if ($seq_num != 1 && $seq_num != 2) {
	return;
    }

    my $selected_stem = 'stem' . $seq_num;

    if ($self->{$selected_stem . '_upstream_end'} >= $self->{$selected_stem . '_downstream_start'}) {
	$self->{'is_' . $selected_stem . '_region_cont'} = 1;
    }
    else {
	$self->{'is_' . $selected_stem . '_region_cont'} = 0;
    }
}

sub append_seq_to_left_end {
    my ($self, $left_seq, $seq_num, $stream_opt) = @_;

    if ($seq_num != 1 && $seq_num != 2) {
	return;
    }

    my $increased_seq_len = length($left_seq);

    $stream_opt = lc($stream_opt);
    my $selected_stem = 'stem' . $seq_num;
    my $base_seq_key = $selected_stem . '_' . $stream_opt . '_base_seq';
    my $base_seq = $self->{$base_seq_key};

    if ($stream_opt eq 'upstream') {
	$self->{$base_seq_key} = $left_seq . $base_seq;

	$self->{$selected_stem . '_upstream_last_pseudo_bp_pos'} += $increased_seq_len;
	$self->{$selected_stem . '_downstream_first_pseudo_bp_pos'} += $increased_seq_len;
	$self->{$selected_stem . '_upstream_len'} += $increased_seq_len;

	my $shifted_pseudo_base_pairs = new rna_stem_align::vector_i2();
	my $pseudo_stem_pair_key = 'stem' . $seq_num . '_pseudo_base_pairs';
	my $pseudo_base_pairs = $self->{$pseudo_stem_pair_key};
	for (my $i = 0; $i < $pseudo_base_pairs->size(); $i++) {
	    my $pseudo_stem_pair = $pseudo_base_pairs->get($i);
	    my $shifted_pseudo_stem_pair = new rna_stem_align::vector_i();
	    $shifted_pseudo_stem_pair->push($pseudo_stem_pair->get(0) + $increased_seq_len);
	    $shifted_pseudo_stem_pair->push($pseudo_stem_pair->get(1) + $increased_seq_len);
	    $shifted_pseudo_base_pairs->push($shifted_pseudo_stem_pair);
	}

	$self->{$pseudo_stem_pair_key} = $shifted_pseudo_base_pairs;
    }
    elsif ($stream_opt eq 'downstream') {
	$self->{$base_seq_key} = $left_seq . $base_seq;

	$self->{$selected_stem . '_downstream_first_pseudo_bp_pos'} += $increased_seq_len;
	$self->{$selected_stem . '_downstream_start'} -= $increased_seq_len;
	_check_region_cont($self, $seq_num);

	my $shifted_pseudo_base_pairs = new rna_stem_align::vector_i2();
	my $pseudo_stem_pair_key = $selected_stem . '_pseudo_base_pairs';
	my $pseudo_base_pairs = $self->{$pseudo_stem_pair_key};
	for (my $i = 0; $i < $pseudo_base_pairs->size(); $i++) {
	    my $pseudo_stem_pair = $pseudo_base_pairs->get($i);
	    my $shifted_pseudo_stem_pair = new rna_stem_align::vector_i();
	    $shifted_pseudo_stem_pair->push($pseudo_stem_pair->get(0));
	    $shifted_pseudo_stem_pair->push($pseudo_stem_pair->get(1) + $increased_seq_len);
	    $shifted_pseudo_base_pairs->push($shifted_pseudo_stem_pair);
	}

	$self->{$pseudo_stem_pair_key} = $shifted_pseudo_base_pairs;
    }
}

sub append_seq_to_right_end {
    my ($self, $right_seq, $seq_num, $stream_opt) = @_;

    if ($seq_num != 1 && $seq_num != 2) {
	return;
    }

    my $increased_seq_len = length($right_seq);

    $stream_opt = lc($stream_opt);
    my $selected_stem = 'stem' . $seq_num;
    my $base_seq_key = $selected_stem . '_' . $stream_opt . '_base_seq';
    my $base_seq = $self->{$base_seq_key};

    if ($stream_opt eq 'upstream') {
	$self->{$base_seq_key} = $base_seq . $right_seq;

	$self->{$selected_stem . '_downstream_first_pseudo_bp_pos'} += $increased_seq_len;
	$self->{$selected_stem . '_upstream_len'} += $increased_seq_len;
	$self->{$selected_stem . '_upstream_end'} += $increased_seq_len;
	_check_region_cont($self, $seq_num);

	my $shifted_pseudo_base_pairs = new rna_stem_align::vector_i2();
	my $pseudo_stem_pair_key = $selected_stem . '_pseudo_base_pairs';
	my $pseudo_base_pairs = $self->{$pseudo_stem_pair_key};
	for (my $i = 0; $i < $pseudo_base_pairs->size(); $i++) {
	    my $pseudo_stem_pair = $pseudo_base_pairs->get($i);
	    my $shifted_pseudo_stem_pair = new rna_stem_align::vector_i();
	    $shifted_pseudo_stem_pair->push($pseudo_stem_pair->get(0));
	    $shifted_pseudo_stem_pair->push($pseudo_stem_pair->get(1) + $increased_seq_len);
	    $shifted_pseudo_base_pairs->push($shifted_pseudo_stem_pair);
	}

	$self->{$pseudo_stem_pair_key} = $shifted_pseudo_base_pairs;
    }
    elsif ($stream_opt eq 'downstream') {
	$self->{$base_seq_key} = $base_seq . $right_seq;
    }
}

sub get_base_seq {
    my ($self, $seq_num, $stream_opt) = @_;

    return $self->{'stem' . $seq_num . '_' . $stream_opt . '_base_seq'};
}

sub get_region_align_bound {
    my ($self, $seq_num) = @_;

    my ($upstream_align_bound, $downstream_align_bound);

    if ($self->{'is_stem' . $seq_num . '_region_cont'}) {
	$upstream_align_bound = $self->{'stem' . $seq_num . '_downstream_first_pseudo_bp_pos'} - 1;
	$downstream_align_bound = $self->{'stem' . $seq_num . '_upstream_last_pseudo_bp_pos'} + 1;
    }
    else {
	$downstream_align_bound = $self->{'stem' . $seq_num . '_upstream_len'};
	$upstream_align_bound = $downstream_align_bound - 1;
    }

    return $upstream_align_bound, $downstream_align_bound;
}

sub get_pseudo_base_pairs {
    my ($self, $seq_num) = @_;

    return $self->{'stem' . $seq_num . '_pseudo_base_pairs'};
}

1;
