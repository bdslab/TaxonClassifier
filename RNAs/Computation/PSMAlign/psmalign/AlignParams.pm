package AlignParams;

use strict;

sub new {
    my $self = {};

    $self->{non_nest_motif_match_threshold} = DefaultAlignParams->NON_NEST_MOTIF_MATCH_THRESHOLD;
    $self->{init_cost_ratio_select_range_end} = DefaultAlignParams->INIT_COST_RATIO_SELECT_RANGE_END;
#    $self->{max_select_range_cost_ratio} = DefaultAlignParams->MAX_SELECT_RANGE_COST_RATIO;
    $self->{select_range_size} = DefaultAlignParams->SELECT_RANGE_SIZE;
    $self->{cost_ratio_selection_factor} = DefaultAlignParams->COST_RATIO_SELECTION_FACTOR;
    $self->{k} = DefaultAlignParams->K;
    $self->{is_no_progressive_stem_align} = 0;
    $self->{is_progressive_stem_align_only} = 0;

    bless $self;

    return $self;
}

sub set_non_nest_motif_match_threshold {
    my ($self, $non_nest_motif_match_threshold) = @_;

    $self->{non_nest_motif_match_threshold} = $non_nest_motif_match_threshold;
}

=comment
sub set_max_select_range_cost_ratio {
    my ($self, $max_select_range_cost_ratio) = @_;

    $self->{max_select_range_cost_ratio} = $max_select_range_cost_ratio;
}
=cut

sub set_k {
    my ($self, $k) = @_;

    $self->{k} = $k;
}

sub set_no_progressive_stem_align {
    my $self = shift;

    $self->{is_no_progressive_stem_align} = 1;
}

sub set_progressive_stem_align_only {
    my $self = shift;

    $self->{is_progressive_stem_align_only} = 1;
}

sub get_non_nest_motif_match_threshold {
    my $self = shift;

    return $self->{non_nest_motif_match_threshold};
}

sub get_init_cost_ratio_select_range_end {
    my $self = shift;

    return $self->{init_cost_ratio_select_range_end};
}

=comment
sub get_max_select_range_cost_ratio {
    my $self = shift;

    return $self->{max_select_range_cost_ratio};
}
=cut

sub get_select_range_size {
    my $self = shift;

    return $self->{select_range_size};
}

sub get_cost_ratio_selection_factor {
    my $self = shift;

    return $self->{cost_ratio_selection_factor};
}

sub get_k {
    my $self = shift;

    return $self->{k};
}

sub is_no_progressive_stem_align {
    my $self = shift;

    return $self->{is_no_progressive_stem_align};
}

sub is_progressive_stem_align_only {
    my $self = shift;

    return $self->{is_progressive_stem_align_only};
}

1;
