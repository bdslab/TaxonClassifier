package BracketPairs;
use strict;

my $open_bracket_map = {")" => "(", "]" => "[", "}" => "{", ">" => "<"};

sub is_open_bracket {
    my (undef, $symbol) = @_;

    if ($symbol =~ /^[\(\[{<A-Z]$/) {
	return 1;
    }

    return 0;
}

sub get_open_bracket {
    my (undef, $close_bracket) = @_;

    if ($close_bracket =~ /^[\)\]}>]$/) {
	return $open_bracket_map->{$close_bracket};
    }
    elsif ($close_bracket =~ /^[a-z]$/) {
	return uc $close_bracket;
    }
    else {
	die "Unknown closing bracket\n";
    }
}

1;
