package DefaultAlignParams;

use constant BASE_SEQ_GAP => '-';
use constant STRUCTURE_DOT => '.';

use constant EDGE_GAP_ALIGN_PATTERN => '^([\\' . BASE_SEQ_GAP . '-]*)[A-Za-z\-]+?([\\' . BASE_SEQ_GAP . ']*)$';
use constant MAX_NON_NEST_MOTIF_SIZE_FOR_EXACT_ALIGN => 300;
use constant MAX_DEL_EDGE_SEQ_LEN_FOR_EXACT_ALIGN => 100;
use constant MAX_SEQ_SIZE_FOR_ENUM => 30;

use constant BOND_BREAKING_COST => 1.5;
use constant BOND_BREAKING_AND_BASE_REMOVAL_COST => 1.75;
use constant PAIR_REMOVAL_COST => 2;
use constant BASE_REMOVAL_COST => 1;
use constant BASE_MISMATCH_COST => 1;
#use constant ARC_BREAKING_ONLY_COST => 1.5;
#use constant ARC_BREAKING_WITH_BASE_REMOVAL_COST => 1.75;
#use constant ARC_BREAKING_WITH_BASE_MISMATCH_AND_REMOVAL_COST => 2.75;
use constant PAIR_BASES_MISMATCH_COST => 1;

use constant NON_NEST_MOTIF_MATCH_THRESHOLD => 5;
use constant INIT_COST_RATIO_SELECT_RANGE_END => 0.1;
#use constant MAX_SELECT_RANGE_COST_RATIO => 0.5;
use constant SELECT_RANGE_SIZE => 0.05;
use constant COST_RATIO_SELECTION_FACTOR => 2;

use constant PARTITION_STEM_COUNT_THRESHOLD => 7;
use constant MAX_MIXED_STEM_TYPE_CANDIDATES => 2;

use constant K => 1;

use constant LINE_WIDTH => 150;

1;
