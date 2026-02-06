#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
#if defined(AUTOPAS_USE_OPENMP)
    int _autopas_prefered_num_threads = autopas_get_max_threads();
#endif
}
