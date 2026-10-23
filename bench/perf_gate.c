#include "perf_gate.h"

#include <errno.h>

#ifdef __linux__
#include <sys/prctl.h>
#endif

bool perf_gate(bool enable)
{
#ifdef __linux__
    return prctl(enable ? PR_TASK_PERF_EVENTS_ENABLE : PR_TASK_PERF_EVENTS_DISABLE,
                 0L, 0L, 0L, 0L) == 0;
#else
    (void)enable;
    errno = ENOTSUP;
    return false;
#endif
}
