#ifndef _QP_CONSTS_HPP_
#define _QP_CONSTS_HPP_

#include <limits>

namespace QP
{
    const double inf = 1e-3*std::numeric_limits<double>::max();
    const double eps = std::numeric_limits<double>::epsilon();
}

#endif