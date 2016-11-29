/**
 * @file TODO
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#include "random.h"
#include <random>

//namespace {
//std::random_device rnd;
//}

namespace CATS {
//std::mt19937 rne(rnd());
std::uniform_real_distribution<double> unif01(0.0, 1.0);
std::uniform_real_distribution<double> unifAngle(0.0, 2.0 * M_PI);

#ifdef ENABLE_OPENMP
OpenMPRNG<> rne;
#endif

}

// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
