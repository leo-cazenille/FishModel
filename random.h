/**
 * @file Random Number Generation stuffs
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#ifndef CATS_RANDOM_H
#define CATS_RANDOM_H

#include <random>
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

namespace Fishmodel {
// Random Distributions
extern std::uniform_real_distribution<double> unif01;
extern std::uniform_real_distribution<double> unifAngle;


#ifdef ENABLE_OPENMP
template <class Engine_t = std::mt19937>
class OpenMPRNG {
	std::vector<Engine_t> _engines;
public:
	OpenMPRNG(): _engines() {
		size_t const threads = std::max(1, omp_get_max_threads());
		for(size_t seed = 0; seed < threads; ++seed) {
			_engines.push_back(Engine_t(seed));
		}
	}

	inline Engine_t& operator()() {
		auto const id = omp_get_thread_num();
		return _engines[id];
	}

	template <class Sseq>
	void seed (Sseq& q) {
		for(size_t s = 0; s < _engines.size(); ++s) {
			_engines.seed(q + s);
		}
	}
};

extern OpenMPRNG<> rne;

#else

inline std::mt19937& rne() {
	static std::random_device rnd;
	static std::mt19937 r(rnd());
	return r;
}
#endif

}

#endif
// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
