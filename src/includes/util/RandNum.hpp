#ifndef GRAPHWALKER_RANDNUM_HPP
#define GRAPHWALKER_RANDNUM_HPP

#include <stdint.h>

typedef long long LL;
typedef unsigned int uint;

class RandNum {
public:
	/**
	 * http://xoroshiro.di.unimi.it/#shootout
	 * We use xoroshiro128+, the fastest generator available
	 */

	uint64_t rng_seed0, rng_seed1;

	RandNum(uint64_t seed) {
		for (int i = 0; i < 2; i++) {
			LL z = seed += UINT64_C(0x9E3779B97F4A7C15);
			z = (z ^ z >> 30) * UINT64_C(0xBF58476D1CE4E5B9);
			z = (z ^ z >> 27) * UINT64_C(0x94D049BB133111EB);
			if (i == 0)
				rng_seed0 = z ^ (z >> 31);
			else
				rng_seed1 = z ^ (z >> 31);
		}
	}

	void reInit(uint64_t seed) {
		for (int i = 0; i < 2; i++) {
			LL z = seed += UINT64_C(0x9E3779B97F4A7C15);
			z = (z ^ z >> 30) * UINT64_C(0xBF58476D1CE4E5B9);
			z = (z ^ z >> 27) * UINT64_C(0x94D049BB133111EB);
			if (i == 0)
				rng_seed0 = z ^ (z >> 31);
			else
				rng_seed1 = z ^ (z >> 31);
		}
	}

	static inline uint64_t rotl(const uint64_t x, int k) {
		return (x << k) | (x >> (64 - k));
	}

	uint64_t lRand() {
		const uint64_t s0 = rng_seed0;
		uint64_t s1 = rng_seed1;
		const uint64_t result = s0 + s1;
		s1 ^= s0;
		rng_seed0 = rotl(s0, 55) ^ s1 ^ (s1 << 14);
		rng_seed1 = rotl(s1, 36);
		return result;
	}

	/* generate a double from [0, 1) */
	double dRand() {
		const union un {
			uint64_t i;
			double d;
		} a = {UINT64_C(0x3FF) << 52 | lRand() >> 12};
		return a.d - 1.0;
	}

	/* generate an int from [0, max) */
	int iRand(int max) { return lRand() % max; }

	uint32_t iRand(uint32_t max) { return lRand() % max; }

	/* generate an int from [int, max) */
	int iRand(int min, int max) { return lRand() % (max - min) + min; }
};

#endif //GRAPHWALKER_RANDNUM_HPP