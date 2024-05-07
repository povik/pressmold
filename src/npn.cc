#include <algorithm>
#include <functional>
#include <bit>
#include <set>

#include <assert.h>

#include "npn.h"

truth6 cofactor_masks[6] = {
	0xaaaaaaaaaaaaaaaa,
	0xcccccccccccccccc,
	0xf0f0f0f0f0f0f0f0,
	0xff00ff00ff00ff00,
	0xffff0000ffff0000,
	0xffffffff00000000
};

truth6 npn_semiclass(truth6 m, int ninputs, NPN &npn)
{
	npn = NPN{};

	if (!ninputs)
		return 0;

	int nbits = 1 << ninputs;
	truth6 mask = ((truth6) 2 << (nbits - 1)) - 1; // damn you, C++!
	m &= mask;
	if (std::popcount(m) > nbits / 2) {
		npn.oc = true;
		m ^= mask;
	}

	bool compls[6] = {};
	int popcount[6];
	int order[6];

	for (int i = 0; i < ninputs; i++) {
		int nfactor = std::popcount(m & ~cofactor_masks[i]);
		int pfactor = std::popcount(m & cofactor_masks[i]);

		if (nfactor > pfactor) {
			std::swap(pfactor, nfactor);
			compls[i] = true;
		}

		int j;
		for (j = 0; j < i; j++)
		if (nfactor < popcount[j]) {
			break;
		}

		for (int k = i - 1; k >= j; k--) {
			popcount[k + 1] = popcount[k];
			order[k + 1] = order[k];
		}

		popcount[j] = nfactor;
		order[j] = i;
	}

	truth6 sc = 0;

	for (int idx1 = 0; idx1 < nbits; idx1++) {
		if (!(m & (truth6) 1 << idx1))
			continue;

		int idx2 = 0;
		for (int j = 0; j < ninputs; j++) {
			int pj = order[j];
			if ((!!(idx1 & 1 << pj)) ^ compls[pj])
				idx2 |= 1 << j;
		}

		sc |= (truth6) 1 << idx2;
	}

	for (int i = 0; i < ninputs; i++) {
		npn.ic[i] = compls[i];
		npn.p[order[i]] = i;
	}

	return sc;
}

void npn_semiclass_allrepr(truth6 m, int ninputs,
						   std::function<void(truth6, NPN&)> cb)
{
	NPN npn = {};

	if (!ninputs) {
		cb(0, npn);
		return;
	}

	bool bipo = false;
	int nbits = 1 << ninputs;
	truth6 mask = ((truth6) 2 << (nbits - 1)) - 1; // damn you, C++!
	m &= mask;
	if (std::popcount(m) > nbits / 2) {
		npn.oc = true;
		m ^= mask;
	} else if (std::popcount(m) == nbits / 2) {
		bipo = true;
	}

repolarize:
	bool compls[6] = {};
	int popcount[6];
	int order[6];
	int ambimask = nbits - 1;

	for (int i = 0; i < ninputs; i++) {
		int nfactor = std::popcount(m & ~cofactor_masks[i]);
		int pfactor = std::popcount(m & cofactor_masks[i]);

		if (nfactor > pfactor) {
			std::swap(pfactor, nfactor);
			compls[i] = true;
		} else if (nfactor == pfactor) {
			ambimask &= ~(1 << i);
		}

		int j;
		for (j = 0; j < i; j++)
		if (nfactor < popcount[j]) {
			break;
		}

		for (int k = i - 1; k >= j; k--) {
			popcount[k + 1] = popcount[k];
			order[k + 1] = order[k];
		}

		popcount[j] = nfactor;
		order[j] = i;
	}

	int tied[6] = {};

	for (int i = 0; i < ninputs - 1; i++) {
		int mark = i;
		for (; popcount[i] == popcount[i + 1] && i < ninputs - 1;) {
			tied[mark]++; i++;
		}
	}

	for (int ambi = 0; ambi < nbits; ambi = ((ambi | ambimask) + 1) & ~ambimask)
	{
next_round:
		{
			truth6 sc = 0;
			for (int idx1 = 0; idx1 < nbits; idx1++) {
				if (!(m & (truth6) 1 << idx1))
					continue;

				int idx2 = 0;
				for (int j = 0; j < ninputs; j++) {
					int pj = order[j];
					if ((!!(idx1 & 1 << pj)) ^ compls[pj] ^ !!(ambi & 1 << pj))
						idx2 |= 1 << j;
				}

				sc |= (truth6) 1 << idx2;
			}


			for (int i = 0; i < ninputs; i++) {
				npn.ic[i] = compls[i] ^ !!(ambi & 1 << i);
				npn.p[order[i]] = i;
			}

			cb(sc, npn);
		}

		for (int j = ninputs - 2; j >= 0; j--)
		if (tied[j]) {
			if (std::next_permutation(std::begin(order) + j,
									  std::begin(order) + j + tied[j] + 1))
				goto next_round;
		}
	}
	if (bipo) {
		bipo = false;
		m ^= mask;
		npn.oc ^= true;
		goto repolarize;
	}
}

#ifdef NPN_MAIN
int main(int argc, char const *argv[])
{
	int K = atoi(argv[1]);
	printf("Classifying functions of %d inputs\n", K);

	std::set<truth6> seen;
	int unique = 0;
	truth6 mask = ((truth6) 2 << (K - 1)) - 1;

	for (uint64_t i = 0; i < (uint64_t) 1 << (1 << K); i++) {
		NPN npn;

		truth6 sc = npn_semiclass(i, K, npn);

		assert((npn(i) & mask) == (sc & mask));
		assert(((npn.inv())(sc) & mask) == (i & mask));
		assert((npn.inv() * npn).is_identity());
		assert((npn * npn.inv()).is_identity());

		if (seen.count(sc))
			continue;

		unique++;

		bool found = false;
		npn_semiclass_allrepr(i, K, [&](truth6 sc_any, NPN &npn_any) {

			assert((npn_any(i) & mask) == (sc_any & mask));
			assert(((npn_any.inv())(sc_any) & mask) == (i & mask));

			seen.insert(sc_any);
			if (sc == sc_any)
				found = true;
		});
		assert(found);
	}

	printf("Found %d unique classes\n", unique);

	return 0;
}
#endif

#ifdef NPN_MAIN2
int main(int argc, char const *argv[])
{
	printf("Class %llx\n", npn_semiclass(0xff00000000000000, 6));

	return 0;
}
#endif


