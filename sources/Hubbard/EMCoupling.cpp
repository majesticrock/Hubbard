#include "EMCoupling.hpp"

namespace Hubbard {
	void EMCoupling::fillHamiltonian()
	{
		hamilton.fill(global_floating_type{});

		NumericalMomentum<Dimension> k;
		NumericalMomentum<Dimension> q;

		size_t i, j;
		do {
			i = k.getIndex();
			hamilton(i, i) = -2. * k.gamma();
			do
			{
				j = q.getIndex();

				hamilton(i, j) = 3;
			} while (q.iterateFullBZ());
			q.reset();
		} while (k.iterateFullBZ());
	}
}
