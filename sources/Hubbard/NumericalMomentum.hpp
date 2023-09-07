#include <iterator>
#include <numeric>
#include "GlobalDefinitions.hpp"

namespace Hubbard {
	template<unsigned int Dimension>
	struct NumericalMomentum {
		global_floating_type momenta[Dimension];

		inline global_floating_type tau() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), global_floating_type{}, [](global_floating_type current, global_floating_type toAdd) {
				return current + sin(toAdd);
				});
		};
		inline global_floating_type gamma() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), global_floating_type{}, [](global_floating_type current, global_floating_type toAdd) {
				return current + cos(toAdd);
				});
		};
		inline global_floating_type unperturbed_energy() const {
			return -2. * gamma();
		};

		inline global_floating_type& operator[](unsigned int index) {
			assert(index < Dimension);
			return momenta[index];
		};
		inline const global_floating_type& operator[](unsigned int index) const {
			assert(index < Dimension);
			return momenta[index];
		};
	};

	template<unsigned int Dimension>
	inline std::ostream& operator<<(std::ostream& os, const NumericalMomentum<Dimension>& momentum) {
		os << momentum.momenta[0];
		for (size_t d = 1U; d < Dimension; ++d)
		{
			os << " " << momentum.momenta[d];
		}
		return os;
	}
}