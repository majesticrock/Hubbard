#include <iterator>
#include <numeric>

namespace Hubbard {
	template<size_t Dimension>
	struct NumericalMomentum {
		double momenta[Dimension];

		inline double tau() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), double{}, [](double current, double toAdd) {
				return current + sin(toAdd);
				});
		};
		inline double gamma() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), double{}, [](double current, double toAdd) {
				return current + cos(toAdd);
				});
		};
		inline double unperturbed_energy() const {
			return -2. * gamma();
		};

		inline double& operator[](size_t index) {
			assert(index < Dimension);
			return momenta[index];
		};
		inline const double& operator[](size_t index) const {
			assert(index < Dimension);
			return momenta[index];
		};
	};
}