#include <iterator>
#include <numeric>

namespace Hubbard {
	template<unsigned int Dimension>
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

		inline double& operator[](unsigned int index) {
			assert(index < Dimension);
			return momenta[index];
		};
		inline const double& operator[](unsigned int index) const {
			assert(index < Dimension);
			return momenta[index];
		};
	};
}