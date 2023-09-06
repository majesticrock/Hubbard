#include <iterator>
#include <numeric>
#include <array>
#include "GlobalDefinitions.hpp"
#include "Constants.hpp"

namespace Hubbard {
	template<unsigned int Dimension>
	struct NumericalMomentum {
		global_floating_type momenta[Dimension];
		int k[Dimension];

		NumericalMomentum()
		{
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = -BASE_PI;
				k[d] = -Constants::K_DISCRETIZATION;
			}
		};
		template <typename... Args>
    	constexpr NumericalMomentum(Args... args) : k{static_cast<int>(args)...} {
        	static_assert(sizeof...(Args) == Dimension, "Incorrect number of arguments");
        	static_assert(std::conjunction_v<std::is_integral<Args>...>, "All arguments must be integers");
        	
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = k[d] * Constants::PI_DIV_DISCRETIZATION;
			}
   		}

		constexpr explicit NumericalMomentum(const Eigen::Array<int, Dimension, 1>& point_in_bz) {
			for (size_t d = 0U; d < Dimension; ++d)
			{
				k[d] = point_in_bz(d);
				momenta[d] = k[d] * Constants::PI_DIV_DISCRETIZATION;
			}
		};
		constexpr explicit NumericalMomentum(Eigen::Array<int, Dimension, 1>&& point_in_bz) {
			for (size_t d = 0U; d < Dimension; ++d)
			{
				k[d] = std::move(point_in_bz(d));
				momenta[d] = k[d] * Constants::PI_DIV_DISCRETIZATION;
			}
		};

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

		inline NumericalMomentum& operator++() {
			_increment();
			return *this;
		};

		inline bool iterateHalfBZ() {
			_increment();
			return (k[Dimension - 1] < 0);
		};

		inline bool iterateFullBZ() {
			_increment();
			return (k[Dimension - 1] < Constants::K_DISCRETIZATION);
		};

	private:
		template <int _d=0>
		inline void _increment() {
			static_assert(_d < Dimension, "Call to increment in NumericalMomentum provides a too high dimension.");
			if (++k[_d] >= Constants::K_DISCRETIZATION) {
				k[_d] = -Constants::K_DISCRETIZATION;
				if constexpr (_d + 1 < Dimension){
					_increment<_d + 1>();
				}
			}
			momenta[_d] = k[_d] * Constants::PI_DIV_DISCRETIZATION;
		};
	};
}