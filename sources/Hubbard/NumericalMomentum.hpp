#include <iterator>
#include <numeric>
#include <array>
#include "GlobalDefinitions.hpp"
#include "Constants.hpp"

namespace Hubbard {
	template<unsigned int Dimension>
	struct NumericalMomentum {
		global_floating_type momenta[Dimension];
		std::array<int, Dimension> k;

		NumericalMomentum()
		{
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = -BASE_PI;
				k[d] = -Constants::K_DISCRETIZATION;
			}
		};
		constexpr explicit NumericalMomentum(const std::array<int, Dimension>& point_in_bz) : k{ point_in_bz }
		{
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = (k[d] * BASE_PI) / Constants::K_DISCRETIZATION;
			}
		};
		constexpr explicit NumericalMomentum(std::array<int, Dimension>&& point_in_bz) : k{ point_in_bz }
		{
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = (k[d] * BASE_PI) / Constants::K_DISCRETIZATION;
			}
		};

		template <bool shiftByPi>
		constexpr explicit NumericalMomentum(const Eigen::Vector<int, Dimension>& point_in_bz) {
			for (size_t d = 0U; d < Dimension; ++d)
			{
				if constexpr (shiftByPi) {
					k[d] = point_in_bz(d) - Constants::K_DISCRETIZATION;
				}
				else {
					k[d] = point_in_bz(d);
				}
				momenta[d] = (k[d] * BASE_PI) / Constants::K_DISCRETIZATION;
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
			_increment<0>();
			return *this;
		};

		inline bool iterateHalfBZ() {
			_increment<0>();
			return (k[Dimension - 1] < 0);
		};

		inline bool iterateFullBZ() {
			_increment<0>();
			return (k[Dimension - 1] < Constants::K_DISCRETIZATION);
		};

	private:
		template <int _d>
		inline void _increment() {
			static_assert(_d < Dimension, "Call to increment in NumericalMomentum provides a too high dimension.");
			if (++k[_d] >= Constants::K_DISCRETIZATION) {
				k[_d] = -Constants::K_DISCRETIZATION;
				_increment<_d + 1>();
			}
			momenta[_d] = (k[_d] * BASE_PI) / Constants::K_DISCRETIZATION;
		};

		template <>
		inline void _increment<Dimension - 1>() {
			momenta[_d] = (++k[_d] * BASE_PI) / Constants::K_DISCRETIZATION;
		};
	};
}