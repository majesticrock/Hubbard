#pragma once
#include "../GlobalDefinitions.hpp"
#include "../ModelAttributes.hpp"
#include "PhaseHelper.hpp"

namespace Hubbard::Helper {
	struct Plaquette {
		/* Layout:
		*  0 1
		*  2 3
		*/
		std::array<ModelAttributes<global_floating_type>, 4> attributes;
		PhaseHelper* parent{};

		double lowerFirst{};
		double upperFirst{};
		double lowerSecond{};
		double upperSecond{};
		int value_index{};

		inline global_floating_type& operator()(const size_t i, const size_t j) {
			assert(i < 4);
			return attributes[i][j];
		};
		inline const global_floating_type& operator()(const size_t i, const size_t j) const {
			assert(i < 4);
			return attributes[i][j];
		};

		inline bool valueIsFinite(const size_t index) const {
			return abs(attributes[index][value_index]) > 1e-12;
		};
		inline bool containsPhaseBoundary() const
		{
			for (size_t i = 1U; i < 4U; ++i)
			{
				if (valueIsFinite(0U) != valueIsFinite(i)) {
					return true;
				}
			}
			return false;
		};
		inline double getCenterFirst() const noexcept {
			return 0.5 * (lowerFirst + upperFirst);
		};
		inline double getCenterSecond() const noexcept {
			return 0.5 * (lowerSecond + upperSecond);
		};
		inline double size() const noexcept {
			return 0.5 * (abs(upperFirst - lowerFirst) + abs(upperSecond - lowerSecond));
		}
		// Devides *this into 4 equally sized plaquettes
		// The result is appended to <appendTo>
		void devidePlaquette(std::vector<Plaquette>& appendTo);
	};
}

