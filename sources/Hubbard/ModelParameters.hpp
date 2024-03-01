#pragma once
#include <string>
#include <vector>
#include <iostream>
#include "GlobalDefinitions.hpp"

namespace Hubbard {
	class ModelParameters {
	private:
		std::string global_iterator_type{ "N/A" };
		std::string second_iterator_type{ "N/A" };
		coefficient_type global_step{ -1 };
		coefficient_type second_step{ -1 };
		coefficient_type global_it_min{ -1 };
		coefficient_type second_it_min{ -1 };

		void init();
		void incrementer(std::string& s, const coefficient_type step);
	public:
		coefficient_type temperature{ -1 };
		coefficient_type U{ -1 };
		coefficient_type V{ -1 };

		ModelParameters(coefficient_type _temperature, coefficient_type _U, coefficient_type _V, coefficient_type _global_step, coefficient_type _second_step,
			std::string _global_iterator_type, std::string _second_iterator_type);
		ModelParameters(const std::vector<coefficient_type>& params, coefficient_type _global_step, coefficient_type _second_step,
			std::string _global_iterator_type, std::string _second_iterator_type);
		// This is offered as a simpler way of initializing the class.
		ModelParameters() = default;

		// Signed is required as the value can be negative (findSingleBoundary in PhaseHelper)
		coefficient_type setGlobalIterator(int it_num);
		coefficient_type setGlobalIteratorExact(coefficient_type newValue);
		coefficient_type setSecondIterator(int it_num);
		coefficient_type setSecondIteratorExact(coefficient_type newValue);
		void incrementGlobalIterator();
		void incrementSecondIterator();
		inline coefficient_type getSecondStep() const noexcept {
			return second_step;
		};
		inline coefficient_type getGlobal() const {
			if (global_iterator_type == "T") {
				return temperature;
			}
			else if (global_iterator_type == "U") {
				return U;
			}
			else if (global_iterator_type == "V") {
				return V;
			}
			return -128;
		};
		inline coefficient_type getSecond() const {
			if (second_iterator_type == "T") {
				return temperature;
			}
			else if (second_iterator_type == "U") {
				return U;
			}
			else if (second_iterator_type == "V") {
				return V;
			}
			return -128;
		};
		inline void reset() {
			setSecondIterator(0);
			setGlobalIterator(0);
		};
		void printGlobal() const;
		void printSecond() const;
		void printParameters() const;
		std::string getFolderName() const;
	};

	inline std::ostream& operator<<(std::ostream& os, const ModelParameters& mp) {
		os << mp.temperature << "\t" << mp.U << "\t" << mp.V;
		return os;
	}
}