#pragma once
#include <string>

namespace Hubbard {
	class ModelParameters {
	private:
		std::string global_iterator_type{"N/A"};
		std::string second_iterator_type{"N/A"};
		double global_step{-1};
		double second_step{-1};
		double global_it_min{-1};
		double second_it_min{-1};

		void incrementer(std::string& s, const double step);
	public:
		double temperature{-1};
		double U{-1};
		double V{-1};

		ModelParameters(double _temperature, double _U, double _V, double global_step, double second_step,
			std::string _global_iterator_type, std::string _second_iterator_type);
		// This is offered as a simpler way of initializing the class.
		ModelParameters() = default;

		// Signed is required as the value can be negative (findSingleBoundary in PhaseHelper)
		double setGlobalIterator(int it_num);
		double setGlobalIteratorExact(double newValue);
		double setSecondIterator(int it_num);
		double setSecondIteratorExact(double newValue);
		void incrementGlobalIterator();
		void incrementSecondIterator();
		inline double getSecondStep() const noexcept {
			return second_step;
		};
		inline double getGlobal() const {
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
		inline double getSecond() const {
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
		void printParameters() const;
		std::string getFileName() const;
	};
}