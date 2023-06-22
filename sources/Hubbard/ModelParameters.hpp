#pragma once
#include <string>

namespace Hubbard {
	class ModelParameters {
	private:
		std::string global_iterator_type;
		std::string second_iterator_type;
		double global_step;
		double second_step;
		double global_it_min;
		double second_it_min;

		void incrementer(std::string& s, const double step);
	public:
		double temperature;
		double U;
		double V;

		ModelParameters(double _temperature, double _U, double _V, double global_step, double second_step,
			std::string _global_iterator_type, std::string _second_iterator_type);
		// This is just a placeholder - a class initialized this way should not be used. Ever.
		ModelParameters() : global_iterator_type(""), second_iterator_type(""), global_step(-1), second_step(-1),
			global_it_min(-1), second_it_min(-1), temperature(-1), U(-1), V(-1) { };

		double setGlobalIterator(int it_num);
		double setGlobalIteratorExact(double newValue);
		double setSecondIterator(int it_num);
		double setSecondIteratorExact(double newValue);
		void incrementGlobalIterator();
		void incrementSecondIterator();
		inline double getSecondStep() const {
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
		std::string getFileName() const;
	};
}