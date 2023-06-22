#pragma once
#include <iostream>
#include <cmath>

namespace Hubbard {
	struct PhaseDataSet
	{
		bool converged = true;
		double delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta;
		inline double operator[](int i) const {
			switch (i) {
			case 0:
				return delta_cdw;
			case 1:
				return delta_afm;
			case 2:
				return delta_sc;
			case 3:
				return gamma_sc;
			case 4:
				return xi_sc;
			case 5:
				return delta_eta;
			default:
				throw std::invalid_argument("ModelParameters[]: Index out of range.");
			}
		}
		inline bool isFinite(int i, double epsilon = 1e-12) const {
			return (std::abs((*this)[i]) > epsilon);
		}
		inline void print() const {
			std::cout << delta_cdw << "\t" << delta_afm << "\t" << delta_sc << "\t" << delta_eta << "\t" << xi_sc
				<< "\n    Delta_tot = " << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
		};
	};
}