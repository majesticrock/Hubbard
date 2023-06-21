#pragma once

#include "BaseModel.hpp"
#include <Eigen/Sparse>
#include <memory>
#include <map>
#include <cmath>
#include <optional>
#include "Constants.hpp"
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include "../Utility/Resolvent.hpp"

namespace Hubbard {
	typedef Utility::Resolvent<double_prec> Resolvent_L;
	class Model : public BaseModel
	{
	protected:
		double_prec gamma_sc, xi_sc;
		double_prec delta_occupation_up, delta_occupation_down;
		
		inline virtual double_prec xi(double_prec k_x, double_prec k_y) const {
			return cos(k_x) - cos(k_y);
		}

		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const {
			return -2 * (1. + delta_occupation_up) * gamma(k_x, k_y);
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const {
			return -2 * (1. + delta_occupation_down) * gamma(k_x, k_y);
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) = 0;

		// maps an index; [0, N_K) -> [-pi, pi)
		template <typename T>
		inline double_prec index_to_k_vector(const T index) const {
			return (((index * L_PI) / Constants::K_DISCRETIZATION) - L_PI);
		};

	public:
		virtual inline double_prec computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector2i& momentum) const {
			if (coeff.name == "\\epsilon_0") {
				//if (!(momentum.has_value())) throw std::length_error("Calling epsilon(k) without specifying k!");
				return (unperturbed_energy(index_to_k_vector(momentum(0)), index_to_k_vector(momentum(1))) - chemical_potential);
			}
			if (coeff.name == "\\frac{U}{N}") {
				return U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				//if (!(momentum.has_value())) throw std::invalid_argument("Calling V without specifying a momentum!");
				// Eventuell ein Faktor 2?
				return V_OVER_N * (cos(index_to_k_vector(momentum(0))) + cos(index_to_k_vector(momentum(1))));
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};

		struct data_set {
			bool converged = true;
			double_prec delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta;
			inline double_prec operator[](int i) const {
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
			inline bool isFinite(int i, double_prec epsilon = 1e-12) const {
				return (std::abs((*this)[i]) > epsilon);
			}
			inline void print() const{
				std::cout << delta_cdw << "\t" << delta_afm << "\t" << delta_sc << "\t" << delta_eta << "\t" << xi_sc
					<< "\n    Delta_tot = " << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
			};
		};

		Model(const  ModelParameters& _params);
		// reciever is the vector the resulting data will be stored in
		// direction gives the angle between the k-path and the k_x axis in multiples of L_PI
		void getEnergies(std::vector<std::vector<double>>& reciever, double_prec direction);
		// saves all one particle energies to reciever
		void getAllEnergies(std::vector<std::vector<double>>& reciever);
		virtual data_set computePhases(const bool print = false) = 0;

		void computeExpectationValues(std::vector<MatrixCL>& expecs, std::vector<complex_prec>& sum_of_all);
		// Returns the total gap value sqrt(sc^2 + cdw^2 + eta^2)
		inline double_prec getTotalGapValue() const {
			return sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta);
		};
	};
}