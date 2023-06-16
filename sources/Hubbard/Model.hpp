#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <map>
#include <cmath>
#include <optional>
#include "Constants.hpp"
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include "../Utility/Resolvent.hpp"

namespace Hubbard {
	// Defines the working precision of the entire project
	// Change to float, double or long double - so far double produces the best results
	typedef double double_prec;
	constexpr double_prec L_PI = 3.141592653589793238462643383279502884L; /* pi */
	typedef Eigen::Matrix<double_prec, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<double_prec, Eigen::Dynamic> Vector_L;
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<std::complex<double>, Eigen::Dynamic> VectorCL;
	typedef Utility::Resolvent<double> Resolvent_L;

	class Model
	{
	protected:
		size_t TOTAL_BASIS;
		double_prec delta_sc, delta_cdw, delta_afm, delta_eta;
		double_prec gamma_sc, xi_sc;
		double_prec delta_occupation_up, delta_occupation_down;
		double_prec chemical_potential;

		typedef std::complex<double_prec> complex_prec;
		typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> SpinorMatrix;
		const complex_prec I = { 0, 1 };
		SpinorMatrix hamilton;

		double_prec temperature;
		double_prec U;
		double_prec V;
		double_prec U_OVER_N;
		double_prec V_OVER_N;

		inline virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};
		inline virtual double_prec gamma(double_prec k_x, double_prec k_y) const {
			return cos(k_x) + cos(k_y);
		}
		inline virtual double_prec xi(double_prec k_x, double_prec k_y) const {
			return cos(k_x) - cos(k_y);
		}

		inline double_prec unperturbed_energy(double_prec k_x, double_prec k_y) const {
			return -2 * gamma(k_x, k_y);
		};
		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const {
			return unperturbed_energy(k_x, k_y);
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const {
			return unperturbed_energy(k_x, k_y);
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) = 0;
		inline double_prec fermi_dirac(double_prec energy) const {
			if (temperature > 1e-8) {
				return (1. / (1 + exp(energy / temperature)));
			}
			else {
				if (std::abs(energy) < 1e-12) {
					return 0.5;
				}
				return ((energy > 0) ? 0 : 1);
			}
		};
		inline void fillRho(SpinorMatrix& rho, const Eigen::SelfAdjointEigenSolver<SpinorMatrix>& solvedHamilton) const {
			rho.fill(0);
			for (int i = 0; i < rho.rows(); i++)
			{
				rho(i, i) = 1 - fermi_dirac(solvedHamilton.eigenvalues()(i));
			}
			rho = solvedHamilton.eigenvectors() * rho * solvedHamilton.eigenvectors().adjoint();
		};

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

		class ModelParameters {
		private:
			std::string global_iterator_type;
			std::string second_iterator_type;
			double_prec global_step;
			double_prec second_step;
			double_prec global_it_min;
			double_prec second_it_min;

			void incrementer(std::string& s, const double_prec step);
		public:
			double_prec temperature;
			double_prec U;
			double_prec V;

			ModelParameters(double_prec _temperature, double_prec _U, double_prec _V, double_prec global_step, double_prec second_step,
				std::string _global_iterator_type, std::string _second_iterator_type);
			// This is just a placeholder - a class initialized this way should not be used. Ever.
			ModelParameters() : global_iterator_type(""), second_iterator_type(""), global_step(-1), second_step(-1),
				global_it_min(-1), second_it_min(-1), temperature(-1), U(-1), V(-1) { };

			double_prec setGlobalIterator(int it_num);
			double_prec setGlobalIteratorExact(double_prec newValue);
			double_prec setSecondIterator(int it_num);
			double_prec setSecondIteratorExact(double_prec newValue);
			void incrementGlobalIterator();
			void incrementSecondIterator();
			inline double_prec getSecondStep() const {
				return second_step;
			};
			inline double_prec getGlobal() const {
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
			inline double_prec getSecond() const {
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
			void print() const;
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