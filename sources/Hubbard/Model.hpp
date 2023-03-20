#pragma once
#define _USE_MATH_DEFINES

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
	class Model
	{
	protected:
		int BASIS_SIZE;
		int TOTAL_BASIS;
		double delta_sc, delta_cdw, delta_eta;
		Eigen::MatrixXd hamilton;
		double temperature;
		double U;

		double chemical_potential;
		// Might not be necessary, depends on the V dependence
		virtual void computeChemicalPotential();

		// maps an index; [0, N_K) -> [-pi, pi)
		template <typename T>
		inline double index_to_k_vector(const T index) const {
			return (((index * M_PI) / Constants::K_DISCRETIZATION) - M_PI);
		};

		/*
		* 0 - number operator
		* 1 - cdw
		* 2 - sc
		* 3 - eta
		*/
		std::vector<Eigen::MatrixXd> expecs;
		double sum_of_all[4] = { 0, 0, 0, 0 };
		// Quartic expecs - contracted using Wick's theorem (exact on MF level)
		std::vector<Eigen::MatrixXd> quartic;

		Eigen::MatrixXd M, N;
		int number_of_basis_terms;
		int start_basis_at;

		const std::map<std::string, int> wick_map = { {"n", 0}, {"g", 1}, {"f", 2}, {"\\eta", 3} };
		std::vector<std::vector<SymbolicOperators::WickTerm>> wicks_M, wicks_N;

		inline double fermi_dirac(double energy) const {
			//energy += chemical_potential;
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
		inline double unperturbed_energy(double k_x, double k_y) const {
			return -2 * (cos(k_x) + cos(k_y));// - chemical_potential;
		};
		virtual void fillHamiltonian(double k_x, double k_y) = 0;

		// Computes the respective x or y component from a given input index
		inline int x(int idx) const {
			return idx / (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		inline int y(int idx) const {
			return idx % (2 * Constants::K_DISCRETIZATION) - Constants::K_DISCRETIZATION;
		};
		inline int equal_up_to_Q(const Eigen::Vector2i& l, const Eigen::Vector2i& r) const {
			if (l == r) return 0;

			if (l(0) == r(0) + Constants::K_DISCRETIZATION || l(0) == r(0) - Constants::K_DISCRETIZATION) {
				if (l(1) == r(1) + Constants::K_DISCRETIZATION || l(1) == r(1) - Constants::K_DISCRETIZATION) {
					return 1;
				}
			}
			return -1;
		};

		// returns a value in [0, N_K), note that N_K = 2*constants::k_disc
		inline void clean_factor_2pi(Eigen::Vector2i& toClean) const {
			// + Q is required for the modulo operation later
			// as well as referencing, which works on indizes from 0 to [2pi] and not from [-pi] to [pi]
			for (int i = 0; i < 2; i++)
			{
				toClean(i) += Constants::K_DISCRETIZATION;
				if (toClean(i) < 0) {
					toClean(i) = ((2 * Constants::K_DISCRETIZATION) - std::abs(toClean(i) % (2 * Constants::K_DISCRETIZATION))) % (2 * Constants::K_DISCRETIZATION);
				}
				else {
					toClean(i) = (toClean(i) % (2 * Constants::K_DISCRETIZATION)) % (2 * Constants::K_DISCRETIZATION);
				}
			}
		};
		inline Eigen::Vector2i computeMomentum(const SymbolicOperators::Momentum& momentum, const std::vector<Eigen::Vector2i>& indizes, const std::vector<char>& momenta) const {
			Eigen::Vector2i buffer = { 0,0 };
			for (int i = 0; i < momenta.size(); ++i)
			{
				int mom_idx = momentum.isUsed(momenta[i]);
				if (mom_idx < 0) continue;
				buffer += momentum.momentum_list[mom_idx].first * indizes[i];
			}
			if (momentum.add_Q) {
				buffer(0) += Constants::K_DISCRETIZATION;
				buffer(1) += Constants::K_DISCRETIZATION;
			}
			clean_factor_2pi(buffer);
			return buffer;
		};

		virtual inline double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector2i& momentum) const {
			if (coeff.name == "\\epsilon_0") {
				//if (!(momentum.has_value())) throw std::length_error("Calling epsilon(k) without specifying k!");
				return (unperturbed_energy(index_to_k_vector(momentum(0)), index_to_k_vector(momentum(1))) - chemical_potential);
			}
			if (coeff.name == "\\frac{U}{N}") {
				return (U / BASIS_SIZE);
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};
		double computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const;
		virtual void fill_M_N();
	public:
		class ModelParameters {
		private:
			std::string global_iterator_type;
			std::string second_iterator_type;
			double global_step;
			double second_step;
			double second_it_min;

			void incrementer(std::string& s, const double step);
		public:
			double temperature;
			double U;
			double V;

			ModelParameters(double _temperature, double _U, double _V, double global_step, double second_step,
				std::string _global_iterator_type, std::string _second_iterator_type);
			void setSecondIterator(int it_num);
			void incrementGlobalIterator();
			void incrementSecondIterator();
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
			void printGlobal() const;
		};
		struct data_set {
			double delta_cdw, delta_sc, delta_eta;
			void print() const;
		};

		Model(double _temperature, double _U, int _number_of_basis_terms, int _start_basis_at);
		Model(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);
		// reciever is the vector the resulting data will be stored in
		// direction gives the angle between the k-path and the k_x axis in multiples of M_PI
		void getEnergies(std::vector<std::vector<double>>& reciever, double direction);
		// saves all one particle energies to reciever
		void getAllEnergies(std::vector<std::vector<double>>& reciever);
		virtual data_set computePhases(const bool print = false) = 0;
		// version 2 use the non mean field hamilton for the commutation,
		// but the mean field system to obtain the expectation values
		std::unique_ptr<Utility::Resolvent> computeCollectiveModes(std::vector<std::vector<double>>& reciever);
		// Returns the total gap value sqrt(sc^2 + cdw^2 + eta^2)
		inline double getTotalGapValue() const {
			return sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta);
		};
		void loadWick(const std::string& filename);
	};

	constexpr size_t triangular_number(const size_t i) {
		return (i * (i + 1) / 2);
	};
}