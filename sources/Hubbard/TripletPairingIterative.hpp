#pragma once
#include "HubbardCDW.hpp"

namespace Hubbard {
	class TripletPairingIterative :
		public HubbardCDW
	{
	private:
		typedef Eigen::Vector<complex_prec, 18> ParameterVector;
	protected:
		complex_prec tau_sc, theta_sc;

		inline double_prec tau(double_prec k_x, double_prec k_y) {
			return sin(k_x) + sin(k_y);
		};
		inline double_prec theta(double_prec k_x, double_prec k_y) {
			return sin(k_x) - sin(k_y);
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		void addToParameterSet(const SpinorMatrix& rho, Eigen::Ref<ParameterVector> F, double_prec k_x, double_prec k_y) {
			HubbardCDW::addToParameterSet(rho, F.head<16>(), k_x, k_y);

			F(16) -= tau(k_x, k_y) * rho(6, 2);
			F(17) -= theta(k_x, k_y) * rho(6, 2);
		};
	public:
		TripletPairingIterative(const ModelParameters& _params);
		virtual Model::data_set computePhases(const bool print = false) override;
	};
}