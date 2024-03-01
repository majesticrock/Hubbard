#pragma once
#include "GlobalDefinitions.hpp"
#include "Constants.hpp"
#include "ModelAttributes.hpp"
#include <string>
#include <memory>

namespace Hubbard {
	// maps an index; [0, N_K) -> [-pi, pi)
	inline global_floating_type index_to_k_vector(const int index) {
		return ((index * Constants::PI_DIV_DISCRETIZATION) - BASE_PI);
	};
	template <int Dimension>
	inline Eigen::Vector<global_floating_type, Dimension> index_vector_to_k_vector(const Eigen::Vector<int, Dimension>& index_vec) {
		Eigen::Vector<global_floating_type, Dimension> ret;
		for (size_t i = 0U; i < Dimension; ++i)
		{
			ret(i) = index_to_k_vector(index_vec(i));
		}
		return ret;
	};

	template <typename DataType>
	class BaseModel
	{
	public:
		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using BaseAttributes = ModelAttributes<DataType>;
		using HamiltonSolver = Eigen::SelfAdjointEigenSolver<SpinorMatrix>;
	private:
		inline void init()
		{
			this->hamilton = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);

			computeChemicalPotential();
			this->parameterCoefficients = {
				0.5 * this->U_OVER_N - 4. * this->V_OVER_N, // CDW
				0.5 * this->U_OVER_N, // AFM
				this->U_OVER_N, // SC
				this->V_OVER_N, // Gamma SC
				this->V_OVER_N, // Xi SC
				this->U_OVER_N, // Eta
				this->V_OVER_N, // Occupation Up
				this->V_OVER_N, // Occupation Down
				(0.5 * this->U_OVER_N + 4. * this->V_OVER_N) // Phase seperation
			};
		};
		inline void multiplyParametersByCoefficients(ParameterVector& F) const {
			for (size_t i = 0U; i < this->model_attributes.size(); ++i)
			{
				F(i) *= parameterCoefficients[i];
			}
		};

	protected:
		ModelAttributes<DataType> model_attributes;
		// Stores the coefficients for the parameters (e.g. V/N) with the appropriate index
		std::vector<global_floating_type> parameterCoefficients = std::vector<global_floating_type>(9U);

		SpinorMatrix hamilton;
		SpinorMatrix rho;
		HamiltonSolver hamilton_solver;

		double temperature{};
		double U{};
		double V{};
		double U_OVER_N{ U / Constants::BASIS_SIZE };
		double V_OVER_N{ V / Constants::BASIS_SIZE };
		double chemical_potential{};

		size_t SPINOR_SIZE{ 4U };

		inline void setParameters(ParameterVector& F) {
			_CONST_FLOATING new_weight{ 0.5 };
			for (size_t i = 0U; i < this->model_attributes.size(); ++i)
			{
				this->model_attributes[i] = new_weight * F(i) + (1 - new_weight) * this->model_attributes[i];
				// Numerical noise correction
				if (abs(this->model_attributes[i]) < DEFAULT_PRECISION * 1e-2) this->model_attributes[i] = DataType{};
			}
		};

		virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};

		inline global_floating_type fermi_dirac(global_floating_type energy) const {
			if (temperature > DEFAULT_PRECISION) {
				return (1. / (1. + exp(energy / temperature)));
			}
			else {
				if (abs(energy) < DEFAULT_PRECISION) {
					return global_floating_type{ 0.5 };
				}
				return ((energy > 0) ? global_floating_type{} : global_floating_type{ 1 });
			}
		};
		inline void fillRho() {
			this->hamilton_solver.compute(this->hamilton);
			rho.fill(DataType{});
			for (size_t i = 0U; i < this->SPINOR_SIZE; ++i)
			{
				rho(i, i) = 1 - fermi_dirac(hamilton_solver.eigenvalues()(i));
			}
			rho = hamilton_solver.eigenvectors() * rho * hamilton_solver.eigenvectors().adjoint();
		};
		inline void applyIteration(ParameterVector& F) {
			this->multiplyParametersByCoefficients(F);
			this->setParameters(F);
		}

		// Utility
		// To be used in the standard 4x4 spinor representation - does not work if the representation is changed
		inline global_floating_type get_n_up() const {
			return 1 - this->rho(0, 0).real();
		};
		inline global_floating_type get_n_down() const {
			return this->rho(2, 2).real();
		};
		inline global_floating_type get_n_up_Q() const {
			return 1 - this->rho(1, 1).real();
		};
		inline global_floating_type get_n_down_Q() const {
			return this->rho(3, 3).real();
		};
		inline global_floating_type get_n_up_plus_down() const {
			return get_n_up() + get_n_down();
		};
		inline complex_prec get_g_up() const {
			return -this->rho(1, 0);
		};
		inline complex_prec get_g_down() const {
			return this->rho(2, 3);
		};
		inline complex_prec get_g_up_plus_down() const {
			return get_g_up() + get_g_down();
		};
		inline complex_prec get_f() const {
			return -this->rho(0, 2);
		};
		inline complex_prec get_f_Q() const {
			return -this->rho(1, 3);
		};
		inline complex_prec get_eta() const {
			return -this->rho(0, 3);
		};

	public:
		explicit BaseModel(const ModelParameters& params, SystemType sytemType = SystemUndefined)
			: model_attributes(params, sytemType), temperature(params.temperature), U(params.U), V(params.V)
		{
			init();
		};

		BaseModel(const ModelParameters& params, const size_t _spinor_size, const size_t number_of_attributes)
			: model_attributes(number_of_attributes), parameterCoefficients(_spinor_size), temperature(params.temperature), U(params.U), V(params.V), SPINOR_SIZE(_spinor_size)
		{
			this->hamilton = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);

			computeChemicalPotential();
		};

		BaseModel(const ModelParameters& _params, const ModelAttributes<DataType>& startingValues)
			: model_attributes(startingValues), temperature(_params.temperature), U(_params.U), V(_params.V)
		{
			init();
		};
		virtual ~BaseModel() = default;

		void setNewModelParameters(const ModelParameters& params, SystemType systemType) {
			this->temperature = params.temperature;
			this->U = params.U;
			this->V = params.V;
			this->model_attributes.reset(params, systemType);
		};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) = 0;

		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy debugPolicy = NoWarning) = 0;

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) = 0;

		virtual void computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all) = 0;

		inline auto getTotalGapValue() const {
			return this->model_attributes.getTotalGapValue();
		}
		inline virtual global_floating_type entropyPerSite() = 0;
		inline virtual global_floating_type internalEnergyPerSite() = 0;
		inline global_floating_type freeEnergyPerSite() {
			return this->internalEnergyPerSite() - temperature * this->entropyPerSite();
		};

		inline std::string parametersAsTriplet() const {
			return ("[T U V] = [" + to_string(temperature) + " " + to_string(U) + " " + to_string(V) + "]");
		}

		inline const ModelAttributes<DataType>& getAttributes() const {
			return this->model_attributes;
		}

		inline void set_CDW_SC_ratio(double cdw_ratio) {
			auto buffer = this->model_attributes[0] + this->model_attributes[2];
			this->model_attributes[0] = buffer * cdw_ratio;
			this->model_attributes[2] = buffer * (1. - cdw_ratio);
		};
	};
}