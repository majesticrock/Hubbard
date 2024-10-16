#pragma once
#include "../GlobalDefinitions.hpp"
#include "../Constants.hpp"
#include "../OrderType.hpp"
#include "ModelAttributes.hpp"
#include <string>
#include <memory>
#include <algorithm>

namespace Hubbard::Models {
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

	protected:
		virtual void init() = 0;
		inline void multiplyParametersByCoefficients(ParameterVector& F) const {
			for (size_t i = 0U; i < this->model_attributes.size(); ++i)
			{
				F(i) *= parameterCoefficients[i];
			}
		};

		ModelAttributes<DataType> model_attributes;
		// Stores the coefficients for the parameters (e.g. V/N) with the appropriate index
		std::vector<global_floating_type> parameterCoefficients = std::vector<global_floating_type>(9U);

		SpinorMatrix hamilton;
		SpinorMatrix rho;
		HamiltonSolver hamilton_solver;

	public:
		coefficient_type temperature{};
		coefficient_type U{};
		coefficient_type V{};
	protected:
		coefficient_type U_OVER_N{ U / Constants::BASIS_SIZE };
		coefficient_type V_OVER_N{ V / Constants::BASIS_SIZE };
		coefficient_type chemical_potential{};

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

		virtual void computeChemicalPotential() = 0;

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
			rho.setZero(this->SPINOR_SIZE, this->SPINOR_SIZE);
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
		inline DataType get_g_up() const {
			if constexpr (Utility::is_complex<DataType>()) {
				return -this->rho(1, 0);
			}
			else {
				return -std::real(this->rho(1, 0));
			}
		};
		inline DataType get_g_down() const {
			if constexpr (Utility::is_complex<DataType>()) {
				return this->rho(2, 3);
			}
			else {
				return std::real(this->rho(2, 3));
			}
		};
		inline DataType get_g_up_plus_down() const {
			return get_g_up() + get_g_down();
		};
		inline DataType get_f() const {
			if constexpr (Utility::is_complex<DataType>()) {
				return -this->rho(0, 2);
			}
			else {
				return -std::real(this->rho(0, 2));
			}
		};
		inline DataType get_f_Q() const {
			if constexpr (Utility::is_complex<DataType>()) {
				return -this->rho(1, 3);
			}
			else {
				return -std::real(this->rho(1, 3));
			}
		};
		inline DataType get_eta() const {
			if constexpr (Utility::is_complex<DataType>()) {
				return -this->rho(0, 3);
			}
			else {
				return -std::real(this->rho(0, 3));
			}
		};

	public:
		explicit BaseModel(const ModelParameters& params, SystemType sytemType = SystemUndefined)
			: model_attributes(params, sytemType), temperature(params.temperature), U(params.U), V(params.V)
		{ };

		BaseModel(const ModelParameters& params, const size_t _spinor_size, const size_t number_of_attributes)
			: model_attributes(number_of_attributes), parameterCoefficients(_spinor_size), temperature(params.temperature), U(params.U), V(params.V), SPINOR_SIZE(_spinor_size)
		{
			this->hamilton = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
		};

		BaseModel(const ModelParameters& _params, const ModelAttributes<DataType>& startingValues)
			: model_attributes(startingValues), temperature(_params.temperature), U(_params.U), V(_params.V)
		{};
		virtual ~BaseModel() = default;

		inline OrderType get_order() const {
			OrderType order{};
			if (this->model_attributes.isFinite(0)) order |= OrderType::CDW;
			if (this->model_attributes.isFinite(1)) order |= OrderType::AFM;
			if (this->model_attributes.isFinite(2)) order |= OrderType::SC;
			return order;
		}

		void setNewModelParameters(const ModelParameters& params, SystemType systemType) {
			this->temperature = params.temperature;
			this->U = params.U;
			this->V = params.V;
			this->U_OVER_N = params.U / Constants::BASIS_SIZE;
			this->V_OVER_N = params.V / Constants::BASIS_SIZE;
			this->init();

			this->model_attributes.reset(params, systemType);
			computeChemicalPotential();
		};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) = 0;

		virtual ModelAttributes<global_floating_type> computePhases() = 0;

		// saves all one particle energies to reciever
		virtual void getAllEnergies(std::vector<global_floating_type>& reciever) = 0;

		std::vector<global_floating_type> continuum_boundaries() {
			std::vector<global_floating_type> energy;
			this->getAllEnergies(energy);
			for (auto& en : energy) {
				if (en < 0) en = std::abs(en);
			}
			return { 2. * (*std::ranges::min_element(energy)), 2. * (*std::ranges::max_element(energy)) };
		}

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

		inline std::string info() const {
			return this->parametersAsTriplet();
		}

		inline const ModelAttributes<DataType>& getAttributes() const {
			return this->model_attributes;
		}

		inline void set_CDW_SC_ratio(global_floating_type cdw_ratio) {
			auto buffer = this->model_attributes[0] + this->model_attributes[2];
			this->model_attributes[0] = buffer * cdw_ratio;
			this->model_attributes[2] = buffer * (1. - cdw_ratio);
		};
	};
}