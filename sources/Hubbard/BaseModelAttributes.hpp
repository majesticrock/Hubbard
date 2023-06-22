#pragma once
#include "BaseModel.hpp"

namespace Hubbard {
	template <typename DataType>
	class BaseModelAttributes {
	private:
		void initializeParamters(const ModelParameters& _params) {
			this->delta_cdw = (std::abs(_params.U) + _params.V) * 0.5 + 0.1;
			this->delta_sc = std::abs(_params.U + std::abs(_params.V)) * 0.3 + 0.05;
			if (_params.V > 0) {
				this->delta_sc *= 0;
			}
			else if (_params.V < 0) {
				this->delta_cdw *= 0;
			}
			if (_params.U > 0) {
				this->delta_afm = std::abs(_params.U - std::abs(_params.V)) * 0.5 + 0.1;
			}
			else {
				this->delta_afm = 0;
			}

			this->delta_eta = 0;//_params.U * 0.1;
			this->delta_occupation_up = _params.V * 0.2;
			this->delta_occupation_down = _params.V * 0.2;

			this->gamma_sc = _params.V * 0.05;
			this->xi_sc = std::abs(_params.V) * 0.2;
		}
	protected:
		DataType delta_sc, delta_eta, gamma_sc, xi_sc;
		DataType delta_occupation_up, delta_occupation_down, delta_cdw, delta_afm;

		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const = 0;
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const = 0;
	public:
		BaseModelAttributes(const ModelParameters& _params) {
			initializeParamters(_params);
		};
		// Returns the total gap value sqrt(sc^2 + cdw^2 + eta^2)
		inline double_prec getTotalGapValue() const {
			return sqrt(std::abs(delta_cdw) * std::abs(delta_cdw) + std::abs(delta_sc) * std::abs(delta_sc)
				+ std::abs(delta_eta) * std::abs(delta_eta));
		};
	};

	class BaseModelRealAttributes : public BaseModelAttributes<double_prec> {
	protected:
		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const override {
			return -2. * (1. + this->delta_occupation_up) * gamma(k_x, k_y);
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2. * (1. + this->delta_occupation_down) * gamma(k_x, k_y);
		};

	public:
		BaseModelRealAttributes(const ModelParameters& _params) : BaseModelAttributes(_params) {};
	};

	class BaseModelComplexAttributes : public BaseModelAttributes<complex_prec> {
	protected:
		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const override {
			return -2. * (1. + this->delta_occupation_up.real()) * gamma(k_x, k_y);
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2. * (1. + this->delta_occupation_down.real()) * gamma(k_x, k_y);
		};

	public:
		BaseModelComplexAttributes(const ModelParameters& _params) : BaseModelAttributes(_params) {};
	};
}