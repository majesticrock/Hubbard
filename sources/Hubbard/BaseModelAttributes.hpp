#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include "ModelParameters.hpp"

namespace Hubbard {
	template <typename DataType>
	struct BaseModelAttributes {
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
			this->gamma_occupation_up = _params.V * 0.2;
			this->gamma_occupation_down = _params.V * 0.2;

			this->gamma_sc = _params.V * 0.05;
			this->xi_sc = std::abs(_params.V) * 0.2;
		}
		void initializeMapper() {
			this->parameterMapper = {
				&(this->delta_cdw),
				&(this->delta_afm),
				&(this->delta_sc),
				&(this->gamma_sc),
				&(this->xi_sc),
				&(this->delta_eta),
				&(this->gamma_occupation_up),
				&(this->gamma_occupation_down),
			};
		};
	protected:
		std::vector<DataType*> parameterMapper;
	public:
		DataType delta_sc = 0., delta_eta = 0., gamma_sc = 0., xi_sc = 0.;
		DataType gamma_occupation_up = 0., gamma_occupation_down = 0., delta_cdw = 0., delta_afm = 0.;
		// Maps the parameters (delta_sc etc) to an index

		inline DataType& operator[](int i) {
			return *(parameterMapper[i]);
		};
		inline const DataType& operator[](int i) const {
			return *(parameterMapper[i]);
		};
		inline size_t size() const noexcept {
			return parameterMapper.size();
		}

		bool isOrdered() const {
			for (const auto p : parameterMapper)
			{
				if (std::abs(*p) > 1e-12) {
					return true;
				}
			}
			return false;
		};

		inline virtual double renormalizedEnergy_up(const double GAMMA) const = 0;
		inline virtual double renormalizedEnergy_down(const double GAMMA) const = 0;

		BaseModelAttributes(const ModelParameters& _params) {
			this->initializeMapper();
			this->initializeParamters(_params);
		};
		BaseModelAttributes() {
			this->initializeMapper();
		};
		BaseModelAttributes(const BaseModelAttributes& copy){
			this->initializeMapper();
			for (size_t i = 0; i < parameterMapper.size(); i++)
			{
				*(this->parameterMapper[i]) = *(copy.parameterMapper[i]);
			}
		};
		// Returns the total gap value sqrt(sc^2 + cdw^2 + eta^2)
		inline double getTotalGapValue() const {
			return sqrt(std::abs(delta_cdw) * std::abs(delta_cdw) + std::abs(delta_sc) * std::abs(delta_sc)
				+ std::abs(delta_eta) * std::abs(delta_eta));
		};

		inline bool isFinite(int i, double epsilon = 1e-12) const {
			return (std::abs((*this)[i]) > epsilon);
		}
		inline void print() const {
			for (const auto p : parameterMapper)
			{
				std::cout << *p << "\t";
			}
			
			std::cout << "\n    Delta_tot = " << getTotalGapValue() << std::endl;
		};
	};

	struct BaseModelRealAttributes;

	struct BaseModelComplexAttributes : public BaseModelAttributes<std::complex<double>> {
		inline virtual double renormalizedEnergy_up(const double GAMMA) const override {
			return -(2. + this->gamma_occupation_up.real()) * GAMMA;
		};
		inline virtual double renormalizedEnergy_down(const double GAMMA) const override {
			return -(2. + this->gamma_occupation_down.real()) * GAMMA;
		};
		BaseModelComplexAttributes();
		BaseModelComplexAttributes(const BaseModelComplexAttributes& copy);
		BaseModelComplexAttributes(const ModelParameters& _params);
		BaseModelComplexAttributes(const BaseModelRealAttributes& realValues);
	};

	struct BaseModelRealAttributes : public BaseModelAttributes<double> {
		inline BaseModelRealAttributes& operator+=(const BaseModelRealAttributes& rhs) {
			for (size_t i = 0; i < this->parameterMapper.size(); i++)
			{
				*(this->parameterMapper[i]) += *(rhs.parameterMapper[i]);
			}
			return *this;
		};
		inline BaseModelRealAttributes& operator*=(const double rhs) {
			for (size_t i = 0; i < this->parameterMapper.size(); i++)
			{
				*(this->parameterMapper[i]) *= rhs;
			}
			return *this;
		};
		inline BaseModelRealAttributes& operator/=(const double rhs) {
			for (size_t i = 0; i < this->parameterMapper.size(); i++)
			{
				*(this->parameterMapper[i]) /= rhs;
			}
			return *this;
		};

		inline virtual double renormalizedEnergy_up(const double GAMMA) const override {
			return -(2. + this->gamma_occupation_up) * GAMMA;
		};
		inline virtual double renormalizedEnergy_down(const double GAMMA) const override {
			return -(2. + this->gamma_occupation_down) * GAMMA;
		};

		BaseModelRealAttributes();
		BaseModelRealAttributes(const BaseModelRealAttributes& copy);
		BaseModelRealAttributes(const ModelParameters& _params);

		// This constructor is used for data output in the complex case.
		// We assume previous analytical knowledge about the the real and imaginary parts
		BaseModelRealAttributes(const BaseModelComplexAttributes& complexValues);
	};

	inline BaseModelRealAttributes operator+(BaseModelRealAttributes lhs, const BaseModelRealAttributes& rhs) {
		return lhs += rhs;
	};
	inline BaseModelRealAttributes operator*(BaseModelRealAttributes lhs, double rhs) {
		return lhs *= rhs;
	};
	inline BaseModelRealAttributes operator*(double lhs, BaseModelRealAttributes rhs) {
		return rhs *= lhs;
	};
	inline BaseModelRealAttributes operator/(BaseModelRealAttributes lhs, double rhs) {
		return lhs /= rhs;
	};
	inline BaseModelRealAttributes operator/(double lhs, BaseModelRealAttributes rhs) {
		return rhs /= lhs;
	};
}