#include "BaseModelAttributes.hpp"

namespace Hubbard {
	BaseModelComplexAttributes::BaseModelComplexAttributes()
		: BaseModelAttributes() {}
	BaseModelComplexAttributes::BaseModelComplexAttributes(const BaseModelComplexAttributes& copy)
		: BaseModelAttributes()
	{
		for (size_t i = 0; i < parameterMapper.size(); i++)
		{
			*(this->parameterMapper[i]) = *(copy.parameterMapper[i]);
		}
		this->converged = copy.converged;
	}
	BaseModelComplexAttributes::BaseModelComplexAttributes(const ModelParameters& _params)
		: BaseModelAttributes(_params) {}
	BaseModelComplexAttributes::BaseModelComplexAttributes(const BaseModelRealAttributes& realValues)
		: BaseModelAttributes()
	{
		this->delta_cdw = realValues.delta_cdw;
		this->delta_afm = realValues.delta_afm;
		this->delta_sc = realValues.delta_sc;
		this->gamma_sc = realValues.gamma_sc;
		this->xi_sc = std::complex<double>(0., realValues.xi_sc);
		this->delta_eta = std::complex<double>(0., realValues.delta_eta);
		this->gamma_occupation_up = realValues.gamma_occupation_up;
		this->gamma_occupation_down = realValues.gamma_occupation_down;
		this->converged = realValues.converged;
	}

	BaseModelRealAttributes::BaseModelRealAttributes()
		: BaseModelAttributes() {}
	BaseModelRealAttributes::BaseModelRealAttributes(const BaseModelRealAttributes& copy)
		: BaseModelAttributes()
	{
		for (size_t i = 0; i < parameterMapper.size(); i++)
		{
			*(this->parameterMapper[i]) = *(copy.parameterMapper[i]);
		}
		this->converged = copy.converged;
	}
	BaseModelRealAttributes::BaseModelRealAttributes(const ModelParameters& _params)
		: BaseModelAttributes(_params) {}
	BaseModelRealAttributes::BaseModelRealAttributes(const BaseModelComplexAttributes& complexValues)
		: BaseModelAttributes()
	{
		this->delta_cdw = complexValues.delta_cdw.real();
		this->delta_afm = complexValues.delta_afm.real();
		this->delta_sc = complexValues.delta_sc.real();
		this->gamma_sc = complexValues.gamma_sc.real();
		this->xi_sc = complexValues.xi_sc.imag();
		this->delta_eta = complexValues.delta_eta.imag();
		this->gamma_occupation_up = complexValues.gamma_occupation_up.real();
		this->gamma_occupation_down = complexValues.gamma_occupation_down.real();
		this->converged = complexValues.converged;
	}
}