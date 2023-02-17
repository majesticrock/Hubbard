#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

using Eigen::MatrixXd;

namespace Hubbard {
	HubbardCDW::HubbardCDW(ModelParameters& _params)
		: BasicHubbardModel(_params), V(_params.V)
	{
		this->delta_cdw = abs(U - V) * 0.5;
		this->delta_sc = abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else {
			this->delta_cdw *= 0.25;
		}
		this->delta_eta = 0.01;

		this->hamilton = MatrixXd::Zero(4, 4);
	}
}