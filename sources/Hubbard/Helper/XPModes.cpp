#include "XPModes.hpp"
#include <chrono>

namespace Hubbard::Helper {
	void XPModes::fill_M()
	{
		constexpr int a = hermitian_size - 1;
		constexpr int b = antihermitian_size - 1;

		K_plus.setZero(a * Constants::BASIS_SIZE, a * Constants::BASIS_SIZE);
		K_minus.setZero(b * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
			}
		}
	}

	void XPModes::fillMatrices()
	{
		constexpr int a = hermitian_size - 1;
		constexpr int b = antihermitian_size - 1;

		K_plus.setZero(a * Constants::BASIS_SIZE, a * Constants::BASIS_SIZE);
		K_minus.setZero(b * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);
		L.setZero(a * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
				// N only contains offdiagonal blocks
				else if (i < hermitian_size && j >= hermitian_size) {
					fill_block_N(i, j);
				}
			}
		}

		if ((K_plus - K_plus.adjoint()).norm() > ERROR_MARGIN * K_plus.rows() * K_plus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_plus - K_plus.adjoint()).norm()));
		if ((K_minus - K_minus.adjoint()).norm() > ERROR_MARGIN * K_minus.rows() * K_minus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_minus - K_minus.adjoint()).norm()));
	}

	void XPModes::createStartingStates()
	{
		this->starting_states.reserve(4);
		this->starting_states.push_back({ Vector_L::Zero(K_minus.rows()), Vector_L::Zero(K_plus.rows()), "SC"}); // SC
		this->starting_states.push_back( _parent_algorithm::OnlyAmplitude(K_plus.rows(), "CDW") ); // CDW
		this->starting_states.push_back( _parent_algorithm::OnlyAmplitude(K_plus.rows(), "AFM") ); // AFM
		this->starting_states.push_back( _parent_algorithm::OnlyAmplitude(K_plus.rows(), "AFM_transversal") ); // AFM transversal

		const global_floating_type norm_constant =
#ifdef _EXACT_DOS
			sqrt(1. / ((global_floating_type)Constants::BASIS_SIZE));
#else
			this->usingDOS ? sqrt((2.0 * this->dos_dimension) / Constants::BASIS_SIZE) : sqrt(1. / ((global_floating_type)Constants::BASIS_SIZE));
#endif
		for (int i = 0; i < Constants::BASIS_SIZE; ++i)
		{
			starting_states[0][0](i) = norm_constant;
			starting_states[0][1](i) = norm_constant;
			starting_states[1][1](2 * Constants::BASIS_SIZE + i) = norm_constant;
			starting_states[2][1](2 * Constants::BASIS_SIZE + i) = (i < Constants::BASIS_SIZE / 2) ? norm_constant : -norm_constant;
			starting_states[3][1](3 * Constants::BASIS_SIZE + i) = norm_constant;
		}
	}

	bool XPModes::matrix_is_negative() {
		return _parent_algorithm::dynamic_matrix_is_negative();
	};

	std::vector<ResolventReturnData> XPModes::computeCollectiveModes()
	{
		return _parent_algorithm::computeCollectiveModes(this->usingDOS ? 150 : 2 * Constants::K_DISCRETIZATION);
	}
}