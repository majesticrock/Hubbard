#include "GeneralBasis.hpp"
#include "../MomentumIndexUtility.hpp"

namespace Hubbard::Helper {
	void GeneralBasis::fill_M()
	{
		M = decltype(M)::Zero(TOTAL_BASIS, TOTAL_BASIS);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				fill_block_M(i, j);
			}
		}
	}

	void GeneralBasis::fillMatrices()
	{
		M.setZero(TOTAL_BASIS, TOTAL_BASIS);
		N.setZero(TOTAL_BASIS, TOTAL_BASIS);

		for (int i = 0; i < std::min(9, number_of_basis_terms); ++i)
		{
			for (int j = 0; j < std::min(9, number_of_basis_terms); ++j)
			{
				fillBlock(i, j);
			}
		}
		// There are no matrix elements connecting the transversal magnons with the rest of the system.
		// Therefore, we may seperate them here.
		for (int i = 9; i < number_of_basis_terms; ++i)
		{
			for (int j = 9; j < number_of_basis_terms; ++j)
			{
				fillBlock(i, j);
			}
		}

		std::cout << "||M - M^+|| = " << (M - M.adjoint()).norm() << std::endl;
		std::cout << "||N - N^+|| = " << (N - N.adjoint()).norm() << std::endl;
	}

	void GeneralBasis::createStartingStates()
	{
		/*
		* 0 - SC Higgs
		* 1 - SC Phase
		* 2 - CDW
		* 3 - AFM
		*/
		const global_floating_type norm_constant = this->usingDOS
			? sqrt((2.0 * this->dos_dimension) / Constants::BASIS_SIZE)
			: sqrt(1. / ((global_floating_type)Constants::BASIS_SIZE));

		int n_starting_states = 2; // SC Higgs and Phase
		if(number_of_basis_terms >= 6) n_starting_states += 2; // AFM longitudinal and CDW
		if(number_of_basis_terms >= 10) ++n_starting_states; // AFM transversal
		this->starting_states.resize(n_starting_states, Vector_L::Zero(TOTAL_BASIS));

		for (int k = 0; k < Constants::BASIS_SIZE; k++)
		{
			this->starting_states[0](get_index(0, k)) = norm_constant;
			this->starting_states[0](get_index(1, k)) = norm_constant;

			this->starting_states[1](get_index(0, k)) = norm_constant;
			this->starting_states[1](get_index(1, k)) = -norm_constant;

			if (number_of_basis_terms >= 6) {
				this->starting_states[2](get_index(4, k)) = norm_constant;
				this->starting_states[2](get_index(5, k)) = norm_constant;

				this->starting_states[3](get_index(4, k)) = norm_constant;
				this->starting_states[3](get_index(5, k)) = -norm_constant;
			}
			if(number_of_basis_terms >= 10) {
				this->starting_states[4](get_index(8, k)) = norm_constant;
				this->starting_states[4](get_index(9, k)) = norm_constant;
			}
		}

		this->resolvent_names.reserve(n_starting_states);
		this->resolvent_names.push_back("amplitude_SC");
		this->resolvent_names.push_back("phase_SC");
		if (number_of_basis_terms >= 6) {
			this->resolvent_names.push_back("amplitude_CDW");
			this->resolvent_names.push_back("amplitude_AFM");
		}
		if(number_of_basis_terms >= 10) {
			this->resolvent_names.push_back("amplitude_AFM_transversal");
		}
	}

	void GeneralBasis::printM(int i, int j) const
	{
		if (std::abs(M(j, i)) < DEFAULT_PRECISION) {
			std::cout << 0 << "\t";
		}
		else {
			std::cout << M(j, i) << "\t";
		}
	}

	void GeneralBasis::printMomentumBlocks() const
	{
		std::cout << std::fixed << std::setprecision(4) << std::endl;

		for (int l = 0; l < Constants::K_DISCRETIZATION; l++)
		{
			for (int k = 0; k < Constants::K_DISCRETIZATION; ++k)
			{
				int idx = l * Constants::K_DISCRETIZATION + k;
				int jdx = addQTo(idx);
				std::cout << idx << ": " << gammaFromIndex(idx) << "\n";

				for (size_t i = 0U; i < number_of_basis_terms; ++i)
				{
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + idx * number_of_basis_terms, i + idx * number_of_basis_terms);
					}
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + jdx * number_of_basis_terms, i + idx * number_of_basis_terms);
					}
					std::cout << std::endl;
				}
				for (size_t i = 0U; i < number_of_basis_terms; ++i)
				{
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + idx * number_of_basis_terms, i + jdx * number_of_basis_terms);
					}
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + jdx * number_of_basis_terms, i + jdx * number_of_basis_terms);
					}
					std::cout << std::endl;
				}

				std::cout << std::endl;
			}
		}
	}

	void GeneralBasis::printDOSBlocks() const
	{
		std::cout << std::fixed << std::setprecision(4) << std::endl;

		for (size_t k = 0U; k < Constants::BASIS_SIZE / 2; ++k)
		{
			for (size_t i = 0U; i < number_of_basis_terms; ++i)
			{
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + k * number_of_basis_terms, i + k * number_of_basis_terms);
				}
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms, i + k * number_of_basis_terms);
				}
				std::cout << std::endl;
			}
			for (size_t i = 0U; i < number_of_basis_terms; ++i)
			{
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + k * number_of_basis_terms, i + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms);
				}
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms, i + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms);
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}

	bool GeneralBasis::matrix_is_negative() {
		return _parent_algorithm::dynamic_matrix_is_negative();
	}

	std::vector<ResolventReturnData> GeneralBasis::computeCollectiveModes() {
		return _parent_algorithm::computeCollectiveModes(this->usingDOS ? 150 : 2 * Constants::K_DISCRETIZATION);
	}
}