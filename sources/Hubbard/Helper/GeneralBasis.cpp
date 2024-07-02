#include "GeneralBasis.hpp"
#include <chrono>
#include <algorithm>
#include "../MomentumIndexUtility.hpp"

namespace Hubbard::Helper {
	void GeneralBasis::fill_M()
	{
		M = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);

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
		M = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				fillBlock(i, j);
			}
		}
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

		this->starting_states.resize(NUMBER_OF_GREENSFUNCTIONS, VectorCL::Zero(TOTAL_BASIS));
		for (int i = 0; i < Constants::BASIS_SIZE; i++)
		{
			this->starting_states[0](i* number_of_basis_terms) = norm_constant;
			this->starting_states[0](i* number_of_basis_terms + 1) = norm_constant;

			this->starting_states[1](i* number_of_basis_terms) = norm_constant;
			this->starting_states[1](i* number_of_basis_terms + 1) = -norm_constant;

			if (number_of_basis_terms >= 6) {
				this->starting_states[2](i* number_of_basis_terms + 4) = norm_constant;
				this->starting_states[2](i* number_of_basis_terms + 5) = norm_constant;

				this->starting_states[3](i* number_of_basis_terms + 4) = norm_constant;
				this->starting_states[3](i* number_of_basis_terms + 5) = -norm_constant;
			}
		}

		
	}

	void GeneralBasis::printM(int i, int j) const
	{
		if (std::abs(M(j, i)) < DEFAULT_PRECISION) {
			std::cout << 0 << "\t";
		}
		else {
			std::cout << M(j, i).real() << "\t";
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