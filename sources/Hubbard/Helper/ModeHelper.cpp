#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	void ModeHelper::loadWick(const std::string& filename)
	{
		wicks_M.resize(number_of_basis_terms * number_of_basis_terms);
		wicks_N.resize(number_of_basis_terms * number_of_basis_terms);
		const int name_offset = (start_basis_at < 0) ? 0 : start_basis_at;

		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "M_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_M[j * number_of_basis_terms + i].clear();
					ia >> wicks_M[j * number_of_basis_terms + i];
					ifs.close();
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "N_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_N[j * number_of_basis_terms + i].clear();
					ia >> wicks_N[j * number_of_basis_terms + i];
					ifs.close();
				}
			}
		}
	}

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
	{
		start_basis_at = input.getInt("start_basis_at");
		this->number_of_basis_terms = input.getInt("number_of_basis_terms");
		if (start_basis_at < 0) {
			// We investigate the special x-p-basis
			this->TOTAL_BASIS = Constants::BASIS_SIZE * 8;
			this->number_of_basis_terms = 10;
		}
		else {
			this->TOTAL_BASIS = Constants::BASIS_SIZE * this->number_of_basis_terms;
		}

		loadWick("../commutators/wick_");
	}
}