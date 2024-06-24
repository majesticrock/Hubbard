#pragma once
#include <optional>
#include <array>
#include <cassert>
#include <memory>
#include <Utility/InputFileReader.hpp>
#include "../Models/BaseModel.hpp"

namespace Hubbard::Helper {
	typedef std::vector<global_floating_type> data_vector;
	class PhaseHelper
	{
	private:
		friend struct Plaquette;
		Models::ModelParameters modelParameters;
		coefficient_type FIRST_IT_MIN{};
		coefficient_type FIRST_IT_MAX{};
		coefficient_type SECOND_IT_MIN{};
		coefficient_type SECOND_IT_MAX{};

		int FIRST_IT_STEPS{};
		int SECOND_IT_STEPS{};

		int rank{};
		int numberOfRanks{};
		bool use_broyden{};
		unsigned char _internal_lattice_type{};

		std::unique_ptr<Models::BaseModel<global_floating_type>> getModelType(const Models::ModelParameters& mp, std::optional<Models::ModelAttributes<global_floating_type>> startingValues = std::nullopt);

		Models::ModelAttributes<global_floating_type> computeDataPoint(const Models::ModelParameters& mp,
			std::optional<Models::ModelAttributes<global_floating_type>> startingValues = std::nullopt);
		Models::ModelAttributes<global_floating_type> computeDataPoint_No_AFM_CDW_Fix(const Models::ModelParameters& mp,
			std::optional<Models::ModelAttributes<global_floating_type>> startingValues = std::nullopt);
	public:
		PhaseHelper(Utility::InputFileReader& input, int _rank, int _nRanks);

		void compute_crude(std::vector<data_vector>& data_mapper);
		void findSingleBoundary(const std::vector<data_vector>& origin, data_vector& recieve_data, int value_index);
		void coexistence_AFM_CDW(std::vector<data_vector>& recieve_data);
	};
}