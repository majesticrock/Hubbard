#pragma once
#include <vector>
#include <string>
#include <optional>

#include "../../Utility/InputFileReader.hpp"
#include "../SquareLattice/Model2D.hpp"

namespace Hubbard::Helper {
	typedef std::vector<double> data_vector;
	class PhaseHelper
	{
	private:
		struct Plaquette {
			PhaseHelper* parent = nullptr;
			int value_index = 0;
			double lowerFirst = 0, upperFirst = 0;
			double lowerSecond = 0, upperSecond = 0;
			/* Layout:
			*  0 1
			*  2 3
			*/
			std::vector<BaseModelRealAttributes> attributes;
			inline double& operator()(int i, int j) {
				return attributes[i][j];
			};

			inline bool valueIsFinite(int index) const {
				return std::abs(attributes[index][value_index]) > 1e-10;
			};
			inline bool containsPhaseBoundary() const
			{
				for (size_t i = 1; i < 4; i++)
				{
					if (valueIsFinite(0) != valueIsFinite(i)) {
						return true;
					}
				}
				return false;
			};
			inline double getCenterFirst() const {
				return 0.5 * (lowerFirst + upperFirst);
			};
			inline double getCenterSecond() const {
				return 0.5 * (lowerSecond + upperSecond);
			};
			inline double size() const {
				return 0.5 * (std::abs(upperFirst - lowerFirst) + std::abs(upperSecond - lowerSecond));
			}
			// Devides *this into 4 equally sized plaquettes
			// The result is appended to <appendTo>
			void devidePlaquette(std::vector<Plaquette>& appendTo);

			Plaquette() : attributes(4) {};
		};

		const std::vector<std::string> option_list = { "T", "U", "V" };
		int FIRST_IT_STEPS;
		double FIRST_IT_MIN, FIRST_IT_MAX;
		int SECOND_IT_STEPS;
		double SECOND_IT_MIN, SECOND_IT_MAX;

		int rank, numberOfRanks;

		bool use_broyden;
		ModelParameters modelParameters;
		BaseModelRealAttributes computeDataPoint(const ModelParameters& mp, std::optional<BaseModelRealAttributes> startingValues = std::nullopt);
	public:
		PhaseHelper(Utility::InputFileReader& input, int _rank, int _nRanks);

		void compute_crude(std::vector<data_vector>& data_mapper);
		void findSingleBoundary(const std::vector<data_vector>& origin, data_vector& recieve_data, int value_index, int rank);
	};
}