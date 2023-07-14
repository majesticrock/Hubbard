#pragma once
#include <optional>
#include <array>
#include "../../Utility/InputFileReader.hpp"
#include "../SquareLattice/Model2D.hpp"

namespace Hubbard::Helper {
	typedef std::vector<double> data_vector;
	class PhaseHelper
	{
	private:
		struct Plaquette {
			/* Layout:
			*  0 1
			*  2 3
			*/
			std::array<ModelAttributes<double>, 4> attributes;
			PhaseHelper* parent{};
			
			double lowerFirst{}; 
			double upperFirst{};
			double lowerSecond{};
			double upperSecond{};
			int value_index{};

			inline double& operator()(const size_t i, const size_t j) {
				assert(i < 4);
				return attributes[i][j];
			};
			inline const double& operator()(const size_t i, const size_t j) const {
				assert(i < 4);
				return attributes[i][j];
			};

			inline bool valueIsFinite(const size_t index) const {
				return std::abs(attributes[index][value_index]) > 1e-10;
			};
			inline bool containsPhaseBoundary() const
			{
				for (size_t i = 1U; i < 4U; ++i)
				{
					if (valueIsFinite(0U) != valueIsFinite(i)) {
						return true;
					}
				}
				return false;
			};
			inline double getCenterFirst() const noexcept {
				return 0.5 * (lowerFirst + upperFirst);
			};
			inline double getCenterSecond() const noexcept {
				return 0.5 * (lowerSecond + upperSecond);
			};
			inline double size() const noexcept {
				return 0.5 * (std::abs(upperFirst - lowerFirst) + std::abs(upperSecond - lowerSecond));
			}
			// Devides *this into 4 equally sized plaquettes
			// The result is appended to <appendTo>
			void devidePlaquette(std::vector<Plaquette>& appendTo);
		};

		ModelParameters modelParameters;
		double FIRST_IT_MIN{};
		double FIRST_IT_MAX{};
		double SECOND_IT_MIN{};
		double SECOND_IT_MAX{};

		int FIRST_IT_STEPS{};
		int SECOND_IT_STEPS{};

		int rank{};
		int numberOfRanks{};
		bool use_broyden{};
		
		ModelAttributes<double> computeDataPoint(const ModelParameters& mp, std::optional<ModelAttributes<double>> startingValues = std::nullopt);
	public:
		PhaseHelper(Utility::InputFileReader& input, int _rank, int _nRanks);

		void compute_crude(std::vector<data_vector>& data_mapper);
		void findSingleBoundary(const std::vector<data_vector>& origin, data_vector& recieve_data, int value_index);
	};
}