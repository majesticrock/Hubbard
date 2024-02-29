#include "UnknownBoundaryHandler.hpp"

#ifdef _DEBUG
#define _NO_MPI
#endif
#ifndef _NO_MPI
#include <mpi.h>
#endif

#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include <filesystem>
#include "Hubbard/Constants.hpp"
#include "Utility/OutputConvenience.hpp"
#include "Hubbard/Helper/DetailModelConstructor.hpp"

using data_vector = std::vector<double>;
constexpr size_t MAX_ITERATIONS = 30U;
constexpr double ERROR_MARGIN = 1e-4;
const std::string BASE_FOLDER = "../../data/phases/";

void UnknownBoundaryHandler::execute(Utility::InputFileReader& input) const {
	Hubbard::Helper::DetailModelConstructorSettings::print_mean_field_result = false;

	double SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");

	for (int i = 0; i < Hubbard::Constants::option_list.size(); i++)
	{
		if (input.getString("second_iterator_type") == Hubbard::Constants::option_list[i]) {
			SECOND_IT_MIN = model_params[i];
		}
	}
	Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, 0.0,
		input.getString("global_iterator_type"), input.getString("second_iterator_type"));

	data_vector local_data(FIRST_IT_STEPS);
	data_vector recieve_data(GLOBAL_IT_STEPS);

	double lower = SECOND_IT_MIN;
	double upper = SECOND_IT_MAX;
	double center = 0.5 * (lower + upper);

	std::chrono::steady_clock::time_point begin, end;
	for (int i = 0; i < FIRST_IT_STEPS; ++i) {
		begin = std::chrono::steady_clock::now();
		std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper_lower, modeHelper_upper;
		lower = SECOND_IT_MIN;
		upper = SECOND_IT_MAX;

		modelParameters.setSecondIteratorExact(lower);
		modeHelper_lower = getHelper(input, modelParameters);
		modelParameters.setSecondIteratorExact(upper);
		modeHelper_upper = getHelper(input, modelParameters);

		{
			bool result_lower = modeHelper_lower->matrix_is_negative();
			bool result_upper = modeHelper_upper->matrix_is_negative();
			if (result_lower == result_upper) {
				// There is no phase transition to be found here
				local_data[i] = std::numeric_limits<double>::quiet_NaN();
				modelParameters.incrementGlobalIterator();
				continue;
			}
			else {
				if (result_lower) {
					// we want the lower value to be the one, where we can actually do the computations
					std::swap(lower, upper);
				}
			}
		}

		for (size_t iter = 0U; iter < MAX_ITERATIONS && abs(upper - lower) > ERROR_MARGIN; ++iter) {
			center = 0.5 * (lower + upper);
			modelParameters.setSecondIteratorExact(center);
			modeHelper_upper = getHelper(input, modelParameters);
			if (modeHelper_upper->matrix_is_negative()) {
				upper = center;
			}
			else {
				lower = center;
			}
			if (iter == MAX_ITERATIONS - 1) std::cerr << "Finished loop #" << i << " with the maximum number of iterations. Precision: " << upper - lower << std::endl;
		}

		local_data[i] = center;
		modelParameters.incrementGlobalIterator();
		end = std::chrono::steady_clock::now();
		std::cout << "Rank #" << rank << ": Runtime for iteration #" << i << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}

#ifndef _NO_MPI
	MPI_Gather(local_data.data(), FIRST_IT_STEPS, MPI_DOUBLE, recieve_data.data(), FIRST_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
	recieve_data = local_data;
#endif
	if (rank == 0) {
		data_vector global_iterator_data(GLOBAL_IT_STEPS);
		const double step = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / GLOBAL_IT_STEPS;
		for (int i = 0; i < GLOBAL_IT_STEPS; ++i) {
			if (std::isnan(recieve_data[i])) {
				global_iterator_data[i] = std::numeric_limits<double>::quiet_NaN();
			}
			else {
				global_iterator_data[i] = GLOBAL_IT_LIMS[0] + i * step;
			}
		}
		// We saved those data points, where there is no afm/cdw phase as nan and we remove them now
		for (auto it = recieve_data.begin(); it != recieve_data.end();) {
			if (std::isnan(*it)) {
				it = recieve_data.erase(it);
			}
			else {
				++it;
			}
		}
		for (auto it = global_iterator_data.begin(); it != global_iterator_data.end();) {
			if (std::isnan(*it)) {
				it = global_iterator_data.erase(it);
			}
			else {
				++it;
			}
		}

		std::string output_folder{ getOutputFolder(input) };
		std::vector<std::string> comments;
		comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
			+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

		comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
			+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

		comments.push_back("Used DOS: " + input.getString("use_DOS"));
		comments.push_back("Discretization: " + input.getString("k_discretization"));
		comments.push_back("Lattice type: " + input.getString("lattice_type"));

		std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
		std::filesystem::create_directories(BASE_FOLDER + output_folder);

		Utility::saveData(global_iterator_data, recieve_data, BASE_FOLDER + output_folder + "unkown_boundary.dat.gz", comments);
	}
}