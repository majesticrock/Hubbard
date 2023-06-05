#include "PhaseHelper.hpp"
#include "Hubbard/UsingBroyden.hpp"
#include "Hubbard/HubbardCDW.hpp"

PhaseHelper::PhaseHelper(Utility::InputFileReader& input, int _rank, int _nRanks)
	: rank(_rank), numberOfRanks(_nRanks)
{
	// Setup the parameters T, U, V
	std::vector<double> model_params = input.getDoubleList("model_parameters");
	// Setup the number of steps
	int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
	FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
	int GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
	const std::vector<std::string> option_list = { "T", "U", "V" };
	double FIRST_IT_RANGE = 0;
	FIRST_IT_MIN = 0, FIRST_IT_MAX = 0;
	for (int i = 0; i < option_list.size(); i++)
	{
		if (input.getString("global_iterator_type") == option_list[i]) {
			GLOBAL_IT_LIMS[0] = model_params[i];
			FIRST_IT_RANGE = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / numberOfRanks;
			FIRST_IT_MIN = GLOBAL_IT_LIMS[0] + rank * FIRST_IT_RANGE;
			FIRST_IT_MAX = FIRST_IT_MIN + FIRST_IT_RANGE;
			model_params[i] = FIRST_IT_MIN;
		}
	}

	use_broyden = input.getBool("use_broyden");
	SECOND_IT_STEPS = input.getInt("second_iterator_steps");
	SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");

	for (int i = 0; i < option_list.size(); i++)
	{
		if (input.getString("second_iterator_type") == option_list[i]) {
			SECOND_IT_MIN = model_params[i];
		}
	}
	modelParameters = Hubbard::Model::ModelParameters(model_params[0], model_params[1], model_params[2],
		(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, (SECOND_IT_MAX - SECOND_IT_MIN) / SECOND_IT_STEPS,
		input.getString("global_iterator_type"), input.getString("second_iterator_type"));
}

void PhaseHelper::compute_crude(std::vector<data_vector*>& data_mapper) {
	int NUMBER_OF_GAP_VALUES = data_mapper.size();
	for (int T = 0; T < FIRST_IT_STEPS; T++)
	{
#pragma omp parallel for num_threads(4) schedule(dynamic)
		for (int U = 0; U < SECOND_IT_STEPS; U++)
		{
			Hubbard::Model::ModelParameters local = modelParameters;
			local.setSecondIterator(U);
			Hubbard::Model::data_set ret;
			if (use_broyden) {
				Hubbard::UsingBroyden model(local, 0, 0);
				ret = model.computePhases();
			}
			else {
				Hubbard::HubbardCDW model(local, 0, 0);
				ret = model.computePhases();
			}

			for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
			{
				(*data_mapper[i])[(T * SECOND_IT_STEPS) + U] = ret[i];
			}
		}
		modelParameters.incrementGlobalIterator();
	}
}

void PhaseHelper::findSingleBoundary(const data_vector& origin, data_vector& recieve_data) {
	modelParameters.reset();
	std::vector<Plaquette> plaqs;

	for (size_t i = 1; i < FIRST_IT_STEPS; i++)
	{
		for (size_t j = 1; j < SECOND_IT_STEPS; j++)
		{
			Plaquette plaq;
			plaq.values = {
				origin[i * SECOND_IT_STEPS + j], origin[i * SECOND_IT_STEPS + j - 1],
				origin[(i - 1) * SECOND_IT_STEPS + j], origin[(i - 1) * SECOND_IT_STEPS + j - 1]
			};
			if (!plaq.containsPhaseBoundary()) continue;

			plaq.lowerFirst = modelParameters.setGlobalIterator(i - 1);
			plaq.upperFirst = modelParameters.setGlobalIterator(i);
			plaq.lowerSecond = modelParameters.setSecondIterator(j - 1);
			plaq.upperSecond = modelParameters.setSecondIterator(j);
			plaqs.push_back(plaq);
		}
	}


}

void PhaseHelper::findBoundaries(const std::vector<data_vector*>& data_mapper, data_vector& recieve_data)
{
	int NUMBER_OF_GAP_VALUES = data_mapper.size();
	recieve_data.reserve(2 * NUMBER_OF_GAP_VALUES * FIRST_IT_STEPS);
	for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
	{
		findSingleBoundary(*(data_mapper[i]), recieve_data);
	}
}