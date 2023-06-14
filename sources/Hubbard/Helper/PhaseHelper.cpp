#include "PhaseHelper.hpp"
#include "../UsingBroyden.hpp"
#include "../HubbardCDW.hpp"

namespace Hubbard::Helper {
	void PhaseHelper::Plaquette::devidePlaquette(std::vector<PhaseHelper::Plaquette>& appendTo, int value_index) {
		/*
		*			x---A---x
		*			|   |   |
		*			| 0 | 1 |
		*			|   |   |
		*			B---C---D
		*			|   |   |
		*			| 2 | 3 |
		*			|   |   |
		*	  . 	x---E---x
		*	 /|\
		*     |
		*   First / Second ->
		*
		*	New values are at the position ABCDE (01234 as indizes)
		*/

		const double centerFirst = this->getCenterFirst();
		const double centerSecond = this->getCenterSecond();
		double new_values[5];

		auto mp = parent->modelParameters;
		mp.setGlobalIteratorExact(this->upperFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_values[0] = parent->computeDataPoint(mp)[value_index];

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(this->lowerSecond);
		new_values[1] = parent->computeDataPoint(mp)[value_index];

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_values[2] = parent->computeDataPoint(mp)[value_index];

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(upperSecond);
		new_values[3] = parent->computeDataPoint(mp)[value_index];

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(this->lowerFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_values[4] = parent->computeDataPoint(mp)[value_index];

		// Upper left
		Plaquette new_plaq = *this;
		new_plaq.values = { this->values[0], new_values[0], new_values[1], new_values[2] };

		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.lowerFirst = centerFirst;
			new_plaq.upperSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Upper right
		new_plaq = *this;
		new_plaq.values = { new_values[0], this->values[1], new_values[2], new_values[3] };

		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.lowerFirst = centerFirst;
			new_plaq.lowerSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Lower left
		new_plaq = *this;
		new_plaq.values = { new_values[1], new_values[2], this->values[2], new_values[4] };

		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.upperFirst = centerFirst;
			new_plaq.upperSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Lower right
		new_plaq = *this;
		new_plaq.values = { new_values[2], new_values[3], new_values[4], this->values[3] };

		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.upperFirst = centerFirst;
			new_plaq.lowerSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}
	}

	PhaseHelper::PhaseHelper(Utility::InputFileReader& input, int _rank, int _nRanks)
		: rank(_rank), numberOfRanks(_nRanks)
	{
		// Setup the parameters T, U, V
		std::vector<double> model_params = input.getDoubleList("model_parameters");
		// Setup the number of steps
		int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
		FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
		double GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
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
		modelParameters = Model::ModelParameters(model_params[0], model_params[1], model_params[2],
			(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, (SECOND_IT_MAX - SECOND_IT_MIN) / SECOND_IT_STEPS,
			input.getString("global_iterator_type"), input.getString("second_iterator_type"));
	}

	Model::data_set PhaseHelper::computeDataPoint(const Model::ModelParameters& mp) {
		Model::data_set ret;
		if (use_broyden) {
			UsingBroyden model(mp, 0, 0);
			return model.computePhases();
		}
		HubbardCDW model(mp, 0, 0);
		return model.computePhases();
	}

	void PhaseHelper::compute_crude(std::vector<data_vector*>& data_mapper) {
		int NUMBER_OF_GAP_VALUES = data_mapper.size();
		for (int T = 0; T < FIRST_IT_STEPS; T++)
		{
#pragma omp parallel for num_threads(4) schedule(dynamic)
			for (int U = 0; U < SECOND_IT_STEPS; U++)
			{
				Model::ModelParameters local = modelParameters;
				local.setSecondIterator(U);
				Model::data_set ret = computeDataPoint(local);

				for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
				{
					(*data_mapper[i])[(T * SECOND_IT_STEPS) + U] = ret[i];
				}
			}
			modelParameters.incrementGlobalIterator();
		}
	}

	template <typename T>
	std::vector<T> joinVectors(const std::vector<std::vector<T>>& vectors)
	{
		std::vector<T> result;

		for (const auto& vec : vectors)
		{
			result.insert(result.end(), vec.begin(), vec.end());
		}

		return result;
	}

	void PhaseHelper::findSingleBoundary(const data_vector& origin, data_vector& recieve_data, int value_index, int rank) {
		modelParameters.reset();
		std::vector<Plaquette> plaqs;
		const int rank_offset = FIRST_IT_STEPS * SECOND_IT_STEPS * rank;

		for (size_t i = (rank > 0) ? 0 : 1; i < FIRST_IT_STEPS; i++)
		{
			for (size_t j = 1; j < SECOND_IT_STEPS; j++)
			{
				Plaquette plaq;
				plaq.values = {
					origin[rank_offset + i * SECOND_IT_STEPS + j - 1], origin[rank_offset + i * SECOND_IT_STEPS + j],
					origin[rank_offset + (i - 1) * SECOND_IT_STEPS + j - 1], origin[rank_offset + (i - 1) * SECOND_IT_STEPS + j]
				};
				if (!plaq.containsPhaseBoundary()) continue;

				plaq.parent = this;
				plaq.lowerFirst = modelParameters.setGlobalIterator(i - 1);
				plaq.upperFirst = modelParameters.setGlobalIterator(i);
				plaq.lowerSecond = modelParameters.setSecondIterator(j - 1);
				plaq.upperSecond = modelParameters.setSecondIterator(j);
				plaqs.push_back(plaq);
			}
		}

		while (plaqs.size() > 0 && plaqs.begin()->size() > 5e-4) {
			std::cout << "Plaquette size: " << plaqs.begin()->size() << "\t" << "Current number of Plaquettes: " << plaqs.size() << std::endl;
			const auto N_PLAQUETTES = plaqs.size();

			const int n_omp_threads = 4;
			std::vector<std::vector<Plaquette>> buffer(n_omp_threads);
#pragma omp parallel for num_threads(n_omp_threads)
			for (size_t i = 0; i < N_PLAQUETTES; i++)
			{
				plaqs[i].devidePlaquette(buffer[omp_get_thread_num()], value_index);
			}
			int totalSize = 0;
			for (const auto& vec : buffer)
			{
				totalSize += vec.size();
			}
			plaqs.clear();
			bool removal = false;
			if (totalSize > 200) {
				totalSize /= 2;
				removal = true;
			}
			plaqs.reserve(totalSize);
			for (const auto& vec : buffer)
			{
				if (removal) {
					for (size_t i = 0; i < vec.size(); i += 2)
					{
						plaqs.push_back(vec[i]);
					}
				}
				else {
					plaqs.insert(plaqs.end(), vec.begin(), vec.end());
				}
			}
		}

		recieve_data.reserve(recieve_data.size() + 2 * plaqs.size());
		for (const auto& p : plaqs) {
			recieve_data.push_back(p.getCenterFirst());
			recieve_data.push_back(p.getCenterSecond());
		}
	}
}