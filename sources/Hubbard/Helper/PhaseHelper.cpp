#include "PhaseHelper.hpp"
#include "../SquareLattice/HubbardCDW.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include "../Selfconsistency/BroydenSolver.hpp"
#include "../Selfconsistency/IterativeSolver.hpp"
#include <omp.h>

namespace Hubbard::Helper {
	void PhaseHelper::Plaquette::devidePlaquette(std::vector<PhaseHelper::Plaquette>& appendTo) {
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

		ModelParameters mp{ parent->modelParameters };
		const double centerFirst{ this->getCenterFirst() };
		const double centerSecond{ this->getCenterSecond() };
		ModelAttributes<double> averageParameters{ this->attributes[0] };
		int finiteCount = averageParameters.isOrdered() ? 1 : 0;
		for (const auto& attr : this->attributes)
		{
			if (attr.isOrdered()) {
				++finiteCount;
				averageParameters += attr;
			}
		}
		if (finiteCount != 0) {
			averageParameters /= finiteCount;
		}
		else {
			averageParameters = ModelAttributes<double>(mp);
		}
		std::array<ModelAttributes<double>, 5> new_attributes;

		mp.setGlobalIteratorExact(this->upperFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_attributes[0] = parent->computeDataPoint(mp, averageParameters);

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(this->lowerSecond);
		new_attributes[1] = parent->computeDataPoint(mp, averageParameters);

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_attributes[2] = parent->computeDataPoint(mp, averageParameters);

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(centerFirst);
		mp.setSecondIteratorExact(this->upperSecond);
		new_attributes[3] = parent->computeDataPoint(mp, averageParameters);

		mp = parent->modelParameters;
		mp.setGlobalIteratorExact(this->lowerFirst);
		mp.setSecondIteratorExact(centerSecond);
		new_attributes[4] = parent->computeDataPoint(mp, averageParameters);

		for (const auto& new_attr : new_attributes)
		{
			if (!new_attr.converged) {
				averageParameters.print();
			}
		}

		// Upper left
		Plaquette new_plaq = *this;
		new_plaq.attributes = { this->attributes[0], new_attributes[0], new_attributes[1], new_attributes[2] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.lowerFirst = centerFirst;
			new_plaq.upperSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Upper right
		new_plaq = *this;
		new_plaq.attributes = { new_attributes[0], this->attributes[1], new_attributes[2], new_attributes[3] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.lowerFirst = centerFirst;
			new_plaq.lowerSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Lower left
		new_plaq = *this;
		new_plaq.attributes = { new_attributes[1], new_attributes[2], this->attributes[2], new_attributes[4] };
		if (new_plaq.containsPhaseBoundary()) {
			new_plaq.upperFirst = centerFirst;
			new_plaq.upperSecond = centerSecond;
			appendTo.push_back(new_plaq);
		}

		// Lower right
		new_plaq = *this;
		new_plaq.attributes = { new_attributes[2], new_attributes[3], new_attributes[4], this->attributes[3] };
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

		for (size_t i = 0U; i < option_list.size(); ++i)
		{
			if (input.getString("second_iterator_type") == option_list[i]) {
				SECOND_IT_MIN = model_params[i];
			}
		}
		modelParameters = ModelParameters(model_params[0], model_params[1], model_params[2],
			(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, (SECOND_IT_MAX - SECOND_IT_MIN) / SECOND_IT_STEPS,
			input.getString("global_iterator_type"), input.getString("second_iterator_type"));
	}

	ModelAttributes<double> PhaseHelper::computeDataPoint(const ModelParameters& mp, std::optional<ModelAttributes<double>> startingValues /*= std::nullopt*/) {
		if (startingValues.has_value()) {
			if (use_broyden) {
				SquareLattice::UsingBroyden model(mp, startingValues.value(), 10);
				return model.computePhases();
			}
			SquareLattice::HubbardCDW model(mp, startingValues.value());
			return model.computePhases();
		}
		else {
			if (use_broyden) {
				SquareLattice::UsingBroyden model(mp);
				ModelAttributes<double> result{model.computePhases()};

				if(mp.U < 0 || mp.V < 0) return result;
				// Remember: [0] returns the cdw and [1] the afm gap
				if (std::abs(result[0]) > 1e-12 || std::abs(result[1]) > 1e-12) {
					ModelAttributes<double> copy{ result };
					if (std::abs(result[0]) > 1e-12) {
						copy[1] = result[0];
						copy[0] = 0;
					}
					else {
						copy[0] = result[1];
						copy[1] = 0;
					}

					SquareLattice::UsingBroyden model_copy(mp, copy);
					copy = model_copy.computePhases({false, false});
					if (copy.converged) {
						if (model_copy.freeEnergyPerSite() < model.freeEnergyPerSite()) {
							return copy;
						}
					}
				}

				return result;
			}
			SquareLattice::HubbardCDW model(mp);
			return model.computePhases();
		}
	}

	void PhaseHelper::compute_crude(std::vector<data_vector>& data_mapper) {
		size_t NUMBER_OF_GAP_VALUES = data_mapper.size();
		for (int T = 0; T < FIRST_IT_STEPS; T++)
		{
#pragma omp parallel for num_threads(4) schedule(dynamic)
			for (int U = 0; U < SECOND_IT_STEPS; U++)
			{
				ModelParameters local{ modelParameters };
				local.setSecondIterator(U);
				ModelAttributes<double> ret{ computeDataPoint(local) };

				for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
				{
					data_mapper[i][(T * SECOND_IT_STEPS) + U] = ret[i];
				}
			}
			modelParameters.incrementGlobalIterator();
		}
	}

	void PhaseHelper::findSingleBoundary(const std::vector<data_vector>& origin, data_vector& recieve_data, int value_index, int rank) {
		modelParameters.reset();
		std::vector<Plaquette> plaqs;
		const int rank_offset = FIRST_IT_STEPS * SECOND_IT_STEPS * rank;

		for (size_t i = (rank > 0) ? 0U : 1U; i < FIRST_IT_STEPS; ++i)
		{
			for (size_t j = 1U; j < SECOND_IT_STEPS; ++j)
			{
				Plaquette plaq;
				plaq.value_index = value_index;
				for (size_t l = 0U; l < origin.size(); ++l)
				{
					plaq(0, l) = origin[l][rank_offset + i * SECOND_IT_STEPS + j - 1];
					plaq(1, l) = origin[l][rank_offset + i * SECOND_IT_STEPS + j];
					plaq(2, l) = origin[l][rank_offset + (i - 1) * SECOND_IT_STEPS + j - 1];
					plaq(3, l) = origin[l][rank_offset + (i - 1) * SECOND_IT_STEPS + j];
				}

				if (!plaq.containsPhaseBoundary()) continue;

				plaq.parent = this;
				plaq.lowerFirst = modelParameters.setGlobalIterator(i - 1);
				plaq.upperFirst = modelParameters.setGlobalIterator(i);
				plaq.lowerSecond = modelParameters.setSecondIterator(j - 1);
				plaq.upperSecond = modelParameters.setSecondIterator(j);

				plaqs.push_back(std::move(plaq));
			}
		}

		while (plaqs.size() > 0 && plaqs.begin()->size() > 5e-4) {
			std::cout << "Plaquette size: " << plaqs.begin()->size() << "\t" << "Current number of Plaquettes: " << plaqs.size() << std::endl;
			const auto N_PLAQUETTES = plaqs.size();

			constexpr int n_omp_threads = 8;
			std::vector<std::vector<Plaquette>> buffer(n_omp_threads);
			// omp wants a signed type, so it shall get one
#pragma omp parallel for num_threads(n_omp_threads)
			for (int i = 0; i < N_PLAQUETTES; i++)
			{
				plaqs[i].devidePlaquette(buffer[omp_get_thread_num()]);
			}
			size_t totalSize = 0U;
			for (const auto& vec : buffer)
			{
				totalSize += vec.size();
			}
			plaqs.clear();
			bool removal = false;
			if (totalSize > 100U) {
				totalSize /= 2;
				removal = true;
			}
			plaqs.reserve(totalSize);
			for (const auto& vec : buffer)
			{
				if (removal) {
					for (size_t i = 0U; i < vec.size(); i += 2)
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