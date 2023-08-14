#include "PhaseHelper.hpp"
#include "../SquareLattice/HubbardCDW.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include "../ChainLattice/ChainTripletPairing.hpp"
#include "../DOSModels/BroydenDOS.hpp"
#include "../DensityOfStates/Square.hpp"
#include "../DensityOfStates/SimpleCubic.hpp"

#include "../Selfconsistency/BroydenSolver.hpp"
#include "../Selfconsistency/IterativeSolver.hpp"
#include "../Constants.hpp"
#include <omp.h>

namespace Hubbard::Helper {
	using Constants::option_list;

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
		ModelAttributes<global_floating_type> averageParameters{ this->attributes[0] };
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
			averageParameters = ModelAttributes<global_floating_type>(mp);
		}
		std::array<ModelAttributes<global_floating_type>, 5> new_attributes;

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
		: rank(_rank), numberOfRanks(_nRanks), use_broyden(input.getBool("use_broyden"))
	{
		std::string lattice_type = input.getString("lattice_type");
		if (lattice_type == "chain") {
			_internal_lattice_type = 1;
		}
		else if (lattice_type == "square") {
			_internal_lattice_type = 2;
		}
		else if (lattice_type == "cube") {
			_internal_lattice_type = 3;
		}
		else {
			throw std::invalid_argument("Did not recognise lattice type: " + lattice_type);
		}

		if (input.getBool("use_DOS")) {
			_internal_lattice_type += 128;
		}

		// Setup the parameters T, U, V
		std::vector<double> model_params = input.getDoubleList("model_parameters");
		// Setup the number of steps
		int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
		FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
		double GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
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

	std::unique_ptr<BaseModel<global_floating_type>> PhaseHelper::getModelType(const ModelParameters& mp,
		std::optional<ModelAttributes<global_floating_type>> startingValues /*= std::nullopt*/)
	{
		if (!use_broyden) {
			throw std::invalid_argument("Not using broyden's method but requiring real returns.");
		}
		if (startingValues.has_value()) {
			switch (_internal_lattice_type) {
			case 1:
				throw std::runtime_error("1D Chain: To be implemented!");
			case 2:
				return std::make_unique<SquareLattice::UsingBroyden>(mp, startingValues.value(), 10);
			case 3:
				throw std::runtime_error("3D Cube: To be implemented!");

			case 129:
				throw std::runtime_error("1D Chain: To be implemented!");
			case 130:
				return std::make_unique<DOSModels::BroydenDOS<DensityOfStates::Square>>(mp, startingValues.value(), 10);
			case 131:
				return std::make_unique<DOSModels::BroydenDOS<DensityOfStates::SimpleCubic>>(mp, startingValues.value(), 10);
			default:
				throw std::runtime_error("_internal_lattice_type not properly set " + std::to_string(_internal_lattice_type));
			}
		}
		else {
			switch (_internal_lattice_type) {
			case 1:
				throw std::runtime_error("1D Chain: To be implemented!");
			case 2:
				return std::make_unique<SquareLattice::UsingBroyden>(mp);
			case 3:
				throw std::runtime_error("3D Cube: To be implemented!");

			case 129:
				throw std::runtime_error("1D Chain: To be implemented!");
			case 130:
				return std::make_unique<DOSModels::BroydenDOS<DensityOfStates::Square>>(mp);
			case 131:
				return std::make_unique<DOSModels::BroydenDOS<DensityOfStates::SimpleCubic>>(mp);
			default:
				throw std::runtime_error("_internal_lattice_type not properly set " + std::to_string(_internal_lattice_type));
			}
		}
	}

	ModelAttributes<global_floating_type> PhaseHelper::computeDataPoint(const ModelParameters& mp, std::optional<ModelAttributes<global_floating_type>> startingValues /*= std::nullopt*/) {
		if (startingValues.has_value()) {
			return getModelType(mp, startingValues)->computePhases();
		}
		else {
			auto model_ptr{ getModelType(mp) };
			ModelAttributes<global_floating_type> result{model_ptr->computePhases()};

			if (mp.U < 0 || mp.V < 0) return result;
			// Remember: [0] returns the cdw and [1] the afm gap
			if (std::abs(result[0]) > 1e-12 || std::abs(result[1]) > 1e-12) {
				ModelAttributes<global_floating_type> copy{ result };
				if (std::abs(result[0]) > 1e-12) {
					copy[1] = result[0];
					copy[0] = 0;
				}
				else {
					copy[0] = result[1];
					copy[1] = 0;
				}

				auto model_copy_ptr = getModelType(mp, copy);
				copy = model_copy_ptr->computePhases({ false, false });
				if (copy.converged) {
					if (model_copy_ptr->freeEnergyPerSite() < model_ptr->freeEnergyPerSite()) {
						return copy;
					}
				}
			}

			return result;
		}
	}

	void PhaseHelper::compute_crude(std::vector<data_vector>& data_mapper) {
		size_t NUMBER_OF_GAP_VALUES = data_mapper.size();
		for (int T = 0; T < FIRST_IT_STEPS; T++)
		{
#pragma omp parallel for schedule(dynamic)
			for (int U = 0; U < SECOND_IT_STEPS; U++)
			{
				ModelParameters local{ modelParameters };
				local.setSecondIterator(U);
				ModelAttributes<global_floating_type> ret{ computeDataPoint(local) };

				for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
				{
					data_mapper[i][(T * SECOND_IT_STEPS) + U] = ret[i];
				}
			}
			modelParameters.incrementGlobalIterator();
		}
	}

	void PhaseHelper::findSingleBoundary(const std::vector<data_vector>& origin, data_vector& recieve_data, int value_index) {
		modelParameters.reset();
		std::vector<Plaquette> plaqs;
		const int rank_offset = FIRST_IT_STEPS * SECOND_IT_STEPS * rank;

		for (int i = (rank > 0) ? 0 : 1; i < FIRST_IT_STEPS; ++i)
		{
			for (int j = 1; j < SECOND_IT_STEPS; ++j)
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

			const int n_omp_threads = omp_get_num_threads();
			std::vector<std::vector<Plaquette>> buffer(n_omp_threads);
			// omp wants a signed type, so it shall get one
#pragma omp parallel for
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
			if (totalSize > 80U) {
				totalSize /= 2;
				removal = true;
			}
			plaqs.reserve(totalSize);
			for (const auto& vec : buffer)
			{
				if (removal) {
					for (size_t i = 0U; i < vec.size(); i += 2U)
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