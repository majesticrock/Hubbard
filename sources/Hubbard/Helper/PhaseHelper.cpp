#include "PhaseHelper.hpp"
#include "Plaquette.hpp"
#include "../SquareLattice/HubbardCDW.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include "../ChainLattice/ChainTripletPairing.hpp"
#include "../DOSModels/BroydenDOS.hpp"
#include "../DOSModels/PhaseSeparationDOS.hpp"
#include "../DensityOfStates/Square.hpp"
#include "../DensityOfStates/SimpleCubic.hpp"

#include "../Selfconsistency/BroydenSolver.hpp"
#include "../Selfconsistency/IterativeSolver.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <limits>
#include <mutex>
#include <algorithm>

namespace Hubbard::Helper {
	using Constants::option_list;

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
				return std::make_unique<SquareLattice::UsingBroyden>(mp, startingValues.value(), 10U);
			case 3:
				throw std::runtime_error("3D Cube: To be implemented!");

			case 129:
				throw std::runtime_error("1D Chain: To be implemented!");
			case 130:
				return std::make_unique<DOSModels::BroydenDOS<DensityOfStates::Square>>(mp, startingValues.value(), 10U);
			case 131:
				return std::make_unique<DOSModels::BroydenDOS<DensityOfStates::SimpleCubic>>(mp, startingValues.value(), 10U);
			default:
				throw std::runtime_error("_internal_lattice_type not properly set " + to_string(_internal_lattice_type));
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
				throw std::runtime_error("_internal_lattice_type not properly set " + to_string(_internal_lattice_type));
			}
		}
	}

	ModelAttributes<global_floating_type> PhaseHelper::computeDataPoint(const ModelParameters& mp,
		std::optional<ModelAttributes<global_floating_type>> startingValues /*= std::nullopt*/, PhaseDebuggingPolicy debug_messages /*= WarnNoConvergence*/)
	{
		auto model_ptr{ getModelType(mp, startingValues) };
		ModelAttributes<global_floating_type> result{ model_ptr->computePhases(debug_messages) };

		if (mp.U < 0 || mp.V < 0) return result;
		// Remember: [0] returns the cdw and [1] the afm gap
		if (result.isFinite(0) || result.isFinite(1)) {
			ModelAttributes<global_floating_type> copy{ result };
			if (result.isFinite(0)) {
				copy[1] = result[0];
				copy[0] = 0;
			}
			else {
				copy[0] = result[1];
				copy[1] = 0;
			}

			auto model_copy_ptr = getModelType(mp, copy);
			copy = model_copy_ptr->computePhases(NoWarning);
			if (copy.converged) {
				if (model_copy_ptr->freeEnergyPerSite() < model_ptr->freeEnergyPerSite()) {
					return copy;
				}
			}
		}
		return result;
	}

	ModelAttributes<global_floating_type> PhaseHelper::computeDataPoint_No_AFM_CDW_Fix(const ModelParameters& mp,
		std::optional<ModelAttributes<global_floating_type>> startingValues, PhaseDebuggingPolicy debug_messages)
	{
		return getModelType(mp, startingValues)->computePhases(debug_messages);
	}

	void PhaseHelper::compute_crude(std::vector<data_vector>& data_mapper)
	{
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

		while (plaqs.size() > 0U && plaqs.begin()->size() > 5e-5) {
			std::cout << "Plaquette size: " << plaqs.begin()->size() << "\t" << "Current number of Plaquettes: " << plaqs.size() << std::endl;
			const size_t N_PLAQUETTES = plaqs.size();

			const unsigned int n_omp_threads = omp_get_max_threads();
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
			plaqs.reserve(totalSize);
			for (const auto& vec : buffer)
			{
				plaqs.insert(plaqs.end(), vec.begin(), vec.end());
			}

			std::sort(plaqs.begin(), plaqs.end(), [](Plaquette const& a, Plaquette const& b) {
				if (a.lowerFirst < b.lowerFirst) {
					return true;
				}
				else {
					if (a.lowerFirst > b.lowerFirst) {
						return false;
					}
					// If equal sort by second
					return a.lowerSecond < b.lowerSecond;
				}
				});

			if (totalSize > (8 * this->FIRST_IT_STEPS > 200 ? 8 * this->FIRST_IT_STEPS : 200)) {
				// remove every second element
				const int offset = plaqs.size() % 2 == 0 ? 0 : 1;
				for (auto it = plaqs.begin() + offset; it != plaqs.end(); ++it)
				{
					it = plaqs.erase(it);
				}
			}
		}

		recieve_data.reserve(recieve_data.size() + 2 * plaqs.size());
		for (const auto& p : plaqs) {
			recieve_data.push_back(p.getCenterFirst());
			recieve_data.push_back(p.getCenterSecond());
		}
	}

	void PhaseHelper::coexistence_AFM_CDW(std::vector<data_vector>& recieve_data) {
		if (recieve_data.size() < 3U || recieve_data.at(0).size() < FIRST_IT_STEPS) {
			std::cerr << "You are calling coexistence_AFM_CDW() with an empty data reciever!" << std::endl;
			recieve_data.resize(3U, data_vector(FIRST_IT_STEPS));
		}

		auto base_U = [this](int it) -> double {
			return FIRST_IT_MIN + ((FIRST_IT_MAX - FIRST_IT_MIN) * it) / FIRST_IT_STEPS;
			};
		auto base_V = [](double U) -> double {
			return 0.25 * U;
			};

		constexpr double PRECISION = 1e-6;
		constexpr double EPSILON = 1e-10;
		for (int i = 0; i < FIRST_IT_STEPS; ++i) {
			ModelParameters local{ modelParameters };
			const double U = base_U(i);
			recieve_data[0][i] = U;
			local.setGlobalIteratorExact(U);
			local.setSecondIteratorExact(base_V(U));

			ModelAttributes<global_floating_type> base_gap{ computeDataPoint(local) };
			base_gap[0] = ONE_OVER_SQRT_2 * sqrt(base_gap[0] * base_gap[0] + base_gap[1] * base_gap[1]);
			base_gap[1] = 0;

			// If there is no cdw nor afm order, we save NaN to the array and remove it later
			if (std::abs(base_gap[0]) < EPSILON) {
				recieve_data[0][i] = std::numeric_limits<double>::quiet_NaN();
				recieve_data[1][i] = std::numeric_limits<double>::quiet_NaN();
				recieve_data[2][i] = std::numeric_limits<double>::quiet_NaN();
				continue;
			}

			double h{ 0.5 };
			double a{ base_V(U) };
			double b{ a - h };
			ModelAttributes<global_floating_type> gap_b{ base_gap };

			while (a - b > PRECISION) {
				local.setSecondIteratorExact(b);
				gap_b = computeDataPoint_No_AFM_CDW_Fix(local, base_gap, NoWarning);
				if (std::abs(gap_b[0]) > EPSILON) {
					a = b;
					b -= h;
				}
				else {
					b = 0.5 * (b + a);
					h *= 0.5;
				}
			}
			recieve_data[1][i] = 0.5 * (a + b);

			base_gap[1] = base_gap[0];
			base_gap[0] = 0;
			h = 0.5;
			a = base_V(U);
			b = a + h;
			while (b - a > PRECISION) {
				local.setSecondIteratorExact(b);
				gap_b = computeDataPoint_No_AFM_CDW_Fix(local, base_gap, NoWarning);
				if (std::abs(gap_b[1]) > EPSILON) {
					a = b;
					b += h;
				}
				else {
					b = 0.5 * (b + a);
					h *= 0.5;
				}
			}
			recieve_data[2][i] = 0.5 * (a + b);
		}
	}
}