#include "TestHandler.hpp"
#include <chrono>
#include <iostream>
#include "Hubbard/ModelParameters.hpp"
#include "Hubbard/SquareLattice/UsingBroyden.hpp"
#include "Hubbard/SquareLattice/HubbardCDW.hpp"
#include "Hubbard/DOSModels/BroydenDOS.hpp"
#include "Hubbard/DensityOfStates/Square.hpp"
#include "Hubbard/DensityOfStates/SimpleCubic.hpp"
#include "Utility/GramSchmidt.hpp"
#include "Hubbard/DOSModels/PhaseSeparationDOS.hpp"
#include "Hubbard/EMCoupling.hpp"

using namespace Hubbard;

std::ostream& operator<<(std::ostream& os, const std::vector<double>& data) {
	for (size_t i = 0U; i < data.size(); ++i)
	{
		os << data[i];
		if (i != data.size() - 1) {
			os << ", ";
		}
	}
	return os;
}

#include <Eigen/Sparse>
#include "Utility/PivotToBlockStructure.hpp"
using namespace Eigen;

void TestHandler::execute(Utility::InputFileReader& input) const
{
	std::chrono::steady_clock::time_point test_b = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point test_e;

	if (input.getBool("em_coupling")) {
		Constants::setDiscretization(input.getInt("k_discretization"));
		Constants::setBasis(4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);

		EMCoupling model(modelParameters);
		model.computePhases().print();
	}
	else {
		//ModelAttributes<global_floating_type> startingValues{ 1., 1., 1., 0., 0., 0., 0.1,  0.1, 0. };
		if (input.getBool("use_DOS")) {
			if (input.getString("lattice_type") == "square") {
				DOSModels::BroydenDOS<DensityOfStates::Square> model(modelParameters);
				model.computePhases(WarnNoConvergence).print();
				std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;
				std::cout << "Sum rule: " << model.cdw_in_sc_sum_rule() << std::endl;
			}
			else if (input.getString("lattice_type") == "cube") {
				DOSModels::PhaseSeparationDOS<DensityOfStates::SimpleCubic> model(modelParameters, 1);
				model.computePhases(WarnNoConvergence).print();
				std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;
				std::cout << "Sum rule: " << model.cdw_in_sc_sum_rule() << std::endl;
			}
		}
		else {
			Constants::setDiscretization(input.getInt("k_discretization"));
			Constants::setBasis(4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);

			SquareLattice::UsingBroyden model2(modelParameters);
			test_b = std::chrono::steady_clock::now();
			model2.computePhases({ false, true }).print();
			std::cout << "Free energy = " << model2.freeEnergyPerSite() << std::endl;

			std::cout << "Sum rule: " << model2.cdw_in_sc_sum_rule() << std::endl;
		}
	}

	//------------------------------------------------------------//
	test_e = std::chrono::steady_clock::now();
	std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
}