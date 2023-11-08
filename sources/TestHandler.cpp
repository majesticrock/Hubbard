#include "TestHandler.hpp"
#include <chrono>
#include <iostream>
#include "Hubbard/ModelParameters.hpp"
#include "Hubbard/SquareLattice/UsingBroyden.hpp"
#include "Hubbard/SquareLattice/HubbardCDW.hpp"
#include "Hubbard/DOSModels/BroydenDOS.hpp"
#include "Hubbard/DensityOfStates/Square.hpp"
#include "Hubbard/DensityOfStates/SimpleCubic.hpp"

using namespace Hubbard::DensityOfStates;
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

void TestHandler::execute(Utility::InputFileReader& input) const
{
	std::chrono::steady_clock::time_point test_b = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point test_e;

	//------------------------------------------------------------//
	ModelAttributes<global_floating_type> startingValues{ 1., 1., 1., 0., 0., 0., 0.1,  0.1 };
	if (input.getBool("use_DOS")) {
		if (input.getString("lattice_type") == "square") {
			DOSModels::BroydenDOS<Square> model(modelParameters, startingValues);
			model.computePhases({ true, true }).print();
			std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;

			Square::writeToBinaryFile("test.bin");
			Square::loadFromBinaryFile("test.bin");
		}
		else if (input.getString("lattice_type") == "cube") {
			DOSModels::BroydenDOS<SimpleCubic> model(modelParameters);
			model.computePhases({ true, true }).print();
			std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;
		}
	}
	else {
		Constants::setDiscretization(input.getInt("k_discretization"));
		Constants::setBasis(4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);

		//SquareLattice::HubbardCDW model(modelParameters, startingValues);
		//model.computePhases({ false, true }).print();
		//std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;
		//std::cout << "\n\n\n";
		SquareLattice::UsingBroyden model2(modelParameters, startingValues);
		test_b = std::chrono::steady_clock::now();
		model2.computePhases({ false, true }).print();
		std::cout << "Free energy = " << model2.freeEnergyPerSite() << std::endl;
	}
	//test_b = std::chrono::steady_clock::now();
	//{
	//	DensityOfStates::SimpleCubic sc_dos;
	//	sc_dos.computeValues();
	//	std::vector<Hubbard::DensityOfStates::dos_precision> abscissa(sc_dos.values.size());
	//	for (size_t i = 0U; i < abscissa.size(); ++i)
	//	{
	//		abscissa[i] = sc_dos.abscissa[i].convert_to<Hubbard::DensityOfStates::dos_precision>();
	//	}
	//	Utility::saveData(abscissa, sc_dos.values, "../../data/3d_dos.dat.gz");
	//}
	//{
	//	DensityOfStates::Square square_dos;
	//	square_dos.computeValues();
	//	std::vector<Hubbard::DensityOfStates::dos_precision> abscissa(square_dos.values.size());
	//	std::cout << square_dos.size() << "   " << square_dos.values.size() << "   " << square_dos.abscissa.size() << std::endl;
	//	for (size_t i = 0U; i < abscissa.size(); ++i)
	//	{
	//		abscissa[i] = square_dos.abscissa[i].convert_to<Hubbard::DensityOfStates::dos_precision>();
	//	}
	//	Utility::saveData(abscissa, square_dos.values, "../../data/2d_dos.dat.gz");
	//}
	//test_e = std::chrono::steady_clock::now();
	//std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
	//std::cout << "\n\n" << std::endl;

	//------------------------------------------------------------//
	test_e = std::chrono::steady_clock::now();
	std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
}