#include "TestHandler.hpp"
#include <chrono>
#include <iostream>
#include "Hubbard/ModelParameters.hpp"
#include "Hubbard/SquareLattice/UsingBroyden.hpp"
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
	/*if (input.getString("lattice_type") == "square") {
		DOSModels::BroydenDOS<Square> model3(modelParameters);
		model3.computePhases({ false, true }).print();
		std::cout << "Free energy = " << model3.freeEnergyPerSite() << std::endl;
	}
	else if (input.getString("lattice_type") == "cube") {
		DOSModels::BroydenDOS<SimpleCubic> model3(modelParameters);
		model3.computePhases({ false, true }).print();
		std::cout << "Free energy = " << model3.freeEnergyPerSite() << std::endl;
	}*/
	Constants::setDiscretization(input.getInt("k_discretization"));
	//test_b = std::chrono::steady_clock::now();
	{
		DensityOfStates::SimpleCubic sc_dos;
		sc_dos.computeValues();
		std::vector<Hubbard::DensityOfStates::dos_precision> abscissa(sc_dos.values.size());
		for (size_t i = 0U; i < abscissa.size(); ++i)
		{
			abscissa[i] = sc_dos.abscissa[i].convert_to<Hubbard::DensityOfStates::dos_precision>();
		}
		Utility::saveData(abscissa, sc_dos.values, "../../data/3d_dos.dat.gz");
	}
	{
		DensityOfStates::Square square_dos;
		square_dos.computeValues();
		std::vector<Hubbard::DensityOfStates::dos_precision> abscissa(square_dos.values.size());
		std::cout << square_dos.size() << "   " << square_dos.values.size() << "   " << square_dos.abscissa.size() << std::endl;
		for (size_t i = 0U; i < abscissa.size(); ++i)
		{
			abscissa[i] = square_dos.abscissa[i].convert_to<Hubbard::DensityOfStates::dos_precision>();
		}
		Utility::saveData(abscissa, square_dos.values, "../../data/2d_dos.dat.gz");
	}

	test_e = std::chrono::steady_clock::now();
	std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
	std::cout << "\n\n" << std::endl;

	//------------------------------------------------------------//

	//SquareLattice::HubbardCDW model(modelParameters);
	//model.computePhases({false, true}).print();
	//std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;
	//
	//test_e = std::chrono::steady_clock::now();
	//std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
	//std::cout << "\n\n" << std::endl;

	//return _DEFAULT_EXIT;
	//------------------------------------------------------------//

	SquareLattice::UsingBroyden model2(modelParameters);
	test_b = std::chrono::steady_clock::now();
	model2.computePhases({ false, true }).print();
	std::cout << "Free energy = " << model2.freeEnergyPerSite() << std::endl;

	test_e = std::chrono::steady_clock::now();
	std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
}