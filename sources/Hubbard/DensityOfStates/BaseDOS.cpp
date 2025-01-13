#include "BaseDOS.hpp"
#include <iostream>
#include <numeric>
#include <mrock/utility/BinaryIO.hpp>

namespace Hubbard::DensityOfStates {
	std::vector<dos_precision> BaseDOS::values;
	std::vector<abscissa_t> BaseDOS::abscissa;
	std::vector<dos_precision> BaseDOS::weights;
	bool BaseDOS::computed = false;
	dos_precision BaseDOS::step = 0;

	dos_precision BaseDOS::integrateValues()
	{
		return step * std::reduce(values.begin(), values.end());
	}

	void BaseDOS::printValues()
	{
		for (const auto& val : values) {
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
	void BaseDOS::clearAll()
	{
		values.clear();
		abscissa.clear();
		weights.clear();
		computed = false;
		step = 0;
	}
	void BaseDOS::printValuesAndAbscissa()
	{
		for (size_t i = 0U; i < values.size(); ++i)
		{
			std::cout << abscissa[i] << "; " << values[i] << std::endl;
		}
		std::cout << std::endl;
	}

	void BaseDOS::writeToBinaryFile(const std::string& filename) {
		std::ofstream writer = mrock::utility::BinaryIO::create_writer(filename);
		if (!writer) {
			std::cerr << "Could not open file stream in writeToBinaryFile - " << filename << std::endl;
			return;
		}
		// we need to create a proper lvalue variable
		// In order to read its address and by extent byte representation later on
		size_t vector_size = values.size();
		mrock::utility::BinaryIO::writeVariable(vector_size, writer);
		mrock::utility::BinaryIO::writeVector(values, writer);
		mrock::utility::BinaryIO::writeVector(abscissa, writer);
		mrock::utility::BinaryIO::writeVector(weights, writer);
		writer.close();
	}

	bool BaseDOS::loadFromBinaryFile(const std::string& filename) {
		std::ifstream reader = mrock::utility::BinaryIO::create_reader(filename);
		if (!reader) {
			std::cerr << "Could not open file stream for " << filename << std::endl;
			return false;
		}
		size_t vector_size;
		mrock::utility::BinaryIO::read_to_variable(vector_size, reader);
		values.resize(vector_size);
		weights.resize(vector_size);
		abscissa.resize(vector_size);

		mrock::utility::BinaryIO::readToVector(values, reader);
		mrock::utility::BinaryIO::readToVector(abscissa, reader);
		mrock::utility::BinaryIO::readToVector(weights, reader);

		reader.close();
		if (!reader.good()) {
			std::cerr << "An error occurred while reading the dos data." << std::endl;
		}
		computed = reader.good();
		return computed;
	}
}