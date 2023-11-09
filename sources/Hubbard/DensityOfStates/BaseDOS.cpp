#include "BaseDOS.hpp"
#include <iostream>
#include <numeric>

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
		std::ofstream writer(filename, std::ios::out | std::ios::binary);
		if (!writer) {
			std::cerr << "Could not open file stream in writeToBinaryFile - " << filename << std::endl;
			return;
		}
		// we need to create a proper lvalue variable
		// In order to read its address and by extent byte representation later on
		size_t vector_size = values.size();
		writer.write((char*)&vector_size, sizeof(size_t));
		for (const auto& value : values) {
			writer.write((char*)&value, sizeof(value));
		}

		for (const auto& absc : abscissa) {
			writer.write((char*)&absc, sizeof(absc));
		}

		for (const auto& weight : weights) {
			writer.write((char*)&weight, sizeof(weight));
		}
		writer.close();
	}

	bool BaseDOS::loadFromBinaryFile(const std::string& filename) {
		std::ifstream reader(filename, std::ios::out | std::ios::binary);
		if (!reader) {
			std::cerr << "Could not open file stream for " << filename << std::endl;
			return false;
		}
		size_t vector_size;
		reader.read((char*)&vector_size, sizeof(size_t));

		values.clear();
		values.resize(vector_size);
		weights.clear();
		weights.resize(vector_size);
		abscissa.clear();
		abscissa.resize(vector_size);

		for (size_t i = 0U; i < vector_size; ++i)
		{
			reader.read((char*)&values[i], sizeof(decltype(values)::value_type));
		}
		for (size_t i = 0U; i < vector_size; ++i)
		{
			reader.read((char*)&abscissa[i], sizeof(decltype(abscissa)::value_type));
		}
		for (size_t i = 0U; i < vector_size; ++i)
		{
			reader.read((char*)&weights[i], sizeof(decltype(weights)::value_type));
		}
		reader.close();
		if (!reader.good()) {
			std::cerr << "An error occurred while reading the dos data." << std::endl;
		}
		computed = reader.good();
		return computed;
	}
}