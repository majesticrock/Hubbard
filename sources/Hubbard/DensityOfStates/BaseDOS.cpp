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

	void BaseDOS::writeToBinaryFile(const std::string& filename){
		std::ofstream writer(filename, std::ios::out | std::ios::binary);
		if(!writer){
			std::cerr << "Could not open file stream in writeToBinaryFile - " << filename << std::endl;
			return;
		}
		writer.write((char *) &values.size(), sizeof(size_t));
		for(const auto& value : values){
			writer.write((char *) &value, sizeof(value));
		}

		writer.write((char *) &abscissa.size(), sizeof(size_t));
		for(const auto& absc : abscissa){
			writer.write((char *) &absc, sizeof(absc));
		}

		writer.write((char *) &weights.size(), sizeof(size_t));
		for(const auto& weight : weights){
			writer.write((char *) &weight, sizeof(weight));
		}
		writer.close();
		std::cout << "size = " << values.size() << std::endl;
	}

	bool loadFromBinaryFile(const std::string& filename){
		std::ifstream reader(filename, std::ios::out | std::ios::binary);
		if(!reader){
			std::cerr << "Could not open file stream for " << filename << std::endl;
			return false;
		}
		size_t vector_size;
		reader.read((char *) &vector_size, sizeof(size_t));
		std::cout << "size = " << vector_size << std::endl;
	}
}