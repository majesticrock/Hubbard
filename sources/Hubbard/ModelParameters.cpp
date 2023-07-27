#include "ModelParameters.hpp"
#include <iostream>
#include <sstream>
#include <cmath>

namespace Hubbard {
	void ModelParameters::init()
	{
		if (global_iterator_type == "T") {
			global_it_min = temperature;
		}
		else if (global_iterator_type == "U") {
			global_it_min = U;
		}
		else if (global_iterator_type == "V") {
			global_it_min = V;
		}
		else {
			global_it_min = 0;
		}
		if (second_iterator_type == "T") {
			second_it_min = temperature;
		}
		else if (second_iterator_type == "U") {
			second_it_min = U;
		}
		else if (second_iterator_type == "V") {
			second_it_min = V;
		}
		else {
			second_it_min = 0;
		}
	}

	ModelParameters::ModelParameters(double _temperature, double _U, double _V, double _global_step, double _second_step,
		std::string _global_iterator_type, std::string _second_iterator_type)
		: global_iterator_type(_global_iterator_type), second_iterator_type(_second_iterator_type),
		global_step(_global_step), second_step(_second_step), temperature(_temperature), U(_U), V(_V)
	{
		init();
	}

	ModelParameters::ModelParameters(const std::vector<double>& params, double _global_step, double _second_step, std::string _global_iterator_type, std::string _second_iterator_type)
		: global_iterator_type(_global_iterator_type), second_iterator_type(_second_iterator_type),
		global_step(_global_step), second_step(_second_step), temperature(params[0]), U(params[1]), V(params[2])
	{
		init();
	}

	void ModelParameters::incrementer(std::string& s, const double step)
	{
		if (s == "T") {
			temperature += step;
		}
		else if (s == "U") {
			U += step;
		}
		else if (s == "V") {
			V += step;
		}
	}
	double ModelParameters::setGlobalIterator(int it_num)
	{
		if (global_iterator_type == "T") {
			temperature = global_it_min + it_num * global_step;
		}
		else if (global_iterator_type == "U") {
			U = global_it_min + it_num * global_step;
		}
		else if (global_iterator_type == "V") {
			V = global_it_min + it_num * global_step;
		}
		return getGlobal();
	}
	double ModelParameters::setGlobalIteratorExact(double newValue)
	{
		if (global_iterator_type == "T") {
			temperature = newValue;
		}
		else if (global_iterator_type == "U") {
			U = newValue;
		}
		else if (global_iterator_type == "V") {
			V = newValue;
		}
		return getGlobal();
	}
	double ModelParameters::setSecondIterator(int it_num)
	{
		if (second_iterator_type == "T") {
			temperature = second_it_min + it_num * second_step;
		}
		else if (second_iterator_type == "U") {
			U = second_it_min + it_num * second_step;
		}
		else if (second_iterator_type == "V") {
			V = second_it_min + it_num * second_step;
		}
		return getSecond();
	}
	double Hubbard::ModelParameters::setSecondIteratorExact(double newValue)
	{
		if (second_iterator_type == "T") {
			temperature = newValue;
		}
		else if (second_iterator_type == "U") {
			U = newValue;
		}
		else if (second_iterator_type == "V") {
			V = newValue;
		}
		return getSecond();
	}
	void ModelParameters::incrementGlobalIterator()
	{
		incrementer(global_iterator_type, global_step);
		setSecondIterator(0);
	}
	void ModelParameters::incrementSecondIterator()
	{
		incrementer(second_iterator_type, second_step);
	}
	void ModelParameters::printGlobal() const
	{
		std::cout << global_iterator_type << " = " << getGlobal();
	}
	void ModelParameters::printParameters() const
	{
		std::cout << "[T U V] = [ " << temperature << " " << U << " " << V << " ]" << std::endl;
	}
	std::string ModelParameters::getFileName() const
	{
		auto improved_string = [](double number) -> std::string {
			if (std::floor(number) == number) {
				// If the number is a whole number, format it with one decimal place
				std::ostringstream out;
				out.precision(1);
				out << std::fixed << number;
				return out.str();
			}
			else {
				std::string str = std::to_string(number);
				// Remove trailing zeroes
				str.erase(str.find_last_not_of('0') + 1, std::string::npos);
				str.erase(str.find_last_not_of('.') + 1, std::string::npos);
				return str;
			}
		};

		std::string ret = "T=" + improved_string(temperature);
		ret += "/U=" + improved_string(U);
		ret += "_V=" + improved_string(V);

		ret += "/";
		return ret;
	}
}