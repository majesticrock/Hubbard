#include "Model.hpp"
#include <string>
#include <sstream>

namespace Hubbard {
	Model::ModelParameters::ModelParameters(double_prec _temperature, double_prec _U, double_prec _V, double_prec _global_step, double_prec _second_step,
		std::string _global_iterator_type, std::string _second_iterator_type)
		: global_iterator_type(_global_iterator_type), second_iterator_type(_second_iterator_type),
		global_step(_global_step), second_step(_second_step), temperature(_temperature), U(_U), V(_V)
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

	void Model::ModelParameters::incrementer(std::string& s, const double_prec step)
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
	double_prec Model::ModelParameters::setGlobalIterator(int it_num)
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
	double_prec Model::ModelParameters::setGlobalIteratorExact(double_prec newValue)
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
	double_prec Model::ModelParameters::setSecondIterator(int it_num)
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
	double_prec Hubbard::Model::ModelParameters::setSecondIteratorExact(double_prec newValue)
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
	void Model::ModelParameters::incrementGlobalIterator()
	{
		incrementer(global_iterator_type, global_step);
		setSecondIterator(0);
	}
	void Model::ModelParameters::incrementSecondIterator()
	{
		incrementer(second_iterator_type, second_step);
	}
	void Model::ModelParameters::printGlobal() const
	{
		std::cout << global_iterator_type << " = " << getGlobal();
	}
	std::string Model::ModelParameters::getFileName() const
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
	void Model::data_set::print() const {
		std::cout << delta_cdw << "\t" << delta_afm << "\t" << delta_sc << "\t" << delta_eta
			<< "\n    Delta_tot = " << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
	}
}