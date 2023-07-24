#include "BaseDOS.hpp"
#include <iostream>
#include <numeric>
#include <boost/math/quadrature/gauss.hpp>

namespace Hubbard::DensityOfStates {
	std::vector<double> BaseDOS::values;
	bool BaseDOS::computed = false;
	double BaseDOS::step = 0;
	void BaseDOS::renormalize()
	{
		// Factor of 2 because we only computed half the DOS due to symmetry
		double norm = getNorm();
		std::cout << "DOS-Norm: " << norm << std::endl;
		for (auto& val : values) {
			val /= norm;
		}
	}
	void BaseDOS::printValues()
	{
		for (const auto& val : values) {
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
	double BaseDOS::getNorm()
	{
		//constexpr size_t num_positions = 30U;
		//auto proc = [](double k) {
		//	for (size_t i = 0U; i < num_positions; ++i)
		//	{
		//		if (boost::math::quadrature::gauss<double, num_positions>::abscissa()[i] == k
		//			|| boost::math::quadrature::gauss<double, num_positions>::abscissa()[i] == -k) {
		//			return values[i];
		//		}
		//	}
		//	throw std::runtime_error("Did not find k within the abscissa!");
		//};
		//return 2 * boost::math::quadrature::gauss<double, num_positions>::integrate(proc, 0., 3.);
		

		return step * std::accumulate(values.begin(), values.end(), double{});
	}
}