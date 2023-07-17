#pragma once
#include <Eigen/Dense>
#include "../Utility/Resolvent.hpp"

namespace Hubbard {
	struct PhaseDebuggingPolicy {
		bool printAll{ false };
		bool convergenceWarning{ true };

		PhaseDebuggingPolicy() = default;
		PhaseDebuggingPolicy(bool _printAll, bool _convergenceWarning)
			: printAll{ _printAll }, convergenceWarning(_convergenceWarning) {};
	};

	constexpr double L_PI = 3.141592653589793238462643383279502884L; /* pi */
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<double, Eigen::Dynamic> Vector_L;

	typedef std::complex<double> complex_prec;
	typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<complex_prec, Eigen::Dynamic> VectorCL;
	using SpinorMatrix = MatrixCL;
	typedef Utility::Resolvent<double> Resolvent_L;

	template <const int vector_size = Eigen::Dynamic>
	void printAsRow(Eigen::Vector<double, vector_size>& printer) {
		for (size_t i = 0U; i < printer.size(); ++i)
		{
			std::cout << " \t" << printer(i);
			if ((i + 1U) % 8U == 0U) {
				std::cout << "\n\t    ";
			}
		}
		std::cout << std::endl;
	}

	template <const int vector_size = Eigen::Dynamic>
	void printAsRow(Eigen::Vector<complex_prec, vector_size>& printer) {
		for (size_t i = 0U; i < printer.size(); ++i)
		{
			std::cout << " \t" << printer(i);
			if ((i + 1U) % 4U == 0U) {
				std::cout << "\n\t    ";
			}
		}
		std::cout << std::endl;
	}
	const complex_prec I = { 0, 1 };

	using ComplexParameterVector = Eigen::Vector<complex_prec, Eigen::Dynamic>;

	inline void complexParametersToReal(const ComplexParameterVector& c, Eigen::VectorXd& r) {
		r(0) = c(0).real(); // CDW
		r(1) = c(1).real(); // AFM
		r(2) = c(2).real(); // SC
		r(3) = c(3).real(); // Gamma SC
		r(4) = c(4).imag(); // Xi SC
		r(5) = c(5).imag(); // Eta
		r(6) = c(6).real(); // Gamma Occupation Up
		r(7) = c(7).real(); // Gamma Occupation Down
	};
}