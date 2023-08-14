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

#define _LONG_PRECISION
#ifdef _LONG_PRECISION
#define _MPI_RETURN_TYPE MPI_LONG_DOUBLE
	typedef long double global_floating_type;
#else
#define _MPI_RETURN_TYPE MPI_DOUBLE
	typedef double global_floating_type;
#endif

	constexpr global_floating_type L_PI = 3.141592653589793238462643383279502884L; /* pi */
	constexpr global_floating_type ONE_HALF = 0.5L;
	typedef Eigen::Matrix<global_floating_type, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<global_floating_type, Eigen::Dynamic> Vector_L;

	typedef std::complex<global_floating_type> complex_prec;
	typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<complex_prec, Eigen::Dynamic> VectorCL;
	using SpinorMatrix = MatrixCL;
	typedef Utility::Resolvent<global_floating_type> Resolvent_L;

	template <const int vector_size = Eigen::Dynamic>
	void printAsRow(Eigen::Vector<global_floating_type, vector_size>& printer) {
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

	inline void complexParametersToReal(const ComplexParameterVector& c, Eigen::Vector<global_floating_type, -1>& r) {
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