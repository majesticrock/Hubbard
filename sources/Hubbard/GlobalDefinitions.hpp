#pragma once
#include <Eigen/Dense>
#include <stddef.h>
#include "../Utility/Resolvent.hpp"
#include <boost/math/constants/constants.hpp>

//#define _BOOST_PRECISION
//#define _LONG_PRECISION
#ifdef _BOOST_PRECISION
#include <boost/multiprecision/cpp_bin_float.hpp>
#define SQUARE_QUAD_CUT_OFF -30
#define SC_QUAD_CUT_OFF -15
#else
#define SQUARE_QUAD_CUT_OFF -14
#define SC_QUAD_CUT_OFF -12
#endif

// Does not quite work yet - but the equally spaced variant is just fine
// Evenly spaced might even be required
//#define _EXACT_DOS

#include <cmath>
#include <complex>
#include <string>

namespace Hubbard {
	using std::abs;
	using std::to_string;

	struct PhaseDebuggingPolicy {
		bool printAll{ false };
		bool convergenceWarning{ true };

		constexpr PhaseDebuggingPolicy() = default;
		constexpr PhaseDebuggingPolicy(bool _printAll, bool _convergenceWarning)
			: printAll{ _printAll }, convergenceWarning(_convergenceWarning) {};
	};
	constexpr PhaseDebuggingPolicy WarnNoConvergence{ false, true };
	constexpr PhaseDebuggingPolicy NoWarning{ false, false };
	constexpr PhaseDebuggingPolicy PrintSteps{ true, true };

#ifdef _BOOST_PRECISION
#define _NO_MPI
	typedef boost::multiprecision::cpp_bin_float_100 global_floating_type;
	// At least long double, but maybe larger
	typedef boost::multiprecision::cpp_bin_float_100 long_double_t;
#define _CONST_FLOATING const long_double_t
#define _CONST_LONG_FLOATING const long_double_t
#else
#ifdef _LONG_PRECISION
#define _CONST_FLOATING constexpr long double
#define _MPI_RETURN_TYPE MPI_LONG_DOUBLE
	typedef long double global_floating_type;
#else
#define _CONST_FLOATING constexpr double
#define _MPI_RETURN_TYPE MPI_DOUBLE
	typedef double global_floating_type;
#endif // _LONG_PRECISION
	// At least long double, but maybe larger
	typedef long double long_double_t;
#define _CONST_LONG_FLOATING constexpr long_double_t
#endif // _BOOST_PRECISION

	_CONST_LONG_FLOATING LONG_PI = boost::math::constants::pi<long_double_t>(); // pi
	// tak, v kancelari mam boost 1.71 a ten nedefinuje 1/pi ale 2/pi. je to na hovno ale nelze delat nic
	_CONST_LONG_FLOATING LONG_1_PI = 0.5 * boost::math::constants::two_div_pi<long_double_t>(); // 1 / pi
	_CONST_LONG_FLOATING LONG_PI_2 = boost::math::constants::half_pi<long_double_t>(); // pi / 2
	_CONST_LONG_FLOATING LOG_4 = 2 * boost::math::constants::ln_two<long_double_t>(); // ln(4) = 2 ln(2)
	_CONST_FLOATING ONE_HALF = 0.5L;
	_CONST_FLOATING BASE_PI = boost::math::constants::pi<global_floating_type>(); // pi
	_CONST_FLOATING ONE_OVER_SQRT_2 = boost::math::constants::one_div_root_two<global_floating_type>(); // 1 / sqrt(2);

	typedef Eigen::Matrix<global_floating_type, Eigen::Dynamic, Eigen::Dynamic> Matrix_L;
	typedef Eigen::Vector<global_floating_type, Eigen::Dynamic> Vector_L;

	typedef std::complex<global_floating_type> complex_prec;
	typedef Eigen::Matrix<complex_prec, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<complex_prec, Eigen::Dynamic> VectorCL;
	typedef Eigen::Array<complex_prec, Eigen::Dynamic, Eigen::Dynamic> ValueArray;
	using SpinorMatrix = MatrixCL;

	using ResolventReal = Utility::Resolvent<global_floating_type, false>;
	using ResolventComplex = Utility::Resolvent<global_floating_type, true>;
	using ResolventReturnData = Utility::ResolventDataWrapper<global_floating_type>;

	template <int vector_size = Eigen::Dynamic>
	void printAsRow(Eigen::Vector<global_floating_type, vector_size>& printer) {
		for (size_t i = 0U; i < printer.size(); ++i)
		{
			std::cout << "   " << printer(i);
		}
		std::cout << std::endl;
	}

	template <int vector_size = Eigen::Dynamic>
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

#ifdef _BOOST_PRECISION
	const complex_prec I = { 0, 1 };
#else
	constexpr complex_prec I = { 0, 1 };
#endif

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
		r(8) = c(8).real(); // PS
	};
}