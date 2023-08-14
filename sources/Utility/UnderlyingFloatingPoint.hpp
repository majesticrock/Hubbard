#pragma once
#include <type_traits>
#include <complex>

namespace Utility {
	template <class T>
	struct UnderlyingFloatingPoint {
		using type = T;
	};

	template <class T>
	struct UnderlyingFloatingPoint<std::complex<T>> {
		using type = T;
	};

	template <class T>
	using UnderlyingFloatingPoint_t = typename UnderlyingFloatingPoint<T>::type;
}