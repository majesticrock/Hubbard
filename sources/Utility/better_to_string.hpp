#pragma once
#include <string>
#include <array>
#include <charconv>
#include <iostream>
#include <string_view>
#include <system_error>

namespace Utility {
	template <typename... Args>
	std::string better_to_string(Args&&... format_args)
	{
		std::array<char, 18> str;
		auto result = std::to_chars(str.data(), str.data() + str.size(), std::forward<Args>(format_args)...);

		if (result.ec == std::errc())
			return std::string(str.data(), result.ptr - str.data());
		else
			return std::make_error_code(result.ec).message();
	}

	template <typename FloatingPoint>
	std::string better_to_string_fixed(FloatingPoint number, unsigned short digits)
	{
		std::array<char, 18> str;
		auto result = std::to_chars(str.data(), str.data() + str.size(), number, std::chars_format::fixed, digits);

		if (result.ec == std::errc())
			return std::string(str.data(), result.ptr - str.data());
		else
			return std::make_error_code(result.ec).message();
	}
}