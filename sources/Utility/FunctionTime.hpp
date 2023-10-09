#pragma once
#include <chrono>
#include <iostream>
#include <utility>
#include <functional>
#include <type_traits>

namespace Utility {
	// Prints the runtime of a function to the console and returns the function's return value
	template <class FunctionTypes, class... Args>
	auto function_time_ms(const std::function<FunctionTypes>& func, Args&&... args)
		-> decltype(func(std::forward<Args>(args)...))
	{
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		if constexpr (!(std::is_same_v<decltype(func(std::forward<Args>(args)...)), void>)) {
			auto return_value = func(std::forward<Args>(args)...);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			return return_value;
		}
		else {
			func(std::forward<Args>(args)...);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		}
	}

	// Prints the runtime of a member function of ObjectType to the console and returns the function's return value
	template <class ObjectType, class ReturnType, class... Args>
	ReturnType member_function_time_ms(ObjectType& obj, ReturnType(ObjectType::* function)(Args...), Args&&... args)
	{
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		if constexpr (!(std::is_same_v<decltype(func(std::forward<Args>(args)...)), void>)) {
			auto return_value = (obj.*function)(std::forward<Args>(args)...);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			return return_value;
		}
		else {
			(obj.*function)(std::forward<Args>(args)...);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		}
	}
}