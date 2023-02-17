#pragma once
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

namespace Utility {
	inline std::tm localtime_xp(std::time_t timer)
	{
		std::tm bt{};
#if defined(__unix__)
		localtime_r(&timer, &bt);
#elif defined(_MSC_VER)
		localtime_s(&bt, &timer);
#else
		static std::mutex mtx;
		std::lock_guard<std::mutex> lock(mtx);
		bt = *std::localtime(&timer);
#endif
		return bt;
	}

	// default = "YYYY-MM-DD HH:MM:SS"
	inline std::string time_stamp(const std::string& fmt = "%F %T")
	{
		auto bt = localtime_xp(std::time(0));
		char buf[64];
		return { buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt) };
	}

	template <typename T>
	std::string to_string_with_precision(const T a_value, const int n)
	{
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << a_value;
		return out.str();
	}

	void saveData(std::vector<double>& data, std::string filename);
	void saveData(std::vector<std::vector<double>>& data, std::string filename);
	void saveData(std::vector<std::vector<double>>& data, std::string filename, std::vector<std::string> comments);
	void saveData(std::vector<double>& data, const int linebreak, std::string filename, std::vector<std::string> comments);
}