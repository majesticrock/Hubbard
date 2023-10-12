#pragma once
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>

namespace Utility {
	// Casts a floating point number to a std::string with desired precision n
	template <typename T>
	std::string to_string_with_precision(const T a_value, const int n)
	{
		std::ostringstream out;
		out.precision(n);
		out << std::fixed << a_value;
		return out.str();
	}

	// Returns the local time - the preprocessor statements deal with different OSs
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

	// Returns the current time in a readable format
	// default = "YYYY-MM-DD HH:MM:SS"
	inline std::string time_stamp(const std::string& fmt = "%F %T")
	{
		auto bt = localtime_xp(std::time(0));
		char buf[64];
		return { buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt) };
	}

	// data_type is a type that overloads the << operator - e.g. float or double, but custom types are allowed too
	// outstream_type is a type of output stream, e.g. std::ofstream or std::ostringstream
	template <typename data_type, typename outstream_type>
	class OutputWriter {
	public:
		inline bool checkDataForNaN(const std::vector<data_type>& data) const{
			for (const auto& value : data) {
				if (std::isnan(data_type)) {
					std::cerr << "Atleast one of your data points is NaN!" << std::endl;
					break;
				}
			}
		};
		// Writes the current time stamp and comments to the file
		// If the latter is not provided only the time stamp is written
		void writeComments(outstream_type& out, const std::vector<std::string>& comments = std::vector<std::string>()) const
		{
			out << "# " << time_stamp() << "\n#\n";
			for (int i = 0; i < comments.size(); i++)
			{
				out << "# " << comments[i] << "\n";
			}
		};

		// Appends a line consisting of <data> to <out>
		void appendLine(const std::vector<data_type>& data, outstream_type& out) const
		{
			checkDataForNaN(data);
			for (int i = 0; i < data.size(); i++)
			{
				out << data[i];
				if (i < data.size() - 1) {
					out << " ";
				}
			}
			out << "\n";
		}

		void saveData(const data_type& data, outstream_type& out,
			const std::vector<std::string>& comments = std::vector<std::string>()) const
		{
			writeComments(out, comments);
			out << std::scientific << std::setprecision(10);
			out << data;
		};

		void saveData(const std::vector<data_type>& data, outstream_type& out,
			const std::vector<std::string>& comments = std::vector<std::string>()) const
		{
			checkDataForNaN(data);
			writeComments(out, comments);
			out << std::scientific << std::setprecision(10);
			for (const auto& dat : data) {
				out << dat << "\n";
			}
		};

		void saveData(const std::vector<std::vector<data_type>>& data, outstream_type& out,
			const std::vector<std::string>& comments = std::vector<std::string>()) const
		{
			writeComments(out, comments);
			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				appendLine(data[i], out);
			}
		};
	};
}