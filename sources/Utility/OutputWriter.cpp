#include "OutputWriter.hpp"
#include <fstream>
#include <iostream>
#include <ctime>

namespace Utility {
	void saveData(const std::vector<double>& data, const std::string& filename)
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n\n";

			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				out << data[i];
				out << " ";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
	void saveData(const std::vector<std::vector<double>>& data, const std::string& filename)
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n\n";

			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				for (int s = 0; s < data[0].size(); s++)
				{
					out << data[i][s];
					if (s < data[0].size() - 1) {
						out << " ";
					}
				}
				out << "\n";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
	void saveData(const std::vector<std::vector<double>>& data, const std::string& filename, const std::vector<std::string>& comments)
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n#\n";

			for (int i = 0; i < comments.size(); i++)
			{
				out << "# " << comments[i] << "\n";
			}

			out << "\n" << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				for (int s = 0; s < data[0].size(); s++)
				{
					out << data[i][s];
					if (s < data[0].size() - 1) {
						out << " ";
					}
				}
				out << "\n";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
	void saveData(const std::vector<double>& data, const int linebreak, const std::string& filename, const std::vector<std::string>& comments) {
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << time_stamp() << "\n#\n";

			for (int i = 0; i < comments.size(); i++)
			{
				out << "# " << comments[i] << "\n";
			}

			out << std::scientific << std::setprecision(10);
			for (int i = 0; i < data.size(); i++)
			{
				if (i % linebreak == 0) {
					out << "\n";
				}
				out << data[i];
				out << " ";
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}
}