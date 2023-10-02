#pragma once
#include "OutputWriter.hpp"
#include <iostream>

#ifndef _DEBUG
#define _USE_BOOST // Disable if boost ist not desired
#endif

#ifdef _USE_BOOST
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <sstream>
#endif

namespace Utility {
#ifdef _USE_BOOST
	template <typename data_type>
	void saveData(const data_type& data, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		// create file
		std::ofstream ofile(filename, std::ios_base::out | std::ios_base::binary);
		if (ofile.is_open()) {
			boost::iostreams::filtering_ostream out;
			out.push(boost::iostreams::gzip_compressor());
			out.push(ofile);

			OutputWriter<data_type, std::ostringstream> ow;
			// The extra stringstream is a proper std::ostream (boost::iostreams are not)
			// and can thereby be used with the standard operator<< overloading
			std::ostringstream oss;
			ow.saveData(data, oss, comments);
			out << oss.str();

			out.pop(); // flushes the filter chain and closes the file
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}

	template <typename data_type>
	void saveData(const std::vector<data_type>& data, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		// create file
		std::ofstream ofile(filename, std::ios_base::out | std::ios_base::binary);
		if (ofile.is_open()) {
			boost::iostreams::filtering_ostream out;
			out.push(boost::iostreams::gzip_compressor());
			out.push(ofile);

			OutputWriter<data_type, std::ostringstream> ow;
			// The extra stringstream is a proper std::ostream (boost::iostreams are not)
			// and can thereby be used with the standard operator<< overloading
			std::ostringstream oss;
			ow.saveData(data, oss, comments);
			out << oss.str();

			out.pop(); // flushes the filter chain and closes the file
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}

	template <typename data_type>
	void saveData(const std::vector<std::vector<data_type>>& data, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		// create file
		std::ofstream ofile(filename, std::ios_base::out | std::ios_base::binary);
		if (ofile.is_open()) {
			boost::iostreams::filtering_ostream out;
			out.push(boost::iostreams::gzip_compressor());
			out.push(ofile);

			OutputWriter<data_type, std::ostringstream> ow;
			// The extra stringstream is a proper std::ostream (boost::iostreams are not)
			// and can thereby be used with the standard operator<< overloading
			std::ostringstream oss;
			ow.saveData(data, oss, comments);
			out << oss.str();

			out.pop(); // flushes the filter chain and closes the file
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	}

	template <typename data_type>
	void saveData(const std::vector<data_type>& first, const std::vector<data_type>& second, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		saveData(std::vector<std::vector<data_type>>{ first, second }, filename, comments);
	};

	// This function assumes that the number of elements of <data> is divisible by linebreak
	template <typename data_type>
	void saveData(const std::vector<data_type>& data, size_t linebreak, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		if (linebreak == 0) {
			linebreak = data.size();
		}
		else if (data.size() % linebreak != 0) {
			std::cerr << "The number of data elements is not divisible by linebreak!" << std::endl;
			return;
		}
		// create file
		std::ofstream ofile(filename, std::ios_base::out | std::ios_base::binary);
		if (ofile.is_open()) {
			boost::iostreams::filtering_ostream out;
			out.push(boost::iostreams::gzip_compressor());
			out.push(ofile);

			OutputWriter<data_type, std::ostringstream> ow;
			// The extra stringstream is a proper std::ostream (boost::iostreams are not)
			// and can thereby be used with the standard operator<< overloading
			std::ostringstream oss;
			ow.writeComments(oss, comments);
			for (size_t n = 0; n < data.size(); n += linebreak)
			{
				std::vector<data_type> part_vec(data.begin() + n, data.begin() + n + linebreak);
				ow.appendLine(part_vec, oss);
			}
			out << oss.str();

			out.pop(); // flushes the filter chain and closes the file
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	};
#else // End using boost
	// Provides an easy-to-use method that uses std::ofstream to write <data> to <filename> in plain text
	template <typename data_type>
	void saveData(const std::vector<data_type>& data, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			OutputWriter<data_type, std::ofstream> ow;
			ow.saveData(data, out, comments);
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	};

	// Provides an easy-to-use method that uses std::ofstream to write <data> to <filename> in plain text
	template <typename data_type>
	void saveData(const std::vector<std::vector<data_type>>& data, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		std::ofstream out(filename);
		if (out.is_open()) {
			OutputWriter<data_type, std::ofstream> ow;
			ow.saveData(data, out, comments);
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	};

	// This function assumes that the number of elements of <data> is divisible by linebreak
	template <typename data_type>
	void saveData(const std::vector<data_type>& data, size_t linebreak, const std::string& filename,
		const std::vector<std::string>& comments)
	{
		if (linebreak == 0U) {
			linebreak = data.size();
		}
		else if (data.size() % linebreak != 0U) {
			std::cerr << "The numbe rof data elements is not divisible by linebreak!" << std::endl;
			return;
		}
		std::ofstream out(filename);
		if (out.is_open()) {
			OutputWriter<data_type, std::ofstream> ow;
			ow.writeComments(out, comments);

			for (size_t n = linebreak; n < data.size(); n += linebreak)
			{
				std::vector<data_type> part_vec(data.begin() + n - linebreak, data.begin() + n);
				ow.appendLine(part_vec, out);
			}
		}
		else {
			std::cerr << "Could not open output filestream for file: " << filename << std::endl;
		}
	};

	template <typename data_type>
	void saveData(const std::vector<data_type>& first, const std::vector<data_type>& second, const std::string& filename,
		const std::vector<std::string>& comments = std::vector<std::string>())
	{
		saveData(std::vector<std::vector<data_type>>{ first, second }, filename, comments);
	};
#endif // End not using boost
}