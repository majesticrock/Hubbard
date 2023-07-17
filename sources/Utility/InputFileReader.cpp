#include "InputFileReader.hpp"
#include <iostream>
#include <string>

namespace Utility {
	InputFileReader::InputFileReader(std::string fileName) {
		std::ifstream f;

		f.open(fileName, std::ios::in);

		std::string s;
		while (std::getline(f, s)) {
			if (s[0] == '#') continue;
			size_t pos = s.find(" ");
			std::string name = s.substr(0, pos);
			std::string content = s.substr(pos + 1);

			if (content.back() == 13) {
				content.erase(prev(content.end()));
			}

			names.push_back(name);
			contents.push_back(content);
			used.push_back(false);
		}
	}

	int InputFileReader::find(std::string s) {
		int pos = -1;
		bool done = false;
		for (unsigned int i = 0U; i < names.size() && !done; ++i) {
			if (names[i] == s) {
				pos = i;
				done = true;
			}
		}
		return pos;
	}

	bool InputFileReader::is(std::string name) {
		int i = find(name);
		if (i != -1) {
			used[i] = true;
			return true;
		}
		else {
			return false;
		}
	}

	bool InputFileReader::getBool(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;
		if (contents[i] == "true") {
			return true;
		}
		else if (contents[i] == "false") {
			return false;
		}
		else {
			throw std::invalid_argument("The Parameter " + name + " is not bool; " + contents[i]);
		}
	}

	int InputFileReader::getInt(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;
		return stoi(contents[i]);
	}

	std::vector<int> InputFileReader::getIntList(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;

		std::string s = contents[i];
		std::vector<int> result;
		size_t pos = 0;
		bool done = false;
		while (!done) {
			result.push_back(stoi(s, &pos));
			if (s.size() > pos) {
				s = s.substr(pos + 1);
			}
			else {
				done = true;
			}
		}
		return result;
	}

	double InputFileReader::getDouble(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;
		return stod(contents[i]);
	}

	std::vector<double> InputFileReader::getDoubleList(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;

		std::string s = contents[i];
		std::vector<double> result;
		size_t pos = 0;
		bool done = false;
		while (!done) {
			result.push_back(stod(s, &pos));
			if (s.size() > pos) {
				s = s.substr(pos + 1);
			}
			else {
				done = true;
			}
		}
		return result;
	}

	std::string InputFileReader::getString(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;

		return contents[i];
	}

	std::vector<std::string> InputFileReader::getStringList(std::string name) {
		int i = find(name);
		if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
		used[i] = true;

		std::string s = contents[i];
		std::vector<std::string> result;
		bool done = false;
		while (!done) {
			size_t pos2 = s.find(" ");
			if (pos2 == std::string::npos) {
				std::string content = s;
				result.push_back(content);
				done = true;
			}
			else {
				std::string content = s.substr(0, pos2);
				s = s.substr(pos2 + 1);
				result.push_back(content);
			}
		}
		return result;
	}

	bool InputFileReader::allUsed() {
		bool result = true;
		for (unsigned int i = 0; i < used.size(); i++) {
			if (used[i] == false) result = false;
		}
		return result;
	}

	std::vector<std::string> InputFileReader::listNotUsed() {
		std::vector<std::string> result;
		for (unsigned int i = 0; i < used.size(); i++) {
			if (!used[i]) {
				result.push_back(names[i]);
			}
		}
		return result;
	}
}