#include "InputFileReader.hpp"
#include <iostream>
#include <string>

namespace Utility {
	using namespace std;

	InputFileReader::InputFileReader(string fileName) {
		ifstream f;

		f.open(fileName, ios::in);

		string s;
		while (getline(f, s)) {
			if (s[0] == '#') continue;
			size_t pos = s.find(" ");
			string name = s.substr(0, pos);
			string content = s.substr(pos + 1);

			if (content.back() == 13) {
				content.erase(prev(content.end()));
			}

			names.push_back(name);
			contents.push_back(content);
			used.push_back(false);
		}
	}

	int InputFileReader::find(string s) {
		int pos = -1;
		bool done = false;
		for (unsigned int i = 0; i < names.size() && !done; i++) {
			if (names[i] == s) {
				pos = i;
				done = true;
			}
		}
		return pos;
	}

	bool InputFileReader::is(string name) {
		int i = find(name);
		if (i != -1) {
			used[i] = true;
			return true;
		}
		else {
			return false;
		}
	}

	bool InputFileReader::getBool(string name) {
		int i = find(name);
		if (i == -1) abort();
		used[i] = true;
		if (contents[i] == "true") {
			return true;
		}
		else if (contents[i] == "false") {
			return false;
		}
		else {
			cerr << "Parameter " << name << " not found in Inputfile" << endl;
			abort();
		}
	}

	int InputFileReader::getInt(string name) {
		int i = find(name);
		if (i == -1) abort();
		used[i] = true;
		return stoi(contents[i]);
	}

	vector<int> InputFileReader::getIntList(string name) {
		int i = find(name);
		if (i == -1) {
			cerr << "Parameter " << name << " not found in Inputfile" << endl;
			abort();
		}
		used[i] = true;

		string s = contents[i];
		vector<int> result;
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

	double InputFileReader::getDouble(string name) {
		int i = find(name);
		if (i == -1) {
			cerr << "Parameter " << name << " not found in Inputfile" << endl;
			abort();
		}
		used[i] = true;
		return stod(contents[i]);
	}

	vector<double> InputFileReader::getDoubleList(string name) {
		int i = find(name);
		if (i == -1) {
			cerr << "Parameter " << name << " not found in Inputfile" << endl;
			abort();
		}
		used[i] = true;

		string s = contents[i];
		vector<double> result;
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

	string InputFileReader::getString(string name) {
		int i = find(name);
		if (i == -1) {
			cerr << "Parameter " << name << " not found in Inputfile" << endl;
			abort();
		}
		used[i] = true;

		return contents[i];
	}

	vector<string> InputFileReader::getStringList(string name) {
		int i = find(name);
		if (i == -1) {
			cerr << "Parameter " << name << " not found in Inputfile" << endl;
			abort();
		}
		used[i] = true;

		string s = contents[i];
		vector<string> result;
		bool done = false;
		while (!done) {
			size_t pos2 = s.find(" ");
			if (pos2 == string::npos) {
				string content = s;
				result.push_back(content);
				done = true;
			}
			else {
				string content = s.substr(0, pos2);
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

	vector<string> InputFileReader::listNotUsed() {
		vector<string> result;
		for (unsigned int i = 0; i < used.size(); i++) {
			if (!used[i]) {
				result.push_back(names[i]);
			}
		}
		return result;
	}
}