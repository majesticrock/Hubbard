#ifndef _InputFileReader_hpp_
#define _InputFileReader_hpp_

#include <fstream>
#include <vector>

/**
* InputFileReader
* This is a general class to handle input files.
* Each line contains a parameter name and the content of the parameter seperated by a space.
* The content can be a bool, int, double or string. Or a list of them seperated by spaces.
* The order of the parameters is irrelevant
*
* Example file:
*
* method Classic
* M 1024
* ak_dist Box
* N_Box 100
* J_dist Const
* mode Mode1
* t_p 0 0.1 0.4
*
* Example code:
* InputFileReader input("input.config");
* string method = getString("method");
* cout << method << endl // Output: "Classic"
* int M = getInt("M");
* cout << M << endl // Output: 1024
* etc.
**/
namespace Utility {
	class InputFileReader {
	private:
		// a list of the parameter names
		std::vector<std::string> names;
		// a list of the parameter content
		// this strings have to be decrypted when a getX() method is called
		std::vector<std::string> contents;
		// tracking which parameters are used
		std::vector<bool> used;

		// get the index of the parameter name
		// return -1 if not found
		int find(std::string s);

	public:
		// Read an input file with the name "fileName"
		InputFileReader(std::string fileName);

		// checks wether the parameter name is specified
		bool is(std::string name);
		// get the content as boolean
		bool getBool(std::string name);
		// get the content as integer
		int getInt(std::string name);
		// get the content as list of integers
		std::vector<int> getIntList(std::string name);
		// get the content as double
		double getDouble(std::string name);
		// get the content as list of doubles
		std::vector<double> getDoubleList(std::string name);
		// get the content as string
		std::string getString(std::string name);
		// get the content as list of strings
		std::vector<std::string> getStringList(std::string name);

		// check wether all parameters are used
		// a parameter is used when a is() method is called
		bool allUsed();
		// list the names of all unused parameters
		std::vector<std::string> listNotUsed();
	};
}
#endif