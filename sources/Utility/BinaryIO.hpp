#include <iostream>
#include <fstream>
#include <vector>

namespace Utility::BinaryIO {
	inline std::ifstream createReader(const std::string& filename) {
		return std::ifstream(filename, std::ios::out | std::ios::binary);
	};
	inline std::ofstream createWriter(const std::string& filename) {
		return std::ofstream(filename, std::ios::out | std::ios::binary);
	};

	template<class T>
	inline T& readToVariable(T& destination, std::ifstream& reader) {
		reader.read((char*)&destination, sizeof(T));
		return destination;
	};
	template<class T>
	inline void writeVariable(const T& source, std::ofstream& writer) {
		writer.write((char*)&source, sizeof(T));
	};

	template<class T>
	inline std::vector<T>& readToVector(std::vector<T>& destination, std::ifstream& reader) {
		for (auto& value : destination)
		{
			readToVariable(value, reader);
		}
		return destination;
	};
	template<class T>
	inline void writeVector(const std::vector<T>& source, std::ofstream& writer) {
		for (const auto& value : source) {
			writeVariable(value, writer);
		}
	};

	template<class T>
	std::ifstream readSerializedVector(std::vector<T>& destination, const std::string& filename) {
		auto reader = createReader(filename);
		if (!reader) {
			std::cerr << "Could not open file stream in readSerializedVector - " << filename << std::endl;
			return reader;
		}
		size_t vector_size;
		readToVariable(vector_size, reader);
		destination.resize(vector_size);
		readToVector(destination, reader);
		return reader;
	};
	template<class T>
	std::ofstream serializeVector(const std::vector<T>& source, const std::string& filename) {
		auto writer = createWriter(filename);
		if (!writer) {
			std::cerr << "Could not open file stream in serializeVector - " << filename << std::endl;
			return writer;
		}
		writeVariable(source.size(), writer);
		writeVector(source, writer);
		return writer;
	};
}