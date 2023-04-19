#pragma once
#define _USE_MATH_DEFINES
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "OutputConvenience.hpp"

namespace Utility {
	template <typename T>
	struct ResolventData {
		std::vector<T> a_i;
		std::vector<T> b_i;
	};

	template <typename T>
	inline std::ostream& operator<<( std::ostream& os, const ResolventData<T>& data )
	{
		for (const auto& elem : data.a_i) {
			os << elem << " ";
		}
		os << "0 \n";
		for (const auto& elem : data.b_i) {
			os << elem << " ";
		}
		os << "\n";
		
  		return os;
	}

	// choose the floating point precision, i.e. float, double or long double
	template <typename T>
	class Resolvent
	{
	private:
		typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_T;
		typedef Eigen::Vector<T, Eigen::Dynamic> vector_T;
		typedef Eigen::Vector<std::complex<T>, Eigen::Dynamic> vector_cT;
		typedef Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> matrix_cT;

		vector_cT startingState;
		typedef ResolventData<T> resolvent_data;
		std::vector<resolvent_data> data;

		constexpr unsigned int findSmallestValue(const vector_T& diagonal) const {
			int position = 0;
			for (int i = 1; i < diagonal.size(); i++)
			{
				if (diagonal(position) > diagonal(i)) {
					position = i;
				}
			}
			return position;
		};
		size_t noEigenvalueChangeAt;
	public:
		// Sets the starting state
		inline void setStartingState(const vector_cT& state) {
			this->startingState = state;
		};
		Resolvent(const vector_cT& _StargingState) : startingState(_StargingState), noEigenvalueChangeAt(0) {};
		Resolvent() : noEigenvalueChangeAt(0) {};

		// Computes the resolvent's parameters a_i and b_i
		void compute(const matrix_T& toSolve, const matrix_T& symplectic, int maxIter, T errorMargin = 1e-10)
		{
			size_t matrixSize = toSolve.rows();
			matrix_T identity(matrixSize, matrixSize);
			identity.setIdentity();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_T currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_T> basisVectors;
			vector_T first = vector_T::Zero(matrixSize); // corresponds to |q_0>
			vector_T second = this->startingState.real().template cast<T>(); // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.dot(symplectic * second));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<T> deltas, gammas;
			gammas.push_back(1);

			vector_T eigenDelta(1);
			vector_T eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_T> diagonalize;
			vector_T diagonal; //stores the diagonal elements in a vector
			T oldEigenValue = 0, newEigenValue = 0;
			unsigned int position = 0;
			long iterNum = 0;
			bool goOn = true;
			vector_T buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				deltas.push_back(basisVectors.back().dot(symplectic * buffer));
				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				T norm_squared = currentSolution.dot(symplectic * currentSolution);
				assertm(norm_squared > 0, ("Norm in loop is complex!" + std::to_string(norm_squared)));

				gammas.push_back(sqrt(norm_squared));
				basisVectors.push_back(currentSolution / gammas.back());
				iterNum++;

				// construct the tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues().real();

				if (diagonalize.eigenvalues().imag().norm() > 1e-8) {
					std::cerr << "Atleast one eigenvalue is complex!" << std::endl;
				}
				position = findSmallestValue(diagonal);
				newEigenValue = diagonal(position);

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (std::abs(gammas.back()) < 1e-7) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (std::abs(newEigenValue - oldEigenValue) / std::abs(oldEigenValue) < errorMargin) {
						//goOn = false;
						if (!noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (long i = 0; i < deltas.size(); i++)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1] * gammas[i + 1]);
			}
			this->data.push_back(res);
		};

		// Same as compute, but for complex matrices
		void compute_complex(const matrix_cT& toSolve, const matrix_cT& symplectic, int maxIter, T errorMargin = 1e-10)
		{
			size_t matrixSize = toSolve.rows();
			matrix_cT identity(matrixSize, matrixSize);
			identity.setIdentity();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_cT currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_cT> basisVectors;
			vector_cT first = vector_cT::Zero(matrixSize); // corresponds to |q_0>
			vector_cT second = this->startingState.template cast<std::complex<T>>(); // corresponds to |q_1>
			resolvent_data res;
			auto norm_buffer = second.dot(symplectic * second);
			assertm(std::abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			res.b_i.push_back(std::abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<T> deltas, gammas;
			gammas.push_back(1);

			vector_T eigenDelta(1);
			vector_T eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_cT> diagonalize;
			vector_T diagonal; //stores the diagonal elements in a vector
			T oldEigenValue = 0, newEigenValue = 0;
			unsigned int position = 0;
			long iterNum = 0;
			bool goOn = true;
			vector_cT buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(symplectic * buffer);
				assertm(std::abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
				deltas.push_back(norm_buffer.real());

				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				assertm(std::abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				gammas.push_back(std::abs(norm_buffer));
				basisVectors.push_back(currentSolution / gammas.back());

				iterNum++;

				// construct tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues().real();

				if (diagonalize.eigenvalues().imag().norm() > 1e-8) {
					std::cerr << "Atleast one eigenvalue is complex!" << std::endl;
				}
				position = findSmallestValue(diagonal);
				newEigenValue = diagonal(position);

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (std::abs(gammas.back()) < 1e-8) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (std::abs(newEigenValue - oldEigenValue) / std::abs(oldEigenValue) < errorMargin) {
						//goOn = false;
						if (!noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (long i = 0; i < deltas.size(); i++)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1]);
			}
			this->data.push_back(res);
		};

		// Computes the resolvent for a Hermitian problem (i.e. the symplectic matrix is the identity)
		void compute(const matrix_T& toSolve, int maxIter, T errorMargin = 1e-10)
		{
			size_t matrixSize = toSolve.rows();
			matrix_T identity(matrixSize, matrixSize);
			identity.setIdentity();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_T currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_T> basisVectors;
			vector_T first = vector_T::Zero(matrixSize); // corresponds to |q_0>
			vector_T second = this->startingState.real(); // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.dot(second));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<T> deltas, gammas;
			gammas.push_back(1);

			vector_T eigenDelta(1);
			vector_T eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_T> diagonalize;
			vector_T diagonal; //stores the diagonal elements in a vector
			T oldEigenValue = 0, newEigenValue = 0;
			unsigned int position = 0;
			long iterNum = 0;
			bool goOn = true;
			vector_T buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				deltas.push_back(basisVectors.back().dot(buffer));
				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				T norm_squared = currentSolution.dot(currentSolution);
				assertm(norm_squared > 0, ("Norm in loop is complex!" + std::to_string(norm_squared)));

				gammas.push_back(sqrt(norm_squared));
				basisVectors.push_back(currentSolution / gammas.back());
				iterNum++;

				// construct the tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues().real();

				if (diagonalize.eigenvalues().imag().norm() > 1e-8) {
					std::cerr << "Atleast one eigenvalue is complex!" << std::endl;
				}
				position = findSmallestValue(diagonal);
				newEigenValue = diagonal(position);

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (std::abs(gammas.back()) < 1e-7) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (std::abs(newEigenValue - oldEigenValue) / std::abs(oldEigenValue) < errorMargin) {
						//goOn = false;
						if (!noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (long i = 0; i < deltas.size(); i++)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1] * gammas[i + 1]);
			}
			this->data.push_back(res);
		};

		// Computes the resolvent directly from M and N. This might be more stable for complex matrices
		void computeFromNM(const matrix_cT& toSolve, const matrix_cT& symplectic, const matrix_cT& N, int maxIter, T errorMargin = 1e-10)
		{
			auto matrixSize = toSolve.rows();
			matrix_cT identity(matrixSize, matrixSize);
			identity.setIdentity();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_cT currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_cT> basisVectors;
			vector_cT first = vector_cT::Zero(matrixSize); // corresponds to |q_0>
			vector_cT second = this->startingState.template cast<std::complex<T>>(); // corresponds to |q_1>
			resolvent_data res;
			auto norm_buffer = second.dot(symplectic * second);
			assertm(std::abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			res.b_i.push_back(std::abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<T> deltas, gammas;
			gammas.push_back(1);

			vector_T eigenDelta(1);
			vector_T eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_cT> diagonalize;
			vector_T diagonal; //stores the diagonal elements in a vector
			T oldEigenValue = 0, newEigenValue = 0;
			unsigned int position = 0;
			long iterNum = 0;
			bool goOn = true;
			vector_cT buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(N * basisVectors.back());
				assertm(std::abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
				deltas.push_back(norm_buffer.real());

				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				assertm(std::abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				gammas.push_back(std::abs(norm_buffer));
				basisVectors.push_back(currentSolution / gammas.back());

				iterNum++;

				// construct tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues().real();

				if (diagonalize.eigenvalues().imag().norm() > 1e-8) {
					std::cerr << "Atleast one eigenvalue is complex!" << std::endl;
				}
				position = findSmallestValue(diagonal);
				newEigenValue = diagonal(position);

				// breaking conditions
				if (iterNum >= maxIter) {
					goOn = false;
				}
				if (iterNum >= toSolve.rows()) {
					goOn = false;
				}
				if (std::abs(gammas.back()) < 1e-8) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (std::abs(newEigenValue - oldEigenValue) / std::abs(oldEigenValue) < errorMargin) {
						//goOn = false;
						if (!noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (long i = 0; i < deltas.size(); i++)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1]);
			}
			this->data.push_back(res);
		};

		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void writeDataToFile(const std::string& filename) const
		{
			saveData_boost(this->data, filename + ".dat.gz");
			saveData(this->data, filename + ".txt");
			return;

			std::cout << "Total Lanczos iterations: " << this->data[0].a_i.size() << "   Point of no change at: " << noEigenvalueChangeAt << std::endl;
			std::ofstream out(filename);
			if (out.is_open()) {
				out << "# " << Utility::time_stamp() << "\n# No Change Iteration: " << noEigenvalueChangeAt << "\n\n" << std::scientific << std::setprecision(12);
				for (const auto& set : this->data) {
					for (const auto& elem : set.a_i) {
						out << elem << " ";
					}
					out << "0 \n";
					for (const auto& elem : set.b_i) {
						out << elem << " ";
					}
					out << "\n";
				}
			}
			else {
				std::cerr << "Could not open output filestream! " << filename << std::endl;
			}
		};
	};
}