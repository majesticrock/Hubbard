#pragma once

#include <Eigen/Dense>
#include <iostream>

namespace Utility {
	namespace NumericalSolver {
		template <typename T>
		class Lanczos {
		private:
			typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_T;
			typedef Eigen::Vector<T, Eigen::Dynamic> vector_T;

			vector_T startingState;
			vector_T eigenDelta;
			vector_T eigenGamma;
			matrix_T basis_trafo;
			vector_T eigenValues;
		public:
			// Sets the starting state
			inline void setStartingState(const vector_T& state) {
				this->startingState = state;
				this->startingState.normalize();
			};

			Lanczos(const vector_T& _StargingState) : startingState(_StargingState) {
				this->startingState.normalize();
			};
			Lanczos() {};

			const vector_T& getEigenValues() const {
				return eigenValues;
			}

			matrix_T pseudoInverse(const T epsilon = 10e-12) const {
				Eigen::SelfAdjointEigenSolver<matrix_T> solver;
				solver.computeFromTridiagonal(eigenDelta, eigenGamma);

				vector_T& evs = const_cast<vector_T&>(solver.eigenvalues());
				for (size_t i = 0; i < evs.size(); i++)
				{
					if (abs(evs(i)) < epsilon) {
						evs(i) = 0;
					}
					else {
						evs(i) = 1 / evs(i);
					}
				}

				matrix_T complete_trafo = basis_trafo * solver.eigenvalues();
				return complete_trafo * evs.asDiagonal() * complete_trafo.adjoint();
			}

			matrix_T pseudoInverseSquareRoot(const T epsilon = 10e-12) const {
				Eigen::SelfAdjointEigenSolver<matrix_T> solver;
				solver.computeFromTridiagonal(eigenDelta, eigenGamma);

				vector_T& evs = const_cast<vector_T&>(solver.eigenvalues());
				for (size_t i = 0; i < evs.size(); i++)
				{
					if (abs(evs(i)) < epsilon) {
						evs(i) = 0;
					}
					else {
						evs(i) = 1 / sqrt(evs(i));
					}
				}

				matrix_T complete_trafo = basis_trafo * solver.eigenvalues();
				return complete_trafo * evs.asDiagonal() * complete_trafo.adjoint();
			}

			void compute(const matrix_T& toSolve, int maxIter, T errorMargin = 1e-12)
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
				vector_T second = this->startingState; // corresponds to |q_1>
				second.normalize();

				basisVectors.push_back(first);
				basisVectors.push_back(second);

				std::vector<T> deltas, gammas;
				gammas.push_back(1);

				eigenDelta = vector_T::Zero(1);
				eigenGamma = vector_T::Zero(1);

				vector_T diagonal; //stores the diagonal elements in a vector
				long iterNum = 0;
				bool goOn = true;
				vector_T buffer;
				while (goOn) {
					// algorithm
					buffer = toSolve * basisVectors.back();
					deltas.push_back(basisVectors.back().dot(buffer));
					currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
					T norm_squared = currentSolution.dot(currentSolution);

					gammas.push_back(sqrt(norm_squared));
					basisVectors.push_back(currentSolution / gammas.back());
					iterNum++;

					eigenDelta.conservativeResize(iterNum);
					eigenGamma.conservativeResize(iterNum - 1);
					eigenDelta(iterNum - 1) = deltas[iterNum - 1];
					if (iterNum > 1) {
						eigenGamma(iterNum - 2) = gammas[iterNum - 1];
					}

					// breaking conditions
					if (iterNum >= maxIter) {
						goOn = false;
					}
					if (iterNum >= toSolve.rows()) {
						goOn = false;
					}
					if (abs(gammas.back()) < 1e-7) {
						goOn = false;
					}
				}

				std::cout << "Number of iterations=" << iterNum << "\n";
				Eigen::SelfAdjointEigenSolver<matrix_T> diagonalize;
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, false);
				this->eigenValues = diagonalize.eigenvalues();

				// Create an Eigen::MatrixXd with the desired size
				basis_trafo = matrix_T::Zero(basisVectors.size() - 1, basisVectors[1].size());

				// Copy the data from the std::vector<Eigen::VectorXd> to the Eigen::MatrixXd
				for (int i = 1; i < basisVectors.size(); i++) {
					basis_trafo.row(i - 1) = basisVectors[i];
				}
			};
		};
	}
}