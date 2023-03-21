#define _USE_MATH_DEFINES
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include "Resolvent.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "OutputWriter.hpp"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

namespace Utility {
	void Resolvent::compute(const MatrixXd& toSolve, const MatrixXd& symplectic, int maxIter, double errorMargin/*=1e-10*/)
	{
		size_t matrixSize = toSolve.rows();
		MatrixXd identity(matrixSize, matrixSize);
		identity.setIdentity();

		if (toSolve.rows() != toSolve.cols()) {
			std::cerr << "Matrix is not square!" << std::endl;
			throw;
		}

		VectorXd currentSolution(matrixSize); // corresponds to |q_(i+1)>
		// First filling
		std::vector<VectorXd> basisVectors;
		VectorXd first = VectorXd::Zero(matrixSize); // corresponds to |q_0>
		VectorXd second = this->startingState.real(); // corresponds to |q_1>
		resolvent_data res;
		res.b_i.push_back(second.dot(symplectic * second));

		second /= sqrt(res.b_i.back());
		basisVectors.push_back(first);
		basisVectors.push_back(second);

		std::vector<double> deltas, gammas;
		gammas.push_back(1);

		Eigen::VectorXd eigenDelta(1);
		Eigen::VectorXd eigenGamma(1);

		Eigen::SelfAdjointEigenSolver<MatrixXd> diagonalize;
		VectorXd diagonal; //stores the diagonal elements in a vector
		double oldEigenValue = 0, newEigenValue = 0;
		unsigned int position = 0;
		long iterNum = 0;
		bool goOn = true;
		VectorXd buffer;
		while (goOn) {
			// algorithm
			buffer = toSolve * basisVectors.back();
			deltas.push_back(basisVectors.back().dot(symplectic * buffer));
			currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
			double norm_squared = currentSolution.dot(symplectic * currentSolution);
			assertm(norm_squared > 0, "Norm in loop is complex!");

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
					if(noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
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
	}

	void Resolvent::compute(const MatrixXcd& toSolve, const MatrixXcd& symplectic, int maxIter, double errorMargin/*=1e-10*/)
	{
		size_t matrixSize = toSolve.rows();
		MatrixXcd identity(matrixSize, matrixSize);
		identity.setIdentity();

		if (toSolve.rows() != toSolve.cols()) {
			std::cerr << "Matrix is not square!" << std::endl;
			throw;
		}

		VectorXcd currentSolution(matrixSize); // corresponds to |q_(i+1)>
		// First filling
		std::vector<VectorXcd> basisVectors;
		VectorXcd first = VectorXcd::Zero(matrixSize); // corresponds to |q_0>
		VectorXcd second = this->startingState; // corresponds to |q_1>
		resolvent_data res;
		auto norm_buffer = second.dot(symplectic * second);
		assertm(std::abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
		res.b_i.push_back(std::abs(norm_buffer));

		second /= sqrt(res.b_i.back());
		basisVectors.push_back(first);
		basisVectors.push_back(second);

		std::vector<double> deltas, gammas;
		gammas.push_back(1);

		Eigen::VectorXd eigenDelta(1);
		Eigen::VectorXd eigenGamma(1);

		Eigen::SelfAdjointEigenSolver<MatrixXd> diagonalize;
		VectorXd diagonal; //stores the diagonal elements in a vector
		double oldEigenValue = 0, newEigenValue = 0;
		unsigned int position = 0;
		long iterNum = 0;
		bool goOn = true;
		VectorXcd buffer;
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
					if(noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
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
	}

	void Resolvent::computeFromNM(const Eigen::MatrixXcd& toSolve, const Eigen::MatrixXcd& symplectic, const Eigen::MatrixXcd& N, int maxIter, double errorMargin)
	{
		auto matrixSize = toSolve.rows();
		MatrixXcd identity(matrixSize, matrixSize);
		identity.setIdentity();

		if (toSolve.rows() != toSolve.cols()) {
			std::cerr << "Matrix is not square!" << std::endl;
			throw;
		}

		VectorXcd currentSolution(matrixSize); // corresponds to |q_(i+1)>
		// First filling
		std::vector<VectorXcd> basisVectors;
		VectorXcd first = VectorXcd::Zero(matrixSize); // corresponds to |q_0>
		VectorXcd second = this->startingState; // corresponds to |q_1>
		resolvent_data res;
		auto norm_buffer = second.dot(symplectic * second);
		assertm(std::abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
		res.b_i.push_back(std::abs(norm_buffer));

		second /= sqrt(res.b_i.back());
		basisVectors.push_back(first);
		basisVectors.push_back(second);

		std::vector<double> deltas, gammas;
		gammas.push_back(1);

		Eigen::VectorXd eigenDelta(1);
		Eigen::VectorXd eigenGamma(1);

		Eigen::SelfAdjointEigenSolver<MatrixXd> diagonalize;
		VectorXd diagonal; //stores the diagonal elements in a vector
		double oldEigenValue = 0, newEigenValue = 0;
		unsigned int position = 0;
		long iterNum = 0;
		bool goOn = true;
		VectorXcd buffer;
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
					if(noEigenvalueChangeAt) noEigenvalueChangeAt = iterNum;
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
	}

	void Resolvent::writeDataToFile(std::string filename) const
	{
		std::cout << "Total Lanczos iterations: " << this->data[0].a_i.size() << "   Point of no change at: " << noEigenvalueChangeAt << std::endl;
		std::ofstream out(filename);
		if (out.is_open()) {
			out << "# " << Utility::time_stamp() << "\n# No Change Iteration: " << noEigenvalueChangeAt << "\n\n" << std::scientific << std::setprecision(12);
			for (auto& set : this->data) {
				for (auto& elem : set.a_i) {
					out << elem << " ";
				}
				out << "0 \n";
				for (auto& elem : set.b_i) {
					out << elem << " ";
				}
				out << "\n";
			}
		}
		else {
			std::cerr << "Could not open output filestream! " << filename << std::endl;
		}
	}
}