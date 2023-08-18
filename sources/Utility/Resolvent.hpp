#pragma once
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include "OutputConvenience.hpp"
#include <Eigen/Dense>

namespace Utility {
	template <typename T>
	struct ResolventData {
		std::vector<T> a_i;
		std::vector<T> b_i;
	};

	template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const ResolventData<T>& data)
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
	template <typename RealType>
	inline size_t findSmallestValue(const Eigen::Vector<RealType, -1>& diagonal) {
		size_t position = 0;
		for (size_t i = 1U; i < diagonal.size(); ++i)
		{
			if (diagonal(position) > diagonal(i)) {
				position = i;
			}
		}
		return position;
	};

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
		size_t noEigenvalueChangeAt;
	public:
		// Sets the starting state
		inline void setStartingState(const vector_cT& state) {
			this->startingState = state;
		};
		const vector_cT& getStartingState() const {
			return this->startingState;
		}
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
			size_t position = 0U;
			size_t iterNum = 0U;
			bool goOn = true;
			vector_T buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				deltas.push_back(basisVectors.back().dot(symplectic * buffer));
				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				T norm_squared = currentSolution.dot(symplectic * currentSolution);
				assertm(norm_squared > 0, ("Norm in loop is complex!" + to_string(norm_squared)));

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
				if (abs(gammas.back()) < 1e-7) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (abs(newEigenValue - oldEigenValue) / abs(oldEigenValue) < errorMargin) {
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
			data.push_back(res);
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
			assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<T> deltas, gammas;
			gammas.push_back(1);

			vector_T eigenDelta(1);
			vector_T eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_cT> diagonalize;
			vector_T diagonal; //stores the diagonal elements in a vector
			T oldEigenValue = 0., newEigenValue = 0.;
			size_t position = 0U;
			size_t iterNum = 0U;
			bool goOn = true;
			vector_cT buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(symplectic * buffer);
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
				deltas.push_back(norm_buffer.real());

				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				gammas.push_back(abs(norm_buffer));
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
				if (abs(gammas.back()) < 1e-8) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (abs(newEigenValue - oldEigenValue) / abs(oldEigenValue) < errorMargin) {
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
			data.push_back(res);
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
			res.b_i.push_back(second.squaredNorm());

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
			size_t position = 0U;
			size_t iterNum = 0U;
			bool goOn = true;
			vector_T buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				deltas.push_back(basisVectors.back().dot(buffer));
				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				T norm_squared = currentSolution.squaredNorm();

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
				if (abs(gammas.back()) < 1e-7) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (abs(newEigenValue - oldEigenValue) / abs(oldEigenValue) < errorMargin) {
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
			data.push_back(res);
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
			assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			res.b_i.push_back(abs(norm_buffer));

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
			size_t position = 0U;
			size_t iterNum = 0U;
			bool goOn = true;
			vector_cT buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(N * basisVectors.back());
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
				deltas.push_back(norm_buffer.real());

				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				gammas.push_back(abs(norm_buffer));
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
				if (abs(gammas.back()) < 1e-8) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (abs(newEigenValue - oldEigenValue) / abs(oldEigenValue) < errorMargin) {
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
			data.push_back(res);
		};

		// Computes the resolvent for a Hermitian problem (i.e. the symplectic matrix is the identity)
		void compute(const matrix_cT& toSolve, int maxIter, T errorMargin = 1e-10)
		{
			size_t matrixSize = toSolve.rows();
			matrix_T identity(matrixSize, matrixSize);
			identity.setIdentity();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_cT currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_cT> basisVectors;
			vector_cT first = vector_cT::Zero(matrixSize); // corresponds to |q_0>
			vector_cT second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

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
			size_t position = 0U;
			size_t iterNum = 0U;
			bool goOn = true;
			vector_cT buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				// This has to be real, as <x|H|x> is always real if H=H^+
				deltas.push_back(basisVectors.back().dot(buffer).real());
				currentSolution = (buffer - ((deltas.back() * identity) * basisVectors.back())) - (gammas.back() * basisVectors.at(iterNum));
				T norm_squared = currentSolution.squaredNorm();

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
				if (abs(gammas.back()) < 1e-7) {
					goOn = false;
				}
				if (oldEigenValue != 0.0) {
					if (abs(newEigenValue - oldEigenValue) / abs(oldEigenValue) < errorMargin) {
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
			data.push_back(res);
		};

		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void writeDataToFile(const std::string& filename) const
		{
			std::cout << "Total Lanczos iterations: " << data[0].a_i.size() << "   Point of no change at: " << noEigenvalueChangeAt << std::endl;
			saveData(data, filename + ".dat.gz");
		};
	};
}