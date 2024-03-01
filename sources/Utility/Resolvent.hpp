#pragma once
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include "OutputConvenience.hpp"
#include "GramSchmidt.hpp"
#include <type_traits>
#include <Eigen/Dense>
#include <cmath>
#include <optional>

namespace Utility {
	using std::abs;

	template <typename T>
	struct ResolventData {
		std::vector<T> a_i;
		std::vector<T> b_i;
	};

	template<class RealType>
	struct ResolventDataWrapper {
		typedef ResolventData<RealType> resolvent_data;
		std::vector<resolvent_data> data;
		size_t noEigenvalueChangeAt{};

		void push_back(resolvent_data&& data_point) {
			data.push_back(std::move(data_point));
		};
		void push_back(const resolvent_data& data_point) {
			data.push_back(data_point);
		};
		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void writeDataToFile(const std::string& filename, const std::optional<std::vector<std::string>>& comments = std::nullopt) const
		{
			std::cout << "Total Lanczos iterations: " << data[0].a_i.size() << "   Point of no change at: " << noEigenvalueChangeAt << std::endl;
			for (const auto& res_data : data) {
				if (checkDataForNaN(res_data.a_i)) std::cerr << "Resolvent a_i" << std::endl;
				if (checkDataForNaN(res_data.b_i)) std::cerr << "Resolvent b_i" << std::endl;
			}
			if (comments.has_value()) {
				saveData(data, filename + ".dat.gz", comments.value());
			}
			else {
				saveData(data, filename + ".dat.gz");
			}
		};
	};

	template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const ResolventData<T>& data)
	{
		for (const auto& elem : data.a_i) {
			os << elem << " ";
		}
		os << "\n";
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
	template <class RealType, bool isComplex>
	class Resolvent
	{
	private:
		using ComputationType = std::conditional_t<isComplex, std::complex<RealType>, RealType>;
		using error_type = std::conditional_t< sizeof(RealType) >= sizeof(double), double, RealType >;

		typedef Eigen::Matrix<ComputationType, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
		typedef Eigen::Vector<ComputationType, Eigen::Dynamic> vector_t;
		typedef Eigen::Vector<RealType, Eigen::Dynamic> RealVector;

		typedef ResolventData<RealType> resolvent_data;

		vector_t startingState;
		ResolventDataWrapper<RealType> data;
	public:
		// Sets the starting state
		inline void setStartingState(const vector_t& state) {
			this->startingState = state;
		};
		const vector_t& getStartingState() const {
			return this->startingState;
		}
		Resolvent(const vector_t& _StargingState) : startingState(_StargingState) {};
		Resolvent() {};

		// Computes the resolvent's parameters a_i and b_i
		// Symplectic needs to be atleast positive semidefinite!
		void compute(const matrix_t& toSolve, const matrix_t& symplectic, int maxIter, error_type errorMargin = 1e-10)
		{
			size_t matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t first = vector_t::Zero(matrixSize); // corresponds to |q_0>
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			ComputationType norm_buffer = second.dot(symplectic * second);
			if constexpr (isComplex) {
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			}
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);

			RealVector eigenDelta(1);
			RealVector eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_t> diagonalize;
			RealVector diagonal; //stores the diagonal elements in a vector
			RealType oldEigenValue{};
			RealType newEigenValue{};
			size_t position{};
			size_t iterNum{};
			bool goOn = true;
			vector_t buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(symplectic * buffer);
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
					deltas.push_back(norm_buffer.real());
				}
				else {
					deltas.push_back(norm_buffer);
				}

				currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				}
				gammas.push_back(abs(norm_buffer));
				basisVectors.push_back(currentSolution / gammas.back());

				++iterNum;

				// construct tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues();

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
						if (!data.noEigenvalueChangeAt) data.noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (long i = 0; i < deltas.size(); i++)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1]);
			}
			data.push_back(std::move(res));
		};

		// Computes the resolvent for a Hermitian problem (i.e. the symplectic matrix is the identity)
		void compute(const matrix_t& toSolve, int maxIter, error_type errorMargin = 1e-10)
		{
			const size_t matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t first = vector_t::Zero(matrixSize); // corresponds to |q_0>
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);

			RealVector eigenDelta(1);
			RealVector eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_t> diagonalize;
			RealVector diagonal; //stores the diagonal elements in a vector
			RealType oldEigenValue{};
			RealType newEigenValue{};
			size_t position{};
			size_t iterNum{};
			bool goOn = true;
			vector_t buffer;

			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					deltas.push_back(basisVectors.back().dot(buffer).real());
				}
				else {
					deltas.push_back(basisVectors.back().dot(buffer));
				}
				currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				gammas.push_back(currentSolution.norm());
				basisVectors.push_back(currentSolution / gammas.back());
				++iterNum;

				// construct the tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues();

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
						if (!data.noEigenvalueChangeAt) data.noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (size_t i = 0U; i < deltas.size(); ++i)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1] * gammas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));
		};

		// Computes the resolvent directly from M and N. This might be more stable for complex matrices
		void computeFromNM(const matrix_t& toSolve, const matrix_t& symplectic, const matrix_t& N, int maxIter, error_type errorMargin = 1e-10)
		{
			auto matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t first = vector_t::Zero(matrixSize); // corresponds to |q_0>
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			ComputationType norm_buffer = second.dot(symplectic * second);
			if constexpr (isComplex) {
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			}
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(first);
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);

			RealVector eigenDelta(1);
			RealVector eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_t> diagonalize;
			RealVector diagonal; //stores the diagonal elements in a vector
			RealType oldEigenValue{};
			RealType newEigenValue{};
			size_t position{};
			size_t  iterNum{};
			bool goOn = true;
			vector_t buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				norm_buffer = basisVectors.back().dot(N * basisVectors.back());
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
					deltas.push_back(norm_buffer.real());
				}
				else {
					deltas.push_back(norm_buffer);
				}

				currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				}
				gammas.push_back(abs(norm_buffer));
				basisVectors.push_back(currentSolution / gammas.back());

				++iterNum;

				// construct tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues();

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
						if (!data.noEigenvalueChangeAt) data.noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (size_t i = 0U; i < deltas.size(); ++i)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1]);
			}
			data.push_back(std::move(res));
		};

		// Computes the resolvent for a Hermitian problem (i.e. the symplectic matrix is the identity)
		// Additionally, this function orthogonalizes the Krylov basis each step
		void computeWithReorthogonalization(const matrix_t& toSolve, int maxIter, error_type errorMargin = 1e-10)
		{
			const size_t matrixSize = toSolve.rows();

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			vector_t currentSolution(matrixSize); // corresponds to |q_(i+1)>
			// First filling
			std::vector<vector_t> basisVectors;
			vector_t second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basisVectors.push_back(second);

			std::vector<RealType> deltas, gammas;
			gammas.push_back(1);

			RealVector eigenDelta(1);
			RealVector eigenGamma(1);

			Eigen::SelfAdjointEigenSolver<matrix_t> diagonalize;
			RealVector diagonal; //stores the diagonal elements in a vector
			RealType oldEigenValue{};
			RealType newEigenValue{};
			size_t position{};
			size_t iterNum{};
			bool goOn = true;
			vector_t buffer;

			while (goOn) {
				// algorithm
				buffer = toSolve * basisVectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					deltas.push_back(basisVectors.back().dot(buffer).real());
				}
				else {
					deltas.push_back(basisVectors.back().dot(buffer));
				}
				if (iterNum > 0U) {
					currentSolution = (buffer - (deltas.back() * basisVectors.back())) - (gammas.back() * basisVectors[iterNum]);
				}
				else {
					currentSolution = (buffer - (deltas.back() * basisVectors.back()));
				}
				GramSchmidt<ComputationType>::orthogonalizeSingleVector(currentSolution, basisVectors);
				gammas.push_back(currentSolution.norm());
				basisVectors.push_back(currentSolution / gammas.back());
				++iterNum;

				// construct the tridiagonal matrix, diagonalize it and find the lowest eigenvalue
				eigenDelta.conservativeResize(iterNum);
				eigenGamma.conservativeResize(iterNum - 1);
				eigenDelta(iterNum - 1) = deltas[iterNum - 1];
				if (iterNum > 1) {
					eigenGamma(iterNum - 2) = gammas[iterNum - 1];
				}
				diagonalize.computeFromTridiagonal(eigenDelta, eigenGamma, 0);
				diagonal = diagonalize.eigenvalues();

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
						if (!data.noEigenvalueChangeAt) data.noEigenvalueChangeAt = iterNum;
					}
				}
				oldEigenValue = newEigenValue;
			}
			for (size_t i = 0U; i < deltas.size(); ++i)
			{
				res.a_i.push_back(deltas[i]);
				res.b_i.push_back(gammas[i + 1] * gammas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));
		};

		const ResolventDataWrapper<RealType>& getData() const {
			return data;
		};

		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void writeDataToFile(const std::string& filename) const
		{
			data.writeDataToFile(filename);
		};
	};
}