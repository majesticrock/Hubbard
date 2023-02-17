#pragma once
#include <Eigen/Dense>
#include <vector>

namespace Utility {
	class Resolvent
	{
	private:
		struct resolvent_data {
			std::vector<double> a_i;
			std::vector<double> b_i;
		};
		Eigen::VectorXcd startingState;
		std::vector<resolvent_data> data;
		constexpr unsigned int findSmallestValue(const Eigen::VectorXd& diagonal) const {
			int position = 0;
			for (int i = 1; i < diagonal.size(); i++)
			{
				if (diagonal(position) > diagonal(i)) {
					position = i;
				}
			}
			return position;
		};
	public:
		// Sets the starting state
		inline void setStartingState(const Eigen::VectorXcd& state) {
			this->startingState = state;
		};
		Resolvent(const Eigen::VectorXcd& _StargingState) : startingState(_StargingState) {};
		Resolvent() {};

		// Computes the resolvent's parameters a_i and b_i
		void compute(const Eigen::MatrixXd& toSolve, const Eigen::MatrixXd& symplectic, int maxIter, double errorMargin = 1e-10);
		void compute(const Eigen::MatrixXcd& toSolve, const Eigen::MatrixXcd& symplectic, int maxIter, double errorMargin = 1e-10);
		void computeFromNM(const Eigen::MatrixXcd& toSolve, const Eigen::MatrixXcd& symplectic, const Eigen::MatrixXcd& N, int maxIter, double errorMargin = 1e-10);
		// Prints the computed data to <filename>
		// Asummes that the data has been computed before...
		void writeDataToFile(std::string filename) const;
	};
}