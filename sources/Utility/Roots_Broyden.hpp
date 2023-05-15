#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <limits>

namespace Utility {
	namespace NumericalSolver {
		namespace Roots {
			template<typename data_type, int t_vector_size>
			class Broyden {
				using MatrixType = Eigen::Matrix<data_type, t_vector_size, t_vector_size>;
				using VectorType = Eigen::Vector<data_type, t_vector_size>;

				static void estimate_jacobian(MatrixType& J_new, const MatrixType& J_old,
					const VectorType& delta_x, const VectorType& delta_F)
				{
					J_new = (delta_x - J_old * delta_F) * (delta_x.transpose() * J_old);
					J_new /= (delta_x.dot(J_old * delta_F));
					J_new += J_old;
				};
			public:
				static bool compute(std::function<void(const VectorType&, VectorType&)>& func,
					VectorType& x0, const int MAX_ITER = 200)
				{
					size_t DIM = x0.rows();
					// You may play around with EPS_X and EPS_F to your desire
					// EPS_X is the minimum distance between x_i and x_i+1
					// EPS_F is the minimum f(x)
					constexpr double EPS_F = std::numeric_limits<double>::denorm_min();
					constexpr double EPS_X = 1e-12;
					double diff_x = 100, diff_F = 100;
					int iter_num = 0;

					VectorType F_old = Eigen::VectorXd::Zero(DIM),
						F_new = Eigen::VectorXd::Zero(DIM),
						delta_x = Eigen::VectorXd::Zero(DIM),
						delta_F = Eigen::VectorXd::Zero(DIM);
					MatrixType J_old = Eigen::MatrixXd::Zero(DIM, DIM),
						J_new = Eigen::MatrixXd::Identity(DIM, DIM);
					func(x0, F_new);

					while (diff_x > EPS_X && diff_F > EPS_F && iter_num++ <= MAX_ITER && F_new.squaredNorm() > EPS_F) {
						delta_x = -J_new * F_new;
						x0 += delta_x;
						diff_x = delta_x.norm();

						F_old = F_new;
						func(x0, F_new);
						delta_F = F_new - F_old;
						diff_F = delta_F.norm();

						J_old = J_new;
						estimate_jacobian(J_new, J_old, delta_x, delta_F);
					}
					// This method returns true if convergence is achieved, in this case if |F(x_final)| < 1e-10
					return (F_new.norm() < 1e-10);
				};
			};
		}
	}
}