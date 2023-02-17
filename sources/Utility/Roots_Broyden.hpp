#pragma once
#include <Eigen/Dense>

namespace Utility {
	namespace Roots {
		namespace Broyden {
			void estimate_jacobian(Eigen::MatrixXd& J_new, const Eigen::MatrixXd& J_old, const Eigen::VectorXd& delta_x, const Eigen::VectorXd& delta_F);

			bool compute(std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)>& func, Eigen::VectorXd& x0, const int MAX_ITER = 200);
		}
	}
}