#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
    typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixCL;
	typedef Eigen::Vector<std::complex<double>, Eigen::Dynamic> VectorCL;

    class GeneralBasis : public ModeHelper
    {
    private:
        void fillBlock(int i, int j);
    protected:
        MatrixCL M, N;

        virtual void fillMatrices() override;
    public:
        GeneralBasis(Utility::InputFileReader& input) : ModeHelper(input) { };

        virtual std::unique_ptr<std::vector<Resolvent_L>> computeCollectiveModes(std::vector<std::vector<double>>& reciever) override;
    };
}