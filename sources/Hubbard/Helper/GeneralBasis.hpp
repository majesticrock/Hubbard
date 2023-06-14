#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
    class GeneralBasis : public ModeHelper
    {
    protected:
        Matrix_L M, N;
    public:
        virtual void fillMatrices() override;
    };
}