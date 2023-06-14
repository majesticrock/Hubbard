#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
    class XPModes : public ModeHelper
    {
    protected:
        Matrix_L K_plus, K_minus, L;
    public:
        virtual void fillMatrices() override;
    };
}