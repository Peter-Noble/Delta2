#pragma once
#include "unit_test.h"

// The matcher class
template<typename real>
class VectorApproxEqual : public Catch::MatcherBase<Eigen::Matrix<real, 3, 1>> {
    Eigen::Matrix<real, 3, 1> m_comp;
    real e;
public:
    VectorApproxEqual(Eigen::Matrix<real, 3, 1> comp) : m_comp(comp) {
        // e = std::numeric_limits<real>::epsilon() * 100;
        e = 1e-3;
    }

    VectorApproxEqual& margin(real m) {
        e = m;
        return *this;
    }

    VectorApproxEqual(Eigen::Matrix<real, 3, 1> comp, real eps) : m_comp(comp), e(eps) { }

    // Performs the test for this matcher
    bool match(Eigen::Matrix<real, 3, 1> const& comp) const override {
        Eigen::Matrix<real, 3, 1> diff = m_comp - comp;
        return fabs(diff[0]) < e && fabs(diff[1]) < e && fabs(diff[2]) < e;
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "is approximately equal to vector: { " << m_comp[0] << ", " << m_comp[1] << ", " << m_comp[2] << "}";
        return ss.str();
    }
};

// The builder function
template<typename real>
VectorApproxEqual<real> VectorEqual(Eigen::Matrix<real, 3, 1> comp) {
    return VectorApproxEqual<real>(comp);
}
