#pragma once

#include "unit_test.h"
#include <Eigen/Dense>

// The matcher class
template<typename real, int size>
class MatrixApproxEqual : public Catch::MatcherBase<Eigen::Matrix<real, size, size>> {
    Eigen::Matrix<real, size, size> m_comp;
    real e;
public:
    MatrixApproxEqual(Eigen::Matrix<real, size, size> comp) : m_comp(comp) {
        e = 1e-3;
    }

    MatrixApproxEqual(Eigen::Matrix<real, size, size> comp, real eps) : m_comp(comp), e(eps) { }

    MatrixApproxEqual& margin(real m) {
        e = m;
        return *this;
    }

    // Performs the test for this matcher
    bool match(Eigen::Matrix<real, size, size> const& comp) const override {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (std::abs(m_comp(i, j) - comp(i, j)) > e) {
                    return false;
                }
            }
        }
        return true;
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "is approximately equal to Matrix: " << std::endl;
        ss << m_comp << std::endl;
        return ss.str();
    }
};

// The builder function
template<typename real, int size>
MatrixApproxEqual<real, size> MatrixEqual(Eigen::Matrix<real, size, size> comp) {
    return MatrixApproxEqual<real, size>(comp);
}
