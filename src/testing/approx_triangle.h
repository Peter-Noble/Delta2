#pragma once

#include "../common/triangle.h"
#include "unit_test.h"
#include "approx_vector.h"

// The matcher class
template<typename real>
class TriangleApproxEqual : public Catch::MatcherBase<Delta2::common::Triangle<real>> {
    Delta2::common::Triangle<real> m_comp;
    real e;
public:
    TriangleApproxEqual(Delta2::common::Triangle<real> comp) : m_comp(comp) {
        e = 1e-3;
    }

    TriangleApproxEqual(Delta2::common::Triangle<real> comp, real eps) : m_comp(comp), e(eps) { }

    TriangleApproxEqual& margin(real m) {
        e = m;
        return *this;
    }

    // Performs the test for this matcher
    bool match(Delta2::common::Triangle<real> const& comp) const override {
        VectorApproxEqual<real> a(comp.A, e);
        VectorApproxEqual<real> b(comp.B, e);
        VectorApproxEqual<real> c(comp.C, e);
        return a.match(m_comp.A) && b.match(m_comp.B) && c.match(m_comp.C);
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "is approximately equal to triangle: {{ " << m_comp[0].x() << ", " << m_comp[0].y() << ", " << m_comp[0].z() << "}, {" << m_comp[1].x() << ", " << m_comp[1].y() << ", " << m_comp[1].z() << "}, {" << m_comp[2].x() << ", " << m_comp[2].y() << ", " << m_comp[2].z() << "}}";
        return ss.str();
    }
};

// The builder function
template<typename real>
TriangleApproxEqual<real> TriangleEqual(Delta2::common::Triangle<real> comp) {
    return TriangleApproxEqual<real>(comp);
}
