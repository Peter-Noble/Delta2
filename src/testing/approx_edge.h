#pragma once

#include "../common/edge.h"
#include "unit_test.h"
#include "approx_vector.h"

// The matcher class
template<typename real>
class EdgeApproxEqual : public Catch::MatcherBase<Delta2::common::Edge<real>> {
    Delta2::common::Edge<real> m_comp;
    real e;
public:
    EdgeApproxEqual(Delta2::common::Edge<real> comp) : m_comp(comp) {
        e = 1e-3;
    }

    EdgeApproxEqual(Delta2::common::Edge<real> comp, real eps) : m_comp(comp), e(eps) { }

    EdgeApproxEqual& margin(real m) {
        e = m;
        return *this;
    }

    // Performs the test for this matcher
    bool match(Delta2::common::Edge<real> const& comp) const override {
        VectorApproxEqual<real> a(comp.A, e);
        VectorApproxEqual<real> b(comp.B, e);
        return a.match(m_comp.A) && b.match(m_comp.B);
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "is approximately equal to Edge: {{ " << m_comp.A.x() << ", " << m_comp.A.y() << ", " << m_comp.A.z() << "}, {" << m_comp.B.x() << ", " << m_comp.B.y() << ", " << m_comp.B.z() << "}}";
        return ss.str();
    }
};

// The builder function
template<typename real>
EdgeApproxEqual<real> EdgeEqual(Delta2::common::Edge<real> comp) {
    return EdgeApproxEqual<real>(comp);
}
