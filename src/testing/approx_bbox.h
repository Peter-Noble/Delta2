#pragma once

#include "unit_test.h"
#include "../embree_common/math/bbox.h"

// The matcher class
class BBoxApproxEqual : public Catch::MatcherBase<embree::BBox3fa> {
    embree::BBox3fa m_comp;
    float e;
public:
    BBoxApproxEqual(embree::BBox3fa comp) : m_comp(comp) {
        e = 1e-3;
    }

    BBoxApproxEqual(embree::BBox3fa comp, float eps) : m_comp(comp), e(eps) { }

    BBoxApproxEqual& margin(float m) {
        e = m;
        return *this;
    }

    // Performs the test for this matcher
    bool match(embree::BBox3fa const& comp) const override {
        bool m = true;
        m &= std::abs(comp.lower.x - m_comp.lower.x) < e;
        m &= std::abs(comp.lower.y - m_comp.lower.y) < e;
        m &= std::abs(comp.lower.z - m_comp.lower.z) < e;
        m &= std::abs(comp.upper.x - m_comp.upper.x) < e;
        m &= std::abs(comp.upper.y - m_comp.upper.y) < e;
        m &= std::abs(comp.upper.z - m_comp.upper.z) < e;
        return m;
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "is approximately equal to BBox: " << m_comp;
        return ss.str();
    }
};

// The builder function
BBoxApproxEqual BBoxEqual(embree::BBox3fa comp) {
    return BBoxApproxEqual(comp);
}
