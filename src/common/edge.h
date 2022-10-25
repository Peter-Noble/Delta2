#pragma once

#include <Eigen/Dense>
#include "basic_utils.h"

namespace Delta2 {
    namespace common {
        template<typename real>
        class Edge {
            public:
                Edge(Eigen::Vector<real, 3> a, Eigen::Vector<real, 3> b) : A(a), B(b) {};
                Edge() {};
                Eigen::Vector<real, 3> A, B;

                template<typename ToT>
                Edge<ToT> cast() const {
                    return Edge<ToT>(A.template cast<ToT>(), B.template cast<ToT>());
                }
                real dist(const Edge<real>& E) const {
                    Eigen::Vector<real, 3> segDC = E.B-E.A;
                    real lineDirSqrMag = segDC.dot(segDC);
                    Eigen::Vector<real, 3> segAC = A-E.A;
                    Eigen::Vector<real, 3> segBC = B-E.A;
                    Eigen::Vector<real, 3> inPlaneA = A-((segAC.dot(segDC)/lineDirSqrMag)*segDC);
                    Eigen::Vector<real, 3> inPlaneB = B-((segBC.dot(segDC)/lineDirSqrMag)*segDC);
                    Eigen::Vector<real, 3> inPlaneBA = inPlaneB-inPlaneA;
                    real t = (E.A-inPlaneA).dot(inPlaneBA)/inPlaneBA.dot(inPlaneBA);
                    t = (inPlaneA != inPlaneB) ? t : 0.0; // Zero's t if parallel
                    Eigen::Vector<real, 3> segABtoLineCD = lerp(clamp01(t));

                    Eigen::Vector<real, 3> segCDtoSegAB = E.constrainToSegment(segABtoLineCD);
                    Eigen::Vector<real, 3> segABtoSegCD = constrainToSegment(segCDtoSegAB);

                    Eigen::Vector<real, 3> Q_out = segCDtoSegAB;
                    Eigen::Vector<real, 3> P_out = segABtoSegCD;

                    return (Q_out - P_out).norm();
                }

                real closest(const Edge<real>& E, Eigen::Vector<real, 3> &P_out, Eigen::Vector<real, 3> &Q_out) const {
                    Eigen::Vector<real, 3> segDC = E.B-E.A;
                    real lineDirSqrMag = segDC.dot(segDC);
                    Eigen::Vector<real, 3> segAC = A-E.A;
                    Eigen::Vector<real, 3> segBC = B-E.A;
                    Eigen::Vector<real, 3> inPlaneA = A-((segAC.dot(segDC)/lineDirSqrMag)*segDC);
                    Eigen::Vector<real, 3> inPlaneB = B-((segBC.dot(segDC)/lineDirSqrMag)*segDC);
                    Eigen::Vector<real, 3> inPlaneBA = inPlaneB-inPlaneA;
                    real t = (E.A-inPlaneA).dot(inPlaneBA)/inPlaneBA.dot(inPlaneBA);
                    t = (inPlaneA != inPlaneB) ? t : 0.0; // Zero's t if parallel
                    Eigen::Vector<real, 3> segABtoLineCD = lerp(clamp01(t));

                    Eigen::Vector<real, 3> segCDtoSegAB = E.constrainToSegment(segABtoLineCD);
                    Eigen::Vector<real, 3> segABtoSegCD = constrainToSegment(segCDtoSegAB);

                    Q_out = segCDtoSegAB;
                    P_out = segABtoSegCD;

                    return (Q_out - P_out).norm();
                }

                Eigen::Vector<real, 3> constrainToSegment(const Eigen::Vector<real, 3>& position) const {
                    Eigen::Vector<real, 3> ba = B - A;
                    real t = (position - A).dot(ba) / ba.dot(ba);
                    return lerp(clamp01(t));
                }

                Eigen::Vector<real, 3> lerp(real t) const {
                    return Delta2::common::lerp(A, B, t);
                }

                real project(Eigen::Vector<real, 3> V) const {
                    Eigen::Vector<real, 3> AB = B - A;
                    return (V - A).dot(AB) / AB.dot(AB);
                }

                Eigen::Vector<real, 3> min() const {
                    return Eigen::Vector<real, 3>(std::min(A.x(), B.x()), std::min(A.y(), B.y()), std::min(A.z(), B.z()));
                }

                Eigen::Vector<real, 3> min(const Edge<real>& other) const {
                    Eigen::Vector<real, 3> other_min = other.min();
                    Eigen::Vector<real, 3> this_min = min();
                    return Eigen::Vector<real, 3>(std::min(this_min.x(), other_min.x()), std::min(this_min.y(), other_min.y()), std::min(this_min.z(), other_min.z()));
                }

                Eigen::Vector<real, 3> max() const {
                    return Eigen::Vector<real, 3>(std::max(A.x(), B.x()), std::max(A.y(), B.y()), std::max(A.z(), B.z()));
                }

                Eigen::Vector<real, 3> max(const Edge<real>& other) const {
                    Eigen::Vector<real, 3> other_max = other.max();
                    Eigen::Vector<real, 3> this_max = max();
                    return Eigen::Vector<real, 3>(std::max(this_max.x(), other_max.x()), std::max(this_max.y(), other_max.y()), std::max(this_max.z(), other_max.z()));
                }
        };
        template class Edge<float>;
        template class Edge<double>;
    }
}