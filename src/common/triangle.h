#pragma once

#include <Eigen/Dense>
#include "edge.h"
#include "basic_utils.h"

namespace Delta2 {
    namespace common {
        enum class TriangleZone {
            V0, V1, V2, E01, E12, E20, region0
        };
        template<typename real>
        class Triangle {
            public:
            Eigen::Vector<real, 3> A, B, C;
            Triangle(const Eigen::Vector<real, 3>& Ap, const Eigen::Vector<real, 3>& Bp, const Eigen::Vector<real, 3>& Cp) {
                A = Ap;
                B = Bp;
                C = Cp;
            }
            real area() const {
                return 0.5 * ((A - C).cross(B - C)).norm();
            };
            Triangle(const Triangle<real>& T) : A(T.A), B(T.B), C(T.C) {};
            Triangle() {};
            real volume(real eps) const {
                // Volume of eps radius sphere swept over this triangle
                auto AB = A - B;
                auto AC = A - C;
                auto BC = B - C;
                real ABnorm = AB.norm();
                real ACnorm = AC.norm();
                real BCnorm = BC.norm();
                real sphere = 4.0 / 3.0 * 3.14159 * eps * eps * eps;
                real edges = 1.0 / 2.0 * 3.14159 * (ABnorm + BCnorm + ACnorm) * eps * eps;
                real midx = AC.y() * BC.z() - AC.z() * BC.y();
                real midy = AC.z() * BC.x() - AC.x() * BC.z();
                real midz = AC.x() * BC.y() - AC.y() * BC.x();
                real mid = sqrt(midx * midx + midy * midy + midz * midz) * 2 * eps;

                return sphere + edges + mid;
            };
            Eigen::Vector<real, 3> normal() const {
                return ((B - A).cross(C - A)).normalized();
            }
            real diameter() const {
                return std::max((A - B).norm(), std::max((B - C).norm(), (C - A).norm()));
            }
            Eigen::Vector<real, 3> operator[](int i) const {
                switch (i) {
                    case 0:
                        return A;
                    case 1:
                        return B;
                    case 2:
                        return C;
                    default:
                        throw std::invalid_argument("Index to triangle must be in range [0,2].");
                }
            }
            Edge<real> AB() const {
                return Edge<real>(A, B);
            }
            Edge<real> BC() const {
                return Edge<real>(B, C);
            }
            Edge<real> CA() const {
                return Edge<real>(C, A);
            }
            template<typename ToT>
            Triangle<ToT> cast() const {
                return Triangle<ToT>(A.template cast<ToT>(), B.template cast<ToT>(), C.template cast<ToT>());
            }
            bool operator==(const Triangle<real>& rhs)
            {
                return A == rhs.A && B == rhs.B && C == rhs.C;
            }
            bbox<real> bbox() const {
                real lower_x = std::min(std::min(A.x(), B.x()), C.x());
                real lower_y = std::min(std::min(A.y(), B.y()), C.y());
                real lower_z = std::min(std::min(A.z(), B.z()), C.z());
                real upper_x = std::max(std::max(A.x(), B.x()), C.x());
                real upper_y = std::max(std::max(A.y(), B.y()), C.y());
                real upper_z = std::max(std::max(A.z(), B.z()), C.z());
                Eigen::Vector<real, 3> lower = {lower_x, lower_y, lower_z};
                Eigen::Vector<real, 3> upper = {upper_x, upper_y, upper_z};
                return {lower, upper};
            }

            std::tuple<real, real, real> calcBarycentric(const Eigen::Vector<real, 3>& P) const {
                Eigen::Vector<real, 3> v0 = B - A, v1 = C - A, v2 = P - A;
                float d00 = v0.dot(v0);
                float d01 = v0.dot(v1);
                float d11 = v1.dot(v1);
                float d20 = v2.dot(v0);
                float d21 = v2.dot(v1);
                float denom = d00 * d11 - d01 * d01;
                float v = (d11 * d20 - d01 * d21) / denom;
                float w = (d00 * d21 - d01 * d20) / denom;
                float u = 1.0f - v - w;

                return {u, v, w};
            }

            Eigen::Vector<real, 3> intersectSegment(const Edge<real>& E, real& t, bool& t_inf, bool& isIntersecting) const {
                const Eigen::Vector<real, 3>& start = E.A;
                const Eigen::Vector<real, 3>& end = E.B;
                Eigen::Vector<real, 3> AP = A - start;
                Eigen::Vector<real, 3> n = normal();
                Eigen::Vector<real, 3> D = end - start;
                real p1 = AP.dot(n);
                real p2 = D.dot(n);
                t_inf = false;
                if (std::fabs(p2) <= std::numeric_limits<real>::epsilon()) {
                    t = 0.0;
                    t_inf = true;
                    isIntersecting = false;
                    return {0, 0, 0};
                }
                
                t = p1 / p2;
                if (0.0 < t && t < 1.0) {
                    auto [u, v, w] = calcBarycentric(start + D * t);

                    isIntersecting = 0.0 < u && u < 1.0 && 0.0 < v && v < 1.0 && u + v < 1.0; // Intersects the triangle
                    return start * (1-t) + end * t;
                } else {
                    isIntersecting = false;
                    return {0, 0, 0};
                }
            }

            Eigen::Vector<real, 3> projectToPlane(const Eigen::Vector<real, 3>& pt, bool &isInTriangle) const {
                Eigen::Vector<real, 3> AP = pt - A;
                Eigen::Vector<real, 3> n = normal();
                
                Eigen::Vector<real, 3> onPlane = pt - AP.dot(n) * n;
                auto [u, v, w] = calcBarycentric(onPlane);
                isInTriangle = 0 < u && u < 1 && 0 < v && v < 1 && u + v < 1;
                return onPlane;
            }

            TriangleZone projectPointToZone(const Eigen::Vector<real, 3>& P) const {
                TriangleZone z;

                Eigen::Vector<real, 3> AX = P - A;
                Eigen::Vector<real, 3> AB = B - A;
                Eigen::Vector<real, 3> AC = C- A;
                real a00 = AB.dot(AB);
                real a01 = AB.dot(AC);
                real a11 = AC.dot(AC);
                real b0 = -AX.dot(AB);
                real b1 = -AX.dot(AC);
                real det = a00 * a11 - a01 * a01;
                real t0 = a01 * b1 - a11 * b0;
                real t1 = a01 * b0 - a00 * b1;

                real t0_ret = t0 / det;
                real t1_ret = t1 / det;

                if (t0 + t1 <= det) {
                    if (t0 < 0.0) {
                        if (t1 < 0.0) {  // region 4
                            if (b0 < 0.0) {
                                if (-b0 >= a00) {  // V1
                                    z = TriangleZone::V1;
                                }
                                else {  // E01
                                    z = TriangleZone::E01;
                                }
                            }
                            else {
                                if (b1 >= 0.0) {  // V0
                                    z = TriangleZone::V0;
                                }
                                else if (-b1 >= a11) {  // V2
                                    z = TriangleZone::V2;
                                }
                                else {  // E20
                                    z = TriangleZone::E20;
                                }
                            }
                        }
                        else {  // region 3
                            if (b1 >= 0.0) {  // V0
                                z = TriangleZone::V0;
                            }
                            else if (-b1 >= a11) {  // V2
                                z = TriangleZone::V2;
                            }
                            else {  // E20
                                z = TriangleZone::E20;
                            }
                        }
                    }
                    else if (t1 < 0.0) {  // region 5
                        if (b0 >= 0.0) { // V0
                            z = TriangleZone::V0;
                        }
                        else if (-b0 >= a00) {  // V1
                            z = TriangleZone::V1;
                        }
                        else {  // E01
                            z = TriangleZone::E01;
                        }
                    }
                    else {  // region 0, interior
                        z = TriangleZone::region0;
                    }
                }
                else {
                    if (t0 < 0.0) {  // region 2
                        real tmp0 = a01 + b0;
                        real tmp1 = a11 + b1;
                        if (tmp1 > tmp0) {
                            real numer = tmp1 - tmp0;
                            real denom = a00 - 2 * a01 + a11;
                            if (numer >= denom) {  // V1
                                z = TriangleZone::V1;
                            }
                            else {  // E12
                                z = TriangleZone::E12;
                            }
                        }
                        else {
                            if (tmp1 <= 0.0) {  // V2
                                z = TriangleZone::V2;
                            }
                            else if (b1 >= 0.0) {  // V0
                                z = TriangleZone::V0;
                            }
                            else {  // E20
                                z = TriangleZone::E20;
                            }
                        }
                    }
                    else if (t1 < 0.0) {  // region 6
                        real tmp0 = a01 + b1;
                        real tmp1 = a00 + b0;
                        if (tmp1 > tmp0) {
                            real numer = tmp1 - tmp0;
                            real denom = a00 - 2 * a01 + a11;
                            if (numer >= denom) {  // V2
                                z = TriangleZone::V2;
                            }
                            else {  // E12
                                z = TriangleZone::E12;
                            }
                        }
                        else {
                            if (tmp1 <= 0.0) {  // V1
                                z = TriangleZone::V1;
                            }
                            else if (b0 >= 0.0) {  // V0
                                z = TriangleZone::V0;
                            }
                            else {  // E01
                                z = TriangleZone::E01;
                            }
                        }
                    }
                    else {  // region 1
                        real numer = a11 + b1 - a01 - b0;
                        if (numer <= 0.0) {  // V2
                            z = TriangleZone::V2;
                        }
                        else {
                            real denom = a00 - 2 * a01 + a11;
                            if (numer >= denom) {  // V1
                                z = TriangleZone::V1;
                            }
                            else {  // E12
                                z = TriangleZone::E12;
                            }
                        }
                    }
                }
                return z;
            }
            real distToPoint(const Eigen::Vector<real, 3> X) const {
                Eigen::Vector<real, 3> AX = X - A;
                Eigen::Vector<real, 3> AB = B - A;
                Eigen::Vector<real, 3> AC = C - A;
                real a00 = AB.dot(AB);
                real a01 = AB.dot(AC);
                real a11 = AC.dot(AC);
                real b0 = -AX.dot(AB);
                real b1 = -AX.dot(AC);
                real det = a00 * a11 - a01 * a01;
                real t0 = a01 * b1 - a11 * b0;
                real t1 = a01 * b0 - a00 * b1;

                if (t0 + t1 <= det)
                {
                    if (t0 < 0.0)
                    {
                        if (t1 < 0.0)  // region 4
                        {
                            if (b0 < 0.0)
                            {
                                t1 = 0.0;
                                if (-b0 >= a00)  // V1
                                {
                                    t0 = 1.0;
                                }
                                else  // E01
                                {
                                    t0 = -b0 / a00;
                                }
                            }
                            else
                            {
                                t0 = 0.0;
                                if (b1 >= 0.0)  // V0
                                {
                                    t1 = 0.0;
                                }
                                else if (-b1 >= a11)  // V2
                                {
                                    t1 = 1.0;
                                }
                                else  // E20
                                {
                                    t1 = -b1 / a11;
                                }
                            }
                        }
                        else  // region 3
                        {
                            t0 = 0.0;
                            if (b1 >= 0.0)  // V0
                            {
                                t1 = 0.0;
                            }
                            else if (-b1 >= a11)  // V2
                            {
                                t1 = 1.0;
                            }
                            else  // E20
                            {
                                t1 = -b1 / a11;
                            }
                        }
                    }
                    else if (t1 < 0.0)  // region 5
                    {
                        t1 = 0.0;
                        if (b0 >= 0.0)  // V0
                        {
                            t0 = 0.0;
                        }
                        else if (-b0 >= a00)  // V1
                        {
                            t0 = 1.0;
                        }
                        else  // E01
                        {
                            t0 = -b0 / a00;
                        }
                    }
                    else  // region 0, interior
                    {
                        real invDet = 1.0 / det;
                        t0 *= invDet;
                        t1 *= invDet;
                    }
                }
                else
                {
                    real tmp0, tmp1, numer, denom;

                    if (t0 < 0.0)  // region 2
                    {
                        tmp0 = a01 + b0;
                        tmp1 = a11 + b1;
                        if (tmp1 > tmp0)
                        {
                            numer = tmp1 - tmp0;
                            denom = a00 - 2.0*a01 + a11;
                            if (numer >= denom)  // V1
                            {
                                t0 = 1.0;
                                t1 = 0.0;
                            }
                            else  // E12
                            {
                                t0 = numer / denom;
                                t1 = 1.0 - t0;
                            }
                        }
                        else
                        {
                            t0 = 0.0;
                            if (tmp1 <= 0.0)  // V2
                            {
                                t1 = 1.0;
                            }
                            else if (b1 >= 0.0)  // V0
                            {
                                t1 = 0.0;
                            }
                            else  // E20
                            {
                                t1 = -b1 / a11;
                            }
                        }
                    }
                    else if (t1 < 0.0)  // region 6
                    {
                        tmp0 = a01 + b1;
                        tmp1 = a00 + b0;
                        if (tmp1 > tmp0)
                        {
                            numer = tmp1 - tmp0;
                            denom = a00 - 2.0*a01 + a11;
                            if (numer >= denom)  // V2
                            {
                                t1 = 1.0;
                                t0 = 0.0;
                            }
                            else  // E12
                            {
                                t1 = numer / denom;
                                t0 = 1.0 - t1;
                            }
                        }
                        else
                        {
                            t1 = 0.0;
                            if (tmp1 <= 0.0)  // V1
                            {
                                t0 = 1.0;
                            }
                            else if (b0 >= 0.0)  // V0
                            {
                                t0 = 0.0;
                            }
                            else  // E01
                            {
                                t0 = -b0 / a00;
                            }
                        }
                    }
                    else  // region 1
                    {
                        numer = a11 + b1 - a01 - b0;
                        if (numer <= 0.0)  // V2
                        {
                            t0 = 0.0;
                            t1 = 1.0;
                        }
                        else
                        {
                            denom = a00 - 2.0*a01 + a11;
                            if (numer >= denom)  // V1
                            {
                                t0 = 1.0;
                                t1 = 0.0;
                            }
                            else  // 12
                            {
                                t0 = numer / denom;
                                t1 = 1.0 - t0;
                            }
                        }
                    }
                }
                Eigen::Vector<real, 3> diff = X - (A + t0 * AB + t1 * AC);
                real sqrDistance = diff.dot(diff);
                return diff.norm();
            }
            Triangle<real> transformed(const Eigen::Matrix<real, 4, 4>& trans) const {
                Eigen::Vector<real, 4> a, b, c;

                a = A.homogeneous();
                b = B.homogeneous();
                c = C.homogeneous();

                Eigen::Vector<real, 4> ra, rb, rc;
                ra = trans * a;
                rb = trans * b;
                rc = trans * c;

                Triangle<real> result(ra.template head<3>() / ra.w(), rb.template head<3>() / rb.w(), rc.template head<3>() / rc.w());

                return result;
            }
            void toEigen(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
                V = Eigen::MatrixXd({
                    {A.x(), A.y(), A.z()},
                    {B.x(), B.y(), B.z()},
                    {C.x(), C.y(), C.z()}
                });
                F = Eigen::MatrixXi({
                    {0, 1, 2}
                });
            }
            /* Return the normal and origin of a plane that this triangle lies on */
            std::pair<Eigen::Vector<real, 3>, Eigen::Vector<real, 3>> toNormalRep() const {
                return std::make_pair((B-A).cross(C-A).normalized(), A);
            }
        };
        template class Triangle<float>;
        template class Triangle<double>;
    }
}
