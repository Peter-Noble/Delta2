#pragma once

#include <Eigen/Dense>
#include "../common/triangle.h"
#include "../model/collision_state.h"
#include "contact.h"
#include "../common/viewer.h"
#include "../globals.h"

#include "tbb/tbb.h"

#include <ittnotify.h>

namespace Delta2 {
    namespace collision {
        template<typename real>
        real sphereSphereFirstIntersect(const Eigen::Vector<real, 3>& a_pos, const Eigen::Vector<real, 3>& a_vel, real a_rad, const Eigen::Vector<real, 3>& b_pos, const Eigen::Vector<real, 3>& b_vel, real b_rad)
        {
            Eigen::Vector<real, 3> ab = a_pos - b_pos; // vector between the centers of each sphere
            Eigen::Vector<real, 3> rel_vel = a_vel - b_vel; // relative velocity between spheres
            real r = a_rad + b_rad;
        
            real c = ab.dot(ab) - r*r; // if negative, they overlap
            if (c < 0.0) // if true, they already overlap
            {
                return 0.0;
            }
        
            real a = rel_vel.dot(rel_vel);
        
            real b = rel_vel.dot(ab);
            if (b >= 0.0)
                return -1.0; // does not move towards each other
        
            real d = b*b - a*c;
            if (d < 0.0)
                return -1.0; // no real roots ... no collision
        
            return (-b - sqrt(d)) / a;
        }
        template float sphereSphereFirstIntersect(const Eigen::Vector3f&, const Eigen::Vector3f&, float, const Eigen::Vector3f&, const Eigen::Vector3f&, float);
        template double sphereSphereFirstIntersect(const Eigen::Vector3d&, const Eigen::Vector3d&, double, const Eigen::Vector3d&, const Eigen::Vector3d&, double);

        // Only check point to triangle surface (ie exclude point to edge because that will be included in edge-edge)?
        template<typename real>
        real pointTriCCDLinear(const Eigen::Vector<real, 3> &pt, const Eigen::Matrix<real, 4, 4> &pt_T_start, const Eigen::Matrix<real, 4, 4> &pt_T_end, const Delta2::common::Triangle<real> &tri, const Eigen::Matrix<real, 4, 4> &tri_T_start, const Eigen::Matrix<real, 4, 4> &tri_T_end, Eigen::Vector<real, 3> &P_out, Eigen::Vector<real, 3> &Q_out, real &TOC_out) {
            Eigen::Vector<real, 4> pt_start_temp = tri_T_start.inverse() * pt_T_start * pt.homogeneous(); // in local triangle space
            Eigen::Vector<real, 3> pt_start = pt_start_temp.template head<3>() / pt_start_temp.w();
            Eigen::Vector<real, 4> pt_end_temp = tri_T_end.inverse() * pt_T_end * pt.homogeneous();
            Eigen::Vector<real, 3> pt_end = pt_end_temp.template head<3>() / pt_end_temp.w();

            common::Edge e(pt_start, pt_end);

            TOC_out = 1.0;

            bool isIntersecting;
            real t;
            bool t_inf;
            Eigen::Vector<real, 3> intersect = tri.intersectSegment(e, t, t_inf, isIntersecting);

            if (0 < t && t < 1 && isIntersecting) {
                Eigen::Matrix<real, 4, 4> t_T = (1 - t) * tri_T_start + t * tri_T_end;
                P_out = Delta2::common::transform(intersect, t_T);
                Q_out = Delta2::common::transform(intersect, t_T);
                TOC_out = t;
                return 0;
            } else if (t_inf) {
                t = 0;
                Eigen::Vector<real, 3> closest = e.lerp(t);
                bool isClosestAboveTriangle;
                Eigen::Vector<real, 3> onPlane = tri.projectToPlane(closest, isClosestAboveTriangle);
                if (isClosestAboveTriangle) {
                    P_out = Delta2::common::transform(closest, tri_T_start);
                    Q_out = Delta2::common::transform(onPlane, tri_T_start);
                    TOC_out = t;
                    return (P_out - Q_out).norm();
                }

                t = 1;
                closest = e.lerp(t);
                onPlane = tri.projectToPlane(closest, isClosestAboveTriangle);
                if (isClosestAboveTriangle) {
                    P_out = Delta2::common::transform(closest, tri_T_end);
                    Q_out = Delta2::common::transform(onPlane, tri_T_end);
                    TOC_out = t;
                    return (P_out - Q_out).norm();
                }
            } else  {
                t = Delta2::common::clamp01(t);
                Eigen::Matrix<real, 4, 4> t_T = (1 - t) * tri_T_start + t * tri_T_end;
                Eigen::Vector<real, 3> closest = e.lerp(t);
                bool isClosestAboveTriangle;
                Eigen::Vector<real, 3> onPlane = tri.projectToPlane(closest, isClosestAboveTriangle);
                if (isClosestAboveTriangle) {
                    P_out = Delta2::common::transform(closest, t_T);
                    Q_out = Delta2::common::transform(onPlane, t_T);
                    TOC_out = t;
                    return (P_out - Q_out).norm();
                }
            }
            return std::numeric_limits<real>::infinity();
        }
        template float pointTriCCDLinear(const Eigen::Vector3f&, const Eigen::Matrix4f&, const Eigen::Matrix4f&, const Delta2::common::Triangle<float>&, const Eigen::Matrix4f&, const Eigen::Matrix4f&, Eigen::Vector3f&, Eigen::Vector3f&, float&);
        template double pointTriCCDLinear(const Eigen::Vector3d&, const Eigen::Matrix4d&, const Eigen::Matrix4d&, const Delta2::common::Triangle<double>&, const Eigen::Matrix4d&, const Eigen::Matrix4d&, Eigen::Vector3d&, Eigen::Vector3d&, double&);

        /**
         * Return the shortest distance between two line segments as their points are linearly interpolated through the time interval 0-1
         * (a_first_start, a_second_start) - line a at the start of the time interval
         * P_out is the position on line a at the closest
         * Q_out is the position on line b at the closest
         * TOC_out is the time (0-1) the closest pass occurs
         */
        template<typename real>
        real edgeEdgeCCDLinear(const Delta2::common::Edge<real>& a_start, const Delta2::common::Edge<real>& a_end, const Delta2::common::Edge<real>& b_start, const Delta2::common::Edge<real>& b_end, Eigen::Vector<real, 3> &P_out, Eigen::Vector<real, 3> &Q_out, real &TOC) {
            // The start and end position of a line are taken.  A crude approximation is made to the volume that it traverses (tetrahedron) and then tested against.

            // TODO by doing all this at the tri-tri level (or even better cluster-cluster level) we could eliminate a whole bunch of duplicate checks.

            //https://zalo.github.io/blog/closest-point-between-segments/
            Eigen::Vector<real, 3> origin_start = a_start.A;
            Eigen::Vector<real, 3> origin_end = a_end.A;

            Eigen::Vector<real, 3> n_start = (a_start.A - a_start.B).normalized();
            Eigen::Vector<real, 3> n_end = (a_end.A - a_end.B).normalized();

            Eigen::Quaternion<real> rot = Eigen::Quaternion<real>().setFromTwoVectors(n_end, n_start);
            Eigen::Matrix<real, 3, 3> R_end_to_start = rot.toRotationMatrix();
            // optimisation could be done here as this matrix will be the same for all edge pairs between a given pair of triangles.  Is that true?

            Eigen::Vector<real, 3> b_first_start_local = b_start.A;
            Eigen::Vector<real, 3> b_second_start_local = b_start.B;
            Eigen::Vector<real, 3> b_first_end_local = R_end_to_start * (b_end.A - origin_end)  + origin_start;
            Eigen::Vector<real, 3> b_second_end_local = R_end_to_start * (b_end.B - origin_end) + origin_start;

            Eigen::Vector<real, 3> P_closest, Q_closest;
            real dist_min = std::numeric_limits<real>::infinity();
            real t_min = std::numeric_limits<real>::infinity();

            common::Edge<real> E;

            Eigen::Vector<real, 3> P_t0, Q_t0;
            E = common::Edge(b_first_start_local, b_second_start_local);
            real dist_t0 = a_start.closest(E, P_t0, Q_t0);
            if (dist_t0 < dist_min || (dist_t0 == dist_min && 0 < t_min)) {
                t_min = 0;
                dist_min = dist_t0;
                P_closest = Eigen::Quaternion<real>::Identity().slerp(0, rot).toRotationMatrix().inverse() * P_t0;
                Q_closest = Eigen::Quaternion<real>::Identity().slerp(0, rot).toRotationMatrix().inverse() * Q_t0;
            }

            Eigen::Vector<real, 3> P_t1, Q_t1;
            E = common::Edge(b_first_end_local, b_second_end_local);
            real dist_t1 = a_start.closest(E, P_t1, Q_t1);
            if (dist_t1 < dist_min || (dist_t1 == dist_min && 1 < t_min)) {
                t_min = 1;
                dist_min = dist_t1;
                P_closest = P_t1;
                Q_closest = Q_t1;
            }

            Eigen::Vector<real, 3> P_first, Q_first;
            E = common::Edge(b_first_start_local, b_first_end_local);
            real dist_first = a_start.closest(E, P_first, Q_first);
            Eigen::Vector<real, 3> D = b_first_end_local - b_first_start_local;
            real t_first = D.dot(Q_first - b_first_start_local) / D.dot(D);
            if (dist_first < dist_min || (dist_first == dist_min && t_first < t_min)) {
                t_min = t_first;
                dist_min = dist_first;
                P_closest = P_first;
                Q_closest = Q_first;
            }

            Eigen::Vector<real, 3> P_second, Q_second;
            E = common::Edge(b_second_start_local, b_second_end_local);
            real dist_second = a_start.closest(E, P_second, Q_second);
            D = b_second_end_local - b_second_start_local;
            real t_second = D.dot(Q_second - b_second_start_local) / D.dot(D);
            if (dist_second < dist_min || (dist_second == dist_min && t_second < t_min)) {
                t_min = t_second;
                dist_min = dist_second;
                P_closest = P_second;
                Q_closest = Q_second;
            }

            // If there is very little relative rotation between the two triangles then the cross edges will effectively be on the plane with the triangles and these
            //	extra checks serve no real purpose.
            Eigen::Vector<real, 3> P_t0cross, Q_t0cross;
            E = common::Edge(b_first_start_local, b_second_end_local);
            real dist_t0cross = a_start.closest(E, P_t0cross, Q_t0cross);
            D = b_second_end_local - b_first_start_local;
            real t_t0cross = D.dot(Q_t0cross - b_first_start_local) / D.dot(D);
            if (dist_t0cross < dist_min || (dist_t0cross == dist_min && t_t0cross < t_min)) {
                t_min = t_t0cross;
                dist_min = dist_t0cross;
                P_closest = P_t0cross;
                Q_closest = Q_t0cross;
            }

            Eigen::Vector<real, 3> P_t1cross, Q_t1cross;
            E = common::Edge(b_first_start_local, b_second_start_local);
            real dist_t1cross = a_start.closest(E, P_t1cross, Q_t1cross);
            D = b_second_start_local - b_first_start_local;
            real t_t1cross = D.dot(Q_t1cross - b_first_start_local) / D.dot(D);
            if (dist_t1cross < dist_min || (dist_t1cross == dist_min && t_t1cross < t_min)) {
                t_min = t_t1cross;
                dist_min = dist_t1cross;
                P_closest = P_t1cross;
                Q_closest = Q_t1cross;
            }

            // Check points against faces
            Delta2::common::Triangle<real> tri_lower_a = Delta2::common::Triangle<real>(b_first_start_local, b_first_end_local, b_second_start_local);
            bool is_intersecting_lower_a;
            real t_lower_a;
            bool t_inf;
            Eigen::Vector<real, 3> intersect_lower_a = tri_lower_a.intersectSegment(a_start, t_lower_a, t_inf, is_intersecting_lower_a);
            if (is_intersecting_lower_a) {
                auto [u, v, w] = tri_lower_a.calcBarycentric(intersect_lower_a);
                real t = 1 - u;
                if (t < t_min) {
                    dist_min = 0.0;
                    t_min = t;
                    P_closest = intersect_lower_a;
                    Q_closest = intersect_lower_a;
                }
            }

            bool pt_to_face;
            Eigen::Vector<real, 3> proj_to_face;
            proj_to_face = tri_lower_a.projectToPlane(a_start.A, pt_to_face);
            if (pt_to_face && (proj_to_face - a_start.A).norm() < dist_min) {
                auto [u, v, _] = tri_lower_a.calcBarycentric(proj_to_face);
                real t = 1 - u;
                if (t < t_min) {
                    dist_min = (proj_to_face - a_start.A).norm();
                    t_min = t;
                    P_closest = a_start.A;
                    Q_closest = proj_to_face;
                }
            }
            proj_to_face = tri_lower_a.projectToPlane(a_start.B, pt_to_face);
            if (pt_to_face && (proj_to_face - a_start.B).norm() < dist_min) {
                auto [u, v, _] = tri_lower_a.calcBarycentric(proj_to_face);
                real t = 1 - u;
                if (t < t_min) {
                    dist_min = (proj_to_face - a_start.B).norm();
                    t_min = t;
                    P_closest = a_start.B;
                    Q_closest = proj_to_face;
                }
            }

            Delta2::common::Triangle<real> tri_lower_b = Delta2::common::Triangle<real>(b_first_start_local, b_second_end_local, b_second_start_local);
            bool is_intersecting_lower_b;
            real t_lower_b;
            Eigen::Vector<real, 3> intersect_lower_b = tri_lower_b.intersectSegment(a_start, t_lower_b, t_inf, is_intersecting_lower_b);
            if (is_intersecting_lower_b) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(intersect_lower_b);
                if (v < t_min) {
                    dist_min = 0.0;
                    t_min = v;
                    P_closest = intersect_lower_b;
                    Q_closest = P_closest;
                }
            }

            proj_to_face = tri_lower_b.projectToPlane(a_start.A, pt_to_face);
            if (pt_to_face) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(proj_to_face);
                if (v < t_min && (proj_to_face - a_start.A).norm() < dist_min) {
                    dist_min = (proj_to_face - a_start.A).norm();
                    t_min = v;
                    P_closest = a_start.A;
                    Q_closest = proj_to_face;
                }
            }
            proj_to_face = tri_lower_b.projectToPlane(a_start.B, pt_to_face);
            if (pt_to_face && (proj_to_face - a_start.B).norm() < dist_min) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(proj_to_face);
                if (v < t_min) {
                    dist_min = (proj_to_face - a_start.B).norm();
                    t_min = v;
                    P_closest = a_start.B;
                    Q_closest = proj_to_face;
                }
            }
            
            // If there is very little relative rotation then these second two will produce almost exactly the same result as the first => skip?
            /*Delta2::common::Triangle<real> tri_upper_a = Delta2::common::Triangle<real>(b_start.A, b_end.A, b_start_local.B);
            bool is_intersecting_upper_a;
            real t_upper_a;
            Eigen::Vector<real, 3> intersect_upper_a = tri_upper_a.intersectSegment(a_start, t_upper_a, t_inf, is_intersecting_upper_a);
            if (is_intersecting_upper_a) {
                auto [u, v, _] = tri_upper_a.calcBarycentric(intersect_upper_a);
                if (v < t_min) {
                    dist_min = 0.0;
                    t_min = u;
                    P_closest = intersect_upper_a;
                    Q_closest = intersect_upper_a;
                }
            }

            Delta2::common::Triangle<real> tri_upper_b = Delta2::common::Triangle<real>(b_end_local.B, b_end.A, b_second_start_local);
            bool is_intersecting_upper_b;
            real t_upper_b;
            Eigen::Vector<real, 3> intersect_upper_b = tri_upper_b.intersectSegment(a_start, t_upper_b, t_inf, is_intersecting_upper_b);
            if (is_intersecting_upper_b) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(intersect_upper_b);
                if (v < t_min) {
                    dist_min = 0.0;
                    t_min = v;
                    P_closest = intersect_upper_b;
                    Q_closest = intersect_upper_b;
                }
            }
            TODO pt-> faces hasn't been done for these two.  See above.
            */

            Eigen::Vector<real, 3> origin_now = (1-t_min) * origin_start + t_min * origin_end;
            P_out = Eigen::Quaternion<real>::Identity().slerp(t_min, rot).toRotationMatrix().inverse() * (P_closest - origin_start) + origin_now;
            Q_out = Eigen::Quaternion<real>::Identity().slerp(t_min, rot).toRotationMatrix().inverse() * (Q_closest - origin_start) + origin_now;
            TOC = t_min;
            return dist_min;
        }
        template float edgeEdgeCCDLinear(const Delta2::common::Edge<float>&, const Delta2::common::Edge<float>&, const Delta2::common::Edge<float>&, const Delta2::common::Edge<float>&, Eigen::Vector3f&, Eigen::Vector3f&, float&);
        template double edgeEdgeCCDLinear(const Delta2::common::Edge<double>&, const Delta2::common::Edge<double>&, const Delta2::common::Edge<double>&, const Delta2::common::Edge<double>&, Eigen::Vector3d&, Eigen::Vector3d&, double&);

        /**
         * Return the shortest distance between two line segments as their points are linearly interpolated through the time interval 0-1
         * (a_first_start, a_second_start) - line a at the start of the time interval
         * P_out is the position on line a at the closest
         * Q_out is the position on line b at the closest
         * TOC_out is the time (0-1) the closest pass occurs
         */
        template<typename real>
        real edgeEdgeCCDLinear(const Delta2::common::Edge<real>& a_orig, const Eigen::Matrix<real, 4, 4>& a_start, const Eigen::Matrix<real, 4, 4>& a_end, const Delta2::common::Edge<real>& b, const Eigen::Matrix<real, 4, 4>& b_start, const Eigen::Matrix<real, 4, 4>& b_end, Eigen::Vector<real, 3> &P_out, Eigen::Vector<real, 3> &Q_out, real &TOC) {
            // The start and end position of a line are taken.  A crude approximation is made to the volume that it traverses (tetrahedron) and then tested against.

            // TODO by doing all this at the tri-tri level (or even better cluster-cluster level) we could eliminate a whole bunch of duplicate checks.

            //https://zalo.github.io/blog/closest-point-between-segments/

            Eigen::Vector<real, 3> a_first_start = common::transform(a_orig.A, a_start);
            Eigen::Vector<real, 3> a_second_start = common::transform(a_orig.B, a_start);
            Eigen::Vector<real, 3> a_first_end = common::transform(a_orig.A, a_end);
            Eigen::Vector<real, 3> a_second_end = common::transform(a_orig.B, a_end);
            common::Edge<real> a(a_first_start, a_second_start);

            Eigen::Matrix<real, 4, 4> b_end_diff_T = a_start * a_end.inverse() * b_end;
            // Eigen::Matrix<real, 4, 4> b_end_diff_T = b_end * a_end.inverse() * a_start;
            Eigen::Vector<real, 3> b_first_start_local = common::transform(b.A, b_start);
            Eigen::Vector<real, 3> b_second_start_local = common::transform(b.B, b_start);
            Eigen::Vector<real, 3> b_first_end_local = common::transform(b.A, b_end_diff_T);
            Eigen::Vector<real, 3> b_second_end_local = common::transform(b.B, b_end_diff_T);

            Eigen::Vector<real, 3> P_closest, Q_closest;
            real dist_min = std::numeric_limits<real>::infinity();
            real t_min = std::numeric_limits<real>::infinity();

            common::Edge<real> E;

            Eigen::Vector<real, 3> P_t0, Q_t0;
            E = common::Edge(b_first_start_local, b_second_start_local);
            real dist_t0 = a.closest(E, P_t0, Q_t0);
            if (dist_t0 < dist_min || (dist_t0 == dist_min && 0 < t_min)) {
                t_min = 0;
                dist_min = dist_t0;
                P_closest = P_t0;
                Q_closest = Q_t0;
            }

            Eigen::Vector<real, 3> P_t1, Q_t1;
            E = common::Edge(b_first_end_local, b_second_end_local);
            real dist_t1 = a.closest(E, P_t1, Q_t1);
            if (dist_t1 < dist_min || (dist_t1 == dist_min && 1 < t_min)) {
                t_min = 1;
                dist_min = dist_t1;
                P_closest = P_t1;
                Q_closest = Q_t1;
            }

            Eigen::Vector<real, 3> P_first, Q_first;
            E = common::Edge(b_first_start_local, b_first_end_local);
            real dist_first = a.closest(E, P_first, Q_first);
            Eigen::Vector<real, 3> D = b_first_end_local - b_first_start_local;
            real t_first = D.dot(Q_first - b_first_start_local) / D.dot(D);
            if (dist_first < dist_min || (dist_first == dist_min && t_first < t_min)) {
                t_min = t_first;
                dist_min = dist_first;
                P_closest = P_first;
                Q_closest = Q_first;
            }

            Eigen::Vector<real, 3> P_second, Q_second;
            E = common::Edge(b_second_start_local, b_second_end_local);
            real dist_second = a.closest(E, P_second, Q_second);
            D = b_second_end_local - b_second_start_local;
            real t_second = D.dot(Q_second - b_second_start_local) / D.dot(D);
            if (dist_second < dist_min || (dist_second == dist_min && t_second < t_min)) {
                t_min = t_second;
                dist_min = dist_second;
                P_closest = P_second;
                Q_closest = Q_second;
            }

            // If there is very little relative rotation between the two triangles then the cross edges will effectively be on the plane with the triangles and these
            //	extra checks serve no real purpose.
            Eigen::Vector<real, 3> P_t0cross, Q_t0cross;
            E = common::Edge(b_first_start_local, b_second_end_local);
            real dist_t0cross = a.closest(E, P_t0cross, Q_t0cross);
            D = b_second_end_local - b_first_start_local;
            real t_t0cross = D.dot(Q_t0cross - b_first_start_local) / D.dot(D);
            if (dist_t0cross < dist_min || (dist_t0cross == dist_min && t_t0cross < t_min)) {
                t_min = t_t0cross;
                dist_min = dist_t0cross;
                P_closest = P_t0cross;
                Q_closest = Q_t0cross;
            }

            Eigen::Vector<real, 3> P_t1cross, Q_t1cross;
            E = common::Edge(b_first_start_local, b_second_start_local);
            real dist_t1cross = a.closest(E, P_t1cross, Q_t1cross);
            D = b_second_start_local - b_first_start_local;
            real t_t1cross = D.dot(Q_t1cross - b_first_start_local) / D.dot(D);
            if (dist_t1cross < dist_min || (dist_t1cross == dist_min && t_t1cross < t_min)) {
                t_min = t_t1cross;
                dist_min = dist_t1cross;
                P_closest = P_t1cross;
                Q_closest = Q_t1cross;
            }

            // Check points against faces
            Delta2::common::Triangle<real> tri_lower_a = Delta2::common::Triangle<real>(b_first_start_local, b_first_end_local, b_second_start_local);
            bool is_intersecting_lower_a;
            real t_lower_a;
            bool t_inf;
            Eigen::Vector<real, 3> intersect_lower_a = tri_lower_a.intersectSegment(a, t_lower_a, t_inf, is_intersecting_lower_a);
            if (is_intersecting_lower_a) {
                auto [u, v, w] = tri_lower_a.calcBarycentric(intersect_lower_a);
                real t = 1 - u;
                if (t < t_min) {
                    dist_min = 0.0;
                    t_min = t;
                    P_closest = intersect_lower_a;
                    Q_closest = intersect_lower_a;
                }
            }

            bool pt_to_face;
            Eigen::Vector<real, 3> proj_to_face;
            proj_to_face = tri_lower_a.projectToPlane(a.A, pt_to_face);
            if (pt_to_face && (proj_to_face - a.A).norm() < dist_min) {
                auto [u, v, _] = tri_lower_a.calcBarycentric(proj_to_face);
                real t = 1 - u;
                if (t < t_min) {
                    dist_min = (proj_to_face - a.A).norm();
                    t_min = t;
                    P_closest = a.A;
                    Q_closest = proj_to_face;
                }
            }
            proj_to_face = tri_lower_a.projectToPlane(a.B, pt_to_face);
            if (pt_to_face && (proj_to_face - a.B).norm() < dist_min) {
                auto [u, v, _] = tri_lower_a.calcBarycentric(proj_to_face);
                real t = 1 - u;
                if (t < t_min) {
                    dist_min = (proj_to_face - a.B).norm();
                    t_min = t;
                    P_closest = a.B;
                    Q_closest = proj_to_face;
                }
            }

            Delta2::common::Triangle<real> tri_lower_b = Delta2::common::Triangle<real>(b_first_start_local, b_second_end_local, b_second_start_local);
            bool is_intersecting_lower_b;
            real t_lower_b;
            Eigen::Vector<real, 3> intersect_lower_b = tri_lower_b.intersectSegment(a, t_lower_b, t_inf, is_intersecting_lower_b);
            if (is_intersecting_lower_b) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(intersect_lower_b);
                if (v < t_min) {
                    dist_min = 0.0;
                    t_min = v;
                    P_closest = intersect_lower_b;
                    Q_closest = P_closest;
                }
            }

            proj_to_face = tri_lower_b.projectToPlane(a.A, pt_to_face);
            if (pt_to_face) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(proj_to_face);
                if (v < t_min && (proj_to_face - a.A).norm() < dist_min) {
                    dist_min = (proj_to_face - a.A).norm();
                    t_min = v;
                    P_closest = a.A;
                    Q_closest = proj_to_face;
                }
            }
            proj_to_face = tri_lower_b.projectToPlane(a.B, pt_to_face);
            if (pt_to_face && (proj_to_face - a.B).norm() < dist_min) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(proj_to_face);
                if (v < t_min) {
                    dist_min = (proj_to_face - a.B).norm();
                    t_min = v;
                    P_closest = a.B;
                    Q_closest = proj_to_face;
                }
            }
            
            // If there is very little relative rotation then these second two will produce almost exactly the same result as the first => skip?
            /*Delta2::common::Triangle<real> tri_upper_a = Delta2::common::Triangle<real>(b_start.A, b_end.A, b_start_local.B);
            bool is_intersecting_upper_a;
            real t_upper_a;
            Eigen::Vector<real, 3> intersect_upper_a = tri_upper_a.intersectSegment(a, t_upper_a, t_inf, is_intersecting_upper_a);
            if (is_intersecting_upper_a) {
                auto [u, v, _] = tri_upper_a.calcBarycentric(intersect_upper_a);
                if (v < t_min) {
                    dist_min = 0.0;
                    t_min = u;
                    P_closest = intersect_upper_a;
                    Q_closest = intersect_upper_a;
                }
            }

            Delta2::common::Triangle<real> tri_upper_b = Delta2::common::Triangle<real>(b_end_local.B, b_end.A, b_second_start_local);
            bool is_intersecting_upper_b;
            real t_upper_b;
            Eigen::Vector<real, 3> intersect_upper_b = tri_upper_b.intersectSegment(a, t_upper_b, t_inf, is_intersecting_upper_b);
            if (is_intersecting_upper_b) {
                auto [u, v, _] = tri_lower_b.calcBarycentric(intersect_upper_b);
                if (v < t_min) {
                    dist_min = 0.0;
                    t_min = v;
                    P_closest = intersect_upper_b;
                    Q_closest = intersect_upper_b;
                }
            }
            TODO pt-> faces hasn't been done for these two.  See above.
            */

            Eigen::Vector<real, 3> P_start = P_closest;
            Eigen::Matrix<real, 4, 4> a_start_to_end_T = a_end * a_start.inverse();
            Eigen::Vector<real, 3> P_end = Delta2::common::transform(P_closest, a_start_to_end_T);

            Eigen::Vector<real, 3> Q_start = Q_closest;
            Eigen::Vector<real, 3> Q_end = Delta2::common::transform(P_closest, a_start_to_end_T);

            Eigen::Quaternion<real> a_rot_start(a_start.template block<3,3>(0,0));
            Eigen::Vector<real, 3> a_trans_start(a_start.template block<3,1>(0,3));
            Eigen::Quaternion<real> a_rot_end(a_end.template block<3,3>(0,0));
            Eigen::Vector<real, 3> a_trans_end(a_end.template block<3,1>(0,3));
            Eigen::Quaternion<real> b_rot_start(b_start.template block<3,3>(0,0));
            Eigen::Vector<real, 3> b_trans_start(b_start.template block<3,1>(0,3));
            Eigen::Quaternion<real> b_rot_end(b_end.template block<3,3>(0,0));
            Eigen::Vector<real, 3> b_trans_end(b_end.template block<3,1>(0,3));

            Eigen::Quaternion<real> a_rot_now = a_rot_start.slerp(1.0-t_min, a_rot_end);
            Eigen::Vector<real, 3> a_trans_now = common::lerp(a_trans_start, a_trans_end, t_min);
            // Eigen::Quaternion<real> b_rot_now = b_rot_start.slerp(1.0-t_min, b_rot_end);
            // Eigen::Vector<real, 3> b_trans_now = common::lerp(b_trans_start, b_trans_end, t_min);

            P_out = a_rot_now.toRotationMatrix() * a_rot_start.toRotationMatrix().inverse() * (P_closest - a_trans_start) + a_trans_now;
            Q_out = a_rot_now.toRotationMatrix() * a_rot_start.toRotationMatrix().inverse() * (P_closest - a_trans_start) + a_trans_now;
            TOC = t_min;
            return dist_min;
        }
        template float edgeEdgeCCDLinear(const Delta2::common::Edge<float>&, const Eigen::Matrix<float, 4, 4>&, const Eigen::Matrix<float, 4, 4>&, const Delta2::common::Edge<float>&, const Eigen::Matrix<float, 4, 4>&, const Eigen::Matrix<float, 4, 4>&, Eigen::Vector<float, 3>&, Eigen::Vector<float, 3>&, float&);
        template double edgeEdgeCCDLinear(const Delta2::common::Edge<double>&, const Eigen::Matrix<double, 4, 4>&, const Eigen::Matrix<double, 4, 4>&, const Delta2::common::Edge<double>&, const Eigen::Matrix<double, 4, 4>&, const Eigen::Matrix<double, 4, 4>&, Eigen::Vector<double, 3>&, Eigen::Vector<double, 3>&, double&);

        template<typename real>
        bool testAxis(real a_lower, real a_upper, real b_lower, real b_upper, real max_dist) {
            if (b_lower > a_upper + max_dist) {
                return false;
            }
            if (b_upper + max_dist < a_lower) {
                return false;
            }
            return true;
        }

        double edgeEdgeCCD(const common::Edge<double>& a_start, const common::Edge<double>& a_end, const common::Edge<double>& b_start, const common::Edge<double>& b_end, double search_dist, Eigen::Vector<double, 3>& P_out, Eigen::Vector<double, 3>& Q_out, double& toc_out);
        float edgeEdgeCCD(const common::Edge<float>& a_start, const common::Edge<float>& a_end, const common::Edge<float>& b_start, const common::Edge<float>& b_end, float search_dist, Eigen::Vector<float, 3>& P_out, Eigen::Vector<float, 3>& Q_out, float& toc_out);

        /**
         * Return the shortest distance between two triangles as their points are linearly interpolated through the time interval 0-1
         */
        template<typename real>
        real triTriCCDLinear(const Delta2::common::Triangle<real> &a, const Eigen::Matrix<real, 4, 4>& a_T_start, const Eigen::Matrix<real, 4, 4> &a_T_end, const Delta2::common::Triangle<real> &b, const Eigen::Matrix<real, 4, 4> &b_T_start, const Eigen::Matrix<real, 4, 4> &b_T_end, Eigen::Vector<real, 3> &P_out, Eigen::Vector<real, 3> &Q_out, real &TOC_out) {
            real dist_min = std::numeric_limits<real>::infinity();
            TOC_out = 1.0;

            Delta2::common::Triangle<real> a_start = a.transformed(a_T_start);
            Delta2::common::Triangle<real> a_end = a.transformed(a_T_end);
            Delta2::common::Triangle<real> b_start = b.transformed(b_T_start);
            Delta2::common::Triangle<real> b_end = b.transformed(b_T_end);

            Eigen::Vector<real, 3> a_start_A, a_start_B, a_end_A, a_end_B, b_start_A, b_start_B, b_end_A, b_end_B;

            for (int i = 0; i < 3; i++) {
                Eigen::Vector<real, 3> pt_a, pt_b;
                if (i == 0) {
                    pt_a = a.A;
                    pt_b = b.A;
                } else if (i == 1) {
                    pt_a = a.B;
                    pt_b = b.B;
                } else {
                    pt_a = a.C;
                    pt_b = b.C;
                }

                Eigen::Vector<real, 3> P_closest_a_to_b, Q_closest_a_to_b;
                real dist_toc_a_to_b;
                real dist_a_to_b = pointTriCCDLinear(pt_a, a_T_start, a_T_end, b, b_T_start, b_T_end, P_closest_a_to_b, Q_closest_a_to_b, dist_toc_a_to_b);

                if (dist_a_to_b < dist_min || (dist_a_to_b == dist_min && dist_toc_a_to_b < TOC_out)) {
                    P_out = P_closest_a_to_b;
                    Q_out = Q_closest_a_to_b;
                    dist_min = dist_a_to_b;
                    TOC_out = dist_toc_a_to_b;
                }

                Eigen::Vector<real, 3> P_closest_b_to_a, Q_closest_b_to_a;
                real dist_toc_b_to_a;
                real dist_b_to_a = pointTriCCDLinear(pt_b, b_T_start, b_T_end, a, a_T_start, a_T_end, P_closest_b_to_a, Q_closest_b_to_a, dist_toc_b_to_a);

                if (dist_b_to_a < dist_min || (dist_b_to_a == dist_min && dist_toc_b_to_a < TOC_out)) {
                    P_out = P_closest_b_to_a;
                    Q_out = Q_closest_b_to_a;
                    dist_min = dist_b_to_a;
                    TOC_out = dist_toc_b_to_a;
                }

                // Edge edge CCD tests
                if (i == 0) {
                    a_start_A = a_start.A;
                    a_start_B = a_start.B;
                    a_end_A = a_end.A;
                    a_end_B = a_end.B;
                } else if (i == 1) {
                    a_start_A = a_start.B;
                    a_start_B = a_start.C;
                    a_end_A = a_end.B;
                    a_end_B = a_end.C;
                } else {
                    a_start_A = a_start.C;
                    a_start_B = a_start.A;
                    a_end_A = a_end.C;
                    a_end_B = a_end.A;
                }
                for (int j = 0; j < 3; j++) {
                    if (j == 0) {
                        b_start_A = b_start.A;
                        b_start_B = b_start.B;
                        b_end_A = b_end.A;
                        b_end_B = b_end.B;
                    } else if (j == 1) {
                        b_start_A = b_start.B;
                        b_start_B = b_start.C;
                        b_end_A = b_end.B;
                        b_end_B = b_end.C;
                    } else {
                        b_start_A = b_start.C;
                        b_start_B = b_start.A;
                        b_end_A = b_end.C;
                        b_end_B = b_end.A;
                    }

                    Eigen::Vector<real, 3> P, Q;
                    real t;
                    common::Edge<real> a_edge_start(a_start_A, a_start_B);
                    common::Edge<real> b_edge_start(b_start_A, b_start_B);
                    common::Edge<real> a_edge_end(a_end_A, a_end_B);
                    common::Edge<real> b_edge_end(b_end_A, b_end_B);
                    real dist = edgeEdgeCCD(a_edge_start, a_edge_end, b_edge_start, b_edge_end, dist_min, P, Q, t);

                    if (dist < dist_min || (dist == dist_min && t < TOC_out)) {
                        dist_min = dist;
                        P_out = P;
                        Q_out = Q;
                        TOC_out = t;
                    }
                }
            }
            return dist_min;
        }
        template float triTriCCDLinear(const Delta2::common::Triangle<float>&, const Eigen::Matrix4f&, const Eigen::Matrix4f&, const Delta2::common::Triangle<float>&, const Eigen::Matrix4f&, const Eigen::Matrix4f&, Eigen::Vector3f&, Eigen::Vector3f&, float&);
        template double triTriCCDLinear(const Delta2::common::Triangle<double>&, const Eigen::Matrix4d&, const Eigen::Matrix4d&, const Delta2::common::Triangle<double>&, const Eigen::Matrix4d&, const Eigen::Matrix4d&, Eigen::Vector3d&, Eigen::Vector3d&, double&);

        template<int branching, class bucket_real>
		void findContactsBucketContinuousComparison(std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<ContinuousContact<bucket_real>>& hits, bucket_real max_toc, std::vector<int>& pair_used_out) {
            __itt_string_handle* continuous_soup_bucket_task = __itt_string_handle_create("Continuous soup bucket");
            __itt_string_handle* continuous_soup_triangle_task = __itt_string_handle_create("Continuous soup triangle");
            
            float search_dist = a.geo_eps + b.geo_eps;
            pair_used_out.reserve(bucket_pairs.size());
            pair_used_out.clear();

            Eigen::Matrix<bucket_real, 4, 4> a_t_start = a.current_state.getTransformation().cast<bucket_real>();
            Eigen::Matrix<bucket_real, 4, 4> a_t_end = a.future_state.getTransformation().cast<bucket_real>();
            Eigen::Matrix<bucket_real, 4, 4> b_t_start = b.current_state.getTransformation().cast<bucket_real>();
            Eigen::Matrix<bucket_real, 4, 4> b_t_end = b.future_state.getTransformation().cast<bucket_real>();

            std::vector<std::vector<ContinuousContact<bucket_real>>> hits_tmp;
            for (int d = 0; d < bucket_pairs.size(); d++) {
                hits_tmp.push_back({});
            }

            const int parallel_grain_size = 10;
            int max_concurrency = std::min(tbb::info::default_concurrency(), (int)bucket_pairs.size() / parallel_grain_size);
            tbb::task_arena arena(max_concurrency);
            arena.execute([&] {
                tbb::parallel_for(tbb::blocked_range<int>(0, bucket_pairs.size(), parallel_grain_size),
                                [&](tbb::blocked_range<int> r) {
                    for (int d = r.begin(); d < r.end(); d++) {
                        DeferredCompare& dc = bucket_pairs[d];
                        // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, continuous_soup_bucket_task);
                        std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
                        std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

                        const std::vector<double>& a_epss = a_bucket->getEps();
                        const std::vector<double>& a_inner_epss = a_bucket->getInnerEps();
                        const std::vector<double>& b_epss = b_bucket->getEps();
                        const std::vector<double>& b_inner_epss = b_bucket->getInnerEps();

                        // TODO loop over vertices and edges rather than triangles to reduce duplicated checks
                        int a_i = 0;
                        for (common::Triangle<bucket_real>& a_tri : a_bucket->getTriangles()) {
                            double a_eps = a_epss[a_i];
                            double a_inner_eps = a_inner_epss[a_i];
                            common::bbox<bucket_real> a_bbox = a_tri.transformed(a_t_start).bbox().expand(a_tri.transformed(a_t_end).bbox()).expand(a.geo_eps + a_eps);

                            int b_i = 0;
                            for (common::Triangle<bucket_real>& b_tri : b_bucket->getTriangles()) {
                                // __itt_task_begin(globals::itt_handles.detailed_domain, __itt_null, __itt_null, continuous_soup_triangle_task);
                                double b_eps = b_epss[b_i];
                                double b_inner_eps = b_inner_epss[b_i];
                            
                                common::bbox<bucket_real> b_bbox = b_tri.transformed(b_t_start).bbox().expand(b_tri.transformed(b_t_end).bbox()).expand(b.geo_eps + b_eps);

                                if (a_bbox.overlap(b_bbox)) {
                                    Eigen::Vector<bucket_real, 3> P, Q;
                                    bucket_real toc;
                                    bucket_real dist = Delta2::collision::triTriCCDLinear(a_tri, a_t_start, a_t_end, b_tri, b_t_start, b_t_end, P, Q, toc);

                                    if (dist < search_dist && toc <= max_toc && toc < 1.0) {
                                        ContinuousContact<bucket_real> ci(P, Q, a_eps, b_eps, a_inner_eps, b_inner_eps, toc, a, b);
                                        hits_tmp[d].push_back(ci);
                                    }
                                }
                                b_i++;
                                // __itt_task_end(globals::itt_handles.detailed_domain);
                            }
                            a_i++;
                        }
                        // __itt_task_end(globals::itt_handles.detailed_domain);
                    }
                });
            });

            for (int d = 0; d < bucket_pairs.size(); d++) {
                for (int t = 0; t < hits_tmp[d].size(); t++) {
                    hits.push_back(hits_tmp[d][t]);
                }
                pair_used_out.push_back(hits.size());
            }
		};

        // TODO actually needs vectorising
        template<int branching, class real>
		void findContactsBucketConnectedContinuousComparison(std::vector<DeferredCompare>& bucket_pairs, Particle& a, Particle& b, std::vector<ContinuousContact<real>>& hits, double max_toc, std::vector<int>& pair_used_out) {
			int total_vert_tri = 0;
            int total_edge_edge = 0;
            
            real search_dist = a.geo_eps + b.geo_eps;
            pair_used_out.reserve(bucket_pairs.size());
            pair_used_out.clear();

            std::vector<std::vector<ContinuousContact<real>>> hits_tmp;
            for (int d = 0; d < bucket_pairs.size(); d++) {
                hits_tmp.push_back({});
            }

            const int parallel_grain_size = 10;
            int max_concurrency = std::min(tbb::info::default_concurrency(), (int)bucket_pairs.size() / parallel_grain_size);
            tbb::task_arena arena(max_concurrency);
            arena.execute([&] {
                tbb::parallel_for(tbb::blocked_range<int>(0, bucket_pairs.size(), parallel_grain_size),
                                [&](tbb::blocked_range<int> r) {
                    for (int d = r.begin(); d < r.end(); d++) {
                        DeferredCompare& dc = bucket_pairs[d];
                        std::shared_ptr<model::Bucket> a_bucket = a.mesh->getSurrogateTree().getBucket(dc.a.bucket_id);
                        std::shared_ptr<model::Bucket> b_bucket = b.mesh->getSurrogateTree().getBucket(dc.b.bucket_id);

                        Eigen::Matrix<real, 4, 4> a_t_start = a.current_state.getTransformation().cast<real>();
                        Eigen::Matrix<real, 4, 4> a_t_end = a.future_state.getTransformation().cast<real>();
                        Eigen::Matrix<real, 4, 4> b_t_start = b.current_state.getTransformation().cast<real>();
                        Eigen::Matrix<real, 4, 4> b_t_end = b.future_state.getTransformation().cast<real>();

                        std::vector<Eigen::Vector<real, 3>> a_verts = a_bucket->getVertices();
                        std::vector<Eigen::Vector<real, 3>> b_verts = b_bucket->getVertices();

                        std::vector<common::Edge<real>> a_edges = a_bucket->getEdges();
                        std::vector<common::Edge<real>> b_edges = b_bucket->getEdges();

                        // Compare a verts to b tris
                        for (Eigen::Vector<real, 3>& a_vert : a_bucket->getVertices()) {
                            for (common::Triangle<real>& b_tri : b_bucket->getTriangles()) {
                                Eigen::Vector<real, 3> P, Q;
                                real toc;
                                real dist = Delta2::collision::pointTriCCDLinear(a_vert, a_t_start, a_t_end, b_tri, b_t_start, b_t_end, P, Q, toc);

                                if (dist < search_dist && toc <= max_toc && toc < 1.0) {
                                // if (dist < search_dist) {
                                    ContinuousContact<real> ci(P, Q, a.geo_eps, b.geo_eps, 0.0, 0.0, toc, a, b);
                                    hits_tmp[d].push_back(ci);
                                }
                            }
                        }
                        // Compare b verts to a tris
                        for (Eigen::Vector<real, 3>& b_vert : b_bucket->getVertices()) {
                            for (common::Triangle<real>& a_tri : a_bucket->getTriangles()) {
                                Eigen::Vector<real, 3> P, Q;
                                real toc;
                                real dist = Delta2::collision::pointTriCCDLinear(b_vert, b_t_start, b_t_end, a_tri, a_t_start, a_t_end, P, Q, toc);

                                if (dist < search_dist && toc <= max_toc && toc < 1.0) {
                                // if (dist < search_dist) {
                                    ContinuousContact<real> ci(Q, P, a.geo_eps, b.geo_eps, 0.0, 0.0, toc, a, b);
                                    hits_tmp[d].push_back(ci);
                                }
                            }
                        }

                        // Compare a edges to b edges
                        for (common::Edge<real>& a_edge : a_bucket->getEdges()) {
                            Eigen::Vector<real, 3> a_start_A = common::transform(a_edge.A, a_t_start);
                            Eigen::Vector<real, 3> a_start_B = common::transform(a_edge.B, a_t_start);
                            common::Edge<real> a_edge_start(a_start_A, a_start_B);
                            Eigen::Vector<real, 3> a_end_A = common::transform(a_edge.A, a_t_end);
                            Eigen::Vector<real, 3> a_end_B = common::transform(a_edge.B, a_t_end);
                            common::Edge<real> a_edge_end(a_end_A, a_end_B);
                            for (common::Edge<real>& b_edge : b_bucket->getEdges()) {
                                Eigen::Vector<real, 3> b_start_A = common::transform(b_edge.A, b_t_start);
                                Eigen::Vector<real, 3> b_start_B = common::transform(b_edge.B, b_t_start);
                                common::Edge<real> b_edge_start(b_start_A, b_start_B);
                                Eigen::Vector<real, 3> b_end_A = common::transform(b_edge.A, b_t_end);
                                Eigen::Vector<real, 3> b_end_B = common::transform(b_edge.B, b_t_end);
                                common::Edge<real> b_edge_end(b_end_A, b_end_B);
                                
                                Eigen::Vector<real, 3> P, Q;
                                real toc;
                                real dist = Delta2::collision::edgeEdgeCCD(a_edge_start, a_edge_end, b_edge_start, b_edge_end, search_dist, P, Q, toc);

                                if (dist < search_dist && toc <= max_toc && toc < 1.0) {
                                // if (dist < search_dist) {
                                    ContinuousContact<real> ci(Q, P, a.geo_eps, b.geo_eps, 0.0, 0.0, toc, a, b);
                                    hits_tmp[d].push_back(ci);
                                }
                            }
                        }
                    }
                });
            });

            for (int d = 0; d < bucket_pairs.size(); d++) {
                for (int t = 0; t < hits_tmp[d].size(); t++) {
                    hits.push_back(hits_tmp[d][t]);
                }
                pair_used_out.push_back(hits.size());
            }
		};
    }
}
