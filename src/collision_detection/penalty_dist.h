#pragma once

#include <Eigen/Dense>

#include "../common/triangle.h"
#include "../common/aligned.h"
#include "contact.h"

namespace Delta2 {
	namespace collision {
        template<typename real>
        real Heaviside(real x) {
            return x<0 ? 0.0 : 1.0;
        }

        template<typename real, int iterations>
		std::vector<std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>>> _penalty1To1Dists(std::vector<Delta2::common::Triangle<real>>& a_tris, std::vector<Delta2::common::Triangle<real>>& b_tris, const std::vector<real>& max_error, std::vector<bool>& failed_out) {
            std::vector<std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>>> result;
            result.reserve(b_tris.size());
            
            const int TrianglePairs = 8;

            assert(a_tris.size() == b_tris.size());
            int size = a_tris.size();

            ALIGNED(real x_a11[TrianglePairs], 64);
            ALIGNED(real x_a12[TrianglePairs], 64);
            ALIGNED(real x_a13[TrianglePairs], 64);
            ALIGNED(real x_a21[TrianglePairs], 64);
            ALIGNED(real x_a22[TrianglePairs], 64);
            ALIGNED(real x_a23[TrianglePairs], 64);
            ALIGNED(real x_a31[TrianglePairs], 64);
            ALIGNED(real x_a32[TrianglePairs], 64);
            ALIGNED(real x_a33[TrianglePairs], 64);

            ALIGNED(real x_b11[TrianglePairs], 64);
            ALIGNED(real x_b12[TrianglePairs], 64);
            ALIGNED(real x_b13[TrianglePairs], 64);
            ALIGNED(real x_b21[TrianglePairs], 64);
            ALIGNED(real x_b22[TrianglePairs], 64);
            ALIGNED(real x_b23[TrianglePairs], 64);
            ALIGNED(real x_b31[TrianglePairs], 64);
            ALIGNED(real x_b32[TrianglePairs], 64);
            ALIGNED(real x_b33[TrianglePairs], 64);

            ALIGNED(real a[TrianglePairs], 64);
            ALIGNED(real b[TrianglePairs], 64);
            ALIGNED(real c[TrianglePairs], 64);
            ALIGNED(real d[TrianglePairs], 64);

            ALIGNED(real a_last[TrianglePairs], 64);
            ALIGNED(real b_last[TrianglePairs], 64);
            ALIGNED(real c_last[TrianglePairs], 64);
            ALIGNED(real d_last[TrianglePairs], 64);

            ALIGNED(real a_new[TrianglePairs], 64);
            ALIGNED(real b_new[TrianglePairs], 64);
            ALIGNED(real c_new[TrianglePairs], 64);
            ALIGNED(real d_new[TrianglePairs], 64);

            for (int batch = 0; batch < (size + TrianglePairs - 1) / TrianglePairs; batch++) {
                for (int t = 0; t < TrianglePairs; t++) {
                    int ind = std::min(TrianglePairs * batch + t, size - 1);
                    x_a11[t] = a_tris[ind].A.x();
                    x_a12[t] = a_tris[ind].A.y();
                    x_a13[t] = a_tris[ind].A.z();
                    x_a21[t] = a_tris[ind].B.x();
                    x_a22[t] = a_tris[ind].B.y();
                    x_a23[t] = a_tris[ind].B.z();
                    x_a31[t] = a_tris[ind].C.x();
                    x_a32[t] = a_tris[ind].C.y();
                    x_a33[t] = a_tris[ind].C.z();

                    x_b11[t] = b_tris[ind].A.x();
                    x_b12[t] = b_tris[ind].A.y();
                    x_b13[t] = b_tris[ind].A.z();
                    x_b21[t] = b_tris[ind].B.x();
                    x_b22[t] = b_tris[ind].B.y();
                    x_b23[t] = b_tris[ind].B.z();
                    x_b31[t] = b_tris[ind].C.x();
                    x_b32[t] = b_tris[ind].C.y();
                    x_b33[t] = b_tris[ind].C.z();
                }

                for (int i = 0; i < TrianglePairs; i++) {
                    a[i] = 0.33;
                    b[i] = 0.33;
                    c[i] = 0.33;
                    d[i] = 0.33;
                }

                real alpha = 1;
                real C_dontdividebyzero = 1e-5;
                real C_damp_start = 1.0;
                real C_damp_end = 1.0;

                for (int i = 0; i < iterations; i++) {
                    real C_damp = C_damp_start - (C_damp_start - C_damp_end) * i / (iterations - 1);
                    #pragma omp simd
                    for (int t = 0; t < TrianglePairs; t++) {
                        a_last[t] = a[t];
                        b_last[t] = b[t];
                        c_last[t] = c[t];
                        d_last[t] = d[t];

                        a_new[t] = a[t] - 1.0f * ((x_a11[t] - x_a21[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) + (x_a12[t] - x_a22[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) + (x_a13[t] - x_a23[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + (x_a11[t] - x_a21[t]) * (x_a11[t] - x_a21[t]) + (x_a12[t] - x_a22[t]) * (x_a12[t] - x_a22[t]) + (x_a13[t] - x_a23[t]) * (x_a13[t] - x_a23[t]));
                        b_new[t] = b[t] - 1.0f * ((x_a11[t] - x_a31[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) + (x_a12[t] - x_a32[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) + (x_a13[t] - x_a33[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + (x_a11[t] - x_a31[t]) * (x_a11[t] - x_a31[t]) + (x_a12[t] - x_a32[t]) * (x_a12[t] - x_a32[t]) + (x_a13[t] - x_a33[t]) * (x_a13[t] - x_a33[t]));
                        c_new[t] = c[t] - 1.0f * (-1.0f * (x_b11[t] - x_b21[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) - 1.0f * (x_b12[t] - x_b22[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) - 1.0f * (x_b13[t] - x_b23[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b21[t]) * (x_b11[t] - x_b21[t]) + 1.0f * (x_b12[t] - x_b22[t]) * (x_b12[t] - x_b22[t]) + 1.0f * (x_b13[t] - x_b23[t]) * (x_b13[t] - x_b23[t]));
                        d_new[t] = d[t] - 1.0f * (-1.0f * (x_b11[t] - x_b31[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) - 1.0f * (x_b12[t] - x_b32[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) - 1.0f * (x_b13[t] - x_b33[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b31[t]) * (x_b11[t] - x_b31[t]) + 1.0f * (x_b12[t] - x_b32[t]) * (x_b12[t] - x_b32[t]) + 1.0f * (x_b13[t] - x_b33[t]) * (x_b13[t] - x_b33[t]));

                        for (int i = 0; i < 4; i++) {
                            a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
                            b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
                            c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));
                            d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));

                            a_new[t] = a[t] + alpha * (-a[t] * Heaviside(-a[t]) - 0.5f * (a[t] + b[t] - 1.0f) * Heaviside(a[t] + b[t] - 1.0f));
                            b_new[t] = b[t] + alpha * (-b[t] * Heaviside(-b[t]) - 0.5f * (a[t] + b[t] - 1.0f) * Heaviside(a[t] + b[t] - 1.0f));
                            c_new[t] = c[t] + alpha * (-c[t] * Heaviside(-c[t]) - 0.5f * (c[t] + d[t] - 1.0f) * Heaviside(c[t] + d[t] - 1.0f));
                            d_new[t] = d[t] + alpha * (-d[t] * Heaviside(-d[t]) - 0.5f * (c[t] + d[t] - 1.0f) * Heaviside(c[t] + d[t] - 1.0f));
                        }

                        a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
                        b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
                        c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));
                        d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));
                    }
                }

                for (int t = 0; t < TrianglePairs && TrianglePairs * batch + t < size; t++) {
                    int ind = TrianglePairs * batch + t;

                    a[t] = common::clamp01(a[t]);
                    b[t] = common::clamp01(b[t]);
                    c[t] = common::clamp01(c[t]);
                    d[t] = common::clamp01(d[t]);

                    if (1.0f < a[t] + b[t]) {
                        a_new[t] = a[t] / (a[t] + b[t]);
                        b_new[t] = b[t] / (a[t] + b[t]);
                        a[t] = a_new[t];
                        b[t] = b_new[t];
                    }

                    if (1.0f < c[t] + d[t]) {
                        c_new[t] = c[t] / (c[t] + d[t]);
                        d_new[t] = d[t] / (c[t] + d[t]);
                        c[t] = c_new[t];
                        d[t] = d_new[t];
                    }

                    real Ax_last_update = (a[t] - a_last[t]) * (x_a21[t] - x_a11[t]) + (b[t] - b_last[t]) * (x_a31[t] - x_a11[t]);
                    real Ay_last_update = (a[t] - a_last[t]) * (x_a22[t] - x_a12[t]) + (b[t] - b_last[t]) * (x_a32[t] - x_a12[t]);
                    real Az_last_update = (a[t] - a_last[t]) * (x_a23[t] - x_a13[t]) + (b[t] - b_last[t]) * (x_a33[t] - x_a13[t]);

                    real Bx_last_update = (c[t] - c_last[t]) * (x_b21[t] - x_b11[t]) + (d[t] - d_last[t]) * (x_b31[t] - x_b11[t]);
                    real By_last_update = (c[t] - c_last[t]) * (x_b22[t] - x_b12[t]) + (d[t] - d_last[t]) * (x_b32[t] - x_b12[t]);
                    real Bz_last_update = (c[t] - c_last[t]) * (x_b23[t] - x_b13[t]) + (d[t] - d_last[t]) * (x_b33[t] - x_b13[t]);

                    real A_all_update = sqrt(Ax_last_update * Ax_last_update + Ay_last_update * Ay_last_update + Az_last_update * Az_last_update) * 1 / C_damp_end;
                    real B_all_update = sqrt(Bx_last_update * Bx_last_update + By_last_update * By_last_update + Bz_last_update * Bz_last_update) * 1 / C_damp_end;

                    failed_out[ind] = A_all_update > max_error[ind] || B_all_update > max_error[ind];
                    Eigen::Vector<real, 3> P = { x_a11[t] + a[t] * (x_a21[t] - x_a11[t]) + b[t] * (x_a31[t] - x_a11[t]), x_a12[t] + a[t] * (x_a22[t] - x_a12[t]) + b[t] * (x_a32[t] - x_a12[t]), x_a13[t] + a[t] * (x_a23[t] - x_a13[t]) + b[t] * (x_a33[t] - x_a13[t]) };
                    Eigen::Vector<real, 3> Q = { x_b11[t] + c[t] * (x_b21[t] - x_b11[t]) + d[t] * (x_b31[t] - x_b11[t]), x_b12[t] + c[t] * (x_b22[t] - x_b12[t]) + d[t] * (x_b32[t] - x_b12[t]), x_b13[t] + c[t] * (x_b23[t] - x_b13[t]) + d[t] * (x_b33[t] - x_b13[t]) };
                    result.push_back(std::make_tuple((P-Q).norm(), P, Q));
                }
            }
            return result;
        }
        template std::vector<std::tuple<float, Eigen::Vector<float, 3>, Eigen::Vector<float, 3>>> _penalty1To1Dists<float, 8>(std::vector<Delta2::common::Triangle<float>>&, std::vector<Delta2::common::Triangle<float>>&, const std::vector<float>&, std::vector<bool>&);
        template std::vector<std::tuple<double, Eigen::Vector<double, 3>, Eigen::Vector<double, 3>>> _penalty1To1Dists<double, 8>(std::vector<Delta2::common::Triangle<double>>&, std::vector<Delta2::common::Triangle<double>>&, const std::vector<double>&, std::vector<bool>&);

        template<typename real, int iterations>
		std::vector<std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>>> _penaltyInnerDists(Delta2::common::Triangle<real>& a_tri, std::vector<Delta2::common::Triangle<real>>& b_tris, const std::vector<real>& max_error, std::vector<bool>& failed_out) {
            std::vector<std::tuple<real, Eigen::Vector<real, 3>, Eigen::Vector<real, 3>>> result;
            result.reserve(b_tris.size());
            
            const int ProjectCorrectIterations = 4;
            const int TrianglePairs = 8;

            const real x_a11 = a_tri.A.x();
            const real x_a12 = a_tri.A.y();
            const real x_a13 = a_tri.A.z();
            const real x_a21 = a_tri.B.x();
            const real x_a22 = a_tri.B.y();
            const real x_a23 = a_tri.B.z();
            const real x_a31 = a_tri.C.x();
            const real x_a32 = a_tri.C.y();
            const real x_a33 = a_tri.C.z();

            ALIGNED(real x_b11[TrianglePairs], 64);
            ALIGNED(real x_b12[TrianglePairs], 64);
            ALIGNED(real x_b13[TrianglePairs], 64);
            ALIGNED(real x_b21[TrianglePairs], 64);
            ALIGNED(real x_b22[TrianglePairs], 64);
            ALIGNED(real x_b23[TrianglePairs], 64);
            ALIGNED(real x_b31[TrianglePairs], 64);
            ALIGNED(real x_b32[TrianglePairs], 64);
            ALIGNED(real x_b33[TrianglePairs], 64);

            for (int t = 0; t < TrianglePairs; t++) {
                x_b11[t] = b_tris[t % b_tris.size()].A.x();
                x_b12[t] = b_tris[t % b_tris.size()].A.y();
                x_b13[t] = b_tris[t % b_tris.size()].A.z();
                x_b21[t] = b_tris[t % b_tris.size()].B.x();
                x_b22[t] = b_tris[t % b_tris.size()].B.y();
                x_b23[t] = b_tris[t % b_tris.size()].B.z();
                x_b31[t] = b_tris[t % b_tris.size()].C.x();
                x_b32[t] = b_tris[t % b_tris.size()].C.y();
                x_b33[t] = b_tris[t % b_tris.size()].C.z();
            }

            ALIGNED(real a[TrianglePairs], 64);
            ALIGNED(real b[TrianglePairs], 64);
            ALIGNED(real c[TrianglePairs], 64);
            ALIGNED(real d[TrianglePairs], 64);

            ALIGNED(real a_last[TrianglePairs], 64);
            ALIGNED(real b_last[TrianglePairs], 64);
            ALIGNED(real c_last[TrianglePairs], 64);
            ALIGNED(real d_last[TrianglePairs], 64);

            ALIGNED(real a_new[TrianglePairs], 64);
            ALIGNED(real b_new[TrianglePairs], 64);
            ALIGNED(real c_new[TrianglePairs], 64);
            ALIGNED(real d_new[TrianglePairs], 64);

            for (int i = 0; i < TrianglePairs; i++) {
                a[i] = 0.33;
                b[i] = 0.33;
                c[i] = 0.33;
                d[i] = 0.33;
            }

            ALIGNED(real failed[TrianglePairs], 64);

            real alpha = 1;
            real C_dontdividebyzero = 1e-5;
            real C_damp_start = 1.0;
            real C_damp_end = 0.5;

            for (int i = 0; i < iterations; i++) {
                real C_damp = C_damp_start - (C_damp_start - C_damp_end) * i / (iterations - 1);
                #pragma omp simd
                for (int t = 0; t < TrianglePairs; t++) {
                    a_last[t] = a[t];
                    b_last[t] = b[t];
                    c_last[t] = c[t];
                    d_last[t] = d[t];

                    a_new[t] = a[t] - C_damp * ((x_a11 - x_a21) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) + (x_a12 - x_a22) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) + (x_a13 - x_a23) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + (x_a11 - x_a21) * (x_a11 - x_a21) + (x_a12 - x_a22) * (x_a12 - x_a22) + (x_a13 - x_a23) * (x_a13 - x_a23));
                    b_new[t] = b[t] - C_damp * ((x_a11 - x_a31) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) + (x_a12 - x_a32) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) + (x_a13 - x_a33) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + (x_a11 - x_a31) * (x_a11 - x_a31) + (x_a12 - x_a32) * (x_a12 - x_a32) + (x_a13 - x_a33) * (x_a13 - x_a33));
                    c_new[t] = c[t] - C_damp * (-1.0f * (x_b11[t] - x_b21[t]) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) - 1.0f * (x_b12[t] - x_b22[t]) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) - 1.0f * (x_b13[t] - x_b23[t]) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b21[t]) * (x_b11[t] - x_b21[t]) + 1.0f * (x_b12[t] - x_b22[t]) * (x_b12[t] - x_b22[t]) + 1.0f * (x_b13[t] - x_b23[t]) * (x_b13[t] - x_b23[t]));
                    d_new[t] = d[t] - C_damp * (-1.0f * (x_b11[t] - x_b31[t]) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) - 1.0f * (x_b12[t] - x_b32[t]) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) - 1.0f * (x_b13[t] - x_b33[t]) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b31[t]) * (x_b11[t] - x_b31[t]) + 1.0f * (x_b12[t] - x_b32[t]) * (x_b12[t] - x_b32[t]) + 1.0f * (x_b13[t] - x_b33[t]) * (x_b13[t] - x_b33[t]));

                    for (int i = 0; i < ProjectCorrectIterations; i++) {
                        a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
                        b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
                        c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
                        d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));

                        a_new[t] = a[t] + alpha * (-a[t] * Heaviside(-a[t]) - 0.5f * (a[t] + b[t] - 1) * Heaviside(a[t] + b[t] - 1));
                        b_new[t] = b[t] + alpha * (-b[t] * Heaviside(-b[t]) - 0.5f * (a[t] + b[t] - 1) * Heaviside(a[t] + b[t] - 1));
                        c_new[t] = c[t] + alpha * (-c[t] * Heaviside(-c[t]) - 0.5f * (c[t] + d[t] - 1) * Heaviside(c[t] + d[t] - 1));
                        d_new[t] = d[t] + alpha * (-d[t] * Heaviside(-d[t]) - 0.5f * (c[t] + d[t] - 1) * Heaviside(c[t] + d[t] - 1));
                    }
                    a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
                    b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
                    c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
                    d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
                }
            }

            for (int t = 0; t < b_tris.size(); t++) {
                a[t] = common::clamp01(a[t]);
                b[t] = common::clamp01(b[t]);
                c[t] = common::clamp01(c[t]);
                d[t] = common::clamp01(d[t]);

                if (1.0f < a[t] + b[t]) {
                    a[t] = a[t] / (a[t] + b[t]);
                    b[t] = b[t] / (a[t] + b[t]);
                }

                if (1.0f < c[t] + d[t]) {
                    c[t] = c[t] / (c[t] + d[t]);
                    d[t] = d[t] / (c[t] + d[t]);
                }

                real Ax_last_update = (a[t] - a_last[t]) * (x_a21 - x_a11) + (b[t] - b_last[t]) * (x_a31 - x_a11);
                real Ay_last_update = (a[t] - a_last[t]) * (x_a22 - x_a12) + (b[t] - b_last[t]) * (x_a32 - x_a12);
                real Az_last_update = (a[t] - a_last[t]) * (x_a23 - x_a13) + (b[t] - b_last[t]) * (x_a33 - x_a13);

                real Bx_last_update = (c[t] - c_last[t]) * (x_b21[t] - x_b11[t]) + (d[t] - d_last[t]) * (x_b31[t] - x_b11[t]);
                real By_last_update = (c[t] - c_last[t]) * (x_b22[t] - x_b12[t]) + (d[t] - d_last[t]) * (x_b32[t] - x_b12[t]);
                real Bz_last_update = (c[t] - c_last[t]) * (x_b23[t] - x_b13[t]) + (d[t] - d_last[t]) * (x_b33[t] - x_b13[t]);

                real A_all_update = sqrt(Ax_last_update * Ax_last_update + Ay_last_update * Ay_last_update + Az_last_update * Az_last_update) * 1 / C_damp_end;
                real B_all_update = sqrt(Bx_last_update * Bx_last_update + By_last_update * By_last_update + Bz_last_update * Bz_last_update) * 1 / C_damp_end;

                failed[t] = A_all_update > max_error[t] || B_all_update > max_error[t];
            }

            for (int t = 0; t < b_tris.size(); t++) {
                Eigen::Vector<real, 3> P = { x_a11 + a[t] * (x_a21 - x_a11) + b[t] * (x_a31 - x_a11), x_a12 + a[t] * (x_a22 - x_a12) + b[t] * (x_a32 - x_a12), x_a13 + a[t] * (x_a23 - x_a13) + b[t] * (x_a33 - x_a13) };
                Eigen::Vector<real, 3> Q = { x_b11[t] + c[t] * (x_b21[t] - x_b11[t]) + d[t] * (x_b31[t] - x_b11[t]), x_b12[t] + c[t] * (x_b22[t] - x_b12[t]) + d[t] * (x_b32[t] - x_b12[t]), x_b13[t] + c[t] * (x_b23[t] - x_b13[t]) + d[t] * (x_b33[t] - x_b13[t]) };
                failed_out[t] = failed[t];
                result.push_back(std::make_tuple((P-Q).norm(), P, Q));
            }
            return result;
        }
        template std::vector<std::tuple<float, Eigen::Vector<float, 3>, Eigen::Vector<float, 3>>> _penaltyInnerDists<float, 8>(Delta2::common::Triangle<float>&, std::vector<Delta2::common::Triangle<float>>&, const std::vector<float>&, std::vector<bool>&);
        template std::vector<std::tuple<double, Eigen::Vector<double, 3>, Eigen::Vector<double, 3>>> _penaltyInnerDists<double, 8>(Delta2::common::Triangle<double>&, std::vector<Delta2::common::Triangle<double>>&, const std::vector<double>&, std::vector<bool>&);
    }
}
