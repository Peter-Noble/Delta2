#pragma once

#include <Eigen/Dense>

#include "../model/particle.h"

namespace Delta2 {
    namespace collision{
        template<typename real>
        struct Contact {
            Contact(const Eigen::Vector<real, 3>& A, const Eigen::Vector<real, 3>& B, real eps_a, real eps_b, real eps_inner_a, real eps_inner_b, Particle& particle_a, Particle& particle_b) : A(A), B(B), eps_a(eps_a), eps_b(eps_b), eps_inner_a(eps_inner_a), eps_inner_b(eps_inner_b), p_a(&particle_a), p_b(&particle_b) {
                force_mag = 0.0;
            };
            bool operator < (const Contact& other) const {
                double this_proportion = (A - B).norm() / (eps_a + eps_b);
                double other_proportion = (other.A - other.B).norm() / (other.eps_a + other.eps_b);
                return this_proportion < other_proportion;
            }
            Eigen::Vector<real, 3> A;
            Eigen::Vector<real, 3> B;
            real eps_a;
            real eps_b;
            real eps_inner_a;
            real eps_inner_b;
            Particle* p_a;
            Particle* p_b;
            real force_mag;
        };
        template struct Contact<float>;
        template struct Contact<double>;

        template<typename real>
        bool compareContacts(const Contact<real>& a, const Contact<real>& b) {
            double a_eps = a.eps_a + a.eps_b;
            double b_eps = b.eps_a + b.eps_b;
            return a_eps < b_eps;
        }
        template bool compareContacts<float>(const Contact<float>& a, const Contact<float>& b);
        template bool compareContacts<double>(const Contact<double>& a, const Contact<double>& b);

        template <typename real>
        std::vector<Contact<real>> filterContacts(std::vector<Contact<real>>& unfiltered, real extra_filter_margin) {
            std::vector<Contact<real>> filtered;
            std::sort(unfiltered.begin(), unfiltered.end(), compareContacts<real>);
            for (const Contact<real>& ct : unfiltered) {
                // This goes smallest to largest
                bool new_point = true;
                for (Contact<real>& h : filtered) {
                    Eigen::Vector<real, 3> ht_mid = (ct.A + ct.B) / 2;
                    Eigen::Vector<real, 3> h_mid = (h.A + h.B) / 2;
                    Eigen::Vector<real, 3> diff = ht_mid - h_mid;
                    real dist = diff.norm();

                    real search = ct.eps_a + ct.eps_b + extra_filter_margin;

                    if (dist < search) {
                        new_point = false;
                        break;
                    }
                }
                if (new_point) {
                    filtered.push_back(ct);
                }
            }
            return filtered;
        }
        template std::vector<Contact<float>> filterContacts(std::vector<Contact<float>>& unfiltered, float extra_filter_margin);
        template std::vector<Contact<double>> filterContacts(std::vector<Contact<double>>& unfiltered, double extra_filter_margin);
        
        template<typename real>
        struct ContinuousContact {
            ContinuousContact(const Eigen::Vector<real, 3>& A, const Eigen::Vector<real, 3>& B, real eps_a, real eps_b, real eps_inner_a, real eps_inner_b, real t, Particle& particle_a, Particle& particle_b) : A(A), B(B), eps_a(eps_a), eps_b(eps_b), eps_inner_a(eps_inner_a), eps_inner_b(eps_inner_b), toc(t), p_a(&particle_a), p_b(&particle_b) {};
            Eigen::Vector<real, 3> A;
            Eigen::Vector<real, 3> B;
            real eps_a;
            real eps_b;
            real eps_inner_a;
            real eps_inner_b;
            real toc;
            Particle* p_a;
            Particle* p_b;
        };
        template struct ContinuousContact<float>;
        template struct ContinuousContact<double>;

        template<typename real>
        bool compareContinuousContacts(const ContinuousContact<real>& a, const ContinuousContact<real>& b) {
            double a_eps = a.eps_a + a.eps_b;
            double b_eps = b.eps_a + b.eps_b;
            return a_eps < b_eps;
        }
        template bool compareContinuousContacts<float>(const ContinuousContact<float>& a, const ContinuousContact<float>& b);
        template bool compareContinuousContacts<double>(const ContinuousContact<double>& a, const ContinuousContact<double>& b);

        template <typename real>
        std::vector<ContinuousContact<real>> filterContinuousContacts(std::vector<ContinuousContact<real>> unfiltered, real extra_filter_margin) {
            std::vector<ContinuousContact<real>> filtered;
            std::sort(unfiltered.begin(), unfiltered.end(), compareContinuousContacts<real>);
            for (const ContinuousContact<real>& ct : unfiltered) {
                // This goes smallest to largest
                bool new_point = true;
                for (ContinuousContact<real>& h : filtered) {
                    Eigen::Vector<real, 3> ht_mid = (ct.A + ct.B) / 2;
                    Eigen::Vector<real, 3> h_mid = (h.A + h.B) / 2;
                    Eigen::Vector<real, 3> diff = ht_mid - h_mid;
                    real dist = diff.norm();

                    real search = ct.eps_a + ct.eps_b + extra_filter_margin;

                    if (dist < search) {
                        new_point = false;
                        break;
                    }
                }
                if (new_point) {
                    filtered.push_back(ct);
                }
            }
            return filtered;
        }
        template std::vector<ContinuousContact<float>> filterContinuousContacts(std::vector<ContinuousContact<float>> unfiltered, float extra_filter_margin);
        template std::vector<ContinuousContact<double>> filterContinuousContacts(std::vector<ContinuousContact<double>> unfiltered, double extra_filter_margin);
    }
}
