#include "continuous_comparison.h"

#include "../common/triangle.h"

using namespace Delta2;

double collision::edgeEdgeCCD(const common::Edge<double>& a_start, const common::Edge<double>& a_end, const common::Edge<double>& b_start, const common::Edge<double>& b_end, double search_dist, Eigen::Vector<double, 3>& P_out, Eigen::Vector<double, 3>& Q_out, double& toc_out) {
    // https://physics.stackexchange.com/questions/373534/sweeping-collision-detection-between-two-line-segments-moving-in-3d
    // Eigen::Vector<double, 3> a_min = a_start.min(a_end);
    // Eigen::Vector<double, 3> a_max = a_start.max(a_end);
    // Eigen::Vector<double, 3> b_min = b_start.min(b_end);
    // Eigen::Vector<double, 3> b_max = b_start.max(b_end);
    
    // bool bbox = testAxis(a_min.x(), a_max.x(), b_min.x(), b_max.x(), search_dist);
    // bbox &= testAxis(a_min.y(), a_max.y(), b_min.y(), b_max.y(), search_dist);
    // bbox &= testAxis(a_min.z(), a_max.z(), b_min.z(), b_max.z(), search_dist);

    // if (!bbox) {
    //     return std::numeric_limits<double>::infinity();
    // }

    double min_dist = std::numeric_limits<double>::max();
    
    Eigen::Vector<double, 3> n_start = (a_start.B - a_start.A).cross(b_start.A - a_start.A);
    double dot_start = (b_start.B - a_start.A).dot(n_start);
    Eigen::Vector<double, 3> n_end = (a_end.B - a_end.A).cross(b_end.A - a_end.A);
    double dot_end = (b_end.B - a_end.A).dot(n_end);
    if (dot_start * dot_end < 0.0) {
        // All points are on a plane at some point => there MIGHT be an intersection
        double t_lower = 0.0;
        double t_upper = 1.0;

        for (int iter = 0; iter < 20; iter++) {
            double t_mid = (t_lower + t_upper) / 2.0;

            common::Edge<double> a_t_edge(common::lerp(a_start.A, a_end.A, t_mid), common::lerp(a_start.B, a_end.B, t_mid));
            common::Edge<double> b_t_edge(common::lerp(b_start.A, b_end.A, t_mid), common::lerp(b_start.B, b_end.B, t_mid));
            Eigen::Vector<double, 3> n_mid = (a_t_edge.B - a_t_edge.A).cross(b_t_edge.A - a_t_edge.A);
            double dot_mid = (b_t_edge.B - a_t_edge.A).dot(n_mid);

            if (std::abs(dot_mid) < 1e-8) {
                t_lower = t_mid;
                t_upper = t_mid;
                break;
            }
            else if (dot_start * dot_mid < 0.0) {
                t_upper = t_mid;
            }
            else { // dot_end * dot_mid < 0.0
                t_lower = t_mid;
            }
        }

        double t_mid = (t_lower + t_upper) / 2.0;

        common::Edge<double> a_t_edge(common::lerp(a_start.A, a_end.A, t_mid), common::lerp(a_start.B, a_end.B, t_mid));
        common::Edge<double> b_t_edge(common::lerp(b_start.A, b_end.A, t_mid), common::lerp(b_start.B, b_end.B, t_mid));

        min_dist = a_t_edge.closest(b_t_edge, P_out, Q_out);
        toc_out = t_mid;
    }

    if (min_dist > 1e-8) {
        double t = 0.0;

        for (int iter = 0; iter < 10; iter++) {
            common::Edge<double> a_t_edge(common::lerp(a_start.A, a_end.A, t), common::lerp(a_start.B, a_end.B, t));
            common::Edge<double> b_t_edge(common::lerp(b_start.A, b_end.A, t), common::lerp(b_start.B, b_end.B, t));

            Eigen::Vector<double, 3> P, Q;
            
            double dist = a_t_edge.closest(b_t_edge, P, Q);

            if (dist < min_dist) {
                min_dist = dist;
                P_out = P;
                Q_out = Q;
                toc_out = t;
            }

            double pa = a_t_edge.project(P);
            double pb = b_t_edge.project(Q);

            double t_diff = ((-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)))/2 + (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)))/2 + (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)))/2)/sqrt(std::pow(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)), 2) + std::pow(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)), 2) + std::pow(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)), 2));
            double t_diff2 = ((-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(-a_start.B.x()*pa/2 + b_start.B.x()*pb/2 + (1 - pa)*(a_end.A.x() - a_start.A.x())/2 + (b_end.A.x() - b_start.A.x())*(pb - 1)/2) + (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(-a_start.B.y()*pa/2 + b_start.B.y()*pb/2 + (1 - pa)*(a_end.A.y() - a_start.A.y())/2 + (b_end.A.y() - b_start.A.y())*(pb - 1)/2) + (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(-a_start.B.z()*pa/2 + b_start.B.z()*pb/2 + (1 - pa)*(a_end.A.z() - a_start.A.z())/2 + (b_end.A.z() - b_start.A.z())*(pb - 1)/2))/sqrt(std::pow(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)), 2) + std::pow(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)), 2) + std::pow(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)), 2) ) + (-(-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)))/2 - (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)))/2 - (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)))/2)*((-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)))/2 + (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)))/2 + (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)))/2)/std::pow(std::pow(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)), 2)+std::pow(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)), 2) + std::pow(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)), 2), (3.0/2.0));

            double alpha = 1.0 / (1.0 + iter);

            t = t - t_diff / t_diff2;
            t = common::clamp01(t);
        }
    }
    return min_dist;
}

float collision::edgeEdgeCCD(const common::Edge<float>& a_start, const common::Edge<float>& a_end, const common::Edge<float>& b_start, const common::Edge<float>& b_end, float search_dist, Eigen::Vector<float, 3>& P_out, Eigen::Vector<float, 3>& Q_out, float& toc_out) {
    // https://physics.stackexchange.com/questions/373534/sweeping-collision-detection-between-two-line-segments-moving-in-3d
    // Eigen::Vector<float, 3> a_min = a_start.min(a_end);
    // Eigen::Vector<float, 3> a_max = a_start.max(a_end);
    // Eigen::Vector<float, 3> b_min = b_start.min(b_end);
    // Eigen::Vector<float, 3> b_max = b_start.max(b_end);
    
    // bool bbox = testAxis(a_min.x(), a_max.x(), b_min.x(), b_max.x(), search_dist);
    // bbox &= testAxis(a_min.y(), a_max.y(), b_min.y(), b_max.y(), search_dist);
    // bbox &= testAxis(a_min.z(), a_max.z(), b_min.z(), b_max.z(), search_dist);

    // if (!bbox) {
    //     return std::numeric_limits<float>::infinity();
    // }

    float min_dist = std::numeric_limits<float>::max();
    
    Eigen::Vector<float, 3> n_start = (a_start.B - a_start.A).cross(b_start.A - a_start.A);
    float dot_start = (b_start.B - a_start.A).dot(n_start);
    Eigen::Vector<float, 3> n_end = (a_end.B - a_end.A).cross(b_end.A - a_end.A);
    float dot_end = (b_end.B - a_end.A).dot(n_end);
    if (dot_start * dot_end < 0.0) {
        // All points are on a plane at some point => there MIGHT be an intersection
        float t_lower = 0.0;
        float t_upper = 1.0;

        for (int iter = 0; iter < 20; iter++) {
            float t_mid = (t_lower + t_upper) / 2.0;

            common::Edge<float> a_t_edge(common::lerp(a_start.A, a_end.A, t_mid), common::lerp(a_start.B, a_end.B, t_mid));
            common::Edge<float> b_t_edge(common::lerp(b_start.A, b_end.A, t_mid), common::lerp(b_start.B, b_end.B, t_mid));
            Eigen::Vector<float, 3> n_mid = (a_t_edge.B - a_t_edge.A).cross(b_t_edge.A - a_t_edge.A);
            float dot_mid = (b_t_edge.B - a_t_edge.A).dot(n_mid);

            if (std::abs(dot_mid) < 1e-8) {
                t_lower = t_mid;
                t_upper = t_mid;
                break;
            }
            else if (dot_start * dot_mid < 0.0) {
                t_upper = t_mid;
            }
            else { // dot_end * dot_mid < 0.0
                t_lower = t_mid;
            }
        }

        float t_mid = (t_lower + t_upper) / 2.0;

        common::Edge<float> a_t_edge(common::lerp(a_start.A, a_end.A, t_mid), common::lerp(a_start.B, a_end.B, t_mid));
        common::Edge<float> b_t_edge(common::lerp(b_start.A, b_end.A, t_mid), common::lerp(b_start.B, b_end.B, t_mid));

        min_dist = a_t_edge.closest(b_t_edge, P_out, Q_out);
        toc_out = t_mid;
    }

    if (min_dist > 1e-8) {
        float t = 0.0;

        for (int iter = 0; iter < 10; iter++) {
            common::Edge<float> a_t_edge(common::lerp(a_start.A, a_end.A, t), common::lerp(a_start.B, a_end.B, t));
            common::Edge<float> b_t_edge(common::lerp(b_start.A, b_end.A, t), common::lerp(b_start.B, b_end.B, t));

            Eigen::Vector<float, 3> P, Q;
            
            float dist = a_t_edge.closest(b_t_edge, P, Q);

            if (dist < min_dist) {
                min_dist = dist;
                P_out = P;
                Q_out = Q;
                toc_out = t;
            }

            float pa = a_t_edge.project(P);
            float pb = b_t_edge.project(Q);

            float t_diff = ((-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)))/2 + (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)))/2 + (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)))/2)/sqrt(std::pow(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)), 2) + std::pow(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)), 2) + std::pow(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)), 2));
            float t_diff2 = ((-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(-a_start.B.x()*pa/2 + b_start.B.x()*pb/2 + (1 - pa)*(a_end.A.x() - a_start.A.x())/2 + (b_end.A.x() - b_start.A.x())*(pb - 1)/2) + (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(-a_start.B.y()*pa/2 + b_start.B.y()*pb/2 + (1 - pa)*(a_end.A.y() - a_start.A.y())/2 + (b_end.A.y() - b_start.A.y())*(pb - 1)/2) + (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(-a_start.B.z()*pa/2 + b_start.B.z()*pb/2 + (1 - pa)*(a_end.A.z() - a_start.A.z())/2 + (b_end.A.z() - b_start.A.z())*(pb - 1)/2))/sqrt(std::pow(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)), 2) + std::pow(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)), 2) + std::pow(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)), 2) ) + (-(-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)))/2 - (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)))/2 - (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)))/2)*((-2*a_start.B.x()*pa + 2*b_start.B.x()*pb + 2*(1 - pa)*(a_end.A.x() - a_start.A.x()) + 2*(b_end.A.x() - b_start.A.x())*(pb - 1))*(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)))/2 + (-2*a_start.B.y()*pa + 2*b_start.B.y()*pb + 2*(1 - pa)*(a_end.A.y() - a_start.A.y()) + 2*(b_end.A.y() - b_start.A.y())*(pb - 1))*(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)))/2 + (-2*a_start.B.z()*pa + 2*b_start.B.z()*pb + 2*(1 - pa)*(a_end.A.z() - a_start.A.z()) + 2*(b_end.A.z() - b_start.A.z())*(pb - 1))*(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)))/2)/std::pow(std::pow(pa*(a_end.B.x() + a_start.B.x()*(1 - t)) - pb*(b_end.B.x() + b_start.B.x()*(1 - t)) + (1 - pa)*(a_end.A.x()*t + a_start.A.x()*(1 - t)) - (1 - pb)*(b_end.A.x()*t + b_start.A.x()*(1 - t)), 2)+std::pow(pa*(a_end.B.y() + a_start.B.y()*(1 - t)) - pb*(b_end.B.y() + b_start.B.y()*(1 - t)) + (1 - pa)*(a_end.A.y()*t + a_start.A.y()*(1 - t)) - (1 - pb)*(b_end.A.y()*t + b_start.A.y()*(1 - t)), 2) + std::pow(pa*(a_end.B.z() + a_start.B.z()*(1 - t)) - pb*(b_end.B.z() + b_start.B.z()*(1 - t)) + (1 - pa)*(a_end.A.z()*t + a_start.A.z()*(1 - t)) - (1 - pb)*(b_end.A.z()*t + b_start.A.z()*(1 - t)), 2), (3.0/2.0));

            float alpha = 1.0 / (1.0 + iter);

            t = t - t_diff / t_diff2;
            t = common::clamp01(t);
        }
    }
    return min_dist;
}