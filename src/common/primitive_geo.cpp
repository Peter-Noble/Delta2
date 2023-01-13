#include "primitive_geo.h"
#include "igl/PI.h"

namespace Delta2 {
    namespace common {
        void cube(double r, Eigen::Matrix<double, -1, -1>& V, Eigen::MatrixXi& F) {
            V = Eigen::Matrix<double, -1, -1>({
                { r,  r,  r},
                { r,  r, -r},
                { r, -r,  r},
                { r, -r, -r},
                {-r,  r,  r},
                {-r,  r, -r},
                {-r, -r,  r},
                {-r, -r, -r},
            });
            F = Eigen::MatrixXi({
                {0, 2, 1},
                {1, 2, 3},
                {2, 6, 3},
                {3, 6, 7},
                {4, 5, 6},
                {5, 7, 6},
                {4, 0, 1},
                {4, 1, 5},
                {0, 4, 6},
                {0, 6, 2},
                {3, 7, 1},
                {7, 5, 1}
            });
        }

        void sphere(const double& r,
                    const int res,
                    Eigen::MatrixXd& V,
                    Eigen::MatrixXi& F)
        {
            if (res < 3) {
                throw std::invalid_argument("A sphere with fewer than 3 rings doesn't make sense");
            }
            V.resize((res-2)*res+2, 3);
            F.resize((res-3)*2*res+2*res, 3);
            
            V.row(0) << 0.0, 0.0, r;
            for (int j = 1; j < res - 1; j++){
                double z = r * cos(igl::PI * (double)j / (double(res - 1)));
                for (int k = 0; k < res; k++){
                    double x = r * sin(igl::PI * (double)j / (double(res - 1))) * cos(2 * igl::PI * (double)k / (double(res)));
                    double y = r * sin(igl::PI * (double)j / (double(res - 1))) * sin(2 * igl::PI * (double)k / (double(res)));
                    V.row(1+(j-1)*res+k) << x, y, z;
                }
            }
            V.row((res-2)*res+2-1) << 0.0, 0.0, -r;

            for (int j = 0; j < res; j++) {
                F.row(j) << 0, j+1, (j+1)%res+1;
            }

            for (int j = 1; j < res - 2; j++) {
                for (int k = 0; k < res; k++) {
                    int v1 = (j - 1) * res + k + 1;
                    int v2 = j * res + k + 1;
                    int v3 = (j - 1) * res + 1 + (k + 1) % res;
                    int v4 = j * res + 1 + (k + 1) % res;
                    F.row(res + 2*((j-1)*res + k)) << v1, v4, v3;
                    F.row(res + 2*((j-1)*res + k) + 1) << v1, v2, v4;
                }
            }

            int face_offset = res + (res - 3) * res * 2;
            int vert_offset = 1 + (res - 3) * res;
            for (int j = 0; j < res; j++) {
                F.row(face_offset + j) << vert_offset+res, vert_offset+(j+1)%res, vert_offset+j;
            }
        }
        void plane(int width, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
            V.resize(4, 3);
            F.resize(2, 3);
            V = Eigen::MatrixXd({
                {  width / 2.0,  width / 2.0,  0.0 },
                {  width / 2.0, -width / 2.0,  0.0 },
                { -width / 2.0,  width / 2.0,  0.0 },
                { -width / 2.0, -width / 2.0,  0.0 }
            });
            F = Eigen::MatrixXi({
                {0, 2, 1},
                {1, 2, 3}
            });
        }
    }
}
