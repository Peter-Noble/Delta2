#include "mesh.h"

#include "igl/tinyply.h"
#include "igl/readOBJ.h"
#include "igl/copyleft/tetgen/tetrahedralize.h"
#include "igl/volume.h"
#include "inertia.h"

using namespace Delta2;

MeshData::MeshData(std::string path, common::Options& opt, bool isSurface) {
    bool result = igl::readOBJ(path, _vertices, _faces);
    _surrogate = model::SurrogateTree(_vertices.template cast<double>(), _faces, opt);
    initBoudingData();
    if (isSurface) {
        initInertiaDataSurface(0.01);
    } else {
        initInertiaDataVolume();
    }
}

MeshData::MeshData(Eigen::MatrixXd V, Eigen::MatrixXi F, common::Options& opt, bool isSurface) {
    _surrogate = model::SurrogateTree(V.template cast<double>(), F, opt);
    _vertices = V;
    _faces = F;
    initBoudingData();
    if (isSurface) {
        initInertiaDataSurface(0.01);
    } else {
        initInertiaDataVolume();
    }
}

Eigen::Matrix3d MeshData::getUnitInertiaTensor() {
    return _unit_inertia_tensor;
}

Eigen::Matrix3d MeshData::getUnitInertiaTensorInverse() {
    return _unit_inertia_tensor_inverse;
}

double MeshData::getVolume() {
    return _volume;
}

void MeshData::initBoudingData() {
    _AABB = embree::BBox3fa(embree::EmptyTy());
    _boundingRadius = 0;
    for (int v = 0; v < _vertices.rows(); v++) {
        _AABB.extend({float(_vertices(v, 0)), float(_vertices(v, 1)), float(_vertices(v, 2))});
        _boundingRadius = std::max(_boundingRadius, Eigen::Vector3d({_vertices(v, 0), _vertices(v, 1), _vertices(v, 2)}).norm());
    }
}

void MeshData::initInertiaDataVolume() {
    Eigen::MatrixXd tet_vertices; // Tetrahedralized vertices  |TV| (or |V| is this just a copy of V?) by 3 vertex positions
    Eigen::MatrixXi tet_tets; // Tetrahedralized tet indices |TT| by 4 vertex positions
    Eigen::MatrixXi tet_faces; // Tetrahedralized faces |TF| (or |F| is this just a copy of F?) by 3 face indices

    //std::string switches = "pq1.414a0.01";
    std::string switches = "Q";

    igl::copyleft::tetgen::tetrahedralize(_vertices, _faces, switches, tet_vertices, tet_tets, tet_faces);
    Eigen::MatrixXd tet_volumes;
    igl::volume(tet_vertices, tet_tets, tet_volumes);
    _volume = tet_volumes.sum();

    Eigen::Vector3d centre_of_mass = {0.0, 0.0, 0.0};

    for (int t = 0; t < tet_tets.rows(); t++) {
        Eigen::Vector3d a({tet_vertices(tet_tets(t, 0), 0), tet_vertices(tet_tets(t, 0), 1), tet_vertices(tet_tets(t, 0), 2)});
        Eigen::Vector3d b({tet_vertices(tet_tets(t, 1), 0), tet_vertices(tet_tets(t, 1), 1), tet_vertices(tet_tets(t, 1), 2)});
        Eigen::Vector3d c({tet_vertices(tet_tets(t, 2), 0), tet_vertices(tet_tets(t, 2), 1), tet_vertices(tet_tets(t, 2), 2)});
        Eigen::Vector3d d({tet_vertices(tet_tets(t, 3), 0), tet_vertices(tet_tets(t, 3), 1), tet_vertices(tet_tets(t, 3), 2)});
        centre_of_mass += tet_volumes(t) * (a + b + c + d) / 4.0;
    }
    centre_of_mass /= _volume;

    _vertices.rowwise() -= centre_of_mass.transpose();
    tet_vertices.rowwise() -= centre_of_mass.transpose(); // TODO do we really need this or just use _vertices?

    _unit_inertia_tensor.setZero();
    for (int t = 0; t < tet_tets.rows(); t++) {
        Eigen::Vector3d a({tet_vertices(tet_tets(t, 0), 0), tet_vertices(tet_tets(t, 0), 1), tet_vertices(tet_tets(t, 0), 2)});
        Eigen::Vector3d b({tet_vertices(tet_tets(t, 1), 0), tet_vertices(tet_tets(t, 1), 1), tet_vertices(tet_tets(t, 1), 2)});
        Eigen::Vector3d c({tet_vertices(tet_tets(t, 2), 0), tet_vertices(tet_tets(t, 2), 1), tet_vertices(tet_tets(t, 2), 2)});
        Eigen::Vector3d d({tet_vertices(tet_tets(t, 3), 0), tet_vertices(tet_tets(t, 3), 1), tet_vertices(tet_tets(t, 3), 2)});

        _unit_inertia_tensor += model::tetrahedronInertiaTensor(a, b, c, d);
    }

    _unit_inertia_tensor_inverse = _unit_inertia_tensor.inverse();
}

void MeshData::initInertiaDataSurface(double thickness) {
    // https://stackoverflow.com/questions/67078659/how-can-i-calculate-the-inertia-tensor-of-a-hollow-object-defined-by-a-triangle
    _volume = 0;
    Eigen::Vector3d centre_of_mass = {0.0, 0.0, 0.0};
    for (int f = 0; f < _faces.rows(); f++) {
        int i_a = _faces(f, 0);
        Eigen::Vector3d a({_vertices(i_a, 0), _vertices(i_a, 1), _vertices(i_a, 2)});
        int i_b = _faces(f, 1);
        Eigen::Vector3d b({_vertices(i_b, 0), _vertices(i_b, 1), _vertices(i_b, 2)});
        int i_c = _faces(f, 2);
        Eigen::Vector3d c({_vertices(i_c, 0), _vertices(i_c, 1), _vertices(i_c, 2)});

        double t_vol = Delta2::common::Triangle<double>(a, b, c).area() * thickness;
        Eigen::Vector3d avg = (a + b + c) / 3.0;
        _volume += t_vol;
        centre_of_mass += avg * t_vol;
    }

    _vertices.rowwise() -= centre_of_mass.transpose() / _volume;

    _unit_inertia_tensor.setZero();
    for (int f = 0; f < _faces.rows(); f++) {
        int i_a = _faces(f, 0);
        Eigen::Vector3d a({_vertices(i_a, 0), _vertices(i_a, 1), _vertices(i_a, 2)});
        int i_b = _faces(f, 1);
        Eigen::Vector3d b({_vertices(i_b, 0), _vertices(i_b, 1), _vertices(i_b, 2)});
        int i_c = _faces(f, 2);
        Eigen::Vector3d c({_vertices(i_c, 0), _vertices(i_c, 1), _vertices(i_c, 2)});

        Eigen::Matrix3d a_inertia = Delta2::model::pointInertiaTensor(a);
        Eigen::Matrix3d b_inertia = Delta2::model::pointInertiaTensor(b);
        Eigen::Matrix3d c_inertia = Delta2::model::pointInertiaTensor(c);

        double t_vol = Delta2::common::Triangle<double>(a, b, c).area() * thickness;
        _unit_inertia_tensor += t_vol * (a_inertia + b_inertia + c_inertia) / 3.0;
    }
    _unit_inertia_tensor_inverse = _unit_inertia_tensor.inverse();
}

embree::BBox3fa MeshData::getAABB() const {
    return _AABB;
}

embree::BBox3fa MeshData::getAABB(const Eigen::Matrix4d& transform) const {
    Eigen::Matrix4f T = transform.cast<float>();
    embree::BBox3fa result = embree::BBox3fa(embree::EmptyTy());
    for (float x_i = 0; x_i < 2; x_i++) {
        float x = x_i == 0 ? _AABB.lower.x : _AABB.upper.x;
        for (float y_i = 0; y_i < 2; y_i++) {
            float y = y_i == 0 ? _AABB.lower.y : _AABB.upper.y;
            for (float z_i = 0; z_i < 2; z_i++) {
                float z = z_i == 0 ? _AABB.lower.z : _AABB.upper.z;

                Eigen::Vector3f v = common::transform(Eigen::Vector3f({x, y, z}), T);

                result.extend({v.x(), v.y(), v.z()});
            }
        }
    }
    return result;
}

embree::BBox3fa MeshData::getSphereBounds() {
    float br = _boundingRadius;
    return embree::BBox3fa({-br, -br, -br}, {br, br, br});
}

double MeshData::getBoundingRadius() {
    return _boundingRadius;
}

Delta2::model::SurrogateTree& MeshData::getSurrogateTree() {
    return _surrogate;
}

const Eigen::MatrixXd& MeshData::getVertices() {
    return _vertices;
}

const Eigen::MatrixXi& MeshData::getFaces() {
    return _faces;
}

const std::string MeshData::serialise() {
    std::stringstream result;
    result << "V\n" << getVertices() << "\n";
    result << "F\n" << getFaces() << "\n";
    return result.str();
}
