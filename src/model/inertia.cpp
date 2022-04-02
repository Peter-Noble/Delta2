#include "inertia.h"

Eigen::Matrix3d star(const Eigen::Vector3d& a) {
	Eigen::Matrix3d result;
	result(0, 0) = 0;
	result(0, 1) = -a[2];
	result(0, 2) = a[1];

	result(1, 0) = a[2];
	result(1, 1) = 0;
	result(1, 2) = -a[0];

	result(2, 0) = -a[1];
	result(2, 1) = a[0];
	result(2, 2) = 0;

	return result;
}

Eigen::Matrix3d Delta2::model::tetrahedronInertiaTensor(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d) {
    // From: Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates
	//	https://thescipub.com/pdf/10.3844/jmssp.2005.8.11
    double density = -1.0; // Use a unit density here to begin with and then we'll multiply later

	double detJ = (a - d).dot((b - d).cross(c - d));
	double ixx = density * detJ * (
		a[1] * a[1] +
		a[1] * b[1] + b[1] * b[1] + 
		a[1] * c[1] + b[1] * c[1] + c[1] * c[1] +
		a[1] * d[1] + b[1] * d[1] + c[1] * d[1] + d[1] * d[1] +
		a[2] * a[2] +
		a[2] * b[2] + b[2] * b[2] +
		a[2] * c[2] + b[2] * c[2] + c[2] * c[2] +
		a[2] * d[2] + b[2] * d[2] + c[2] * d[2] + d[2] * d[2]) / 60;
	double iyy = density * detJ * (
		a[0] * a[0] +
		a[0] * b[0] + b[0] * b[0] +
		a[0] * c[0] + b[0] * c[0] + c[0] * c[0] +
		a[0] * d[0] + b[0] * d[0] + c[0] * d[0] + d[0] * d[0] +
		a[2] * a[2] +
		a[2] * b[2] + b[2] * b[2] +
		a[2] * c[2] + b[2] * c[2] + c[2] * c[2] +
		a[2] * d[2] + b[2] * d[2] + c[2] * d[2] + d[2] * d[2]) / 60;
	double izz = density * detJ * (
		a[0] * a[0] +
		a[0] * b[0] + b[0] * b[0] +
		a[0] * c[0] + b[0] * c[0] + c[0] * c[0] +
		a[0] * d[0] + b[0] * d[0] + c[0] * d[0] + d[0] * d[0] +
		a[1] * a[1] +
		a[1] * b[1] + b[1] * b[1] +
		a[1] * c[1] + b[1] * c[1] + c[1] * c[1] +
		a[1] * d[1] + b[1] * d[1] + c[1] * d[1] + d[1] * d[1]) / 60;
	double iyz = -density * detJ * (
		2 * a[1] * a[2] + b[1] * a[2] + c[1] * a[2] + d[1] * a[2] +
		a[1] * b[2] + 2 * b[1] * b[2] + c[1] * b[2] + d[1] * b[2] +
		a[1] * c[2] + b[1] * c[2] + 2 * c[1] * c[2] + d[1] * c[2] +
		a[1] * d[2] + b[1] * d[2] + c[1] * d[2] + 2 * d[1] * d[2]) / 120;
	double ixy = -density * detJ * (
		2 * a[0] * a[2] + b[0] * a[2] + c[0] * a[2] + d[0] * a[2] +
		a[0] * b[2] + 2 * b[0] * b[2] + c[0] * b[2] + d[0] * b[2] +
		a[0] * c[2] + b[0] * c[2] + 2 * c[0] * c[2] + d[0] * c[2] +
		a[0] * d[2] + b[0] * d[2] + c[0] * d[2] + 2 * d[0] * d[2]) / 120;
	double ixz = -density * detJ * (
		2 * a[0] * a[1] + b[0] * a[1] + c[0] * a[1] + d[0] * a[1] +
		a[0] * b[1] + 2 * b[0] * b[1] + c[0] * b[1] + d[0] * b[1] +
		a[0] * c[1] + b[0] * c[1] + 2 * c[0] * c[1] + d[0] * c[1] +
		a[0] * d[1] + b[0] * d[1] + c[0] * d[1] + 2 * d[0] * d[1]) / 120;
	Eigen::Matrix3d result;
	result << ixx, ixy, ixz,
		      ixy, iyy, iyz,
		      ixz, iyz, izz;
	return result;
}

Eigen::Matrix3d Delta2::model::pointInertiaTensor(const Eigen::Vector3d& pt) {
	// https://stackoverflow.com/questions/67078659/how-can-i-calculate-the-inertia-tensor-of-a-hollow-object-defined-by-a-triangle
	double x = pt.x();
	double y = pt.y();
	double z = pt.z();
	return Eigen::Matrix3d({
        {y*y + z*z,      -x*y,      -x*z},
        {     -x*y, x*x + z*z,      -y*z},
        {     -x*z,      -y*z, x*x + y*y}
	});
}
