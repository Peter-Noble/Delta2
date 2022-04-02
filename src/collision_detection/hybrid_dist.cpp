// #include "hybrid_dist.h"
// #include "penalty_dist.h"

// void penalty_sympy(triangle<float>& A, std::vector<triangle<float>>& Bs, std::vector<Eigen::Vector3d>& P, std::vector<Eigen::Vector3d>& Q, std::vector<float>& max_error, std::vector<bool>& failed_vec) {
// 	const int Iterations = 4;
// 	const int ProjectCorrectIterations = 4;
// 	const int TrianglePairs = 8;

// 	const float x_a11 = A.xs[0];
// 	const float x_a12 = A.ys[0];
// 	const float x_a13 = A.zs[0];
// 	const float x_a21 = A.xs[1];
// 	const float x_a22 = A.ys[1];
// 	const float x_a23 = A.zs[1];
// 	const float x_a31 = A.xs[2];
// 	const float x_a32 = A.ys[2];
// 	const float x_a33 = A.zs[2];

// 	ALIGNED(float x_b11[TrianglePairs], 64);
// 	ALIGNED(float x_b12[TrianglePairs], 64);
// 	ALIGNED(float x_b13[TrianglePairs], 64);
// 	ALIGNED(float x_b21[TrianglePairs], 64);
// 	ALIGNED(float x_b22[TrianglePairs], 64);
// 	ALIGNED(float x_b23[TrianglePairs], 64);
// 	ALIGNED(float x_b31[TrianglePairs], 64);
// 	ALIGNED(float x_b32[TrianglePairs], 64);
// 	ALIGNED(float x_b33[TrianglePairs], 64);

// 	for (int t = 0; t < TrianglePairs; t++) {
// 		x_b11[t] = Bs[t % Bs.size()].xs[0];
// 		x_b12[t] = Bs[t % Bs.size()].ys[0];
// 		x_b13[t] = Bs[t % Bs.size()].zs[0];
// 		x_b21[t] = Bs[t % Bs.size()].xs[1];
// 		x_b22[t] = Bs[t % Bs.size()].ys[1];
// 		x_b23[t] = Bs[t % Bs.size()].zs[1];
// 		x_b31[t] = Bs[t % Bs.size()].xs[2];
// 		x_b32[t] = Bs[t % Bs.size()].ys[2];
// 		x_b33[t] = Bs[t % Bs.size()].zs[2];
// 	}

// 	ALIGNED(float a[TrianglePairs], 64);
// 	ALIGNED(float b[TrianglePairs], 64);
// 	ALIGNED(float c[TrianglePairs], 64);
// 	ALIGNED(float d[TrianglePairs], 64);

// 	ALIGNED(float a_last[TrianglePairs], 64);
// 	ALIGNED(float b_last[TrianglePairs], 64);
// 	ALIGNED(float c_last[TrianglePairs], 64);
// 	ALIGNED(float d_last[TrianglePairs], 64);

// 	ALIGNED(float a_new[TrianglePairs], 64);
// 	ALIGNED(float b_new[TrianglePairs], 64);
// 	ALIGNED(float c_new[TrianglePairs], 64);
// 	ALIGNED(float d_new[TrianglePairs], 64);

// 	for (int i = 0; i < TrianglePairs; i++) {
// 		a[i] = 0.33;
// 		b[i] = 0.33;
// 		c[i] = 0.33;
// 		d[i] = 0.33;
// 	}

// 	ALIGNED(float failed[TrianglePairs], 64);

// 	float alpha = 1;
// 	float C_dontdividebyzero = 1e-5;
// 	float C_damp_start = 1.0;
// 	float C_damp_end = 0.5;

// 	for (int i = 0; i < Iterations; i++) {
// 		float C_damp = C_damp_start - (C_damp_start - C_damp_end) * i / (Iterations - 1);
// #pragma omp simd
// 		for (int t = 0; t < TrianglePairs; t++) {
// 			a_last[t] = a[t];
// 			b_last[t] = b[t];
// 			c_last[t] = c[t];
// 			d_last[t] = d[t];

// 			a_new[t] = a[t] - C_damp * ((x_a11 - x_a21) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) + (x_a12 - x_a22) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) + (x_a13 - x_a23) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + (x_a11 - x_a21) * (x_a11 - x_a21) + (x_a12 - x_a22) * (x_a12 - x_a22) + (x_a13 - x_a23) * (x_a13 - x_a23));
// 			b_new[t] = b[t] - C_damp * ((x_a11 - x_a31) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) + (x_a12 - x_a32) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) + (x_a13 - x_a33) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + (x_a11 - x_a31) * (x_a11 - x_a31) + (x_a12 - x_a32) * (x_a12 - x_a32) + (x_a13 - x_a33) * (x_a13 - x_a33));
// 			c_new[t] = c[t] - C_damp * (-1.0f * (x_b11[t] - x_b21[t]) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) - 1.0f * (x_b12[t] - x_b22[t]) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) - 1.0f * (x_b13[t] - x_b23[t]) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b21[t]) * (x_b11[t] - x_b21[t]) + 1.0f * (x_b12[t] - x_b22[t]) * (x_b12[t] - x_b22[t]) + 1.0f * (x_b13[t] - x_b23[t]) * (x_b13[t] - x_b23[t]));
// 			d_new[t] = d[t] - C_damp * (-1.0f * (x_b11[t] - x_b31[t]) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) - 1.0f * (x_b12[t] - x_b32[t]) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) - 1.0f * (x_b13[t] - x_b33[t]) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b31[t]) * (x_b11[t] - x_b31[t]) + 1.0f * (x_b12[t] - x_b32[t]) * (x_b12[t] - x_b32[t]) + 1.0f * (x_b13[t] - x_b33[t]) * (x_b13[t] - x_b33[t]));

// 			/*if (t == 0) {
// 				printf("pre\n");
// 				printf("a: %f, b: %f, c: %f, d: %f\n", a_new[t], b_new[t], c_new[t], d_new[t]);
// 			}*/

// 			for (int i = 0; i < ProjectCorrectIterations; i++) {
// 				a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 				b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 				c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
// 				d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));

// 				a_new[t] = a[t] + alpha * (-a[t] * Heaviside(-a[t]) - 0.5f * (a[t] + b[t] - 1) * Heaviside(a[t] + b[t] - 1));
// 				b_new[t] = b[t] + alpha * (-b[t] * Heaviside(-b[t]) - 0.5f * (a[t] + b[t] - 1) * Heaviside(a[t] + b[t] - 1));
// 				c_new[t] = c[t] + alpha * (-c[t] * Heaviside(-c[t]) - 0.5f * (c[t] + d[t] - 1) * Heaviside(c[t] + d[t] - 1));
// 				d_new[t] = d[t] + alpha * (-d[t] * Heaviside(-d[t]) - 0.5f * (c[t] + d[t] - 1) * Heaviside(c[t] + d[t] - 1));
// 			}
// 			a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 			b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 			c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
// 			d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));

// 			P[t] = { x_a11 + a[t] * (x_a21 - x_a11) + b[t] * (x_a31 - x_a11), x_a12 + a[t] * (x_a22 - x_a12) + b[t] * (x_a32 - x_a12), x_a13 + a[t] * (x_a23 - x_a13) + b[t] * (x_a33 - x_a13) };
// 			Q[t] = { x_b11[t] + c[t] * (x_b21[t] - x_b11[t]) + d[t] * (x_b31[t] - x_b11[t]), x_b12[t] + c[t] * (x_b22[t] - x_b12[t]) + d[t] * (x_b32[t] - x_b12[t]), x_b13[t] + c[t] * (x_b23[t] - x_b13[t]) + d[t] * (x_b33[t] - x_b13[t]) };
// 		}
// 	}

// 	for (int t = 0; t < Bs.size(); t++) {
// 		a[t] = clamp(a[t]);
// 		b[t] = clamp(b[t]);
// 		c[t] = clamp(c[t]);
// 		d[t] = clamp(d[t]);

// 		if (1.0f < a[t] + b[t]) {
// 			a[t] = a[t] / (a[t] + b[t]);
// 			b[t] = b[t] / (a[t] + b[t]);
// 		}

// 		if (1.0f < c[t] + d[t]) {
// 			c[t] = c[t] / (c[t] + d[t]);
// 			d[t] = d[t] / (c[t] + d[t]);
// 		}

// 		float Ax_last_update = (a[t] - a_last[t]) * (x_a21 - x_a11) + (b[t] - b_last[t]) * (x_a31 - x_a11);
// 		float Ay_last_update = (a[t] - a_last[t]) * (x_a22 - x_a12) + (b[t] - b_last[t]) * (x_a32 - x_a12);
// 		float Az_last_update = (a[t] - a_last[t]) * (x_a23 - x_a13) + (b[t] - b_last[t]) * (x_a33 - x_a13);

// 		float Bx_last_update = (c[t] - c_last[t]) * (x_b21[t] - x_b11[t]) + (d[t] - d_last[t]) * (x_b31[t] - x_b11[t]);
// 		float By_last_update = (c[t] - c_last[t]) * (x_b22[t] - x_b12[t]) + (d[t] - d_last[t]) * (x_b32[t] - x_b12[t]);
// 		float Bz_last_update = (c[t] - c_last[t]) * (x_b23[t] - x_b13[t]) + (d[t] - d_last[t]) * (x_b33[t] - x_b13[t]);

// 		float A_all_update = sqrt(Ax_last_update * Ax_last_update + Ay_last_update * Ay_last_update + Az_last_update * Az_last_update) * 1 / C_damp_end;
// 		float B_all_update = sqrt(Bx_last_update * Bx_last_update + By_last_update * By_last_update + Bz_last_update * Bz_last_update) * 1 / C_damp_end;
// 		//if (t == 0) {
// 		//	printf("update: %.4f, %.4f\n", A_all_update, B_all_update);
// 		//	printf("max error: %.4f\n", max_error[t]);
// 		//	//printf("a: %f, b: %f, c: %f, d: %f\n", a[t] - a_last[t], b[t] - b_last[t], c[t] - c_last[t], d[t] - d_last[t]);
// 		//	printf("a: %f, b: %f, c: %f, d: %f\n", a[t], b[t], c[t], d[t]);
// 		//}

// 		failed[t] = A_all_update > max_error[t] || B_all_update > max_error[t];
// 	}

// 	for (int t = 0; t < Bs.size(); t++) {
// 		//printf("a: %.4f, b: %.4f, c: %.4f, d: %.4f\n", a[t], b[t], c[t], d[t]);
// 		failed_vec[t] = failed[t];
// 		P[t] = { x_a11 + a[t] * (x_a21 - x_a11) + b[t] * (x_a31 - x_a11), x_a12 + a[t] * (x_a22 - x_a12) + b[t] * (x_a32 - x_a12), x_a13 + a[t] * (x_a23 - x_a13) + b[t] * (x_a33 - x_a13) };
// 		Q[t] = { x_b11[t] + c[t] * (x_b21[t] - x_b11[t]) + d[t] * (x_b31[t] - x_b11[t]), x_b12[t] + c[t] * (x_b22[t] - x_b12[t]) + d[t] * (x_b32[t] - x_b12[t]), x_b13[t] + c[t] * (x_b23[t] - x_b13[t]) + d[t] * (x_b33[t] - x_b13[t]) };
// 	};
// }


// void penalty_sympy_1_to_1_times_n(std::vector<triangle<float>>& As, std::vector<triangle<float>>& Bs, std::vector<Eigen::Vector3d>& P, std::vector<Eigen::Vector3d>& Q, std::vector<float>& max_error, std::vector<bool>& failed_vec) {
// 	const int Iterations = 6;
// 	const int TrianglePairs = 8;

// 	assert(As.size() == Bs.size());
// 	int size = As.size();

// 	ALIGNED(float x_a11[TrianglePairs], 64);
// 	ALIGNED(float x_a12[TrianglePairs], 64);
// 	ALIGNED(float x_a13[TrianglePairs], 64);
// 	ALIGNED(float x_a21[TrianglePairs], 64);
// 	ALIGNED(float x_a22[TrianglePairs], 64);
// 	ALIGNED(float x_a23[TrianglePairs], 64);
// 	ALIGNED(float x_a31[TrianglePairs], 64);
// 	ALIGNED(float x_a32[TrianglePairs], 64);
// 	ALIGNED(float x_a33[TrianglePairs], 64);

// 	ALIGNED(float x_b11[TrianglePairs], 64);
// 	ALIGNED(float x_b12[TrianglePairs], 64);
// 	ALIGNED(float x_b13[TrianglePairs], 64);
// 	ALIGNED(float x_b21[TrianglePairs], 64);
// 	ALIGNED(float x_b22[TrianglePairs], 64);
// 	ALIGNED(float x_b23[TrianglePairs], 64);
// 	ALIGNED(float x_b31[TrianglePairs], 64);
// 	ALIGNED(float x_b32[TrianglePairs], 64);
// 	ALIGNED(float x_b33[TrianglePairs], 64);

// 	ALIGNED(float a[TrianglePairs], 64);
// 	ALIGNED(float b[TrianglePairs], 64);
// 	ALIGNED(float c[TrianglePairs], 64);
// 	ALIGNED(float d[TrianglePairs], 64);

// 	ALIGNED(float a_last[TrianglePairs], 64);
// 	ALIGNED(float b_last[TrianglePairs], 64);
// 	ALIGNED(float c_last[TrianglePairs], 64);
// 	ALIGNED(float d_last[TrianglePairs], 64);

// 	ALIGNED(float a_new[TrianglePairs], 64);
// 	ALIGNED(float b_new[TrianglePairs], 64);
// 	ALIGNED(float c_new[TrianglePairs], 64);
// 	ALIGNED(float d_new[TrianglePairs], 64);

// 	for (int batch = 0; batch < (size + TrianglePairs - 1) / TrianglePairs; batch++) {
// 		for (int t = 0; t < TrianglePairs; t++) {
// 			int ind = std::min(TrianglePairs * batch + t, size - 1);
// 			x_a11[t] = As[ind].xs[0];
// 			x_a12[t] = As[ind].ys[0];
// 			x_a13[t] = As[ind].zs[0];
// 			x_a21[t] = As[ind].xs[1];
// 			x_a22[t] = As[ind].ys[1];
// 			x_a23[t] = As[ind].zs[1];
// 			x_a31[t] = As[ind].xs[2];
// 			x_a32[t] = As[ind].ys[2];
// 			x_a33[t] = As[ind].zs[2];

// 			x_b11[t] = Bs[ind].xs[0];
// 			x_b12[t] = Bs[ind].ys[0];
// 			x_b13[t] = Bs[ind].zs[0];
// 			x_b21[t] = Bs[ind].xs[1];
// 			x_b22[t] = Bs[ind].ys[1];
// 			x_b23[t] = Bs[ind].zs[1];
// 			x_b31[t] = Bs[ind].xs[2];
// 			x_b32[t] = Bs[ind].ys[2];
// 			x_b33[t] = Bs[ind].zs[2];
// 		}

// 		for (int i = 0; i < TrianglePairs; i++) {
// 			a[i] = 0.33;
// 			b[i] = 0.33;
// 			c[i] = 0.33;
// 			d[i] = 0.33;
// 		}

// 		float alpha = 1;
// 		float C_dontdividebyzero = 1e-5;
// 		float C_damp_start = 1.0;
// 		float C_damp_end = 1.0;

// 		for (int i = 0; i < Iterations; i++) {
// 			float C_damp = C_damp_start - (C_damp_start - C_damp_end) * i / (Iterations - 1);
// #pragma omp simd
// 			for (int t = 0; t < TrianglePairs; t++) {
// 				a_last[t] = a[t];
// 				b_last[t] = b[t];
// 				c_last[t] = c[t];
// 				d_last[t] = d[t];

// 				a_new[t] = a[t] - 1.0f * ((x_a11[t] - x_a21[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) + (x_a12[t] - x_a22[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) + (x_a13[t] - x_a23[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + (x_a11[t] - x_a21[t]) * (x_a11[t] - x_a21[t]) + (x_a12[t] - x_a22[t]) * (x_a12[t] - x_a22[t]) + (x_a13[t] - x_a23[t]) * (x_a13[t] - x_a23[t]));
// 				b_new[t] = b[t] - 1.0f * ((x_a11[t] - x_a31[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) + (x_a12[t] - x_a32[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) + (x_a13[t] - x_a33[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + (x_a11[t] - x_a31[t]) * (x_a11[t] - x_a31[t]) + (x_a12[t] - x_a32[t]) * (x_a12[t] - x_a32[t]) + (x_a13[t] - x_a33[t]) * (x_a13[t] - x_a33[t]));
// 				c_new[t] = c[t] - 1.0f * (-1.0f * (x_b11[t] - x_b21[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) - 1.0f * (x_b12[t] - x_b22[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) - 1.0f * (x_b13[t] - x_b23[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b21[t]) * (x_b11[t] - x_b21[t]) + 1.0f * (x_b12[t] - x_b22[t]) * (x_b12[t] - x_b22[t]) + 1.0f * (x_b13[t] - x_b23[t]) * (x_b13[t] - x_b23[t]));
// 				d_new[t] = d[t] - 1.0f * (-1.0f * (x_b11[t] - x_b31[t]) * (a[t] * (x_a11[t] - x_a21[t]) + b[t] * (x_a11[t] - x_a31[t]) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11[t] + x_b11[t]) - 1.0f * (x_b12[t] - x_b32[t]) * (a[t] * (x_a12[t] - x_a22[t]) + b[t] * (x_a12[t] - x_a32[t]) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12[t] + x_b12[t]) - 1.0f * (x_b13[t] - x_b33[t]) * (a[t] * (x_a13[t] - x_a23[t]) + b[t] * (x_a13[t] - x_a33[t]) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13[t] + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b31[t]) * (x_b11[t] - x_b31[t]) + 1.0f * (x_b12[t] - x_b32[t]) * (x_b12[t] - x_b32[t]) + 1.0f * (x_b13[t] - x_b33[t]) * (x_b13[t] - x_b33[t]));

// 				for (int i = 0; i < 4; i++) {
// 					a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
// 					b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
// 					c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));
// 					d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));

// 					a_new[t] = a[t] + alpha * (-a[t] * Heaviside(-a[t]) - 0.5f * (a[t] + b[t] - 1.0f) * Heaviside(a[t] + b[t] - 1.0f));
// 					b_new[t] = b[t] + alpha * (-b[t] * Heaviside(-b[t]) - 0.5f * (a[t] + b[t] - 1.0f) * Heaviside(a[t] + b[t] - 1.0f));
// 					c_new[t] = c[t] + alpha * (-c[t] * Heaviside(-c[t]) - 0.5f * (c[t] + d[t] - 1.0f) * Heaviside(c[t] + d[t] - 1.0f));
// 					d_new[t] = d[t] + alpha * (-d[t] * Heaviside(-d[t]) - 0.5f * (c[t] + d[t] - 1.0f) * Heaviside(c[t] + d[t] - 1.0f));
// 				}

// 				a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
// 				b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1.0f) * Heaviside(a_new[t] + b_new[t] - 1.0f));
// 				c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));
// 				d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1.0f) * Heaviside(c_new[t] + d_new[t] - 1.0f));
// 			}
// 		}

// 		for (int t = 0; t < TrianglePairs && TrianglePairs * batch + t < size; t++) {
// 			int ind = TrianglePairs * batch + t;

// 			a[t] = clamp(a[t]);
// 			b[t] = clamp(b[t]);
// 			c[t] = clamp(c[t]);
// 			d[t] = clamp(d[t]);

// 			if (1.0f < a[t] + b[t]) {
// 				a_new[t] = a[t] / (a[t] + b[t]);
// 				b_new[t] = b[t] / (a[t] + b[t]);
// 				a[t] = a_new[t];
// 				b[t] = b_new[t];
// 			}

// 			if (1.0f < c[t] + d[t]) {
// 				c_new[t] = c[t] / (c[t] + d[t]);
// 				d_new[t] = d[t] / (c[t] + d[t]);
// 				c[t] = c_new[t];
// 				d[t] = d_new[t];
// 			}

// 			float Ax_last_update = (a[t] - a_last[t]) * (x_a21[t] - x_a11[t]) + (b[t] - b_last[t]) * (x_a31[t] - x_a11[t]);
// 			float Ay_last_update = (a[t] - a_last[t]) * (x_a22[t] - x_a12[t]) + (b[t] - b_last[t]) * (x_a32[t] - x_a12[t]);
// 			float Az_last_update = (a[t] - a_last[t]) * (x_a23[t] - x_a13[t]) + (b[t] - b_last[t]) * (x_a33[t] - x_a13[t]);

// 			float Bx_last_update = (c[t] - c_last[t]) * (x_b21[t] - x_b11[t]) + (d[t] - d_last[t]) * (x_b31[t] - x_b11[t]);
// 			float By_last_update = (c[t] - c_last[t]) * (x_b22[t] - x_b12[t]) + (d[t] - d_last[t]) * (x_b32[t] - x_b12[t]);
// 			float Bz_last_update = (c[t] - c_last[t]) * (x_b23[t] - x_b13[t]) + (d[t] - d_last[t]) * (x_b33[t] - x_b13[t]);

// 			float A_all_update = sqrt(Ax_last_update * Ax_last_update + Ay_last_update * Ay_last_update + Az_last_update * Az_last_update) * 1 / C_damp_end;
// 			float B_all_update = sqrt(Bx_last_update * Bx_last_update + By_last_update * By_last_update + Bz_last_update * Bz_last_update) * 1 / C_damp_end;

// 			failed_vec[ind] = A_all_update > max_error[ind] || B_all_update > max_error[ind];
// 			P[ind] = { x_a11[t] + a[t] * (x_a21[t] - x_a11[t]) + b[t] * (x_a31[t] - x_a11[t]), x_a12[t] + a[t] * (x_a22[t] - x_a12[t]) + b[t] * (x_a32[t] - x_a12[t]), x_a13[t] + a[t] * (x_a23[t] - x_a13[t]) + b[t] * (x_a33[t] - x_a13[t]) };
// 			Q[ind] = { x_b11[t] + c[t] * (x_b21[t] - x_b11[t]) + d[t] * (x_b31[t] - x_b11[t]), x_b12[t] + c[t] * (x_b22[t] - x_b12[t]) + d[t] * (x_b32[t] - x_b12[t]), x_b13[t] + c[t] * (x_b23[t] - x_b13[t]) + d[t] * (x_b33[t] - x_b13[t]) };
// 		}
// 	}
// }

// void sympy_bucket(
// 	triangle<float>& a_parent_tri,
// 	float a_parent_eps,
// 	triangle<float>& b_parent_tri,
// 	float b_parent_eps,
// 	std::vector<int>& a_children,
// 	std::vector<int>& b_children,
// 	const std::vector<triangle<double>>& a_tris,
// 	const std::vector<triangle<double>>& b_tris,
// 	Eigen::Matrix4d a_trans,
// 	Eigen::Matrix4d b_trans,
// 	double a_eps,
// 	double b_eps,
// 	std::vector<contact>& result,
// 	int& ttd_total,
// 	int& ttd_penatly_total)
// {
// 	float search_dist = a_eps + b_eps;

// 	const int Iterations = 4;
// 	const int ProjectCorrectIterations = 8;
// 	const int TrianglePairs = 8;

// 	ALIGNED(float x_b11[TrianglePairs], 64);
// 	ALIGNED(float x_b12[TrianglePairs], 64);
// 	ALIGNED(float x_b13[TrianglePairs], 64);
// 	ALIGNED(float x_b21[TrianglePairs], 64);
// 	ALIGNED(float x_b22[TrianglePairs], 64);
// 	ALIGNED(float x_b23[TrianglePairs], 64);
// 	ALIGNED(float x_b31[TrianglePairs], 64);
// 	ALIGNED(float x_b32[TrianglePairs], 64);
// 	ALIGNED(float x_b33[TrianglePairs], 64);

// 	ALIGNED(float a[TrianglePairs], 64);
// 	ALIGNED(float b[TrianglePairs], 64);
// 	ALIGNED(float c[TrianglePairs], 64);
// 	ALIGNED(float d[TrianglePairs], 64);

// 	ALIGNED(float a_last[TrianglePairs], 64);
// 	ALIGNED(float b_last[TrianglePairs], 64);
// 	ALIGNED(float c_last[TrianglePairs], 64);
// 	ALIGNED(float d_last[TrianglePairs], 64);

// 	ALIGNED(float a_new[TrianglePairs], 64);
// 	ALIGNED(float b_new[TrianglePairs], 64);
// 	ALIGNED(float c_new[TrianglePairs], 64);
// 	ALIGNED(float d_new[TrianglePairs], 64);

// 	ALIGNED(float failed[TrianglePairs], 64);

// 	std::vector<triangle<float>> As_failed, Bs_failed;

// 	for (int b_grp = 0; b_grp < (b_children.size() + TrianglePairs - 1) / TrianglePairs; b_grp++) {
// 		int tris_used = 0;
// 		for (int t = 0; t < TrianglePairs; t++) {
// 			int ind = std::min(b_grp * TrianglePairs + t, (int) b_children.size() - 1);
// 			tris_used = (ind) % TrianglePairs + 1;
// 			triangle<float> B = transform(b_tris[b_children[ind]], b_trans).cast<float>();
// 			x_b11[t] = B.xs[0];
// 			x_b12[t] = B.ys[0];
// 			x_b13[t] = B.zs[0];
// 			x_b21[t] = B.xs[1];
// 			x_b22[t] = B.ys[1];
// 			x_b23[t] = B.zs[1];
// 			x_b31[t] = B.xs[2];
// 			x_b32[t] = B.ys[2];
// 			x_b33[t] = B.zs[2];
// 		}

// 		for (int a_child : a_children) {
// 			triangle<float> A = transform(a_tris[a_child], a_trans).cast<float>();

// 			const float x_a11 = A.xs[0];
// 			const float x_a12 = A.ys[0];
// 			const float x_a13 = A.zs[0];
// 			const float x_a21 = A.xs[1];
// 			const float x_a22 = A.ys[1];
// 			const float x_a23 = A.zs[1];
// 			const float x_a31 = A.xs[2];
// 			const float x_a32 = A.ys[2];
// 			const float x_a33 = A.zs[2];

// 			for (int i = 0; i < TrianglePairs; i++) {
// 				a[i] = 0.33;
// 				b[i] = 0.33;
// 				c[i] = 0.33;
// 				d[i] = 0.33;
// 			}

// 			float alpha = 1;
// 			float C_dontdividebyzero = 1e-5;
// 			float C_damp_start = 1.0;
// 			float C_damp_end = 1.0;

// 			//ttd_penatly_total += TrianglePairs;
// 			ttd_penatly_total += tris_used;

// 			for (int i = 0; i < Iterations; i++) {
// 				float C_damp = C_damp_start - (C_damp_start - C_damp_end) * i / (Iterations - 1);
// #pragma omp simd
// 				for (int t = 0; t < TrianglePairs; t++) {
// 					a_last[t] = a[t];
// 					b_last[t] = b[t];
// 					c_last[t] = c[t];
// 					d_last[t] = d[t];

// 					a_new[t] = a[t] - C_damp * ((x_a11 - x_a21) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) + (x_a12 - x_a22) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) + (x_a13 - x_a23) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + (x_a11 - x_a21) * (x_a11 - x_a21) + (x_a12 - x_a22) * (x_a12 - x_a22) + (x_a13 - x_a23) * (x_a13 - x_a23));
// 					b_new[t] = b[t] - C_damp * ((x_a11 - x_a31) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) + (x_a12 - x_a32) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) + (x_a13 - x_a33) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + (x_a11 - x_a31) * (x_a11 - x_a31) + (x_a12 - x_a32) * (x_a12 - x_a32) + (x_a13 - x_a33) * (x_a13 - x_a33));
// 					c_new[t] = c[t] - C_damp * (-1.0f * (x_b11[t] - x_b21[t]) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) - 1.0f * (x_b12[t] - x_b22[t]) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) - 1.0f * (x_b13[t] - x_b23[t]) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b21[t]) * (x_b11[t] - x_b21[t]) + 1.0f * (x_b12[t] - x_b22[t]) * (x_b12[t] - x_b22[t]) + 1.0f * (x_b13[t] - x_b23[t]) * (x_b13[t] - x_b23[t]));
// 					d_new[t] = d[t] - C_damp * (-1.0f * (x_b11[t] - x_b31[t]) * (a[t] * (x_a11 - x_a21) + b[t] * (x_a11 - x_a31) - c[t] * (x_b11[t] - x_b21[t]) - d[t] * (x_b11[t] - x_b31[t]) - x_a11 + x_b11[t]) - 1.0f * (x_b12[t] - x_b32[t]) * (a[t] * (x_a12 - x_a22) + b[t] * (x_a12 - x_a32) - c[t] * (x_b12[t] - x_b22[t]) - d[t] * (x_b12[t] - x_b32[t]) - x_a12 + x_b12[t]) - 1.0f * (x_b13[t] - x_b33[t]) * (a[t] * (x_a13 - x_a23) + b[t] * (x_a13 - x_a33) - c[t] * (x_b13[t] - x_b23[t]) - d[t] * (x_b13[t] - x_b33[t]) - x_a13 + x_b13[t])) / (C_dontdividebyzero + 1.0f * (x_b11[t] - x_b31[t]) * (x_b11[t] - x_b31[t]) + 1.0f * (x_b12[t] - x_b32[t]) * (x_b12[t] - x_b32[t]) + 1.0f * (x_b13[t] - x_b33[t]) * (x_b13[t] - x_b33[t]));

// 					for (int i = 0; i < ProjectCorrectIterations; i++) {
// 						a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 						b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 						c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
// 						d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));

// 						a_new[t] = a[t] + alpha * (-a[t] * Heaviside(-a[t]) - 0.5f * (a[t] + b[t] - 1) * Heaviside(a[t] + b[t] - 1));
// 						b_new[t] = b[t] + alpha * (-b[t] * Heaviside(-b[t]) - 0.5f * (a[t] + b[t] - 1) * Heaviside(a[t] + b[t] - 1));
// 						c_new[t] = c[t] + alpha * (-c[t] * Heaviside(-c[t]) - 0.5f * (c[t] + d[t] - 1) * Heaviside(c[t] + d[t] - 1));
// 						d_new[t] = d[t] + alpha * (-d[t] * Heaviside(-d[t]) - 0.5f * (c[t] + d[t] - 1) * Heaviside(c[t] + d[t] - 1));
// 					}
// 					a[t] = a_new[t] + alpha * (-a_new[t] * Heaviside(-a_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 					b[t] = b_new[t] + alpha * (-b_new[t] * Heaviside(-b_new[t]) - 0.5f * (a_new[t] + b_new[t] - 1) * Heaviside(a_new[t] + b_new[t] - 1));
// 					c[t] = c_new[t] + alpha * (-c_new[t] * Heaviside(-c_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
// 					d[t] = d_new[t] + alpha * (-d_new[t] * Heaviside(-d_new[t]) - 0.5f * (c_new[t] + d_new[t] - 1) * Heaviside(c_new[t] + d_new[t] - 1));
// 				}
// 			}

// 			for (int t = 0; t < TrianglePairs; t++) {
// 				a[t] = clamp(a[t]);
// 				b[t] = clamp(b[t]);
// 				c[t] = clamp(c[t]);
// 				d[t] = clamp(d[t]);

// 				if (1.0f < a[t] + b[t]) {
// 					a[t] = a[t] / (a[t] + b[t]);
// 					b[t] = b[t] / (a[t] + b[t]);
// 				}

// 				if (1.0f < c[t] + d[t]) {
// 					c[t] = c[t] / (c[t] + d[t]);
// 					d[t] = d[t] / (c[t] + d[t]);
// 				}

// 				float Ax_last_update = (a[t] - a_last[t]) * (x_a21 - x_a11) + (b[t] - b_last[t]) * (x_a31 - x_a11);
// 				float Ay_last_update = (a[t] - a_last[t]) * (x_a22 - x_a12) + (b[t] - b_last[t]) * (x_a32 - x_a12);
// 				float Az_last_update = (a[t] - a_last[t]) * (x_a23 - x_a13) + (b[t] - b_last[t]) * (x_a33 - x_a13);

// 				float Bx_last_update = (c[t] - c_last[t]) * (x_b21[t] - x_b11[t]) + (d[t] - d_last[t]) * (x_b31[t] - x_b11[t]);
// 				float By_last_update = (c[t] - c_last[t]) * (x_b22[t] - x_b12[t]) + (d[t] - d_last[t]) * (x_b32[t] - x_b12[t]);
// 				float Bz_last_update = (c[t] - c_last[t]) * (x_b23[t] - x_b13[t]) + (d[t] - d_last[t]) * (x_b33[t] - x_b13[t]);

// 				float A_all_update = sqrt(Ax_last_update * Ax_last_update + Ay_last_update * Ay_last_update + Az_last_update * Az_last_update) * 1 / C_damp_end;
// 				float B_all_update = sqrt(Bx_last_update * Bx_last_update + By_last_update * By_last_update + Bz_last_update * Bz_last_update) * 1 / C_damp_end;

// 				failed[t] = A_all_update > search_dist || B_all_update > search_dist;
// 			}

// 			for (int t = 0; t < TrianglePairs && b_grp * TrianglePairs + t < b_children.size(); t++) {
// 				if (failed[t]) {
// 					As_failed.push_back(A);
// 					int ind = b_grp * TrianglePairs + t;
// 					triangle<float> B = transform(b_tris[b_children[ind]], b_trans).cast<float>();
// 					Bs_failed.push_back(B);
// 				}
// 				else {
// 					Eigen::Vector3d P = { x_a11 + a[t] * (x_a21 - x_a11) + b[t] * (x_a31 - x_a11), x_a12 + a[t] * (x_a22 - x_a12) + b[t] * (x_a32 - x_a12), x_a13 + a[t] * (x_a23 - x_a13) + b[t] * (x_a33 - x_a13) };
// 					Eigen::Vector3d Q = { x_b11[t] + c[t] * (x_b21[t] - x_b11[t]) + d[t] * (x_b31[t] - x_b11[t]), x_b12[t] + c[t] * (x_b22[t] - x_b12[t]) + d[t] * (x_b32[t] - x_b12[t]), x_b13[t] + c[t] * (x_b23[t] - x_b13[t]) + d[t] * (x_b33[t] - x_b13[t]) };
// 					if ((P - Q).norm() < a_eps + b_eps) {
// 						result.push_back(std::make_tuple(P, Q, a_eps, b_eps));
// 					}
// 				}
// 			}
// 		}
// 	}

// 	std::vector<std::tuple<float, Eigen::Vector3d, Eigen::Vector3d>> fast_result = fast_1_to_1_times_n(As_failed, Bs_failed);
// 	ttd_total += As_failed.size();

// 	for (auto& r : fast_result) {
// 		if (std::sqrt(std::get<0>(r)) < a_eps + b_eps) {
// 			result.push_back(std::make_tuple(std::get<1>(r), std::get<2>(r), a_eps, b_eps));
// 		}
// 	}
// }



// void sympy_bucket_deferred(
// 	std::vector<deferred_compare>& buckets,
// 	const std::vector<triangle<double>>& a_tris,
// 	const std::vector<triangle<double>>& b_tris,
// 	Eigen::Matrix4d a_trans,
// 	Eigen::Matrix4d b_trans,
// 	float a_eps,
// 	float b_eps,
// 	std::vector<contactIO>& hits,
// 	int& ttd_total,
// 	int& ttd_penatly_total)
// {
// 	if (buckets.size() == 0) {
// 		return;
// 	}

// 	float search_dist = a_eps + b_eps;

// 	int batch_size = 256;

// 	std::vector<Eigen::Vector3d> Ps, Qs;
// 	Ps.resize(batch_size);
// 	Qs.resize(batch_size);
// 	std::vector<bool> failed;
// 	failed.resize(batch_size);
// 	std::vector<float> max_error;
// 	max_error.resize(batch_size);
// 	for (int i = 0; i < batch_size; i++) {
// 		max_error[i] = search_dist;
// 	}

// 	std::vector<triangle<float>> As, Bs;
// 	As.resize(batch_size);
// 	Bs.resize(batch_size);

// 	bool done = false;

// 	int d = 0;
// 	int a = 0;
// 	int b = 0;

// 	while (!done) {
// 		int a_child = buckets[d].a.children[a];
// 		int b_child = buckets[d].b.children[b];

// 		triangle<float> A = transform(a_tris[a_child], a_trans).cast<float>();
// 		triangle<float> B = transform(b_tris[b_child], b_trans).cast<float>();

// 		int i = 0;

// 		while (i < batch_size && !done) {
// 			As[i] = A;
// 			Bs[i] = B;

// 			b++;
// 			if (b >= buckets[d].b.children.size()) {
// 				b = 0;
// 				a++;
// 				if (a >= buckets[d].a.children.size()) {
// 					a = 0;
// 					d++;
// 					if (d >= buckets.size()) {
// 						done = true;
// 					}
// 				}
// 				if (!done) {
// 					a_child = buckets[d].a.children[a];
// 					A = transform(a_tris[a_child], a_trans).cast<float>();
// 				}
// 			}
// 			if (!done) {
// 				b_child = buckets[d].b.children[b];
// 				B = transform(b_tris[b_child], b_trans).cast<float>();
// 			}
// 			i++;
// 		}
// 		if (i < batch_size) {
// 			As.resize(i);
// 			Bs.resize(i);
// 		}

// 		penalty_sympy_1_to_1_times_n(As, Bs, Ps, Qs, max_error, failed);
// 		ttd_penatly_total += As.size();

// 		std::vector<triangle<float>> As_failed, Bs_failed;

// 		for (int i = 0; i < As.size(); i++) {
// 			if (failed[i]) {
// 				As_failed.push_back(As[i]);
// 				Bs_failed.push_back(Bs[i]);
// 			}
// 			else if ((Ps[i] - Qs[i]).norm() < a_eps + b_eps) {
// 				contactIO ci;
// 				ci.A = Ps[i];
// 				ci.B = Qs[i];
// 				ci.eps_a = ci.eps_inner_a = a_eps;
// 				ci.eps_b = ci.eps_inner_b = b_eps;
// 				hits.push_back(ci);
// 			}
// 		}

// 		std::vector<std::tuple<float, Eigen::Vector3d, Eigen::Vector3d>> result = fast_1_to_1_times_n(As_failed, Bs_failed);
// 		ttd_total += As_failed.size();

// 		for (auto& r : result) {
// 			if (std::sqrt(std::get<0>(r)) < a_eps + b_eps) {
// 				contactIO ci;
// 				ci.A = std::get<1>(r);
// 				ci.B = std::get<2>(r);
// 				ci.eps_a = ci.eps_inner_a = a_eps;
// 				ci.eps_b = ci.eps_inner_b = b_eps;
// 				hits.push_back(ci);
// 			}
// 		}
// 	}
// }
