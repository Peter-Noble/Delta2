//#include <immintrin.h>
#include "fast_dist.h"
#include <iostream>
#include <limits>
#include <array>

#define NOMINMAX
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/avx.h"
#include "contact.h"

//typedef __m128 sseb;
struct sseb
{
	sseb() {};
	sseb(__m128 m) { m128 = m; };
	//sseb(float V[]) { for (int i = 0; i < 4; i++) { v[i] = V[i]; } };
	sseb(double a, double b, double c, double d) { m128 = { float(a), float(b), float(c), float(d) }; };
	// data
	union { __m128 m128; float v[4]; };
};
// operations
const sseb operator &(const sseb& a, const sseb& b) {
	return _mm_and_ps(a.m128, b.m128);
}
sseb& operator&=(sseb& a, const sseb& b) {
	a = _mm_and_ps(a.m128, b.m128);
	return a;
}
const sseb operator |(const sseb& a, const sseb& b) {
	return _mm_or_ps(a.m128, b.m128);
}
const sseb operator ^(const sseb& a, const sseb& b) {
	return _mm_xor_ps(a.m128, b.m128);
}

bool all(const sseb& b) { return _mm_movemask_ps(b.m128) == 0xf; }
bool any(const sseb& b) { return _mm_movemask_ps(b.m128) != 0x0; }
bool none(const sseb& b) { return _mm_movemask_ps(b.m128) == 0x0; }

template<typename T> struct Vec3 {
	// data
	T x, y, z;
};
// operations
template<typename T> Vec3<T> operator+(const Vec3<T>& a, const Vec3 <T>& b) {
	Vec3<T> r;
	r.x = a.x + b.x;
	r.y = a.y + b.y;
	r.z = a.z + b.z;
	return r;
}
template<typename T> Vec3<T> operator-(const Vec3<T>& a, const Vec3 <T>& b) {
	Vec3<T> r;
	r.x = a.x - b.x;
	r.y = a.y - b.y;
	r.z = a.z - b.z;
	return r;
}
template<typename T> Vec3<T> rcp(const Vec3<T>& a) {
	return Vec3<T>( rcp(a.x), rcp(a.y), rcp(a.z));
}
template<typename T> Vec3<T> rsqrt(const Vec3<T>& a) {
	return Vec3<T>(rsqrt(a.x), rsqrt(a.y), rsqrt(a.z));
}
template<typename T> T dot(const Vec3<T>& a, const Vec3<T>& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
template<typename T> Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) {
	Vec3<T> r;
	r.x = a.y * b.z - a.z * b.y;
	r.y = a.z * b.x - a.x * b.z;
	r.z = a.x * b.y - a.y * b.x;
	return r;
}
template<typename T> T length2(const Vec3<T>& a) { return dot(a, a); }
template<typename T> Vec3<T> operator*(const Vec3<T>& a, const T& b) {
	Vec3<T> r;
	r.x = a.x * b;
	r.y = a.y * b;
	r.z = a.z * b;
	return r;
}

struct h {
	int a;
};

typedef sseb simdBool;
typedef sseb simdFloat;
typedef Vec3<simdFloat> simdFloatVec;
typedef std::array<simdFloatVec, 3> simdTriangle_type;
typedef std::array<simdFloatVec, 2> simdLine_type;
typedef simdFloatVec simdPoint_type;


const simdFloat operator *(const simdFloat& a, const simdFloat& b) {
	return _mm_mul_ps(a.m128, b.m128);
}

const simdFloat operator -(const simdFloat& a, const simdFloat& b) {
	return _mm_sub_ps(a.m128, b.m128);
}

const simdFloat operator +(const simdFloat& a, const simdFloat& b) {
	return _mm_add_ps(a.m128, b.m128);
}

const simdFloat operator /(const simdFloat& a, const simdFloat& b) {
	return _mm_div_ps(a.m128, b.m128);
}

const simdBool operator !=(const simdFloat& a, const simdFloat& b) {
	return _mm_cmp_ps(a.m128, b.m128, _CMP_NEQ_OQ);
}

const simdBool operator <=(const simdFloat& a, const simdFloat& b) {
	return _mm_cmp_ps(a.m128, b.m128, _CMP_LE_OQ);
}

const simdBool operator <(const simdFloat& a, const simdFloat& b) {
	return _mm_cmp_ps(a.m128, b.m128, _CMP_LT_OQ);
}

const simdBool operator >=(const simdFloat& a, const simdFloat& b) {
	return _mm_cmp_ps(a.m128, b.m128, _CMP_GE_OQ);
}

const simdBool operator &&(const simdBool& a, const simdBool& b) {
	return _mm_and_ps(a.m128, b.m128);
}

const simdBool And(const simdBool& a, const simdBool& b) {
	return _mm_and_ps(a.m128, b.m128);
}

const simdFloat min(const simdFloat& a, const simdFloat& b) {
	return _mm_min_ps(a.m128, b.m128);
}

const simdFloat min(const simdFloat& a, const simdFloat& b, const simdFloat& c) {
	return min(a, min(b, c));
}

const simdFloat max(const simdFloat& a, const simdFloat& b) {
	return _mm_max_ps(a.m128, b.m128);
}

const simdFloat max(const simdFloat& a, const simdFloat& b, const simdFloat& c) {
	return max(a, max(b, c));
}

template<typename T> T clamp(const T& x, const T& lower, const T& upper) {
	return max(lower, min(x, upper));
}

const simdFloat select(const simdFloat& mask, const simdFloat& a, const simdFloat& b) {
	return _mm_blendv_ps(b.m128, a.m128, mask.m128);
}

const simdFloatVec select(const simdFloat& mask, const simdFloatVec& a, const simdFloatVec& b) {
	simdFloatVec r;
	r.x = select(mask, a.x, b.x);
	r.y = select(mask, a.y, b.y);
	r.z = select(mask, a.z, b.z);
	return r;
}

simdFloat zero = { 0.0, 0.0, 0.0, 0.0 };
simdFloat one = { 1.0, 1.0, 1.0, 1.0 };
simdFloat ulp = { std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::epsilon() };
const simdBool true_mask() {
	return _mm_castsi128_ps(_mm_set1_epi32(-1));
}
const simdBool false_mask() {
	return _mm_setzero_ps();
	// return _mm_castsi128_ps(_mm_set1_epi32(0));
	// return _mm_castsi128_ps(0x0000000000000000);
}

// ie. is ax NOT a separating axis for the 6 points.
inline simdBool simdProject6(const simdFloatVec& ax,
	const simdFloatVec& p1, const simdFloatVec& p2, const simdFloatVec& p3,
	const simdFloatVec& q1, const simdFloatVec& q2, const simdFloatVec& q3) {
	simdFloat P1 = dot(ax, p1);
	simdFloat P2 = dot(ax, p2);
	simdFloat P3 = dot(ax, p3);

	simdFloat Q1 = dot(ax, q1);
	simdFloat Q2 = dot(ax, q2);
	simdFloat Q3 = dot(ax, q3);

	simdFloat mx1 = max(P1, P2, P3);
	simdFloat mn1 = min(P1, P2, P3);
	simdFloat mx2 = max(Q1, Q2, Q3);
	simdFloat mn2 = min(Q1, Q2, Q3);

	return And((mn1 <= mx2), (mn2 <= mx1));
}


//Code is taken from Real Time Collision Detection section 5.1.9 and has been adapted and changed to suit my purposes.
simdFloat simdSegmentSegment2(simdFloatVec& oLine1Point, simdFloatVec& oLine2Point, const simdLine_type& iLine1, const simdLine_type& iLine2) {

	const simdFloatVec dir1 = iLine1[1] - iLine1[0]; // Direction vector of segment S1
	const simdFloatVec dir2 = iLine2[1] - iLine2[0]; // Direction vector of segment S2
	const simdFloatVec r = iLine1[0] - iLine2[0];
	const simdFloat a = dot(dir1, dir1); // Squared length of segment S1, always nonnegative
	simdFloat e = dot(dir2, dir2); // Squared length of segment S2, always nonnegative
	const simdFloat f = dot(dir2, r);
	const simdFloat c = dot(dir1, r);
	const simdFloat b = dot(dir1, dir2);

	// s and t are the parameter values form Line1 and iLine2 respectively.
	simdFloat s, t;
	//The following is always nonnegative.
	simdFloat denom = a * e - b * b;
	// If segments not parallel, compute closest point on L1 to L2, and
	// clamp to segment S1. Else pick arbitrary s (here 0)

	//EVAN: As the previous description says if s can be arbitrary then we take the value given below instead of an if statement.
	//To avoid a nonnegative denominator, we clip it. We know that it always has to be non-negative therefore we clip it with the following value.
	denom = max(denom, simdFloat(ulp));
	s = clamp<simdFloat>((b * f - c * e) / denom, simdFloat(zero), simdFloat(one));
	// Compute point on L2 closest to S1(s) using
	// t = dot((P1+D1*s)-P2,D2) / dot(D2,D2) = (b*s + f) / e
	e = max(e, simdFloat(ulp));
	t = (b * s + f) / e;
	// If t in [0,1] done. Else clamp t, recompute s for the new value
	// of t using s = dot((P2+D2*t)-P1,D1) / dot(D1,D1)= (t*b - c) / a
	// and clamp s to [0, 1]
	const simdFloat newT = clamp(t, simdFloat(zero), simdFloat(one));
	simdBool mask = (newT != t);

	//Now test if all true or none true or some true. Use the select function to choose the respective values.
	s = select(mask, clamp((newT * b - c) / a, simdFloat(zero), simdFloat(one)), s);

	oLine1Point = iLine1[0] + dir1 * s;
	oLine2Point = iLine2[0] + dir2 * newT;
	return length2(oLine1Point - oLine2Point);
}

simdBool simdTriContact(const simdFloatVec& P1, const simdFloatVec& P2, const simdFloatVec& P3, const simdFloatVec& Q1, const simdFloatVec& Q2, const simdFloatVec& Q3) {
	const simdFloatVec p1 = P1 - P1; // TODO ZeroTy();
	const simdFloatVec p2 = P2 - P1;
	const simdFloatVec p3 = P3 - P1;

	const simdFloatVec q1 = Q1 - P1;
	const simdFloatVec q2 = Q2 - P1;
	const simdFloatVec q3 = Q3 - P1;

	const simdFloatVec e1 = P2 - P1;
	const simdFloatVec e2 = P3 - P2;

	const simdFloatVec f1 = Q2 - Q1;
	const simdFloatVec f2 = Q3 - Q2;

	simdBool mask = true_mask();

	const simdFloatVec n1 = cross(e1, e2);    mask = mask & simdProject6(n1, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec m1 = cross(f1, f2);    mask = mask & simdProject6(m1, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec ef11 = cross(e1, f1);  mask = mask & simdProject6(ef11, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec ef12 = cross(e1, f2);  mask = mask & simdProject6(ef12, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec f3 = q1 - q3;
	const simdFloatVec ef13 = cross(e1, f3);  mask = mask & simdProject6(ef13, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec ef21 = cross(e2, f1);  mask = mask & simdProject6(ef21, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec ef22 = cross(e2, f2);  mask = mask & simdProject6(ef22, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec ef23 = cross(e2, f3);  mask = mask & simdProject6(ef23, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec e3 = p1 - p3;
	const simdFloatVec ef31 = cross(e3, f1);  mask = mask & simdProject6(ef31, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec ef32 = cross(e3, f2);  mask = mask & simdProject6(ef32, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec ef33 = cross(e3, f3);  mask = mask & simdProject6(ef33, p1, p2, p3, q1, q2, q3); if (none(mask)) return mask;
	const simdFloatVec g1 = cross(e1, n1);    mask = mask & simdProject6(g1, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec g2 = cross(e2, n1);    mask = mask & simdProject6(g2, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec g3 = cross(e3, n1);    mask = mask & simdProject6(g3, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec h1 = cross(f1, m1);    mask = mask & simdProject6(h1, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec h2 = cross(f2, m1);    mask = mask & simdProject6(h2, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	const simdFloatVec h3 = cross(f3, m1);    mask = mask & simdProject6(h3, p1, p2, p3, q1, q2, q3);   if (none(mask)) return mask;
	return mask;
}

simdBool closestEdgePoints(const simdFloatVec& iTri1Pt, const simdFloatVec& iClosestPtToTri1, const simdFloatVec& iTri2Pt, const simdFloatVec& iClosestPtToTri2, const simdFloatVec& iSepDir) {
	simdFloatVec awayDirection = iTri1Pt - iClosestPtToTri1;
	const simdFloat isDiffDirection = dot(awayDirection, iSepDir);

	awayDirection = iTri2Pt - iClosestPtToTri2;
	const simdFloat isSameDirection = dot(awayDirection, iSepDir);

	return (isDiffDirection <= simdFloat(zero)) & (isSameDirection >= simdFloat(zero));
}

simdFloat closestEdgeToEdge(simdBool& oIsFinished, simdFloatVec& oTriAPoint, simdFloatVec& oTriBPoint, const simdLine_type(&iTriAEdges)[3], const simdLine_type& iTriBEdge, const simdFloatVec& iTriBLastPt) {
	//Test the triangle edge against all three edges of the triangle iTriA.
	simdFloatVec A1p, A2p, A3p, B1p, B2p, B3p, separatingDir;

	const simdFloat A = simdSegmentSegment2(oTriAPoint, oTriBPoint, iTriAEdges[0], iTriBEdge);
	//Test to see if the distances found so far were the closest:
	separatingDir = oTriBPoint - oTriAPoint;
	oIsFinished = oIsFinished | closestEdgePoints(iTriAEdges[1][0], oTriAPoint, iTriBLastPt, oTriBPoint, separatingDir);
	if (all(oIsFinished))
		return A;

	const simdFloat B = simdSegmentSegment2(A2p, B2p, iTriAEdges[1], iTriBEdge);
	separatingDir = B2p - A2p;
	oIsFinished = oIsFinished | closestEdgePoints(iTriAEdges[2][0], A2p, iTriBLastPt, B2p, separatingDir);

	const simdBool AB = A < B;
	const simdFloat ABdist = select(AB, A, B);
	oTriAPoint = select(AB, oTriAPoint, A2p);
	oTriBPoint = select(AB, oTriBPoint, B2p);

	if (all(oIsFinished))
		return ABdist;

	const simdFloat C = simdSegmentSegment2(A3p, B3p, iTriAEdges[2], iTriBEdge);
	separatingDir = B3p - A3p;
	oIsFinished = oIsFinished | closestEdgePoints(iTriAEdges[0][0], A3p, iTriBLastPt, B3p, separatingDir);

	const simdBool ABC = ABdist < C;
	oTriAPoint = select(ABC, oTriAPoint, A3p);
	oTriBPoint = select(ABC, oTriBPoint, B3p);

	return select(ABC, ABdist, C);
}

const simdFloat simdTriPoint2(simdFloatVec& oTriPoint, const simdTriangle_type& iTri, const simdPoint_type& iPoint) {
	const simdFloatVec ab = iTri[1] - iTri[0];
	const simdFloatVec ac = iTri[2] - iTri[0];
	const simdFloatVec ap = iPoint - iTri[0];
	const simdFloat d1 = dot(ab, ap);
	const simdFloat d2 = dot(ac, ap);
	const simdBool mask1 = (d1 <= simdFloat(zero)) & (d2 <= simdFloat(zero));
	oTriPoint = iTri[0];
	simdBool exit(mask1);
	if (all(exit))
		return length2(oTriPoint - iPoint);
	
	const simdFloatVec bp = iPoint - iTri[1];
	const simdFloat d3 = dot(ab, bp);
	const simdFloat d4 = dot(ac, bp);
	const simdBool mask2 = (d3 >= simdFloat(zero)) & (d4 <= d3);
	//Closest point is the point iTri[1]. Update if necessary.
	oTriPoint = select(exit, oTriPoint, select(mask2, iTri[1], oTriPoint));
	exit = exit | mask2;
	if (all(exit))
		return length2(oTriPoint - iPoint);
	
	const simdFloatVec cp = iPoint - iTri[2];
	const simdFloat d5 = dot(ab, cp);
	const simdFloat d6 = dot(ac, cp);
	const simdBool mask3 = (d6 >= simdFloat(zero)) & (d5 <= d6);
	//Closest point is the point iTri[2]. Update if necessary.
	oTriPoint = select(exit, oTriPoint, select(mask3, iTri[2], oTriPoint));
	exit = exit | mask3;
	if (all(exit))
		return length2(oTriPoint - iPoint);
	
	const simdFloat vc = d1 * d4 - d3 * d2;
	const simdBool mask4 = (vc <= simdFloat(zero)) & (d1 >= simdFloat(zero)) & (d3 <= simdFloat(zero));
	const simdFloat v1 = d1 / (d1 - d3);
	const simdFloatVec answer1 = iTri[0] + ab * v1;
	//Closest point is on the line ab. Update if necessary.
	oTriPoint = select(exit, oTriPoint, select(mask4, answer1, oTriPoint));
	exit = exit | mask4;
	if (all(exit))
		return length2(oTriPoint - iPoint);
	
	const simdFloat vb = d5 * d2 - d1 * d6;
	const simdBool mask5 = (vb <= simdFloat(zero)) & (d2 >= simdFloat(zero)) & (d6 <= simdFloat(zero));
	const simdFloat w1 = d2 / (d2 - d6);
	const simdFloatVec answer2 = iTri[0] + ac * w1;
	//Closest point is on the line ac. Update if necessary.
	oTriPoint = select(exit, oTriPoint, select(mask5, answer2, oTriPoint));
	exit = exit | mask5;
	if (all(exit))
		return length2(oTriPoint - iPoint);
	
	const simdFloat va = d3 * d6 - d5 * d4;
	const simdBool mask6 = (va <= simdFloat(zero)) & ((d4 - d3) >= simdFloat(zero)) & ((d5 - d6) >= simdFloat(zero));
	simdFloat w2 = (d4 - d3) / ((d4 - d3) + (d5 - d6));
	const simdFloatVec answer3 = iTri[1] + (iTri[2] - iTri[1]) * w2;
	//Closest point is on the line bc. Update if necessary.
	oTriPoint = select(exit, oTriPoint, select(mask6, answer3, oTriPoint));
	exit = exit | mask6;
	if (all(exit))
		return length2(oTriPoint - iPoint);
	
	const simdFloat denom = simdFloat(one) / (va + vb + vc);
	const simdFloat v2 = vb * denom;
	const simdFloat w3 = vc * denom;
	const simdFloatVec answer4 = iTri[0] + ab * v2 + ac * w3;
	const simdBool mask7 = length2(answer4 - iPoint) < length2(oTriPoint - iPoint);
	//Closest point is inside triangle. Update if necessary.
	oTriPoint = select(exit, oTriPoint, select(mask7, answer4, oTriPoint));
	return length2(oTriPoint - iPoint);
}


simdFloat closestVertToTri(simdFloatVec& oTriAPoint, simdFloatVec& oTriBPoint, const simdTriangle_type& iTriA, const simdTriangle_type& iTriB) {
	simdFloatVec Ap, Bp, Cp;

	const simdFloatVec edge[2] = { iTriA[1] - iTriA[0], iTriA[2] - iTriA[1] };
	simdFloatVec TriNormal = cross(edge[1], edge[0]);
	const simdFloat norm2 = length2(TriNormal);

	/*const simdFloat A = simdTriPoint2(Ap, TriNormal, norm2, iTriA, iTriB[0]);
	const simdFloat B = simdTriPoint2(Bp, TriNormal, norm2, iTriA, iTriB[1]);
	const simdFloat C = simdTriPoint2(Cp, TriNormal, norm2, iTriA, iTriB[2]);*/

	const simdFloat A = simdTriPoint2(Ap, iTriA, iTriB[0]);
	const simdFloat B = simdTriPoint2(Bp, iTriA, iTriB[1]);
	const simdFloat C = simdTriPoint2(Cp, iTriA, iTriB[2]);

	const simdBool AB = A < B;
	const simdFloat ABdist = select(AB, A, B);
	const simdFloatVec ABp = select(AB, Ap, Bp);

	const simdBool ABC = ABdist < C;
	oTriAPoint = select(ABC, ABp, Cp);
	oTriBPoint = select(ABC, select(AB, iTriB[0], iTriB[1]), iTriB[2]);
	return select(ABC, ABdist, C);
}

//Return the four floats (one simdFloat) and two sets of four closest points (the two simdFloatVec) between the two sets of triangles.
simdFloat simdTriTriOpt2(simdFloatVec& oTri1Point, simdFloatVec& oTri2Point, const simdTriangle_type& iTri1, const simdTriangle_type& iTri2) {

	//The three edges of the triangle. Keep orientation consistent.
	const simdLine_type tri1Edges[3] = { { iTri1[1], iTri1[0] }, { iTri1[2], iTri1[1] }, { iTri1[0], iTri1[2] } };
	const simdLine_type tri2Edges[3] = { { iTri2[1], iTri2[0] }, { iTri2[2], iTri2[1] }, { iTri2[0], iTri2[2] } };

	simdFloatVec tri1Vector, tri2Vector;
	simdBool isFinished = false_mask();

	simdFloat minDistsTriTri = closestEdgeToEdge(isFinished, oTri1Point, oTri2Point, tri1Edges, tri2Edges[0], iTri2[2]);
	if (all(isFinished))
		return minDistsTriTri;

	simdFloat tmpMinDist = closestEdgeToEdge(isFinished, tri1Vector, tri2Vector, tri1Edges, tri2Edges[1], iTri2[0]);
	simdBool mask = tmpMinDist < minDistsTriTri;
	minDistsTriTri = select(mask, tmpMinDist, minDistsTriTri);
	oTri1Point = select(mask, tri1Vector, oTri1Point);
	oTri2Point = select(mask, tri2Vector, oTri2Point);
	if (all(isFinished))
		return minDistsTriTri;

	tmpMinDist = closestEdgeToEdge(isFinished, tri1Vector, tri2Vector, tri1Edges, tri2Edges[2], iTri2[1]);
	mask = tmpMinDist < minDistsTriTri;
	minDistsTriTri = select(mask, tmpMinDist, minDistsTriTri);
	oTri1Point = select(mask, tri1Vector, oTri1Point);
	oTri2Point = select(mask, tri2Vector, oTri2Point);
	if (all(isFinished))
		return minDistsTriTri;

	//Now do vertex-triangle distances.
	tmpMinDist = closestVertToTri(tri2Vector, tri1Vector, iTri2, iTri1);
	mask = tmpMinDist < minDistsTriTri;
	oTri1Point = select(mask, tri1Vector, oTri1Point);
	oTri2Point = select(mask, tri2Vector, oTri2Point);
	minDistsTriTri = select(mask, tmpMinDist, minDistsTriTri);

	tmpMinDist = closestVertToTri(tri1Vector, tri2Vector, iTri1, iTri2);
	mask = tmpMinDist < minDistsTriTri;
	oTri1Point = select(mask, tri1Vector, oTri1Point);
	oTri2Point = select(mask, tri2Vector, oTri2Point);

	minDistsTriTri = select(mask, tmpMinDist, minDistsTriTri);
	//We need to rule out the triangles colliding with each other otherwise we can get a distance that is not equal to 0 although
	//the true distance is 0. Hence we use simdTriContact here.

	simdBool colliding = simdTriContact(iTri1[0], iTri1[1], iTri1[2], iTri2[0], iTri2[1], iTri2[2]);
	return select(colliding, simdFloat(zero), minDistsTriTri);
	//return select(colliding, minDistsTriTri, zero);
}

void fastBucket(
	Delta2::common::Triangle<float>& a_parent_tri,
	float a_parent_eps,
	Delta2::common::Triangle<float>& b_parent_tri,
	float b_parent_eps,
	std::vector<int>& a_children,
	std::vector<int>& b_children,
	const std::vector<Delta2::common::Triangle<float>>& a_tris,
	const std::vector<Delta2::common::Triangle<float>>& b_tris,
	Eigen::Matrix4f a_trans,
	Eigen::Matrix4f b_trans,
	float eps_a,
	float eps_b,
	std::vector<Delta2::collision::Contact<float>>& result,
	Delta2::Particle& p_a,
	Delta2::Particle& p_b)
{
	// TODO do computations in a_local space to reduce number of transforms
	const int fast_vec_width = 4;
	const int a_groups = (a_children.size() + fast_vec_width - 1) / fast_vec_width;
	std::vector<simdTriangle_type> __As;
	__As.reserve(a_children.size());
	for (int a = 0; a < a_children.size(); a+=4) {
		Delta2::common::Triangle<float> a0, a1, a2, a3;
		a0 = a_tris[a_children[a]].transformed(a_trans);
		if (a + 1 < a_children.size()) {
			a1 = a_tris[a_children[a + 1]].transformed(a_trans);
		}
		else {
			a1 = a0;
		}
		if (a + 2 < a_children.size()) {
			a2 = a_tris[a_children[a+2]].transformed(a_trans);
		}
		else {
			a2 = a0;
		}
		if (a + 3 < a_children.size()) {
			a3 = a_tris[a_children[a+3]].transformed(a_trans);
		}
		else {
			a3 = a0;
		}
		//printf("1st: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", a0.xs[0], a0.ys[0], a0.zs[0], a0.xs[1], a0.ys[1], a0.zs[1], a0.xs[2], a0.ys[2], a0.zs[2]);

		simdTriangle_type __a;
		for (int vert = 0; vert < 3; vert++) {
			simdPoint_type __vert;
			__vert.x = simdFloat(a0[vert].x(), a1[vert].x(), a2[vert].x(), a3[vert].x());
			__vert.y = simdFloat(a0[vert].y(), a1[vert].y(), a2[vert].y(), a3[vert].y());
			__vert.z = simdFloat(a0[vert].z(), a1[vert].z(), a2[vert].z(), a3[vert].z());
			__a[vert] = __vert;
		}
		__As.push_back(__a);
	}

	const float search_dist = eps_a + eps_b;

	for (int b = 0; b < b_children.size(); b++) {
		Delta2::common::Triangle<float> B = b_tris[b_children[b]].transformed(b_trans);
		//printf("2nd: (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", B.xs[0], B.ys[0], B.zs[0], B.xs[1], B.ys[1], B.zs[1], B.xs[2], B.ys[2], B.zs[2]);
		simdTriangle_type __B;
		for (int i = 0; i < 3; i++) {
			simdPoint_type __vert;
			__vert.x = simdFloat(B[i].x(), B[i].x(), B[i].x(), B[i].x());
			__vert.y = simdFloat(B[i].y(), B[i].y(), B[i].y(), B[i].y());
			__vert.z = simdFloat(B[i].z(), B[i].z(), B[i].z(), B[i].z());
			__B[i] = __vert;
		}
		simdFloatVec oTri1Point, oTri2Point;
		for (int a = 0; a < __As.size(); a++) {
			simdFloat __dists2 = simdTriTriOpt2(oTri1Point, oTri2Point, __As[a], __B);
			// TODO compare any(__dists2 < search) before continuing
			float dists2[4];
			_mm_store_ps(dists2, __dists2.m128);

			float Px[4];  _mm_store_ps(Px, oTri1Point.x.m128);
			float Py[4];  _mm_store_ps(Py, oTri1Point.y.m128);
			float Pz[4];  _mm_store_ps(Pz, oTri1Point.z.m128);
			float Qx[4];  _mm_store_ps(Qx, oTri2Point.x.m128);
			float Qy[4];  _mm_store_ps(Qy, oTri2Point.y.m128);
			float Qz[4];  _mm_store_ps(Qz, oTri2Point.z.m128);

			for (int i = 0; i < 4; i++) {
				if (a * 4 + i < a_children.size()) {
					float dist2 = dists2[i];

					if (dist2 == 0) { // TODO handle zero distances better
						Eigen::Vector3f P = Eigen::Vector3f({ 0, 0, 0 });
						Eigen::Vector3f Q = Eigen::Vector3f({ 0, 0, 0 });
						result.emplace_back(P, Q, eps_a, eps_b, 0.0, 0.0, p_a, p_b);
					}
					else if (dist2 < search_dist * search_dist) {
						Eigen::Vector3f P = Eigen::Vector3f({ Px[i], Py[i], Pz[i] });
						Eigen::Vector3f Q = Eigen::Vector3f({ Qx[i], Qy[i], Qz[i] });
						result.emplace_back(P, Q, eps_a, eps_b, 0.0, 0.0, p_a, p_b);
					}
				}
			}
		}
	}
	// TODO might be beneficial to filter results here and then add that to result to take adanvantage of collision being very local at this point.
}

// Returns <distance, contact point on A, contact point on B>
std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> Delta2::collision::_fastInnerDists(
	Delta2::common::Triangle<float>& a_tri,
	std::vector<Delta2::common::Triangle<float>>& b_tris)
{
	const int fast_vec_width = 4;

	std::vector<std::tuple<float, Eigen::Vector3f, Eigen::Vector3f>> result;
	result.reserve(b_tris.size());

	simdTriangle_type __a;
	for (int i = 0; i < 3; i++) {
		simdPoint_type __vert;
		__vert.x = simdFloat(a_tri[i].x(), a_tri[i].x(), a_tri[i].x(), a_tri[i].x());
		__vert.y = simdFloat(a_tri[i].y(), a_tri[i].y(), a_tri[i].y(), a_tri[i].y());
		__vert.z = simdFloat(a_tri[i].z(), a_tri[i].z(), a_tri[i].z(), a_tri[i].z());
		__a[i] = __vert;
	}

	for (int a = 0; a < b_tris.size(); a += 4) {
		Delta2::common::Triangle<float> a0, a1, a2, a3;
		a0 = b_tris[a];
		if (a + 1 < b_tris.size()) {
			a1 = b_tris[a + 1];
		}
		else {
			a1 = a0;
		}
		if (a + 2 < b_tris.size()) {
			a2 = b_tris[a + 2];
		}
		else {
			a2 = a0;
		}
		if (a + 3 < b_tris.size()) {
			a3 = b_tris[a + 3];
		}
		else {
			a3 = a0;
		}

		simdTriangle_type __b;
		for (int vert = 0; vert < 3; vert++) {
			simdPoint_type __vert;
			__vert.x = simdFloat(a0[vert].x(), a1[vert].x(), a2[vert].x(), a3[vert].x());
			__vert.y = simdFloat(a0[vert].y(), a1[vert].y(), a2[vert].y(), a3[vert].y());
			__vert.z = simdFloat(a0[vert].z(), a1[vert].z(), a2[vert].z(), a3[vert].z());
			__b[vert] = __vert;
		}

		simdFloatVec oTri1Point, oTri2Point;

		simdFloat __dists2 = simdTriTriOpt2(oTri1Point, oTri2Point, __a, __b);

		float dists2[4];
		_mm_store_ps(dists2, __dists2.m128);

		float Px[4];  _mm_store_ps(Px, oTri1Point.x.m128);
		float Py[4];  _mm_store_ps(Py, oTri1Point.y.m128);
		float Pz[4];  _mm_store_ps(Pz, oTri1Point.z.m128);
		float Qx[4];  _mm_store_ps(Qx, oTri2Point.x.m128);
		float Qy[4];  _mm_store_ps(Qy, oTri2Point.y.m128);
		float Qz[4];  _mm_store_ps(Qz, oTri2Point.z.m128);

		int num_results = 4;
		if (a + 4 > b_tris.size()) {
			num_results = b_tris.size() % 4;
		}
		for (int i = 0; i < num_results; i++) {
			Eigen::Vector3f P = Eigen::Vector3f({ Px[i], Py[i], Pz[i] });
			Eigen::Vector3f Q = Eigen::Vector3f({ Qx[i], Qy[i], Qz[i] });
			result.push_back(std::make_tuple(sqrt(dists2[i]), P, Q));
		}
	}
	return result;
}

std::vector<std::tuple<float, Eigen::Vector3d, Eigen::Vector3d>> fast1To1TimesN(
	std::vector<Delta2::common::Triangle<float>>& a_tris,
	std::vector<Delta2::common::Triangle<float>>& b_tris)
{
	const int fast_vec_width = 4;

	std::vector<std::tuple<float, Eigen::Vector3d, Eigen::Vector3d>> result;
	result.reserve(b_tris.size());

	int batches = (a_tris.size() + 4 - 1) / 4;

	for (int batch = 0; batch < batches; batch++) {
		simdTriangle_type __a;
		for (int i = 0; i < 3; i++) {
			simdPoint_type __vert;
			int one = std::min(batch * 4 + 0, int(a_tris.size() - 1));
			int two = std::min(batch * 4 + 1, int(a_tris.size() - 1));
			int three = std::min(batch * 4 + 2, int(a_tris.size() - 1));
			int four = std::min(batch * 4 + 3, int(a_tris.size() - 1));
			__vert.x = simdFloat(a_tris[one][i].x(), a_tris[two][i].x(), a_tris[three][i].x(), a_tris[four][i].x());
			__vert.y = simdFloat(a_tris[one][i].y(), a_tris[two][i].y(), a_tris[three][i].y(), a_tris[four][i].y());
			__vert.z = simdFloat(a_tris[one][i].z(), a_tris[two][i].z(), a_tris[three][i].z(), a_tris[four][i].z());
			__a[i] = __vert;
		}

		simdTriangle_type __b;
		for (int i = 0; i < 3; i++) {
			simdPoint_type __vert;
			int one = std::min(batch * 4 + 0, int(b_tris.size() - 1));
			int two = std::min(batch * 4 + 1, int(b_tris.size() - 1));
			int three = std::min(batch * 4 + 2, int(b_tris.size() - 1));
			int four = std::min(batch * 4 + 3, int(b_tris.size() - 1));
			__vert.x = simdFloat(b_tris[one][i].x(), b_tris[two][i].x(), b_tris[three][i].x(), b_tris[four][i].x());
			__vert.y = simdFloat(b_tris[one][i].y(), b_tris[two][i].y(), b_tris[three][i].y(), b_tris[four][i].y());
			__vert.z = simdFloat(b_tris[one][i].z(), b_tris[two][i].z(), b_tris[three][i].z(), b_tris[four][i].z());
			__b[i] = __vert;
		}

		simdFloatVec oTri1Point, oTri2Point;

		simdFloat __dists2 = simdTriTriOpt2(oTri1Point, oTri2Point, __a, __b);

		float dists2[4];
		_mm_store_ps(dists2, __dists2.m128);

		float Px[4];  _mm_store_ps(Px, oTri1Point.x.m128);
		float Py[4];  _mm_store_ps(Py, oTri1Point.y.m128);
		float Pz[4];  _mm_store_ps(Pz, oTri1Point.z.m128);
		float Qx[4];  _mm_store_ps(Qx, oTri2Point.x.m128);
		float Qy[4];  _mm_store_ps(Qy, oTri2Point.y.m128);
		float Qz[4];  _mm_store_ps(Qz, oTri2Point.z.m128);

		for (int i = 0; i < 4 && 4 * batch + i < a_tris.size(); i++) {
			Eigen::Vector3d P = Eigen::Vector3d({ Px[i], Py[i], Pz[i] });
			Eigen::Vector3d Q = Eigen::Vector3d({ Qx[i], Qy[i], Qz[i] });
			result.push_back(std::make_tuple(dists2[i], P, Q));
		}
	}

	return result;
}

// This probably isn't fast at all as it's having to pad one pair of triangles up to 4.
// TODO An unvectorised version of this might be faster.
float fastIndividualDist(
	Delta2::common::Triangle<float>& a_tri,
	Delta2::common::Triangle<float>& b_tri,
	Eigen::Vector3f& P,
	Eigen::Vector3f& Q)
{
	simdTriangle_type __a;
	for (int i = 0; i < 3; i++) {
		simdPoint_type __vert;
		__vert.x = simdFloat(a_tri[i].x(), a_tri[i].x(), a_tri[i].x(), a_tri[i].x());
		__vert.y = simdFloat(a_tri[i].y(), a_tri[i].y(), a_tri[i].y(), a_tri[i].y());
		__vert.z = simdFloat(a_tri[i].z(), a_tri[i].z(), a_tri[i].z(), a_tri[i].z());
		__a[i] = __vert;
	}
	
	simdTriangle_type __b;
	for (int i = 0; i < 3; i++) {
		simdPoint_type __vert;
		__vert.x = simdFloat(b_tri[i].x(), b_tri[i].x(), b_tri[i].x(), b_tri[i].x());
		__vert.y = simdFloat(b_tri[i].y(), b_tri[i].y(), b_tri[i].y(), b_tri[i].y());
		__vert.z = simdFloat(b_tri[i].z(), b_tri[i].z(), b_tri[i].z(), b_tri[i].z());
		__b[i] = __vert;
	}

	simdFloatVec oTri1Point, oTri2Point;

	simdFloat __dists2 = simdTriTriOpt2(oTri1Point, oTri2Point, __a, __b);

	float dists2[4];
	_mm_store_ps(dists2, __dists2.m128);

	if (dists2[0] <= std::numeric_limits<float>::epsilon()) {
		P = { 0, 0, 0 };
		Q = { 0, 0, 0 };
	}
	else {
		float Px[4];  _mm_store_ps(Px, oTri1Point.x.m128);
		float Py[4];  _mm_store_ps(Py, oTri1Point.y.m128);
		float Pz[4];  _mm_store_ps(Pz, oTri1Point.z.m128);
		float Qx[4];  _mm_store_ps(Qx, oTri2Point.x.m128);
		float Qy[4];  _mm_store_ps(Qy, oTri2Point.y.m128);
		float Qz[4];  _mm_store_ps(Qz, oTri2Point.z.m128);

		P = { Px[0], Py[0], Pz[0] };
		Q = { Qx[0], Qy[0], Qz[0] };
	}

	return dists2[0];
}
