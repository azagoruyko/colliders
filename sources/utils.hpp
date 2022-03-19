#pragma once

#include <maya/MPoint.h>

#define MSTR(v) MString(to_string(v).c_str())

class Plane
{
public:
    MPoint orig;
    MVector normal;

    Plane() {};
    Plane(const Plane& plane)
    {
        orig = plane.orig;
        normal = plane.normal;
    };

    Plane(const MPoint &orig, const MVector &normal) : orig(orig), normal(normal.normal()) {};

    MVector projectVector(const MVector vec) const
    {
        return vec - (vec * normal) * normal;
    };

    double distance(const MPoint& src) const
    {
        return (src - orig) * normal;
    };

    MPoint projectPoint(const MPoint& src) const
    {
        double dist = (src - orig) * normal;
        return src - dist * normal;
    };

    MPoint findLineIntersection(const MPoint linePoint, const MVector lineDirection) const
    {
        MVector lineVector = lineDirection.normal();

        double d = ((orig - linePoint) * normal) / (lineVector * normal);
        return linePoint + d * lineVector;
    };
};

inline MVector maxis(const MMatrix& mat, unsigned int index) { return MVector(mat[index][0], mat[index][1], mat[index][2]); }
inline MVector xaxis(const MMatrix &mat) { return maxis(mat, 0); }
inline MVector yaxis(const MMatrix& mat) { return maxis(mat, 1); }
inline MVector zaxis(const MMatrix& mat) { return maxis(mat, 2); }
inline MPoint taxis(const MMatrix& mat) { return maxis(mat, 3); }

inline MMatrix& set_maxis(MMatrix& mat, const unsigned int a, const MVector& v)
{
    mat[a][0] = v.x;
    mat[a][1] = v.y;
    mat[a][2] = v.z;
    return mat;
}

MPointArray findSphereLineIntersection(const MPoint &linePoint, const MVector &lineDirection, const MPoint &sphereCenter, double sphereRadius)
{
    MVector lineVector = lineDirection.normal();

    double a = pow(lineVector.x, 2) + pow(lineVector.y, 2) + pow(lineVector.z, 2);
    double b = 2 * (lineVector.x * (linePoint.x - sphereCenter.x) + lineVector.y * (linePoint.y - sphereCenter.y) + lineVector.z * (linePoint.z - sphereCenter.z));
    double c = pow(linePoint.x - sphereCenter.x, 2) + pow(linePoint.y - sphereCenter.y, 2) + pow(linePoint.z - sphereCenter.z, 2) - pow(sphereRadius, 2);

    double delta = pow(b, 2) - 4 * a * c;
    if (delta <= 0) // when 0 or 1 intersection found
        return MPointArray();

    double d1 = (-b + sqrt(delta)) / 2.0 * a;
    double d2 = (-b - sqrt(delta)) / 2.0 * a;

    MPoint p1(linePoint.x + lineVector.x * d1,
        linePoint.y + lineVector.y * d1,
        linePoint.z + lineVector.z * d1);

    MPoint p2(linePoint.x + lineVector.x * d2,
        linePoint.y + lineVector.y * d2,
        linePoint.z + lineVector.z * d2);

    MPointArray points;
    points.append(p1);
    points.append(p2);
    return points;
}