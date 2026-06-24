#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MQuaternion.h>
#include <maya/MEulerRotation.h>
#include <maya/MIntArray.h>
#include <maya/MPointArray.h>
#include <maya/MDoubleArray.h>
#include <cmath>

#include "bellColliderSolver.h"
#include "utils.hpp"

using namespace std;

static double wrapParam(double param, double minParam, double maxParam)
{
    double range = maxParam - minParam;
    if (range <= 0.0) return minParam;
    double p = param - minParam;
    p = fmod(p, range);
    if (p < 0.0) p += range;
    return p + minParam;
}

MObject BellColliderSolver::makeBellMesh(const MMatrix& matrix, unsigned int axis, unsigned int numSides, double height, double bottomRadius, double topRadius)
{
    const int numVertices = numSides * 2 + 1;
    const int numPolygons = numSides * 2;

    MPointArray vertexArray;
    MIntArray polygonCounts, polygonConnects;

    vertexArray.append(MPoint(0, 0, 0) * matrix);

    // bottom
    for (int i = 0; i < numSides; i++)
    {
        const double rad = (double)i / numSides * 2 * M_PI;
        const double x = bottomRadius * cos(rad);
        const double z = bottomRadius * sin(rad);
        
        MPoint p;
        switch (axis)
        {
        case 0: p = MPoint(0, x, z); break;
        case 1: p = MPoint(x, 0, z); break;
        case 2: p = MPoint(x, z, 0); break;
        }

        vertexArray.append(p * matrix);

        polygonCounts.append(3);
        polygonConnects.append(0);
        polygonConnects.append(i + 1); // i=0 is center
        polygonConnects.append(i == numSides - 1 ? 1 : i + 2);
    }

    // top
    for (int i = 0; i < numSides; i++)
    {
        const double rad = (double)i / numSides * 2 * M_PI;
        const double x = topRadius * cos(rad);
        const double z = topRadius * sin(rad);

        MPoint p;
        switch (axis)
        {
        case 0: p = MPoint(height, x, z); break;
        case 1: p = MPoint(x, height, z); break;
        case 2: p = MPoint(x, z, height); break;
        }
        vertexArray.append(p * matrix);

        polygonCounts.append(4);
        polygonConnects.append(i + 1); // i=0 is center
        polygonConnects.append(numSides + i + 1);
        polygonConnects.append(i == numSides - 1 ? numSides + 1 : numSides + i + 2);
        polygonConnects.append(i == numSides - 1 ? 1 : i + 2);
    }

    MFnMeshData meshData;
    MObject meshObject = meshData.create();

    MFnMesh meshFn;
    meshFn.create(numVertices, numPolygons, vertexArray, polygonCounts, polygonConnects, meshObject);
    return meshObject;
}

static MObject makeBellCurve(const MPointArray &points, int bellSubdivision, bool use_bottom=false)
{
    const int START = use_bottom ? 1 : bellSubdivision + 1;
    const int END = use_bottom ? bellSubdivision : points.length();

    MPointArray cvs;
    MDoubleArray knots;
    for (int i = START; i < END; i++)
    {
        knots.append(cvs.length());
        cvs.append(points[i]);
    }

    knots.append(knots[knots.length() - 1] + 1);
    cvs.append(cvs[0]);

    MFnNurbsCurveData curveDataFn;
    MObject curveData = curveDataFn.create();

    MFnNurbsCurve curveFn;
    curveFn.create(cvs, knots, 1, MFnNurbsCurve::kPeriodic, false, false, curveData);
    return curveData;
}

void BellColliderSolver::deformPoints(const BellColliderInputs& inputs, const MPointArray& baseBellPoints, const Plane& bellPlane, vector<MPointArray>& bellPointsList)
{
    const MMatrix bellMatrix = inputs.bellMatrix;
    const MMatrix bellMatrixInverse = bellMatrix.inverse();
    const int bellSubdivision = inputs.bellSubdivision;
    const float falloff = inputs.falloff;
    const float collision = inputs.collision;

    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1); // Y axis
    const MVector bellNormal = bellAxis.normal();

    bellPointsList.clear();

    for (const auto& ringMatrix : inputs.ringMatrices)
    {
        const MMatrix ringMatrixInverse = ringMatrix.inverse();

        const MVector ringDirection = maxis(ringMatrix, 1); // Y axis
        const MPoint ring_translate = taxis(ringMatrix);
        const MVector ringNormal = ringDirection.normal();
        const Plane ringPlane(ring_translate, ringNormal);

        const MPoint ring_translate_proj = bellPlane.projectPoint(taxis(ringMatrix));
        const MVector ringDirection_proj = bellPlane.projectVector(ringDirection);

        MPointArray bellPoints = baseBellPoints;

        if (ringDirection_proj.length() > 1e-3)
        {
            const MPointArray hitPoints = findSphereLineIntersection(ring_translate_proj * bellMatrixInverse, ringDirection_proj * bellMatrixInverse, MPoint(0,0,0), 1.001);

            MPoint collisionPointBell, collisionPointRing;
            if (hitPoints.length() > 0)
            {
                collisionPointBell = hitPoints[0] * bellMatrix + bellAxis;

                const double linePointCoeff = ringNormal * bellNormal > 0 ? 1 : -1;

                const MVector ring_proj = ringPlane.projectVector(ringDirection_proj * linePointCoeff);
                double delta = 1.0;
                MVector ring_proj_scaled(0, 0, 0);
                double ring_proj_len = ring_proj.length();
                if (ring_proj_len > 1e-5)
                {
                    double local_len = (ring_proj * ringMatrixInverse).length();
                    if (local_len > 1e-5)
                        delta = ring_proj_len / local_len;
                    ring_proj_scaled = ring_proj.normal() * delta;
                }
                const MPoint linePoint = ring_translate + ring_proj_scaled;

                const MPointArray sphereLinePoints = findSphereLineIntersection(linePoint, ringDirection, ring_translate, (collisionPointBell - ring_translate).length());

                for (int k = 0; k < sphereLinePoints.length(); k++)
                {
                    if ((sphereLinePoints[k] - ring_translate) * ringDirection > 0)
                        collisionPointRing = sphereLinePoints[k];
                }
            }

            double bellAxisLen = bellAxis.length();
            const double collisionDelta = bellAxisLen > 1e-5 ? (bellPlane.distance(collisionPointRing) - bellPlane.distance(collisionPointBell)) / bellAxisLen : 0.0;
            if (collisionDelta < 0)
            {
                MTransformationMatrix rotationMatrixFn;
                rotationMatrixFn.setTranslation(ring_translate, MSpace::kWorld);
                const MMatrix rotateMatrixInverse = rotationMatrixFn.asMatrixInverse();

                const MQuaternion quat(collisionPointBell - ring_translate, collisionPointRing - ring_translate); // rotate X to final point
                rotationMatrixFn.rotateBy(quat, MSpace::kTransform);
                const MMatrix rotateMatrix = rotationMatrixFn.asMatrix();

                const Plane upperBellPlane(bell_translate + bellAxis, bellNormal);

                // bell top deformation
                for (int j = bellSubdivision + 1; j < (int)bellPoints.length(); j++)
                {
                    const MPoint bellPoint_proj = bellPlane.projectPoint(bellPoints[j]);
                    const MVector offset_proj = bellPoint_proj * bellMatrixInverse - ring_translate_proj * bellMatrixInverse;

                    double weight = offset_proj.normal() * (ringDirection_proj * bellMatrixInverse).normal(); // -1..1

                    if (weight > falloff)
                    {
                        double divisor = 1.0 - falloff;
                        weight = divisor > 1e-5 ? (weight - falloff) / divisor : 1.0;

                        const MPoint rp = bellPoints[j] * rotateMatrixInverse * rotateMatrix;

                        // constrain by upper plane
                        MPoint p = rp;
                        if (upperBellPlane.distance(rp) > 0)
                            p = upperBellPlane.projectPoint(rp);

                        bellPoints[j] = p * weight + bellPoints[j] * (1.0 - weight);                               
                    }
                }
            }
        }

        // ring collision
        if (collision > 1e-5)
        {
            for (int j = bellSubdivision + 1; j < (int)bellPoints.length(); j++)
            {                    
                const MVector vec = ringPlane.projectPoint(bellPoints[j]) - ring_translate;
                double vec_len = vec.length();
                if (vec_len > 1e-5)
                {
                    double local_len = (vec * ringMatrixInverse).length();
                    double delta = local_len > 1e-5 ? vec_len / local_len : 1.0;
                    const MVector vec_proj_scaled = vec.normal() * delta; // scale vector

                    if (vec_proj_scaled.length() > vec_len)
                        bellPoints[j] += vec.normal() * (vec_proj_scaled.length() - vec_len) * collision;
                }
            }
        }

        bellPointsList.push_back(bellPoints);
    }
}

void BellColliderSolver::averageDisplacements(int bellSubdivision, const MPointArray& baseBellPoints, const vector<MPointArray>& bellPointsList, MPointArray& outBellPoints)
{
    outBellPoints = baseBellPoints;

    for (size_t i = bellSubdivision + 1; i < baseBellPoints.length(); i++)
    {
        double sum = 0;
        double maxDist = 0;
        for (const auto& bellPoints : bellPointsList)
        {
            const MVector vec = bellPoints[i] - baseBellPoints[i];
            const double d = vec.length();
            sum += pow(d, 2);

            if (d > maxDist)
                maxDist = d;
        }

        if (sum > 0)
        {
            MVector wp;
            for (const auto& bellPoints : bellPointsList)
            {
                const MVector vec = bellPoints[i] - baseBellPoints[i];
                const double d = pow(vec.length(), 2);
                const double w = d / sum;

                wp += vec * w;
            }

            double wp_len = wp.length();
            if (wp_len > 1e-5)
            {
                outBellPoints[i] += wp.normal() * maxDist;
            }
        }
    }
}

MStatus BellColliderSolver::solve(const BellColliderInputs& inputs, BellColliderOutputs& outputs)
{
    const MMatrix bellMatrix = inputs.bellMatrix;
    const int bellSubdivision = inputs.bellSubdivision;
    const float bellBottomRadius = inputs.bellBottomRadius;

    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1); // Y axis
    const MVector bellNormal = bellAxis.normal();
    const Plane bellPlane(bell_translate, bellNormal);

    MObject bellMesh = makeBellMesh(bellMatrix, 1, bellSubdivision, 1, bellBottomRadius, 1);
    MFnMesh bellMeshFn(bellMesh);

    MPointArray baseBellPoints;
    bellMeshFn.getPoints(baseBellPoints);

    vector<MPointArray> bellPointsList;
    deformPoints(inputs, baseBellPoints, bellPlane, bellPointsList);

    MPointArray outBellPoints;
    averageDisplacements(bellSubdivision, baseBellPoints, bellPointsList, outBellPoints);

    bellMeshFn.setPoints(outBellPoints);

    MObject outCurve = makeBellCurve(outBellPoints, bellSubdivision);

    outputs.outputCurveData = outCurve;
    outputs.outputBellMeshData = bellMesh;

    return MS::kSuccess;
}
