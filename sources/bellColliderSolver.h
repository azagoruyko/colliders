#pragma once

#include <maya/MMatrix.h>
#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>
#include <maya/MObject.h>
#include <maya/MStatus.h>
#include <vector>

#include "utils.hpp"

struct BellColliderInputs
{
    MMatrix bellMatrix;
    std::vector<MMatrix> ringMatrices;
    int bellSubdivision = 16;
    int ringSubdivision = 16;
    float bellBottomRadius = 0.8f;
    float falloff = 0.0f;
    float collision = 0.0f;
};

struct BellColliderOutputs
{
    MObject outputCurveData;       // NurbsCurve data object
    MObject outputBellMeshData;    // Deformed bell mesh data object
};

class BellColliderSolver
{
public:
    static MObject makeBellMesh(const MMatrix& matrix, unsigned int axis, unsigned int numSides, double height = 1, double bottomRadius = 1, double topRadius = 1);
    static MStatus solve(const BellColliderInputs& inputs, BellColliderOutputs& outputs);

private:
    static void deformPoints(const BellColliderInputs& inputs, const MPointArray& baseBellPoints, const Plane& bellPlane, std::vector<MPointArray>& bellPointsList);
    static void averageDisplacements(int bellSubdivision, const MPointArray& baseBellPoints, const std::vector<MPointArray>& bellPointsList, MPointArray& outBellPoints);
};
