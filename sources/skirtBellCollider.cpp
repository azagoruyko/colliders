#include <maya/MPlug.h>
#include <maya/MPointArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MIntArray.h>
#include <maya/MMatrix.h>
#include <maya/MVector.h>
#include <maya/MPoint.h>
#include <maya/MGlobal.h>
#include <maya/MQuaternion.h>

#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MRampAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MUIDrawManager.h>
#include <maya/MFnMatrixData.h>
#include <maya/MColorArray.h>
#include <maya/MFloatPointArray.h>

#include <vector>
#include <cmath>

#include "skirtBellCollider.h"
#include "bellColliderSolver.h"
#include "utils.hpp"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MTypeId SkirtBellCollider::typeId(1274436);
MString SkirtBellCollider::typeName("skirtBellCollider");

MString SkirtBellCollider::drawDbClassification = "drawdb/geometry/skirtBellCollider";
MString SkirtBellCollider::drawRegistrantId = "collidersPlugin";

MObject SkirtBellCollider::attr_bellMatrix;
MObject SkirtBellCollider::attr_leftHipMatrix;
MObject SkirtBellCollider::attr_leftKneeMatrix;
MObject SkirtBellCollider::attr_leftHeelMatrix;
MObject SkirtBellCollider::attr_rightHipMatrix;
MObject SkirtBellCollider::attr_rightKneeMatrix;
MObject SkirtBellCollider::attr_rightHeelMatrix;

MObject SkirtBellCollider::attr_skirtType;
MObject SkirtBellCollider::attr_height;
MObject SkirtBellCollider::attr_ringScale;
MObject SkirtBellCollider::attr_bellScale;
MObject SkirtBellCollider::attr_bellSubdivision;
MObject SkirtBellCollider::attr_ringSubdivision;
MObject SkirtBellCollider::attr_falloff;
MObject SkirtBellCollider::attr_collision;
MObject SkirtBellCollider::attr_bellScaleRamp;
MObject SkirtBellCollider::attr_leftRingAxis;
MObject SkirtBellCollider::attr_rightRingAxis;
MObject SkirtBellCollider::attr_bellAxis;

MObject SkirtBellCollider::attr_outputSurface;

static MPointArray getCurvePoints(const MPointArray& points, int bellSubdivision, bool use_bottom)
{
    MPointArray cvs;
    int start = use_bottom ? 1 : bellSubdivision + 1;
    int count = bellSubdivision;
    for (int i = 0; i < count; i++)
    {
        cvs.append(points[start + i]);
    }
    // Overlap the first 3 CVs for cubic periodic continuity
    cvs.append(cvs[0]);
    cvs.append(cvs[1]);
    cvs.append(cvs[2]);
    return cvs;
}

static MVector getAxis(const MMatrix& m, short idx)
{
    // 0=X, 1=Y, 2=Z, 3=-X, 4=-Y, 5=-Z
    const MVector v = maxis(m, (unsigned int)(idx % 3));
    return idx >= 3 ? -v : v;
}

static MMatrix createRingMatrix(const MMatrix& jointMatrix, const MVector& scale, short axisIndex = 0, const MPoint* targetPos = nullptr)
{
    const MPoint pos = taxis(jointMatrix);

    // Use the selected joint axis as the ring normal (Y row in the ring matrix).
    const MVector rawY = getAxis(jointMatrix, axisIndex);
    MVector Y = rawY.normal();

    // Extract the remaining joint axes to preserve the joint's actual rotation/twist
    MVector jointX = maxis(jointMatrix, (unsigned int)((axisIndex + 1) % 3));
    MVector uX = jointX.normal();

    // If targetPos is provided, rotate the axes using a quaternion to align Y with the direction vector
    if (targetPos)
    {
        const MVector V = *targetPos - pos;
        if (V.length() > 1e-6)
        {
            MQuaternion Q = Y.rotateTo(V);
            MMatrix rotMatrix = Q.asMatrix();
            Y = Y * rotMatrix;
            uX = uX * rotMatrix;
        }
    }

    // Project uX to be perpendicular to Y
    MVector X = (uX - (uX * Y) * Y).normal();
    MVector Z = (Y ^ X).normal();

    X *= scale.x;
    Y *= scale.y;
    Z *= scale.z;

    double m[4][4] = {
        {X.x, X.y, X.z, 0.0},
        {Y.x, Y.y, Y.z, 0.0},
        {Z.x, Z.y, Z.z, 0.0},
        {pos.x, pos.y, pos.z, 1.0}
    };
    return MMatrix(m);
}

static void drawCylinder(MHWRender::MUIDrawManager& drawManager, const MObject& mesh)
{
    if (mesh.isNull())
        return;

    const MFnMesh meshFn(mesh);
    const int numSides = (meshFn.numVertices() - 1) / 2;

    MPointArray points;
    meshFn.getPoints(points);

    for (int i = 0; i < numSides; i++)
    {
        drawManager.line(points[i + 1], i == numSides - 1 ? points[1] : points[i + 2]); // bottom
        drawManager.line(points[numSides + i + 1], i == numSides - 1 ? points[numSides + 1] : points[numSides + i + 2]); // top
        drawManager.line(points[i + 1], points[numSides + i + 1]); // edge
    }
}

static void drawMesh(MHWRender::MUIDrawManager& drawManager, const MObject& mesh, const MColor &color)
{
    MIntArray triangleCount, triangleIndices;
    MFnMesh meshFn(mesh);

    MPointArray points;
    meshFn.getTriangles(triangleCount, triangleIndices);
    meshFn.getPoints(points);

    MFloatPointArray positions(triangleIndices.length());
    MColorArray colors(triangleIndices.length(), color);

    for (int i = 0; i < triangleIndices.length(); i++)
    {
        positions[i] = points[triangleIndices[i]];
    }

    drawManager.mesh(MUIDrawManager::kTriangles, positions, NULL, &colors);
}

void SkirtBellCollider::postConstructor()
{
    MObject thisObj = thisMObject();
    MRampAttribute rampAttr(thisObj, attr_bellScaleRamp);

    MFloatArray positions;
    MFloatArray values;
    MIntArray interpolations;

    positions.append(0.0f);
    values.append(1.0f);
    interpolations.append(MRampAttribute::kLinear);

    positions.append(1.0f);
    values.append(1.0f);
    interpolations.append(MRampAttribute::kLinear);

    rampAttr.addEntries(positions, values, interpolations);
}

MStatus SkirtBellCollider::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnTypedAttribute tAttr;
    MStatus stat;

    // Joint matrices
    attr_bellMatrix = mAttr.create("bellMatrix", "bellMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_bellMatrix);

    attr_leftHipMatrix = mAttr.create("leftHipMatrix", "leftHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftHipMatrix);

    attr_leftKneeMatrix = mAttr.create("leftKneeMatrix", "leftKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftKneeMatrix);

    attr_leftHeelMatrix = mAttr.create("leftHeelMatrix", "leftHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_leftHeelMatrix);

    attr_rightHipMatrix = mAttr.create("rightHipMatrix", "rightHipMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightHipMatrix);

    attr_rightKneeMatrix = mAttr.create("rightKneeMatrix", "rightKneeMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightKneeMatrix);

    attr_rightHeelMatrix = mAttr.create("rightHeelMatrix", "rightHeelMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_rightHeelMatrix);

    // Skirt Type
    attr_skirtType = eAttr.create("skirtType", "skirtType", 1);
    eAttr.addField("Short", 0);
    eAttr.addField("Long", 1);
    eAttr.setKeyable(true);
    addAttribute(attr_skirtType);


    // Height
    attr_height = nAttr.create("height", "height", MFnNumericData::kFloat, 1.0f);
    nAttr.setMin(0.01f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_height);

    // Ring Scale
    attr_ringScale = nAttr.create("ringScale", "ringScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.5, 1.0, 0.5);
    nAttr.setKeyable(true);
    addAttribute(attr_ringScale);

    // Bell Scale
    attr_bellScale = nAttr.create("bellScale", "bellScale", MFnNumericData::k3Double);
    nAttr.setDefault(0.8, 1.0, 0.8);
    nAttr.setKeyable(true);
    addAttribute(attr_bellScale);

    // Bell Subdivision
    attr_bellSubdivision = nAttr.create("bellSubdivision", "bellSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(3);
    nAttr.setKeyable(true);
    addAttribute(attr_bellSubdivision);

    // Ring Subdivision
    attr_ringSubdivision = nAttr.create("ringSubdivision", "ringSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(3);
    nAttr.setKeyable(true);
    addAttribute(attr_ringSubdivision);

    // Falloff
    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0.0f);
    nAttr.setMin(-1.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    // Collision
    attr_collision = nAttr.create("collision", "collision", MFnNumericData::kFloat, 1.0f);
    nAttr.setMin(0.0f);
    nAttr.setMax(1.0f);
    nAttr.setKeyable(true);
    addAttribute(attr_collision);

    // Bell Scale Ramp (Curve)
    attr_bellScaleRamp = MRampAttribute::createCurveRamp("bellScaleRamp", "bellScaleRamp");
    addAttribute(attr_bellScaleRamp);

    // Left Leg Axis
    attr_leftRingAxis = eAttr.create("leftRingAxis", "leftRingAxis", 0);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_leftRingAxis);

    // Right Leg Axis
    attr_rightRingAxis = eAttr.create("rightRingAxis", "rightRingAxis", 3);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_rightRingAxis);

    // Skirt Axis
    attr_bellAxis = eAttr.create("bellAxis", "bellAxis", 1);
    eAttr.addField("X", 0);
    eAttr.addField("Y", 1);
    eAttr.addField("Z", 2);
    eAttr.addField("-X", 3);
    eAttr.addField("-Y", 4);
    eAttr.addField("-Z", 5);
    eAttr.setKeyable(true);
    addAttribute(attr_bellAxis);

    // Output NURBS Surface (Visible to users)
    attr_outputSurface = tAttr.create("outputSurface", "outputSurface", MFnData::kNurbsSurface);
    tAttr.setWritable(false);
    tAttr.setStorable(false);
    addAttribute(attr_outputSurface);

    // Set up attribute affects relationships
    attributeAffects(attr_bellMatrix, attr_outputSurface);
    attributeAffects(attr_leftHipMatrix, attr_outputSurface);
    attributeAffects(attr_leftKneeMatrix, attr_outputSurface);
    attributeAffects(attr_leftHeelMatrix, attr_outputSurface);
    attributeAffects(attr_rightHipMatrix, attr_outputSurface);
    attributeAffects(attr_rightKneeMatrix, attr_outputSurface);
    attributeAffects(attr_rightHeelMatrix, attr_outputSurface);
    attributeAffects(attr_skirtType, attr_outputSurface);
    attributeAffects(attr_height, attr_outputSurface);
    attributeAffects(attr_ringScale, attr_outputSurface);
    attributeAffects(attr_bellScale, attr_outputSurface);
    attributeAffects(attr_bellSubdivision, attr_outputSurface);
    attributeAffects(attr_ringSubdivision, attr_outputSurface);
    attributeAffects(attr_falloff, attr_outputSurface);
    attributeAffects(attr_collision, attr_outputSurface);
    attributeAffects(attr_bellScaleRamp, attr_outputSurface);
    attributeAffects(attr_leftRingAxis, attr_outputSurface);
    attributeAffects(attr_rightRingAxis, attr_outputSurface);
    attributeAffects(attr_bellAxis, attr_outputSurface);

    return MS::kSuccess;
}

MStatus SkirtBellCollider::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    if (plug != attr_outputSurface)
        return MS::kUnknownParameter;

    MStatus stat;

    // Get input values
    const MMatrix inputBellMatrix = dataBlock.inputValue(attr_bellMatrix).asMatrix();
    const MMatrix leftHipMatrix = dataBlock.inputValue(attr_leftHipMatrix).asMatrix();
    const MMatrix leftKneeMatrix = dataBlock.inputValue(attr_leftKneeMatrix).asMatrix();
    const MMatrix leftHeelMatrix = dataBlock.inputValue(attr_leftHeelMatrix).asMatrix();
    const MMatrix rightHipMatrix = dataBlock.inputValue(attr_rightHipMatrix).asMatrix();
    const MMatrix rightKneeMatrix = dataBlock.inputValue(attr_rightKneeMatrix).asMatrix();
    const MMatrix rightHeelMatrix = dataBlock.inputValue(attr_rightHeelMatrix).asMatrix();

    const short skirtType = dataBlock.inputValue(attr_skirtType).asShort();
    const float height = dataBlock.inputValue(attr_height).asFloat();
    const MVector ringScale = dataBlock.inputValue(attr_ringScale).asVector();
    const MVector bellScale = dataBlock.inputValue(attr_bellScale).asVector();
    const int bellSubdivision = dataBlock.inputValue(attr_bellSubdivision).asInt();
    const int ringSubdivision = dataBlock.inputValue(attr_ringSubdivision).asInt();
    const float falloff = dataBlock.inputValue(attr_falloff).asFloat();
    const float collision = dataBlock.inputValue(attr_collision).asFloat();
    const short leftRingAxis = dataBlock.inputValue(attr_leftRingAxis).asShort();
    const short rightRingAxis = dataBlock.inputValue(attr_rightRingAxis).asShort();
    const short bellAxis = dataBlock.inputValue(attr_bellAxis).asShort();

    // Ramp attribute
    MRampAttribute rampAttr(thisMObject(), attr_bellScaleRamp);

    // Midpoints
    const MPoint LH = taxis(leftHipMatrix);
    const MPoint LK = taxis(leftKneeMatrix);
    const MPoint LHe = taxis(leftHeelMatrix);
    const MPoint RH = taxis(rightHipMatrix);
    const MPoint RK = taxis(rightKneeMatrix);
    const MPoint RHe = taxis(rightHeelMatrix);

    const MPoint H = (LH + RH) * 0.5;
    const MPoint K = (LK + RK) * 0.5;
    const MPoint He = (LHe + RHe) * 0.5;

    const MPoint W = taxis(inputBellMatrix);

    // Skirt height calculation using rigid bone lengths to prevent stretching during knee bends
    const double L_thigh = ((LK - LH).length() + (RK - RH).length()) * 0.5;
    const double L_calf = ((LHe - LK).length() + (RHe - RK).length()) * 0.5;

    const double d_hip = (H - W).length();
    const double d_knee = d_hip + L_thigh;
    const double d_heel = d_knee + L_calf;
    const double d_mid = (d_hip + d_knee) * 0.5;

    double h_param = (double)height;
    if (h_param < 0.0) h_param = 0.0;
    if (h_param > 1.0) h_param = 1.0;

    double h_val = 0.0;
    if (skirtType == 1) // Long skirt: bottom interpolates from Knee to Heel
    {
        h_val = d_knee + (d_heel - d_knee) * h_param;
    }
    else // Short skirt: bottom interpolates from Hip to Knee
    {
        h_val = d_hip + (d_knee - d_hip) * h_param;
    }

    const double raw_defaultHeight = (skirtType == 1) ? d_heel : d_knee;
    const double defaultHeight = raw_defaultHeight < 1e-5 ? 1.0 : raw_defaultHeight;
    const double s = h_val / defaultHeight;

    // Get direction vector based on the selected bell matrix axis
    const MVector raw_dir_vector = getAxis(inputBellMatrix, bellAxis);
    const MVector dir_vector = raw_dir_vector.length() < 1e-4 ? MVector(0, -1, 0) : raw_dir_vector.normal();

    // Start point of the skirt aligned rigidly along the chosen waist axis
    const MPoint P_start = W;

    const int N = (skirtType == 0) ? 2 : 3;

    // Construct collision rings
    const MVector thighScale(ringScale.x, L_thigh * ringScale.y, ringScale.z);
    const MMatrix leftHipToKnee = createRingMatrix(leftHipMatrix, thighScale, leftRingAxis, &LK);
    const MMatrix rightHipToKnee = createRingMatrix(rightHipMatrix, thighScale, rightRingAxis, &RK);

    MMatrix leftHipToHeel;
    MMatrix rightHipToHeel;
    if (skirtType == 1)
    {
        const double leftLegLen = L_thigh + L_calf;
        const double rightLegLen = L_thigh + L_calf;
        const MVector leftHeelScale(ringScale.x, leftLegLen * ringScale.y, ringScale.z);
        const MVector rightHeelScale(ringScale.x, rightLegLen * ringScale.y, ringScale.z);
        leftHipToHeel = createRingMatrix(leftHipMatrix, leftHeelScale, leftRingAxis, &LHe);
        rightHipToHeel = createRingMatrix(rightHipMatrix, rightHeelScale, rightRingAxis, &RHe);
    }

    // Setup level distances along the skirt axis
    vector<double> levelDistances;
    levelDistances.push_back(0.0);
    if (skirtType == 1) // Long skirt: height controls the heel level only
    {
        levelDistances.push_back(d_mid);
        levelDistances.push_back(d_knee);
        levelDistances.push_back(h_val);
    }
    else // Short skirt: height scales all levels
    {
        levelDistances.push_back(d_mid * s);
        levelDistances.push_back(d_knee * s);
    }

    vector<MPointArray> rows(N + 1);
    MPointArray controlPoints;
    MDoubleArray uKnots;
    MDoubleArray vKnots;

    // Setup uKnots: uniform for kPeriodic form (degree 3)
    int numKnotsInU = bellSubdivision + 5;
    for (int j = 0; j < numKnotsInU; j++)
    {
        uKnots.append((double)(j - 2));
    }

    const double h_safe = h_val < 1e-5 ? 1e-5 : h_val;

    // Solve for each bell
    for (int i = 0; i < N; i++)
    {
        const double dist = levelDistances[i+1] - levelDistances[i];
        const double safeDist = dist < 1e-4 ? 1e-4 : dist;
        const MPoint P_bell = P_start + dir_vector * levelDistances[i];

        // Build bellMatrix orientation (bulletproof against zero/unconnected inputs)
        const MVector dir_y = dir_vector;

        const MVector raw_waist_X = xaxis(inputBellMatrix);
        const MVector waist_X = raw_waist_X.length() < 1e-4 ? MVector(1, 0, 0) : raw_waist_X.normal();

        const MVector raw_X = (waist_X - (waist_X * dir_y) * dir_y);
        const MVector X = [&]() {
            if (raw_X.length() < 1e-4) {
                const MVector raw_waist_Z = zaxis(inputBellMatrix);
                const MVector waist_Z = raw_waist_Z.length() < 1e-4 ? MVector(0, 0, 1) : raw_waist_Z.normal();
                const MVector crossed = dir_y ^ waist_Z;
                if (crossed.length() < 1e-4) {
                    return MVector(1, 0, 0);
                }
                return crossed.normal();
            }
            return raw_X.normal();
        }();

        const MVector raw_Z = X ^ dir_y;
        const MVector Z = raw_Z.length() < 1e-4 ? MVector(0, 0, 1) : raw_Z.normal();

        // Query Ramp values using actual distance ratios
        const float t_bottom = (float)(levelDistances[i] / h_safe);
        const float t_top = (float)(levelDistances[i+1] / h_safe);
        float raw_scale_bottom = 1.0f;
        float raw_scale_top = 1.0f;
        rampAttr.getValueAtPosition(t_bottom, raw_scale_bottom);
        rampAttr.getValueAtPosition(t_top, raw_scale_top);

        const float scale_bottom = raw_scale_bottom;
        const float scale_top = raw_scale_top < 1e-5f ? 1e-5f : raw_scale_top;

        // Scale axes using bellScale attributes
        const MVector scaled_X = X * (bellScale.x * scale_top);
        const MVector Y = dir_y * (safeDist * bellScale.y);
        const MVector scaled_Z = Z * (bellScale.z * scale_top);

        const double m[4][4] = {
            {scaled_X.x, scaled_X.y, scaled_X.z, 0.0},
            {Y.x, Y.y, Y.z, 0.0},
            {scaled_Z.x, scaled_Z.y, scaled_Z.z, 0.0},
            {P_bell.x, P_bell.y, P_bell.z, 1.0}
        };
        const MMatrix bellMatrix(m);

        vector<MMatrix> bellRings;
        if (i == 0 || i == 1)
        {
            bellRings.push_back(leftHipToKnee);
            bellRings.push_back(rightHipToKnee);
            if (skirtType == 1)
            {
                bellRings.push_back(leftHipToHeel);
                bellRings.push_back(rightHipToHeel);
            }
        }
        else if (i == 2 && skirtType == 1)
        {
            bellRings.push_back(leftHipToHeel);
            bellRings.push_back(rightHipToHeel);
        }

        // Setup solver inputs
        BellColliderInputs inputs;
        inputs.bellMatrix = bellMatrix;
        inputs.ringMatrices = bellRings;
        inputs.bellSubdivision = bellSubdivision;
        inputs.ringSubdivision = ringSubdivision;
        inputs.bellBottomRadius = scale_bottom / scale_top;
        inputs.falloff = falloff;
        inputs.collision = collision;

        // Solve
        BellColliderOutputs outputs;
        const MStatus solveStat = BellColliderSolver::solve(inputs, outputs);
        if (solveStat != MS::kSuccess)
        {
            MGlobal::displayError("SkirtBellCollider: Solver failed at bell index " + MString(to_string(i).c_str()) + " with: " + solveStat.errorString());
            return solveStat;
        }

        // Extract points from output mesh
        MPointArray meshPoints;
        const MFnMesh meshFn(outputs.outputBellMeshData);
        meshFn.getPoints(meshPoints);

        // For the first bell, extract the bottom ring points
        if (i == 0)
        {
            rows[0] = getCurvePoints(meshPoints, bellSubdivision, true);
        }

        // Extract top ring points
        rows[i + 1] = getCurvePoints(meshPoints, bellSubdivision, false);
    }

    // Populate vKnots matching actual distance levels
    for (int i = 0; i <= N; i++)
    {
        vKnots.append(levelDistances[i]);
    }

    // Populate controlPoints using the transposed indexing layout required by Maya:
    // index = (numCVsInV * uIndex) + vIndex
    int numU = bellSubdivision + 3;
    int numV = N + 1;
    controlPoints.setLength(numU * numV);
    for (int u = 0; u < numU; u++)
    {
        for (int v = 0; v < numV; v++)
        {
            controlPoints.set(rows[v][numU - 1 - u], u * numV + v);
        }
    }



    // Create NURBS surface
    MFnNurbsSurfaceData surfaceDataFn;
    MObject surfaceData = surfaceDataFn.create(&stat);
    if (stat != MS::kSuccess)
    {
        MGlobal::displayError("SkirtBellCollider: Failed to create NURBS surface data object: " + stat.errorString());
        return stat;
    }

    MFnNurbsSurface surfaceFn;
    surfaceFn.create(
        controlPoints,
        uKnots,
        vKnots,
        3, // uDegree (cubic periodic)
        1, // vDegree (linear open)
        MFnNurbsSurface::kPeriodic, // uForm
        MFnNurbsSurface::kOpen, // vForm
        false, // createRational
        surfaceData,
        &stat
    );
    if (stat != MS::kSuccess)
    {
        MGlobal::displayError("SkirtBellCollider: MFnNurbsSurface::create failed: " + stat.errorString());
        return stat;
    }

    // Set to output plug
    MDataHandle outputHandle = dataBlock.outputValue(attr_outputSurface, &stat);
    if (stat != MS::kSuccess)
        return stat;

    outputHandle.setMObject(surfaceData);
    dataBlock.setClean(plug);

    return MS::kSuccess;
}

MUserData* SkirtBellColliderDrawOverride::prepareForDraw(
    const MDagPath& objPath,
    const MDagPath& cameraPath,
    const MHWRender::MFrameContext& frameContext,
    MUserData* oldData)
{
    MStatus stat;
    MObject obj = objPath.node(&stat);
    if (stat != MS::kSuccess)
        return NULL;

    auto* data = dynamic_cast<SkirtBellColliderDrawData*>(oldData);
    if (!data)
        data = new SkirtBellColliderDrawData();

    // Clear old draw data
    data->drawData.bellCurves.clear();
    data->drawData.ringMeshList.clear();
    data->drawData.ringMatrices.clear();

    // Extract attributes
    MPlug leftHipMatrixPlug(obj, SkirtBellCollider::attr_leftHipMatrix);
    MPlug leftKneeMatrixPlug(obj, SkirtBellCollider::attr_leftKneeMatrix);
    MPlug leftHeelMatrixPlug(obj, SkirtBellCollider::attr_leftHeelMatrix);
    MPlug rightHipMatrixPlug(obj, SkirtBellCollider::attr_rightHipMatrix);
    MPlug rightKneeMatrixPlug(obj, SkirtBellCollider::attr_rightKneeMatrix);
    MPlug rightHeelMatrixPlug(obj, SkirtBellCollider::attr_rightHeelMatrix);
    MPlug skirtTypePlug(obj, SkirtBellCollider::attr_skirtType);
    MPlug ringScalePlug(obj, SkirtBellCollider::attr_ringScale);
    MPlug ringSubdivisionPlug(obj, SkirtBellCollider::attr_ringSubdivision);
    MPlug leftRingAxisPlug(obj, SkirtBellCollider::attr_leftRingAxis);
    MPlug rightRingAxisPlug(obj, SkirtBellCollider::attr_rightRingAxis);

    MMatrix leftHipMatrix, leftKneeMatrix, leftHeelMatrix;
    MMatrix rightHipMatrix, rightKneeMatrix, rightHeelMatrix;
    MObject tempObj;

    if (leftHipMatrixPlug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
        leftHipMatrix = MFnMatrixData(tempObj).matrix();
    if (leftKneeMatrixPlug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
        leftKneeMatrix = MFnMatrixData(tempObj).matrix();
    if (leftHeelMatrixPlug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
        leftHeelMatrix = MFnMatrixData(tempObj).matrix();

    if (rightHipMatrixPlug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
        rightHipMatrix = MFnMatrixData(tempObj).matrix();
    if (rightKneeMatrixPlug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
        rightKneeMatrix = MFnMatrixData(tempObj).matrix();
    if (rightHeelMatrixPlug.getValue(tempObj) == MS::kSuccess && !tempObj.isNull())
        rightHeelMatrix = MFnMatrixData(tempObj).matrix();

    short skirtType = 1;
    skirtTypePlug.getValue(skirtType);

    MVector ringScale(0.5, 1.0, 0.5);
    if (ringScalePlug.numChildren() == 3)
    {
        double sx = 0.5, sy = 1.0, sz = 0.5;
        ringScalePlug.child(0).getValue(sx);
        ringScalePlug.child(1).getValue(sy);
        ringScalePlug.child(2).getValue(sz);
        ringScale = MVector(sx, sy, sz);
    }

    int ringSubdivision = 16;
    ringSubdivisionPlug.getValue(ringSubdivision);
    if (ringSubdivision < 3) ringSubdivision = 3;
    data->drawData.ringSubdivision = ringSubdivision;

    short leftRingAxis = 0;
    leftRingAxisPlug.getValue(leftRingAxis);
    short rightRingAxis = 0;
    rightRingAxisPlug.getValue(rightRingAxis);

    MPoint LH = taxis(leftHipMatrix);
    MPoint LK = taxis(leftKneeMatrix);
    MPoint LHe = taxis(leftHeelMatrix);
    MPoint RH = taxis(rightHipMatrix);
    MPoint RK = taxis(rightKneeMatrix);
    MPoint RHe = taxis(rightHeelMatrix);

    // Construct ring matrices exactly like in compute()
    double leftThigh = (LK - LH).length();
    double rightThigh = (RK - RH).length();

    const MVector leftThighScale(ringScale.x, leftThigh * ringScale.y, ringScale.z);
    const MVector rightThighScale(ringScale.x, rightThigh * ringScale.y, ringScale.z);

    MMatrix leftHipToKnee = createRingMatrix(leftHipMatrix, leftThighScale, leftRingAxis, &LK);
    MMatrix rightHipToKnee = createRingMatrix(rightHipMatrix, rightThighScale, rightRingAxis, &RK);

    data->drawData.ringMatrices.push_back(leftHipToKnee);
    data->drawData.ringMatrices.push_back(rightHipToKnee);

    if (skirtType == 1)
    {
        double leftLegLen = leftThigh + (LHe - LK).length();
        double rightLegLen = rightThigh + (RHe - RK).length();
        const MVector leftHeelScale(ringScale.x, leftLegLen * ringScale.y, ringScale.z);
        const MVector rightHeelScale(ringScale.x, rightLegLen * ringScale.y, ringScale.z);
        MMatrix leftHipToHeel = createRingMatrix(leftHipMatrix, leftHeelScale, leftRingAxis, &LHe);
        MMatrix rightHipToHeel = createRingMatrix(rightHipMatrix, rightHeelScale, rightRingAxis, &RHe);

        data->drawData.ringMatrices.push_back(leftHipToHeel);
        data->drawData.ringMatrices.push_back(rightHipToHeel);
    }

    // Build the meshes for the collision rings using the solver
    for (const auto& matrix : data->drawData.ringMatrices)
    {
        data->drawData.ringMeshList.push_back(BellColliderSolver::makeBellMesh(matrix, 1, ringSubdivision, 1));
    }

    // Extract the deformed bell curves from the output NURBS surface
    MPlug outputSurfacePlug(obj, SkirtBellCollider::attr_outputSurface);
    MObject surfaceData;
    if (outputSurfacePlug.getValue(surfaceData) == MS::kSuccess && !surfaceData.isNull())
    {
        MFnNurbsSurface surfaceFn(surfaceData, &stat);
        if (stat == MS::kSuccess)
        {
            unsigned int numU = surfaceFn.numCVsInU();
            unsigned int numV = surfaceFn.numCVsInV();
            MPointArray cvs;
            if (surfaceFn.getCVs(cvs, MSpace::kObject) == MS::kSuccess && cvs.length() == numU * numV)
            {
                // Each index v in 0..numV-1 represents one bell curve ring along the V direction.
                for (unsigned int v = 0; v < numV; v++)
                {
                    MPointArray curveLoop;
                    for (unsigned int u = 0; u < numU; u++)
                    {
                        unsigned int idx = u * numV + v;
                        if (idx < cvs.length())
                        {
                            curveLoop.append(cvs[idx]);
                        }
                    }
                    if (curveLoop.length() > 0)
                    {
                        data->drawData.bellCurves.push_back(curveLoop);
                    }
                }
            }
        }
    }

    // Default transparent cyan color for drawing collider rings
    data->drawData.color = MColor(0.0f, 0.6f, 1.0f, 0.25f);

    return data;
}

void SkirtBellColliderDrawOverride::addUIDrawables(
    const MDagPath& objPath,
    MHWRender::MUIDrawManager& drawManager,
    const MHWRender::MFrameContext& frameContext,
    const MUserData* data)
{
    auto* skirtDrawData = dynamic_cast<const SkirtBellColliderDrawData*>(data);
    if (!skirtDrawData)
        return;

    const auto& drawData = skirtDrawData->drawData;

    drawManager.beginDrawable();

    // 1. Draw collision rings (shaded transparent cylinders and dark wireframe outlines)
    for (const auto& ringMesh : drawData.ringMeshList)
    {
        if (!ringMesh.isNull())
        {
            drawMesh(drawManager, ringMesh, drawData.color);
            drawManager.setColor(MColor(0.0f, 0.1f, 0.2f, 1.0f));
            drawCylinder(drawManager, ringMesh);
        }
    }

    // 2. Draw all bell curves
    // Using a vibrant yellow/orange to make them clearly visible and premium looking
    drawManager.setColor(MColor(1.0f, 0.75f, 0.0f, 1.0f));
    drawManager.setLineWidth(2.0f);

    for (const auto& curve : drawData.bellCurves)
    {
        if (curve.length() > 1)
        {
            for (unsigned int i = 0; i < curve.length() - 1; i++)
            {
                drawManager.line(curve[i], curve[i + 1]);
            }
        }
    }

    drawManager.endDrawable();
}
