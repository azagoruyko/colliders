#include <set>
#include <map>

#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MIntArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnDependencyNode.h>

#include <maya/MTransformationMatrix.h>
#include <maya/MQuaternion.h>
#include <maya/MEulerRotation.h>

#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MMeshIntersector.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MFnMatrixData.h>

#include <maya/MArrayDataBuilder.h>
#include <maya/MArrayDataHandle.h>

#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnEnumAttribute.h>

#include <tbb/parallel_for.h>

#include "bellCollider.h"
#include "utils.hpp"

using namespace std;

MTypeId BellCollider::typeId(1274434);

MObject BellCollider::attr_bellMatrix;
MObject BellCollider::attr_ringMatrix;
MObject BellCollider::attr_bellSubdivision;
MObject BellCollider::attr_ringSubdivision;
MObject BellCollider::attr_bellBottomRadius;
MObject BellCollider::attr_falloff;
MObject BellCollider::attr_collision;
MObject BellCollider::attr_drawColor;
MObject BellCollider::attr_drawOpacity;
MObject BellCollider::attr_outputCurve;
MObject BellCollider::attr_outputBellMesh;

MString BellCollider::drawDbClassification = "drawdb/geometry/bellCollider";
MString BellCollider::drawRegistrantId = "collidersPlugin";

MStatus BellCollider::compute(const MPlug &plug, MDataBlock &dataBlock)
{
    if (plug != attr_outputCurve && plug != attr_outputBellMesh)
        return MS::kUnknownParameter;

    // Extract inputs
    BellColliderInputs inputs;
    inputs.bellMatrix = dataBlock.inputValue(attr_bellMatrix).asMatrix();
    
    auto ringMatrixHandle = dataBlock.inputArrayValue(attr_ringMatrix);
    if (ringMatrixHandle.elementCount() == 0)
        return MS::kFailure;

    for (int i = 0; i < ringMatrixHandle.elementCount(); i++)
    {
        ringMatrixHandle.jumpToElement(i);
        inputs.ringMatrices.push_back(ringMatrixHandle.inputValue().asMatrix());
    }

    inputs.bellSubdivision = dataBlock.inputValue(attr_bellSubdivision).asInt();
    inputs.ringSubdivision = dataBlock.inputValue(attr_ringSubdivision).asInt();
    inputs.bellBottomRadius = dataBlock.inputValue(attr_bellBottomRadius).asFloat();
    inputs.falloff = dataBlock.inputValue(attr_falloff).asFloat();
    inputs.collision = dataBlock.inputValue(attr_collision).asFloat();

    // Call solver
    BellColliderOutputs outputs;
    MStatus stat = BellColliderSolver::solve(inputs, outputs);
    if (stat != MS::kSuccess)
        return stat;

    dataBlock.outputValue(attr_outputCurve).setMObject(outputs.outputCurveData);
    dataBlock.outputValue(attr_outputBellMesh).setMObject(outputs.outputBellMeshData);

    dataBlock.setClean(attr_outputCurve);
    dataBlock.setClean(attr_outputBellMesh);

    return MS::kSuccess;
}

MStatus BellCollider::initialize()
{
    MFnNumericAttribute nAttr;
    MFnMatrixAttribute mAttr;
    MFnEnumAttribute eAttr;
    MFnTypedAttribute tAttr;
    
    attr_bellMatrix = mAttr.create("bellMatrix", "bellMatrix");
    mAttr.setHidden(true);
    addAttribute(attr_bellMatrix);

    attr_ringMatrix = mAttr.create("ringMatrix", "ringMatrix");
    mAttr.setArray(true);
    mAttr.setHidden(true);
    addAttribute(attr_ringMatrix);

    attr_bellSubdivision = nAttr.create("bellSubdivision", "bellSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(3);
    nAttr.setKeyable(true);
    addAttribute(attr_bellSubdivision);

    attr_ringSubdivision = nAttr.create("ringSubdivision", "ringSubdivision", MFnNumericData::kInt, 16);
    nAttr.setMin(3);
    nAttr.setKeyable(true);
    addAttribute(attr_ringSubdivision);

    attr_bellBottomRadius = nAttr.create("bellBottomRadius", "bellBottomRadius", MFnNumericData::kFloat, 0.8);
    nAttr.setMin(0);
    nAttr.setKeyable(true);
    addAttribute(attr_bellBottomRadius);

    attr_falloff = nAttr.create("falloff", "falloff", MFnNumericData::kFloat, 0);
    nAttr.setMin(-1);
    nAttr.setMax(1);
    nAttr.setKeyable(true);
    addAttribute(attr_falloff);

    attr_collision = nAttr.create("collision", "collision", MFnNumericData::kFloat, 0);
    nAttr.setMin(0);
    nAttr.setMax(1);
    nAttr.setKeyable(true);
    addAttribute(attr_collision);

    attr_drawColor = nAttr.create("drawColor", "drawColor", MFnNumericData::k3Double);
    nAttr.setDefault(0.0, 0.01, 0.11);
    nAttr.setMin(0, 0, 0);
    nAttr.setMax(1, 1, 1);
    nAttr.setKeyable(true);
    addAttribute(attr_drawColor);

    attr_drawOpacity = nAttr.create("drawOpacity", "drawOpacity", MFnNumericData::kFloat, 0.3);
    nAttr.setMin(0);
    nAttr.setMax(1);
    nAttr.setKeyable(true);
    addAttribute(attr_drawOpacity);

    attr_outputCurve = tAttr.create("outputCurve", "outputCurve", MFnData::kNurbsCurve);
    tAttr.setHidden(true);
    addAttribute(attr_outputCurve);

    attr_outputBellMesh = tAttr.create("outputBellMesh", "outputBellMesh", MFnData::kMesh);
    tAttr.setHidden(true);
    addAttribute(attr_outputBellMesh);

    attributeAffects(attr_bellMatrix, attr_outputCurve);
    attributeAffects(attr_ringMatrix, attr_outputCurve);
    attributeAffects(attr_bellSubdivision, attr_outputCurve);
    attributeAffects(attr_ringSubdivision, attr_outputCurve);
    attributeAffects(attr_bellBottomRadius, attr_outputCurve);
    attributeAffects(attr_falloff, attr_outputCurve);
    attributeAffects(attr_collision, attr_outputCurve);
    attributeAffects(attr_drawColor, attr_outputCurve);
    attributeAffects(attr_drawOpacity, attr_outputCurve);

    attributeAffects(attr_bellMatrix, attr_outputBellMesh);
    attributeAffects(attr_ringMatrix, attr_outputBellMesh);
    attributeAffects(attr_bellSubdivision, attr_outputBellMesh);
    attributeAffects(attr_ringSubdivision, attr_outputBellMesh);
    attributeAffects(attr_bellBottomRadius, attr_outputBellMesh);
    attributeAffects(attr_falloff, attr_outputBellMesh);
    attributeAffects(attr_collision, attr_outputBellMesh);
    attributeAffects(attr_drawColor, attr_outputBellMesh);
    attributeAffects(attr_drawOpacity, attr_outputBellMesh);

    return MS::kSuccess;
}

void drawCylinder(MHWRender::MUIDrawManager& drawManager, const MObject& mesh)
{
    if (mesh.isNull())
        return;

    const MFnMesh meshFn(mesh);
    const int numSides = (meshFn.numVertices() - 1) / 2;

    MPointArray points;
    meshFn.getPoints(points);

    for (int i = 0; i < numSides; i++)
    {
        drawManager.line(points[i + 1], i == numSides - 1 ? points[1] : points[i + 2]); // buttom
        drawManager.line(points[numSides + i + 1], i == numSides - 1 ? points[numSides + 1] : points[numSides + i + 2]); // top
        drawManager.line(points[i + 1], points[numSides + i + 1]); // edge
    }
}

void drawMesh(MHWRender::MUIDrawManager& drawManager, const MObject& mesh, const MColor &color)
{
    MIntArray triangleCount, triangleIndices;
    MFnMesh meshFn(mesh);

    MPointArray points;
    meshFn.getTriangles(triangleCount, triangleIndices);
    meshFn.getPoints(points);

    //MFloatVectorArray meshNormals;
    //meshFn.getVertexNormals(true, meshNormals);

    MFloatPointArray positions(triangleIndices.length());
    MColorArray colors(triangleIndices.length(), color);
    //MFloatVectorArray normals;

    for (int i = 0; i < triangleIndices.length(); i++)
    {
        positions[i] = points[triangleIndices[i]];
        //normals.append(meshNormals[triangleIndices[i]]);
    }

    drawManager.mesh(MUIDrawManager::kTriangles, positions, NULL, &colors);
}

MUserData* BellColliderDrawOverride::prepareForDraw(
    const MDagPath& objPath, 
    const MDagPath& cameraPath, 
    const MHWRender::MFrameContext& frameContext, 
    MUserData* oldData)
{
    MStatus stat;
    MObject obj = objPath.node(&stat);
    if (stat != MS::kSuccess)
        return NULL;

    auto* data = dynamic_cast<BellColliderDrawData*>(oldData);
    if (!data)
        data = new BellColliderDrawData();

    MPlug bellMatrixPlug(obj, BellCollider::attr_bellMatrix);
    MPlug ringMatrixPlug(obj, BellCollider::attr_ringMatrix);
    MPlug bellSubdivisionPlug(obj, BellCollider::attr_bellSubdivision);
    MPlug ringSubdivisionPlug(obj, BellCollider::attr_ringSubdivision);
    MPlug bellBottomRadiusPlug(obj, BellCollider::attr_bellBottomRadius);
    MPlug drawColorPlug(obj, BellCollider::attr_drawColor);
    MPlug drawOpacityPlug(obj, BellCollider::attr_drawOpacity);
    MPlug outputBellMeshPlug(obj, BellCollider::attr_outputBellMesh);

    MMatrix bellMatrix;
    MObject matrixObj;
    if (bellMatrixPlug.getValue(matrixObj) == MS::kSuccess)
        bellMatrix = MFnMatrixData(matrixObj).matrix();

    const MMatrix bellMatrixInverse = bellMatrix.inverse();

    int bellSubdivision = 16;
    bellSubdivisionPlug.getValue(bellSubdivision);

    int ringSubdivision = 16;
    ringSubdivisionPlug.getValue(ringSubdivision);

    float bellBottomRadius = 0.8f;
    bellBottomRadiusPlug.getValue(bellBottomRadius);

    double rVal = 0.0, gVal = 0.01, bVal = 0.11;
    if (drawColorPlug.numChildren() == 3)
    {
        drawColorPlug.child(0).getValue(rVal);
        drawColorPlug.child(1).getValue(gVal);
        drawColorPlug.child(2).getValue(bVal);
    }

    float drawOpacity = 0.3f;
    drawOpacityPlug.getValue(drawOpacity);

    MObject bellMesh;
    outputBellMeshPlug.getValue(bellMesh);

    data->drawData.color = MColor(rVal, gVal, bVal, drawOpacity);
    data->drawData.bellCenter = taxis(bellMatrix);
    data->drawData.bellMesh = bellMesh;
    data->drawData.ringMeshList.clear();
    data->drawData.collisionPointBellList.clear();
    data->drawData.collisionPointRingList.clear();
    data->drawData.ringDirectionList.clear();
    data->drawData.ringPositionList.clear();

    const MPoint bell_translate = taxis(bellMatrix);
    const MVector bellAxis = maxis(bellMatrix, 1); // Y axis
    const MVector bellNormal = bellAxis.normal();
    const Plane bellPlane(bell_translate, bellNormal);

    unsigned int numRings = ringMatrixPlug.numElements();
    for (unsigned int i = 0; i < numRings; i++)
    {
        MPlug ringMatrixElementPlug = ringMatrixPlug.elementByPhysicalIndex(i);
        MObject ringMatrixObj;
        if (ringMatrixElementPlug.getValue(ringMatrixObj) != MS::kSuccess)
            continue;
        MMatrix ringMatrix = MFnMatrixData(ringMatrixObj).matrix();
        const MMatrix ringMatrixInverse = ringMatrix.inverse();

        const MVector ringDirection = maxis(ringMatrix, 1); // Y axis
        const MPoint ring_translate = taxis(ringMatrix);
        const MVector ringNormal = ringDirection.normal();
        const Plane ringPlane(ring_translate, ringNormal);

        const MPoint ring_translate_proj = bellPlane.projectPoint(taxis(ringMatrix));
        const MVector ringDirection_proj = bellPlane.projectVector(ringDirection);

        data->drawData.ringDirectionList.push_back(ringDirection_proj);
        data->drawData.ringPositionList.push_back(taxis(ringMatrix));

        if (ringDirection_proj.length() > 1e-3)
        {
            const MPointArray hitPoints = findSphereLineIntersection(ring_translate_proj * bellMatrixInverse, ringDirection_proj * bellMatrixInverse, MPoint(0,0,0), 1.001);

            if (hitPoints.length() > 0)
            {
                MPoint collisionPointBell = hitPoints[0] * bellMatrix + bellAxis;

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

                MPoint collisionPointRing;
                for (int k = 0; k < sphereLinePoints.length(); k++)
                {
                    if ((sphereLinePoints[k] - ring_translate) * ringDirection > 0)
                        collisionPointRing = sphereLinePoints[k];
                }

                data->drawData.collisionPointRingList.push_back(linePoint);
                data->drawData.collisionPointRingList.push_back(collisionPointRing);
                data->drawData.collisionPointBellList.push_back(collisionPointBell);
            }
        }

        data->drawData.ringMeshList.push_back(BellColliderSolver::makeBellMesh(ringMatrix, 1, ringSubdivision, 1));
    }

    return data;
}

void BellColliderDrawOverride::addUIDrawables(
    const MDagPath& objPath, 
    MHWRender::MUIDrawManager& drawManager, 
    const MHWRender::MFrameContext& frameContext, 
    const MUserData* data)
{
    auto* bellColliderData = dynamic_cast<const BellColliderDrawData*>(data);

    if (bellColliderData)
    {
        drawManager.beginDrawable();

        const auto& drawData = bellColliderData->drawData;
        if (!drawData.bellMesh.isNull())
        {
            drawMesh(drawManager, drawData.bellMesh, drawData.color);
            drawManager.setColor(MColor(0, 0, 0));
            drawCylinder(drawManager, drawData.bellMesh); // wireframe

            for (const auto& m : drawData.ringMeshList)
            {
                drawMesh(drawManager, m, drawData.color * 0.5);
                drawManager.setColor(MColor(0, 0, 0));
                drawCylinder(drawManager, m); // wireframe
            }

            drawManager.setColor(MColor(0, 0, 0));
            drawManager.setPointSize(5);

            for (const auto& p : drawData.collisionPointBellList)
                drawManager.point(p);

            for (const auto& p : drawData.collisionPointRingList)
                drawManager.point(p);

            for (size_t i = 0; i < drawData.ringDirectionList.size(); i++)
                drawManager.line(drawData.ringPositionList[i], drawData.ringPositionList[i] + drawData.ringDirectionList[i]);
        }

        drawManager.endDrawable();
    }
}