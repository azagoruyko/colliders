#pragma once

#include <maya/MPxLocatorNode.h>
#include <maya/MPxDrawOverride.h>
#include <maya/MDrawRegistry.h>
#include <maya/MPointArray.h>
#include <maya/MStatus.h>
#include <maya/MTypeId.h>
#include <maya/MMatrix.h>
#include <maya/MColor.h>
#include <maya/MObject.h>
#include <vector>

struct SkirtDrawData
{
    std::vector<MPointArray> bellCurves;
    std::vector<MObject> ringMeshList;
    std::vector<MMatrix> ringMatrices;
    int ringSubdivision = 16;
    MColor color;
};

class SkirtBellCollider : public MPxLocatorNode
{
public:
    static MTypeId typeId;
    static MString typeName;

    static MString drawDbClassification;
    static MString drawRegistrantId;

    // Inputs
    static MObject attr_bellMatrix;
    static MObject attr_leftHipMatrix;
    static MObject attr_leftKneeMatrix;
    static MObject attr_leftHeelMatrix;
    static MObject attr_rightHipMatrix;
    static MObject attr_rightKneeMatrix;
    static MObject attr_rightHeelMatrix;

    static MObject attr_skirtType;
    static MObject attr_height;
    static MObject attr_ringScale;
    static MObject attr_bellScale;
    static MObject attr_bellSubdivision;
    static MObject attr_ringSubdivision;
    static MObject attr_falloff;
    static MObject attr_collision;
    static MObject attr_tightness;
    static MObject attr_bellScaleRamp;
    static MObject attr_leftRingAxis;
    static MObject attr_rightRingAxis;
    static MObject attr_bellAxis;

    // Outputs
    static MObject attr_outputSurface;

    SkirtBellCollider() : MPxLocatorNode() {}
    virtual ~SkirtBellCollider() override {}

    static void* creator() { return new SkirtBellCollider(); }
    static MStatus initialize();
    virtual void postConstructor() override;

    virtual MStatus compute(const MPlug& plug, MDataBlock& dataBlock) override;
};

class SkirtBellColliderDrawData : public MUserData
{
public:
    SkirtBellColliderDrawData() : MUserData() {}
    SkirtDrawData drawData;
};

class SkirtBellColliderDrawOverride : public MHWRender::MPxDrawOverride
{
public:
    static MHWRender::MPxDrawOverride* creator(const MObject& obj) { return new SkirtBellColliderDrawOverride(obj); }

    MHWRender::DrawAPI supportedDrawAPIs() const { return MHWRender::kOpenGL | MHWRender::kDirectX11 | MHWRender::kOpenGLCoreProfile; }

    MUserData* prepareForDraw(const MDagPath& objPath, const MDagPath& cameraPath, const MHWRender::MFrameContext& frameContext, MUserData* oldData) override;

    virtual bool hasUIDrawables() const override { return true; }
    virtual void addUIDrawables(const MDagPath& objPath, MHWRender::MUIDrawManager& drawManager, const MHWRender::MFrameContext& frameContext, const MUserData* data) override;

private:
    SkirtBellColliderDrawOverride(const MObject& obj) : MHWRender::MPxDrawOverride(obj, NULL, true) {} // alwaysDirty is true
};
