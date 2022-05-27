#include <maya/MPxLocatorNode.h>
#include <maya/MPxDrawOverride.h>
#include <maya/MDrawRegistry.h>
#include <maya/MPointArray.h>

#include <vector>
#include <string>

using namespace std;

struct DrawData
{
	MObject bellMesh;
	vector<MObject> ringMeshList;
	vector<MPoint> collisionPointBellList;
	vector<MPoint> collisionPointRingList;
	vector<MVector> ringDirectionList;
	vector<MPoint> ringPositionList;
	MPoint bellCenter;
	MColor color;
};

class BellCollider : public MPxLocatorNode
{
public:
	static MTypeId typeId;

	static MString drawDbClassification;
	static MString drawRegistrantId;

	static MObject attr_bellMatrix;
	static MObject attr_ringMatrix;
	static MObject attr_bellSubdivision;
	static MObject attr_ringSubdivision;
	static MObject attr_bellBottomRadius;
	static MObject attr_falloff;
	static MObject attr_positionCount;
	static MObject attr_drawColor;
	static MObject attr_drawOpacity;
	static MObject attr_outputPositions;
	static MObject attr_outputRotations;
	static MObject attr_outputCurve;

	DrawData drawData;

	static void* creator() { return new BellCollider(); }
	static MStatus initialize();

	MStatus compute(const MPlug&, MDataBlock&);

	void drawUI(MHWRender::MUIDrawManager&);

private:
};

class BellColliderDrawData : public MUserData
{
public:
	BellColliderDrawData() : MUserData(false) {} // deleteAfterUse = false
	BellCollider* bellCollider{nullptr};
};

class BellColliderDrawOverride : public MHWRender::MPxDrawOverride
{
public:
	static MHWRender::MPxDrawOverride* creator(const MObject& obj) { return new BellColliderDrawOverride(obj); }

	MHWRender::DrawAPI supportedDrawAPIs() const { return MHWRender::kOpenGL | MHWRender::kDirectX11 | MHWRender::kOpenGLCoreProfile; }

	MUserData* prepareForDraw(const MDagPath& objPath, const MDagPath& cameraPath, const MHWRender::MFrameContext& frameContext, MUserData* oldData);

	virtual bool hasUIDrawables() const { return true; }
	virtual void addUIDrawables(const MDagPath& objPath, MHWRender::MUIDrawManager& drawManager, const MHWRender::MFrameContext& frameContext, const MUserData* data);

	//bool isBounded(const MDagPath &objPath, const MDagPath &cameraPath) const	{return true;}
	//MBoundingBox boundingBox(const MDagPath &objPath, const MDagPath &cameraPath) const {return MBoundingBox();}

private:
	BellColliderDrawOverride(const MObject& obj) : MHWRender::MPxDrawOverride(obj, NULL, true) {} // alwaysDirty is true
};