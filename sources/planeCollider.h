#include <maya/MPxLocatorNode.h>
#include <maya/MPxDrawOverride.h>
#include <maya/MDrawRegistry.h>
#include <maya/MPointArray.h>

#include <vector>
#include <string>

using namespace std;

struct PlaneDrawData
{
	MPoint planeCenter;
	MVector planeNormal;
	MColor color;
	double size;
};

class PlaneCollider : public MPxLocatorNode
{
public:
	static MTypeId typeId;

	static MString drawDbClassification;
	static MString drawRegistrantId;

	static MObject attr_planeMatrix;
	static MObject attr_normalAxis;
	static MObject attr_inputPosition;
	static MObject attr_drawColor;
	static MObject attr_drawOpacity;
	static MObject attr_outputPosition;

	static void* creator() { return new PlaneCollider(); }
	static MStatus initialize();

	MStatus compute(const MPlug&, MDataBlock&);

private:
};

class PlaneColliderDrawData : public MUserData
{
public:
	PlaneColliderDrawData() : MUserData() {}
	PlaneDrawData drawData;
};

class PlaneColliderDrawOverride : public MHWRender::MPxDrawOverride
{
public:
	static MHWRender::MPxDrawOverride* creator(const MObject& obj) { return new PlaneColliderDrawOverride(obj); }

	MHWRender::DrawAPI supportedDrawAPIs() const { return MHWRender::kOpenGL | MHWRender::kDirectX11 | MHWRender::kOpenGLCoreProfile; }

	MUserData* prepareForDraw(const MDagPath& objPath, const MDagPath& cameraPath, const MHWRender::MFrameContext& frameContext, MUserData* oldData);

	virtual bool hasUIDrawables() const { return true; }
	virtual void addUIDrawables(const MDagPath& objPath, MHWRender::MUIDrawManager& drawManager, const MHWRender::MFrameContext& frameContext, const MUserData* data);

private:
	PlaneColliderDrawOverride(const MObject& obj) : MHWRender::MPxDrawOverride(obj, NULL, true) {} // alwaysDirty is true
};