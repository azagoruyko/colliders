#include <maya/MFnPlugin.h>

#include "bellCollider.h"
#include "planeCollider.h"

MStatus initializePlugin(MObject plugin)
{
	MStatus stat;

	MFnPlugin pluginFn(plugin);
	stat = pluginFn.registerNode("bellCollider", BellCollider::typeId, BellCollider::creator, BellCollider::initialize, MPxNode::kLocatorNode, &BellCollider::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::registerDrawOverrideCreator(BellCollider::drawDbClassification, BellCollider::drawRegistrantId, BellColliderDrawOverride::creator);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.registerNode("planeCollider", PlaneCollider::typeId, PlaneCollider::creator, PlaneCollider::initialize, MPxNode::kLocatorNode, &PlaneCollider::drawDbClassification);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::registerDrawOverrideCreator(PlaneCollider::drawDbClassification, PlaneCollider::drawRegistrantId, PlaneColliderDrawOverride::creator);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	return MS::kSuccess;
}

MStatus uninitializePlugin(MObject plugin)
{
	MStatus stat;

	MFnPlugin pluginFn(plugin);

	stat = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator(BellCollider::drawDbClassification, BellCollider::drawRegistrantId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(BellCollider::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator(PlaneCollider::drawDbClassification, PlaneCollider::drawRegistrantId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	stat = pluginFn.deregisterNode(PlaneCollider::typeId);
	CHECK_MSTATUS_AND_RETURN_IT(stat);

	return MS::kSuccess;
}