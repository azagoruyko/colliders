import math
import maya.cmds as cmds

WINDOW_ID = "bellColliderUI"
WINDOW_TITLE = "Bell Collider Creator"

BELL_NODE_TYPE = "bellCollider"
SKIRT_NODE_TYPE = "skirtBellCollider"


# ---------------------------------------------------------------------------
# Node-creation backend
# ---------------------------------------------------------------------------

def createBellCollider(numRings=1, prefix=""):
    """
    Auto-create a bell locator + numRings ring locators, wire them
    into a new bellCollider node. All parameters use node defaults.
    """
    bellLoc = cmds.spaceLocator(name=f"{prefix}skirt_bell_locator")[0]

    ringLocs = []
    for i in range(numRings):
        loc = cmds.spaceLocator(name=f"{prefix}skirt_ring_locator{i + 1}")[0]
        cmds.setAttr(loc + ".s", 0.4, 1.5, 0.4)
        ringLocs.append(loc)

    node = cmds.createNode(BELL_NODE_TYPE, name=f"{prefix}bellCollider")
    nodeXform = cmds.listRelatives(node, parent=True, type="transform")[0]
    cmds.rename(nodeXform, f"{prefix}bellCollider_transform")
    cmds.connectAttr(f"{bellLoc}.worldMatrix[0]", f"{node}.bellMatrix", force=True)
    for idx, rl in enumerate(ringLocs):
        cmds.connectAttr(f"{rl}.worldMatrix[0]", f"{node}.ringMatrix[{idx}]", force=True)

    # Lock t/r/s on the bellCollider's parent transform
    nodeXform = cmds.listRelatives(node, parent=True)[0]
    for attr in ("tx", "ty", "tz", "rx", "ry", "rz", "sx", "sy", "sz"):
        cmds.setAttr(f"{nodeXform}.{attr}", lock=True)

    # Create a curve and connect outputCurve so it's visible in the viewport
    curveXform = cmds.createNode("transform", name=f"{prefix}bellCurve")
    curveShape = cmds.createNode("nurbsCurve", name=f"{prefix}bellCurveShape", parent=curveXform)
    cmds.connectAttr(f"{node}.outputCurve", f"{curveShape}.create", force=True)

    cmds.select(node)
    return node, bellLoc, ringLocs, curveXform


def _findAxis(startObj, endObj):
    if not startObj or not endObj:
        return 1
    pwm = cmds.xform(startObj, query=True, worldSpace=True, matrix=True)
    p1 = cmds.xform(startObj, query=True, worldSpace=True, translation=True)
    p2 = cmds.xform(endObj, query=True, worldSpace=True, translation=True)
    
    vx = p2[0] - p1[0]
    vy = p2[1] - p1[1]
    vz = p2[2] - p1[2]
    l = math.sqrt(vx*vx + vy*vy + vz*vz)
    if l < 1e-5:
        return 1
    vx /= l
    vy /= l
    vz /= l
    
    px = (pwm[0], pwm[1], pwm[2])
    py = (pwm[4], pwm[5], pwm[6])
    pz = (pwm[8], pwm[9], pwm[10])
    
    dx = vx*px[0] + vy*px[1] + vz*px[2]
    dy = vx*py[0] + vy*py[1] + vz*py[2]
    dz = vx*pz[0] + vy*pz[1] + vz*pz[2]
    
    axes = [(0, dx), (1, dy), (2, dz), (3, -dx), (4, -dy), (5, -dz)]
    axes.sort(key=lambda x: x[1], reverse=True)
    return axes[0][0]


def createSkirtBellCollider(
    prefix="",
    leftHipObj=None,
    leftKneeObj=None,
    leftHeelObj=None,
    rightHipObj=None,
    rightKneeObj=None,
    rightHeelObj=None,
    waistObj=None,
    skirtType=1,
):
    """Create a SkirtBellCollider. Bell locator is placed at the midpoint
    between the two hip objects (world space)."""

    lp = cmds.xform(leftHipObj, query=True, worldSpace=True, translation=True)
    rp = cmds.xform(rightHipObj, query=True, worldSpace=True, translation=True)
    mid = [(lp[i] + rp[i]) * 0.5 for i in range(3)]

    dist = math.sqrt((lp[0]-rp[0])**2 + (lp[1]-rp[1])**2 + (lp[2]-rp[2])**2)
    if dist < 1e-4:
        dist = 1.0

    bellLoc = cmds.spaceLocator(name=f"{prefix}skirt_locator")[0]
    cmds.xform(bellLoc, worldSpace=True, translation=mid)
    cmds.xform(bellLoc, worldSpace=True, rotation=(180, 0, 0))
    
    if waistObj:
        cmds.parentConstraint(waistObj, bellLoc, maintainOffset=True)

    node = cmds.createNode(SKIRT_NODE_TYPE, name=f"{prefix}skirtBellCollider")
    nodeXform = cmds.listRelatives(node, parent=True, type="transform")[0]
    cmds.rename(nodeXform, f"{prefix}skirtBellCollider_transform")

    cmds.setAttr(f"{node}.skirtType", skirtType)
    cmds.setAttr(f"{node}.bellScale", dist, 1.0, dist, type="double3")
    cmds.setAttr(f"{node}.ringScale", dist / 2.0, 1.0, dist / 2.0, type="double3")
    
    cmds.setAttr(f"{node}.leftRingAxis", _findAxis(leftHipObj, leftKneeObj))
    cmds.setAttr(f"{node}.rightRingAxis", _findAxis(rightHipObj, rightKneeObj))
    cmds.setAttr(f"{node}.bellAxis", 1)  # We set rotation to (180,0,0), so Y is down

    jointPairs = [
        (bellLoc, "bellMatrix"),
        (leftHipObj, "leftHipMatrix"),
        (leftKneeObj, "leftKneeMatrix"),
        (leftHeelObj, "leftHeelMatrix"),
        (rightHipObj, "rightHipMatrix"),
        (rightKneeObj, "rightKneeMatrix"),
        (rightHeelObj, "rightHeelMatrix"),
    ]
    for obj, attr in jointPairs:
        if obj:
            cmds.connectAttr(f"{obj}.worldMatrix[0]", f"{node}.{attr}", force=True)

    # Lock t/r/s on the skirtBellCollider's parent transform
    nodeXform = cmds.listRelatives(node, parent=True)[0]
    for attr in ("tx", "ty", "tz", "rx", "ry", "rz", "sx", "sy", "sz"):
        cmds.setAttr(f"{nodeXform}.{attr}", lock=True)

    # Rebuild the surface with history
    surfXform = cmds.createNode("transform", name=f"{prefix}skirt_surface")
    surfShape = cmds.createNode("nurbsSurface", name=f"{prefix}skirt_surfaceShape", parent=surfXform)
    
    rebuildNode = cmds.createNode("rebuildSurface")
    cmds.setAttr(f"{rebuildNode}.direction", 1) # 0=U, 1=V, 2=U and V
    cmds.setAttr(f"{rebuildNode}.spansU", 1)
    cmds.setAttr(f"{rebuildNode}.spansV", 4)
    
    cmds.connectAttr(f"{node}.outputSurface", f"{rebuildNode}.inputSurface", force=True)
    cmds.connectAttr(f"{rebuildNode}.outputSurface", f"{surfShape}.create", force=True)
    cmds.sets(surfShape, edit=True, forceElement="initialShadingGroup")

    cmds.select(node)
    return node, bellLoc, surfXform


def attachJointsToSurface(surface, u_num, v_num, prefix=""):
    """Creates a grid of joints attached to the given NURBS surface using pointOnSurfaceInfo and fourByFourMatrix."""
    shapes = cmds.listRelatives(surface, shapes=True, path=True)
    if not shapes:
        cmds.warning(f"No shape found for surface {surface}")
        return []
    surfShape = shapes[0]
    
    joints = []
    
    # Group all attachments under a single transform for neatness
    grp = cmds.createNode('transform', name=f"{prefix}attachedJoints_group")
    
    min_u = cmds.getAttr(f"{surfShape}.minValueU")
    max_u = cmds.getAttr(f"{surfShape}.maxValueU")
    min_v = cmds.getAttr(f"{surfShape}.minValueV")
    max_v = cmds.getAttr(f"{surfShape}.maxValueV")
    
    for i in range(u_num):
        # Periodic in U, so [0, 1) to avoid overlap at seam
        u_ratio = i / float(u_num) if u_num > 0 else 0.5
        u_val = min_u + u_ratio * (max_u - min_u)
        for j in range(v_num):
            # Open in V, so [0, 1] exactly
            v_ratio = j / float(max(1, v_num - 1)) if v_num > 1 else 0.5
            v_val = min_v + v_ratio * (max_v - min_v)
            
            # PointOnSurfaceInfo
            posi = cmds.createNode('pointOnSurfaceInfo', name=f"{prefix}skirt_{i}_{j}_posi")
            cmds.connectAttr(f"{surfShape}.worldSpace[0]", f"{posi}.inputSurface", force=True)
            cmds.setAttr(f"{posi}.parameterU", u_val)
            cmds.setAttr(f"{posi}.parameterV", v_val)
            
            # Nodes for matrix assembly
            fourByFour = cmds.createNode("fourByFourMatrix", name=f"{prefix}skirt_{i}_{j}_fourByFourMatrix")
            decompose = cmds.createNode("decomposeMatrix", name=f"{prefix}skirt_{i}_{j}_decomposeMatrix")
            inv = cmds.createNode("multMatrix", name=f"{prefix}skirt_{i}_{j}_multMatrix")
            
            # Transform to hold the joint
            jntGrp = cmds.createNode('transform', name=f"{prefix}skirt_{i}_{j}_transform", parent=grp)
            
            # Feed posi vectors into fourByFourMatrix
            cmds.connectAttr(f"{posi}.positionX", f"{fourByFour}.in30", force=True)
            cmds.connectAttr(f"{posi}.positionY", f"{fourByFour}.in31", force=True)
            cmds.connectAttr(f"{posi}.positionZ", f"{fourByFour}.in32", force=True)
            
            cmds.connectAttr(f"{posi}.normalX", f"{fourByFour}.in00", force=True)
            cmds.connectAttr(f"{posi}.normalY", f"{fourByFour}.in01", force=True)
            cmds.connectAttr(f"{posi}.normalZ", f"{fourByFour}.in02", force=True)
            
            cmds.connectAttr(f"{posi}.tangentUx", f"{fourByFour}.in10", force=True)
            cmds.connectAttr(f"{posi}.tangentUy", f"{fourByFour}.in11", force=True)
            cmds.connectAttr(f"{posi}.tangentUz", f"{fourByFour}.in12", force=True)
            
            cmds.connectAttr(f"{posi}.tangentVx", f"{fourByFour}.in20", force=True)
            cmds.connectAttr(f"{posi}.tangentVy", f"{fourByFour}.in21", force=True)
            cmds.connectAttr(f"{posi}.tangentVz", f"{fourByFour}.in22", force=True)
            
            # Multiply world matrix by parentInverseMatrix to get local
            cmds.connectAttr(f"{fourByFour}.output", f"{inv}.matrixIn[0]", force=True)
            cmds.connectAttr(f"{jntGrp}.parentInverseMatrix[0]", f"{inv}.matrixIn[1]", force=True)
            cmds.connectAttr(f"{inv}.matrixSum", f"{decompose}.inputMatrix", force=True)
            
            # Drive the transform
            cmds.connectAttr(f"{decompose}.outputTranslate", f"{jntGrp}.translate", force=True)
            cmds.connectAttr(f"{decompose}.outputRotate", f"{jntGrp}.rotate", force=True)
            
            # Actual joint
            cmds.select(clear=True)
            jnt = cmds.joint(name=f"{prefix}skirt_{i}_{j}_joint")
            cmds.parent(jnt, jntGrp, relative=True)
            cmds.setAttr(f"{jnt}.radius", 0.5)
            joints.append(jnt)
            
    return joints


# ---------------------------------------------------------------------------
# Shared UI helpers
# ---------------------------------------------------------------------------




def _pickIntoField(tf):
    sel = cmds.ls(selection=True, shortNames=True)
    if sel:
        cmds.textField(tf, edit=True, text=sel[0])
    else:
        cmds.warning("Nothing selected.")


def _objRow(label):
    """Compact label + text field + Pick button row. Returns textField name."""
    cmds.rowLayout(numberOfColumns=3, columnWidth3=(60, 130, 30),
                   adjustableColumn=2,
                   columnAlign3=("right", "left", "left"),
                   columnAttach3=("right", "both", "left"))
    cmds.text(label=label, align="right")
    tf = cmds.textField(placeholderText="<none>", height=20)
    cmds.button(label="<", height=20, width=25,
                command=lambda *_, f=tf: _pickIntoField(f))
    cmds.setParent("..")
    return tf


def _prefixRow():
    """Prefix label + text field. Returns textField name."""
    cmds.rowLayout(numberOfColumns=2, columnWidth2=(60, 160),
                   columnAlign2=("right", "left"))
    cmds.text(label="Prefix", align="right")
    tf = cmds.textField(placeholderText="e.g. hero", text="test_", height=20)
    cmds.setParent("..")
    return tf

# ---------------------------------------------------------------------------
# Tab 1 — Bell Collider
# ---------------------------------------------------------------------------

class BellColliderTab:
    def __init__(self, parent):
        self._build(parent)

    def _build(self, parent):
        cmds.setParent(parent)
        cmds.columnLayout(adjustableColumn=True, rowSpacing=2, columnOffset=("both", 6))
        cmds.separator(height=4, style="none")

        self._prefix = _prefixRow()

        cmds.rowLayout(numberOfColumns=2, columnWidth2=(80, 60), columnAlign2=("right", "left"))
        cmds.text(label="Rings", align="right")
        self._ringCount = cmds.intField(value=1, minValue=1, height=20)
        cmds.setParent("..")

        cmds.separator(height=6, style="none")
        cmds.button(label="Create Bell Collider", height=26,
                    backgroundColor=(0.24, 0.37, 0.60),
                    command=lambda *_: self._create())
        cmds.separator(height=4, style="none")

    def _create(self):
        prefix = cmds.textField(self._prefix, query=True, text=True).strip()
        n = cmds.intField(self._ringCount, query=True, value=True)
        node, bellLoc, ringLocs, curveXform = createBellCollider(numRings=n, prefix=prefix)
        cmds.inViewMessage(
            amg=(f'<hl>BellCollider</hl> "<b>{node}</b>" created '
                 f'with {n} ring{"s" if n != 1 else ""}.'),
            pos="topCenter", fade=True
        )


# ---------------------------------------------------------------------------
# Tab 2 — Skirt Bell Collider
# ---------------------------------------------------------------------------

class SkirtBellColliderTab:
    def __init__(self, parent):
        self._f = {}
        self._build(parent)

    def _build(self, parent):
        cmds.setParent(parent)
        cmds.columnLayout(adjustableColumn=True, rowSpacing=2, columnOffset=("both", 6))
        cmds.separator(height=4, style="none")

        self._f["prefix"] = _prefixRow()

        cmds.separator(height=6, style="in")
        cmds.text(label=" Joint Connections", align="left", font="boldLabelFont", height=16)
        self._f["waist"] = _objRow("Waist")
        
        cmds.rowColumnLayout(numberOfColumns=2, columnWidth=[(1, 220), (2, 220)], columnSpacing=[(2, 10)])
        for rKey, rLabel, lKey, lLabel in [
            ("rightHip", "R Hip", "leftHip", "L Hip"),
            ("rightKnee", "R Knee", "leftKnee", "L Knee"),
            ("rightHeel", "R Heel", "leftHeel", "L Heel"),
        ]:
            self._f[rKey] = _objRow(rLabel)
            self._f[lKey] = _objRow(lLabel)
        cmds.setParent("..")

        cmds.separator(height=6, style="none")
        cmds.button(label="Create Skirt Bell Collider", height=26,
                    backgroundColor=(0.24, 0.50, 0.38),
                    command=lambda *_: self._create())
        cmds.separator(height=4, style="none")

        # --- Attach Joints Section ---
        cmds.separator(height=6, style="in")
        cmds.text(label=" Attach Joints to Surface", align="left", font="boldLabelFont", height=16)

        self._f["attachSurf"] = _objRow("Surface")

        cmds.rowLayout(numberOfColumns=4, columnWidth4=(60, 50, 50, 50), columnAlign4=("right", "left", "right", "left"))
        cmds.text(label="U Joints", align="right")
        self._f["uJoints"] = cmds.intField(value=7, minValue=1, height=20)
        cmds.text(label="V Joints", align="right")
        self._f["vJoints"] = cmds.intField(value=5, minValue=1, height=20)
        cmds.setParent("..")

        cmds.separator(height=6, style="none")
        cmds.button(label="Attach Joints", height=26,
                    backgroundColor=(0.38, 0.24, 0.50),
                    command=lambda *_: self._attachJoints())
        cmds.separator(height=4, style="none")

    def _create(self):
        f = self._f
        prefix = cmds.textField(f["prefix"], query=True, text=True).strip()

        required = ["leftHip", "leftKnee", "rightHip", "rightKnee"]
        missing = [k for k in required if not (cmds.textField(f[k], query=True, text=True).strip() or None)]
        if missing:
            cmds.warning("Skirt Bell Collider: missing — " + ", ".join(missing))
            return

        leftHeel = cmds.textField(f["leftHeel"], query=True, text=True).strip() or None
        rightHeel = cmds.textField(f["rightHeel"], query=True, text=True).strip() or None
        heelsPresent = leftHeel and rightHeel
        skirtType = 1 if heelsPresent else 0  # Long if heels provided, Short otherwise

        node, bellLoc, surfXform = createSkirtBellCollider(
            prefix=prefix,
            waistObj=cmds.textField(f["waist"], query=True, text=True).strip() or None,
            leftHipObj=cmds.textField(f["leftHip"], query=True, text=True).strip() or None,
            leftKneeObj=cmds.textField(f["leftKnee"], query=True, text=True).strip() or None,
            leftHeelObj=leftHeel,
            rightHipObj=cmds.textField(f["rightHip"], query=True, text=True).strip() or None,
            rightKneeObj=cmds.textField(f["rightKnee"], query=True, text=True).strip() or None,
            rightHeelObj=rightHeel,
            skirtType=skirtType,
        )
        cmds.inViewMessage(
            amg=f'<hl>SkirtBellCollider</hl> "<b>{node}</b>" created '
                f'({"Long" if heelsPresent else "Short"}).',
            pos="topCenter", fade=True
        )

    def _attachJoints(self):
        f = self._f
        prefix = cmds.textField(f["prefix"], query=True, text=True).strip()
        surface = cmds.textField(f["attachSurf"], query=True, text=True).strip() or None
        
        if not surface:
            cmds.warning("Please specify a target surface.")
            return
            
        u_num = cmds.intField(f["uJoints"], query=True, value=True)
        v_num = cmds.intField(f["vJoints"], query=True, value=True)
        
        joints = attachJointsToSurface(surface, u_num, v_num, prefix=prefix)
        if joints:
            cmds.inViewMessage(
                amg=f'<hl>Attach Joints</hl> Created {len(joints)} joints on "<b>{surface}</b>".',
                pos="topCenter", fade=True
            )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def show():
    """Open (or re-open) the Bell Collider Creator window."""
    if cmds.window(WINDOW_ID, exists=True):
        cmds.deleteUI(WINDOW_ID, window=True)

    win = cmds.window(WINDOW_ID, title=WINDOW_TITLE,
                      widthHeight=(480, 300), sizeable=True, resizeToFitChildren=True)

    cmds.scrollLayout(childResizable=True)
    tabs = cmds.tabLayout(innerMarginWidth=2, innerMarginHeight=2)

    tab1 = cmds.columnLayout(adjustableColumn=True)
    BellColliderTab(tab1)
    cmds.setParent(tabs)

    tab2 = cmds.columnLayout(adjustableColumn=True)
    SkirtBellColliderTab(tab2)
    cmds.setParent(tabs)

    cmds.tabLayout(tabs, edit=True,
                   tabLabel=[(tab1, "Bell Collider"),
                             (tab2, "Skirt Bell Collider")])

    cmds.showWindow(win)


if __name__ == "__main__":
    show()
