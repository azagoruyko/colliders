# Colliders

My implementation of Bell Collider. It supports multiple rings and viewport 2.0 debug drawing. 

![bell](https://user-images.githubusercontent.com/9614751/170646159-575afcd8-76b3-47cb-ad13-65deb6a56942.PNG)

Youtube: https://www.youtube.com/watch?v=89Dbd-n8EzY

## How to run
#### Compile C++ plugin for Maya version you want.
You need Visual Studio and CMake to do this. 

#### Python script to create Bell Collider node with some bell rings.
```python
import pymel.core as pm

name = "test"
numRings = 2
numOutJoints = 0

grp = pm.createNode("transform", n=name+"_bellCollider_group")
grp.addAttr("scaleFactor", min=0.001, dv=1, k=True)
grp.scaleFactor >> grp.sx
grp.scaleFactor >> grp.sy
grp.scaleFactor >> grp.sz
grp.sx.set(l=True, k=False)
grp.sy.set(l=True, k=False)
grp.sz.set(l=True, k=False)

bellCollider = pm.createNode("bellCollider", n=name+"_bellCollider")
bellColliderParent = bellCollider.getParent()
bellColliderParent.rename(name+"_bellCollider_transform")
bellColliderParent.t.lock()
bellColliderParent.r.lock()
bellColliderParent.s.lock()
grp | bellColliderParent

bell = pm.createNode("transform", n=name+"_bell_transform", p=grp)
bell.m >> bellCollider.bellMatrix

for i in range(numRings):
    ring = pm.createNode("transform", n=name+"_ring_"+str(i+1)+"_transform", p=grp)
    ring.s.set([0.2, 2, 0.2])
    ring.m >> bellCollider.ringMatrix[i]

curve = pm.createNode("nurbsCurve", n=name+"_bellCollider_curveShape")
curve.getParent().rename(name+"_bellCollider_curve")
grp | curve.getParent()
bellCollider.outputCurve >> curve.create

bellMesh = pm.createNode("mesh", n=name+"_bellCollider_mesh")
bellMesh.getParent().rename(name+"_bellCollider_mesh")
grp | bellMesh.getParent()
bellCollider.outputBellMesh >> bellMesh.inMesh

bellCollider.positionCount.set(numOutJoints)
for i in range(numOutJoints):
    tr = pm.createNode("transform", n=name+"_"+str(i+1)+"_transform", p=grp)
    bellCollider.outputPositions[i] >> tr.t
    bellCollider.outputRotations[i] >> tr.r
```
