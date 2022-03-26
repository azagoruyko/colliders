# Colliders

My implementation of Bell Collider. It supports multiple rings and viewport 2.0 debug drawing. It uses mesh intersection approach.

![bell](https://user-images.githubusercontent.com/9614751/159114116-ec6ac52e-00fd-41e1-8475-8e068ed21467.PNG)

Youtube: https://www.youtube.com/watch?v=89Dbd-n8EzY

## How to run
#### Compile C++ plugin for Maya version you want.
You need Visual Studio and CMake to do this. You can find precompiled plugins for Maya in *plugins* folder.

#### Python script to create Bell Collider node with some bell rings.
```python
import pymel.core as pm

name = "test"
numRings = 2
numOutJoints = 6

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
grp | curve
bellCollider.outputCurve >> curve.create

if numOutJoints > 0:
    bellCollider.positionCount.set(numOutJoints)
    for i in range(numOutJoints):
        tr = pm.createNode("transform", n=name+"_"+str(i+1)+"_transform", p=grp)
        bellCollider.outputPositions[i] >> tr.t
        bellCollider.outputRotations[i] >> tr.r
```        

## Thanks and gracias

If you want to support me you're welcome!

197VXFacHv8VsYR1fQcs65sugBBQm1FSa5 BTC<br>
0x27c561D37a9163ffeAFB9FC9999f088aE4544F3B ETH<br>
0x27c561D37a9163ffeAFB9FC9999f088aE4544F3B BSC BEP20<br>
bnb1jfguvl4dwcg3ldx3vnv9qm23xve237sj5wyz26  BTCB BEP2<br>
TWBKt2EYHS7uArz6Gbsz3bnjg293Jzdh6m  TRX USDT TRC20<br>
2N2o6EviF3yhJZ5tSreHshJ77FTkBNz1X5EZMkXN18Z9 SOLANA<br>
