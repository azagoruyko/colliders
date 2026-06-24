# Colliders

My implementation of Bell Collider and Skirt Bell Collider for long/short skirts. 

![bell](https://user-images.githubusercontent.com/9614751/170646159-575afcd8-76b3-47cb-ad13-65deb6a56942.PNG)

Youtube: https://www.youtube.com/watch?v=89Dbd-n8EzY

## Features
- **Bell Collider**: Multi-ring bell collider.
- **Skirt Bell Collider**: Mixture of Bell Colliders specifically for long/shirt skirts in a single node.
- **Simple UI**

## Installation
#### Compile C++ plugin for Maya version you want.
You need Visual Studio and CMake to do this. 

#### Using Python
Use the script `python/colliders.py` to run the UI.

```python
import sys
sys.path.insert(0, r"<path_to_colliders>/python")
import colliders
colliders.show()
```
