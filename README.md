# Colliders

My implementation of `Bell Collider` as well as custom `Skirt Bell Collider` for skirts.

<img width="977" height="550" alt="image" src="https://github.com/user-attachments/assets/922bf5db-6b6c-40f4-8ed8-0e4748e0cf30" />

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
