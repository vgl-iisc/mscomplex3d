# ugly hack to make package structure consistent across operating systems

try:
    from .pyms3d_core.pyms3d_core import *
except ImportError:
    from .pyms3d_core import *  # this loads all symbols from your compiled module