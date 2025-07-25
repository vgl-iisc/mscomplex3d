import pickle
import pytest
import os
import sys

sys.path.append("./build/pyms3d/Debug")
import pyms3d_core as pyms3d

def test_pickling_unpickling():
    pyms3d.get_hw_info()
    msc = pyms3d.MsComplexPyms3D()

    msc.compute_bin("Hydrogen_128x128x128.raw", (128, 128, 128))
    
    dump = pickle.dumps(msc)    
    msc_new = pickle.loads(dump)

    assert msc.num_cps() == msc_new.num_cps()

    for i in range(msc.num_cps()):
        if not msc.cp_func(i) == msc_new.cp_func(i):
            pytest.fail(f"cp function value doesn't match at cp {i}")

test_pickling_unpickling()