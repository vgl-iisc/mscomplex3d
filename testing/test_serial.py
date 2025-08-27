from multiprocessing import Pool, cpu_count
from math import ceil

import pickle
import pytest
import pyms3d_core as pyms3d

def msc_worker(msc, start, end):
    num = 0

    for i in range(start, end):
        num += msc.cp_func(i)
        
    return num

def test_pickling_unpickling():
    pyms3d.select_device()
    msc = pyms3d.MsComplex()

    msc.compute_bin("testing/Hydrogen_128x128x128.raw", (128, 128, 128))
    
    dump = pickle.dumps(msc)    
    msc_new = pickle.loads(dump)

    assert msc.num_cps() == msc_new.num_cps()

    errors = []
    for i in range(msc.num_cps()):
        if not msc.cp_func(i) == msc_new.cp_func(i):
            errors.append(f"cp function value doesn't match at cp {i}")

    if len(errors) > 0:
        pytest.fail("\n".join(errors))

def test_multiprocessing_cp_sum():
    pyms3d.select_device()
    msc = pyms3d.MsComplex()

    msc.compute_bin("testing/Hydrogen_128x128x128.raw", (128, 128, 128))

    n_processes = cpu_count()
    N = msc.num_cps()
    task_size = ceil(N / n_processes)

    ranges = [[start, start + task_size] for start in range(0, N, task_size)]
    ranges[-1][-1] = N
    ranges = list(map(tuple, ranges))

    tasks = Pool(n_processes)

    res = tasks.starmap(msc_worker, [(msc,) + ranges[i] for i in range(n_processes)])

    total = sum(res)

    all_funcs = msc.cps_func()

    assert abs(total - all_funcs.sum()) <= 1e-3, "incorrect sum"