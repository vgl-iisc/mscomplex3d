from multiprocessing import Pool, cpu_count
from math import ceil
import os
import sys

sys.path.append("build/pyms3d/Debug")
import pyms3d_core as pyms3d

def msc_worker(msc, start, end):
    num = 0
    
    # total = end - start

    for i in range(start, end):
        num += msc.cp_func(i)
        
        # if (i - start) % (total // 10) == 0:
        #     print(f"{start}: {i - start}/{total}")

    return num

def test_multiprocessing_cp_sum():
    pyms3d.get_hw_info()
    msc = pyms3d.MsComplexPyms3D()

    msc.compute_bin("Hydrogen_128x128x128.raw", (128, 128, 128))

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

    assert abs(total - all_funcs.sum()) <= 1.0