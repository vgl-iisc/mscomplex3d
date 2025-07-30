from timeit import default_timer
import MSDebugFileCreator
import sys

sys.path.append("../../build/pyms3d/Debug")

# Now you can import your custom module
import pyms3d_core as pyms3d

DataFile1 = "../../testing/Hydrogen_128x128x128.raw"
Dim1      = (128,128,128)

def print_num_cps(msc):
    print(", ".join(list(map(lambda arr: str(len(arr)), [msc.cps(i) for i in range(4)]))))


pyms3d.get_hw_info(0)

msc_gpu = pyms3d.MsComplexPyms3D()

start = default_timer()
msc_gpu.compute_bin(DataFile1,Dim1)
end = default_timer()

print(f"gpu computation time: {end - start}s")
print("gpu pre-simplification")
print_num_cps(msc_gpu)
state_gpu = msc_gpu.dbg_serialize()
MSDebugFileCreator.createDebugFile(msc_gpu, "gpu_raw")

start = default_timer()
msc_gpu.simplify_pers(thresh=0.05)
end = default_timer()

print(f"gpu simplification time: {end - start}s")
print("gpu post-simplification")
print_num_cps(msc_gpu)
MSDebugFileCreator.createDebugFile(msc_gpu, "gpu")

print("done with gpu, cpu time")
sys.stdout.flush()

pyms3d.get_hw_info(1)
msc_cpu = pyms3d.MsComplexPyms3D()

start = default_timer()
msc_cpu.compute_bin(DataFile1,Dim1)
end = default_timer()

print(f"cpu computation time: {end - start}s")
print("cpu pre-simplification")
print_num_cps(msc_cpu)
state_cpu = msc_cpu.dbg_serialize()
MSDebugFileCreator.createDebugFile(msc_cpu, "cpu_raw")

start = default_timer()
msc_cpu.simplify_pers(thresh=0.05)
end = default_timer()

print(f"cpu simplification time: {end - start}s")
print("cpu post-simplification")
print_num_cps(msc_cpu)

MSDebugFileCreator.createDebugFile(msc_cpu, "cpu")

print("done")
sys.stdout.flush()

print("serialised string test")

def diff_ns_dict(dictA, dictB, prefix=""):
    children = []
    
    dev = []

    assert dictA.keys() == dictB.keys()
    
    for k in dictA.keys():
        valA = dictA[k]
        valB = dictB[k]

        if isinstance(valA, dict):
            assert isinstance(valB, dict)
            children.append(k)
            continue

        if valA != valB:
            dev.append(f"deviation in {prefix}{k}")
        
    for k in children:
        dev.extend(diff_ns_dict(dictA[k], dictB[k], k + "/"))

    return dev

deviation = diff_ns_dict(state_gpu, state_cpu)
print("\n".join(deviation))
assert len(deviation) == 0