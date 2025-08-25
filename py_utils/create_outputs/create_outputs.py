import pyms3d_core as pyms3d
import pickle
import os

indir = "."
outdir = "out"

# simplification levels
levels = sorted([0.05])

os.makedirs(outdir, exist_ok=True)

pyms3d.get_hw_info()

for file in os.listdir(indir):
    if not file.endswith(".raw"):
        continue

    name, _ = os.path.splitext(file)
    dims = tuple(map(int, name.split("_")[-1].split("x")))
    
    msc = pyms3d.MsComplexPyms3D()
    msc.compute_bin(file, dims)

    base_state = msc.dbg_serialize()
    
    with open(f"stored_outputs/{name}.pkl", "wb") as f:
        pickle.dump(base_state, f)

    for lvl in levels:
        msc.simplify_pers(lvl)
        simp_state = msc.dbg_serialize()

        with open(f"stored_outputs/{name}_{lvl}.pkl", "wb") as f:
            pickle.dump(simp_state, f)