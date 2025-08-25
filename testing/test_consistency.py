import pyms3d_core as pyms3d
from diff_state_dict import diff_ns_dict

import pickle
import os

def test_gpu_cpu_agreement():
    df = "testing/Hydrogen_128x128x128.raw"
    dim      = (128,128,128)

    pyms3d.get_hw_info(0)
    msc_gpu = pyms3d.MsComplexPyms3D()
    msc_gpu.compute_bin(df, dim)
    state_gpu = msc_gpu.dbg_serialize()

    pyms3d.get_hw_info(1)
    msc_cpu = pyms3d.MsComplexPyms3D()
    msc_cpu.compute_bin(df, dim)
    state_cpu = msc_cpu.dbg_serialize()

    deviation = diff_ns_dict(state_gpu, state_cpu)
    assert len(deviation) == 0, "deviations b/w cpu & gpu:\n" + "\n".join(deviation)

def test_output_match():
    pyms3d.get_hw_info()

    outputs = os.listdir("testing/stored_outputs")

    name_dim_levs = {}
    for out in outputs:
        if not out.endswith(".pkl"):
            continue
        
        name, _ = os.path.splitext(out)
        parts = name.split("_")

        simp = None if len(parts) < 3 else float(parts[-1])
        dims = tuple(map(int, parts[1].split("x")))

        basename = "_".join(parts) if simp is None else "_".join(parts[:-1])

        if basename not in name_dim_levs:
            name_dim_levs[basename] = (dims, [] if simp is None else [simp])
        elif simp is not None:
            name_dim_levs[basename][1].append(simp)

    for name, (dim, levs) in name_dim_levs.items():
        msc = pyms3d.MsComplexPyms3D()
        msc.compute_bin(f"testing/{name}.raw", dim)

        state_dict = msc.dbg_serialize()

        with open(f"testing/stored_outputs/{name}.pkl", "rb") as f:
            target_dict = pickle.load(f)

        dev = diff_ns_dict(state_dict, target_dict)

        assert len(dev) == 0, f"deviations in unsimplified msc for {name}:\n" + "\n".join(dev)

        levs_s = sorted(levs)

        for lvl in levs_s:
            msc.simplify_pers(lvl)

            state_dict = msc.dbg_serialize()

            with open(f"testing/stored_outputs/{name}_{lvl}.pkl", "rb") as f:
                target_dict = pickle.load(f)

            dev = diff_ns_dict(state_dict, target_dict)

            assert len(dev) == 0, f"deviations in {lvl}-simplified msc for {name}:\n" + "\n".join(dev)