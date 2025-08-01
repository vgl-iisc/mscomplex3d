from multiprocessing import Pool
import numpy as np
import requests
import argparse
import shutil
import json
import os

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("path", help="path to the datasets.json file", default="datasets.json")
    parser.add_argument("-o", "--outdir", help="output directory for datasets", default="datasets", required=True)
    parser.add_argument("--maxsize", help="only download files of this size (MB) or less", type=int, default=500)
    parser.add_argument("--maxtotalsize", help="stop downloading files if the directory grows larger than this size (MB)", type=int, default=1024*10)
    parser.add_argument("--no-float32", help="don't cast raw files to float32", action="store_false", dest="float32")
    parser.add_argument("--no-skip64", help="don't skip files with 64-bit elements", action="store_false", dest="skip64")

    return parser.parse_args()

def get_element_size(metadata):
    dtype = metadata["type"]

    i = len(dtype) - 1
    
    while dtype[i].isdigit():
        i -= 1

    return int(dtype[i+1:])

def download_worker(id, tasks, outdir):
    N = len(tasks)
    
    for i, (name, url) in enumerate(tasks):
        resp = requests.get(url, timeout=600)
        file = resp.content

        with open(os.path.join(outdir, name + ".raw"), "wb") as f:
            f.write(file)

        print(f"worker {id}: {i+1}/{N}")

def download_parallel(download_tasks, outdir):
    workers = 4

    work = [(i, download_tasks[i::workers], outdir) for i in range(workers)]
    
    pool = Pool(workers)
    pool.starmap(download_worker, work)

def transform_f32(download_dir, out_dir, meta_by_name):
    for download in os.listdir(download_dir):
        if not download.endswith(".raw"):
            continue

        name = download.split(".")[0]
        dims = meta_by_name[name]["size"]

        arr = np.fromfile(os.path.join(download_dir, download), dtype=meta_by_name[name]["type"]).reshape(dims)
        arr = arr.astype(np.float32)

        arr.tofile(os.path.join(out_dir, name + "_" + "_".join(map(str, dims)) + "_f32.raw"))
        
def main():
    args = parse_args()

    float32 = args.float32
    skip64 = args.skip64

    max_size = args.maxsize
    max_total_size = args.maxtotalsize
    out_dir = args.outdir

    os.makedirs(out_dir, exist_ok=True)

    path = args.path
    
    with open(path, "r", encoding="utf8") as f:
        datasets = json.load(f)
    
    dir_size = 0
    download_tasks = []
    meta_by_name = {}

    for _, meta in datasets.items():
        # size in bytes
        el_size = get_element_size(meta) // 8
        name = meta["name"]

        if el_size > 4 and skip64:
            print(f"skipping: {name} of type {meta['type']} sized {el_size}B")
            continue

        dims = tuple(meta["size"])

        assert len(dims) == 3

        # file in MB
        file_size = (el_size * dims[0] * dims[1] * dims[2]) / (1024**2)

        if file_size > max_size:
            print(f"skipping: {name} is of size {file_size} > {max_size} MB")
            continue

        if dir_size + file_size > max_total_size:
            print(f"skipping: {name} is of size {file_size}, which would make directory {dir_size + file_size} > {max_total_size} MB")
            continue

        download_tasks.append((name, meta["url"]))

        dir_size += file_size
        print(f"downloading: {name} of size {file_size}, dir_size: {dir_size} / {max_total_size}")
        meta_by_name[name] = meta

    print(f"starting downloads")

    download_dir = out_dir
    if float32:
        download_dir = os.path.join(out_dir, "tmp")
        os.makedirs(download_dir, exist_ok=True)

    download_parallel(download_tasks, download_dir)

    if float32:
        transform_f32(download_dir, out_dir, meta_by_name)
        shutil.rmtree(download_dir)
        

if __name__ == "__main__":
    main()