import MSOutput
import os

def createDebugFile(msc, suffix=""):
    
    cellids = msc.cps_cellid()

    # Collect all the data into a list of tuples
    data = []
    for cpInd in [0, 1, 2, 3]:
        for cp in range(len(msc.cps(cpInd))):
            cp_value = msc.cps(cpInd)[cp]
            cell_id = cellids[cp_value]
            func_value = msc.cp_func(cp_value)
            data.append((cell_id, cpInd, func_value))

    # Sort the data based on the third element (msc.cp_func values)
    data.sort(key=lambda x: x[2])

    filetext = ""
    for cell_id, cpInd, func_value in data:
        filetext += f"{cell_id}\t{cpInd}\t{func_value}\n"

    data = []
    for cpInd in [0, 1, 2, 3]:
        for cp in range(len(msc.cps(cpInd))):
            cp_value1 = msc.cps(cpInd)[cp]
            cell_id = cellids[cp_value1]  # 1x3 vector
            for c in msc.asc(cp_value1):
                asc_cell_id = cellids[c[0]]  # 1x3 vector
                # Convert the 1x3 vectors to string for easier sorting and writing
                data.append((tuple(cell_id), tuple(asc_cell_id)))

    # Sort the data based on the first vector (cell_id)
    data.sort(key=lambda x: (x[0], x[1]))


    filetext+="\n"
    for cell_id, asc_cell_id in data:
        # Format each line with proper vector representation
        cell_id_str = f"[{cell_id[0]}, {cell_id[1]}, {cell_id[2]}]"
        asc_cell_id_str = f"[{asc_cell_id[0]}, {asc_cell_id[1]}, {asc_cell_id[2]}]"
        filetext += f"asc {cell_id_str}\t-->\t{asc_cell_id_str}\n"

    ###############################

    data = []
    for cpInd in [0, 1, 2, 3]:
        for cp in range(len(msc.cps(cpInd))):
            cp_value1 = msc.cps(cpInd)[cp]
            cell_id = cellids[cp_value1]  # 1x3 vector
            for c in msc.des(cp_value1):
                asc_cell_id = cellids[c[0]]  # 1x3 vector
                # Convert the 1x3 vectors to string for easier sorting and writing
                data.append((tuple(cell_id), tuple(asc_cell_id)))

    # Sort the data based on the first vector (cell_id)
    data.sort(key=lambda x: (x[0], x[1]))

    filetext+="\n"
    for cell_id, asc_cell_id in data:
        # Format each line with proper vector representation
        cell_id_str = f"[{cell_id[0]}, {cell_id[1]}, {cell_id[2]}]"
        asc_cell_id_str = f"[{asc_cell_id[0]}, {asc_cell_id[1]}, {asc_cell_id[2]}]"
        filetext += f"des {cell_id_str}\t-->\t{asc_cell_id_str}\n"


    suffix = f"_{suffix}" if suffix is not None else ""
    pathname = lambda i: f"scalar{suffix}_{i}.txt" if i is not None else f"scalar{suffix}.txt"
    
    i = None
    if os.path.exists(pathname(i)):
        i = 0
        while os.path.exists(pathname(i)):
            i += 1        

    MSOutput.WriteToOutputFile(filetext, pathname(i))