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
