import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns
import matplotlib as mpl

def load_fort46(data_path):
    """Load result from fort.46 file with EIRENE output on EIRENE grid.

    Returns
    -------
    out : dict
        Dictionary for each loaded file containing a subdictionary with keys for each loaded field from each file.
    """
    out = {}

    out = {}
    # load each of these files into dictionary structures
    with open(data_path, "r") as f:
        contents = f.readlines()

    out["NTRII"] = contents[0].split()[0]

    NATM = out["NATM"] = int(contents[1].split()[0][0])
    NMOL = out["NMOL"] = int(contents[1].split()[1][0])
    NION = out["NION"] = int(contents[1].split()[2][0])

    ii = 2
    out["species"] = []
    while not contents[ii].startswith("*eirene"):
        if not contents[ii].split()[0].isnumeric():
            out["species"].append(contents[ii].split()[0])
        ii += 1

    while ii < len(contents):
        if contents[ii].startswith("*eirene"):
            key = "".join(contents[ii].split()[3:-3])
            dim = int(contents[ii].split()[-1])
            out[key] = []
        elif not contents[ii].split()[0][0].isalpha():
            [out[key].append(float(val))
             for val in contents[ii].strip().split()]
        ii += 1

    for key in out:

        # set to the right shape, depending on whether quantity is atomic, molecular or ionic
        if key.endswith("a"):
            num = NATM
        elif key.endswith("m"):
            num = NMOL
        elif key.endswith("i"):
            num = NION
        else:
            num = 1
        if key in ["species", "NTRII"]:
            continue
        out[key] = np.array(out[key]).reshape(-1, num, order="F")

    return out

baserun_path = "/home/bianyu/nt/triangle"
i = load_fort46(baserun_path+r'/fort.46')