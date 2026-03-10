import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns
import matplotlib as mpl
import load_eirene_mesh
import load_fort46

baserun_path = "/home/bianyu/nt/triangle"
nodes = np.fromfile(
    baserun_path+r'/fort.33', sep=" "
)
n = int(nodes[0])
xnodes = nodes[1: n + 1] / 100  # cm -->m
ynodes = nodes[n + 1:] / 100

# EIRENE triangulation
triangles = np.loadtxt(
    baserun_path+r'/fort.34',
    skiprows=1,
    usecols=(1, 2, 3),
)


'''
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))
for i in range(n):
    ax.plot([xnodes[triangles[i, 1]],ynodes[triangles[i, 1]]], [xnodes[triangles[i, 2]],ynodes[triangles[i, 2]]])
ax.show()
'''