import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns
import matplotlib as mpl
from load_fort46 import load_fort46

baserun_path = "/home/bianyu/nt/triangle"
i = load_fort46(baserun_path+r'/fort.46')
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
vals = i['pdena']
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))
norm1 = mpl.colors.LogNorm()
cntr = ax.tripcolor(xnodes,ynodes,triangles-1,cmap=mpl.cm.jet, facecolors=vals.flatten(),norm=norm1)
fig.colorbar(cntr,ax=ax)
plt.axis('equal')
plt.show()
'''
for i in range(n):
    ax.plot([xnodes[triangles[i, 1]],ynodes[triangles[i, 1]]], [xnodes[triangles[i, 2]],ynodes[triangles[i, 2]]])
ax.show()'''