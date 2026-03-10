#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@File    :   plot_mesh_test.py
@Time    :   2023/11/27 10:20:22
@Author  :   Xuele Zhao
@Contact :   zhaoxuele@mail.dlut.edu.cn
@Site    :   http://plasmaphys.dlut.edu.cn/
'''

import re
import numpy as np
import matplotlib.pyplot as plt

def load_eirene_mesh(baserun_path):
        """Load EIRENE nodes from the fort.33 file and triangulation from the fort.34 file"""
        nodes = np.fromfile(
            baserun_path+r'/fort.33', sep=" "
        )
        n = int(nodes[0])
        xnodes = nodes[1 : n + 1] / 100  # cm -->m
        ynodes = nodes[n + 1 :] / 100

        # EIRENE triangulation
        triangles = np.loadtxt(
            baserun_path+r'/fort.34',
            skiprows=1,
            usecols=(1, 2, 3),
        )
        tag = np.loadtxt(
            baserun_path+r'/fort.34',
            skiprows=1,
            usecols=(0),
        )
        return xnodes, ynodes, triangles - 1,tag-1

if __name__ == "__main__":
    fig=plt.figure(1)
    ax=fig.gca()
    case_path=r'/home/bianyu/nt/triangle'
    
    xnodes,ynodes,triangles,tags=load_eirene_mesh(case_path)
    ax.triplot(xnodes,ynodes,triangles,color='g',lw=0.2)
    center_x=[]
    center_y=[]
    for i,test in enumerate(triangles):
        x_total=0
        y_total=0
        for j in test:
            x_total=x_total+xnodes[int(j)]
            y_total=y_total+ynodes[int(j)]
            ax.annotate(int(j), (xnodes[int(j)], ynodes[int(j)]), color='black', weight='bold', ha='center', va='center')
        center_x=x_total/3
        center_y=y_total/3
        #ax.annotate(int(tags[i]+1), (center_x, center_y), color='k', weight='bold',  ha='center', va='center')
    ax.axis('equal')
    plt.show()


