import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.cm
import copy
import os


config = {'font.family': 'serif', 'mathtext.fontset': 'stix', 'font.size': 20,
          'font.weight':'bold','legend.frameon':False,'figure.subplot.wspace':0,
          'legend.loc':'best','legend.fontsize':20,'figure.subplot.hspace':0,
          'ytick.minor.visible':True,'xtick.minor.visible':True,
          'xtick.minor.top':True,'xtick.minor.bottom':True,
          'ytick.minor.left':True,'ytick.minor.right':True,
          'ytick.major.left':True,'ytick.major.right':True,
          'xtick.major.top':True,'xtick.major.bottom':True,
          'ytick.labelright':False,'ytick.labelleft':True,
          'xtick.labelbottom':True,'xtick.labeltop':False,
          'xtick.bottom':True,'xtick.top':True,
          'ytick.left':True,'ytick.right':True,
          'xtick.labelsize':20,'ytick.labelsize':20,
          'figure.autolayout':True,'xtick.direction':'in','ytick.direction':'in',
          'figure.figsize':(10, 10),'axes.formatter.use_mathtext':True,
          'axes.facecolor':'None','axes.formatter.use_mathtext':True,
          'axes.formatter.limits':(-3,3)

          }
rcParams.update(config)

mass_D = 3.34e-27

class SOLPS:
    def __init__(self, path, IfT):
        self.path = path
        self.IfT = IfT

    def SOLPSread(self):
        self.n_D_0 = np.loadtxt(fname=self.path + "2D_data/nDatom_2D.data", skiprows=1)
        self.n_D_0 = np.transpose(self.n_D_0[:, 3].reshape(38, 98))
        self.T_D_0 = np.loadtxt(fname=self.path + "2D_data/tDatom_2D.data", skiprows=1)
        self.T_D_0 = np.transpose(self.T_D_0[:, 3].reshape(38, 98))
        self.n_D_1 = np.loadtxt(fname=self.path + "2D_data/nDion_2D.data", skiprows=1)
        self.n_D_1 = np.transpose(self.n_D_1[:, 3].reshape(38, 98))
        self.n_e_read = np.loadtxt(fname=self.path + "2D_data/ne_2D.dat", skiprows=1)
        self.n_e = np.transpose(self.n_e_read[:, 3].reshape(38, 98))
        self.T_e = np.loadtxt(fname=self.path + "2D_data/te_2D.dat", skiprows=1)
        self.T_e = np.transpose(self.T_e[:, 3].reshape(38, 98))
        self.T_i = np.loadtxt(fname=self.path + "2D_data/ti_2D.dat", skiprows=1)
        self.T_i = np.transpose(self.T_i[:, 3].reshape(38, 98))
        self.Vol = np.loadtxt(fname=self.path + "2D_data/vol_2D.dat", skiprows=1)
        self.Vol = np.transpose(self.Vol[:, 3].reshape(38, 98))
        self.Flux_i_recyc_l = np.loadtxt(self.path + "profiles_data/recycled_neutral_flux_D_l.data", skiprows=1)
        self.Flux_i_recyc_r = np.loadtxt(self.path + "profiles_data/recycled_neutral_flux_D_r.data", skiprows=1)
        self.Ua_Dion_2D = np.loadtxt(self.path + "2D_data/ua_Dion_2D.dat", skiprows=1)
        self.Ua_Dion_2D = np.transpose(self.Ua_Dion_2D[:, 3].reshape(38, 98))
        self.sign_Ua = np.ones((98, 38))
        if (self.IfT):
            self.n_T_0 = np.loadtxt(fname=self.path + "2D_data/nTatom_2D.data", skiprows=1)
            self.n_T_0 = np.transpose(self.n_T_0[:, 3].reshape(38, 98))
            self.n_T_1 = np.loadtxt(fname=self.path + "2D_data/nTion_2D.data", skiprows=1)
            self.n_T_1 = np.transpose(self.n_T_1[:, 3].reshape(38, 98))
            self.n_e = np.loadtxt(fname=self.path + "2D_data/ne_2D.dat", skiprows=1)
            self.n_e = np.transpose(self.n_e[:, 3].reshape(38, 98))
            self.T_e = np.loadtxt(fname=self.path + "2D_data/te_2D.dat", skiprows=1)
            self.T_e = np.transpose(self.T_e[:, 3].reshape(38, 98))

        
        self.Te_1D_Inner = np.transpose(self.T_e[1, 1:37].reshape(1,36))
        self.Ti_1D_Inner = np.transpose(self.T_i[1, 1:37].reshape(1,36))
        self.Te_1D_Outer = np.transpose(self.T_e[96, 1:37].reshape(1,36))
        self.Ti_1D_Outer = np.transpose(self.T_i[96, 1:37].reshape(1,36))

        self.Te_matrix_1 = np.zeros([98,38])
        self.Te_matrix_2 = np.zeros([98,38])
        self.Te_matrix_3 = np.zeros([98,38])
        self.Te_matrix_4 = np.zeros([98,38])
        self.Te_matrix_5 = np.zeros([98,38])
        self.Te_matrix_6 = np.zeros([98,38])
        self.Te_matrix = [self.Te_matrix_1, self.Te_matrix_2, self.Te_matrix_3, self.Te_matrix_4, self.Te_matrix_5, self.Te_matrix_6]

        self.Mu = self.n_D_1 * np.abs(self.Ua_Dion_2D) * 3.34e-27 * self.Vol
    
    def Te_matrix_cal(self):
        for i in range(98):
            for j in range(38):
                if (self.T_e[i, j] < 0.5): 
                    self.Te_matrix_1[i, j] = 1
                elif (self.T_e[i, j] < 1):
                    self.Te_matrix_2[i, j] = 1
                elif (self.T_e[i, j] < 2):
                    self.Te_matrix_3[i, j] = 1
                elif (self.T_e[i, j] < 5):
                    self.Te_matrix_4[i, j] = 1
                elif (self.T_e[i, j] < 10):
                    self.Te_matrix_5[i, j] = 1
                else:
                    self.Te_matrix_6[i, j] = 1

class Neutral:
    def __init__(self, path, IfT):
        self.path = path
        self.IfT = IfT
        self.k_array=0

    def Neuread(self, K_Rec = 1):
        self.K_Rec = K_Rec
        self.n_D_0 = np.loadtxt(fname=self.path + "data/n_D_0")
        self.n_D2_0 = np.loadtxt(fname=self.path + "data/n_D2_0")
        self.n_D2_1 = np.loadtxt(fname=self.path + "data/n_D2_1")
        self.T_D_0 = np.loadtxt(fname=self.path + "data/T_D_0")
        self.T_D2_0 = np.loadtxt(fname=self.path + "data/T_D2_0")
        
        if(self.k_array):
            self.n_D_0_array = Tn_array(39, self.path + "data/n_D_0_")

        self.recycling_D = np.loadtxt(fname=self.path + "data/recycling_D.txt")
        self.SE_recycling = self.recycling_D[:,0] * self.recycling_D[:,1] * 1.6e-19
        self.Sn_Ion_D_0 = np.loadtxt(fname=self.path + "data/Sn_Ion_D_0")
        self.Sn_CX_D_0 = np.loadtxt(fname=self.path + "data/Sn_CX_D_0")
        if (self.K_Rec):
            self.Sn_Rec_D_1 = np.loadtxt(fname=self.path + "data/Sn_Rec_D_1")
        else:
            self.Sn_Rec_D_1 = np.zeros([98,38])
        self.Sn_Ion_D2_0 = np.loadtxt(fname=self.path + "data/Sn_Ion_D2_0")
        self.Sn_CX_D2_0 = np.loadtxt(fname=self.path + "data/Sn_CX_D2_0")
        self.Sn_Ela_D2_0 = np.loadtxt(fname=self.path + "data/Sn_Ela_D2_0")
        self.Sn_Diss1_D2_0 = np.loadtxt(fname=self.path + "data/Sn_Diss1_D2_0")
        self.Sn_Diss2_D2_0 = np.loadtxt(fname=self.path + "data/Sn_Diss2_D2_0")
        self.Sn_DS1_D2_1 = np.loadtxt(fname=self.path + "data/Sn_DS1_D2_1")
        self.Sn_DS2_D2_1 = np.loadtxt(fname=self.path + "data/Sn_DS2_D2_1")
        self.Sn_DS3_D2_1 = np.loadtxt(fname=self.path + "data/Sn_DS3_D2_1")
        self.Sn_DS4_D2_1 = np.loadtxt(fname=self.path + "data/Sn_DS4_D2_1")

        self.Pra_Ion_D_0 = np.loadtxt(fname = self.path + "data/Pra_Ion_D_0")
        self.Pra_CX_D_0 = np.loadtxt(fname = self.path + "data/Pra_CX_D_0")
        if (self.K_Rec):
            self.Pra_Rec_D_1 = np.loadtxt(fname = self.path + "data/Pra_Rec_D_1")
        else:
            self.Pra_Rec_D_1 = np.zeros([98,38])
        self.SE_Ion_D_0 = np.loadtxt(fname = self.path + "data/SE_Ion_D_0")
        self.SE_CX_D_0 = np.loadtxt(fname = self.path + "data/SE_CX_D_0")
        if (self.K_Rec):
            self.SE_Rec_D_1 = np.loadtxt(fname = self.path + "data/SE_Rec_D_1")
        else:
            self.Pra_Rec_D_1 = np.zeros([98,38])
        self.SE_CX_D2_0 = np.loadtxt(fname = self.path + "data/SE_CX_D2_0")
        self.SE_Ela_D2_0 = np.loadtxt(fname = self.path + "data/SE_Ela_D2_0")

        self.Smu_Ion_D_0 = np.loadtxt(fname = self.path + "data/Smu_Ion_D_0")
        if (self.K_Rec):
            self.Smu_Rec_D_1 = np.loadtxt(fname = self.path + "data/Smu_Rec_D_1")
        else:
            self.Smu_Rec_D_1 = np.zeros([98,38])
        self.Smu_CX_D_0 = np.loadtxt(fname = self.path + "data/Smu_CX_D_0")
        self.Smu1_CX_D_0 = np.loadtxt(fname = self.path + "data/Smu1_CX_D_0")
        self.Smu2_CX_D_0 = np.loadtxt(fname = self.path + "data/Smu2_CX_D_0")
        self.Smu_CX_D2_0 = np.loadtxt(fname = self.path + "data/Smu_CX_D2_0")
        self.Smu_Ela_D2_0 = np.loadtxt(fname = self.path + "data/Smu_Ela_D2_0")
        

        self.Bx = np.loadtxt(fname = self.path + "data/Bx")
        self.By = np.loadtxt(fname = self.path + "data/By")
        self.Bz = np.loadtxt(fname = self.path + "data/Bz")
        self.B = np.sqrt(self.Bx * self.Bx + self.By * self.By + self.Bz * self.Bz)

        self.Vx_D_0 = np.loadtxt(fname = self.path + "data/D_0_Vx")
        self.Vy_D_0 = np.loadtxt(fname = self.path + "data/D_0_Vy")

        self.Vx_D_1 = np.loadtxt(fname = self.path + "data/Vx_D_1")
        self.Vy_D_1 = np.loadtxt(fname = self.path + "data/Vy_D_1")
        self.Vz_D_1 = np.loadtxt(fname = self.path + "data/Vz_D_1")
        self.V_D_1 = V_Magnitude(self.Vx_D_1, self.Vy_D_1, self.Vz_D_1, self.Bx, self.By, self.Bz)

        self.Vx_CX_Ion_Be = np.loadtxt(fname = self.path + "data/Vx_D_1_CX_Be")
        self.Vx_CX_Ion_Af = np.loadtxt(fname = self.path + "data/Vx_D_1_CX_Af")
        self.Vx_CX_Neu_Be = np.loadtxt(fname = self.path + "data/Vx_D_0_CX_Be")
        self.Vx_CX_Neu_Af = np.loadtxt(fname = self.path + "data/Vx_D_0_CX_Af")
        self.Vy_CX_Ion_Be = np.loadtxt(fname = self.path + "data/Vy_D_1_CX_Be")
        self.Vy_CX_Ion_Af = np.loadtxt(fname = self.path + "data/Vy_D_1_CX_Af")
        self.Vy_CX_Neu_Be = np.loadtxt(fname = self.path + "data/Vy_D_0_CX_Be")
        self.Vy_CX_Neu_Af = np.loadtxt(fname = self.path + "data/Vy_D_0_CX_Af")
        self.Vz_CX_Ion_Be = np.loadtxt(fname = self.path + "data/Vz_D_1_CX_Be")
        self.Vz_CX_Ion_Af = np.loadtxt(fname = self.path + "data/Vz_D_1_CX_Af")
        self.Vz_CX_Neu_Be = np.loadtxt(fname = self.path + "data/Vz_D_0_CX_Be")
        self.Vz_CX_Neu_Af = np.loadtxt(fname = self.path + "data/Vz_D_0_CX_Af")
        self.V_CX_Ion_Be = V_Magnitude(self.Vx_CX_Ion_Be, self.Vy_CX_Ion_Be, self.Vz_CX_Ion_Be, self.Bx, self.By, self.Bz)
        self.V_CX_Ion_Af = V_Magnitude(self.Vx_CX_Ion_Af, self.Vy_CX_Ion_Af, self.Vz_CX_Ion_Af, self.Bx, self.By, self.Bz)
        self.V_CX_Neu_Be = V_Magnitude(self.Vx_CX_Neu_Be, self.Vy_CX_Neu_Be, self.Vz_CX_Neu_Be, self.Bx, self.By, self.Bz)
        self.V_CX_Neu_Af = V_Magnitude(self.Vx_CX_Neu_Af, self.Vy_CX_Neu_Af, self.Vz_CX_Neu_Af, self.Bx, self.By, self.Bz)

        self.Vx_CX_Neu_Delta = self.Vx_CX_Neu_Af - self.Vx_CX_Neu_Be
        self.Vy_CX_Neu_Delta = self.Vy_CX_Neu_Af - self.Vy_CX_Neu_Be
        self.Vz_CX_Neu_Delta = self.Vz_CX_Neu_Af - self.Vz_CX_Neu_Be
        self.V_CX_Neu_Delta = V_Magnitude(self.Vx_CX_Neu_Af, self.Vy_CX_Neu_Af, self.Vz_CX_Neu_Af, self.Bx, self.By, self.Bz)
        self.Vx_CX_Ion_Delta = self.Vx_CX_Ion_Af - self.Vx_CX_Ion_Be
        self.Vy_CX_Ion_Delta = self.Vy_CX_Ion_Af - self.Vy_CX_Ion_Be
        self.Vz_CX_Ion_Delta = self.Vz_CX_Ion_Af - self.Vz_CX_Ion_Be
        self.V_CX_Ion_Delta = V_Magnitude(self.Vx_CX_Ion_Af, self.Vy_CX_Ion_Af, self.Vz_CX_Ion_Af, self.Bx, self.By, self.Bz)


        self.Vx_Ua = np.loadtxt(fname = self.path + "data/Vx_Ua")
        self.Vy_Ua = np.loadtxt(fname = self.path + "data/Vy_Ua")
        self.Vz_Ua = np.loadtxt(fname = self.path + "data/Vz_Ua")
        self.V_Ua = V_Magnitude(self.Vx_Ua, self.Vy_Ua, self.Vz_Ua, self.Bx, self.By, self.Bz)

        
        self.Smu_Ion_D_0_Ua = np.zeros([98, 38])
        self.Smu_Rec_D_1_Ua = np.zeros([98, 38])
        self.Smu_CX_D_0_Ua = np.zeros([98, 38])
        self.Smu_CX_D2_0_Ua = np.zeros([98, 38])
        self.Smu_Ela_D2_0_Ua = np.zeros([98, 38])

        if (self.IfT):
            self.n_T_0 = np.loadtxt(fname=self.path + "data/n_T_0")
            self.n_T2_0 = np.loadtxt(fname=self.path + "data/n_T2_0")
            self.T_T_0 = np.loadtxt(fname=self.path + "data/T_T_0")
            self.T_T2_0 = np.loadtxt(fname=self.path + "data/T_T2_0")
        
        
        self.Sn_D_0 = -self.Sn_Ion_D_0 + self.Sn_Rec_D_1 + 2 * self.Sn_Diss1_D2_0 + self.Sn_Diss2_D2_0 + self.Sn_CX_D2_0 + self.Sn_DS2_D2_1 + 2 * self.Sn_DS3_D2_1
        self.Sn_D_1 = self.Sn_Ion_D_0 - self.Sn_Rec_D_1 + self.Sn_Diss2_D2_0 - self.Sn_CX_D2_0 + 2 * self.Sn_DS1_D2_1 + self.Sn_DS2_D2_1
        if(self.k_array):
            self.Sn_D_0_array = [np.sum(-self.Sn_Ion_D_0[1:24,:]), np.sum(self.Sn_Rec_D_1[1:24,:]), 2 * np.sum(self.Sn_Diss1_D2_0[1:24,:]), np.sum(self.Sn_Diss2_D2_0[1:24,:]), np.sum(self.Sn_CX_D2_0[1:24,:]), np.sum(self.Sn_DS2_D2_1[1:24,:]), 2 * np.sum(self.Sn_DS3_D2_1[1:24,:])]
            self.Sn_D_1_array = [np.sum(self.Sn_Ion_D_0[1:24,:]), np.sum(-self.Sn_Rec_D_1[1:24,:]), np.sum(self.Sn_Diss2_D2_0[1:24,:]), np.sum(-self.Sn_CX_D2_0[1:24,:]), 2 * np.sum(self.Sn_DS1_D2_1[1:24,:]), np.sum(self.Sn_DS2_D2_1[1:24,:])]

        self.Sn_Ion_D_0_1D_Inner = np.transpose(self.Sn_Ion_D_0[1, 1:37].reshape(1,36))
        self.Sn_CX_D_0_1D_Inner = np.transpose(self.Sn_CX_D_0[1, 1:37].reshape(1,36))
        self.Sn_Ion_D_0_1D_Outer = np.transpose(self.Sn_Ion_D_0[96, 1:37].reshape(1,36))
        self.Sn_CX_D_0_1D_Outer = np.transpose(self.Sn_CX_D_0[96, 1:37].reshape(1,36))

    def Smu_array(self, SOLPScase):
        self.Smu_Ion_D_0_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_Ion_D_0)
        self.Smu_Rec_D_1_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_Rec_D_1)
        self.Smu_CX_D_0_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_CX_D_0)
        self.Smu_CX_D2_0_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_CX_D2_0)
        self.Smu_Ela_D2_0_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_Ela_D2_0)
        self.Smu_Ion_D_0_Ua_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_Ion_D_0_Ua)
        self.Smu_Rec_D_1_Ua_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_Rec_D_1_Ua)
        self.Smu_CX_D_0_Ua_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_CX_D_0_Ua)
        self.Smu_CX_D2_0_Ua_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_CX_D2_0_Ua)
        self.Smu_Ela_D2_0_Ua_array = S_Te_array(SOLPScase.Te_matrix, self.Smu_Ela_D2_0_Ua)
        self.Te_array = S_Te_array(SOLPScase.Te_matrix, SOLPScase.T_e)


    def SOLPSMu(self, SOLPScase):
        self.Smu_Ion_D_0_Ua = self.Sn_Ion_D_0 * mass_D * SOLPScase.Ua_Dion_2D
        self.Smu_Rec_D_1_Ua = -self.Sn_Rec_D_1 * mass_D * SOLPScase.Ua_Dion_2D
        self.Smu_CX_D_0_Ua = -self.Sn_CX_D_0 * mass_D * SOLPScase.Ua_Dion_2D
        self.Smu_CX_D2_0_Ua = -self.Sn_CX_D2_0 * mass_D * SOLPScase.Ua_Dion_2D
        self.Smu_Ela_D2_0_Ua = -self.Sn_Ela_D2_0 * mass_D * SOLPScase.Ua_Dion_2D

    def Smu_sum(self):
        self.Smu_D = self.Smu_CX_D_0 + self.Smu_Ion_D_0 + self.Smu_Rec_D_1
        self.Smu_D2 = self.Smu_CX_D2_0 + self.Smu_Ela_D2_0
        self.Smu = self.Smu_D + self.Smu_D2
           
    def Smu_fix(self, EAST_SOLPS):
        for i in range(98):
            for j in range(38):
                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_Ion_D_0[i, j] > 0):
                    self.Smu_Ion_D_0[i, j] = -1 * abs(self.Smu_Ion_D_0[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu_Ion_D_0[i, j] = abs(self.Smu_Ion_D_0[i, j]) * EAST_SOLPS.Vol[i, j]

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_Rec_D_1[i, j] > 0):
                    self.Smu_Rec_D_1[i, j] = -1 * abs(self.Smu_Rec_D_1[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu_Rec_D_1[i, j] = abs(self.Smu_Rec_D_1[i, j]) * EAST_SOLPS.Vol[i, j]

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_CX_D_0[i, j] > 0):
                    self.Smu_CX_D_0[i, j] = -1 * abs(self.Smu_CX_D_0[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu_CX_D_0[i, j] = abs(self.Smu_CX_D_0[i, j]) * EAST_SOLPS.Vol[i, j]

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu1_CX_D_0[i, j] > 0):
                    self.Smu1_CX_D_0[i, j] = -1 * abs(self.Smu1_CX_D_0[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu1_CX_D_0[i, j] = abs(self.Smu1_CX_D_0[i, j]) * EAST_SOLPS.Vol[i, j]
                    
                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu2_CX_D_0[i, j] > 0):
                    self.Smu2_CX_D_0[i, j] = -1 * abs(self.Smu2_CX_D_0[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu2_CX_D_0[i, j] = abs(self.Smu2_CX_D_0[i, j]) * EAST_SOLPS.Vol[i, j]

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_CX_D2_0[i, j] > 0):
                    self.Smu_CX_D2_0[i, j] = -1 * abs(self.Smu_CX_D2_0[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu_CX_D2_0[i, j] = abs(self.Smu_CX_D2_0[i, j]) * EAST_SOLPS.Vol[i, j]

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_Ela_D2_0[i, j] > 0):
                    self.Smu_Ela_D2_0[i, j] = -1 * abs(self.Smu_Ela_D2_0[i, j]) * EAST_SOLPS.Vol[i, j]
                else:
                    self.Smu_Ela_D2_0[i, j] = abs(self.Smu_Ela_D2_0[i, j]) * EAST_SOLPS.Vol[i, j]
                    


                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_Ion_D_0_Ua[i, j] > 0):
                    self.Smu_Ion_D_0_Ua[i, j] = -1 * abs(self.Smu_Ion_D_0_Ua[i, j])
                else:
                    self.Smu_Ion_D_0_Ua[i, j] = abs(self.Smu_Ion_D_0_Ua[i, j])

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_Rec_D_1_Ua[i, j] > 0):
                    self.Smu_Rec_D_1_Ua[i, j] = -1 * abs(self.Smu_Rec_D_1_Ua[i, j])
                else:
                    self.Smu_Rec_D_1_Ua[i, j] = abs(self.Smu_Rec_D_1_Ua[i, j])

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_CX_D_0_Ua[i, j] > 0):
                    self.Smu_CX_D_0_Ua[i, j] = -1 * abs(self.Smu_CX_D_0_Ua[i, j])
                else:
                    self.Smu_CX_D_0_Ua[i, j] = abs(self.Smu_CX_D_0_Ua[i, j])

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_CX_D2_0_Ua[i, j] > 0):
                    self.Smu_CX_D2_0_Ua[i, j] = -1 * abs(self.Smu_CX_D2_0_Ua[i, j])
                else:
                    self.Smu_CX_D2_0_Ua[i, j] = abs(self.Smu_CX_D2_0_Ua[i, j])

                if (EAST_SOLPS.Ua_Dion_2D[i, j] * self.Smu_Ela_D2_0_Ua[i, j] > 0):
                    self.Smu_Ela_D2_0_Ua[i, j] = -1 * abs(self.Smu_Ela_D2_0_Ua[i, j])
                else:
                    self.Smu_Ela_D2_0_Ua[i, j] = abs(self.Smu_Ela_D2_0_Ua[i, j])

    def set_K_Rec(self):
        self.K_Rec
def axset(ax_, p1, wall, ax_title, xlim, ylim, clim):
    ax_.add_collection(p1)
    ax_.autoscale_view()
    ax_.plot(wall[:, 0], wall[:, 1], color="b", lw=0.1)
    ax_.axis('equal')
    ax_.set_title(ax_title)
    ax_.set_xlim(xlim)
    ax_.set_ylim(ylim)
    p1.set_clim(clim)

def Divertor(value):
    return np.sum(value[1:24,1:37])+np.sum(value[74:97,1:37])
def OuterDivertor(value):
    return np.sum(value[74:97,1:37])
def InnerDivertor(value):
    return np.sum(value[1:24,1:37])

def V_norm(U,V):
    return [U / np.sqrt(U*U+V*V), V / np.sqrt(U*U+V*V)] * 2

def V_Magnitude(Vx, Vy, Vz, Bx, By, Bz):
    return (Vx * Bx + Vy * By + Vz * Bz)/1

def plot_scan(value, Value_label, file_name, Legend, cclim, iflog = 1):
    with PdfPages(save_figures_path + file_name) as pdf:
        patches = []
        ixs = 1
        com_nx = 96
        com_ny = 36
        ixe = com_nx+1
        iys = 1
        iye = com_ny+1
        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                rcol = com_rm[ix, iy, [1, 2, 4, 3]]
                zcol = com_zm[ix, iy, [1, 2, 4, 3]]
                rcol.shape = (4, 1)
                zcol.shape = (4, 1)
                polygon = Polygon(np.column_stack((rcol, zcol)), True)
                patches.append(polygon)

        vals1 = np.zeros((ixe-ixs)*(iye-iys))
        vals2 = np.zeros((ixe-ixs)*(iye-iys))
        vals3 = np.zeros((ixe-ixs)*(iye-iys))

        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                k = (ix-ixs)+(ixe-ixs)*(iy-iys)
                vals1[k] = value[0][ix, iy]
                vals2[k] = value[1][ix, iy]
                vals3[k] = value[2][ix, iy]
        if (iflog):
            p1 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
            p1.set_array(np.array(vals1))
            p2 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
            p2.set_array(np.array(vals2))
            p3 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
            p3.set_array(np.array(vals3))
        else:
            p1 = PatchCollection(patches, cmap=mpl.cm.jet)
            p1.set_array(np.array(vals1))
            p2 = PatchCollection(patches, cmap=mpl.cm.jet)
            p2.set_array(np.array(vals2))
            p3 = PatchCollection(patches, cmap=mpl.cm.jet)
            p3.set_array(np.array(vals3))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))
        axset(ax[0], p1, wall, Legend[0] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[1], p2, wall, Legend[1] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[2], p3, wall, Legend[2] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        fig.colorbar(p1, ax=ax[0])
        fig.colorbar(p2, ax=ax[1])
        fig.colorbar(p3, ax=ax[2])
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        plt.close()
        
        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                k = (ix-ixs)+(ixe-ixs)*(iy-iys)
                vals1[k] = value[3][ix, iy]
                vals2[k] = value[4][ix, iy]
                vals3[k] = value[8][ix, iy]
        if (iflog):
            p1 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
            p1.set_array(np.array(vals1))
            p2 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
            p2.set_array(np.array(vals2))
            p3 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
            p3.set_array(np.array(vals3))
        else:
            p1 = PatchCollection(patches, cmap=mpl.cm.jet)
            p1.set_array(np.array(vals1))
            p2 = PatchCollection(patches, cmap=mpl.cm.jet)
            p2.set_array(np.array(vals2))
            p3 = PatchCollection(patches, cmap=mpl.cm.jet)
            p3.set_array(np.array(vals3))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))
        axset(ax[0], p1, wall, Legend[3] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[1], p2, wall, Legend[4] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[2], p3, wall, Legend[8] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        fig.colorbar(p1, ax=ax[0])
        fig.colorbar(p2, ax=ax[1])
        fig.colorbar(p3, ax=ax[2])
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        plt.close()
        

        # Target n_D plot
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(7, 15))
        if (iflog):
            ax[0].semilogy(x, value[0][63, 1:37], lw=0.1)
            #ax[0].semilogy(x, value[1][63, 1:37], lw=0.1)
            #ax[0].semilogy(x, value[2][63, 1:37], lw=0.1)
            #ax[0].semilogy(x, value[3][63, 1:37], lw=0.1)
            ax[0].semilogy(x, value[4][63, 1:37], lw=0.1)
            ax[0].semilogy(x, value[8][63, 1:37], lw=0.1)
        else:
            ax[0].plot(x, value[0][63, 1:37], lw=0.1)
            #ax[0].plot(x, value[1][63, 1:37], lw=0.1)
            #ax[0].plot(x, value[2][63, 1:37], lw=0.1)
            #ax[0].plot(x, value[3][63, 1:37], lw=0.1)
            ax[0].plot(x, value[4][63, 1:37], lw=0.1)
            ax[0].plot(x, value[8][63, 1:37], lw=0.1)
        #ax[0].semilogy(x, value[5][63, 1:37], lw=0.1)
        #ax[0].semilogy(x, value[6][63, 1:37], lw=0.1)
        #ax[0].semilogy(x, value[7][63, 1:37], lw=0.1)
        ax[0].set_title(Value_label)
        ax[0].set_xlabel("Radial Index")

        if (iflog):
            ax[1].semilogy(x_l[1:37], value[0][1, 1:37], lw=0.1)
            #ax[1].semilogy(x_l[1:37], value[1][1, 1:37], lw=0.1)
            #ax[1].semilogy(x_l[1:37], value[2][1, 1:37], lw=0.1)
            #ax[1].semilogy(x_l[1:37], value[3][1, 1:37], lw=0.1)
            ax[1].semilogy(x_l[1:37], value[4][1, 1:37], lw=0.1)
            ax[1].semilogy(x_l[1:37], value[8][1, 1:37], lw=0.1)
        else:
            ax[1].plot(x_l[1:37], value[0][1, 1:37], lw=0.1)
            #ax[1].plot(x_l[1:37], value[1][1, 1:37], lw=0.1)
            #ax[1].plot(x_l[1:37], value[2][1, 1:37], lw=0.1)
            #ax[1].plot(x_l[1:37], value[3][1, 1:37], lw=0.1)
            ax[1].plot(x_l[1:37], value[4][1, 1:37], lw=0.1)
            ax[1].plot(x_l[1:37], value[8][1, 1:37], lw=0.1)
        #ax[1].semilogy(x_l[1:37], value[5][1, 1:37], lw=0.1)
        #ax[1].semilogy(x_l[1:37], value[6][1, 1:37], lw=0.1)
        #ax[1].semilogy(x_l[1:37], value[7][1, 1:37], lw=0.1)
        ax[1].set_title(Value_label)
        ax[1].set_xlabel("${r-r}_{sep}$ at Inner Target(cm)")

        if (iflog):
            ax[2].semilogy(x_r[1:37], value[0][96, 1:37],lw=0.1)
            #ax[2].semilogy(x_r[1:37], value[1][96, 1:37],lw=0.1)
            #ax[2].semilogy(x_r[1:37], value[2][96, 1:37],lw=0.1)
            #ax[2].semilogy(x_r[1:37], value[3][96, 1:37],lw=0.1)
            ax[2].semilogy(x_r[1:37], value[4][96, 1:37],lw=0.1)
            ax[2].semilogy(x_r[1:37], value[8][96, 1:37],lw=0.1)
        else:
            ax[2].plot(x_r[1:37], value[0][96, 1:37],lw=0.1)
            #ax[2].plot(x_r[1:37], value[1][96, 1:37],lw=0.1)
            #ax[2].plot(x_r[1:37], value[2][96, 1:37],lw=0.1)
            #ax[2].plot(x_r[1:37], value[3][96, 1:37],lw=0.1)
            ax[2].plot(x_r[1:37], value[4][96, 1:37],lw=0.1)
            ax[2].plot(x_r[1:37], value[8][96, 1:37],lw=0.1)
        #ax[2].semilogy(x_r[1:37], value[5][96, 1:37],lw=0.1)
        #ax[2].semilogy(x_r[1:37], value[6][96, 1:37],lw=0.1)
        #ax[2].semilogy(x_r[1:37], value[7][96, 1:37],lw=0.1)
        ax[2].set_title(Value_label)
        ax[2].set_xlabel("${r-r}_{sep}$ at Outer Target(cm)")
        ax[2].legend(["3.2e19", "5e19", "5.8e19"], loc='best')
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        #plt.show()
        plt.close()

        fig,ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 7))
        ax[0].semilogy(dp_l[0:9], value[0][1:10,20])
        ax[0].semilogy(dp_l[0:9], value[4][1:10,20])
        ax[0].semilogy(dp_l[0:9], value[8][1:10,20])
        ax[1].semilogy(dp_r[0:9], value[0][88:97,20][::-1])
        ax[1].semilogy(dp_r[0:9], value[4][88:97,20][::-1])
        ax[1].semilogy(dp_r[0:9], value[8][88:97,20][::-1])
        ax[0].set_title(Value_label)
        ax[0].legend(["3.2e19", "5e19", "5.8e19"], loc='best')
        ax[0].set_xlabel("Distance to Inner Target(m)")
        ax[1].set_xlabel("Distance to Outer Target(m)")
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        #plt.show()
        plt.close()

def row_scan(value1, value2, Value_magntude, Value_label, file_name, Legend, cclim):
    with PdfPages(save_figures_path + file_name) as pdf:
        patches = []
        ixs = 1
        com_nx = 96
        com_ny = 36
        ixe = com_nx+1
        iys = 1
        iye = com_ny+1
        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                rcol = com_rm[ix, iy, [1, 2, 4, 3]]
                zcol = com_zm[ix, iy, [1, 2, 4, 3]]
                rcol.shape = (4, 1)
                zcol.shape = (4, 1)
                polygon = Polygon(np.column_stack((rcol, zcol)), True)
                patches.append(polygon)

        vals1 = np.zeros((ixe-ixs)*(iye-iys))
        vals2 = np.zeros((ixe-ixs)*(iye-iys))
        vals3 = np.zeros((ixe-ixs)*(iye-iys))

        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                k = (ix-ixs)+(ixe-ixs)*(iy-iys)
                vals1[k] = Value_magntude[0][ix, iy]
                vals2[k] = Value_magntude[1][ix, iy]
                vals3[k] = Value_magntude[2][ix, iy]

        p1 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p1.set_array(np.array(vals1))
        p2 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p2.set_array(np.array(vals2))
        p3 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p3.set_array(np.array(vals3))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))
        axset(ax[0], p1, wall, Legend[0] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[1], p2, wall, Legend[1] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[2], p3, wall, Legend[2] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        for i in range(3):
            ax[i].quiver(location_x, location_y, V_norm(value1[i], value2[i])[0], V_norm(value1[i], value2[i])[1])
        fig.colorbar(p1, ax=ax[0])
        fig.colorbar(p2, ax=ax[1])
        fig.colorbar(p3, ax=ax[2])
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        plt.close()
        
        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                k = (ix-ixs)+(ixe-ixs)*(iy-iys)
                vals1[k] = Value_magntude[3][ix, iy]
                vals2[k] = Value_magntude[4][ix, iy]
                vals3[k] = Value_magntude[5][ix, iy]
        p1 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p1.set_array(np.array(vals1))
        p2 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p2.set_array(np.array(vals2))
        p3 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p3.set_array(np.array(vals3))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))
        axset(ax[0], p1, wall, Legend[3] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[1], p2, wall, Legend[4] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[2], p3, wall, Legend[8] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        case = [3,4,8]
        for j in range(3):
            ax[j].quiver(location_x, location_y, V_norm(value1[case[j]], value2[case[j]])[0], V_norm(value1[case[j]], value2[case[j]])[1])
        fig.colorbar(p1, ax=ax[0])
        fig.colorbar(p2, ax=ax[1])
        fig.colorbar(p3, ax=ax[2])
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        plt.close()

def plot_array(value, Value_label, file_name, Legend, cclim):
    with PdfPages(save_figures_path + file_name) as pdf:
        patches = []
        ixs = 1
        com_nx = 96
        com_ny = 36
        ixe = com_nx+1
        iys = 1
        iye = com_ny+1
        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                rcol = com_rm[ix, iy, [1, 2, 4, 3]]
                zcol = com_zm[ix, iy, [1, 2, 4, 3]]
                rcol.shape = (4, 1)
                zcol.shape = (4, 1)
                polygon = Polygon(np.column_stack((rcol, zcol)), True)
                patches.append(polygon)

        vals1 = np.zeros((ixe-ixs)*(iye-iys))
        vals2 = np.zeros((ixe-ixs)*(iye-iys))
        vals3 = np.zeros((ixe-ixs)*(iye-iys))

        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                k = (ix-ixs)+(ixe-ixs)*(iy-iys)
                vals1[k] = value[0][ix, iy]
                vals2[k] = value[1][ix, iy]
                vals3[k] = value[2][ix, iy]

        p1 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p1.set_array(np.array(vals1))
        p2 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p2.set_array(np.array(vals2))
        p3 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p3.set_array(np.array(vals3))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))
        axset(ax[0], p1, wall, Legend[0] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[1], p2, wall, Legend[1] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[2], p3, wall, Legend[2] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        fig.colorbar(p1, ax=ax[0])
        fig.colorbar(p2, ax=ax[1])
        fig.colorbar(p3, ax=ax[2])
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        plt.close()
        
        for iy in np.arange(iys, iye):
            for ix in np.arange(ixs, ixe):
                k = (ix-ixs)+(ixe-ixs)*(iy-iys)
                vals1[k] = value[3][ix, iy]
                vals2[k] = value[4][ix, iy]
                vals3[k] = value[5][ix, iy]
        p1 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p1.set_array(np.array(vals1))
        p2 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p2.set_array(np.array(vals2))
        p3 = PatchCollection(patches, cmap=mpl.cm.jet, norm=colors.LogNorm())
        p3.set_array(np.array(vals3))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))
        axset(ax[0], p1, wall, Legend[3] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[1], p2, wall, Legend[4] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        axset(ax[2], p3, wall, Legend[5] + "  " + Value_label, xlim_EAST, ylim_EAST, cclim)
        fig.colorbar(p1, ax=ax[0])
        fig.colorbar(p2, ax=ax[1])
        fig.colorbar(p3, ax=ax[2])
        pdf.savefig()    # saves the current figure into a pdf page
        plt.draw()
        plt.close()