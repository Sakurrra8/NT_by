import matplotlib as mpl 
import seaborn as sns

def load_eirene_mesh(self, baserun_path):
    """Load EIRENE nodes from the fort.33 file and triangulation from the fort.34 file"""
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

    return xnodes, ynodes, triangles - 1

 def plot2d_eirene(self,case_path,vals,ax=None,title='',label="",DDN=False,lb=None,ub=None,cb_scale='log',clevels=None,cdict=sns.color_palette('hls',8),cmap=mpl.cm.jet,cbar_show=True,cax=None,xlabel="R (m)",ylabel="Z (m)", **kwargs):
        """Method to plot 2D fields from EIRENE.

        Parameters
        ----------
        vals : array (self.triangles)
            Data array for an EIRENE variable of interest.
        ax : matplotlib Axes instance
            Axes on which to plot. If left to None, a new figure is created.
        scale : str
            Choice of 'linear','log' and 'symlog' for matplotlib.colors.
        label : str
            Label to set on the colorbar. No label by default.
        lb : float
            Lower bound for colorbar. If left to None, the minimum value in `vals` is used.
        ub : float
            Upper bound for colorbar. If left to None, the maximum value in `vals` is used.
        replace_zero : boolean
            If True (default), replace all zeros in 'vals' with minimum value in 'vals'
        kwargs
            Additional keyword arguments passed to the `tripcolor` function.

        """
        # Avoid zeros that may derive from low Monte Carlo statistics
        # np.nan_to_num(vals,copy=False)

        
        vals = vals.flatten()  
        vals1=vals[~np.isnan(vals)][~np.isinf(vals)]
        if cb_scale=='log':
            vals[vals==0]= np.min(vals[vals!=0])
        if lb is None:
            lb = np.nanmin(vals)
        if ub is None:
            ub = np.nanmax(vals)

        if lb <=0 and cb_scale=='log':
            lb=1e-40
        if clevels is None and cb_scale=="manual":
            if isinstance(cdict,list):
                y_step=int((ub-lb)/len(cdict))
                clevels=list(np.arange(lb,ub,y_step).astype('float'))
            else:
                return ('LinearSegmentedColormap has no len(),please set clevels')
        else:
            clevels=clevels

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 11))

        

        # given quantity is on EIRENE triangulation
        if cb_scale=='linear':
            norm = mpl.colors.Normalize(vmin=lb, vmax=ub)
        elif cb_scale=='log':
            norm = mpl.colors.LogNorm(vmin=lb,vmax=ub)
            
        elif cb_scale=='symlog':
            norm = mpl.colors.SymLogNorm(linthresh=ub/10.,base=10, lincb_scale=0.5, vmin=lb,vmax=ub)
        elif cb_scale=='manual':
            if isinstance(cdict,list):
                cmap=mpl.colors.ListedColormap(cdict)
            else:
                cmap=cdict
            norm=mpl.colors.BoundaryNorm(clevels, cmap.N)
            
            
        else:
            raise ValueError('Unrecognized cb_scale parameter')
        xnodes,ynodes,triangles = self.load_eirene_mesh(case_path)
        cntr = ax.tripcolor(xnodes,ynodes,triangles,cmap=cmap,facecolors=vals.flatten(),norm=norm,**kwargs,)
        ax.axis("scaled")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        b2fgmtry_path=case_path+r'/'+self.__b2fgmtry_data_path
        if Path(b2fgmtry_path).exists():
            gem_data=self.read_gem_data(b2fgmtry_path)
            crx_data=gem_data['crx'].astype(np.float64)
            cry_data=gem_data['cry'].astype(np.float64)
            NX=gem_data['nx']+2
            NY=gem_data['ny']+2
            xx=crx_data.transpose(2,1,0)
            yy=cry_data.transpose(2,1,0)
            topcut=gem_data['topcut']
            leftcut=gem_data['leftcut']
            rightcut=gem_data['rightcut']
            bottomcut=gem_data['bottomcut']
        if not DDN:
            ax.plot(crx_data[0,topcut+1,:],cry_data[0,topcut+1,:],crx_data[1,topcut+1,-1],cry_data[1,topcut+1,-1],color='purple',lw=1.0)
        else:
            print('The separatrix can not been plotted because index is wrong')
        if cbar_show:
            if cax is None:
                cbar = plt.colorbar(cntr, ax=ax,label=label)  # format='%.3g'
            else:
                cbar = plt.colorbar(cntr, cax=cax,label=label)  # format='%.3g'
        ax.set_title(title)
        return cntr