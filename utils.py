import pandas as pd
import numpy as np
import string
from matplotlib.transforms import ScaledTranslation
from matplotlib import pyplot as plt

def print_setup():
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["savefig.dpi"] = 300

def computeRelativeGrowthRate(df,switch_time,norm_time=None,make_plot=True,cutoff = 0.002):
    # by default, use 50 minutes before media switch to compute average growth rate
    if norm_time is None:
        norm_time = max(0,switch_time-50)

    mean_gr_df = df[(df['Time (min)']>=norm_time) & (df['Time (min)']<=switch_time)].groupby(['Trap','Strain'],as_index=False)['GrowthRate'].mean()
    mean_gr_df=mean_gr_df.set_index('Trap')
    
    if make_plot:
        mean_gr_df.hist(by="Strain",column="GrowthRate",sharex=True,sharey=True)

    #discard outliers
    mean_gr_df.loc[(mean_gr_df['GrowthRate']<cutoff), 'GrowthRate'] = np.nan
    for strain in mean_gr_df['Strain'].unique():
        print(f"{strain}: discarded {mean_gr_df.loc[mean_gr_df['Strain']==strain,'GrowthRate'].isna().sum()} traps")

    #compute normalized growth rate
    df['Relative Growth Rate'] = df.apply(lambda x: x['GrowthRate'] / mean_gr_df.loc[x['Trap'],'GrowthRate'], axis=1)
    df['GrowthRate']= df.apply(lambda x: np.nan if np.isnan(x['Relative Growth Rate']) else x['GrowthRate'],axis=1)
    return df


def add_antibiotic_info(ax,switch_timepoints,labels,color="black",lw=1,yfrac=0.9):
    ylim = ax.get_ylim()
    fontweight="bold" if lw>1 else "normal"
    for t,l in zip(switch_timepoints,labels):
        ax.vlines(t,ylim[0],ylim[1],colors='black',linestyles="dashed",color=color,lw=lw)
        ax.text(t,ylim[0]+yfrac*(ylim[1]-ylim[0])," + "+l,color=color,fontweight=fontweight)
    ax.set_ylim(ylim)

def italic_strain(strain):
    parts = strain.split(",")
    if len(parts)==2:
        label = rf"$\it{{{parts[0]}}}$, {parts[1]}"
    else:
        label = rf"$\it{{{parts[0]}}}$"
    return label

def make_italic_legend(ax):
    h,l=ax.get_legend_handles_labels()
    ax.legend(h,[italic_strain(li) for li in l])

def label_subplots(fig):
    axs = fig.get_axes()
    labels = list(string.ascii_lowercase[:len(axs)])
    for label, ax in zip(labels,axs):
    # Use ScaledTranslation to put the label
    # - at the top left corner (axes fraction (0, 1)),
    # - offset 20 pixels left and 7 pixels up (offset points (-20, +7)),
    # i.e., just outside the axes.
        ax.text(
            0.0, 1.0, label, transform=(
                ax.transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
            fontsize='large', va='bottom', weight="bold")#, fontfamily='serifbold')
            
def convertTimeToKymoCoords(t,im,num_traps,dt,t_start=0):
    trap_width = im.shape[1]/num_traps
    w = trap_width/2
    t = np.array(t)
    ticks = t*trap_width/dt+w
    ticklabels = t_start+t
    return ticks,ticklabels
    
def plotKymograph(ax,im,num_traps,dt,t=None,t_start=0,labelstep=50):
    ax.imshow(im,cmap="gray")
    if t is None:
        trap_width = im.shape[1]/num_traps
        w = trap_width/2
        tmax = (num_traps-1)*dt
        t = np.arange(0,tmax,labelstep)
    ticks, ticklabels = convertTimeToKymoCoords(t,im,num_traps,dt,t_start)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.set_yticks([])
    ax.set_xlabel("Time (min)") 

def frame2time(frames,dt):
    if type(frames)==list:
        frames = [(f-1)*dt for f in frames]
    else:
        frames = (frames-1)*dt
    return frames

def decodeTraps(traps,num_traps):
    if not type(traps)==list:
        traps = list(traps)
    return [(t//num_traps+1, t % num_traps) for t in traps]


