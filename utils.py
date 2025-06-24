import pandas as pd
import numpy as np
import string
from matplotlib.transforms import ScaledTranslation
from matplotlib import pyplot as plt

def print_setup():
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["savefig.dpi"] = 300

def computeRelativeGrowthRate_v0(df,switch_frame,norm_frame=[]):
    if norm_frame==[]:
        norm_frame=max(1,switch_frame-40)

    growth_rate_before_switch_df = df[(df['Frame']>=norm_frame) & (df['Frame']<=switch_frame)].groupby(['Trap','Strain'],as_index=False)['GrowthRate'].mean()

    #discard outliers
    #means = growth_rate_before_switch_df.groupby('Strain')['GrowthRate'].mean()
    #stdevs = growth_rate_before_switch_df.groupby('Strain')['GrowthRate'].std()
    #cutoffs = means-3*stdevs
    cutoff=0.002
    #print(cutoffs)
    for strain in growth_rate_before_switch_df['Strain'].unique():
        growth_rate_before_switch_df.loc[(growth_rate_before_switch_df['GrowthRate']<cutoff) & 
                                        (growth_rate_before_switch_df['Strain']==strain), 'GrowthRate'] = np.nan
        #print(f"{strain}: discarded {growth_rate_before_switch_df.loc[growth_rate_before_switch_df['Strain']==strain,'GrowthRate'].isna().sum()} traps")


    #compute normalized growth rate
    growth_rate_before_switch_df=growth_rate_before_switch_df.set_index('Trap')
    #print(growth_rate_before_switch_df[10])
    df['Relative Growth Rate'] = df.apply(lambda x: x['GrowthRate'] / growth_rate_before_switch_df.loc[x['Trap'],'GrowthRate'], axis=1)
    return df


def computeRelativeGrowthRate(df,switch_frame,norm_frame=None,make_plot=False,cutoff = 0.002):
    if norm_frame is None:
        growth_rate_before_switch_df = pd.DataFrame()
        if make_plot:
            fig,axs=plt.subplots(1,len(df['Strain'].unique()),sharex=True,sharey=True,layout="constrained")
        else:
            axs = np.arange(len(set(df['Strain'])))
            
        for (i,s) in enumerate(set(df['Strain'])):
            df_s=df[df['Strain']==s]
            traps = df_s['Trap'].unique()
            num_frames = np.arange(20, min(70,switch_frame))
            norm_factors=np.empty((len(traps),len(num_frames)),)*np.nan
            for (j,num_frame) in enumerate(num_frames):
                growth_rate_before_switch_df_j = df_s[(df_s['Frame']>=switch_frame-num_frame) & (df_s['Frame']<=switch_frame)].groupby('Trap')['GrowthRate'].mean()
                norm_factors[:,j] = growth_rate_before_switch_df_j.values
            norm_factors_df = pd.DataFrame(norm_factors,columns=num_frames,index=traps)
            rel_norm_factors_df = norm_factors_df.apply(lambda x: x/x.iloc[-1],axis=1)
            meanvals = rel_norm_factors_df.mean()
            stdvals = rel_norm_factors_df.std()
            ind=np.argmax((np.abs(meanvals-1)<0.001) & (stdvals<0.01))
            #print(meanvals)
            #print(ind)

            print(f"Num frames for growth rate normalization of {s}: {num_frames[ind]}")
            growth_rate_before_switch_df_s = norm_factors_df[num_frames[ind]].to_frame()
            growth_rate_before_switch_df_s.rename(columns={num_frames[ind]:"GrowthRate"},inplace=True)
            growth_rate_before_switch_df_s.loc[:,'Strain'] = s
            growth_rate_before_switch_df = pd.concat([growth_rate_before_switch_df,growth_rate_before_switch_df_s])
            if make_plot:
                if len(set(df['Strain']))>1:
                    ax=axs[i]
                else:
                    ax=axs
                ax.plot(num_frames,meanvals)
                ax.fill_between(num_frames,meanvals+stdvals,meanvals-stdvals,alpha=0.25)
                ax.set_title(s)
        if make_plot:
            plt.show()

        growth_rate_before_switch_df = growth_rate_before_switch_df.reset_index().rename({"index":"Trap"},axis="columns")
    else:
        growth_rate_before_switch_df = df[(df['Frame']>=norm_frame) & (df['Frame']<=switch_frame)].groupby(['Trap','Strain'],as_index=False)['GrowthRate'].mean()

    #discard outliers
    #means = growth_rate_before_switch_df.groupby('Strain')['GrowthRate'].mean()
    #stdevs = growth_rate_before_switch_df.groupby('Strain')['GrowthRate'].std()
    #cutoffs = means-3*stdevs
    if make_plot:
        growth_rate_before_switch_df.hist(by="Strain",column="GrowthRate",sharex=True,sharey=True)
    #print(cutoffs)
    for strain in growth_rate_before_switch_df['Strain'].unique():
        growth_rate_before_switch_df.loc[(growth_rate_before_switch_df['GrowthRate']<cutoff) & 
                                        (growth_rate_before_switch_df['Strain']==strain), 'GrowthRate'] = np.nan
        print(f"{strain}: discarded {growth_rate_before_switch_df.loc[growth_rate_before_switch_df['Strain']==strain,'GrowthRate'].isna().sum()} traps")


    #compute normalized growth rate
    growth_rate_before_switch_df=growth_rate_before_switch_df.set_index('Trap')
    #print(growth_rate_before_switch_df[10])
    df['Relative Growth Rate'] = df.apply(lambda x: x['GrowthRate'] / growth_rate_before_switch_df.loc[x['Trap'],'GrowthRate'], axis=1)
    df['GrowthRate']= df.apply(lambda x: np.nan if np.isnan(x['Relative Growth Rate']) else x['GrowthRate'],axis=1)
    return df


def add_antibiotic_info(ax,switch_timepoints,labels,color="black",lw=1):
    ylim = ax.get_ylim()
    fontweight="bold" if lw>1 else "normal"
    for t,l in zip(switch_timepoints,labels):
        ax.vlines(t,ylim[0],ylim[1],colors='black',linestyles="dashed",color=color,lw=lw)
        ax.text(t,ylim[0]+0.9*(ylim[1]-ylim[0])," + "+l,color=color,fontweight=fontweight)
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
