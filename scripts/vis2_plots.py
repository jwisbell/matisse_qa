"""
Author: Jacob Isbell 
Date: October 1, 2024

Script which contains functions to plot visibility and correlated flux data.
"""
import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import baseline_idx_from_stapair

mpl.rcParams['font.family']='serif'
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['xtick.top']=True
mpl.rcParams['ytick.direction']='in'
mpl.rcParams['ytick.right']=True


normpa = mpl.colors.Normalize(vmin=-180,vmax=180)
normvis = mpl.colors.Normalize(vmin=0,vmax=1.1)


def plot_vis(data_dict, output_dir:str = '/.', save_fig:bool = False):
    """
    data_dict: dictionary containing the following keys
    'vis':{"cflux":[],"cflux_err":[],"u":[], "v":[],"wl_vis":[],"vis2":[], "vis2_err":[],'bcd':[],'vis2_sta':[]

    Used to plot the correlated flux and visibility
    """

    fig1, axarr1 = plt.subplots(3,3, figsize=(8.5,8.5))
    fig2, axarr2 = plt.subplots(3,3, figsize=(8.5,8.5))
    fig3, axarr3 = plt.subplots(4,3, figsize=(8.5,8.5))

    bl_lengths = [[],[],[],[],[],[]]
    bl_names = ['UT3-UT4', 'UT1-UT2','UT2-UT3','UT2-UT4','UT1-UT3','UT1-UT4' ]
    for i in range(len(data_dict['cflux'])):
        u = data_dict['u'][i]
        v = data_dict['v'][i]
        bl = np.sqrt(u**2 + v**2)
        sta = data_dict['vis2_sta'][i]
        idx = baseline_idx_from_stapair(sta)
        wls = data_dict['wl_vis'][i]
        #name = baseline_name_from_stapair(sta)
        
        bl_lengths[idx].append(bl)
        pa = np.degrees(np.arctan2(v,u))
        im = axarr1.flatten()[idx].scatter(  data_dict['wl_vis'][i]*1e6, data_dict['cflux'][i], s=2, c=[pa]*len(data_dict['wl_vis'][i]), norm=normpa, cmap='twilight',zorder=1)
        axarr1.flatten()[idx].errorbar(  data_dict['wl_vis'][i]*1e6, data_dict['cflux'][i], yerr=data_dict['cflux_err'][i],zorder=0,color='k',ls='none',marker='.',alpha=0.1)
        #axarr1.flatten()[i%6].set_ylim([0,0.5])
        axarr1.flatten()[idx].set_ylim([0,1.1])

        axarr1.flatten()[6].scatter(bl / wls / 1e6, data_dict['cflux'][i],c=[pa]*len(data_dict['wl_vis'][i]), norm=normpa, cmap='twilight',zorder=1, s=2)
        axarr1.flatten()[6].errorbar(bl / wls / 1e6, data_dict['cflux'][i], yerr=data_dict['cflux_err'][i],zorder=0,color='k',ls='none',marker='.',alpha=0.1)
        axarr1.flatten()[6].set_ylim([0,1.1])

        im_uv1 = axarr1.flatten()[7].scatter(u/wls[::4]/1e6, v/wls[::4]/1e6, c=data_dict['cflux'][i][::4], norm=normvis, cmap='rainbow',zorder=1,alpha=0.5,s=2)
        axarr1.flatten()[7].scatter(-u/wls[::4]/1e6, -v/wls[::4]/1e6, c=data_dict['cflux'][i][::4], norm=normvis, cmap='rainbow',zorder=1,alpha=0.5,s=2)
        axarr1.flatten()[7].set_xlim([45,-45])
        axarr1.flatten()[7].set_ylim([-45,45])
        axarr1.flatten()[7].set_xlabel('u [Mlambda]')
        axarr1.flatten()[7].set_ylabel('v [Mlambda]')
        axarr1.flatten()[7].set_aspect('equal')
        

    bcd_sorted = [[[],[],[],[],[],[]], [[],[],[],[],[],[]], [[],[],[],[],[],[]], [[],[],[],[],[],[]],[[],[],[],[],[],[]],[[],[],[],[],[],[]]]
    for i in range(len(data_dict['vis2'])):
        u = data_dict['u'][i]
        v = data_dict['v'][i]
        bl = np.sqrt(u**2 + v**2)
        sta = data_dict['vis2_sta'][i]
        idx = baseline_idx_from_stapair(sta)

        bl_lengths[idx].append(bl)
        pa = np.degrees(np.arctan2(v,u))

        im2 = axarr2.flatten()[idx].scatter(  data_dict['wl_vis'][i]*1e6, data_dict['vis2'][i], s=2, c=[pa]*len(data_dict['wl_vis'][i]), norm=normpa, cmap='twilight',zorder=1)
        axarr2.flatten()[idx].errorbar(  data_dict['wl_vis'][i]*1e6, data_dict['vis2'][i], yerr=data_dict['vis2_err'][i],zorder=0,color='k',ls='none',marker='.',alpha=0.1)
        axarr2.flatten()[idx].set_ylim([0,1.1])

        axarr2.flatten()[6].scatter(bl / wls / 1e6, data_dict['vis2'][i],c=[pa]*len(data_dict['wl_vis'][i]), norm=normpa, cmap='twilight',zorder=1,s=2)
        axarr2.flatten()[6].errorbar(bl / wls / 1e6, data_dict['vis2'][i], yerr=data_dict['vis2_err'][i],zorder=0,color='k',ls='none',marker='.',alpha=0.1)
        axarr2.flatten()[6].set_ylim([0,1.1])

        im2_uv1 = axarr2.flatten()[7].scatter(u/wls[::4]/1e6, v/wls[::4]/1e6, c=data_dict['vis2'][i][::4], norm=normvis, cmap='rainbow',zorder=1,alpha=0.5,s=2)
        axarr2.flatten()[7].scatter(-u/wls[::4]/1e6, -v/wls[::4]/1e6, c=data_dict['vis2'][i][::4], norm=normvis, cmap='rainbow',zorder=1,alpha=0.5,s=2)
        axarr2.flatten()[7].set_xlim([45,-45])
        axarr2.flatten()[7].set_ylim([-45,45])
        axarr2.flatten()[7].set_xlabel('u [Mlambda]')
        axarr2.flatten()[7].set_ylabel('v [Mlambda]')
        axarr2.flatten()[7].set_aspect('equal')

        bcd = data_dict['bcd'][i]
        if bcd == 'oo':
            bcd_sorted[0][idx].append( data_dict['vis2'][i] )
        elif bcd == 'ii':
            bcd_sorted[1][idx].append( data_dict['vis2'][i] )
        elif bcd == 'oi':
            bcd_sorted[2][idx].append( data_dict['vis2'][i] )
        elif bcd == 'io':
            bcd_sorted[3][idx].append( data_dict['vis2'][i] )
        elif bcd == 'oo_phot':
            bcd_sorted[4][idx].append( data_dict['vis2'][i] )
        elif bcd == 'ii_phot':
            bcd_sorted[5][idx].append( data_dict['vis2'][i] )

        
    for i in range(6):
        label = 'in-in'
        if i != 0:
            label = None 
        axarr3.flatten()[i].errorbar( np.median(bcd_sorted[0][i],0), np.median(bcd_sorted[1][i],0), yerr=np.std(bcd_sorted[1][i],0), c='r',zorder=1, marker='o',ls='none',alpha=0.5,label=label)
        label = 'out-in'
        if i != 0:
            label = None 
        axarr3.flatten()[i].errorbar( np.median(bcd_sorted[0][i],0), np.median(bcd_sorted[2][i],0), yerr=np.std(bcd_sorted[2][i],0),c='b',zorder=1, marker='s',ls='none',alpha=0.5,label=label)
        label = 'in-out'
        if i != 0:
            label = None 
        axarr3.flatten()[i].errorbar( np.median(bcd_sorted[0][i],0), np.median(bcd_sorted[3][i],0), yerr=np.std(bcd_sorted[3][i],0), c='g',zorder=1, marker='^',ls='none',alpha=0.5,label=label)
        

        axarr3.flatten()[i+6].errorbar( np.median(bcd_sorted[4][i],0), np.median(bcd_sorted[5][i],0), yerr=np.std(bcd_sorted[4][i],0), c='r',zorder=1, marker='o',ls='none',alpha=0.5)
        
        
        axarr3.flatten()[i].plot([0,1.1],[0,1.1],c='k',ls='--',zorder=0)
        axarr3.flatten()[i+6].plot([0,1.1],[0,1.1],c='k',ls='--',zorder=0)
        axarr3.flatten()[i].set_ylim([0,1.1])
        axarr3.flatten()[i].set_xlim([0,1.1])
        axarr3.flatten()[i].set_title('UNCHOPPED - '+bl_names[i])
        
        axarr3.flatten()[i+6].set_ylim([0,1.1])
        axarr3.flatten()[i+6].set_xlim([0,1.1])
        axarr3.flatten()[i+6].set_title('CHOPPED - '+bl_names[i])
    axarr3.flatten()[0].legend()
    axarr3[3,0].set_xlabel('Out-Out Vis2')
    axarr3[3,0].set_ylabel('Other BCD Vis2')
    


    axarr1.flatten()[-1].axis('off')
    #axarr1.flatten()[-2].axis('off')
    #axarr1.flatten()[-3].axis('off')

    axarr2.flatten()[-1].axis('off')
    #axarr2.flatten()[-2].axis('off')
    #axarr2.flatten()[-3].axis('off')

    for i in range(6):
        axarr1.flatten()[i].set_title(f"{bl_names[i]}\nMean Proj. BL: {np.mean(bl_lengths[i]):.1f} m")
        axarr2.flatten()[i].set_title(f"{bl_names[i]}\nMean Proj. BL: {np.mean(bl_lengths[i]):.1f} m")

    #axarr1[1,0].set_ylabel("Flux [cts]")
    axarr1[2,0].set_ylabel("Coherent Vis [cts or None]")
    axarr1[1,0].set_xlabel("Wavelength [micron]")
    axarr1[2,0].set_xlabel("Spatial Frequency [MLambda]")
    axarr2[2,0].set_ylabel("Incoherent Vis2")
    axarr2[1,0].set_xlabel("Wavelength [micron]")
    axarr1[2,0].set_xlabel("Spatial Frequency [MLambda]")

    #axarr2.suptitle.set_text(f'')

    plt.colorbar(im, ax=axarr1.flatten()[-1], orientation='horizontal', label="Position Angle [deg]")
    plt.colorbar(im_uv1, ax=axarr1.flatten()[-2], orientation='vertical', label="Vis")
    plt.colorbar(im2, ax=axarr2.flatten()[-1], orientation='horizontal', label="Position Angle [deg]")
    plt.colorbar(im2_uv1, ax=axarr2.flatten()[-2], orientation='vertical', label="Vis")

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()

    if output_dir is not None and save_fig:
        fig1.savefig(f"{output_dir}/coher_vis.png")
        fig2.savefig(f"{output_dir}/incoher_vis.png")
        fig3.savefig(f"{output_dir}/bcd_compare_vis.png")

    
    if not save_fig:
        plt.show()
    plt.close("all")


    return None 