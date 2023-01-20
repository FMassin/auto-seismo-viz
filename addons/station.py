#!/usr/bin/env python
from matplotlib.pyplot import figure
from matplotlib.gridspec import  GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator
from numpy import meshgrid,logspace,log10,abs,median
from obspy.imaging.cm import pqlx
from obspy.signal.tf_misfit import cwt
from matplotlib.pyplot import colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import patheffects 
from addons.core import eewmap
from obspy.core.event.catalog import Catalog


golden = (1 + 5 **0.5) / 2

def myaxstyle(ax):
    ax.set_facecolor('w')
    ax.tick_params(right=True, top=True,
                left=True, bottom=True,
                which='both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(visible=True, which='major', color='gray', linestyle='dashdot', zorder=-9999)
    ax.grid(visible=True, which='minor', color='beige',  ls='-', zorder=-9999)  

def dataaxes(fig,nsubplots, gs, labelright, labelleft):
    gs01 = GridSpecFromSubplotSpec(nsubplots, 1, height_ratios=[1 for n in range(nsubplots)], hspace=0, subplot_spec=gs)
    axes = []
    for i in range(nsubplots):
        axes += [fig.add_subplot(gs01[i])]
        axes[-1].tick_params(labelright=labelright, 
                    labelleft=labelleft, 
                    labeltop=i==0, 
                    labelbottom=i+1==nsubplots, 
                    direction="in", 
                    which='both')
        myaxstyle(axes[-1]) 
    return axes

def stationdataaxes(stationwf,nsubplots):
    fig = figure(figsize=(10*golden,10))
    
    gs0 = GridSpec(1, 2, width_ratios=(1,golden), wspace=0, figure=fig)
    gs00 = GridSpecFromSubplotSpec(2, 1, height_ratios=(golden,1), hspace=0, subplot_spec=gs0[0])

    ax_map = fig.add_subplot(gs00[1])
    ax_map.tick_params(direction="in", 
                    which='both')
    myaxstyle(ax_map)

    axes_spec = dataaxes(fig, nsubplots, gs0[1], labelright=True, labelleft=False)
    axes_data = dataaxes(fig, nsubplots, gs00[0], labelright=False, labelleft=True)
    axes_data[-1].tick_params(labelbottom=False)
    
    return ax_map, axes_data, axes_spec
    


def scalogram(tr,ax,output):
    npts = tr.stats.npts
    dt = tr.stats.delta
    t = tr.times("matplotlib") #linspace(0, dt * npts, npts)
    f_min = 1/(npts*dt/8)
    f_max = 1/dt/4
    
    s = cwt(tr.data-median(tr.data), dt, 8, f_min, f_max)

    x, y = meshgrid(t,
                    logspace(log10(f_min), 
                             log10(f_max), 
                             s.shape[0]))
    
    im = ax.pcolormesh(x, y, abs(s), cmap=pqlx)

    ax.set_ylabel("Frequency [Hz]",position='right')
    ax.set_yscale('log')
    ax.set_ylim(f_min, f_max)

    cbaxes = inset_axes(ax, width="1%", height="90%", loc=2) 
    cbaxes.set_yticklabels([ ], path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])  # vertically oriented colorbar
    colorbar(im, cax=cbaxes, orientation='vertical')
    cbaxes.tick_params(labelright=True, 
                labelleft=False, 
                labeltop=False, 
                labelbottom=False,
                which='both')
    cbaxes.yaxis.set_minor_locator(AutoMinorLocator())

    lax = ax.figure.add_axes(ax.get_position())
    lax.patch.set_facecolor('None')
    lax.tick_params(right=False, top=False,
                    left=False, bottom=False,
                    labelright=False, labeltop=False,
                    labelleft=False, labelbottom=False,
                    which='both')
    for pos in ['right', 'top', 'bottom', 'left']:
        lax.spines[pos].set_visible(False)
    return lax


def plotstationdata(streams,event,inventory, 
                    outputs=['raw'],#acc','vel'], 
                    units={'acc':'m.s$^{-2}$',
                           'vel':'m/s',
                           'disp':'m',
                           'raw':'counts'}):
    stations=[]
    figs=[]
    for i,trace in enumerate(streams['raw']):
        if len(stations)>6:
            break
        station = '%s.%s'%(trace.stats.network,trace.stats.station)
        if station in stations:
            continue
        stations += [station]
        stationwf = {output:streams[output].select(station=trace.stats.station,channel='*[Z,E,N,1,2,3]') for output in streams if output in outputs}
        ncha = len(stationwf[outputs[0]])
        noutput = len(stationwf)

        ax_map, axes_data, axes_spec = stationdataaxes(stationwf, ncha*noutput)

        for iout,(output, wf) in enumerate(stationwf.items()):                 
            for itr,tr in enumerate(wf):
                a = axes_data[itr+iout*ncha]
                a.plot(tr.times("matplotlib"), tr.data)
                a.xaxis_date()
                a.set_xlim([min(tr.times("matplotlib")),max(tr.times("matplotlib"))])
                a.set_ylabel('%s'%( units[output]))
                a.legend([],[],title='%s. %s'%(output.capitalize(),tr.id),labelspacing=-.2)

        for iout,(output, wf) in enumerate(stationwf.items()):                 
            for itr,tr in enumerate(wf):
                a=axes_spec[itr+iout*ncha]
                lax = scalogram(tr,a,output)
                lax.legend([],[],title='%s. %s'%(output.capitalize(),tr.id),labelspacing=-.2)
                a.xaxis_date()
                a.yaxis.set_label_position("right")
        
        eewmap({'event':event,
                'inventory':inventory.select(network=trace.stats.network,
                                             station=trace.stats.station)},
               title=False,
               legend=False,
               mapbounds=False,
               others=[Catalog(events=[event])],
               ax=ax_map)
        #ax_map.legend([],[],title=trace.id,labelspacing=-.2)
        a.figure.basename = '%s.%s'%(station,'.'.join(outputs))
        figs += [a.figure]
    
    return figs