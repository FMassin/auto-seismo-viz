#!/usr/bin/env python

from matplotlib.pyplot import figure,colorbar
from matplotlib.gridspec import  GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator, AutoLocator, FormatStrFormatter
from matplotlib import patheffects , colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates


from numpy import meshgrid,logspace,log10,abs,median, percentile, nan
from numpy.ma import masked_array

from obspy.imaging.cm import pqlx
from obspy.signal.tf_misfit import cwt

from addons.core import eewmap
from obspy.core.event.catalog import Catalog


golden = (1 + 5 **0.5) / 2

def myaxstyle(ax):
    ax.set_facecolor('w')
    ax.tick_params(right=True, top=True,
                left=True, bottom=True,
                which='both')

    ax.set_yticklabels([ ], path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])  
    ax.set_xticklabels([ ], path_effects=[patheffects.withStroke(linewidth=4, foreground="w")]) 
    ax.set_ylabel('', path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])  
    ax.set_xlabel('', path_effects=[patheffects.withStroke(linewidth=4, foreground="w")]) 
 
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_major_locator(AutoLocator())

    ax.grid(visible=True, which='major', color='gray', linestyle='dashdot', zorder=-999)
    ax.grid(visible=True, which='minor', color='beige',  ls='-', zorder=-9999) 

def dataaxes(fig, nsubplots, gs, labelright, labelleft):
    gs01 = GridSpecFromSubplotSpec(nsubplots, 1, 
                                   height_ratios=[1 for n in range(nsubplots)], 
                                   hspace=0, 
                                   subplot_spec=gs)
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

    for ax in axes[1:]:
        ax.sharex(axes[0])
        ax.sharey(axes[0])

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
    


def scalogram(tr,ax):
    npts = tr.stats.npts
    dt = tr.stats.delta
    t = tr.times("matplotlib") #linspace(0, dt * npts, npts)
    f_min = max([0.5, 1/(npts*dt/8)])
    f_max = min([50, 1/dt/4])
    
    s = abs(cwt(tr.data-median(tr.data), dt, 8, f_min, f_max))

    x, y = meshgrid(t,
                    logspace(log10(f_min), 
                             log10(f_max), 
                             s.shape[0]))
    
    im = ax.pcolormesh(x, y, s, 
                       cmap=pqlx, 
                       norm=colors.Normalize(vmin=s.min()+(s.max()-s.min())/99, 
                                             vmax=s.max()),
                       alpha=0.8, 
                       zorder=999)

    ax.set_ylabel("Frequency [Hz]",position='right')
    ax.set_yscale('log')
    #ax.set_ylim(f_min, f_max)

    cbaxes = inset_axes(ax, width="1%", height="90%", loc=2) 
    cbaxes.set_yticklabels([ ], 
                           path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])  # vertically oriented colorbar
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

def allstationsdataaxes(output,nsubplots):

    fig = figure(figsize=(10,10*golden))
    
    gs0 = GridSpec(1, 2,   
                width_ratios=(golden, 1), 
                wspace=0, 
                figure=fig)
    gs00 = GridSpecFromSubplotSpec(2, 1,    
                                height_ratios=(1, golden), 
                                hspace=0, 
                                subplot_spec=gs0[1])

    ax_map = fig.add_subplot(gs00[0])
    ax_map.tick_params(direction="in", 
                       which='both')
    myaxstyle(ax_map)

    fig.axes_spec = dataaxes(fig, nsubplots, gs0[0], labelright=False, labelleft=True)
    fig.axes_data = dataaxes(fig, nsubplots, gs00[1], labelright=True, labelleft=False)
    fig.basename = '%s'%(output)
    return fig

def plotallstationsdata(streams,event,inventory, 
                        outputs=['vel'],#'raw',acc'], 
                        units={'acc':'m.s$^{-2}$',
                                'vel':'m/s',
                                'disp':'m',
                                'raw':'counts'},
                        max_num_station=16):
    
    print(streams['raw'].__str__(extended=True))

    stations = []
    print(event.preferred_origin())
    for i,trace in enumerate(streams['raw']):
        station = '%s.%s'%(trace.stats.network,trace.stats.station)
        ok=False
        for a in event.preferred_origin().arrivals:
            p=a.pick_id.get_referred_object()
            if station == p.waveform_id.network_code+'.'+p.waveform_id.station_code:
                ok=True
                break
        if not ok:
            continue
        if len(stations)>=max_num_station:
            break
        if station in stations:
            continue
        stations += [station]
    print(stations)
    figs = [allstationsdataaxes(output,len(stations)) for output in outputs]
    done = []
    for station in stations:

        ax_index = stations.index(station)

        selectopt = {'network':station.split('.')[0],
                     'station':station.split('.')[1],
                     'channel':'*Z'}
        stationwf = {output:streams[output].select(**selectopt) for output in streams if output in outputs}
        print(station)
        print(stationwf)
        ncha = len(stationwf[outputs[0]]) #  ALLOW MULTIPLE INSTRUMENTS SAME STATION? 

        for iout,(output, wf) in enumerate(stationwf.items()):    
            if iout+1>len(figs):
                break          
            for itr,tr in enumerate(wf): #  ALLOW MULTIPLE INSTRUMENTS SAME STATION? 
                a = figs[iout].axes_data[ax_index]
                a.plot(tr.times("matplotlib"), tr.data)
                a.set_xlim([min(tr.times("matplotlib")),max(tr.times("matplotlib"))])
                a.set_ylabel('%s'%( units[output]))
                a.legend([],[],
                        title='%s. %s'%(output.capitalize(),tr.id),
                        title_fontsize='x-small',
                        labelspacing=-.2)
                a.yaxis.set_label_position("right")
                a.tick_params(labelright=True, 
                                labelleft=False, 
                                labeltop=(ax_index==0), 
                                labelbottom=(ax_index+1)==len(figs[0].axes_spec),
                                which='both')
                
                locator = mdates.AutoDateLocator()
                formatter = mdates.ConciseDateFormatter(locator)
                a.xaxis.set_major_locator(locator)
                a.xaxis.set_major_formatter(formatter)

                a.yaxis.set_minor_locator(AutoMinorLocator())
                a.yaxis.set_major_locator(AutoLocator())
                a.yaxis.set_major_formatter(FormatStrFormatter("$%g$"))

                a=figs[iout].axes_spec[ax_index]                
                if True:
                    lax = scalogram(tr,a)      
                else:
                    a.plot(tr.times("matplotlib"), tr.data)          
                a.legend([],[],
                        title='%s. %s'%(output.capitalize(),tr.id),
                        title_fontsize='x-small',
                        labelspacing=-.2)
                a.yaxis.set_label_position("left")
                a.tick_params(labelright=False, 
                                labelleft=True, 
                                labeltop=(ax_index==0), 
                                labelbottom=(ax_index+1)==len(figs[0].axes_spec),
                                which='both')
                
                locator = mdates.AutoDateLocator()
                formatter = mdates.ConciseDateFormatter(locator)
                a.xaxis.set_major_locator(locator)
                a.xaxis.set_major_formatter(formatter)

                break
    
    return figs