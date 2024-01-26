#!/usr/bin/env python

from matplotlib.pyplot import figure,colorbar,rcParams
from matplotlib.gridspec import  GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator, AutoLocator, EngFormatter, LogLocator
from matplotlib import patheffects , colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates
from matplotlib.markers import MarkerStyle, CapStyle


from numpy import meshgrid,logspace,log10,abs,median, percentile, nan,linspace,argmin,argsort
from numpy.ma import masked_array

from scipy.interpolate import interp2d

from obspy.imaging.cm import pqlx
from obspy.signal.tf_misfit import cwt
from obspy.geodetics.base import gps2dist_azimuth

from fnmatch import fnmatch
from cartopy.crs import  Geodetic


from addons.core import eewmap
from obspy.core.event.catalog import Catalog
from addons.bmap import bmap


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

def dataaxes(fig, nsubplots, gs, labelright, labelleft,sharey=False):
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
        if sharey:
            ax.sharey(axes[0])

    return axes

def stationdataaxes(stationwf,nsubplots,map):
    fig = figure(figsize=(10*golden,10))
    
    gs0 = GridSpec(1, 2, width_ratios=(1,golden), wspace=0, figure=fig)

    if map:
        gs00 = GridSpecFromSubplotSpec(2, 1, height_ratios=(golden,1), hspace=0, subplot_spec=gs0[0])

        ax_map = fig.add_subplot(gs00[1])
        ax_map.tick_params(direction="in", 
                        which='both')
        myaxstyle(ax_map)
        axes_data = dataaxes(fig, nsubplots, gs00[0], labelright=False, labelleft=True)

    else :
        ax_map = None    
        axes_data = dataaxes(fig, nsubplots, gs0[0], labelright=False, labelleft=True)


    axes_spec = dataaxes(fig, nsubplots, gs0[1], labelright=True, labelleft=False, sharey=True)
    axes_data[-1].tick_params(labelbottom=False)
    
    return ax_map, axes_data, axes_spec
    


def scalogram(trace,
              ax,
              output,
              nfreqs=64,
              ntimes=512*2
              ):
    
    tr = trace.copy()
    
    npts = tr.stats.npts
    dt = tr.stats.delta
    t = tr.times("matplotlib") #linspace(0, dt * npts, npts)
    f_min = max([0.5, 1/(npts*dt/8)])
    f_max = min([50, 1/dt/4])
    
    #tr.detrend('polynomial',order=4)
    tr.detrend()
    tr.filter('highpass',freq=f_min)
    
    
    s = abs(cwt(tr.data-median(tr.data), dt, 8, f_min, f_max, nfreqs))

    x, y = meshgrid(t,
                    logspace(log10(f_min), 
                             log10(f_max), 
                             s.shape[0]))

    # scipy interp. cubic
    f = interp2d(x[0,:], 
                 y[:,0], 
                 s, 
                 kind='cubic'
                 )
    xnew = linspace(min(x[0,:]),max(x[0,:]),ntimes)
    ynew = y[:,0]
    Sn = f(xnew,ynew)
    Xn, Yn = meshgrid(xnew, ynew)

    nanthreshold = Sn.min()+(Sn.max()-Sn.min())/128
    if not hasattr(ax.figure,'scalogram_lims'):
        ax.figure.scalogram_lims = [nanthreshold,Sn.max()]
    else:
        ax.figure.scalogram_lims[0] = min([nanthreshold,ax.figure.scalogram_lims[0]])
        ax.figure.scalogram_lims[1] = max([Sn.max(),ax.figure.scalogram_lims[1]])
    #print(ax.figure.scalogram_lims, (Sn.min(),Sn.max()))

    
    Sn[Sn<nanthreshold] = nan

    ax.scalogram = ax.pcolormesh(Xn, Yn, Sn, #x, y, s, #
                       cmap=pqlx,
                       norm=colors.LogNorm(*ax.figure.scalogram_lims),
                       #alpha=0.7,
                       rasterized=True, 
                       zorder=999)

    ax.set_ylabel("Frequency (Hz)",position='right')
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(EngFormatter(places=1, sep="\N{THIN SPACE}"))

    cbaxes = inset_axes(ax, 
                        width="1%", 
                        height="90%", 
                        loc='center left') 
    cbaxes.set_yticklabels([ ],
                           fontsize='xx-small', 
                           path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])  # vertically oriented colorbar
    colorbar(ax.scalogram, 
             cax=cbaxes, 
             orientation='vertical')
    cbaxes.tick_params(labelright=True, 
                labelleft=False, 
                labeltop=False, 
                labelbottom=False,
                which='both')
    cbaxes.yaxis.set_minor_locator(LogLocator())
    cbaxes.yaxis.set_major_formatter(EngFormatter(unit=output, 
                                                  places=1, 
                                                  sep="\N{THIN SPACE}"))
    ax.scalogram_cbar = cbaxes
    if False:
        lax = ax.figure.add_axes(ax.get_position())
        lax.patch.set_facecolor('None')
        lax.tick_params(right=False, top=False,
                        left=False, bottom=False,
                        labelright=False, labeltop=False,
                        labelleft=False, labelbottom=False,
                        which='both')
        for pos in ['right', 'top', 'bottom', 'left']:
            lax.spines[pos].set_visible(False)
    return ax


def plotstationdata(streams,event,inventory, 
                    outputs=['vel'],#raw'],#acc','vel'], 
                    units={'acc':'m.s$^{-2}$',
                           'vel':'m/s',
                           'disp':'m',
                           'raw':'counts'},
                    max_num_station=4,
                    map=False):
    stations=[]
    figs=[]
    for i,trace in enumerate(streams['raw']):
        if len(stations)>=max_num_station:
            break
        station = '%s.%s'%(trace.stats.network,trace.stats.station)
        if station in stations:
            continue
        stations += [station]
        stationwf = {output:streams[output].select(station=trace.stats.station,channel='*[Z,E,N,1,2,3]') for output in streams if output in outputs}
        ncha = len(stationwf[outputs[0]])
        noutput = len(stationwf)

        ax_map, axes_data, axes_spec = stationdataaxes(stationwf, ncha*noutput, map)

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
                lax = scalogram(tr,a,units[output])
                leg = lax.legend([],[],
                           title='%s. %s'%(output.capitalize(),tr.id),
                           labelspacing=-.2,
                           )
                leg.set_zorder(99999)
                a.xaxis_date()
                a.yaxis.set_label_position("right")
        if map:
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

def allstationsdataaxes(output,nsubplots,map=False):

    fig = figure(figsize=(10,10*golden))
    fig.basename = '%s'%(output)
    
    gs0 = GridSpec(1, 2,   
                width_ratios=(golden, 1), 
                wspace=0, 
                figure=fig)
    if map:
        gs00 = GridSpecFromSubplotSpec(2, 1,    
                                    height_ratios=(1, golden), 
                                    hspace=0, 
                                    subplot_spec=gs0[1])
        fig.ax_map = fig.add_subplot(gs00[0])
        fig.ax_map.tick_params(direction="in", 
                        which='both')
        myaxstyle(fig.ax_map)
        fig.axes_data = dataaxes(fig, nsubplots, gs00[1], labelright=True, labelleft=False)
    
    else:
        fig.axes_data = dataaxes(fig, nsubplots, gs0[1], labelright=True, labelleft=False)

    fig.axes_spec = dataaxes(fig, nsubplots, gs0[0], labelright=False, labelleft=True, sharey=True)
    return fig

def plottimes(a,
              vlinelabels,
              vlinetimes,
              alllines=False,
              y=0
              ):

    for i,label in enumerate(vlinelabels):
        for t in vlinetimes[i]:
            if 'rg' in vlinelabels[i] or alllines:
                a.axvline(t,
                        label=label,
                        color='C%d'%range(1,10)[vlinelabels.index(vlinelabels[i])],
                        zorder=999999999,
                        linewidth=1,
                        path_effects=[patheffects.withStroke(linewidth=2, foreground="w")]
                        )
            else:
                a.plot(t,y,
                        marker=4,
                        label=label,
                        markeredgecolor='w',
                        markeredgewidth=1,
                        color='C%d'%range(1,10)[vlinelabels.index(vlinelabels[i])]
                        )
            label=None
            
def plotallstationsdata(streams,event,inventory, 
                        outputs=['vel'],#'raw',acc'], 
                        units={'acc':'m.s$^{-2}$',
                                'vel':'m/s',
                                'disp':'m',
                                'raw':'counts'},
                        max_num_station=4,#16
                        sharec=True,
                        ot=False,
                        ew=True,#False,
                        map=False,
                        all=True,
                        orderbydistancefrom=None,#(14.73878,-91.56802),#Santiaguuito volcano
                        whitelist='*',#GI.STG*,CH.*,8D.*',
                        blacklist='*.BH?,*.HD?,*.BD?,*GI.STG2.*,GI.STG8.*,GI.STG5.*.HH?,GI.STG11.*,GI.STG*SH?'
                        ):
    
    # EMPTY SPECTROGRAM
    # fix gain EHZ (399635995.8 c/m/s)

    blacklist = blacklist.split(',')
    whitelist = whitelist.split(',')
    #print(streams['raw'].__str__(extended=True))
    #print(event.preferred_origin())

    instruments = []
    instrumentcoordinates = []
    if all:
        for n in inventory:
            if len(instruments)==max_num_station:
                break
            for s in n:
                if len(instruments)==max_num_station:
                    break
                for c in s:
                    if len(instruments)==max_num_station:
                        break
                    mseedid = '%s.%s.%s.%s'%(n.code,s.code,c.location_code,c.code)
                    if mseedid[:-1] in instruments:
                        continue

                    whitelisted = False
                    for white in whitelist:
                        if fnmatch(mseedid, white):
                            whitelisted = True
                    if not whitelisted:
                        print(mseedid,'is not whitelisted')
                        continue

                    blacklisted = False
                    for black in blacklist :
                        if fnmatch(mseedid, black):
                            blacklisted = True
                    if blacklisted:
                        print(mseedid,'is blacklisted')
                        continue

                    instruments += [mseedid[:-1]]
                    instrumentcoordinates += [(c.latitude,c.longitude)]
    else:
        for output in streams:
            if output not in outputs:
                continue
            for i,trace in enumerate(streams[output]):

                mseedid = '%s.%s.%s.%s'%(trace.stats.network,trace.stats.station,trace.stats.location,trace.stats.channel)
                if True in [fnmatch(mseedid, black) for black in blacklist ]:
                    print(mseedid,'is blacklisted')
                    continue

                if True not in [fnmatch(mseedid, white) for white in whitelist ]:
                    print(mseedid,'is not whitelisted')
                    continue

                print(trace)
                instrument = '%s.%s.%s.%s'%(trace.stats.network,
                                    trace.stats.station,
                                    trace.stats.location,
                                    trace.stats.channel[:-1])
                ok=False
                for o in event.origins:
                    for a in o.arrivals:
                        p=a.pick_id.get_referred_object()
                        if '.'.join(instrument.split('.')[:-2]) == '%s.%s'%(p.waveform_id.network_code,p.waveform_id.station_code):
                            ok=True
                            break
                if not ok:
                    continue
                if len(instruments)>=max_num_station:
                    break
                if instrument in instruments:
                    continue
                instruments += [instrument]

    if orderbydistancefrom is not None:
        distances = [gps2dist_azimuth(*orderbydistancefrom,*c)[0] for c in instrumentcoordinates]
        instruments = [instruments[i] for i in argsort(distances)]

        
        
    #print(instruments)
    figs = [allstationsdataaxes(output,len(instruments),map=map) for output in outputs]
    
    if map:
        for fig in figs:
            la = [s.latitude for n in inventory for s in n]
            lo = [s.longitude for n in inventory for s in n]
            typ = ['HH']
            fig.bmap = bmap(bounds=[min(lo), max(lo), min(la), max(la)],
                right=True,
                top=True,
                ax=fig.ax_map,
                zoom=15)
            marker = MarkerStyle('2')
            marker._capstyle = CapStyle.round
            fig.bmap.plot(lo,la,
                        linestyle='',
                        marker = marker,
                        markeredgewidth = 2,
                        path_effects = [patheffects.withStroke(linewidth=4,foreground="w")],
                        transform=Geodetic())

    #return figs
    # Preferred origin
    preforig = event.preferred_origin()
    if preforig is None:
        preforig = event.origins[-1]
    for o in event.origins:
        pref = ['','pref.'][str(event.preferred_origin().resource_id) == str(o.resource_id)]
    
    # Earliest origin
    first = argmin([o.creation_info.creation_time for o in event.origins])
    first = event.origins[first]

    vlinelabels = []
    vlinetimes = []
    for instrument in instruments:

        ax_index = instruments.index(instrument)
        wfid = instrument.split('.')
        selectopt = {'network':wfid[0],
                     'station':wfid[1],
                     'location':wfid[2],
                     'channel':'%sZ' % wfid[3]
                     }
        #print(selectopt)
        stationwf = {output:streams[output].select(**selectopt) for output in streams if output in outputs}
        print(instrument)
        #print(stationwf[outputs[0]])
        ncha = len(stationwf[outputs[0]]) #  ALLOW MULTIPLE INSTRUMENTS SAME STATION? 

        for iout,(output, wf) in enumerate(stationwf.items()):    
            if iout+1>len(figs):
                break
            
            for i,label in enumerate(vlinelabels):
                vlinetimes[i] = []                   

            if ot and len(pref):
                label = '$Org._{T.}^{pref.}$'
                vlinelabels += [label] 
                vlinetimes += [[pref.time.matplotlib_date]] 
            
            if ew and len(first):
                label = '$Org._{CT.}^{early.}$'
                vlinelabels += [label] 
                vlinetimes += [[first.creation_info.creation_time.matplotlib_date]] 

            for o in event.origins:

                for a in o.arrivals:
                    p=a.pick_id.get_referred_object()
                    if '.'.join(instrument.split('.')[:-2]) == '%s.%s'%(p.waveform_id.network_code,p.waveform_id.station_code):
                        label = '$Pick_{%s}^{%s.%s}$' % ( p.phase_hint, p.evaluation_mode[:4], pref ) 
                        if label not in vlinelabels:
                            vlinelabels += [label] 
                            vlinetimes += [[]] 
                        vlinetimes[vlinelabels.index(label)] += [p.time.matplotlib_date]

            for itr,tr in enumerate(wf): #  ALLOW MULTIPLE INSTRUMENTS SAME STATION? 
                
                a = figs[iout].axes_data[ax_index]

                a.plot(tr.times("matplotlib"), tr.data)

                plottimes(a,vlinelabels,vlinetimes,y=median(tr.data))

                if ax_index == 0:
                    leg = a.legend(bbox_to_anchor=(0.98, 1.05), 
                                   loc='lower right',
                                   ncol=3, 
                                   fontsize='x-small',
                                   title='%s on %s'%(event.event_type.capitalize(), str(preforig.time)[:16]),
                                   title_fontsize='x-small',
                                   #mode="expand", 
                                   borderaxespad=0.)
                    a.add_artist(leg)

                a.set_xlim([min(tr.times("matplotlib")),max(tr.times("matplotlib"))])
                
                a.legend([],[],
                        title='%s. %s'%(output.capitalize(),tr.id),
                        title_fontsize='xx-small',
                        labelspacing=-.2)
                a.yaxis.set_label_position("right")
                if ax_index == int(len(figs[iout].axes_data)/2):
                    a.set_ylabel('%s. (%s)'%( output.capitalize(), units[output] ))
                a.tick_params(labelright=True, 
                                labelleft=False, 
                                labeltop=False, #( or not map), 
                                labelbottom=(ax_index+1)==len(figs[0].axes_spec),
                                which='both')
                a.spines['top'].set_visible(ax_index==0)
                a.spines['bottom'].set_visible(ax_index+1)==len(figs[0].axes_spec)
                a.tick_params(right=True, 
                              top=(ax_index==0),
                              left=True, 
                              bottom=(ax_index+1)==len(figs[0].axes_spec),
                              which='both')

                
                locator = mdates.AutoDateLocator()
                formatter = mdates.ConciseDateFormatter(locator)
                a.xaxis.set_major_locator(locator)
                a.xaxis.set_major_formatter(formatter)

                a.yaxis.set_minor_locator(AutoMinorLocator())
                a.yaxis.set_major_locator(AutoLocator())
                a.yaxis.set_major_formatter(EngFormatter(places=1, sep="\N{THIN SPACE}"))

                a=figs[iout].axes_spec[ax_index]          
                lax = scalogram(tr,a,units[output])

                plottimes(a,vlinelabels,vlinetimes,alllines=True)

                leg = a.legend([],[],
                        title='%s. %s'%(output.capitalize(),tr.id),
                        title_fontsize='x-small',
                        labelspacing=-.2,
                        )
                leg.set_zorder(99999)
                a.yaxis.set_label_position("left")
                if ax_index == int(len(figs[iout].axes_data)/2):
                    a.set_ylabel("Frequency (Hz)",position='right')
                else:
                    a.set_ylabel("")
                    if sharec:
                        a.scalogram_cbar.remove()

                a.tick_params(labelright=False, 
                                labelleft=True, 
                                labeltop=(ax_index==0), 
                                labelbottom=(ax_index+1)==len(figs[0].axes_spec),
                                which='both')
                
                a.spines['top'].set_visible(ax_index==0)
                a.spines['bottom'].set_visible(ax_index+1)==len(figs[0].axes_spec)
                a.tick_params(right=True, 
                              top=(ax_index==0),
                              left=True, 
                              bottom=(ax_index+1)==len(figs[0].axes_spec),
                              which='both')

                locator = mdates.AutoDateLocator()
                formatter = mdates.ConciseDateFormatter(locator)
                a.xaxis.set_major_locator(locator)
                a.xaxis.set_major_formatter(formatter)

                a.grid(visible=True, which='major', color='gray', linestyle='dashdot', zorder=-999)
                a.grid(visible=True, which='minor', color='beige',  ls='-', zorder=-9999) 

                break

    # Make all scalogram same color limits
    if sharec:
        for f,fig in enumerate(figs):
            for a in fig.axes_spec:
                try:
                    a.scalogram.set_clim(*fig.scalogram_lims)
                    a.scalogram_cbar.yaxis.set_major_formatter(EngFormatter(unit=units[outputs[f]], 
                                                                places=1, 
                                                                sep="\N{THIN SPACE}"))
                except:
                    pass
    
    return figs