#!/usr/bin/env python

from matplotlib.pyplot import figure,colorbar,rcParams
from matplotlib.gridspec import  GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator, AutoLocator, EngFormatter, LogLocator
from matplotlib import patheffects , colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.dates as mdates
from matplotlib.markers import MarkerStyle, CapStyle
from obspy.core.util import Enum


from numpy import meshgrid,cumsum,logspace,argmax,log10,abs,median, std,percentile, nan,linspace,argmin,argsort,nanmin,nanmax,sort,average
from numpy.ma import masked_array

from scipy.interpolate import interp2d,interp1d
from scipy.ndimage import center_of_mass

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
                    labeltop=0,
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

def moving_average(a, n=3):
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def scalogram(trace,
              ax,
              units,
              output,
              nfreqs=64,
              ntimes=512*2,
              cb=True,
              spec_pick=False,
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
        ax.figure.scalogram_lims[0] = min([nanthreshold,ax.figure.scalogram_lims[0]])*.99
        ax.figure.scalogram_lims[1] = max([Sn.max(),ax.figure.scalogram_lims[1]])
    #print(ax.figure.scalogram_lims, (Sn.min(),Sn.max()))

    if hasattr(ax.figure,'scalogram_cmax') and ax.figure.scalogram_cmax is not None:
        ax.figure.scalogram_lims[1] = ax.figure.scalogram_cmax
    if hasattr(ax.figure,'scalogram_cmin') and ax.figure.scalogram_cmin is not None:
        ax.figure.scalogram_lims[0] = ax.figure.scalogram_cmin*.99
        nanthreshold = ax.figure.scalogram_lims[0]

    if spec_pick:
        average_freqs = average(Yn,
                                weights=Sn,
                                axis=0)

        average_freq_amps = [f(Xn[0,iamp],amp)[0] for iamp,amp in enumerate(average_freqs)]

        lowpass_average_freq_amps = moving_average(average_freq_amps, n=5)
        maxindex = argmax(lowpass_average_freq_amps)

        ax.freq_max_amp = average_freq_amps[maxindex]
        ax.freq_max_time = Xn[0,maxindex]
        ax.freq_max = average_freqs[argmax(average_freq_amps)]
        print('Max freq: %s m/s at %s Hz on %s'%(ax.freq_max_amp,
                                                 ax.freq_max,
                                                 ax.freq_max_time ))


    Sn[Sn<=nanthreshold] = nan

    ax.scalogram = ax.pcolormesh(Xn, Yn, Sn, #x, y, s, #
                       cmap=pqlx,
                       norm=colors.LogNorm(*ax.figure.scalogram_lims),
                       #alpha=0.7,
                       rasterized=True,
                       zorder=999)
    if False:
        ax.scatter(Xn[0,:],
            average_freqs,
            s=[amp/average_freq_max_amp*80 for amp in average_freq_amps],
            c='None',
            marker='X',
            linewidths=0.5,
            edgecolors='w',
                    zorder=9999)

        ax.scatter(ax.freq_max_time,
                    ax.freq_max,
                    c=ax.freq_max_amp,
                    s=80,
                    marker='X',
                    norm=colors.LogNorm(*ax.figure.scalogram_lims),
                    linewidths=0.8,
                    edgecolors='k',
                    zorder=9999)

    ax.set_ylabel("Frequency (Hz)",position='right')
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(EngFormatter(places=0, sep="\N{THIN SPACE}"))

    if cb:
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
        cbaxes.yaxis.set_major_formatter(EngFormatter(unit=units,
                                                    places=0,
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
              annotation=False,
              zorder=999999999,
              ):

    for i,label in enumerate(vlinelabels):
        for t in vlinetimes[i]:

            lw=[3.5,2.5,1.5,2.5,1,1,1,1][i]
            color='C%d'%range(1,10)[i]
            ls=['-','--',':','-.',':',':'][i]

            if not annotation:#'rg' in vlinelabels[i] or alllines:
                a.axvline(t,
                        label=None,
                        color=color,
                        ls=ls,
                        linewidth=lw,
                        zorder=zorder,
                        dash_capstyle='round',
                        path_effects=[patheffects.withStroke(linewidth=lw+1.5, foreground="w")]
                        )
                a.axvline(t,
                        label=None,
                        color='w',
                        ls=ls,
                        linewidth=lw,
                        alpha=0.5,
                        zorder=9999,
                        dash_capstyle='round'
                        )
                a.axvline(t,
                        label=None,
                        color=color,
                        ls=ls,
                        linewidth=lw,
                        alpha=0.5,
                        zorder=9999,
                        dash_capstyle='round'
                        )
            else:

                for zorder in [-999,None]:
                    a.annotate('',
                                xy=(t, a.get_ylim()[1]),
                                xycoords='data',
                                xytext=(0, 1), textcoords='offset points',
                                arrowprops=dict(color=['w',color][zorder==None],
                                                headwidth=lw**.7*6,
                                                headlength=lw**.7*6,
                                                zorder=zorder,
                                                linewidth=[lw*2,None][zorder==None]
                                                ),
                                horizontalalignment='center',
                                verticalalignment='bottom',
                                zorder=zorder
                                )
                if False:
                    a.plot(t,y,
                            marker=4,
                            label=None,
                            markeredgecolor='w',
                            markeredgewidth=1,
                            color='C%d'%range(1,10)[vlinelabels.index(vlinelabels[i])]
                            )

def plottimelabels(a,
                   vlinelabels,
                   title,
                   alllines=False,
                   y=0,
                   offset=1+.03,
                   ):
    labeldone = []
    for i,label in enumerate(vlinelabels):
        if label in labeldone:
            continue
        lw=[3.5,2.5,1.5,2.5,1,1,1,1][i]
        color='C%d'%range(1,10)[i]
        ls=['-','--',':','-.',':',':'][i]
        a.axvline(nan,
                    label=label,
                    color=color,
                    ls=ls,
                    linewidth=lw,
                    zorder=999999999,
                    )
        labeldone += [label]

    a.plot(nan,nan,'vk',linewidth=0,label='Alert')

    fig = a.figure
    if hasattr(fig,'freq_max') and fig.freq_max > 0:
        a.scatter(nan,
                    nan,
                    s=90,
                    c='None',
                    marker='X',
                    linewidths=0.8,
                    edgecolors='k',
                    label=r'F$_{max}^{%.1f Hz}$ %.1f µm/s'%(fig.freq_max, fig.freq_max_amp*1000000))

    leg = a.legend(bbox_to_anchor=(1, offset),
                    loc='lower right',
                    ncol=len(vlinelabels)+2,
                    #fontsize='x-small',
                    title=title,
                    #title_fontsize='x-small',
                    #mode="expand",
                    borderaxespad=0.)
    a.add_artist(leg)

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
                        datalim=True,
                        orderbydistancefrom=None,#(14.73878,-91.56802),#Santiaguuito volcano
                        whitelist='*',#GI.STG*,CH.*,8D.*',
                        blacklist='*.BH?,*.HD?,*.BD?,*GI.STG2.*,GI.STG8.*,GI.STG5.*.HH?,GI.STG*SH?,GI.STG0.*',
                        nfreqs = 64,
                        ntimes = 512*2,
                        cmax=None,#0.00005,
                        cmin=None,#0.00000001,
                        pickorder=None,
                        spec_pick=False,
                        **args):

    if spec_pick is not None:
        spec_pick=int(spec_pick)

    if cmin is not None:
        cmin = float(cmin)
    if cmax is not None:
        cmax = float(cmax)

    nfreqs = int(nfreqs)
    ntimes = int(ntimes)

    if isinstance(all,str):
        all = bool(all)

    if isinstance(orderbydistancefrom,str):
        orderbydistancefrom = tuple([float(x) for x in orderbydistancefrom.split(',')])
        print(orderbydistancefrom)

    # EMPTY SPECTROGRAM
    # fix gain EHZ (399635995.8 c/m/s)

    blacklist = blacklist.split(',')
    whitelist = whitelist.split(',')
    #print(streams['raw'].__str__(extended=True))
    #print(event.preferred_origin())

    instruments = []
    instrumentcoordinates = []

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
                print(mseedid,'is not picked')
                continue

            if instrument in instruments:
                print(mseedid,'is already in')
                continue
            instruments += [instrument]
            c = inventory.select(network=trace.stats.network,
                                 station=trace.stats.station,
                                 location=trace.stats.location,
                                 channel=trace.stats.channel,
                                 time=trace.stats.starttime)[-1][-1][-1]
            instrumentcoordinates += [(c.latitude,c.longitude)]

    if all:
        for n in inventory:
            for s in n:
                for c in s:
                    mseedid = '%s.%s.%s.%s'%(n.code,s.code,c.location_code,c.code)
                    if mseedid[:-1] in instruments:
                        print(mseedid,'is already in')
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

    if orderbydistancefrom is not None:
        distances = [gps2dist_azimuth(*orderbydistancefrom,*c)[0] for c in instrumentcoordinates]
        instruments = [instruments[i] for i in argsort(distances)]

    if len(instruments) > max_num_station:
        distances = distances[:max_num_station]
        instruments = instruments[:max_num_station]

    #print(instruments)
    figs = [allstationsdataaxes(output,len(instruments),map=map) for output in outputs]

    for fig in figs:
        if cmin is not None:
            fig.scalogram_cmin = cmin
        if cmax is not None:
            fig.scalogram_cmax = cmax

    # Make all time limits equal to data span
    if datalim:
        alltimes = [t for output in outputs for tr in streams[output] for t in tr.times("matplotlib")]
        timelim = [ nanmin(alltimes),
                    nanmax(alltimes) ]
        print('Fixing time limits to: ',timelim)
        for f,fig in enumerate(figs):
            fig.cbdone = False
            for a in fig.axes_spec:
                a.set_xlim(timelim)
                a.freq_max_time = -9999

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
    #for o in event.origins:
    #    pref = ['','pref.'][str(preforig.resource_id) == str(o.resource_id)]

    # Earliest origin
    first = argmin([o.creation_info.creation_time for o in event.origins])
    first = event.origins[first]

    vlinelabels = []
    vlinetimes = []
    orgsymboltimes = {}
    orgsymbolpicks = {}
    if pickorder is not None:
        for p in pickorder.split(','):
           vlinelabels+=[p]
           vlinetimes+=[[]]
        orgsymboltimes = {k:[] for k in pickorder.split(',')}
        orgsymbolpicks = {k:[] for k in pickorder.split(',')}

    #if ot and len(pref):
    #    label = '$Org._{T.}^{pref.}$'
    #    vlinelabels += [label]
    #    vlinetimes += [[pref.time.matplotlib_date]]

    #if ew and len(first):
    #    label = '$Org._{CT.}^{early.}$'
    #    vlinelabels += [label]
    #    vlinetimes += [[first.creation_info.creation_time.matplotlib_date]]

    eventtype = 'Automatic'
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
        #print(instrument)
        #print(stationwf[outputs[0]])
        ncha = len(stationwf[outputs[0]]) #  ALLOW MULTIPLE INSTRUMENTS SAME STATION?

        for iout,(output, wf) in enumerate(stationwf.items()):
            if iout+1>len(figs):
                break

            if not len(wf):
                for a in [ figs[iout].axes_data[ax_index]]:#, figs[iout].axes_spec[ax_index] ]:

                    plotlabel = instrument.split('.')
                    plotlabel = '%s.%s'%(plotlabel[1],plotlabel[-1])

                    leg = a.legend([],[],
                                title=plotlabel,#'%s. %s'%(output.capitalize(),instrument),
                                #title_fontsize='xx-small',
                                labelspacing=-.2)
                    a.add_artist(leg)

                if ax_index == int(len(figs[iout].axes_data)/2):
                    figs[iout].axes_spec[ax_index].yaxis.set_label_position("left")
                    figs[iout].axes_spec[ax_index].set_ylabel("Frequency (Hz)",position='right')
                    figs[iout].axes_data[ax_index].yaxis.set_label_position("right")
                    figs[iout].axes_data[ax_index].set_ylabel('%s. (%s)'%( output.capitalize(), units[output] ))

            for i,label in enumerate(vlinelabels):
                if 'Org' == vlinelabels[i][1:4]:
                    continue
                vlinetimes[i] = []

            ismanual=False
            for o in event.origins:
                for a in o.arrivals:
                    p=a.pick_id.get_referred_object()
                    if '.'.join(instrument.split('.')[:-2]) == '%s.%s'%(p.waveform_id.network_code,p.waveform_id.station_code):
#                        label = '$Pick_{%s\ %s.}^{%s.%s}$' % ( p.phase_hint, str(o.method_id).split('/')[-1][:4], p.evaluation_mode[:4], pref )
                        label = '%s' % ( str(p.method_id).split('/')[-1][:4] )
                        if 'Fixe' in label or 'None' in label :
                            continue
                        if label not in vlinelabels:
                            vlinelabels += [label]
                            vlinetimes += [[]]
                        vlinetimes[vlinelabels.index(label)] += [p.time.matplotlib_date]

                        if 'Manu' in label:
                            ismanual=True
                            eventtype = 'Manual'

                        if p in orgsymbolpicks[label]:
                            continue
                        if label not in orgsymboltimes :
                            orgsymboltimes[label] = [p.time]
                            orgsymbolpicks[label] = [p]
                        else :
                            orgsymboltimes[label] += [p.time]
                            orgsymbolpicks[label] += [p]

            for itr,tr in enumerate(wf): #  ALLOW MULTIPLE INSTRUMENTS SAME STATION?

                #'%s. %s'%(output.capitalize(),tr.id)
                plotlabel = tr.id.split('.')
                plotlabel = '%s.%s'%(plotlabel[1],plotlabel[-1])

                # Seismogram
                a = figs[iout].axes_data[ax_index]

                plottimes(a,vlinelabels,vlinetimes,zorder=None)

                a.plot(tr.times("matplotlib"), tr.data)


                a.set_xlim([min(tr.times("matplotlib")),max(tr.times("matplotlib"))])

                leg = a.legend([],[],
                        title=plotlabel,
                        #title_fontsize='xx-small',
                        labelspacing=-.2)
                a.add_artist(leg)

                a.yaxis.set_label_position("right")
                if ax_index == int(len(figs[iout].axes_data)/2):
                    a.set_ylabel('%s. (%s)'%( output.capitalize(), units[output] ))

                a.yaxis.set_minor_locator(AutoMinorLocator())
                a.yaxis.set_major_locator(AutoLocator())
                a.yaxis.set_major_formatter(EngFormatter(places=0, sep="\N{THIN SPACE}"))


                # Scalogram
                a=figs[iout].axes_spec[ax_index]
                lax = scalogram(tr,
                                a,
                                units[output],
                                output,
                                nfreqs = nfreqs,
                                ntimes = ntimes,
                                cb=sharec==False,
                                spec_pick=spec_pick*ismanual)

                plottimes(a,vlinelabels,vlinetimes)

                if False:
                    leg = a.legend([],[],
                            title=plotlabel,
                            #title_fontsize='x-small',
                            labelspacing=-.2,
                            )
                    leg.set_zorder(99999)
                a.yaxis.set_label_position("left")
                if ax_index == int(len(figs[iout].axes_data)/2):
                    a.set_ylabel("Frequency (Hz)",position='right')
                else:
                    a.set_ylabel("")

                if sharec and not figs[iout].cbdone:

                    figs[iout].cbdone = True

                    nplots = len(figs[iout].axes_spec)
                    cbaxes = figs[iout].axes_spec[0].inset_axes([0.0,  1+.03*nplots,
                                                                 0.45, 0.015*nplots])

                    cbaxes.set_xticklabels([ ],
                                        #fontsize='xx-small',
                                        path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])  # vertically oriented colorbar
                    colorbar(a.scalogram,
                            cax=cbaxes,
                            orientation='horizontal')

                    cbaxes.tick_params(labeltop=True,
                                    labelbottom=False,
                                    top=True,
                                    bottom=False,
                                    which='both')

                    cbaxes.xaxis.tick_top()
                    cbaxes.xaxis.set_label_position('top')
                    cbaxes.set_xlabel('%s. (%s)'%( output.capitalize(), units[output] ),
                                      position='top',
                                      path_effects=[patheffects.withStroke(linewidth=4, foreground="w")])

                    cbaxes.xaxis.set_major_locator(LogLocator(base=10))
                    cbaxes.xaxis.set_major_formatter(EngFormatter(places=0,
                                                                  sep="\N{THIN SPACE}"))
                    a.scalogram_cbar = cbaxes

                # Seismogram
                a = figs[iout].axes_data[ax_index]
                a.set_ylim([median(tr.data)-fig.scalogram_cmax,median(tr.data)+fig.scalogram_cmax])
                a.set_xlim(figs[iout].axes_spec[ax_index].get_xlim())

                break

    if spec_pick:
        for fig in figs:

            index_last_max_freq = argmax([ a.freq_max_time for a in fig.axes_spec ])
            ax = fig.axes_spec[index_last_max_freq]

            if ax.freq_max_time < 0:
                continue

            fig.freq_max_time = ax.freq_max_time
            fig.freq_max = ax.freq_max
            fig.freq_max_amp = ax.freq_max_amp

            ax.scatter(ax.freq_max_time,
                       ax.freq_max,
                       c=ax.freq_max_amp,
                       s=90,
                       marker='X',
                       norm=colors.LogNorm(*ax.figure.scalogram_lims),
                       linewidths=0.8,
                       edgecolors='k',
                       zorder=999999)

            axes_data = fig.axes_data[index_last_max_freq]
            axes_data.scatter(ax.freq_max_time,
                       average(axes_data.get_ylim()),
                       c='None',#ax.freq_max_amp,
                       s=90,
                       marker='X',
                       #norm=colors.LogNorm(*ax.figure.scalogram_lims),
                       linewidths=0.8,
                       edgecolors='k',
                       zorder=999999)

    if preforig.time is None:
        preforig.time  = min([sort(orgsymboltimes[k])[0] for k in orgsymboltimes])

    orgsymboltimes_list = []
    for k in vlinelabels:
        if len(orgsymboltimes[k])>1 :
            orgsymboltimes_list += [[sort(orgsymboltimes[k])[1].matplotlib_date]]
        elif 'Flow' in k and len(orgsymboltimes[k]):
            orgsymboltimes_list += [[sort(orgsymboltimes[k])[0].matplotlib_date]]
        else:
            orgsymboltimes_list += [[]]

    for fig in figs:

        if event.event_type is None:
            event.event_type = "other event"

        for a in [fig.axes_spec[0], fig.axes_data[0]]:
            plottimes(a,
                      vlinelabels,
                      orgsymboltimes_list,
                    annotation=True)

        a = fig.axes_data[0]
        nplots = len(fig.axes_data)
        plottimelabels(a,
                       vlinelabels,
                       '%s event on %s'%(eventtype,#event.event_type.capitalize(),
                                         str(preforig.time)[:16]),
                        offset=1+.03*nplots)


    # Make all spec and seis plot tick and axis lines
    for f,fig in enumerate(figs):
        for grpindex,axgroup in enumerate([fig.axes_spec, fig.axes_data]): # order matters
            for ax_index,a in enumerate(axgroup):

                a.tick_params(right=grpindex!=0,
                              top=ax_index==0,
                              left=True,
                              bottom=(ax_index+1)==len(axgroup),
                              direction='out',
                              which='both')

                a.spines['bottom'].set_visible(ax_index+1==len(axgroup))
                a.spines['top'].set_visible(ax_index==0)

                locator = mdates.AutoDateLocator()
                formatter = mdates.ConciseDateFormatter(locator)

                a.xaxis.set_major_locator(locator)
                a.xaxis.set_major_formatter(formatter)

                a.grid(visible=True, which='major', color='gray', linestyle='dashdot', zorder=-999)
                a.grid(visible=True, which='minor', color='beige',  ls='-', zorder=-9999)


    return figs
