#!/usr/bin/env python
import imp,obspy,numpy
import matplotlib,matplotlib.animation
import obspy.signal.freqattributes,scipy.interpolate
from IPython.display import HTML
from datetime import datetime
from addons.core import scfinderauthor, nicemap, haversine, ipe_allen2012_hyp, getradius


g = 9.80665
golden = (1 + 5 **0.5) / 2

def eventmapboundary(latitudes,# = [origin.latitude, tr.stats.coordinates['latitude']],
                     longitudes,# = [origin.longitude, tr.stats.coordinates['longitude']]
                     radius=110000,
                     ):

    
    center = [numpy.median(longitudes),numpy.median(latitudes)]    
    longitude_radius = 0
    km=0
    while km<=radius:
        km,az,baz=obspy.geodetics.base.gps2dist_azimuth(center[1],center[0],center[1],
                                                        center[0]+longitude_radius)
        longitude_radius+=0.001

    latitude_radius = 0
    
    km=0
    while km<=radius:
        km,az,baz=obspy.geodetics.base.gps2dist_azimuth(center[1],center[0],
                                                       center[1]+latitude_radius,
                                                       center[0])
        latitude_radius+=0.001
        
    return [[center[0]+longitude_radius/golden-longitude_radius*golden,
             center[0]+longitude_radius/golden+longitude_radius*golden],
            [center[1]-latitude_radius,
             center[1]+latitude_radius]]

def default(ax): 

    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))

    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="white")]

    [l.set_path_effects(path_effects) for l in ax.get_xticklabels()]
    [l.set_path_effects(path_effects) for l in ax.get_yticklabels()]
    
def map_event(event,
              disp,
              stream,
              radius=110000,
              frames=10,
              vs=3.7,
              vp=3.7*2.7,
              xpixels=1500,
              label='Ground motion acceleration (%sg, observed)'%"%",#Intensity (Allen et al., 2012)',#
              figsize=[15,9],
              resolution='h',
              duration=60,
              **kwargs):
    
    
    magnitude = event.preferred_magnitude()
    if magnitude is None:
        print('no pref mag')
        event.preferred_magnitude = event.magnitudes[0].resource_id.get_referred_object
        print(event.preferred_magnitude())
    origin = event.preferred_origin()
    tr=disp.select(channel='*Z')[1]
    origin_time=origin.time
    
    mapbounds = eventmapboundary(latitudes=[origin.latitude],# = [origin.latitude, tr.stats.coordinates['latitude']],
                    longitudes=[origin.longitude],# = [origin.longitude, tr.stats.coordinates['longitude']]
                    radius=radius,
                    )

    fig,ax,bmap = nicemap(mapbounds=mapbounds,
                          labels=[1,0,0,1],
                          xpixels=xpixels,
                          resolution=resolution,
                          figsize=figsize,
                          **kwargs
                          )
    fig.tight_layout(rect=(0, 0, 1, .97))
    fig.event_creation_info = event.creation_info
    fig.creation_time = obspy.UTCDateTime.now()

    xlims=ax.get_xlim()
    ylims=ax.get_ylim()
    patchheight=1/15 
    patchpad=1/100
    ax.add_patch(matplotlib.pyplot.Rectangle((xlims[0]+(xlims[1]-xlims[0])/golden,
                                              ylims[1]+(ylims[1]-ylims[0])*patchpad),
                                             (xlims[1]-xlims[0])-(xlims[1]-xlims[0])/golden, 
                                             (ylims[1]-ylims[0])*patchheight,
                                             facecolor='w',
                                             zorder=999,
                                             clip_on=False,
                                             linewidth = 0))
    ax.add_patch(matplotlib.pyplot.Rectangle((xlims[0]+(xlims[1]-xlims[0])/golden,
                                              ylims[0]-(ylims[1]-ylims[0])*patchpad-(ylims[1]-ylims[0])*patchheight),
                                             (xlims[1]-xlims[0])-(xlims[1]-xlims[0])/golden, 
                                             (ylims[1]-ylims[0])*patchheight,
                                             facecolor='w',
                                             zorder=999,
                                             clip_on=False,
                                             linewidth = 0))
    
    region=', '.join([desc.text for desc in event.event_descriptions if 'region' in desc.type])
    desc='%s - %s: M$_{%s}$%.1f, %.1fkm deep'
    desc=desc%(region,
               str(event.preferred_origin().time)[:19],
               event.preferred_magnitude().magnitude_type[1:],
               event.preferred_magnitude().mag,
               event.preferred_origin().depth/1000)
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                  foreground="white")]
    ax.set_title(desc, 
                 ha='left',
                 va='bottom',
                 pad=2,
                 x=.0, y=1,
                 color='k',
                 fontweight='bold',
                 path_effects=path_effects)    
    
    fig.bmap=bmap
    fig.vs=vs
    fig.vp=vp
    fig.frames=frames
    fig.duration=duration
    fig.secperframe = fig.duration/fig.frames
    fig.deglatatlon = haversine(origin.longitude,
                                   origin.latitude-.5,
                                   origin.longitude,
                                   origin.latitude+.5)/1000
    
    
    

    fig.insets = []
    pos = ax.get_position()

    ## PSA subplot
    w = pos.width*(golden-int(golden))/golden
    h = pos.height*(golden-int(golden))/golden
    position = [pos.x1-w, pos.y0, w, h]
    ax = fig.add_axes(position,
                      facecolor=[0,0,0,.5],
                      zorder=999)
    default(ax)
    ax.yaxis.set_label_position('right')
    ax.tick_params(right=True, 
                   top=False,
                   left=True, 
                   bottom=True,
                   labelbottom=True, 
                   labeltop=False, 
                   labelleft=False, 
                   labelright=True,
                   which='both')
    fig.insets += [ax]
    
    ## GM-distance subplot
    h = pos.height*int(golden)/golden
    position = [pos.x1-w, pos.y1-h, w, h]
    axmmi = fig.add_axes(position,
                       facecolor=[1,1,1,.5],
                            zorder=999)
    axmmi.xaxis.set_label_position('top') 
    axmmi.twin = axmmi.twinx()
    axmmi.twin.set_zorder(axmmi.get_zorder() + 1)

    default(axmmi)
    axmmi.tick_params(right=False, 
                   top=True,
                   left=True, 
                   bottom=False,
                   labelbottom=False, 
                   labeltop=True, 
                   labelleft=True, 
                   labelright=False,
                    which='both')

    axmmi.twin.tick_params(right=True, 
                   top=False,
                   left=False, 
                   bottom=False,
                   labelbottom=False, 
                   labeltop=False, 
                   labelleft=False, 
                   labelright=False,
                   which='both')
    axmmi.twin.set_yscale('log')
    axmmi.set_xlim([0,#event.preferred_origin().depth/1000, 
                 ((duration*vs*2)**2-(event.preferred_origin().depth/1000)**2)**.5])
    axmmi.spines['left'].set_color('C0')
    axmmi.tick_params(axis='y', colors='C0',which="both")
    axmmi.yaxis.label.set_color('C0')
    fig.insets += [axmmi]

    ## Colorbar
    cw = pos.width*1/128
    position = [pos.x1-cw/2, pos.y1-h, cw, h]
    fig.cax = fig.add_axes(position,
                           facecolor=[0,0,0,0],
                           sharey=axmmi.twin,
                           zorder=9999)
    fig.cax.tick_params(right=False, 
                   top=False,
                   left=False, 
                   bottom=True,
                   labelbottom=False, 
                   labeltop=False, 
                   labelleft=False, 
                   labelright=True,
                   which='both')
    default(fig.cax)

    
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                  foreground="white")]
    opts={'x':origin.longitude,
          'y':origin.latitude,
          'marker':'*',
          'markeredgecolor':'k',
          'markerfacecolor':'None',
          'markersize':20,
          'latlon':True,
          'zorder':999,
          'path_effects':path_effects}
    fig.preferred_magnitude = [magnitude,fig.bmap.plot(**opts)[0]]
    
    inputs=[origin.longitude,
            origin.latitude,
            0 / fig.deglatatlon,
            64]
    opts={'fill':False,
          'edgecolor':'k'}
    copts={'linewidth':2,
          'path_effects':path_effects}
    fig.swaves = []

    print(stream)
    cmap = matplotlib.cm.nipy_spectral
    normopts = {'vmin': min([numpy.percentile(abs(tr.data[tr.data!=0]/0.0981),8) for tr in stream]), 
                'vmax': max([max(abs(tr.data[tr.data!=0]/0.0981)) for tr in stream])*1.1
                }
    norm = matplotlib.colors.LogNorm(**normopts)
    scalarcolormap = matplotlib.cm.ScalarMappable(norm=norm, 
                                                  cmap=cmap)
    cb = fig.colorbar(scalarcolormap,
                      cax=fig.cax,
                      ticklocation='right', 
                      extend='both',
                      label=label)
    default(cb.ax)
    cb.ax.minorticks_on()
    #cb.ax.set_yticks(range(13))
    #cb.ax.set_yticklabels(['','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII'])
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="white")]  
    minpath_effects = [matplotlib.patheffects.withStroke(linewidth=0,
                                                         foreground="white")]  
    fig.cax.set_ylabel(label,
                       path_effects=path_effects,
                       fontweight='bold')
    axmmixlim = axmmi.twin.get_xlim()
    for index in range(fig.frames):
        frametime = index*fig.secperframe
        linestyle='solid'
        for v in [fig.vs, fig.vp]:
            hypd=frametime*v
            if hypd<=abs(origin.depth)/1000:
                continue
            inputs[2]=((hypd**2-(abs(origin.depth)/1000)**2)**.5)/fig.deglatatlon
            if inputs[2]<=0:
                continue
            fig.swaves+=[[frametime,
                        frametime+fig.secperframe,
                        fig.bmap.tissot(*inputs,linestyle=linestyle,**opts,**copts)]]
            fig.swaves+=[[frametime,
                        frametime+fig.secperframe,
                        fig.insets[1].axvline(hypd,color='k',linestyle=linestyle,**copts)]]
            linestyle = 'dashed'


            for tr in stream.slice(endtime=origin_time+frametime):
                
                tropts={'x':tr.stats.coordinates['longitude'],
                    'y':tr.stats.coordinates['latitude'],
                    'linestyle':'None',
                    'latlon':True,
                    'marker':'^',
                    'markersize':8,
                    'markeredgewidth':.5,
                    'markeredgecolor':'k',
                    'color':scalarcolormap.to_rgba(max(abs(tr.data/0.0981))),#fig.acc2mmi(),
                    'path_effects':path_effects}
                fig.swaves[-1] += [fig.bmap.plot(**tropts,zorder=99999)[0]]
                
                tropts['markersize']=6
                tropts.pop("latlon")
                tropts.pop("x")
                tropts.pop("y")
                tropts.pop("path_effects")
                fig.swaves[-1] += [axmmi.twin.plot(tr.stats.distance/1000,
                                            max(abs(tr.data/0.0981)),#fig.acc2mmi(),
                                            **tropts,
                                            zorder=9999)[0]]
                tropts['markersize'] = 0
                tropts['linestyle'] = '-'
                tropts['linewidth'] = 1
                tropts['color'] = 'w'
                fig.swaves[-1] += [axmmi.twin.plot([tr.stats.distance/1000,axmmixlim[1]],
                                            [max(abs(tr.data/0.0981)),max(abs(tr.data/0.0981))],#fig.acc2mmi(),
                                            **tropts,
                                            path_effects=None,
                                            zorder=9998)[0]]
    return fig

def samemagupdate(m,mm,dt,origin_time):
    test = m.creation_info.creation_time - origin_time
    if (mm.magnitude_type == m.magnitude_type and 
            m.creation_info.author == mm.creation_info.author and
            test>dt):
        return test
    return 99999

def add_magnitudes(event,
                   bmap,
                   axmmi,
                   stream,
                   label='Ground motion acceleration (%sg, observed)'%"%",#Intensity (Allen et al., 2012)',#
                   lineauthors=[],
                   types = {'Mfd':'C1','MVS':'C2'},
                   latextypes = {'Mfd':'orange','MVS':'green'},
                   blacklist = ['scvsmag2']):

    fig = bmap.ax.figure
    fig.magnitude_mmi=[]
    pref_origin = event.preferred_origin()
    pref_magnitude = event.preferred_magnitude()
    origin_time=pref_origin.time
   
    ###############################################################
    lineauthors = scfinderauthor(pref_origin,lineauthors=lineauthors)
    ###############################################################
        
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="white")]
    fig.magnitudes = []
    fig.title_info=[]
    magnitudes=[]
    for m in event.magnitudes:
        if m.creation_info.author is None:
            continue
        if (m.creation_info.author.split('@')[0] in blacklist or 
            m.magnitude_type in blacklist or
            m.magnitude_type not in types or
            (m.magnitude_type == 'Mfd' and 
             m.creation_info.author.split('@')[0] not in lineauthors)):
            continue
        dt=m.creation_info.creation_time - origin_time
        fig.magnitudes+=[[dt]]  
        magnitudes+=[m]
        
    station_distances=[]
    for tr in stream:
        station_distances+=[tr.stats.distance/1000]
    station_distances=numpy.linspace(min(station_distances),max(station_distances+list(axmmi.get_xlim())),16)
    
    for i,m in enumerate(magnitudes):
        origin = m.origin_id.get_referred_object()
        dt=m.creation_info.creation_time - origin_time
        inputs=[m,dt,origin_time]
        update=min([ samemagupdate(mm,*inputs) for mm in magnitudes])
        fig.magnitudes[i]+=[update]
        #desc='$\color{%s}M_{%s}%.f$'
        #desc=desc%(latextypes[m.magnitude_type],m.magnitude_type,m.mag)
        desc='M$_{%s}$%.1f'
        desc=desc%(m.magnitude_type[1:],m.mag)
        desc={'color':types[m.magnitude_type],'label':desc}
        fig.title_info+=[[dt,update,desc]]
        
        opts={'x':origin.longitude,
              'y':origin.latitude,
              'marker':'*',
              'color':types[m.magnitude_type],
              #'markeredgecolor':types[m.magnitude_type],
              #'markerfacecolor':'none',
              'markersize':20*(m.mag**2)/(pref_magnitude.mag**2),
              'latlon':True,
              'zorder':9999,
              'path_effects':path_effects}
        fig.magnitudes[i]+=[fig.bmap.plot(**opts)[0]]
        
        minputs=[numpy.array([origin.depth/1000]),
                 numpy.array([m.mag])]
        mmi=list(ipe_allen2012_hyp(*minputs))
        distances=[abs(origin.depth)/1000]
        for tr in stream.slice(endtime=m.creation_info.creation_time):
            minputs=[origin.longitude,origin.latitude,
                     tr.stats.coordinates['longitude'],
                     tr.stats.coordinates['latitude']]
            epd=haversine(*minputs)/1000
            hyp=(epd**2+(abs(origin.depth)/1000+tr.stats.coordinates['elevation']/1000)**2)**.5
            distances+=[hyp]

        minputs=[numpy.linspace(min(distances),max(distances+list(axmmi.get_xlim())),16),
                 numpy.array([m.mag for n in range(16)])]
        minputs[1]= ipe_allen2012_hyp(*minputs)
        minputs[0]= station_distances
        mopts={'color':types[m.magnitude_type],
               #'linestyle':'dashed',
               'linewidth':2,
               'path_effects':path_effects}
        fig.magnitudes[i]+=[axmmi.plot(*minputs,**mopts)[0]]
        minputs[1]= [minputs[1][0]]+[v for v in minputs[1]]+[minputs[1][-1]] 
        minputs[0]= [0]+[v for v in minputs[0]]+[0]
        mopts['linewidth']=1
        mopts['alpha']=.5
        mopts.pop("path_effects")
        fig.magnitudes[i]+=[axmmi.plot(*minputs,path_effects=None,**mopts)[0]]

            
        inputs=[origin.longitude,
                origin.latitude,
                0 / fig.deglatatlon,
                64]
        opts={'fill':False,
              'edgecolor':types[m.magnitude_type],
              #'linestyle':'dashed',
              'linewidth':2,
              'path_effects':path_effects}  
        for index in range(fig.frames):
            frametime = index*fig.secperframe
            if frametime<dt: 
                continue
            if frametime>=update:
                continue
            linestyle = None
            for v in [fig.vs]:#, fig.vp]:
                hypd=frametime*v
                if hypd<=origin.depth/1000:
                    continue
                inputs[2]=((hypd**2-(origin.depth/1000)**2)**.5)/fig.deglatatlon
                if inputs[2]<=0:
                    continue
                fig.swaves+=[[frametime,
                              frametime+fig.secperframe,
                              fig.bmap.tissot(*inputs,
                                              linestyle=linestyle,
                                              **opts)]]
                linestyle = 'dashed'

     

    ## Align MMI and acc
    #http://shakemap.rm.ingv.it/shake/17769831/download/intensity.jpg
    #http://www.seismo.ethz.ch/static/shakemap/shake/KP201211290852/intensity.html
    #modAcc = [ 0.06/10,0.06/2, 0.2, .8, 2, 4.8, 12, 29, 70, 171]
    #modMMI = [       0,     1, 2.5,  4, 5,   6,  7,  8,  9,  10]
    int2acc = [[ 1, 0.06/2], 
               [ 4,  .8   ], 
               [ 5, 2     ]]

    mmilim = axmmi.get_ylim()
    acclim = axmmi.twin.get_ylim()
    ha = (numpy.log10(acclim[1])-numpy.log10(int2acc[0][1]))/(numpy.log10(acclim[1])-numpy.log10(acclim[0]))
    bottom = -1*(mmilim[1]-int2acc[0][0]-ha*mmilim[1])/ha 
    ## linked mmi & acc : axmmi.set_ylim(bottom=bottom)
    
    mmilim = axmmi.get_ylim()
    acclim = axmmi.twin.get_ylim()
    ## linked mmi & acc : axmmi.set_ylim(top=max([fig.acc2mmi(acclim[1]), mmilim[1]]))
    ## linked mmi & acc : axmmi.twin.set_ylim(top=max([fig.mmi2acc(mmilim[1]), acclim[1]]))

    ## Fill MMI plot
    mmilim = axmmi.get_ylim()
    if mmilim[0]<fig.mmifillmin:
        axmmi.fill_between(axmmi.get_xlim(),
                        [fig.mmifillmin,fig.mmifillmin],## linked mmi & acc :mmilim[0],mmilim[0]],
                        fig.mmifillmin,
                        **fig.mmifillopt)    
        ## linked mmi & acc : axmmi.set_ylim(mmilim)    
    
def add_psa(tr,ax,time, tr1=None, tr2=None, maxperiod=5):  
    path_effects = [matplotlib.patheffects.withStroke(linewidth=10,
                                                      foreground="white")]  
    ax.set_xlabel('Period (s)',path_effects=path_effects,fontweight='bold',)
    ax.set_ylabel('PSA (%g, observed)',path_effects=path_effects,fontweight='bold',)
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="k")]
    ax.set_title('%s %.3g km ep.$_{\ }$'%(tr.id[:-1],tr.stats.distance/1000), 
                 ha='right',
                 va='top',
                 pad=-2,
                 x=1, y=1,
                 color='w',
                 fontweight='bold',
                path_effects=path_effects)
    ax.set_xscale('log')
    minor = matplotlib.ticker.LogLocator(base = 10.0, 
                                         subs = numpy.arange(1.0, 10.0) * 0.1,
                                         numticks = 10)
    ax.xaxis.set_minor_locator(minor)

    opts = [tr.data,tr.stats.delta]
    frequencies = numpy.logspace(numpy.log10(1/maxperiod),numpy.log10(tr.stats.sampling_rate/2),32)
    ax.set_xlim([1/max(frequencies),1/min(frequencies)])
    ax.set_ylim(bottom=0)
    psa = []
    for freq in frequencies:
        psa += [obspy.signal.freqattributes.peak_ground_motion(*opts,freq)[0]/g*100]
    
    if tr1 is not None:
        opts[0]=tr1.data
        for i,freq in enumerate(frequencies):
            psa[i] = max([obspy.signal.freqattributes.peak_ground_motion(*opts,freq)[0]/g*100, psa[i]])
    if tr2 is not None:
        opts[0]=tr2.data
        for i,freq in enumerate(frequencies):
            psa[i] = max([obspy.signal.freqattributes.peak_ground_motion(*opts,freq)[0]/g*100, psa[i]])

    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="white")]
    inputs=[1/frequencies,psa]
    opts={'color':'C0',
          'path_effects':path_effects}
    fopts={'color':'C0',
           'alpha':.5,
           'linewidth':0}
    ax.figure.psa = [(tr,
                      tr,
                      ax.plot(*inputs,**opts)[0],
                      ax.fill_between(*inputs,**fopts))]
    
    pinputs = [tr.data,tr.stats.delta]    
    opts={'color':'C0',
          'alpha':1,
          'linestyle':'dashed',
          'path_effects':path_effects}
    ax.figure.psa_evolution =[]
    for index in range(ax.figure.frames):
        frametime = index*ax.figure.secperframe 
        pinputs[0] = tr.slice(starttime=time-1/min(frequencies),endtime=time+frametime).data
        psa = []
        for freq in frequencies:
            if len(pinputs[0]):
                psa += [obspy.signal.freqattributes.peak_ground_motion(*pinputs,freq)[0]/g*100]
            else:
                psa += [0]
        
        if tr2 is not None:
            pinputs[0] = tr2.slice(endtime=time+frametime).data
            for i,freq in enumerate(frequencies):
                if len(pinputs[0]):
                    psa[i] = max([obspy.signal.freqattributes.peak_ground_motion(*pinputs,freq)[0]/g*100, psa[i]])
        if tr1 is not None:
            pinputs[0] = tr1.slice(endtime=time+frametime).data
            for i,freq in enumerate(frequencies):
                if len(pinputs[0]):
                    psa[i] = max([obspy.signal.freqattributes.peak_ground_motion(*pinputs,freq)[0]/g*100, psa[i]])

        inputs=[1/frequencies,psa]
        ax.set_ylim(top=max([ax.get_ylim()[1],max(psa)]))
        ax.figure.psa_evolution += [[frametime,
                                     frametime+ax.figure.secperframe,
                                     ax.plot(*inputs,**opts)[0]]]
                   
    path_effects = [matplotlib.patheffects.withStroke(linewidth=7,
                                                      foreground="k")]
    opts={'x':tr.stats.coordinates['longitude'],
          'y':tr.stats.coordinates['latitude'],
          'latlon':True,
          'marker':'^',
          'markersize':8,
          'markeredgewidth':.5,
          'markeredgecolor':'k',
          'color':'w',
          'path_effects':path_effects}
    ax.figure.psa_stations = [(tr,tr,ax.figure.bmap.plot(**opts)[0])]
        
def add_mmi(stream,event,ax):
    fig = ax.figure
    pref_origin = event.preferred_origin()
    pref_magnitude = event.preferred_magnitude()
    
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="w")]
    ax.set_title('Ground motion distribution$^{\ }$.', 
                 ha='right',
                 va='top', 
                 pad=-2,
                 x=1, y=1,
                 color='k',
                 fontweight='bold',
                 path_effects=path_effects)
    ax.set_xlabel('Distance (km, hyp.)',fontweight='bold')#,path_effects=path_effects)
    ax.set_ylabel('Ground motion intensity (MMI, modelled)',fontweight='bold',path_effects=path_effects)
    ax.set_yticks(range(13))
    ax.set_yticklabels(['','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII'])
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      foreground="w")]
    
    inputs=[numpy.array([pref_origin.depth/1000]),
            numpy.array([pref_magnitude.mag])]
    mmi=list(ipe_allen2012_hyp(*inputs))
    distances=[pref_origin.depth/1000]
    lon=[pref_origin.longitude]
    lat=[pref_origin.latitude]
    ids=['Epicentre']
    modAcc=[]
    for tr in stream:
        lon+=[tr.stats.coordinates['longitude']]
        lat+=[tr.stats.coordinates['latitude']]
        ids+=[tr.id]
        distances+=[tr.stats.distance/1000]
        modAcc+=[max(abs(tr.data/0.0981))]

    modMMI=list(ipe_allen2012_hyp(numpy.array(distances[1:]),
                                                                numpy.array([pref_magnitude.mag for n in distances[1:]]),
                                                                ))
    
    #http://shakemap.rm.ingv.it/shake/17769831/download/intensity.jpg
    #http://www.seismo.ethz.ch/static/shakemap/shake/KP201211290852/intensity.html
    modAcc += [ 0.06/10,0.06/2, 0.2, .8, 2, 4.8, 12, 29, 70, 171]
    modMMI += [       0,     1, 2.5,  4, 5,   6,  7,  8,  9,  10]
    fig.mmi2acc = scipy.interpolate.interp1d(modMMI,modAcc,fill_value="extrapolate")#,kind='quadratic')
    fig.acc2mmi = scipy.interpolate.interp1d(modAcc,modMMI,fill_value="extrapolate")#,kind='quadratic')
    #print('Interpolation data:')
    #print('x:',modMMI)
    #print('y:',modAcc)
    #print('Interpolation consistency check:')
    #print('y:',fig.mmi2acc_fmodel(modMMI))

    inputs=[numpy.linspace(min(distances),max(distances+list(ax.get_xlim())),16),
            numpy.array([pref_magnitude.mag for n in range(16)])]
    inputs[1] = ipe_allen2012_hyp(*inputs)
    opts={'color':'C0',
          'path_effects':path_effects}
    fopts={'color':'C0',
           'alpha':.5,
           'linewidth':0}
    ax.figure.mmi = [(tr,tr,
                      ax.plot(*inputs,**opts)[0],
                      ax.fill_betweenx(*inputs[::-1],0,**fopts))]
    ax.figure.mmifillmin = min(inputs[1])
    ax.figure.mmifillopt = fopts
    path_effects = [matplotlib.patheffects.withStroke(linewidth=4,
                                                      foreground="w")]
    opts={'x':lon[1:],
          'y':lat[1:],
          'linestyle':'None',
          'latlon':True,
          'marker':'^',
          'markersize':8,
          'markeredgewidth':.5,
          'markeredgecolor':'k',
          'color':'w',
          'path_effects':path_effects}
    ax.figure.mmi_stations = [(tr,tr,ax.figure.bmap.plot(**opts)[0])]

def init_func():
    print('hello')

def set_path_effects(artist,alpha):
    path_effects = [matplotlib.patheffects.withStroke(linewidth=5,
                                                      alpha=alpha,
                                                      foreground="white")]
    artist.set_path_effects(path_effects)

def set_alpha(artist,dt,et,frametime):
    if artist.get_alpha() is None:
        artist.set_alpha(1)
    if dt>frametime:
        artist.set_alpha(0)
        if artist.get_path_effects() is not None:
            set_path_effects(artist,0)
    elif  dt<=frametime and et>frametime:#artist.get_alpha()==0:
        artist.set_alpha(1)
        if artist.get_path_effects() is not None:
            set_path_effects(artist,1)
    else:
        if artist.get_path_effects() is not None:
            set_path_effects(artist,(artist.get_alpha()/2)**4+0.000001)
        artist.set_alpha((artist.get_alpha()/2)**4+0.000001)
        
def set_title(info,dt,et,frametime,ax):
    inputs=[[numpy.nan,numpy.nan],[numpy.nan,numpy.nan]]
    handles=[]
    labels=[]
    try:
        handles, labels = ax.my_legend_handles_labels
    except:
        pass
    old = labels
    done=False
    for l,lab in enumerate(labels):
        if 'wave' in lab:
            done=True
            handles[l]=ax.plot(*inputs,color='k')[0]
            labels[l] = 'S-wave at origin time +%ds'%frametime
            #handles+=ax.plot(*inputs,color='k',linestyle='dashed')
            #labels+=['P-wave']
    if not done:
        handles+=ax.plot(*inputs,color='k')
        labels+=['S-wave at origin time +%ds'%frametime]
        #handles+=ax.plot(*inputs,color='k',linestyle='dashed')
        #labels+=['P-wave']
        done=True
    
    if  dt<=frametime and et>=frametime: 
        done=False
        for l,lab in enumerate(labels):
            if info['label'][4:7] in lab:
                #print('overwrite',info['label'],'onto',lab)
                handles[l]=ax.plot(*inputs,**info)[0]
                labels[l]=info['label']
                done=True
        if not done:
            #print('add',info['label'],'to',labels)
            handles+=ax.plot(*inputs,**info)
            labels+=[info['label']]
            done=True
            
    
    ax.my_legend_handles_labels=(handles, labels)
    ax.legend(handles,
              labels,
              title='Data from %s updated at %s'%(ax.figure.event_creation_info.agency_id,str(ax.figure.creation_time)[:19]),
              bbox_to_anchor=(0, 1, 1, 0), 
              loc="upper left", ncol=len(labels))
            
def refresh(index,fig):
    frametime = index*fig.secperframe
    print('%4gs (%d/%d)'%(frametime,index, fig.frames))
    
    for t in fig.title_info:
        set_title(t[-1],t[0],t[1],frametime,fig.bmap.ax)
            
    for m in fig.magnitudes:
        for artist in m[2:]:
            set_alpha(artist,m[0],m[1],frametime)
    
    for s in fig.swaves:
        for artist in s[2:]:
            set_alpha(artist,s[0],s[1],frametime)
            
    for p in fig.psa_evolution:
        for artist in p[2:]:
            set_alpha(artist,p[0],p[1],frametime)

def animate(event,
            acc,
            disp,
            inventory,
            duration=30,
            frames=30*6,
            **opts):

    if frames is None:
        frames=duration*6

    radius = getradius(event.preferred_origin(),inventory)    

    fig = map_event(event,
                    disp,
                    acc.select(channel='*b')+acc.select(channel='*X'),
                    #epsg=3857,
                    frames=frames,
                    duration=duration,
                    radius=radius,
                    **opts)
    add_mmi(acc.select(channel='*b')+acc.select(channel='*X'),
                event,
                fig.insets[1])

    add_psa(disp[0],
        fig.insets[0],
        event.preferred_origin().time,
        tr1=disp[1],
        tr2=disp[2])

    lineauthors=scfinderauthor(event.preferred_origin())
    add_magnitudes(event,
                   fig.bmap,
                   fig.insets[1],
                   acc.select(channel='*Z'),
                   lineauthors=lineauthors)

    animopts = {'init_func':init_func,
                'frames':fig.frames,
                'fargs':[fig], 
                'interval':150,#fig.secperframe/1000, 
                'blit':False,
                'repeat':False}
    dat = [fig, refresh]

    return matplotlib.animation.FuncAnimation(*dat,**animopts) #anim = 

    ##HTML/(anim.to_html5_video())
    #anim.save('%s/%s.mp4'%(archive,shorteventid),dpi=300)
    #toopen='file://%s/%s.mp4'%(archive,shorteventid)
    #print(toopen)


