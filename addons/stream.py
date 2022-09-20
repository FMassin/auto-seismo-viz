#!/usr/bin/env python
from gettext import Catalog
import matplotlib.pyplot,matplotlib.colors,matplotlib.patheffects,numpy
from obspy.taup import TauPyModel
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def tracedistance(tr,
                  inventory,
                  origin,
                  optline=None,
                  onedeg=111.19492664455875):
    #print(tr)
    #print(inventory)
    try: #if waveform_id
        station=inventory.select(network=tr.network_code,
                                   station=tr.station_code)[0][0]
    except: #if trace
        try:
            station=inventory.select(network=tr.stats.network,
                                     station=tr.stats.station,
                                     location=tr.stats.location)[0][0]
        except:
            try:
                station=inventory.select(network=tr.stats.network,
                                        station=tr.stats.station)[0][0]
            except:
                return 9999,9999
    epr = obspy_addons.haversine(origin.longitude,
                                   origin.latitude,
                                   station.longitude,
                                   station.latitude)/1000
    if optline is not None:
        epr = obspy_addons.DistancePointLine(station.longitude,
                                                 station.latitude,
                                                 **optline)
        epr = epr*onedeg

    r = (epr**2+(origin.depth/1000+station.elevation/1000)**2)**.5
    return epr, r

def get_velmodel_correction(event,
                            inventory,
                            lim=999,
                            model=TauPyModel(model='iasp91')):

    time = event.preferred_origin().time
    depth = event.preferred_origin().depth/1000
    logs=[]
    modelresidual={}#print(inventory)
    for a in event.preferred_origin().arrivals:
        phase=a.phase[0].lower()
        p=a.pick_id.get_referred_object()
        if p is None:
            print(a,'misses related pick',a.pick_id)
            continue
        obstt = p.time - time
        epr,r=tracedistance(p.waveform_id,
                            inventory,
                            event.preferred_origin())
        if epr>lim:
            #print(a)
            continue
        arrivals = model.get_travel_times(depth,
                                          distance_in_degree=epr/111,
                                          phase_list=['tt'+phase],
                                          receiver_depth_in_km=0)
        modtt=numpy.nanmin([numpy.nan]+[a.time for a in arrivals])
        if phase in modelresidual:
            modelresidual[phase]+=[obstt - modtt]
        else:
            modelresidual[phase]=[obstt - modtt]
        logs += ['%s %s:%.1fs-%.1fs @ %dkm'%(p.waveform_id.station_code,
                                       a.phase,
                                       obstt,modtt,
                                       epr)]
    modelresidual = {k:numpy.percentile(modelresidual[k],50)for k in modelresidual}
    for phase in 'ps':
        if phase not in modelresidual:
            modelresidual[phase]=0
    print("model residual:",modelresidual)
    return modelresidual,model

def ploteqdata(self, #eqdata[output].select(channel='*b')
               catalog,#eqdata['catalog']
               inventory,#eqdata['inventory']
               model = TauPyModel(model='iasp91'),
               gain=1,
               minsnr=9,
               cmap='nipy_spectral',#tab20b'
               mtypes={'MVS':'C2','Mfd':'C1'},
               ax=None,
               lim=30,
               nopgalabel=False,
               nostationlabel=False,
               ):
    if ax is None:
        ax = matplotlib.pyplot.figure().gca()
    if hasattr(catalog,'events'):
        event = catalog[0]
    else:
        event = catalog

    time = event.preferred_origin().time
    depth = event.preferred_origin().depth/1000
    modelresidual,model=get_velmodel_correction(event,inventory,lim,model)
    tmax=0
    vmax=[]
    vmin=[]
    for tr in self:
        imax=numpy.nanargmax(tr.data)
        tmax=numpy.nanmax([tmax, tr.times(type="utcdatetime")[imax]])
        vmax+=list(tr.data)
    for tr in self:
        vmin+=list(tr.data[:99])
    print(vmin,vmax)
    vmin=numpy.nanpercentile(vmin,84)
    vmax=numpy.nanpercentile(vmax,100)*1.1
    print(vmin,vmax)
    norm=matplotlib.colors.LogNorm(vmin=vmin*gain,
                                   vmax=vmax*gain)
    cmap=matplotlib.pyplot.get_cmap(cmap)
    eprs=[]
    rs=[]
        
    skipped=[]
    withStroke = matplotlib.patheffects.withStroke
    for tr in self:
        if max(tr.data)<vmin*minsnr:
            continue
        epr,r=tracedistance(tr,
                      inventory,
                      event.preferred_origin())
        if epr>=lim:
            skipped+=[tr.id]
            continue
        eprs+=[epr]
        rs+=[r]
        imax=numpy.argmax(tr.data)
        times=tr.times(type="utcdatetime")[:imax+1]
        d=[epr for v in times]
        y=[numpy.nanmax(tr.data[:i+1])*gain for i,v in enumerate(tr.data[:imax+1])]
        s=[1 for v in times]
        t=[v-time for v in times]
        skip=False
        for test in [t,y,d,s]:
            if True in numpy.isnan(numpy.array(t)):
                print('skipping')
                skip=True
        if skip:
            continue
        order = numpy.argsort(t)[::-1]
        t = [t[i] for i in order]
        y = [y[i] for i in order]
        for shift in numpy.arange(tr.stats.delta/2,tr.stats.delta,0.1):
            h=ax.scatter(t+shift,d,s,y,marker='_',#'.',
                        cmap=cmap,
                        zorder=4+y[0],
                        norm=norm)
        x=t[0]
        if True:#x<max(ax.get_xlim()):
            label_path_effects=[withStroke(linewidth=5,
                                           foreground=cmap(norm(y[0]))),
                                withStroke(linewidth=3,
                                           foreground="w")]
            data_label='  %.2g m/s$^2$ %s'%(y[0],tr.id)
            if nostationlabel:
                data_label='  %.2g m/s$^2$'%(y[0])
            if nopgalabel:
                data_label='  %s'%(tr.id)
            if not nopgalabel or not nostationlabel:
                ax.text(x,
                        epr,#r,tr.data[imax],
                        data_label,#),#c[imax]
                        ha='left',
                        va='center',
                        clip_on=True,
                        size='xx-small',
                        fontweight='semibold',
                        path_effects=label_path_effects,
                        zorder=5+y[0],
                        #color='w',
                        #bbox={'fc': '0.8', 'pad': 2},
                        #rotation=90
                        )
    if len(skipped):
        print('Skipped further than %d km: '%(lim)+', '.join(skipped))
    opts={'markersize':15**.5,'color':'r','markeredgewidth':0}
    #ax.plot(numpy.nan,1,label='PGA$_{obs}$',**opts)
    
    cax = inset_axes(ax,
                    width="50%",  # width = 50% of parent_bbox width
                    height="3%",  # height : 3%
                    loc='upper right',
                    borderpad=0)
    h=matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
    cb=matplotlib.pyplot.colorbar(h,
                                  cax=cax,
                                  orientation="horizontal")
    ticks=[]
    ticklabels=[]
    #http://shakemap.rm.ingv.it/shake/17769831/download/intensity.jpg
    #http://www.seismo.ethz.ch/static/shakemap/shake/KP201211290852/intensity.html
    theMMIs     = [ 0.06/10,0.06/2, 0.2,     .8,  2,  4.8, 12,   29,    70,  171]
    theMMIcodes = ['',     'I',    'II-III','IV','V','VI','VII','VIII','IX','X']
    if False:
        ticks = [ i/100*9.81 for i in theMMIs]
    else:
        for index,i in enumerate(theMMIs):
            if index==0:
                continue
            bond = numpy.mean(theMMIs[index:index+1])*4 #-1, numpy.exp(numpy.log())
            bond = 9.81*bond/100
            if bond<vmin*gain or bond>vmax*gain:
                continue
            ticks+=[bond]
            ticklabels+=[theMMIcodes[index]]
        if len(ticklabels):
            ticklabels[0] = ticklabels[0]+'\nIntensity'
    cax.set_xscale('log')
    cax.set_xticks(ticks)
    cb.set_ticks(ticks)
    path_effects=[matplotlib.patheffects.withStroke(linewidth=2,
                                                    foreground="w")]
    cax.set_xticklabels(ticklabels,#theMMIcodes,#
                        fontsize='small',
                        horizontalalignment='left',
                        path_effects=path_effects)
    #[l.set_path_effects(path_effects) for l in cax.get_xticklabels()]
    #cax.tick_params(axis='x')
    cax.minorticks_off()
    cax.set_xlim(vmin*gain, vmax*gain)
    
    cbb=cax.twiny()
    cbb.set_xlabel('Ground motion acceleration (m/s$^2$)',fontsize='small')
    cbb.set_xscale('log')
    cbb.set_xlim(vmin*gain, vmax*gain)
    cbb.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    cbb.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    cbb.tick_params(axis='x',labelsize='small')
        
    P=[]
    S=[]
    ttdepth=depth
    if ttdepth<0:
        print("depth is negative")
        ttdepth=0.01
    PSdist=numpy.linspace(0.01,max(eprs),64)
    for epr in PSdist:#sort(eprs):
        arrivals = model.get_travel_times(ttdepth,
                                          distance_in_degree=epr/111,
                                          phase_list=['ttp'],
                                          receiver_depth_in_km=0)
        P+=[modelresidual['p']+numpy.nanmin([ a.time for a in arrivals ])]
        arrivals = model.get_travel_times(ttdepth,
                                          distance_in_degree=epr/111,
                                          phase_list=['tts'],
                                          receiver_depth_in_km=0)
        S+=[modelresidual['s']+numpy.nanmin([ a.time for a in arrivals ])]

    path_effects=[matplotlib.patheffects.withStroke(alpha=.5,
                                                    linewidth=3,
                                                    foreground="w")]
    opts={'color':'k','alpha':.5,'zorder':-9,
          'path_effects':path_effects}
    ax.plot(P,PSdist,':',label='P-waves',**opts)#numpy.sort(eprs)
    ax.plot(S,PSdist,'--',label='S-waves',**opts)#numpy.sort(eprs)
    ax.set_ylim(bottom=0)
    
    ax.set_ylim(ax.get_ylim())
    ax.set_xlim(ax.get_xlim())
    ealiest={}
    for m in event.magnitudes:
        mdt = m.creation_info.creation_time - event.preferred_origin().time
        if m.magnitude_type not in mtypes:
            continue
        if mdt<0 or mdt>max(ax.get_xlim()) :
            continue
        ax.axvline(mdt,
                   color=mtypes[m.magnitude_type],
                   alpha=0.1,
                   linewidth=0.5,
                   zorder=-2)
        label='M$_{%s}$'%m.magnitude_type[1:]
        if label not in ealiest:
            pass
        if label in ealiest and ealiest[label][0]<mdt:
            continue
        ealiest[label]=[mdt,mtypes[m.magnitude_type]]
    for label in ealiest:
        ax.fill_betweenx(ax.get_ylim(),
                         [ealiest[label][0],ealiest[label][0]],
                         x2=max(ax.get_xlim()),
                         label=label,
                         #color=ealiest[label][1],
                         edgecolor=ealiest[label][1],
                         facecolor='None',
                         linewidth=1,
                         #alpha=0.1,
                         zorder=-3)
            
    ax.invert_yaxis()
    ax.tick_params(axis='both',labelsize='small')
    axx=ax.twinx()
    axx.tick_params(axis='both',labelsize='small')
    axx.invert_yaxis()
    axx.set_ylabel('Hypocentral distance (km)',fontsize='small')
    axx.set_ylim(ax.get_ylim())
    fmt = lambda x, pos: "%1g"%round((x**2+(event.preferred_origin().depth/1000)**2)**.5,1)
    axx.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))
    
    #ax.set_yscale('log')
    ax.set_ylabel('Epicentral distance (km)',fontsize='small')
    ax.grid(zorder=-9)
    ax.set_xlabel('Time after origin (seconds)',fontsize='small')
    ax.set_xlim(left=0)
    if False:
        leg=ax.legend(ncol=2,
                      loc="lower left",#"lower right",
                      bbox_to_anchor=(0, 1, 1, 0),#(0, 1/2, 1, 0)
                      fontsize='small')
    else:
        #ax.plot(P,numpy.sort(eprs),':',label='P-waves',**opts)
        #ax.plot(S,numpy.sort(eprs),'--',label='S-waves',**opts)
        ealiest['P-wave'] = [max(P), '0.5', 0]
        ealiest['S-wave'] = [max(S), '0.5', numpy.diff(axx.set_ylim())*0.02]
        for label in ealiest:
            label_path_effects=[matplotlib.patheffects.withStroke(linewidth=5,
                                                #alpha=0.8,
                                                foreground=ealiest[label][1]),
              matplotlib.patheffects.withStroke(linewidth=3,
                                                #alpha=0.8,
                                                foreground="w")]
            y=max(axx.set_ylim())
            va='bottom'
            if 'wave'in label:
                pass#y=min(axx.set_ylim())#-ealiest[label][2]
                #va='center'
            axx.text(ealiest[label][0],
                    y,
                    label,
                    #ha='center',
                    va=va,
                    size='xx-small',
                    fontweight='semibold',
                    path_effects=label_path_effects,
                    zorder=9,
                    #color=ealiest[label][1],
                    #bbox={'fc': '0.8', 'pad': 2},
                    #rotation=90
                    )
                  
    for a in [ax,axx]:
        for axis in [a.xaxis,a.yaxis]:
            axis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    
    if hasattr(event,'event_descriptions'):
        desc = ','.join([d.text for d in event.event_descriptions])
    else:
        geocode = geocoder.osm([event.preferred_origin().latitude,
                                event.preferred_origin().longitude],
                                method='reverse').json
        desc = ','.join(geocode['raw']['display_name'].split(',')[:-2])
    event_type=event.event_type
    if event_type is None or event_type in ['None']:
        event_type='event'
    title='M$_{%s}$%.1f %s (%s)\n'
    if event.preferred_magnitude() is not None:
        title = title % (event.preferred_magnitude().magnitude_type[0:],
                        event.preferred_magnitude().mag,
                        event_type,
                        desc)
    else:
        title = title % ('',
                        -9.9,
                        event_type,
                        desc)

                     
    descs = ', '.join([desc.text for desc in event.event_descriptions if 'region' in desc.type])
    t = '%s\nM$_{%s}$%.1f, %s, %.1fkm deep'
    if event.preferred_magnitude() is not None:
        tt=(str(event.preferred_origin().time)[:19],
            event.preferred_magnitude().magnitude_type[1:],
            event.preferred_magnitude().mag,
            descs,
            event.preferred_origin().depth / 1000.)
    else:
        tt=(str(event.preferred_origin().time)[:19],
            '',
            -9.9,
            descs,
            event.preferred_origin().depth / 1000.)

    ax.set_title(t % tt,
                 loc='left',
                 fontsize='small'
                 )
    return ax.figure