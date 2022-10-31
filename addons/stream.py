#!/usr/bin/env python
from gettext import Catalog
import matplotlib.pyplot,matplotlib.colors,matplotlib.patheffects,numpy
from obspy.taup import TauPyModel
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from obspy.core.stream import Stream, Trace
from obspy.core import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
gold= (1 + 5 ** 0.5) / 2.

def lineMagnitude(x1, y1, x2, y2):
    lineMagnitude = numpy.sqrt(numpy.power((x2 - x1), 2)+ numpy.power((y2 - y1), 2))
    return lineMagnitude

#Calc minimum distance from a point and a line segment (i.e. consecutive vertices in a polyline).
def DistancePointLine(px, py, x1, y1, x2, y2):
    #http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/source.vba
    LineMag = lineMagnitude(x1, y1, x2, y2)

    if LineMag < 0.00000001:
        DistancePointLine = numpy.sqrt(numpy.power((px - x1), 2)+ numpy.power((py - y1), 2)) # 9999
        return DistancePointLine

    u1 = (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u = u1 / (LineMag * LineMag)

    if (u < 0.00001) or (u > 1):
        #// closest point does not fall within the line segment, take the shorter distance
        #// to an endpoint
        ix = lineMagnitude(px, py, x1, y1)
        iy = lineMagnitude(px, py, x2, y2)
        if ix > iy:
            DistancePointLine = iy
        else:
            DistancePointLine = ix
    else:
        # Intersecting point is on the line, use the formula
        ix = x1 + u * (x2 - x1)
        iy = y1 + u * (y2 - y1)
        DistancePointLine = lineMagnitude(px, py, ix, iy)

    return DistancePointLine

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
    epr = gps2dist_azimuth(origin.latitude,origin.longitude,
                           station.latitude,station.longitude)[0]/1000
    if optline is not None:
        epr = DistancePointLine(station.longitude,
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
    #print("model residual:",modelresidual)
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
               combine='both',
               nopgalabel=False,
               nostationlabel=False,
               ):
    if ax is None:
        w, h = matplotlib.figure.figaspect(gold) # gold aspect
        ax = matplotlib.pyplot.figure(figsize=(w*gold, h*gold)).add_axes([0.1, 0.1, 0.7, 0.7]) # bigger by factor of gold

        #size = ax.figure.get_size_inches()
        #sax.figure.set_size_inches([size[0],size[0]*gold])
        #ax.set_position(ax.get_position().shrunk_to_aspect(gold))
        #axx.set_position(axx.get_position().shrunk_to_aspect(gold))

    if hasattr(catalog,'events'):
        event = catalog[0]
    else:
        event = catalog
    axx=ax.twinx()
    #ax.sharey(axx)

    time = event.preferred_origin().time
    depth = event.preferred_origin().depth/1000
    modelresidual,model=get_velmodel_correction(event,inventory,lim,model)
    tmax=0
    vmax=[]
    vmin=[]
    for tr in self:
        if numpy.nanmax(tr.data) > 0 :
            imax=numpy.nanargmax(tr.data)
            tmax=numpy.nanmax([tmax, tr.times(type="utcdatetime")[imax]])
            vmax+=[numpy.nanmax(tr.data)]
    for tr in self:
        if numpy.nanmax(tr.data[:99]) > 0 :
            vmin+=[numpy.nanmax(tr.data[:99])]

    vmin=numpy.nanpercentile(vmin,84)
    vmax=numpy.nanpercentile(vmax,100)*1.1
    norm=matplotlib.colors.LogNorm(vmin=vmin*gain,
                                   vmax=vmax*gain)
    print(vmin,'*',gain,vmax,'*',gain)
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
                if True:
                    axx.annotate(data_label, 
                                  xy=(x, epr), xycoords="data",
                                  ha='left',va='center',
                                  size='xx-small',
                                  fontweight='semibold',
                                  path_effects=label_path_effects,
                                  zorder=5+y[0]
                                  #bbox=dict(boxstyle="round", fc="w")
                                  )
                else:
                    ax.text(x,
                        epr,#r,tr.data[imax],
                        data_label,#),#c[imax]
                        ha='left',
                        va='center',
                        clip_on=False,
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
                    width="%d%%"%((1-1/gold)*100),  
                    height="1%",  # height : 3%
                    loc='upper right',
                    borderpad=-0.25)
    h=matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap)
    cb=matplotlib.pyplot.colorbar(h,
                                  cax=cax,
                                  extend='both',
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
#    cax.minorticks_off()
    cax.set_xlim(vmin*gain, vmax*gain)
    cax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs=numpy.arange(2, 10) * .1, numticks=100))
    
    cbb=cax.twiny()
    cbb.set_xlabel('Ground motion acceleration (m/s$^2$)',fontsize='small')
    cbb.set_xscale('log')
    cbb.set_xlim(vmin*gain, vmax*gain)
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
    axx.tick_params(axis='both',labelsize='small')
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

def select(stream,
           startafter=None,
           endbefore=None,
           startbefore=None,
           endafter=None,
           **kargs):
    out = stream.select(**kargs)
    if startafter is not None:
        out = [tr for tr in out if tr.stats.starttime>=startafter]
    if endbefore is not None:
        out = [tr for tr in out if tr.stats.endtime<=endbefore]
    if startbefore is not None:
        out = [tr for tr in out if tr.stats.starttime<startbefore]
    if endafter is not None:
        out = [tr for tr in out if tr.stats.endtime>endafter]
    return Stream(out)

def overlap(reftrace,
            trace):
    starttime = max([reftrace.stats.starttime , trace.stats.starttime])
    endtime =  min([reftrace.stats.endtime , trace.stats.endtime])
    if starttime>=endtime:
        print('WARNING! No overlaping data in')
        print(starttime,endtime)
        return [],[],starttime,0
    iref = (reftrace.times(reftime=starttime)>=0)&(reftrace.times(reftime=endtime)<=0)
    itr = (trace.times(reftime=starttime)>=0)&(trace.times(reftime=endtime)<=0)
    return iref, itr, starttime, sum(iref)

def combinehoriz(EN,
                 horizontal_code='h'):
    starttime = max([tr.stats.starttime for tr in EN])
    endtime = min([tr.stats.endtime for tr in  EN])
    if endtime<=starttime:
        longer = numpy.argmax([tr.stats.npts for tr in EN])
        return EN[longer:longer],EN[longer]
    EN = EN.slice(starttime,endtime)
    npts = min([tr.stats.npts for tr in EN])
    for tr in EN:
        tr.data=tr.data[:npts]
    iref, itr, starttime, npts = overlap(EN[0],EN[1])
    tr = Trace(header=EN[0].stats) #EN[0].copy()
    tr.stats.channel = tr.stats.channel[:-1]+horizontal_code
    tr.stats.npts = npts
    tr.stats.startime = starttime
    if isinstance(EN[0].data[0] , UTCDateTime):
        tr.data = numpy.max([EN[0].data[iref],
                             EN[1].data[itr]],
                            axis=0)
    else:
        tmp = EN[0].data[iref]**2
        l=min([len(EN[1].data[itr]),len(EN[0].data[iref])])
        tmp[:l-1] += EN[1].data[itr][:l-1]**2
        tr.data = tmp**.5
    return EN, tr

def combinechannels(stream = Stream(),
                    combine = 'all',
                    horizontal_code = 'b',
                    tridimentional_code = 't',
                    difference_code = 'd',
                    ref = Stream(),
                    max_code = 'm',
                    percentile=False,
                    verbose=False,
                    limit=9999999):
    dontdo=[]
    tridim = Stream()
    horiz = Stream()
    diff = Stream()
    for trace in stream:
        if trace.id[:-1] not in dontdo:
            dontdo.append(trace.id[:-1])
            if combine in 'differences':
                if len(diff)>limit:
                    print('REACHED LIMIT!!!')
                    return diff
                #reftrace=ref.select(id=trace.id)
                reftrace = select(ref,
                                  startafter = trace.stats.starttime,
                                  endbefore = trace.stats.endtime,
                                  id=trace.id)
                if (not len(reftrace) and
                    len(trace.stats.channel)==2):
                    E = select(ref,
                               startafter = trace.stats.starttime,
                               endbefore = trace.stats.endtime,
                               id=trace.id+'[E,2]')
                    N = select(ref,
                              startafter = trace.stats.starttime,
                              endbefore = trace.stats.endtime,
                              id=trace.id+'[N,3]')
                    #E = ref.select(id=trace.id+'[E,2]').slice(trace.stats.starttime,trace.stats.endtime)
                    #N = ref.select(id=trace.id+'[N,3]').slice(trace.stats.starttime,trace.stats.endtime)
                    try:
                        EN, reftrace = combinehoriz(Stream([E[0], N[0]]), horizontal_code='')
                    except:
                        if verbose:
                            print('WARNING! No reference streams found in ref for ')
                            print(trace)
                        continue
                else:
                    try:
                        reftrace=reftrace[0]
                    except:
                        if verbose:
                            print('WARNING! No reference stream found in ref for ')
                            print(trace)
                        continue
                if not len(reftrace.data) :
                    if verbose:
                        print('WARNING! No reference trace found in ref for ')
                        print(trace)
                    continue
                        
                iref, itr, starttime, npts = overlap(reftrace,trace)
                tr = Trace(header=trace.stats) #reftrace.copy()
                tr.stats.starttime = starttime
                tr.stats.npts = npts
                tr.stats.channel = tr.stats.channel+difference_code
                
                if isinstance(trace.data[0] , UTCDateTime):
                    tr.data = numpy.asarray([ reftrace.data[iref][i]-d for i,d in enumerate(trace.data[itr])])
                else:
                    tr.data = reftrace.data[iref]-trace.data[itr]
                    if percentile:
                        tr.data[reftrace.data<numpy.sort(reftrace.data)[int(len(reftrace.data)*percentile)]]=numpy.nan
                diff += tr
            else:
                Z = stream.select(id=trace.id[:-1]+'[Z,1]').slice(trace.stats.starttime,trace.stats.endtime)
                E = stream.select(id=trace.id[:-1]+'[E,2]').slice(trace.stats.starttime,trace.stats.endtime)
                N = stream.select(id=trace.id[:-1]+'[N,3]').slice(trace.stats.starttime,trace.stats.endtime)

                if combine in 'horizontal2dimensionaltwodimensionalboth' and E and N:
                    EN,tr = combinehoriz(E+N,
                                         horizontal_code=horizontal_code)
                    horiz += EN
                    horiz += tr

                if Z and E and N and  combine in 'all3dimensionaltridimensionalboth' :
                    starttimes = [tr.stats.starttime for tr in E+N+Z]
                    endtimes = [tr.stats.endtime for tr in E+N+Z]
                    if max(starttimes)>min(endtimes):
                        longer = numpy.argmax([tr.stats.npts for tr in E+N+Z])
                        tridim += (E+N+Z)[longer]
                        tr = (E+N+Z)[longer].copy()
                    else:
                        ZEN = (Z+E+N).slice(max(starttimes),min(endtimes))
                        npts = [tr.stats.npts for tr in ZEN]
                        tridim += ZEN
                        tr = ZEN[0].copy()
                    tr.stats.channel = tr.stats.channel[:-1]+tridimentional_code
                    tr.data = (ZEN[0].data[:min(npts)]**2+ZEN[1].data[:min(npts)]**2+ZEN[2].data[:min(npts)]**2)**.5
                    tridim += tr

    if combine in 'both':
        return tridim , horiz
    elif combine in 'all3dimensionaltridimensional':
        return tridim
    elif combine in 'differences':
        return diff
    return horiz