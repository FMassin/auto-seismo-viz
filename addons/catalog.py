#!/usr/bin/env python
from obspy.core.event.catalog import Catalog
from obspy.taup import TauPyModel
from obspy.imaging.beachball import beach
import addons.core
import numpy
import matplotlib.pyplot,matplotlib.figure

def distfilter(self=Catalog(),
               dmax=None,
               dmin=0.,
               x1=None,
               y1=None,
               x2=None,
               y2=None,
               z1=None,
               out=False):
    """
        Filter the events in self (obspy:Catalog) with specified distance range form specified point or line.

    """

    distances=list()
    bad=list()
    good=list()
    for e in self.events:
        d=numpy.nan
        for o in e.origins: # IS THIS ALWAY THE RIGHT ORDER ? SHOULD IT BE SORTED BY CREATION TIME?
            if str(e.preferred_origin_id) == str(o.resource_id):
                break
        if (x1 and y1) or (x2 and y2):
            d = 110.*addons.core.DistancePointLine(o.longitude, o.latitude, x1, y1, x2, y2)
        elif o.quality.minimum_distance :
            d = o.quality.minimum_distance*110.
        else:
            print('warning: no minimum_distance in origin, excluded')
            print(o)
            d=999999999999.
        if z1 is not None:
            d = numpy.sqrt(d**2. + (o.depth/1000-z1)**2)
        distances.append(d)

    for i,d in enumerate(distances):
        for o in self.events[i].origins: # IS THIS ALWAY THE RIGHT ORDER ? SHOULD IT BE SORTED BY CREATION TIME?
            if str(self.events[i].preferred_origin_id) == str(o.resource_id):
                break
        if str(self.events[i].event_type) == 'not existing' or str(o.evaluation_mode)=='automatic':
            pass
        else:
            if d<dmin or d>dmax:
                bad.append(self.events[i])
            else :
                good.append(self.events[i])


    if out in ['d','dist','distance']:
        return distances
    else:
        if out:
            return Catalog(events=good), Catalog(events=bad) #inside, outside
        else:
            return Catalog(events=good)

def get(self, lst, 
        att=None, 
        types=[] , 
        full=False, 
        pref=False, 
        last=False, 
        first=False, 
        nan=True, 
        ct=False, 
        subatt=None):
    """
        1 use
        types = [ 'bests' [, 'all'] [, 'lasts' ] [, 'firsts' ] [, 'nans'] ]

        instead of
        full=False, pref=False, last=False, first=False, nan=True

        2 use
        atts = [ level0, level1, etc ]
        instead
        lst = level0 , att=level1


    """
    out = list()

    for e in self.events:
        patt = None

        if hasattr(e,'preferred_'+lst[:-1]):
            try :
                e['preferred_'+lst[:-1]][att]
                patt = 'preferred_'+lst[:-1]
            except:
                pass
        elif hasattr(e,'preferred_'+lst):
            try :
                e['preferred_'+lst][att]
                patt = 'preferred_'+lst
            except:
                pass

        if hasattr(e,lst) and ((set(['a','all']) & set(types)) or (full and not last and not first and not pref)):
            for o in e[lst]:
                if hasattr(o,att):
                    out.append(o[att])
                elif nan:
                    out.append(numpy.nan)

        elif patt and ((set(['p','pref']) & set(types)) or (pref or (not last and not first))):
            out.append(e[patt][att])

        elif hasattr(e,lst) and len(e[lst]) and ((set(['f','first']) & set(types)) or (first and not pref)):
            val = e[lst][0][att]
            mem=9999999999999999999999999999
            for elt in e[lst]:
                if hasattr(elt,att) and hasattr(elt,'creation_info'):
                    if hasattr(elt.creation_info,'creation_time'):
                        if elt.creation_info.creation_time:
                            if elt.creation_info.creation_time < mem:
                                mem = elt.creation_info.creation_time
                                val = elt[att]
                                if subatt and hasattr(elt[att],subatt) :
                                    val = elt[att][subatt]
            if ct:
                out.append(mem)
            else:
                out.append(val)

        elif hasattr(e,lst) and len(e[lst]) and ((set(['l','last']) & set(types)) or (last or not pref)):
            val = e[lst][-1][att]
            mem=0
            for elt in e[lst]:
                if hasattr(elt,att) and hasattr(elt,'creation_info'):
                    if hasattr(elt.creation_info,'creation_time'):
                        if elt.creation_info.creation_time :
                            if elt.creation_info.creation_time > mem:
                                mem = elt.creation_info.creation_time
                                val = elt[att]
                                if subatt and hasattr(elt[att],subatt) :
                                    val = elt[att][subatt]
            if ct:
                out.append(mem)
            else:
                out.append(val)
        else:
            if nan or set(['n','nan']) & set(types) : #:
                out.append(numpy.nan)

    return out

def map_events(self=Catalog(),
               bmap=None,
               fig=None,
               titletext='',
               eqcolorfield = 'depth',
               colorbar=None,
               fontsize=8,
               cmap = None,
               vmax=40.,
               vmin=0.,
               prospective_inventory=None,
               latencies=None,
               flat_delay=0.,
               VpVsRatio=1.75,
               model_correction='iasp91',
               Vs=3.428,
               extra=None,
               extramarker='1',
               extraname='',
               fp=False,
               magnitude_types=['MVS','Mfd'],
               force=None):
    """
        Map event locations on given basemap and colorcoded with depth or times or delays.

    """


    cf=[]
    mags = get(self,'magnitudes','mag',['b'] )
    times = get(self, 'origins','time', types=['b'])
    stndrdth = ['','st','nd','rd','th']
    for i in range(10):
        stndrdth.append('th')

    if prospective_inventory:

        if eqcolorfield[0] == 'P':
            model = TauPyModel(model=model_correction)
        stations_longitudes, stations_latitudes = addons.core.search(prospective_inventory,
                                                                  fields=['longitude','latitude'],
                                                                  levels=['networks','stations'])

    if len(mags) >0:
        longitudes = [e.preferred_origin().longitude for e in self]
        latitudes = [e.preferred_origin().latitude for e in self]
        #longitudes = get(self, 'origins','longitude', types=['b'])
        #latitudes = get(self, 'origins','latitude', types=['b'])


        if eqcolorfield == 'depth':
            eqcolor = get(self, 'origins',eqcolorfield, types=['b'])
            eqcolor = numpy.asarray(eqcolor)/-1000.
            eqcolorlabel = 'Event depth (km) and magnitude'

        elif eqcolorfield[0] in ['t', 'd', 'e', 'P'] and (eqcolorfield[-1] in [str(int(n)) for n in range(9)] or eqcolorfield[1] in ['M', 'O']):
            eqcolorlabel = 'Travel time to %s$^{th}$ station \n(s'%(eqcolorfield[1:])
            if eqcolorfield[1] in ['M', 'O']:
                if eqcolorfield[-1] == 'l':
                    eqcolorlabel = 'last %s \n(s after Ot)'%(eqcolorfield[-2].upper())
                else :
                    eqcolorlabel = '%s$^{%s}$ %s \n(s after Ot)'%(eqcolorfield[-1],stndrdth[int(eqcolorfield[-1])],eqcolorfield[1:-1])


            if eqcolorfield[0] in ['d']:
                if vmin is None:
                    vmax=130
                eqcolorlabel = 'Error at '+eqcolorlabel
                eqcolorlabel = eqcolorlabel.replace('s after Ot','km from O$_{pref.}$')
            elif eqcolorfield[0] in ['e']:
                if vmax is None:
                    vmax=0.5
                if vmin is None:
                    vmin=-1.0
                eqcolorlabel = 'Error at '+eqcolorlabel
                eqcolorlabel = eqcolorlabel.replace('s after Ot','M unit from M$_{pref.}$')
            elif  latencies is None and ( 'P' in eqcolorfield ):
                eqcolorlabel = 'Travel time for '+eqcolorlabel
            else:
                eqcolorlabel = 'Time at '+eqcolorlabel

            eqcolorlabel = eqcolorlabel.replace(' O (',' origin (')
            eqcolorlabel = eqcolorlabel.replace(' M (',' magnitude (')
            eqcolorlabel = eqcolorlabel.replace(' P (',' P trigger (')
            eqcolorlabel = eqcolorlabel.replace(' OP (',' P trigger (')
            eqcolorlabel = eqcolorlabel.replace(' Op (',' P trigger (')

            ntht=list()
            if prospective_inventory:
                print('Using travel time modeling')
            else:
                print('Using earthquake location travel time')

            for ie,e in enumerate(self.events):
                t=[]
                metric=[]
                origin=e.preferred_origin_id.get_referred_object()
                try:
                    magnitude=e.preferred_magnitude_id.get_referred_object()
                except:
                    print(e)
                    print('no pref mag')
                    magnitude=e.magnitudes[-1]

                if eqcolorfield[1] in ['O']:
                    o=e.origins[0]
                    for o in e.origins:
                        if eqcolorfield[2] not in ['p','P','s','S']:
                            if eqcolorfield[2:-1] in ['',str(o.method_id),str(o.creation_info.author)] : #or eqcolorfield[1] in ['M']:
                                if o.creation_info.creation_time-origin.time>0:
                                    t.append(o.creation_info.creation_time-origin.time)

                        else:
                            if eqcolorfield[-1] == 'l':
                                t.append(numpy.nan)
                            for a in o.arrivals:
                                #if a.phase[0] in [ eqcolorfield[2].upper() , eqcolorfield[2].lower() ]:
                                p = a.pick_id.get_referred_object()
                                if p.time - origin.time>0 and p.creation_info.creation_time - origin.time>0:
                                    if eqcolorfield[2] in ['s','p']:
                                        if eqcolorfield[-1] == 'l' :#and t[-1]<9999999999999999.:
                                            t[-1] = numpy.nanmax([p.creation_info.creation_time - origin.time, t[-1]])
                                        else:
                                            t.append(p.creation_info.creation_time - origin.time)
                                    else:
                                        if eqcolorfield[-1] == 'l' :#and t[-1]<9999999999999999.:
                                            t[-1] = numpy.nanmax([p.time - origin.time, t[-1]])
                                        else:
                                            t.append(p.time - origin.time)

                                        mlatency=0
                                        if latencies:
                                            mlatency = numpy.median(latencies[p.waveform_id.network_code+'.'+p.waveform_id.station_code])
                                        t[-1] += mlatency+flat_delay

                                #if eqcolorfield[-1] =='l':
                                #    pass#    break
                                #elif len([tmp for tmp in t if tmp<vmax]) > int(eqcolorfield[-1])-1:
                                #    break
                        #if eqcolorfield[-1] == 'l':
                        #    print(t)
                elif eqcolorfield[1] in ['M']:
                    for m in e.magnitudes:
                        if m.magnitude_type == 'M':
                            continue
                        #test = re.sub('.*rigin#','', str(m.resource_id ))
                        #test = re.sub('#.*','', str(test))
                        if eqcolorfield[1:-1] in str(m.magnitude_type) or str(m.magnitude_type) in eqcolorfield[1:-1] : # and str( test )  in  str( o.resource_id ) :
                            if eqcolorfield[0] in ['d']:
                                o = m.origin_id.get_referred_object()
                                ep_d = addons.core.haversine(origin.longitude,
                                                              origin.latitude,
                                                              o.longitude,
                                                              o.latitude)
                                hyp_d =  numpy.sqrt((ep_d)**2+(origin.depth - o.depth)**2)
                                metric.append(hyp_d/1000.)
                            elif eqcolorfield[0] in ['e']:
                                try:
                                    metric.append(m.mag-magnitude.mag)
                                except:
                                    metric.append(99999)
                            if m.creation_info.creation_time is None:
                                t.append(99999)
                            elif m.creation_info.creation_time-origin.time>0:
                                t.append(m.creation_info.creation_time-origin.time)
                                #break

                elif eqcolorfield[0] in ['p','P'] :
                    d=[9999999999999999 for d in range(100)]
                    if not prospective_inventory:
                        for a in origin.arrivals:
                            if a.phase[0] in [ eqcolorfield[0].upper() , eqcolorfield[0].lower() ]:
                                p = a.pick_id.get_referred_object()
                                if eqcolorfield[0] in ['s','p']:
                                    t.append(p.creation_info.creation_time - origin.time)
                                else:
                                    t.append(p.time - origin.time)
                                    if latencies:
                                        mlatency = numpy.median(latencies[p.waveform_id.network_code+'.'+p.waveform_id.station_code])
                                        t[-1] += mlatency+flat_delay

                            if eqcolorfield[-1] =='l':
                                pass#    break
                            elif len([tmp for tmp in t if tmp<vmax])>int(eqcolorfield[-1])-1:
                                break


                    else:
                        for istation,lonstation in enumerate(stations_longitudes):
                            latstation =  stations_latitudes[istation]
                            ep_d =  numpy.sqrt((lonstation-origin.longitude)**2+(latstation-origin.latitude)**2)
                            d.append(ep_d)
                        d = numpy.sort(d)
                        for ep_d in  d[:int(eqcolorfield[1:])]:
                            arrivals = model.get_travel_times(origin.depth/1000.,
                                                              distance_in_degree=ep_d,
                                                              phase_list=['ttp'],
                                                              receiver_depth_in_km=0.0)
                            try:
                                t.append( numpy.nanmin([ a.time for a in arrivals ]))
                            except:
                                print('No phase for depth',origin.depth/1000.,'and distance',ep_d)
                                pass
                for tmp in range(100):
                    t.append(99999999999999.)
                if eqcolorfield[0] in ['d','e']:
                    for tmp in range(100):
                        metric.append(99999999999999.)

                if eqcolorfield[0] in ['d','e']:
                    t = [metric[i] for i in numpy.argsort(t)]
                else:
                    t = numpy.sort(t)
                if eqcolorfield[-1] =='l':
                    ntht.append( numpy.nanmin(t) )
                else:
                    ntht.append( t[int(eqcolorfield[-1])-1] )


            ntht = numpy.asarray(ntht)
            if eqcolorfield[0] not in ['e']:
                ntht[ntht<0.]=0.
            eqcolor = ntht

        elif eqcolorfield[0] == 'I' and eqcolorfield[-1] in [str(int(n)) for n in range(9)]:

            if eqcolorfield[1] in ['o','O','m','M','p', 'P']:
                toadd=eqcolorfield[1:-1]
                if toadd in ['P']:
                    toadd='trigger'
                eqcolorlabel = 'MMI at %s$^{%s}$ %s '%(eqcolorfield[-1],stndrdth[int(eqcolorfield[-1])],toadd)
            else:
                eqcolorlabel = 'MMI at %s$^{%s}$ travel time'%(eqcolorfield[-1],stndrdth[int(eqcolorfield[-1])])
            if latencies:
                eqcolorlabel += ' (delays incl.)'
            if eqcolorfield[1]=='0':
                eqcolorlabel = 'Epicentral MMI'

            eqcolorlabel = eqcolorlabel.replace(' O ',' origin ')
            eqcolorlabel = eqcolorlabel.replace(' M ',' magnitude ')
            eqcolorlabel = eqcolorlabel.replace(' P ',' P trigger ')
            eqcolorlabel = eqcolorlabel.replace(' p ',' P trigger ')


            r=list()
            error=list()
            for ie,e in enumerate(self.events):
                d=[9999999999999999. for d in range(100)]
                preferred_origin = e.preferred_origin_id.get_referred_object()
                if eqcolorfield[1] in ['M', 'm']:
                    for m in e.magnitudes:
                        if eqcolorfield[1:-1] in str(m.magnitude_type) or str(m.magnitude_type) in eqcolorfield[1:-1] :
                            if m.creation_info.creation_time-preferred_origin.time <0 :
                                print('WARNING magnitude at', preferred_origin.time, 'created on', m.creation_info.creation_time )
                            else:
                                d.append((m.creation_info.creation_time-preferred_origin.time)*Vs)

                elif eqcolorfield[1] in ['O', 'o']:
                    for o in e.origins:
                        if o.creation_info.creation_time-preferred_origin.time <0 :
                            print('WARNING origin at', preferred_origin.time, 'created on',  o.creation_info.creation_time)
                        else:
                            d.append((o.creation_info.creation_time-preferred_origin.time)*Vs)

                else:
                    if not prospective_inventory:
                        for a in preferred_origin.arrivals:
                            if 'p' in str(a.phase).lower():
                                p = a.pick_id.get_referred_object()
                                for testpick in e.picks:
                                   if (testpick.waveform_id == p.waveform_id and
                                       testpick.phase_hint == p.phase_hint and
                                       testpick.creation_info.creation_time<p.creation_info.creation_time ):
                                       p=testpick

                                if eqcolorfield[1] in ['p', 'P']:
                                    if max([p.creation_info.creation_time, p.time]) < preferred_origin.time:
                                        print('WARNING', 'P at',p.time, ' created on', p.creation_info.creation_time, 'while origin at', preferred_origin.time)
                                    else:
                                        d.append((max([p.creation_info.creation_time, p.time])-preferred_origin.time)*Vs)
                                else:
                                    testd = numpy.sqrt((a.distance*110.)**2+(preferred_origin.depth/1000.)**2)
                                    if testd <= preferred_origin.depth/1000. :
                                        print('WARNING','P at ',testd,'km while origin at ',preferred_origin.depth/1000. ,'km deep' )
                                    d.append( testd/VpVsRatio )

                                if latencies:
                                    mlatency = numpy.median(latencies[p.waveform_id.network_code+'.'+p.waveform_id.station_code])
                                    d[-1] += (mlatency+flat_delay)*Vs
                    else:
                        for istation,lonstation in enumerate(stations_longitudes):
                            latstation =  stations_latitudes[istation]
                            ep_d = addons.core.haversine(lonstation,
                                                        latstation,
                                                        o.longitude,
                                                        o.latitude)
                            d.append( numpy.sqrt((ep_d/1000.)**2+(o.depth/1000.)**2)/VpVsRatio)

                d = numpy.sort(d)
                if eqcolorfield[1]=='0':
                    r.append( preferred_origin.depth/1000. )
                else:
                    r.append( d[int(eqcolorfield[-1])-1] )
            r = numpy.asarray(r)
            r[r<1.]=1.
            eqcolor = obspy_addons.ipe_allen2012_hyp(r,
                                                     numpy.asarray(mags))
            #eqcolor[numpy.isnan(eqcolor)] = numpy.nanmin(eqcolor[eqcolor>.5])
            #eqcolor[eqcolor<1.] =  numpy.nan #numpy.nanmin(eqcolor[eqcolor>.5])
            cmap='nipy_spectral'
            vmin=1.
            vmax=8.5
        mags = []
        for e in self:
            try:
                mags+=[e.preferred_magnitude().mag ]
            except:
                mags+=[e.magnitudes[-1].mag ]
        sizes, sizesscale, labelsscale, x = eventsize( mags = mags,force=force)

        longitudes,latitudes=bmap(longitudes,latitudes)
        bmap.scatter(longitudes,
                         latitudes,
                         sizes,
                         edgecolor='w',
                         lw=1.5,
                        zorder=992)
        bmap.scatter(longitudes,
                         latitudes,
                         sizes,
                         edgecolor='k',
                         lw=1,
                        zorder=992)
        bmap.scatter(longitudes,
                     latitudes,
                     sizes,
                     color='w',
                     edgecolor='none',
                     lw=0.,
                     zorder=992)

        sortorder =  numpy.argsort(sizes)#eqcolor)
        if eqcolorfield[0] in [ 't' , 'P' , 'p' ]:
            sortorder =  sortorder[::-1]

        norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
        if eqcolorfield[0] not in 'd':
            eqcolor[eqcolor>vmax]=numpy.nan
        else:
            pass#norm = matplotlib.colors.LogNorm(vmin=1,vmax=vmax)
        #eqcolor[eqcolor<vmin]=numpy.nan
        cf = bmap.scatter([longitudes[i] for i in sortorder],
                          [latitudes[i] for i in sortorder],
                          [sizes[i] for i in sortorder],
                          [eqcolor[i] for i in sortorder ],
                          edgecolor='None',
                          cmap=cmap,
                          norm = norm,
                          zorder=995)


        nextra=''
        if extra:
            nextra=0
            for i in numpy.argsort(eqcolor):
                test = [((e[10]-longitudes[i])**2+(e[9]-latitudes[i])**2+((obspy.UTCDateTime(e[0],e[1],e[2],e[3],e[4],e[5])-times[i])/1200)**2)**.5 for e in extra]


                if numpy.nanmin(test)<.1 :

                    t = extra[numpy.argmin(test)]
                    nextra+=1
                    bmap.scatter(x=longitudes[i],
                                  y=latitudes[i],
                                  s=sizes[i]*1.1,
                                  c='w',
                                  marker=extramarker,
                                  latlon=True,
                                  zorder=997)
                    bmap.ax.text(bmap(longitudes[i], latitudes[i])[0]+t[12][0]/2,
                                  bmap(longitudes[i], latitudes[i])[1]+t[12][1],
                                  s='%s'%(t[-1]),#,t[1],t[2]),
                                  va='top',
                                  ha=t[11],
                                  alpha=.5,
                                  fontsize='x-small',
                                  path_effects=[
                                      matplotlib.patheffects.withStroke(linewidth=3,
                                                       foreground="white")],
                                 zorder=1000)
            extraname='and %s %s ' % (str(nextra), extraname)

        if fp:
            cmap=cf.get_cmap()

            #client = Client("USGS")
            for i in numpy.argsort(sizes)[max([-100,-1*len(sizes)]):]:
                if sizes[i]>=fp:
                    id = str(self.events[i].resource_id).split('id=')[-1].split('&')[0]
                    #try:
                    e = self.events[i] #client.get_events(eventid=id,producttype='moment-tensor').events[0]
                    o = e.preferred_origin()
                    for f in e.focal_mechanisms[::-1]:
                        t=f.moment_tensor.tensor
                        xy=bmap(longitudes[i],latitudes[i])
                        try:
                            b = beach([t.m_rr,t.m_tt,t.m_pp,t.m_rt,t.m_rp,t.m_tp],
                                      xy=(xy[0],xy[1]),
                                      width=.028*sizes[i]**.525,
                                      linewidth=0,
                                      alpha=1,
                                      facecolor=cmap(norm(eqcolor[i])),
                                      #edgecolor=cmap(norm(eqcolor[i]))
                                      )
                            b.set_zorder(999)
                            bmap.ax.add_collection(b)
                            break
                        except:
                            pass
                        #except:
                        #pass

        titletext += '\n%s events %s(%s' % (len(times), extraname, str(min(times))[:10])
        if str(max(times))[:10] > str(min(times))[:10]:
            titletext += ' to %s)' % str(max(times))[:10]
        else:
            titletext += '%s)' % titletext

        if colorbar and len(self.events)>0 :
            fig.cb = obspy_addons.nicecolorbar(cf,
                                  axcb = bmap.ax,
                                  label = eqcolorlabel,
                                  data = eqcolor,
                                  fontsize=fontsize)
            fig.cb.ax.scatter(numpy.repeat(.1,len(sizesscale)),
                              x*(vmax-vmin)+vmin,
                              s=sizesscale,
                              facecolor='None',
                              linewidth=3,
                              edgecolor='w')
            fig.cb.ax.scatter(numpy.repeat(.1,len(sizesscale)),
                              x*(vmax-vmin)+vmin,
                              s=sizesscale,
                              facecolor='None',
                              edgecolor='k')
            for i,v in enumerate(labelsscale):
                fig.cb.ax.text(-.2,x[i]*(vmax-vmin)+vmin,'M'+str(v), va='center',ha='right', rotation='vertical',fontsize=fontsize)
            if eqcolorfield == 'depth':
                fig.cb.ax.set_yticklabels([ str(int(l)).split('−')[-1] for l in fig.cb.ax.get_yticks()])

            if eqcolorfield[0] == 'I' and eqcolorfield[-1] in [str(int(n)) for n in range(100)]:
                fig.cb.ax.set_yticklabels(obspy_addons.num2roman(list(range(1,15))))

    return titletext

def performance_timelines(event,
                          stations=None,
                          vs=3.7,
                          magnitude_types=['Mfd','MVS'],
                          authors=['scvsmag2','scvsmag'],
                          lineauthors=None, #'scfinder'],
                          type='soil',
                          #xlim=59,
                          lettering_offset=0):

    for author in authors:
            if author in [m.creation_info.author.split('@')[0] for m in event.magnitudes if m.creation_info.author is not None]:
                print("Old author list:",authors)
                authors= [author]
                print("New author list:",authors)
                break
            
    authors += addons.core.scfinderauthor(event.preferred_origin(),
                                          lineauthors=lineauthors)
    w,h=matplotlib.figure.figaspect(1.)
    f=matplotlib.pyplot.figure(figsize=(w*1.4,h*1.3))
    initaxes = f.subplots(4,1,sharex=True)
    axes = [initaxes[2], initaxes[1], initaxes[0], initaxes[3]]
    
    axesdata=[{'x':{},
               'y':{},
               'yep':{},
               'yem':{},
               'label':{},
               'title':r'M$_{type}$',
               'ylabel':'Magnitude'},
              {'x':{},
               'y':{},
               'yep':{},
               'yem':{},
               'label':{},
               'title':r'M$_{type}$ (last loc.)',
               'ylabel':r'Location $\delta$ (km)'},
              {'x':{},
               'y':{},
               'yep':{},
               'yem':{},
               'label':{},
               'title':r'M$_{type}$',# (last #)',
               'ylabel':r'Stations #'},
              {'x':{},
               'y':{},
               'yep':{},
               'yem':{},
               'label':{},
               'title':r'M$_{type}$ (intensity)',
               'ylabel':'Max. pred.\nintensity'}]



    preferred_magnitude=event.preferred_magnitude()
    preferred_origin=event.preferred_origin()
    mPG={}
    magnitudes = []
    for t in magnitude_types:
        for m in event.magnitudes[::-1]:
            if m.magnitude_type==t and m.creation_info.author.split('@')[0] in authors:
                    magnitudes +=[m]
    for m in magnitudes:
        if (m.magnitude_type in magnitude_types and
            m.creation_info.author.split('@')[0] in authors):
            
            
            for itmp,tmp in enumerate([preferred_magnitude, m]):
                key = (tmp.magnitude_type,tmp.creation_info.author)
                o=tmp.origin_id.get_referred_object()
                tmpeewdelay = [m.creation_info.creation_time-preferred_origin.time]
                
                if key not in axesdata[0]['x'].keys() and tmp.magnitude_type not in magnitude_types:
                    tmpeewdelay = list(numpy.linspace(int(preferred_origin.depth/1000/vs),tmpeewdelay[0],64))
                    
                for eewdelay in tmpeewdelay:
                    if key not in axesdata[0]['x'].keys():
                        axesdata[0]['x'][key]=[]
                        axesdata[0]['y'][key]=[]
                        axesdata[0]['yep'][key]=[]
                        axesdata[0]['yem'][key]=[]
                        axesdata[1]['x'][key]=[]
                        axesdata[1]['y'][key]=[]
                        axesdata[1]['yep'][key]=[]
                        axesdata[1]['yem'][key]=[]
                        axesdata[3]['x'][key]=[]
                        axesdata[3]['y'][key]=[]
                        axesdata[3]['yep'][key]=[]
                        axesdata[3]['yem'][key]=[]
                        mPG[key]=0

                    #############################################
                    axesdata[0]['x'][key].append(eewdelay)
                    axesdata[0]['y'][key].append(tmp.mag)
                    mag_errors_uncertainty=tmp.mag_errors['uncertainty'] or \
                                           numpy.percentile(numpy.abs([sm.residual for sm in tmp.station_magnitude_contributions]+[0]),34) or \
                                           0
                    axesdata[0]['yep'][key].append(mag_errors_uncertainty)
                    axesdata[0]['yem'][key].append(mag_errors_uncertainty)

                    if eewdelay==min(axesdata[0]['x'][key]) or (eewdelay==0 and tmp.magnitude_type not in magnitude_types):
                        label = r'M$_{%s}$ %.1f'%(tmp.magnitude_type[1:],
                                                    tmp.mag)
                        axesdata[0]['label'][key]=label

                    #############################################
                    d = addons.core.haversine(preferred_origin.longitude,
                                                   preferred_origin.latitude,
                                                   o.longitude,
                                                   o.latitude)
                    d = (d**2+(preferred_origin.depth-o.depth)**2)**.5
                    axesdata[1]['y'][key].append(d/1000)

                    de=((o.latitude_errors['uncertainty']or 5)**2+
                        (o.longitude_errors['uncertainty']or 5)**2+
                        ((o.depth_errors['uncertainty']or 5000)/1000)**2)**.5
                    #de+=110*((preferred_origin.latitude_errors['uncertainty']or 0)**2+
                    #        (preferred_origin.longitude_errors['uncertainty']or 0)**2+
                    #        ((preferred_origin.depth_errors['uncertainty']or 0)/110000)**2)**.5
                    axesdata[1]['yep'][key].append(de)
                    axesdata[1]['yem'][key].append(de)
                    axesdata[1]['x'][key].append(eewdelay)
                    if eewdelay==min(axesdata[1]['x'][key]) or  (eewdelay==0 and tmp.magnitude_type not in magnitude_types):
                        label = r'%s$_{%s}$ $^{N%.2f°}_{E%.2f°}$%.0f$^{bsl}_{km}$'%(tmp.magnitude_type[0],
                                                             tmp.magnitude_type[1:],
                                                             o.latitude,
                                                             o.longitude,
                                                             o.depth/1000. )
                        axesdata[1]['label'][key]=label

                    #############################################
                    memlogPGcm=-99999
                    if stations is None:
                        eewdistance=eewdelay*vs-d/1000
                        memlogPGcm=9
                    else:
                        for network in stations:
                            for station in network:
                                distance=addons.core.haversine(o.longitude,
                                                       o.latitude,
                                                       station.longitude,
                                                       station.latitude)

                                distance = ((distance**2+(-station.elevation-o.depth)**2)**.5)/1000.
                                            
                                logPGcm = addons.core.ipe_allen2012_hyp(numpy.asarray([max([.1,abs(eewdistance)])]),
                                                                        numpy.asarray([tmp.mag]),
                                                                        s=-1)[0]
                                          
                                if distance >= eewdelay*vs and logPGcm>memlogPGcm:
                                    memlogPGcm=logPGcm
                                    eewdistance=distance

                    if memlogPGcm>-99999:
                        axesdata[3]['x'][key].append(eewdelay)
                        logPGcm = addons.core.ipe_allen2012_hyp(numpy.asarray([max([.1,abs(eewdistance)])]),
                                                                numpy.asarray([tmp.mag]),
                                                                s=-1)
                        if numpy.isnan(logPGcm[0]):
                            print('WARNING: NaN PGA at %g km'%(eewdistance))

                        PG=logPGcm[0]
                        axesdata[3]['y'][key].append(PG)
                        if PG>mPG[key] or eewdelay==0:
                            mPG[key]=PG
                            if eewdelay==min(axesdata[3]['x'][key]) or (eewdelay==0 and tmp.magnitude_type not in magnitude_types):
                                label = r'%s$_{%s}$ %.0e m/s$^2$'%(tmp.magnitude_type[0],
                                                               tmp.magnitude_type[1:],
                                                               PG)
                                axesdata[3]['label'][key]=label

                        if mag_errors_uncertainty<2.:
                            logPGcm = addons.core.ipe_allen2012_hyp(numpy.asarray([max([.1,abs(eewdistance-de)])]),
                                                                            numpy.asarray([tmp.mag+mag_errors_uncertainty]),
                                                                            s=-1)
                            if numpy.isnan(logPGcm[0]):
                                print('WARNING: Nan PGA for %g = %g - %g'%(eewdistance-de,eewdistance,de))
                            axesdata[3]['yep'][key].append(logPGcm[0] - PG)
                            
                            logPGcm = addons.core.ipe_allen2012_hyp(numpy.asarray([max([.1,abs(eewdistance+de)])]),
                                                                    numpy.asarray([tmp.mag-mag_errors_uncertainty]),
                                                                    s=-1)
                            if numpy.isnan(logPGcm[0]):
                                print('WARNING: Nan PGA %g = %g + %g'%(eewdistance-de,eewdistance,de))
                            axesdata[3]['yem'][key].append(PG - logPGcm[0])
                        else:

                            axesdata[3]['yep'][key].append(0)
                            axesdata[3]['yem'][key].append(0)

                    #############################################

                    if not hasattr(o.quality,'used_phase_count'):
                        continue
                    if key not in axesdata[2]['x'].keys():
                        axesdata[2]['x'][key]=[]
                        axesdata[2]['y'][key]=[]
                        axesdata[2]['yep'][key]=[]
                        axesdata[2]['yem'][key]=[]
                        
                    n = o.quality.used_phase_count #tmp.station_count
                    narr=n
                    if len(preferred_origin.arrivals)>0:
                        narr=0
                        for arr in preferred_origin.arrivals:
                            pick = arr.pick_id.get_referred_object()
                            if pick is not None:
                                test=pick.time-preferred_origin.time
                                if ('P' in arr.phase and
                                    test<=eewdelay):#tmp.creation_info.creation_time):
                                    narr+=1
                    else:
                        pass
                        # use best origin available?
                    if itmp:
                        axesdata[2]['y'][key].append(n)
                    else:
                        axesdata[2]['y'][key].append(narr)

                    axesdata[2]['yep'][key].append(0)#max([0,narr-n]))
                    axesdata[2]['yem'][key].append(0)#max([0,n-narr]))
                    axesdata[2]['x'][key].append(eewdelay)
                    if eewdelay==min(axesdata[2]['x'][key]) or (eewdelay==0 and tmp.magnitude_type not in magnitude_types):
                        author=''
                        if tmp.magnitude_type in magnitude_types:
                            author = tmp.creation_info.author.split('@')[0]

                        label = r'M$_{%s}^{%s}$'%(tmp.magnitude_type[1:],author# %.0f (%s.)
                                                       #n,
                                                       #key[0][0]
                                                       )
                        axesdata[2]['label'][key]=label


    maxt=0
    maxy=0


    for i,data in enumerate(axesdata):
        for key in data['y']:
            for e in range(9):
                data['label'][key] = data['label'][key].replace('e-%02.0f'%e,'$e^{-%.0f}$'%e)
                data['label'][key] = data['label'][key].replace('e%02.0f'%e,'$e^{%.0f}$'%e)
            data['label'][key] = data['label'][key].replace('N-','S')
            data['label'][key] = data['label'][key].replace('E-','W')
            indexes=numpy.argsort(data['x'][key])
            color='C0'
            if key[0] in magnitude_types:
                color = 'C%d'%(magnitude_types.index(key[0])+1)
            axes[i].fill_between([data['x'][key][i] for i in indexes],
                                 [data['y'][key][i]-data['yem'][key][i] for i in indexes],
                                 [data['y'][key][i]+data['yep'][key][i] for i in indexes],
                                 step = 'post',
                                 facecolor='w',
                                 alpha=0.7,
                                 zorder=9)
            axes[i].fill_between([data['x'][key][i] for i in indexes],
                                 [data['y'][key][i]-data['yem'][key][i] for i in indexes],
                                 [data['y'][key][i]+data['yep'][key][i] for i in indexes],
                                 step = 'post',
                                 facecolor=color,
                                 zorder=10,
                                 alpha=0.3)
            axes[i].step([data['x'][key][i] for i in indexes],
                         [data['y'][key][i] for i in indexes],
                         where='post',
                         color=color,
                         zorder=11,
                         label=data['label'][key])

            axes[i].set_ylabel(data['ylabel'],fontsize='small')
                      
            maxt=max([maxt,max(data['x'][key])])
            maxy=max([maxy,max(data['y'][key])])
        
        if i==2:
            legend=axes[i].legend(#title=data['title'],
                           prop={'size': 'small'},
                           ncol=3,
                           loc='lower right',
                           bbox_to_anchor=(0, 1, 1, 0)
                           )
            legend.remove()
            axes[i].add_artist(legend)

        axes[i].legend([],[],
        title='abcdefghijklmnopqrst'[list(initaxes).index(axes[i])+lettering_offset],
                  title_fontsize='small',
                  loc='best',
                  #fontsize=0,
                  labelspacing=-.2)


    x=[]
    for arr in preferred_origin.arrivals:
        pick = arr.pick_id.get_referred_object()
        if pick is not None:
            if 'P' in arr.phase:
                x.append(pick.time-preferred_origin.time)
    axes[2].step(numpy.sort(x),
                 range(1,len(x)+1),
                 color='k',
                 alpha=.1,
                 where='post',
                 zorder=-9
                    )
    x=[]
    y=[]
    ylim=axes[3].get_ylim()

    axes[3].set_ylim(bottom=ylim[0])        
    for ax in axes:
        #ax.grid()
        ax.set_facecolor('w')
        ax.tick_params(right=True, top=True,
                       left=True, bottom=True,
                       which='both')
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.grid(b=True, which='major', color='gray', linestyle='dashdot', zorder=-9999)
        ax.grid(b=True, which='minor', color='beige',  ls='-', zorder=-9999)  
        
    for a in [2,1]:
        axes[a].set_yscale('log')
        lim=axes[a].get_ylim()
        axes[a].set_yticks([0.001,0.01,0.1,1,10,100,1000])
        axes[a].set_yticklabels([0.001,0.01,0.1,1,10,100,1000])
        axes[a].set_ylim(lim)
        y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = numpy.arange(1.0, 10.0) * 0.1, numticks = 10)
        axes[a].yaxis.set_minor_locator(y_minor)
        
    lim=axes[2].get_ylim()
    axes[2].set_yticks([4,6,10,100,1000])
    axes[2].set_yticklabels(["$^{_4}$","$^{_6}$",10,100,1000])
    axes[2].yaxis.set_minor_formatter('')
    axes[2].set_ylim([.91,max(lim)])
    
            
    lim=axes[3].get_ylim()
    axes[3].set_yticks([1,2,3,4,5,6,7,8,9,10,11,12])
    axes[3].set_yticklabels(['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII'])
    axes[3].set_ylim(lim)
    top=max([max(axesdata[3]['y'][key][1:]+data['yep'][key][1:]) for key in axesdata[3]['y']])
    bot=numpy.nanmin([numpy.nanmin([axesdata[3]['y'][key][i]-data['yem'][key][i] for i in range(len(axesdata[3]['y'][key]))]) for key in axesdata[3]['y']])
    axes[3].set_ylim([bot-.3,top+.3])
    #axes[3].set_ylim(bottom=1)

    lab='Time after origin (second)'
    ax.set_xlabel(lab,fontsize='small')

    if max(axes[0].get_ylim())>9.5:
        axes[0].set_ylim(top=9.5)
    if min(axes[0].get_ylim())<4:
        pass#axes[0].set_ylim(bottom=4)
    if max(axes[3].get_ylim())>40:
        axes[3].set_ylim(top=40)
    #if min(axes[3].get_ylim())<1:
    #    axes[3].set_ylim(bottom=1)

    xlim = numpy.nanmax([numpy.nanmax(axesdata[0]['x'][key]) for key in axesdata[0]['x']])
    #print([axesdata[0]['x'][key] for key in axesdata[0]['x']])
    #print(xlim)
    axes[3].set_xlim([1,xlim+xlim/9])
    #axes[3].set_xlim(left=1)#,xlim])
    
    descs = ', '.join([desc.text for desc in event.event_descriptions if 'region' in desc.type])
    t = '%s\nM$_{%s}$%.1f, %s, %.1fkm deep'
    tt=(str(event.preferred_origin().time)[:19],
        event.preferred_magnitude().magnitude_type[1:],
        event.preferred_magnitude().mag,
        descs,
        event.preferred_origin().depth / 1000.)
    axes[2].set_title(t % tt,
                      loc='left',
                      fontsize='small')

    f.tight_layout()
    f.subplots_adjust(hspace = 0)
    f.align_ylabels(axes)
    f.tight_layout(pad=0)
    return f