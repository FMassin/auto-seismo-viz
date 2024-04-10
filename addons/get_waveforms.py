#!/usr/bin/env python
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel
import numpy


def locations2degrees_origin2channel(origin,channel):

    lat1 = origin.latitude
    lat2 = channel.latitude
    lon1 = origin.longitude
    lon2 = channel.longitude

    return locations2degrees(lat1, lon1, lat2, lon2)

def fixrespunit(inventory,debug=False):

    for n,net in enumerate(inventory):
        for s,stat in enumerate(net):
            for c,cha in enumerate(stat):

                if (sensitivity:=cha.response.instrument_sensitivity) is not None:

                    inventory[n][s][c].response.instrument_sensitivity.input_units = sensitivity.input_units.replace('m/s','M/S')

                    for a,stage in enumerate(cha.response.response_stages):
                        if (units:=stage.input_units) is not None :#and "m/s" in stage.input_units:
                            
                            inventory[n][s][c].response.response_stages[a].input_units = units.replace('m/s','M/S')
    
    return inventory #anyways happens in place ???

def attach_distance(self,
                    inventory,
                    origin,
                    debug=False):

    for trace in self.select():
        station = inventory.select(network=trace.stats.network,
                                   station=trace.stats.station)
        if len(station.networks)==0:
            if debug:
                print('Cannot find',trace)
            self.remove(trace)
            continue

        station = station.networks[0].stations[0]
        if debug:
            print('Found %s.%s for %s'%(trace.stats.network,trace.stats.station, str(trace)))
        trace.stats.coordinates = {'latitude': station.latitude,
                                   'longitude':station.longitude,
                                   'elevation':station.elevation}
        if origin.latitude is None or origin.longitude is None:
            trace.stats.distance = 0.
            continue
        distance = gps2dist_azimuth(station.latitude,
                                    station.longitude,
                                    origin.latitude,
                                    origin.longitude)[0]
        distance = ((distance**2+(trace.stats.coordinates['elevation']*-1)**2.)**.5)
        if not hasattr(trace.stats, 'distance'):
            trace.stats.distance = 0.
        trace.stats.distance = distance
    
    return self.sort(keys=['distance'])

def event2stations_bulk(self,
                           origin,
                           endafter=40,
                           model = TauPyModel(),
                           debug=False,
                           location='*',
                           channel='SN*,SH*,EN*,EH*,HN*,HH*,HG*'):
    allthes = []
    for mseedid in self.get_contents()['channels']:
        nsc = self.select(*mseedid.split('.'))
        s=(nsc[0].code, nsc[0][0].code)       
        epr = locations2degrees_origin2channel(origin,nsc[0][0][0])
        receiver_depth_in_km=nsc[0][0][0].elevation/-1000
        origin_depth = origin.depth
        if origin_depth<0:
            origin_depth=1
        arrivals = model.get_travel_times(origin_depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['tts'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        allthes += [numpy.nanmin([ a.time for a in arrivals ])]
    firsts=numpy.nanmin(allthes)

    bulk = []
    done=[]
    if debug:
        print('Requesting:')
    for mseedid in self.get_contents()['channels']:

        nsc = self.select(*mseedid.split('.'))
        s=(nsc[0].code, nsc[0][0].code)
        if s in done:
            continue
        done+=[s]
        
        epr = locations2degrees_origin2channel(origin,nsc[0][0][0])
        receiver_depth_in_km=nsc[0][0][0].elevation/-1000
        origin_depth = origin.depth
        if origin_depth<0:
            origin_depth=1
        arrivals = model.get_travel_times(origin_depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['ttp'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        p = numpy.nanmin([ a.time for a in arrivals ])

        origin_depth = origin.depth
        if origin_depth<0:
            origin_depth=1
        arrivals = model.get_travel_times(origin_depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['tts'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        s = numpy.nanmin([ a.time for a in arrivals ])

        # make selection rectangular
        if s > firsts+endafter:
            if debug:
                print("%s.%s is too far: %f > %f+%f"%(nsc[0].code, nsc[0][0].code,s,firsts,endafter))
            continue

        bulk += [(nsc[0].code, 
                  nsc[0][0].code, 
                  location, 
                  channel, 
                  origin.time+p/2, 
                  origin.time+s+(s-p))]
        if debug:
            print("%s.%s.%s.%s [%s, %s]"%bulk[-1],'Ot+%.3g deg'%epr)
    return bulk

def inventory2waveforms_bulk(self,
                          origin,
                          model = TauPyModel(),
                          debug=False,
                          location='*', 
                          padafter=0,
                          padbefore=0,
                          channel='SN*,SH*,EN*,EH*,HN*,HH*,HG*'):
    bulk = []
    done=[]
    if debug:
        print('Requesting:')
    for mseedid in self.get_contents()['channels']:

        nsc = self.select(*mseedid.split('.'))
        s=(nsc[0].code, nsc[0][0].code)
        if s in done:
            continue
        done+=[s]
        
        epr = locations2degrees_origin2channel(origin,nsc[0][0][0])
        receiver_depth_in_km=nsc[0][0][0].elevation/-1000

        origin_depth = origin.depth
        if origin_depth<0:
            origin_depth=1
        arrivals = model.get_travel_times(origin_depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['ttp'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        p = numpy.nanmin([ a.time for a in arrivals ])

        origin_depth = origin.depth
        if origin_depth<0:
            origin_depth=1
        arrivals = model.get_travel_times(origin_depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['tts'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        s = numpy.nanmin([ a.time for a in arrivals ])
        bulk += [(nsc[0].code, 
                  nsc[0][0].code, 
                  location, 
                  channel, 
                  origin.time+p/2-padbefore, 
                  origin.time+s+(s-p)+padafter)]
        if debug:
            print("%s.%s.%s.%s [%s, %s]"%bulk[-1],'Ot+%.3g deg'%epr)
    return bulk

def remove_response(eventstreams,
                    outputs=['acc'],#['acc','vel','disp']
                    pre_filt=[0.1, 0.33, 45, 50],
                    water_level=60,
                    **detrend_options):
    # remove responses
    
    if 'type' not in detrend_options:
        detrend_options['type']='polynomial'
    if  detrend_options['type']=='polynomial' and 'order' not in detrend_options:
        detrend_options['order']=3


    for output in outputs:#['acc','vel','disp']:

        eventstreams[output] = eventstreams['raw'].copy()
        eventstreams[output].output = output
        eventstreams[output].correction_method='remove_response'
        eventstreams[output].detrend(**detrend_options)#, plot=True)
        for i,tr in enumerate(eventstreams[output]):
            try:
                eventstreams[output][i].remove_response(pre_filt=pre_filt,
                                                    water_level=water_level,
                                                    output=output) # plot=True,
            except:
                print('Cannot remove response:')
                print(tr)
    return eventstreams

def remove_sensitivity(eventstreams,
                       outputs=['acc','vel','disp'],#['acc','vel','disp']
                       filters={'pre':None,#{'type':'highpass','freq':0.075},
                                'acc':None,
                                'vel':None,
                                'disp':{'type':'highpass','freq':1/3.}},
                       integrate_method='cumtrapz',
                       differentiate_method='gradient',
                       sceewenv=False,
                       **detrend_options):
    # remove sensitivity
    
    
    if 'type' not in detrend_options:
        detrend_options['type']='polynomial'
    if  detrend_options['type']=='polynomial' and 'order' not in detrend_options:
        detrend_options['order']=4
    
    if sceewenv:
        tmp = eventstreams['raw'].copy()
        tmp.detrend(**detrend_options)
        tmp.filter(type='lowpass',freq=1/60.)
    
    for output in outputs:#['acc','vel','disp']:

        eventstreams[output] = eventstreams['raw'].copy()
        eventstreams[output].output = output
        eventstreams[output].correction_method='remove_sensitivity'
        if sceewenv:
            for tri,tr in enumerate(eventstreams[output]):
                tr.data = tr.data*1.
                tr.data[:len(tr.data)-2] -= tmp[tri].data[:len(tr.data)-2]
        
        if filters['pre']:
            eventstreams[output].detrend(**detrend_options)#, plot=True)
            eventstreams[output].filter(**filters['pre'])
        eventstreams[output].remove_sensitivity()

        for trace in eventstreams[output]:
            try :
                trace.stats.response.get_paz()
            except:
                print('Cannot get response for')
                print(trace)
                continue
            if 'M/S**2' == trace.stats.response.get_paz().input_units:
                if output in ['acc']:
                    pass
                elif output in ['vel']:
                    eventstreams[output].detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.integrate(method=integrate_method)
                elif output in ['disp']:
                    eventstreams[output].detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.integrate(method=integrate_method)
                    trace.taper(.05,side='left')
                    trace.detrend(**detrend_options)#, plot=True)
                    trace.integrate(method=integrate_method)

            elif 'M/S' == trace.stats.response.get_paz().input_units:
                if output in ['acc']:
                    eventstreams[output].detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.differentiate(method=differentiate_method)
                elif output in ['vel']:
                    pass
                elif output in ['disp']:
                    eventstreams[output].detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.integrate(method=integrate_method)

            elif 'M' == trace.stats.response.get_paz().input_units:
                if output in ['acc']:
                    eventstreams[output].detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.integrate(method=integrate_method)
                elif output in ['vel']:
                    eventstreams[output].detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.integrate(method=integrate_method)
                    trace.detrend(**detrend_options)#, plot=True)
                    trace.taper(.05,side='left')
                    trace.integrate(method=integrate_method)
                elif output in ['disp']:
                    pass
            else:
                print('WARNING: unknown units for trace:')
                print(trace)

        #eventstreams[output].detrend(**detrend_options)#, plot=True)
        if filters[output]:
            for tr in eventstreams[output]:
                try:
                    tr.filter(**filters[output])
                except:
                    print("WARNING: CANNOT FILTER THIS:")
                    print(tr)
                    print(filters[output])
            
    return eventstreams

def clean_inventorystream(inventory,stream,cleaninv=False):
    #Remove traces without response, cha, sta, net without data from inventory
    rmtr = []
    for tr in stream:
        response = None
        if tr.stats.location not in ['PV','PA','PD','EV','EA','ED']:
            try:
                response = inventory.get_response(tr.id, tr.stats.starttime)
            except:
                rmtr += [tr]
                continue
        else:
            try:
                response = inventory.select(network=tr.stats.network,
                                            station=tr.stats.station,
                                            channel=tr.stats.channel.replace('X','*'),
                                            time=tr.stats.starttime)
                mseedid='%s.%s.%s.%s'%(tr.stats.network,
                                       tr.stats.station,
                                       response[-1][-1][-1].location_code,
                                       response[-1][-1][-1].code)
                response = inventory.get_response(mseedid, tr.stats.starttime)
            except:
                rmtr += [tr]
                continue

        if response is None:
            rmtr += [tr]
        elif response.instrument_sensitivity is None:
            rmtr += [tr]
        elif response.instrument_sensitivity.value is None:
            rmtr += [tr]

    for tr in rmtr:
        print('Removing (no response)',tr)
        stream.remove(tr)
    
    if cleaninv:
        rmcha=[]
        for n, net in enumerate(inventory):
            for s,sta in enumerate(net):
                for c,cha in enumerate(sta):
                    mseedid = {'network':net.code,
                            'station':sta.code,
                            #'location':cha.location_code,
                            'channel':cha.code.replace('X','*')}
                    if not len(stream.select(**mseedid)):
                        rmcha+=[mseedid]
                
        for mseedid in rmcha:
            print('Removing',' '.join(["%s: %s"%(k,mseedid[k]) for k in mseedid]))
            tmp = inventory.remove(**mseedid)
            inventory = tmp

    
    print('Stream left:')
    print(stream)
    print('Inventory left:')
    print(inventory)
    
    return inventory,stream

def get_events_waveforms(self,
                         catalog,
                         inventory_client=None,
                         maxStraveltime=20.,
                         location='*', 
                         channel='SN*,SH*,EN*,EH*,HN*,HH*,HG*',
                         model=TauPyModel('iasp91'),
                         quality=None,
                         minimumlength=None,
                         longestonly=None,
                         debug=True,
                         correction_method = remove_sensitivity,
                         padafter=0,
                         padbefore=0,
                         **kwargs):
    
    if inventory_client is None:
        inventory_client=self

    bulkargs = {'quality':quality,
                'minimumlength':minimumlength,
                'longestonly':longestonly}

    inv2bulkargs = {'model':model,
                    'location':location,
                    'channel':channel,
                    'debug':debug}
    
    eventstreams=[]
    eventinventories=[]
    for event in catalog:        

        if False: # Load files if available
            if 'xml' not in files and 'mseed' not in files:
                files='data/%s.quakeml,data/%s.%s.mseed,data/%s.stationxml'(files,files,'%s',files)

            if isinstance(files,str):
                files=files.split(',')

            from obspy import read_events, read_inventory, read
            
            if debug:
                print('Expecting files=<Ev. catalog>,<ev#1 stream (%s for "raw" and "acc")>,<ev#1 inventory>,<ev#2 stream>,<ev#2 inventory>...')
            
            catalog = read_events(files[0])
            eventstreams = [{o:read(f%o) for o in ['raw','acc','vel','disp'] if os.path.exists(f%o) } for f in files[1::2] ]
            eventinventories = [read_inventory(f) for f in files[2::2]]

            return cleandata(catalog,eventstreams,eventinventories)

        origin = event.preferred_origin()

        # Get all stations without resp
        inventory = inventory_client.get_stations(level='channel',
                                                  channel=channel,
                                                  startbefore=origin.time,
                                                  endafter=origin.time+maxStraveltime*9)
        # Get stations with S within first S + maxStraveltime
        bulk = event2stations_bulk(inventory,
                                   origin,
                                   endafter=maxStraveltime,
                                   **inv2bulkargs)
        inventory = inventory_client.get_stations_bulk(bulk,
                                                    level='response',
                                                    **kwargs)

        # Fix response with bad units
        inventory = fixrespunit(inventory,debug=debug)
        eventinventories += [inventory]

        # Get data from ttp/2 to tts+tts-ttp
        bulk = inventory2waveforms_bulk(inventory,
                                        origin,
                                        padafter=padafter,
                                        padbefore=padbefore,
                                        **inv2bulkargs)
        stream = self.get_waveforms_bulk(bulk,**bulkargs)

        # Remove empty (meta)data 
        inventory,stream = clean_inventorystream(inventory,stream)

        # Improve waveforms attributes
        if debug:
            print('Attaching response')
        stream.attach_response(inventory)
        
        if debug:
            print('Attaching trace.stats.distance')
        stream = attach_distance(stream,inventory,origin,debug=debug)
        
        if correction_method:
            if debug:
                print('Correcting amplitudes with',correction_method)
            eventstream = correction_method({'raw':stream})
            if debug:
                print('Correction of amplitudes done.')

        eventstreams += [eventstream]

    return eventstreams,eventinventories
