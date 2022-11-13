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
                for a,stage in enumerate(cha.response.response_stages):
                    input_units = cha.response.instrument_sensitivity.input_units.replace('m/s','M/S')
                    inventory[n][s][c].response.instrument_sensitivity.input_units=input_units
                    if stage.input_units is not None :#and "m/s" in stage.input_units:
                        input_units = stage.input_units.replace('m/s','M/S')
                        inventory[n][s][c].response.response_stages[a].input_units = input_units
    
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
        arrivals = model.get_travel_times(origin.depth/1000,
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
        arrivals = model.get_travel_times(origin.depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['ttp'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        p = numpy.nanmin([ a.time for a in arrivals ])

        arrivals = model.get_travel_times(origin.depth/1000,
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
        arrivals = model.get_travel_times(origin.depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['ttp'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        p = numpy.nanmin([ a.time for a in arrivals ])

        arrivals = model.get_travel_times(origin.depth/1000,
                                            distance_in_degree=epr,
                                            phase_list=['tts'],
                                            #receiver_depth_in_km=receiver_depth_in_km
                                            )
        s = numpy.nanmin([ a.time for a in arrivals ])
        bulk += [(nsc[0].code, 
                  nsc[0][0].code, 
                  location, 
                  channel, 
                  origin.time+p/2, 
                  origin.time+s+(s-p))]
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
            eventstreams[output].filter(**filters[output])            
            
    return eventstreams

def clean_inventorystream(inventory,stream):
    #Remove traces without response, cha, sta, net without data from inventory
    rmtr = []
    for tr in stream:
        try:
            inventory.get_response(tr.id, tr.stats.starttime)
        except:
            rmtr += [tr]
    for tr in rmtr:
        print('Removing (no response)',tr)
        stream.remove(tr)

    rmcha=[]
    for n, net in enumerate(inventory):
        for s,sta in enumerate(net):
            for c,cha in enumerate(sta):
                mseedid = {'network':net.code,
                           'station':sta.code,
                           'location':cha.location_code,
                           'channel':cha.code}
                if not len(stream.select(**mseedid)):
                    #print('Cannot find data for',mseedid)
                    rmcha+=[mseedid]
            
    for mseedid in rmcha:
        print('Removing',' '.join(["%s: %s"%(k,mseedid[k]) for k in mseedid]))
        tmp = inventory.remove(**mseedid)
        inventory = tmp
    
    return inventory,stream

def get_events_waveforms(self,
                         catalog,
                         inventory_client=None,
                         maxStraveltime=40.,
                         location='*', 
                         channel='SN*,SH*,EN*,EH*,HN*,HH*,HG*',
                         model=TauPyModel('iasp91'),
                         quality=None,
                         minimumlength=None,
                         longestonly=None,
                         debug=True,
                         correction_method = remove_sensitivity,
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