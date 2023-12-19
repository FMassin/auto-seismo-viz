#!/usr/bin/env python
from addons.get_waveforms import remove_sensitivity,attach_distance,clean_inventorystream
from addons.get_events import match_events
from obspy import read, read_inventory, read_events
import os 

def main(catalog=None,
         inventory=None,
         stream=None,
         refcatalog=None,
         debug=True,
         channel='SN*,SH*,EN*,EH*,HN*,HH*,HG*',
         correction_method = remove_sensitivity,
         ):
    
    eventstreams=[]
    eventinventories=[]

    if refcatalog is None:
        catalog = read_events(catalog)
    else:
        ## Get additional solution form another server if any provided
        toadd = read_events(catalog)
        catalog = read_events(refcatalog)
        if debug:
            print('Adding more event info...')#with',client)
        catalog = match_events(catalog,
                                toadd=toadd,
                                lastmagpref=False,
                                includeallmagnitudes=True,
                                includeallorigins=True,
                                includearrivals=False,
                                includecomments=True,
                                filename='/tmp/xml',
                                format='sc3ml',
                                v=debug,
                                x=False)[0]
        
    inventory = read_inventory(inventory)#.select(channel=channel)
    stream = read(stream)#.select(channel=channel)                    

    for event in catalog:
        origin = event.preferred_origin()

        if ((accenv:=stream.select(location='EA')) and
            (velenv:=stream.select(location='EV')) and
            (dispenv:=stream.select(location='ED')) and 
            (acc:=stream.select(location='PA')) and
            (vel:=stream.select(location='PV')) and
            (disp:=stream.select(location='PD'))):

            print('This mseed has been preprocessed by sceewenv (it will be missing single vertical stations)')

            for l in ['EA','EV','ED','PA','PV','PD']:
                for tr in stream.select(location=l):
                    stream.remove(tr)  

            acc += stream.select(channel='HG*')
            acc += stream.select(channel='HN*')
            acc += stream.select(channel='SN*')
            acc += stream.select(channel='EN*')

            vel += stream.select(channel='HH*')
            vel += stream.select(channel='EH*')
            vel += stream.select(channel='SH*')
            vel += stream.select(channel='BH*')

            eventstreams += [{'raw':stream,'acc':acc,'vel':vel,'disp':disp}]
            eventinventories += [inventory]

        else:
            print('Correcting mseed')

            # Remove empty (meta)data 
            eventinventory,eventstream = clean_inventorystream(inventory,stream)

            # Improve waveforms attributes
            if debug:
                print('Attaching response')
            eventstream.attach_response(eventinventory)
            
            if debug:
                print('Attaching trace.stats.distance')
            eventstream = attach_distance(eventstream,eventinventory,origin,debug=debug)
            
            if correction_method:
                if debug:
                    print('Correcting amplitudes with',correction_method)
                eventstreamcorrect = correction_method({'raw':eventstream})
                if debug:
                    print('Correction of amplitudes done.')

            eventstreams += [eventstreamcorrect]
            eventinventories += [eventinventory]
    
    print(zip(catalog,eventstreams,eventinventories))
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):
        
        shorteventid = str(event.resource_id).split('/')[-1]
        ## Make output dirs
        if not os.path.exists("data"):
            os.makedirs("data")
    
        ## Save event data
        if event.preferred_origin() is None:
            try:
                event.preferred_origin_id = event.origins[-1].resourceid
            except:
                print('NO PREF ORG!')
        if event.preferred_magnitude() is None:
            try:
                event.preferred_magnitude_id = event.magnitudes[-1].resourceid
            except:
                print('NO PREF MAG!')

        event.write('data/%s.quakeml'%shorteventid,
                    format='quakeml')
        print('data/%s.quakeml'%shorteventid)

        ## Save event metadata
        inventory.write('data/%s.stationxml'%shorteventid,
                        format='STATIONXML')
        print('data/%s.stationxml'%shorteventid)
        
        ## Save event seismic data
        for output, stream in streams.items():
            # {'raw': ..., 'acc': ...,'vel': ..., 'disp': ...}
            stream.write('data/%s.%s.mseed'%(shorteventid,
                                             output),
                         format='mseed')
            print('data/%s.%s.mseed'%(shorteventid,output))