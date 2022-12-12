#!/usr/bin/env python
from obspy.clients.fdsn import Client 
from obspy.core import UTCDateTime
import os
from addons.get_events import get_events_updatedafter
from addons.get_waveforms import get_events_waveforms,attach_distance,clean_inventorystream
from addons.stream import ploteqdata,combinechannels
from addons.core import eewmap
from addons.catalog import performance_timelines
from addons.event import animate
from matplotlib.patheffects import withStroke
from matplotlib.text import Text
path_effects=[withStroke(linewidth=2,foreground="w")]

def cleandata(catalog,eventstreams,eventinventories):
    for e,event in enumerate(catalog):
        for output, stream in eventstreams[e].items():

            ## Removing bad channels from plots
            for tr in eventstreams[e][output].select(location='99'):
                eventstreams[e][output].remove(tr)
            
            # Remove empty (meta)data 
            eventinventories[e],eventstreams[e][output] = clean_inventorystream(eventinventories[e],eventstreams[e][output])
            
            # Improve waveforms attributes
            eventstreams[e][output].attach_response(eventinventories[e])
            eventstreams[e][output] = attach_distance(eventstreams[e][output],eventinventories[e],event.preferred_origin())

            ## Combine horizontal data
            tridim , horiz = combinechannels(eventstreams[e][output], combine='both')
            eventstreams[e][output] += horiz.select(channel='*b')

    return catalog,eventstreams,eventinventories

def event_data(catalog_uri='USGS',
         inventory_url=None,
         stream_url=None,
         files=None,
         ndays=4,
         nseconds=2*60*60, # Runs twice in case tried every hours 
         user=None, 
         password=None, 
         debug=False, 
         timeout=120*4, 
         channel='SN*,SH*,EN*,EH*,HN*,HH*,HG*',
         eida_token=None,
         **get_events_wargs):

    if files is not None :
        if isinstance(files,str):
            files=files.split(',')

        from obspy import read_events, read_inventory, read
        
        if debug:
            print('Expecting files=<Ev. catalog>,<ev#1 stream (%s for "raw" and "acc")>,<ev#1 inventory>,<ev#2 stream>,<ev#2 inventory>...')
        
        catalog = read_events(files[0])
        eventstreams = [{o:read(f%o) for o in ['raw','acc','vel','disp'] if os.path.exists(f%o) } for f in files[1::2] ]
        eventinventories = [read_inventory(f) for f in files[2::2]]

        return cleandata(catalog,eventstreams,eventinventories)


    if not isinstance(catalog_uri,list):
        catalog_uri = catalog_uri.split(',')

    if inventory_url is None and stream_url is not None:
        inventory_url = stream_url
    elif inventory_url is None :
        inventory_url = catalog_uri[-1]

    if  stream_url is None and  inventory_url is not None:
        stream_url = inventory_url
    elif stream_url is None :
        stream_url = catalog_uri[-1]

    clientargs = {'user':user, 
                  'password':password, 
                  'debug':debug, 
                  'timeout':timeout, 
                  'eida_token':eida_token}

    inventory_client = Client(inventory_url,
                            **clientargs)
    stream_client = Client(stream_url,
                            **clientargs)    

    catalog_clients=[]
    for url in catalog_uri:
        catalog_clients.append(Client(base_url=url,**clientargs))

    if ndays is not None :
        ndays = float(ndays)
    if nseconds is not None :
        nseconds = float(nseconds)
    updated_after = UTCDateTime.now()-nseconds
    if 'starttime' not in get_events_wargs:
        get_events_wargs['starttime'] = UTCDateTime.now()-ndays*24*60*60
    if 'longitude' in get_events_wargs :
        get_events_wargs['longitude'] = float(get_events_wargs['longitude'])
    if 'latitude' in get_events_wargs :
        get_events_wargs['latitude'] = float(get_events_wargs['latitude'])
    if 'maxradius' in get_events_wargs :
        get_events_wargs['maxradius'] = float(get_events_wargs['maxradius'])

    catalog = get_events_updatedafter(catalog_clients,
                                      updated_after=updated_after,
                                      **get_events_wargs)

    eventstreams,eventinventories = get_events_waveforms(stream_client,
                                                         catalog,
                                                         channel=channel,
                                                         inventory_client=inventory_client)
    
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):
        
        shorteventid = str(event.resource_id).split('/')[-1]
        ## Make output dirs
        if not os.path.exists("data"):
            os.makedirs("data")
    
        ## Save event data
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
        
    return cleandata(catalog,eventstreams,eventinventories)

def event_plots(catalog,
                eventstreams,
                eventinventories, 
                plots=True, 
                animates=False, 
                test=False,
                **args):

    saveopt = {'dpi':512,
               'facecolor':'none',
               'transparent':False}
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):

        shorteventid = str(event.resource_id).split('/')[-1]

        if test:
            pass#figs = plotstationdata(streams,event,inventory)
            
        if not plots==False and not plots==0 and not plots=='0' :

            ## Plot data
            fig = ploteqdata(streams['acc'].select(channel='*b'),event,inventory,lim=999)
            fig.tight_layout()
            for ax in fig.axes:
                for t in ax.findobj(Text):
                    if not t.get_path_effects():
                        t.set(path_effects=path_effects) 
            fig.savefig('data/%s_data.png'%shorteventid,bbox_inches=None,**saveopt)
            print('data/%s_data.png'%shorteventid)

            ## Map results
            fig = eewmap({'event':event,
                        'inventory':inventory},
                        reference=False,
                        stationgroups={})
            
            for ax in fig.axes:
                for t in ax.findobj(Text):
                    if not t.get_path_effects():
                        t.set(path_effects=path_effects) 
            fig.savefig('data/%s_map.png'%shorteventid,bbox_inches='tight',**saveopt)
            print('data/%s_map.png'%shorteventid)


            ## Plot results timeline
            fig = performance_timelines(event)
            for ax in fig.axes:
                for t in ax.findobj(Text):
                    if not t.get_path_effects():
                        t.set(path_effects=path_effects) 
            fig.savefig('data/%s_timeline.png'%shorteventid,bbox_inches='tight',**saveopt)
            print('data/%s_timeline.png'%shorteventid)

        if not animates==False and not animates==0 and not animates=='0' :
        
            ## Animate data and results
            anim = animate(event,
                        streams['acc'],
                        streams['disp'],
                        inventory)
            anim.save('data/%s_anim.mp4'%shorteventid)#,dpi=300)
            print('data/%s_anim.mp4'%shorteventid)

def main(**args):

        # Getting data
        catalog,eventstreams,eventinventories = event_data(**args)
        
        # Creating plots
        event_plots(catalog,eventstreams,eventinventories,**args)