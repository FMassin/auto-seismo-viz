#!/usr/bin/env python
from obspy.clients.fdsn import Client 
from obspy.core import UTCDateTime
from sys import argv
from re import sub
import os,glob
from addons.get_events import get_events_updatedafter
from addons.get_waveforms import get_events_waveforms,attach_distance
from addons.stream import ploteqdata,combinechannels
from addons.core import eewmap
from addons.catalog import performance_timelines
from addons.event import animate

def event_data(catalog_uri='USGS',
         inventory_url=None,
         stream_url=None,
         files=None,
         ndays=4,
         nseconds=60*60, # Updates 3 times at most in case it's ran every 15 min
         user=None, 
         password=None, 
         debug=False, 
         timeout=120, 
         eida_token=None,
         **get_events_wargs):

    if files is not None :
        if isinstance(files,str):
            files=files.split(',')

        from obspy import read_events, read_inventory, read
        
        if debug:
            print('Expecting files=<Ev. catalog>,<ev#1 stream (%s for "raw" and "acc")>,<ev#1 inventory>,<ev#2 stream>,<ev#2 inventory>...')
        
        catalog = read_events(files[0])
        eventstreams = [{o:read(f%o) for o in ['raw','acc','disp'] } for f in files[1::2] ]
        eventinventories = [read_inventory(f) for f in files[2::2]]
        return catalog,eventstreams,eventinventories


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
    # To do : Remove cha, sta, net without data from inventory
    eventstreams,eventinventories = get_events_waveforms(stream_client,
                                                         catalog,
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
        
    return catalog,eventstreams,eventinventories

def event_plots(catalog,eventstreams,eventinventories):
    
    saveopt = {'dpi':512,
               'facecolor':'none',
               'transparent':False}
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):

        shorteventid = str(event.resource_id).split('/')[-1]

        ## Post processing seismic data
        for output, stream in streams.items():

            ## Removing bad channels from plots
            for tr in streams[output].select(location='99'):
                streams[output].remove(tr)
            
            # Improve waveforms attributes
            stream.attach_response(inventory)
            streams[output] = attach_distance(streams[output],inventory,event.preferred_origin())

            ## Combine horizontal data
            tridim , horiz = combinechannels(streams[output], combine='both')
            streams[output] += horiz.select(channel='*b')

        if True:
            ## Plot data
            fig = ploteqdata(streams['acc'].select(channel='*b'),event,inventory,lim=999)
            fig.tight_layout()
            fig.savefig('data/%s_data.png'%shorteventid,**saveopt)
            print('data/%s_data.png'%shorteventid)


            ## Map results
            fig = eewmap({'event':event,
                        'inventory':inventory},
                        reference=False,
                        stationgroups={})
            fig.savefig('data/%s_map.png'%shorteventid,bbox_inches='tight',**saveopt)
            print('data/%s_map.png'%shorteventid)


            ## Plot results timeline
            # To do : Fix issue with time axis
            fig = performance_timelines(event)
            fig.savefig('data/%s_timeline.png'%shorteventid,bbox_inches='tight',**saveopt)
            print('data/%s_timeline.png'%shorteventid)


        else:
            ## Animate data and results
            anim = animate(event,
                        streams['acc'],
                        streams['disp'],
                        inventory)
            anim.save('data/%s_anim.mp4'%shorteventid)#,dpi=300)
            print('data/%s_anim.mp4'%shorteventid)

if __name__ == '__main__':

    args={ arg.split('=')[0]:arg.split('=')[1] for arg in argv[1:] if "=" in arg}
    tmp='\n'.join(['%s: %s'%(arg,args[arg]) for arg in args])
    tmp=sub('catalog_uri: .*\n',  'catalog_uri: *****\n',tmp)
    tmp=sub('stream_url: .*\n',   'stream_url: *****\n',tmp)
    tmp=sub('inventory_url: .*\n','inventory_url: *****\n',tmp)
    tmp=sub('user:.*\n',       'user: *****\n',      tmp)
    tmp=sub('password:.*\n',   'password: *****\n',  tmp)
    print('Arguments provided:')
    print(tmp)

    # Creating plots
    catalog,eventstreams,eventinventories = event_data(**args)
    
    # Creating plots
    event_plots(catalog,eventstreams,eventinventories)

    # To do:
    # - add agency
    # - webhook notif
    # - html or panel report
    # - station and event tables
    # - pip install
    # - BETTER OPTIONS:
    #from getopt import getopt, GetoptError
    #py2 = sys.version_info < (3,)
    #pipedinputs = sys.stdin if py2 else sys.stdin.buffer
    #... = None
    #
    #try:
    #    opts, args = getopt(sys.argv[1:], "x:...:",
    #                        ["stdout", "...=", ])
    #except GetoptError:
    #    usage()
    #    return 1
    #
    #out_channel = None
    #
    #
    #for flag, arg in opts:
    #    if flag in ("-c", "--stdout"):
    #        out_channel = sys.stdout if py2 else sys.stdout.buffer
    #    elif ... in ("-...", "--..."):
    #        ... = arg
    #    else:
    #        usage()
    #        if flag in ("-h", "--help"):
    #            return 0
    #        return 1
    #
    #if ...
    #    sys.exit(1)

