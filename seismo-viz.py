#!/usr/bin/env python
from obspy.clients.fdsn import Client 
from obspy.core import UTCDateTime
from sys import argv
from re import sub
import os

from addons.get_events import get_events_updatedafter

def main(catalog_uri='USGS',
         inventory_url=None,
         stream_url=None,
         ndays=4,
         nseconds=15*60*60,
         user=None, 
         password=None, 
         debug=False, 
         timeout=120, 
         eida_token=None,
         **get_events_wargs):

    if not isinstance(catalog_uri,list):
        catalog_uri = catalog_uri.split(',')

    if inventory_url is None and stream_url is not None:
        inventory_url = stream_url
    elif inventory_url is None :
        inventory_url = catalog_uri[0]

    if  stream_url is None and  inventory_url is not None:
        stream_url = inventory_url
    elif stream_url is None :
        stream_url = catalog_uri[0]

    if ndays is not None :
        ndays = float(ndays)

    if nseconds is not None :
        nseconds = float(nseconds)
    catalog_clients=[]
    for url in catalog_uri:
        catalog_clients.append(Client(base_url=url,
                                user=user, 
                                password=password, 
                                debug=debug, 
                                timeout=timeout, 
                                eida_token=eida_token))

    updated_after = UTCDateTime.now()-nseconds
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

    

    for event in catalog:
        
        ## Make output dirs
        if not os.path.exists("data"):
            os.makedirs("data")
            
        ## Save Event xml
        event.write('data/%s.quakeml'%str(event.resource_id).split('/')[-1],
                    format='quakeml')

        ## Get seismic (meta)data
        #

        ## Plot data
        #

        ## Map results
        #

        ## Plot results timeline
        #

        ## Animate data and results
        #

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

    main(**args)