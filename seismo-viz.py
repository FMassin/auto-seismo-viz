#!/usr/bin/env python
from obspy.core.event.catalog import Catalog
from obspy.clients.fdsn import Client 
from obspy.core import UTCDateTime
from obspy.clients.fdsn.header import FDSNNoDataException
from sys import argv

def main(fdsnws_uri='USGS',
         ndays=4,
         nseconds=15*60*60,
         areas={'global':{'minmagnitude':7}},
         **kwargs):

    nsecongs_ago = UTCDateTime.now()-nseconds
    myclient = Client(fdsnws_uri,**kwargs)

    ## The catalog in last N days
    catalog_last_ndays = Catalog()
    for i, (area, options) in enumerate(areas.items()):
        options['starttime'] = UTCDateTime.now()-ndays*24*60*60
        try:
            catalog_last_ndays += myclient.get_events(**options)
        except FDSNNoDataException:
            print('No event found in', fdsnws_uri, 'with', options)

    print('Found in the last %d days:'%ndays)    
    print(catalog_last_ndays.__str__(print_all=True))

    ## The catalog in last N days updated in last N seconds
    catalog_updated_last_nseconds = Catalog()
    for e in catalog_last_ndays:
        for m in e.magnitudes:
            if m.creation_info.creation_time > nsecongs_ago:
                eventid = str(e.resource_id).split("?")[-1]
                if '&' in eventid:
                    eventid=[p.split('=')[-1] for p in eventid.split('&') if 'eventid' in p][0]

                try:
                    allsolutions_noarrivals = myclient.get_events(eventid=eventid,
                                                                  includeallmagnitudes=True,
                                                                  includeallorigins=True,
                                                                  includearrivals=False)
                    prefsolution_witharrivals = myclient.get_events(eventid=eventid,
                                                                    includeallmagnitudes=False,
                                                                    includeallorigins=False,
                                                                    includearrivals=True)
                    print('Event %s should be updated'%eventid)
                except:
                    try:
                        allsolutions_noarrivals = myclient.get_events(eventid=eventid,
                                                                    includeallmagnitudes=True,
                                                                    includeallorigins=True)
                        prefsolution_witharrivals = allsolutions_noarrivals
                        print('Cannot use includearrivals but event %s should be updated'%eventid)
                    except FDSNNoDataException:
                        print('No event found in', fdsnws_uri, 'with eventid =',eventid)
                        break

                catalog_updated_last_nseconds += prefsolution_witharrivals[0]
                catalog_updated_last_nseconds[-1].origins = allsolutions_noarrivals[0].origins
                catalog_updated_last_nseconds[-1].magnitudes = allsolutions_noarrivals[0].magnitudes
                break

    print('Updated in the last %d seconds:'%nseconds)    
    print(catalog_updated_last_nseconds.__str__(print_all=True))

if __name__ == '__main__':

    args={ arg.split('=')[0]:arg.split('=')[1] for arg in argv if "=" in arg}
    print('Arguments provided:')
    print('\n'.join(['%s: %s'%(arg,args[arg]) for arg in args]))

    main(**args)