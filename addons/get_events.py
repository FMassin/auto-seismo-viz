#!/usr/bin/env python
from obspy import read_events
from obspy.core.event.catalog import Catalog
from obspy.clients.fdsn.header import FDSNNoDataException
import numpy
import addons.core

def match_events(self=Catalog(),
                 toadd=Catalog(),
                 client=None,
                 v=False,
                 x=True,
                 lastmagpref=True,
                 **get_events_args):
    """
        Matchs the events in self (obspy:Catalog) with specified catalog.

        Matchs with catalog available in webservice if client is specified.
    """

    tofind = toadd.copy()
    matchs = Catalog()
    extras = Catalog()
    matchs.description = 'Intersections of '+str(self.description)+' and '+str(tofind.description)
    missed = Catalog()
    missed.description = 'Part of '+str(self.description)+' not in '+str(tofind.description)
    listfound=list()
    
    if v:
        if client:
            print('Reference client')
            print(client)

    for e in self.events :

        o = e.preferred_origin() 
        if o is None:
            o = e.origins[-1]


        found=False
        memdl=99999999
        memdt=99999999
        if v:
            print('---- event',e.short_str(),'...')

        if len(listfound) < len(self.events):

            filter = ["time >= "+str(o.time-120),
                      "time <= "+str(o.time+120),
                      "latitude >= "+str(o.latitude-2.5),
                      "latitude <= "+str(o.latitude+2.5),
                      "longitude >= "+str(o.longitude-2.5),
                      "longitude <= "+str(o.longitude+2.5)]

            if client :
                eq_specs = {'starttime':str(o.time-120),
                            'endtime':str(o.time+120),
                            'minlatitude':max([-90,o.latitude-2.5]),
                            'maxlatitude':min([90, o.latitude+2.5]),
                            'minlongitude':max([-180,o.longitude-2.5]),
                            'maxlongitude':min([180,o.longitude+2.5])}
                try:
                    tofind = client.get_events( **eq_specs)
                    if v:
                        print('Event in reference client')
                        print(tofind)
                except:
                    if v:
                        print('No event in reference client')
                    continue
                if len(tofind.events)==0:
                    try:
                        tofind = read_events(get_events_args['filename'],
                                                   format=get_events_args['format'].upper())
                    except:
                        if v:
                            print('can t read filename')
                            print(get_events_args['filename'])
                            return
                        pass
            if v:
                print(tofind)
                print(filter)
            for candidate in tofind.filter( *filter ).events:
                if v:
                    print(candidate)
                candidateo = candidate.preferred_origin()
                if o is  None:
                    candidateo = candidate.origins[-1]

                if candidateo.time not in listfound:

                    dt = abs(o.time-candidateo.time)
                    dl = addons.core.haversine(o.longitude,
                                                    o.latitude,
                                                    candidateo.longitude,
                                                    candidateo.latitude)
                    if v:
                        print(dt,dl)
                    if (dt < 10 and dl <=50000 and
                        (dl<memdl or dt<memdt)):

                        found = True
                        memi = str(candidateo.time)
                        memdl = dl
                        memdt = dt
                        meme = candidate
                        if v:
                            print('fits nicely input catalog ',tofind.description,':\n  ', candidate.short_str())

                        break

                    elif (dt < 40 and dl <=100000 and
                          (dl<memdl or dt<memdt)):

                        found = True
                        memi = str(candidateo.time)
                        memdl = dl
                        memdt = dt
                        meme = candidate
                        if v:
                            print('fits input catalog ',tofind.description,':\n  ', candidate.short_str())

                    elif (dt < 80 and dl <200000 and
                          (dl<memdl or dt<memdt)):

                        found = True
                        memi = str(candidateo.time)
                        memdl = dl
                        memdt = dt
                        meme = candidate
                        if v:
                            print('poorly fits input catalog ',tofind.description,':\n  ', candidate.short_str())


        if not found:
            if v:
                print('does not exist in current catalog')
            missed.events.append(e)

        elif found:
            if v:
                print('merged with ', meme.short_str())
            #matchs.events.append(e)
            o = meme.preferred_origin()
            if o is  None:
                o = meme.origins[-1]

            
            if client is not None:
                meme = client.get_events(starttime=str(o.time-1),
                                        endtime=str(o.time+1),
                                        **get_events_args )

                if 'filename' in get_events_args.keys():
                    f = open(get_events_args['filename'],'r')
                    filedata = f.read()
                    f.close()

                    newdata = filedata.replace('0.10" version="0.10"','0.9" version="0.9"')
                    newdata = filedata.replace('0.11" version="0.11"','0.9" version="0.9"')
                    newdata = filedata.replace('0.12" version="0.12"','0.9" version="0.9"')

                    f = open(get_events_args['filename'],'w')
                    f.write(newdata)
                    f.close()

                    meme = read_events(get_events_args['filename'],
                                            format=get_events_args['format'].upper())

                meme=meme[0]

            matchs.events.append(meme)
            listfound.append(memi)

            for prefatt in [ 'preferred_origin_id', 'preferred_magnitude_id' ]:
                if hasattr(e, prefatt):
                    matchs.events[-1][prefatt] = e[prefatt]

            #for prefatt in [ 'preferred_origin_id', 'preferred_magnitude_id' ]:
            #    if hasattr(meme, prefatt):
            #        matchs.events[-1][prefatt] = meme[prefatt]

            print(e)
            for listatt in [ 'origins' , 'magnitudes', 'picks', 'station_magnitudes','focal_mechanisms' ]:
            #    if hasattr(meme, listatt):
            #        matchs.events[-1][listatt].extend( meme[listatt] )

                if hasattr(e, listatt):
                    for i,att in enumerate(e[listatt] ):
                        matchs.events[-1][listatt].insert(i,att)

            if matchs.events[-1].preferred_magnitude() is None:
                print('No pref mag')
                print(matchs.events[-1])
                continue

            if not hasattr(matchs.events[-1].preferred_magnitude(), 'origin_id'):
                matchs.events[-1].preferred_magnitude().origin_id=matchs.events[-1].preferred_origin().resource_id

            if matchs.events[-1].preferred_magnitude().origin_id is None:
                matchs.events[-1].preferred_magnitude().origin_id=matchs.events[-1].preferred_origin().resource_id

            if False:
                matchs.events[-1].preferred_origin=matchs.events[-1].preferred_origin_id.get_referred_object
                matchs.events[-1].preferred_magnitude=matchs.events[-1].preferred_magnitude_id.get_referred_object
                if matchs.events[-1].preferred_focal_mechanism_id is not None:
                    matchs.events[-1].preferred_focal_mechanism=matchs.events[-1].preferred_focal_mechanism_id.get_referred_object

            #tofind.events.remove(meme)
    if lastmagpref:
        for e in matchs:
            e.preferred_origin_id = e.origins[-1].resource_id
            e.preferred_magnitude_id = e.magnitudes[-1].resource_id
    if x :
        matchs, extras, trash = match_events(self=tofind,
                                                          toadd=self,
                                                          client=client,
                                                          v=v,
                                                          x=False)

    return matchs, missed, extras

def get_events_updatedafter(self,
                            updated_after=None,
                            debug=True,
                            **kwargs):

    if debug:
        print('//////// get_events_updatedafter ////////')
    refclient = self
    if isinstance(self,list):
        if debug:
            print('List of clients provided')
        refclient=self[0]

    ## The catalog in last N days
    catalog_last_ndays = Catalog()
    if 'eventid' in kwargs:
        for badkey in ['starttime','endtime','limit']:
            if  badkey in  kwargs:
                kwargs.pop(badkey)
    print(kwargs)
    try:
        catalog_last_ndays += refclient.get_events(**kwargs)
    except FDSNNoDataException:
        if debug:
            print('No event found with', kwargs)

    if debug:
        print('Found:')    
        print(catalog_last_ndays.__str__(print_all=True))

    ## The catalog in last N days updated in last N seconds
    catalog_updated_last_nseconds = Catalog()
    for e in catalog_last_ndays:
        for m in e.origins:
            if m.creation_info.creation_time > updated_after :
                eventid=str(e.resource_id)
                if not 'smi:' == eventid[:4]:
                    eventid = str(e.resource_id).split("?")[-1]
                    if '&' in eventid:
                        eventid=[p.split('=')[-1] for p in eventid.split('&') if 'eventid' in p][0]
                    
                    eventid = eventid.split("/")[-1]
                usgs = False
                try:
                    allsolutions_noarrivals = refclient.get_events(eventid=eventid,
                                                                  includeallmagnitudes=True,
                                                                  includeallorigins=True,
                                                                  includearrivals=False)
                    prefsolution_witharrivals = refclient.get_events(eventid=eventid,
                                                                    includeallmagnitudes=False,
                                                                    includeallorigins=False,
                                                                    includearrivals=True)
                
                except FDSNNoDataException:
                    try:
                        allsolutions_noarrivals = refclient.get_events(eventid=eventid.split("/")[-1],
                                                                    includeallmagnitudes=True,
                                                                    includeallorigins=True,
                                                                    includearrivals=False)
                        prefsolution_witharrivals = refclient.get_events(eventid=eventid.split("/")[-1],
                                                                        includeallmagnitudes=False,
                                                                        includeallorigins=False,
                                                                        includearrivals=True)
                    except FDSNNoDataException:
                        try:
                            eventid = 'us'+eventid
                            print('Trying USGS style with no origin or magnitude specification, and eventid=', eventid)
                            allsolutions_noarrivals = refclient.get_events(eventid=eventid,
                                                                        #includeallmagnitudes=True,
                                                                        #includeallorigins=True,
                                                                        #includearrivals=False
                                                                        )
                            prefsolution_witharrivals = allsolutions_noarrivals
                            usgs=True
                        except FDSNNoDataException:
                            if debug:
                                print('No event found with eventid=', eventid, 'or', eventid[2:])
                            break

                except:
                    if debug:
                        print('Cannot get eventid=', eventid,'... Trying without arrivals')
                    
                    try:
                        allsolutions_noarrivals = refclient.get_events(eventid=eventid,
                                                                    includeallmagnitudes=True,
                                                                    includeallorigins=True)
                        prefsolution_witharrivals = allsolutions_noarrivals
                        if debug:
                            print('Cannot use includearrivals')

                    except FDSNNoDataException:
                        if debug:
                            print('No event found with eventid =',eventid,'(includearrivals not used)')
                        break
                if debug:
                    print('Event %s should be updated'%eventid)

                catalog_updated_last_nseconds += prefsolution_witharrivals[0]
                if not usgs:
                    catalog_updated_last_nseconds[-1].origins = allsolutions_noarrivals[0].origins
                    catalog_updated_last_nseconds[-1].magnitudes = allsolutions_noarrivals[0].magnitudes

                if not isinstance(self,list):
                    break

                ## Get additional solution form another server if any provided
                for client in self[1:]:
                    if debug:
                        print('Adding more event info...')#with',client)
                    catalog_updated_last_nseconds = match_events(catalog_updated_last_nseconds,
                                                                client=client,
                                                                lastmagpref=False,
                                                                includeallmagnitudes=True,
                                                                includeallorigins=True,
                                                                includearrivals=False,
                                                                includecomments=True,
                                                                filename='/tmp/xml',
                                                                format='sc3ml',
                                                                v=debug,
                                                                x=False)[0]
                break

    if debug:
        print('Updated after %s:'%updated_after)    
        print(catalog_updated_last_nseconds.__str__(print_all=True))

    return catalog_updated_last_nseconds