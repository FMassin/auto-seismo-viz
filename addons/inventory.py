#!/usr/bin/env python
from obspy.core.inventory.inventory import Inventory
import matplotlib.colors
import numpy
import addons.core

def get_best_instrument(self,
                        instruments_markers,
                        preforder=['HG','HN','HH','EN','SN','BH','EH','SH','HL'],
                        rank=-1):

    channels = self.get_contents()['channels']
    found=[]
    for instrument_type in preforder :
        if (instrument_type in instruments_markers.keys() and
            instrument_type in [str(cs.split('.')[-1][:2]) for cs in channels] and
            instrument_type not in found):
            found += [instrument_type]
            #return instrument_type
    if len(found):
        return '-'.join(found[:rank])
    else:
        print('None preferred in',channels)
        return 'none'

def distfilter(self=Inventory([],''),
               dmax=None,
               dmin=0.,
               x1=None,
               y1=None,
               x2=None,
               y2=None,
               z1=None,
               save=True,
               out=False):
    if not x1 or not y1 :
        y1 = numpy.mean([s.latitude for n in self.networks for s in n.stations])
        x1 = numpy.mean([s.longitude for n in self.networks for s in n.stations])
    if not y1 or not y2:
        y2 = y1+0.01
        x2 = x1+0.01

    distances=list()
    for n in self.networks:
        d=numpy.nan
        for s in n.stations:
            if (x1 and y1) or (x2 and y2):
                d = 110.*addons.core.DistancePointLine(s.longitude, s.latitude, x1, y1, x2, y2)

            if z1 is not None:
                d = numpy.sqrt(d**2. + (s.elevation/1000-z1)**2)

            distances.append(d)

            if d>=dmin and d<=dmax:
                s.code+='_inside'
            else:
                s.code+='_outside'

    inside = self.select(station='*_inside')
    outside = self.select(station='*_outside')

    for n in self.networks:
        for s in n.stations:
            s.code = s.code.replace("_inside", '')
            s.code = s.code.replace("_outside", '')

    for n in inside.networks:
        for s in n.stations:
            s.code = s.code.replace('_inside', '')

    for n in outside.networks:
        for s in n.stations:
            s.code = s.code.replace('_outside', '')
    if out in ['d','dist','distance']:
        return distances
    else:
        if out:
            return inside, outside
        else:
            return inside

def make_datacaption(self,
                     dim,
                     instruments_markers,
                     instruments_captions,
                     orientations_markers,
                     orientations_captions,
                     samplerates_markers,
                     samplerates_captions,
                     stationgroups=None,
                     rank=-1,
                     legend_sum=True):
    data = list()
    captions=list()
    counts={}
    for n in  self.networks:
        for s in n.stations:
            if dim == 'networks':
                data.append(n.code)
                captions.append(n.code)
            elif dim == 'locations':
                data.append(s.code)
                captions.append(s.code)
            elif dim == 'elevation':
                data.append(self.get_coordinates(s)['elevation'])
                captions.append(self.get_coordinates(s)['elevation'])
            elif dim == 'local_depth':
                data.append(self.get_coordinates(s)['local_depth'])
                captions.append(self.get_coordinates(s)['local_depth'])
            elif dim == 'instruments':
                best_instrument=get_best_instrument(s,instruments_markers,rank=rank)
                if best_instrument[:2] in instruments_markers:
                    data.append(instruments_markers[best_instrument[:2]])
                cap = []
                for bi in best_instrument.split('-'):
                    if bi not in instruments_captions:
                        pass#cap+=['']
                    elif not instruments_captions[bi] in cap:
                        cap+=[instruments_captions[bi]]
                if len(cap)>1:
                    pass#print(n.code,s.code, cap)
                captions.append('-'.join(cap))
            elif dim == 'orientations':
                best_orientation=get_best_orientation(s,orientations_markers)
                data.append(orientations_markers[best_orientation])
                captions.append(orientations_captions[best_orientation])
            elif dim == 'sample_rates':
                best_samplerate=get_best_samplerate(s,samplerates_markers)
                data.append(samplerates_markers[best_samplerate])
                captions.append(samplerates_captions[best_samplerate])
            
            if stationgroups is not None:
                for key in stationgroups.keys():
                    for code in stationgroups[key]:
                        if not key in counts.keys():
                            counts[key]={}
                        if not captions[-1] in counts[key].keys():
                            counts[key][captions[-1]] = 0
                        if code == ('%s.%s'%(n.code,s.code))[:len(code)]:
                            counts[key][captions[-1]] += 1

    for icap,cap in enumerate(captions):
        if legend_sum:
            detail = ', '.join(['%s %s'%(counts[key][cap], key) for key in counts.keys() if counts[key][cap]])
        else:
            detail = ', '.join(['%s'%(key) for key in counts.keys() if counts[key][cap]])
        if len(detail):
            captions[icap] += ' (%s)'%(detail)
        if 'ther' in detail:
            print(icap, cap)

    if dim == 'instruments':
        for instrument_type in instruments_markers.keys():
            data.append(instruments_markers[instrument_type])
    elif dim == 'orientations':
        for orientation_type in orientations_markers.keys():
            data.append(orientations_markers[orientation_type])

    return data,captions

def codes2nums(data,
               used=list()):
    #print(data)
    try :
        data[0]*1.
        nums = data
    except:
        nums=list()
        for code in data:
            if code not in used:
                used.append(code)
            nums.append(used.index(code))
    #print(nums)
    return nums

def scatter6d(self,
              longitudes,
              latitudes,
              sizes,
              colors,
              markers,
              captions,
              captions_dim,
              legend_sum=True,
              **kwargs):

    longitudes, latitudes = self(longitudes, latitudes) 

    if not isinstance(colors[0], str)  :
        norm = matplotlib.colors.Normalize()
        norm.autoscale(colors)
        cm = matplotlib.pyplot.get_cmap('nipy_spectral')

    used=list()
    for i,m in enumerate(markers):
        if isinstance(colors[i],str):
            colors[i] = colors[i][:2]
        l=None
        if captions[i] not in used and captions_dim not in ['none']:
            n = sum([1 for c in captions if captions[i]==c])
            l = '%s' % (captions[i])
            if legend_sum:
                l = '%s %s' % (n, captions[i])
            used.append(captions[i])

        if isinstance(colors[i], str)  :
            try:
                self.ax.scatter(x=longitudes[i],
                             y=latitudes[i],
                             s=sizes[i],
                             c=colors[i],
                             marker=markers[i],
                             label=l,
                             **kwargs)
            except:
                print('FAIL! What is wrong in self.scatter inputs below???')
                print('x=',longitudes[i])
                print('y=',latitudes[i])
                print('s=',sizes[i])
                print('c=',colors[i])
                print('marker=',markers[i])
                print('label=',l)
                print('**kwargs=',**kwargs)
        else:
            self.ax.scatter(x=longitudes[i],
                         y=latitudes[i],
                         s=sizes[i],
                         c=cm(norm(colors[i])),
                         marker=markers[i],
                         label=l,
                         **kwargs)

def map_stations(self=Inventory([],''),
                 bmap=None,
                 fig=None,
                 stationgroups = {'NU':['NU.'],
                                  #'CAM':['G.','GI.','SV.','HN.','NU.','OV.','TC.'],
                                  'ECOS':['CH','C4']},
                 colors='instruments',
                 markers='instruments',
                 sizes='none',
                 rank=1,
                 titletext='',
                 markersize=30,
                 fontsize=8,
                 instruments_markers= {'SP':'<','SH':'<', 'EH':'<','HL':'v',
                                       'BB':'v','BB-SP':'v','SP-BB':'v','HH':'v', 'BH':'v',
                                       'SM':'^','SM-SP':'^','SP-SM':'^','BB-SM':'^','SM-BB':'^',
                                       'HG':'^',   'HN':'^',   'EN':'^',   'SN':'^',
                                       'none':'*'},
                 instruments_captions={'SH':'SP','EH':'SP',
                                       'BH':'BB','HH':'BB','HL':'BB',
                                       'HG':'SM','HN':'SM','EN':'SM','SN':'SM',
                                       'none':'Other'},
                 orientations_markers={'N':'^','2':'s','Z':'P','none':'*'},
                 orientations_captions={'N':'triax.','2':'hori.','Z':'vert.','none':'other'},
                 samplerates_markers={100:'^',40:'s','none':'*'},
                 samplerates_captions={100:'Hi rate',40:'Low rate','none':'other'},
                 stations_colors=None,
                 stations_markers=None,
                 stations_sizes=None,
                 stations_colordata=None,
                 stations_markerdata=None,
                 stations_sizedata=None,
                 filled_markers = ('^', 'v', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','o'),
                 dateintitle=False,
                 stationintitle=False,
                 legend_sum=True,
                 mapbounds=None
                ):

    stations_longitudes, stations_latitudes = addons.core.search(self,
                                                                  fields=['longitude','latitude'],
                                                                  levels=['networks','stations'])
    if not stations_sizedata:
        if sizes in ['none']:
            stations_sizedata = numpy.zeros(len(stations_longitudes))+markersize
            stations_sizecaptions = numpy.zeros(len(stations_longitudes))+markersize
        else:
            stations_sizedata, \
                stations_sizecaptions=make_datacaption(self,
                                                        sizes,
                                                        instruments_markers,
                                                        instruments_captions,
                                                        orientations_markers,
                                                        orientations_captions,
                                                        samplerates_markers,
                                                       samplerates_captions,
                                                       stationgroups=stationgroups)


    if not stations_colordata:
        if colors in ['none']:
            stations_colordata = numpy.zeros(len(stations_longitudes))
            stations_colorcaptions = numpy.zeros(len(stations_longitudes))
        else:
            stations_colordata, \
                stations_colorcaptions=make_datacaption(self,
                                                        colors,
                                                        instruments_markers,
                                                        instruments_captions,
                                                        orientations_markers,
                                                        orientations_captions,
                                                        samplerates_markers,
                                                        samplerates_captions,
                                                        stationgroups=stationgroups,
                                                        rank=rank,                         legend_sum=legend_sum)

    if not stations_markerdata:
        if markers in ['none']:
            stations_markerdata = numpy.repeat(filled_markers[0],len(stations_longitudes))
            stations_markercaptions = numpy.repeat(filled_markers[0],len(stations_longitudes))
        else:
            stations_markerdata, \
                stations_markercaptions =make_datacaption(self,
                                                          markers,
                                                          instruments_markers,
                                                          instruments_captions,
                                                          orientations_markers,
                                                          orientations_captions,
                                                          samplerates_markers,
                                                          samplerates_captions,
                                                          stationgroups=stationgroups,
                                                           rank=rank,                         legend_sum=legend_sum)

    #make stations_colors colors and stations_markers markers
    if not stations_sizes:
        stations_sizes= codes2nums(stations_sizedata)
    if not stations_colors:
        stations_colors = codes2nums(stations_colordata,
                                                  used = [instruments_markers[k] for k in instruments_markers.keys()] )
        indexes = numpy.sort(list(set(stations_colors)))
        compress = {}
        for index,colorindex in enumerate(indexes):
            compress[colorindex]=index
        for index,colorindex in enumerate(stations_colors):
            stations_colors[index]=compress[colorindex]
        for index,colorindex in enumerate(stations_colors):
            stations_colors[index]
        for index,colorindex in enumerate(stations_colors):
            if colorindex > 9:
                stations_colors[index] = 'k'
            else:
                stations_colors[index] = 'C'+str(colorindex)
        #stations_colors = ['C'+str(colorindex) for colorindex in stations_colors]

        
    if not stations_markers:
        stations_markers = [ filled_markers[min([d,len(filled_markers)])] for d in codes2nums(stations_markerdata) ]
        for k,m in enumerate(stations_markers):
            cap=stations_colorcaptions[k].split(' ')[0]
            if cap in instruments_markers:
                stations_markers[k]=instruments_markers[cap]
    
    
    if mapbounds is not None:
        print('all stations:',len(stations_longitudes))

        stations_longitudes = numpy.asarray(stations_longitudes)
        stations_latitudes = numpy.asarray(stations_latitudes)
        stations_sizes = numpy.asarray(stations_sizes)
        stations_colors = numpy.asarray(stations_colors)
        stations_markers = numpy.asarray(stations_markers)
        stations_sizecaptions = numpy.asarray(stations_sizecaptions)
        stations_markercaptions = numpy.asarray(stations_markercaptions)
        stations_colorcaptions = numpy.asarray(stations_colorcaptions)
        if False:
            if markers is not None:
                markers = numpy.asarray(markers)
            if sizes is not None:
                colors = numpy.asarray(colors)
            if sizes is not None:
                sizes = numpy.asarray(sizes)

        mask = (stations_longitudes>=min(mapbounds[0])) * (stations_longitudes<=max(mapbounds[0])) * (stations_latitudes>=min(mapbounds[1])) * (stations_latitudes<=max(mapbounds[1]))

        stations_longitudes = list(stations_longitudes[mask])
        stations_latitudes = list(stations_latitudes[mask])
        stations_sizes = list(stations_sizes[mask])
        stations_colors = list(stations_colors[:len(mask)][mask])
        stations_markers = list(stations_markers[mask])
        stations_sizecaptions = list(stations_sizecaptions[mask])
        stations_markercaptions = list(stations_markercaptions[mask])
        stations_colorcaptions = list(stations_colorcaptions[mask])
        if False:
            if markers is not None:
                markers = list(markers[mask])
            if colors is not None:
                colors = list(colors[mask])
            if sizes is not None:
                sizes = list(sizes[mask])
        print('Filtered stations:',len(stations_longitudes))

    #print(stations_latitudes)
    #print(stations_sizes)
    #print(stations_colors)
    #print(stations_markers)
    #print(stations_sizecaptions)
     
    scatter6d(bmap,
              stations_longitudes,
              stations_latitudes,
              stations_sizes,
              stations_colors,
              stations_markers,
              stations_sizecaptions,
              sizes,
              facecolor='w',
              edgecolor='w',
              linewidths=2.5,
              zorder=97,
                            legend_sum=legend_sum
              )
   
    scatter6d(bmap,
              stations_longitudes,
              stations_latitudes,
              stations_sizes,
              stations_colors,
              stations_markers,
              stations_markercaptions,
              markers,
              facecolor='k',
              edgecolor='k',
              linewidths=.5,
              zorder=98,
                            legend_sum=legend_sum
              )
    scatter6d(bmap,
              stations_longitudes,
              stations_latitudes,
              stations_sizes,
              stations_colors,
              stations_markers,
              stations_colorcaptions,
              colors,
              edgecolor='none',
              linewidths=0,
              zorder=99,
                            legend_sum=legend_sum
              )
    times, = addons.core.search(self, fields=['start_date'], levels=['networks','stations'])
    if len(times)>0 and stationintitle:
        titletext+= '\n%s stations' % (len(times))
        if dateintitle:
            titletext+= '(%s to %s)' % (str(min(times))[:10], str(max(times))[:10])

    return titletext
