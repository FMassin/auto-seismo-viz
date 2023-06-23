#!/usr/bin/env python
from obspy.geodetics.base import gps2dist_azimuth
from mpl_toolkits.basemap import Basemap
from obspy.core.event.catalog import Catalog
from obspy.core.inventory.inventory import Inventory
from obspy.imaging.source import beach
import glob,geocoder,numpy
import matplotlib.patheffects,matplotlib.gridspec,matplotlib.pyplot,matplotlib.cm
import scipy.interpolate
from addons.catalog import map_events,get, legend_title
from addons.catalog import distfilter as catalog_distfilter
from addons.inventory import map_stations
from addons.inventory import distfilter as inventory_distfilter
from addons.stream import get_velmodel_correction
gold= (1 + 5 ** 0.5) / 2.


def ipe_allen2012_hyp(r,
                      m,
                      a = 2.085,
                      b = 1.428,#.913,#1.06,
                      c = -1.402,#-1.107,#-0.0010,
                      d = 0.078,#.813,#-3.37,
                      s = 1,
                      m1=-0.209,
                      m2=2.042):
    rm = m1+m2*numpy.exp(m-5)
    I = a + b*m + c*numpy.log(numpy.sqrt(r**2 + rm**2))+s
    for i,ri in enumerate(r):
        try:
            for j,rj in enumerate(ri):
                if rj<50:
                    I[i,j] = a + b*m[i,j] + c*numpy.log(numpy.sqrt(r[i,j]**2 + rm[i,j]**2))+d*numpy.log(r[i,j]/50)+s
        except:
            if ri<50:
                I[i] = a + b*m[i] + c*numpy.log(numpy.sqrt(r[i]**2 + rm[i]**2))+d*numpy.log(r[i]/50)+s
    return I

def lineMagnitude(x1, y1, x2, y2):
    lineMagnitude = numpy.sqrt(numpy.power((x2 - x1), 2)+ numpy.power((y2 - y1), 2))
    return lineMagnitude

#Calc minimum distance from a point and a line segment (i.e. consecutive vertices in a polyline).
def DistancePointLine(px, py, x1, y1, x2, y2):
    #http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/source.vba
    LineMag = lineMagnitude(x1, y1, x2, y2)

    if LineMag < 0.00000001:
        DistancePointLine = numpy.sqrt(numpy.power((px - x1), 2)+ numpy.power((py - y1), 2)) # 9999
        return DistancePointLine

    u1 = (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
    u = u1 / (LineMag * LineMag)

    if (u < 0.00001) or (u > 1):
        #// closest point does not fall within the line segment, take the shorter distance
        #// to an endpoint
        ix = lineMagnitude(px, py, x1, y1)
        iy = lineMagnitude(px, py, x2, y2)
        if ix > iy:
            DistancePointLine = iy
        else:
            DistancePointLine = ix
    else:
        # Intersecting point is on the line, use the formula
        ix = x1 + u * (x2 - x1)
        iy = y1 + u * (y2 - y1)
        DistancePointLine = lineMagnitude(px, py, ix, iy)

    return DistancePointLine

def scfinderauthor(origin, lineauthors=None):

    if lineauthors is not None and len(lineauthors):
        print('WARNING: User defined',lineauthors)
        return lineauthors
    lineauthors = ['scfinder']
    geocode = geocoder.osm([origin.latitude,
                            origin.longitude],
                           method='reverse').json
    print("==============\ngeocode")
    print(geocode)
    
    if geocode is not None and geocode['country'] in ['Schweiz/Suisse/Svizzera/Svizra','Switzerland','France']:
        print(geocode)
        if ( geocode['region'] in ['Valais/Wallis', 'Valais'] or
            geocode['county'] in ['Sierre','Haute-Savoie' ] or
            'Glarus' in geocode['address']):
            lineauthors += ['scfdalpine']
        else:
            lineauthors += ['scfdforela']
    else:
        bm = Basemap()
        if origin.depth / 1000 > 70:
            lineauthors += ['scfd85sym']
        elif bm.is_land(origin.longitude, origin.latitude):
            lineauthors += ['scfdcrust']
        else:
            lineauthors += ['scfd20asym']
    return lineauthors

def plot_eewsourcepoints(event,
                         bmap,
                         authors=['scvsmag','scvsmag2'],
                         magnitude_types=['MVS'],
                         colors={'MVS':'C2','Mfd':'C1'}):

    bmap.plot(event.preferred_origin().longitude,
              event.preferred_origin().latitude,
              '*',
              label=r'M$_{%s}$%.1f'%(event.preferred_magnitude().magnitude_type[1:],
                                     event.preferred_magnitude().mag),
              markeredgecolor='k',
              markerfacecolor='None',
              latlon=True,
              markersize=20,
              zorder=9999,
              path_effects=[matplotlib.patheffects.withStroke(linewidth=4,
                                                              foreground="white")])

    magnitudes = []
    for m in event.magnitudes:
        if (m.creation_info.author is not None and
            m.magnitude_type in magnitude_types and
            m.creation_info.author.split('@')[0] in authors ):
            magnitudes+=[m]
    if len(magnitudes)==0:
        return

    indexes = numpy.argsort([m.creation_info.creation_time for m in magnitudes])
    magnitudes = [magnitudes[i] for i in indexes]
    alphanorm=(magnitudes[0].creation_info.creation_time - event.preferred_origin().time)
    label='M$_{%s}$%.1f\n%.1fs'%(magnitudes[0].magnitude_type[1:],
                                  magnitudes[0].mag,
                                  magnitudes[0].creation_info.creation_time - event.preferred_origin().time)

    for k,m in enumerate(magnitudes):
        alpha = m.creation_info.creation_time - event.preferred_origin().time
        alpha = alpha - alphanorm
        alpha = min([1,max([0.05,1-alpha/(alphanorm*2)])])#**.3
        l,=bmap.plot(m.origin_id.get_referred_object().longitude,
                      m.origin_id.get_referred_object().latitude,
                      '*',
                      label=label,
                      markerfacecolor="None",
                      markeredgecolor=colors[m.magnitude_type],
                      latlon=True,
                      markersize=20*(m.mag**2)/(event.preferred_magnitude().mag**2),
                      alpha=alpha,
                      zorder=99-k,
                      path_effects=[matplotlib.patheffects.withStroke(linewidth=3,
                                                                      foreground=[1,1,1,alpha])])
                                                                                                       
        l.set_solid_capstyle('round')
        label=None

    bmap.ax.legend()

def plot_eewsourcelines(event,
                        bmap,
                        drawline=True,
                        authors=['scfd85sym'],
                        magnitude_types=['Mfd'],
                         colors={'MVS':'C2','Mfd':'C1'}):

    magnitudes=[]
    for m in event.magnitudes:
        if (m.creation_info.author is not None and
            m.magnitude_type in magnitude_types and
            m.creation_info.author.split('@')[0] in authors ):
            magnitudes+=[m]
    if len(magnitudes)==0:
        print('WARNING: No source line found')
        return
    indexes = numpy.argsort([m.creation_info.creation_time for m in magnitudes])
    magnitudes = [magnitudes[i] for i in indexes]
    alphanorm =  (magnitudes[0].creation_info.creation_time - event.preferred_origin().time)
    label='M$_{%s}$%.1f\n%.1fs'%(magnitudes[0].magnitude_type[1:],
                                  magnitudes[0].mag,
                                  magnitudes[0].creation_info.creation_time - event.preferred_origin().time)
    for k,m in enumerate(magnitudes):
        length=None
        for c in m.comments:
            if 'rupture-length' == str(c.resource_id).split('/')[-1]:
                length = float(c.text)
            if 'rupture-strike' == str(c.resource_id).split('/')[-1]:
                strike = float(c.text)
            if 'likelyhood' == str(c.resource_id).split('/')[-1]:
                likelyhood = float(c.text)
        
        alpha = m.creation_info.creation_time - event.preferred_origin().time
        alpha = alpha - alphanorm
        alpha = (max([0.05,1-alpha/(alphanorm*2)]))#**.3
        path_effects=[matplotlib.patheffects.withStroke(linewidth=4,
                                                      foreground=[1,1,1,alpha])]
        
        if length is not None and drawline:
            elon = m.origin_id.get_referred_object().longitude+\
            length/110.*\
            numpy.sin(numpy.deg2rad(strike))
            elat = m.origin_id.get_referred_object().latitude+\
            length/110.*\
            numpy.cos(numpy.deg2rad(strike))
            l,=bmap.plot([m.origin_id.get_referred_object().longitude,elon],
                         [m.origin_id.get_referred_object().latitude,elat],
                         '-',
                         label=label,
                         color=colors[m.magnitude_type],
                         latlon=True,
                         zorder=999-k,
                         alpha=alpha,
                        path_effects=path_effects
                        )
            l.set_solid_capstyle('round')
        label=None
        if True:
            l,=bmap.plot(m.origin_id.get_referred_object().longitude,
                      m.origin_id.get_referred_object().latitude,
                      '*',
                      label=label,
                      markerfacecolor="None",
                      markeredgecolor=colors[m.magnitude_type],
                      latlon=True,
                      markersize=20*(m.mag**2)/(event.preferred_magnitude().mag**2),
                      alpha=alpha,
                      zorder=99-k,
                      path_effects=path_effects)
            l.set_solid_capstyle('round')

    bmap.ax.legend()

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    #lon1 = numpy.deg2rad(lon1)
    #lat1 = numpy.deg2rad(lat1)
    #lon2 = numpy.deg2rad(lon2)
    #lat2 = numpy.deg2rad(lat2)
    lon1, lat1, lon2, lat2 = map(numpy.deg2rad, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = numpy.sin(dlat/2)**2 + numpy.cos(lat1) * numpy.cos(lat2) * numpy.sin(dlon/2)**2
    c = 2 * numpy.arcsin(numpy.sqrt(a))
    r = 6371000 # Radius of earth in meters. Use 3956 for miles
    return c * r

def nicemapscale(self):
    m=self
    fig=self.ax.get_figure()
    diaginch = (sum(m.ax.get_position().size**2.))**.5*(sum(fig.get_size_inches()**2.))**.5
    diagkm = haversine(min(m.boundarylons),
                               min(m.boundarylats),
                               max(m.boundarylons),
                               max(m.boundarylats))/1000.
    widthinch = m.ax.get_position().size[0]*fig.get_size_inches()[0]
    widthkm = haversine(min(m.boundarylons),
                       numpy.mean(m.boundarylats),
                       max(m.boundarylons),
                       numpy.mean(m.boundarylats))/1000.
    heightinch = m.ax.get_position().size[1]*fig.get_size_inches()[1]
    heightkm = haversine(numpy.mean(m.boundarylons),
                         min(m.boundarylats),
                         numpy.mean(m.boundarylons),
                         max(m.boundarylats))/1000.

    m.ax.set_aspect((heightkm/widthkm) / (heightinch/widthinch))

    diagdeg = ((max(m.boundarylons)-min(m.boundarylons))**2.+(max(m.boundarylats)-min(m.boundarylats))**2)**.5
    steps = numpy.asarray([ 0.001,0.0025,0.005,0.01,0.025,0.05,0.1,0.25,.5,1,2.5,5,10,25,50,100,250,500,1000])
    mainscalekm = steps[numpy.argmin((diagkm/(1*diaginch)-steps)**2)]
    mainscaledeg = mainscalekm * diagdeg/diagkm
    m.plot(numpy.mean(m.boundarylons)+[-mainscaledeg/2., mainscaledeg/2.],
           [min(m.boundarylats),min(m.boundarylats)],
           'w',
           linewidth=8,
           latlon=True,
           solid_capstyle='butt',
           alpha=.7,
           zorder=999999)
    for i in range(3):
        scalekm = steps[numpy.argmin((diagkm/(1*diaginch)-steps)**2)-i]
        scaledeg = scalekm * diagdeg/diagkm
        if mainscalekm/scalekm<6.:
            n=0
            c='k'
            while scalekm*(n+1)<=mainscalekm:
                n+=1
                m.plot(numpy.mean(m.boundarylons)-mainscaledeg/2.+[scaledeg*(n-1), scaledeg*n],
                       [min(m.boundarylats),min(m.boundarylats)],
                       c,
                       linewidth=5,
                       latlon=True,
                       solid_capstyle='butt',
                       zorder=999999)
                if c=='k':
                    c='w'
                else:
                    c='k'
    xy=m(numpy.mean(m.boundarylons), min(m.boundarylats)+mainscaledeg/12)
    m.ax.text(xy[0],xy[1],
              s='%s km'%(mainscalekm),
              va='bottom',
              ha='center',
              fontsize='xx-small',
              alpha=.7,
              path_effects=[matplotlib.patheffects.withStroke(linewidth=3,
                                                              foreground="white")],
              zorder=999999)

def search(self,
           fields=['longitude','latitude','elevation'],
           levels=['networks','stations']):

    out = list()
    for i in fields:
        out.append(list())

    if not levels:
        for i,a in enumerate(fields):
            if hasattr(self,a):
                out[i].append(getattr(self,a))
    else:
        for level0 in getattr(self,levels[0]) :
            if len(levels)>1 and hasattr(level0,levels[1]):

                for level1 in getattr(level0,levels[1]):
                    if len(levels)>2 and hasattr(level1,levels[2]):

                        for level2 in getattr(level1,levels[2]):
                            if len(levels)>3 and hasattr(level2,levels[3]):

                                for level3 in getattr(level2,levels[3]):
                                    if len(levels)>4 and hasattr(level3,levels[4]):
                                        print('WARNING: Cannot go further than level 4')

                                    elif len(levels)>4 and not hasattr(level3,levels[4]):
                                        for i,a in enumerate(fields):
                                            out[i].append(numpy.nan)
                                    else:
                                        for i,a in enumerate(fields):
                                            if hasattr(level3,a):
                                                out[i].append(getattr(level3,a))
                                            else:
                                                out[i].append(numpy.nan)

                            elif len(levels)>3 and not hasattr(level2,levels[3]):
                                for i,a in enumerate(fields):
                                    out[i].append(numpy.nan)
                            else:
                                for i,a in enumerate(fields):
                                    if hasattr(level2,a):
                                        out[i].append(getattr(level2,a))
                                    else:
                                        out[i].append(numpy.nan)

                    elif len(levels)>2 and not hasattr(level1,levels[2]):
                        for i,a in enumerate(fields):
                            out[i].append(numpy.nan)
                    else:
                        for i,a in enumerate(fields):
                            if hasattr(level1,a):
                                out[i].append(getattr(level1,a))
                            else:
                                out[i].append(numpy.nan)

            elif len(levels)>1 and not hasattr(level0,levels[1]):
                for i,a in enumerate(fields):
                    out[i].append(numpy.nan)
            else:
                for i,a in enumerate(fields):
                    if hasattr(level0,a):
                        out[i].append(getattr(level0,a))
                    else:
                        out[i].append(numpy.nan)

    return out

def get_aspectratio(catalog=Catalog(),
                    inventory=Inventory([],''),
                    aspectratio=gold,
                    figsize=4):

    lons, lats = search(inventory, fields=['longitude','latitude'], levels=['networks','stations'])
    
    for e in catalog:
        for o in e.origins:
            lats.append(o.latitude)
            lons.append(o.longitude)

    if max(lons)-min(lons) < max(lats)-min(lats):

        if isinstance(figsize,list):
            figsize=figsize[0]
        figsize=(figsize, figsize*aspectratio)
        aspectratio = 1/aspectratio
    else:
        if isinstance(figsize,list):
            figsize=figsize[1]
        figsize=(figsize*aspectratio,figsize)
    return aspectratio, lons, lats, figsize

def nicemap(catalog=Catalog(),
            inventory=Inventory([],''),
            aspectratio=gold,
            f=None,
            ax=None,
            alpha=1.,
            xpixels=400,
            resolution='l',
            fontsize='small',
            figsize=[6*gold,6],
            labels=[1,0,0,1],
            shift=0.,
            dpi=300,
            reference=True,
            dark=False,
            mapbounds=None,
            arcgis=True,
            epsg = 3857, #anywhere?
            #epsg=5520, # Switzerland?
            #epsg=4326 # CAM?
            scale=True,
            server='http://server.arcgisonline.com/ArcGIS',
            service='Ocean/World_Ocean_Base',
            shapefile='/Users/fred/Documents/Data/misc/gem-global-active-faults/shapefile/gem_active_faults_harmonized',
            shapecolors='slip_type',
            default_encoding='iso-8859-15'):

    if mapbounds:
        lons = mapbounds[0]
        lats = mapbounds[1]
    else:
        aspectratio, lons, lats, figsize = get_aspectratio(catalog=catalog,
                        inventory=inventory,
                        aspectratio=aspectratio,
                        figsize=figsize)
        
        y = lats.copy()
        x = lons.copy()

        if max(lons)-min(lons) > aspectratio*(max(lats)-min(lats)):
            # Too large
            m = sum(y)/len(y)
            d = max(x)-min(x)
            lats.append(m + d/(2*aspectratio))
            lats.append(m - d/(2*aspectratio))
        else:
            # Too high
            m = sum(x)/len(x)
            d = max(y)-min(y)
            lons.append(m + d*aspectratio/2)
            lons.append(m - d*aspectratio/2)

    if ax:
        f = ax.get_figure()
    elif f:
        ax = f.add_subplot(111)
    else:
        f = matplotlib.pyplot.figure(figsize=figsize)
        ax = f.add_subplot(111)

    projection=None#'merc'
    centerlon = (max(lons)+min(lons))/2
    centerlat = (max(lats)+min(lats))/2
    if centerlon>-100 and centerlon<-70 and centerlat>0 and centerlat<30:
        print('Central america')
        #epsg=4326 # CAM
        projection=None

    bmap = Basemap(llcrnrlon=min(lons)-(max(lons)-min(lons))*.05,
                   llcrnrlat=min(lats)-(max(lats)-min(lats))*.05,
                   urcrnrlon=min([179.99,max(lons)+(max(lons)-min(lons))*.05]),
                   urcrnrlat=max(lats)+(max(lats)-min(lats))*.05,
                   epsg=epsg,#5520,#4326,
                   projection=projection,
                   resolution=resolution,
                   ax=ax,
                   #lon_0=min(lons)-(max(lons)-min(lons))*.1,
                   #o_lat_p=90,
                   #o_lon_p=45
                   )
    if arcgis:
        if not isinstance(server,list):
            server = [server]
            service = [service]
        for i,s in enumerate(server):
            try:
                im2 = bmap.arcgisimage(service=service[i],
                                    server=server[i],
                                    xpixels=xpixels,
                                    dpi=dpi,
                                    zorder=-999999,
                                        verbose= True
                                    )
            except:
                im2 = bmap.arcgisimage(service=service[i],
                                    server=server[i],
                                    dpi=dpi,
                                    zorder=-999999,
                                        verbose= True
                                    )
            im2.set_alpha(alpha)
        if True:
            if dark:
                im3 = bmap.arcgisimage(service='Elevation/World_Hillshade_Dark',
                                   xpixels=xpixels,
                                   dpi=dpi,
                                   zorder=-99999
                                   )
                #im3.set_alpha(.9)
            else:
                try:
                    im3 = bmap.arcgisimage(service='Elevation/World_Hillshade',
                                    xpixels=xpixels,
                                    dpi=dpi,
                                    zorder=-99999
                                    )
                except:
                    im3 = bmap.arcgisimage(service='Elevation/World_Hillshade',
                                    dpi=dpi,
                                    zorder=-99999
                                    )

            data=im3.get_array()
            data[:,:,3] = 1-(data[:,:,0]*data[:,:,1]*data[:,:,2])
            if numpy.nanmax(data[:,:,:2])>1:
                data[:,:,3] = (1-(data[:,:,0]/255*data[:,:,1]/255*data[:,:,2]/255))*255
            im3.set_array(data)
        if reference:
            try:
                im1 = bmap.arcgisimage(service='Reference/World_Boundaries_and_Places_Alternate',
                                    xpixels=int(xpixels/2),
                                    dpi=dpi,
                                    zorder=-9999)
            except:
                im1 = bmap.arcgisimage(service='Reference/World_Boundaries_and_Places_Alternate',
                                    dpi=dpi,
                                    zorder=-9999)
            im1.set_alpha(.8)

    
    steps=numpy.asarray([ 0.001,0.0025,0.005,0.01,0.025,0.05,0.1,0.25,.5,1,2.5,5,10,25,50,100,250,500,1000])
    if aspectratio>1:
        step = steps[numpy.argmin(abs((max(lons)-min(lons))/5 - steps))]
    else:
        step = steps[numpy.argmin(abs((max(lats)-min(lats))/5 - steps))]


    bmap.resolution=resolution
    bmap.drawmeridians(meridians=numpy.arange(int(min(lons))-1,int(max(lons))+2,step),
                         labels=labels,
                         linewidth=.25,
                         color = 'w',
                         rotation=20,
                         fontsize=fontsize)
    bmap.drawparallels(circles=numpy.arange(int(min(lats))-1,int(max(lats))+2,step),
                         labels=labels,
                         linewidth=.25,
                         rotation=20,
                         color = 'w',
                         fontsize=fontsize)
    try:
        l = bmap.drawcoastlines(linewidth=1.,
                              color = 'w')
        l.set_alpha(alpha)
        l=bmap.drawcoastlines(linewidth=.5)
        l.set_alpha(alpha)
    except:
        print('WARNING: No coastline')
    #l=bmap.drawstates(linewidth=.5,
    #                  color = 'w')
    #l.set_alpha(alpha)

    l=bmap.drawcountries(linewidth=1.5,
                         color='w')
    l.set_alpha(alpha)
    l=bmap.drawcountries(linewidth=1)
    l.set_alpha(alpha)
    
    #l=bmap.drawstates(linewidth=.5,
    #                  color='w')
    #l.set_alpha(alpha)
    #l=bmap.drawstates(linewidth=.25)
    #l.set_alpha(alpha)
    
    #shapefile='/Users/fmassin/Documents/Data/misc/tectonicplates/PB2002_plates'
    import matplotlib._color_data as mcd
    if glob.glob('%s*'%shapefile):
        bmap.readshapefile(shapefile,
                           name=shapefile.split('/')[-1],
                           color='black',
                           default_encoding=default_encoding,
                           drawbounds = False)
        infos=bmap.__dict__.get('%s_info'%shapefile.split('/')[-1])
        shapes=bmap.__dict__.get(shapefile.split('/')[-1])
        slip_types = list(set([info[shapecolors] for info in infos]))
        colornames = [name for name in mcd.CSS4_COLORS if "xkcd:" + name in mcd.XKCD_COLORS]
        caption=[]
        for info, shape in zip(infos,shapes):
            nameindex = slip_types.index(info[shapecolors])
            color = mcd.XKCD_COLORS["xkcd:"+colornames[nameindex]]
            if 'ubduction' not in info['slip_type']:
                color = 'k'
            x, y = zip(*shape)
            mask = [i for i,x in enumerate(x)if x>=min(lons)-2
                                             and x<=max(lons)+2
                                             and y[i]>=min(lats)-2
                                             and y[i]<=max(lats)+2]
            if not len(mask):
                continue
            caption+=['%s : %s'%(colornames[nameindex],info['slip_type'])]
            pe=[matplotlib.patheffects.withStroke(linewidth=2.5,
                                                  alpha=.1,
                                                  foreground=color),
                        matplotlib.patheffects.withStroke(linewidth=1.8,
                                                  alpha=.1,
                                                  foreground=color),
                        matplotlib.patheffects.withStroke(linewidth=1.2,
                                                  alpha=.1,
                                                  foreground=color)]

            l,=bmap.plot([x[i] for i in mask],
                      [y[i] for i in mask],
                      marker=None,
                      color='white',
                      linewidth=0,
                      path_effects=pe)
            l.set_solid_capstyle('round')
            l,=bmap.plot([x[i] for i in mask],
                      [y[i] for i in mask],
                      marker=None,
                      linewidth=.2,
                      color='white',
                      zorder=3)
            l.set_solid_capstyle('round')
        #print(', '.join(list(set(caption)))+'.')

    if scale:
        nicemapscale(bmap)

    return f,ax,bmap

def map_all(self=None,
            others=[],
            xpixels=900,
            resolution='h',
            fontsize=8,
            vmin=None,
            vmax=None,
            rank=1,
            legend_sum=True,
            alpha=1,
            aspectratio=gold,
            markers='none',#networks',
            colors='instruments',
            eventcolors='depth',
            label=False,
            showlegend=True,
            colorbar=True,
            title=True,
            labels=[1,0,0,1],
            titleaddons='',
            stationgroups = {'NU':['NU.'],
                             #'CAM':['G.','GI.','SV.','HN.','NU.','OV.','TC.'],
                            'ECOS':['CH','C4']},
            titlereplacement=None,
            insetfilter=None,
            fig=None,
            ax=None,
            reference=None,
            dark=False,
            mapbounds=None,
            prospective_inventory=None,
            latencies=None,
            extracatalog=None,
            extramarkercatalog='1',
            extranamecatalog='',
            markersize=30,
            forquiver=None,
            arcgis=True,
            titletext_alpha=1,
            fp=False,
            scale=True,
            epsg=3857,
            server='http://server.arcgisonline.com/ArcGIS',
            service='Ocean/World_Ocean_Base',
            force=None,
            loc=1,
            **kwargs):
    #Inventory.plot(self, projection='global', resolution='l',continent_fill_color='0.9', water_fill_color='1.0', marker="v",s ize=15**2, label=True, color='#b15928', color_per_network=False, colormap="Paired", legend="upper left", time=None, show=True, outfile=None, method=None, fig=None, **kwargs)

    inventory=Inventory([],'')
    catalog=Catalog()
    others.append(self)

    for other in others:
        if isinstance(other, Catalog):
            catalog += other
        elif isinstance(other, Inventory):
            inventory += other


    fig, ax, bmap = nicemap(catalog=catalog,
                            inventory=inventory,
                            alpha=alpha,
                            aspectratio=aspectratio,
                            xpixels=xpixels,
                            dark=dark,
                            resolution=resolution,
                            fontsize=fontsize,
                            f=fig,
                            ax=ax,
                            epsg=epsg,
                            labels=labels,
                            shift=1/4.,
                            reference=reference,
                            mapbounds=mapbounds,
                            arcgis=arcgis,
                            scale=scale,
                            server=server,
                            service=service)
    titletext = ''
    if False:
        titletext = map_events(catalog,
                                    bmap=bmap,
                                    fig=fig,
                                    eqcolorfield = eventcolors,
                                    titletext=titleaddons,
                                    colorbar=colorbar,
                                    fontsize=fontsize,
                                    prospective_inventory=prospective_inventory,
                                    latencies=latencies,
                                    extra=extracatalog,
                                    extramarker=extramarkercatalog,
                                    extraname=extranamecatalog,
                                    fp=fp,
                                    vmin=vmin,
                                    vmax=vmax,
                                    force=force)

    titletext = map_stations(inventory,
                                       bmap=bmap,
                                       fig=fig,
                                       colors=colors,
                                       markersize=markersize,
                                       rank=rank,
                                       legend_sum=legend_sum,
                                       markers=markers,
                                       titletext=titletext,
                                       fontsize=fontsize,
                                       stationgroups=stationgroups)

    if title:
        #sticker(titletext, bmap.ax, x=0, y=1, ha='left', va='bottom')
        if titlereplacement is not None:
            titletext = titlereplacement
        bmap.ax.annotate(titletext,
                         xy=(0, 1),
                         xycoords='axes fraction',
                         horizontalalignment='left',
                         verticalalignment='bottom',
                         alpha=titletext_alpha)
    if showlegend:
        fig.legend = bmap.ax.legend(prop={'size':fontsize}, loc=loc)#,ncol=2


    if forquiver:
        if not 'latlon' in forquiver:
            forquiver['latlon']=True
        if not 'scale_units' in forquiver:
            forquiver['scale_units']='inches'
        if not 'color' in forquiver:
            forquiver['color']='gray'

        bmap.quiver(**forquiver)

    if insetfilter:

        insetcatalog = catalog_distfilter(catalog,
                                                **insetfilter)
        insetinventory = inventory_distfilter(inventory,
                                                    **insetfilter)

        mainaspectratio, lons, lats, figsize = get_aspectratio(catalog=catalog,
                                                           inventory=inventory,
                                                           aspectratio=aspectratio)
        insetaspectratio, lons, lats, figsize = get_aspectratio(catalog=insetcatalog,
                                                           inventory=insetinventory,
                                                           aspectratio=aspectratio)
        nl=2
        nc=2
        pos=1
        if mainaspectratio>1:
            if insetaspectratio>1:
                #print('HORIZONTAL & horizontal')
                width_ratios=[1, aspectratio-1]
                height_ratios = [aspectratio-1, 1]
            else:
                #print('HORIZONTAL & vertical')
                width_ratios=[1*.86, aspectratio-1]
                height_ratios = [1]
                nl=1
        else:
            if insetaspectratio>1:
                #print('VERTICAL & horizontal')
                nc=1
                pos=0
                width_ratios=[1,]
                height_ratios = [insetaspectratio-1,1*.85]
            else:
                #print('VERTICAL & vertical')
                width_ratios=[1, insetaspectratio]
                height_ratios = [1/insetaspectratio-1,1]

        if colorbar and len(catalog.events)>0 :
            nc+=1
            width_ratios.append(.05)
            #print('cb qspace')


        gs = matplotlib.gridspec.GridSpec(nl, nc,
                                          width_ratios=width_ratios,
                                          height_ratios=height_ratios)
        axinset = fig.add_subplot(gs[pos])
        fig, fig.axinset, fig.bmapinset = map_all(others=[insetinventory,
                                                          ],#insetcatalog],
                                                 label=label,
                                                 xpixels=xpixels/3,
                                                 resolution=resolution,
                                                 fontsize=5,
                                                 labels=[0,0,0,0],#[0,1,1,0],
                                                 markers=markers,
                                                 colors=colors,
                                                 alpha=alpha,
                                                 showlegend=False,
                                                 aspectratio=aspectratio,
                                                 ax=axinset,
                                                 title=False,
                                                 colorbar=False,
                                                  scale=False,
                                                  arcgis=arcgis,
                                                 **kwargs)
        insetbounds = [fig.bmapinset.boundarylons,
                       fig.bmapinset.boundarylats]
        insetbounds[0].append(fig.bmapinset.boundarylons[0])
        insetbounds[1].append(fig.bmapinset.boundarylats[0])
        bmap.plot(insetbounds[0],
                  insetbounds[1],
                  'red')
        for spine in fig.bmapinset.ax.spines.values():
            spine.set_edgecolor('red')

    fig.bmap=bmap
    return fig, ax, bmap

def fixtissot(polygon,bmap,rad,lon):
    xy = polygon.get_xy()
    newxy=[]
    for i in range(len(xy)):
        test=bmap(xy[i][0],xy[i][1],inverse=True)
        if False:#test[0]<0 :
            newxy+=[list(bmap(180,test[1]))]
            continue
        newxy+=[xy[i]]
    polygon.set_xy(newxy)

def plot_focmech(event,lineauthors,ax,color='C1'):
    pos = ax.get_position()
    axins = ax.figure.add_axes([pos.x1-pos.width/6,pos.y1-pos.height/6,pos.width/3,pos.height/3],
                         facecolor='None',
                                   zorder=999,
                                   frame_on=False)
    axins.set_axis_off()
    axins.set(xlim=(-50, 50), ylim=(-50, 50))
    norm = matplotlib.pyplot.Normalize(-1., 1.)
    cmap = matplotlib.cm.get_cmap('bwr')
    opts= {'xy':(0, 0), 'width':50,
           'facecolor':'0.5',
           'edgecolor':'None'
           #'alpha':0.5
          }
    focmec = event.preferred_focal_mechanism()
    if focmec is None:
        #print(event)
        return
    fm = [focmec.nodal_planes.nodal_plane_1.strike,
              focmec.nodal_planes.nodal_plane_1.dip,
              focmec.nodal_planes.nodal_plane_1.rake]
    bball = obspy.imaging.source.beach(fm,**opts)
    axins.add_collection(bball)

    opts['nofill']=True
    opts['edgecolor']='k'
    opts['linewidth']=1
    if hasattr(focmec.moment_tensor,'tensor'):
        mtensor = focmec.moment_tensor.tensor
        mt = [mtensor.m_rr, mtensor.m_tt, mtensor.m_pp,
              mtensor.m_rt, mtensor.m_rp, mtensor.m_tp]
        bball = obspy.imaging.source.beach(mt,**opts)
    else:
        bball = obspy.imaging.source.beach(fm,**opts)
    axins.add_collection(bball)

    opts['edgecolor']=color
    withStroke=matplotlib.patheffects.withStroke
    magnitudes = []
    for m in event.magnitudes:
        if not str(m.creation_info.author).split('@')[0] in lineauthors:
            continue
        if not hasattr(m,'comments'):
            continue
        for c  in m.comments:
            if 'strike' in str(c.resource_id):
                magnitudes+=[m]
                break
        
    creation_delays=[m.creation_info.creation_time-event.preferred_origin().time for m in magnitudes]
    for m in magnitudes:
        for c  in m.comments:
            if 'strike' in str(c.resource_id):
                f = [int(c.text)+90,0,0]
                if f[0]<0:
                    f+=180
                if f[0]>360:
                    f-=180
                delta = m.creation_info.creation_time - event.preferred_origin().time
                alpha = (delta - min(creation_delays)) / (max(creation_delays) - min(creation_delays))
                #print(event.preferred_origin().time, m.creation_info.creation_time, event.magnitudes[0].creation_info.creation_time)
                #print(delta, alpha, 1-alpha**.15)
                opts['alpha'] = 1-alpha**.2
                opts['linewidth'] = (1-alpha**.2)
                
                bball = beach(f,**opts)
                axins.add_collection(bball)
                x = 35*numpy.cos(numpy.deg2rad(f[0]))
                y = -35*numpy.sin(numpy.deg2rad(f[0]))
                axins.text(x, y,'%d s'%delta,
                           fontsize='x-small',
                           zorder = -delta,
                           horizontalalignment='center',
                           verticalalignment='center',
                           path_effects=[withStroke(linewidth=4,
                                                    alpha = opts['alpha'],
                                                    foreground=color),
                                          withStroke(linewidth=2,
                                                     foreground="w")])
        
    axins.set(xlim=(-50, 50), ylim=(-50, 50))
    
def getradius(origin,inventory,maxdist=1.5):

    minlatitude=origin.latitude
    maxlatitude=origin.latitude
    minlongitude=origin.longitude
    maxlongitude=origin.longitude
    for n in inventory:
        for s in n:
            if minlatitude > s.latitude and s.latitude>origin.latitude-maxdist:
                minlatitude = s.latitude
            if maxlatitude < s.latitude and s.latitude<origin.latitude+maxdist:
                maxlatitude = s.latitude
            if minlongitude > s.longitude and s.longitude>origin.longitude-maxdist:
                minlongitude = s.longitude
            if maxlongitude < s.longitude and s.longitude<origin.longitude+maxdist:
                maxlongitude = s.longitude
    
    return max([ gps2dist_azimuth(origin.latitude,origin.longitude,
                                    origin.latitude,minlongitude)[0],
                   gps2dist_azimuth(origin.latitude,origin.longitude,
                                    minlatitude,origin.longitude)[0],
                   gps2dist_azimuth(origin.latitude,origin.longitude,
                                    origin.latitude,maxlongitude)[0],
                   gps2dist_azimuth(origin.latitude,origin.longitude,
                                    maxlatitude,origin.longitude)[0] ])

def eewmap(data,
           vs=3.7,
           drawline=True,
           title=True,
           legend=True,
           lineauthors=None,
           stationgroups=None,
           mapbounds=True,
           **kwargs):

    origin = data['event'].preferred_origin()
    magnitude = data['event'].preferred_magnitude()
    lineauthors = scfinderauthor(origin, lineauthors=lineauthors)

    radius = getradius(origin,data['inventory'])

    possibledelaysteps=[1,2,5,10,20,50,100]
    delaystep=numpy.argmin( [ abs(s-radius/10/vs/1000) for s in possibledelaysteps])
    delaystep=possibledelaysteps[delaystep]

    longitude_radius = 0
    km=0
    while km<=radius:
        km,az,baz=gps2dist_azimuth(origin.latitude,
                                    origin.longitude,
                                    origin.latitude,
                                    origin.longitude+longitude_radius)
        longitude_radius+=0.001
        
    latitude_radius = 0
    km=0
    while km<=radius:
        km,az,baz=gps2dist_azimuth(origin.latitude,
                                    origin.longitude,
                                    origin.latitude+latitude_radius,
                                    origin.longitude)
        latitude_radius+=0.001
    if mapbounds:
        mapbounds=[[origin.longitude-longitude_radius,
                    origin.longitude+longitude_radius],
                    [origin.latitude-latitude_radius,
                    origin.latitude+latitude_radius]]
    else:
        mapbounds=None
    fig,ax,bmap = map_all(data['inventory'],
                          mapbounds=mapbounds,
                          stationgroups=stationgroups,
                          **kwargs)
        
    plot_eewsourcepoints(event=data['event'],
                         bmap=fig.bmap)
    plot_eewsourcelines(event=data['event'],
                        authors=lineauthors,
                        drawline=drawline,
                        bmap=fig.bmap)
    
    magnitudes = []
    for m in data['event'].magnitudes:
        if (m.creation_info.author is not None and
                m.magnitude_type in ['Mfd', 'MVS'] and
                m.creation_info.author.split('@')[0] in lineauthors + ['scvsmag']):
            magnitudes += [m]
    indexes = numpy.argsort([m.creation_info.creation_time for m in magnitudes])
    magnitudes = [magnitudes[i] for i in indexes]
    
    deglatatlon = haversine(origin.longitude,
                                   origin.latitude-.5,
                                   origin.longitude,
                                   origin.latitude+.5)/1000
    deglonatlat = haversine(origin.longitude-.5,
                                   origin.latitude,
                                   origin.longitude+.5,
                                   origin.latitude)/1000

    dblind = (magnitudes[0].creation_info.creation_time - origin.time) * vs
    if dblind > origin.depth / 1000:
        dblind = (dblind ** 2 - (origin.depth / 1000) ** 2) ** .5
    else:
        dblind=0
    #print(magnitudes[0].creation_info.creation_time - origin.time,dblind,deglatatlon)
    
    modelresidual,model=get_velmodel_correction(data['event'],data['inventory'],radius/1000)
    triald=list(numpy.arange(((radius*2/1000)**2*2)**.5/2))
    modtts=[]
    for epr in triald:
        arrivals = model.get_travel_times(origin.depth / 1000,
                                          distance_in_degree=epr/111,
                                          phase_list=['tts'],
                                          receiver_depth_in_km=0)
        modtts+=[modelresidual['s']+numpy.nanmin([numpy.nan]+[a.time for a in arrivals])]
    fmodel = scipy.interpolate.interp1d([0]+modtts,[0]+triald,fill_value="extrapolate")
    #print(modtts, magnitudes[0].creation_info.creation_time - origin.time)
    dblind = fmodel(max([0.01, magnitudes[0].creation_info.creation_time - origin.time ]))
        
    polygon = fig.bmap.tissot(origin.longitude,
                    origin.latitude,
                    dblind / deglatatlon,
                    64,
                    fill=True,
                    color='k',
                    alpha=0.1,
                    edgecolor='k',
                    label='Late alert\n%d km'%dblind,
                    path_effects=[matplotlib.patheffects.withStroke(linewidth=4,
                                                                    foreground="white")])
    
    fixtissot(polygon,fig.bmap,dblind,origin.longitude)
    
    mmis=['I','II','III','IV','V','VI','VII','IIX','IIX','X','XI','XII']
    mmii=-1
    d=9999
    triald=numpy.arange(((radius*5/1000)**2*2)**.5/2)
    mmi = ipe_allen2012_hyp((triald**2+(origin.depth/1000)**2)**.5,
                                                triald*0+magnitude.mag,
                                                s=-1)
    #print(mmi,triald)
    f = scipy.interpolate.interp1d(mmi,triald)
    minmaxmmi=[max([1,int(min(mmi))+1]), min([9,1+int(max(mmi))])]
    mmi = numpy.arange(*minmaxmmi)
    mmis = numpy.asarray([mmis[i-1] for i in mmi])
    #print(mmi)
    triald = f(mmi)
    refmmi = mmis[triald>1]#(dblind/4)]
    d = triald[triald>1]
        
    withStroke=matplotlib.patheffects.withStroke
    for t in range(len(d)):
        label=None
        if t==0:
            label='Final intensity'#+' '+' to '.join(list(set([refmmi[0],refmmi[-1]])))
        polygon=fig.bmap.tissot(origin.longitude,
                        origin.latitude,
                        d[t] / deglatatlon,
                        64,
                        fill=False,
                        edgecolor='k',
                        label=label,
                        linewidth=2,
                        path_effects=[withStroke(linewidth=4,
                                                 foreground="white")])
        fixtissot(polygon,fig.bmap,d[t],origin.longitude)
        
        if t>0:
            pass#continue
        if t==0:
            label='EEW intensity'
        #print(refmmi[t], d[t], deglatatlon)
        polygon=fig.bmap.tissot(magnitudes[0].origin_id.get_referred_object().longitude,
                        magnitudes[0].origin_id.get_referred_object().latitude,
                        d[t] / deglatatlon,
                        64,
                        fill=False,
                        edgecolor={'MVS':'C2','Mfd':'C1'}[magnitudes[0].magnitude_type],
                        linestyle='--',
                        linewidth=1,
                        label=label,#'EEW MMI '+refmmi[t],
                        path_effects=[withStroke(linewidth=4,
                                                 foreground="w")])
        fixtissot(polygon,fig.bmap,d[t],magnitudes[0].origin_id.get_referred_object().longitude)
        x, y = fig.bmap(origin.longitude-d[t]/deglonatlat*.7,
                        origin.latitude+d[t]/deglatatlon*.7)
        fig.bmap.ax.annotate('%s\n$^{Intensity}$'%refmmi[t],
                             xy=(x, y),
                             xycoords='data',
                             xytext=(x, y),
                             textcoords='data',
                         horizontalalignment='center',
                         verticalalignment='center',
                         rotation=45,
                         fontweight='bold',
                         fontsize='x-small',
                         color='k',
                         clip_on=True,
                        path_effects=[withStroke(linewidth=2,
                                                 foreground="w")])

    
    if delaystep is not None and delaystep>0:
        delay=0
        while dblind<radius/1000:#*1.5:
            delay += delaystep
            dblind = fmodel(max([0.01, magnitudes[0].creation_info.creation_time - origin.time + delay]))
            #dblind = (magnitudes[0].creation_info.creation_time - origin.time + delay) * vs
            #if dblind > origin.depth / 1000:
                #dblind = (dblind ** 2 - (origin.depth / 1000) ** 2) ** .5
            #print(delay, dblind, deglatatlon)
            polygon=fig.bmap.tissot(origin.longitude,
                            origin.latitude,
                            dblind / deglatatlon,
                            64,
                            fill=False,
                            edgecolor='k',
                            lw=0.5,
                            path_effects=[matplotlib.patheffects.withStroke(linewidth=1.5,
                                                                            foreground="white")])
            fixtissot(polygon,fig.bmap,dblind,origin.longitude)
            x, y = fig.bmap(origin.longitude+dblind/deglonatlat*0.7,
                            origin.latitude+dblind/deglatatlon*0.7)
            fig.bmap.ax.annotate('%ds\n$^{EEW}$'%delay,
                                 xy=(x, y),
                                 xycoords='data',
                                 xytext=(x, y),
                         clip_on=True,
                                 textcoords='data',
                             horizontalalignment='center',
                             verticalalignment='center',
                     rotation=-45,
                     fontsize='x-small',
                     color='k',
                    path_effects=[withStroke(linewidth=2,
                                             foreground="w")])

    if legend:
        fig.bmap.ax.legend(ncol=1,
                           title=legend_title(data['event'],['Mfd', 'MVS']), 
                           title_fontproperties={'size':'xx-small',
                                                 'weight':'bold'}, 
                           bbox_to_anchor=(1,0), loc="lower left",
                           prop={'size': 'xx-small'})
    else:
        fig.bmap.ax.get_legend().remove()
    descs = ', '.join([desc.text for desc in data['event'].event_descriptions if 'region' in desc.type])
    if title:
        t = '%s\nM$_{%s}$%.1f, %s, %.1fkm deep'
        tt=(str(origin.time)[:19],
            magnitude.magnitude_type[1:],
            magnitude.mag,
            descs,
            origin.depth / 1000.)
        fig.bmap.ax.set_title(t % tt,
                              loc='left',
                              fontsize='small')
    m='EEW performance map:'
    m+=' The authoritative location (black star) is compared with EEW locations.'
    m+=' If available, VS locations are represented with red stars and FinDer with red lines.'
    m+=' The brightest symbols are the first solutions (see delay in legend).'
    if delaystep is not None and delaystep>0:
        m+=' The blind zone and EEW delays every %ds are shown with black circles considering the authoritative location.'%delaystep
    m+=' The blind zone can be compared to MMI contours considering the authoritative location (blue circles).'
    m+=' The smallest MMI contour is also represented with the first EEW solution (dash blue circle).'
    #print(m)
    
    plot_focmech(data['event'],lineauthors,fig.bmap.ax)
    
    return fig
