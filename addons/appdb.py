#!/usr/bin/env python
import sqlite3
from contextlib import closing
from numpy import log,sqrt,exp
import datetime, geopy.distance
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.pyplot import get_cmap, Rectangle
import roman
from matplotlib.ticker import FuncFormatter
from obspy.imaging.cm import pqlx
from matplotlib.ticker import EngFormatter

REPORT_REQUEST = 'SELECT eventnotif.userid, eventnotif.estmdintensity, eventinfo.magnitude, \
                        eventnotif.swavearrival, eventnotif.alertsite, \
                        intensityreports.lon,  intensityreports.lat, \
                        intensityreports.lon, intensityreports.lat,  \
                        eventinfo.longitude, eventinfo.latitude, eventinfo.depth,\
                        intensityreports.intensity, eventinfo.origintime \
                    FROM intensityreports, eventinfo, eventnotif \
                    WHERE eventinfo.eventid = "%s" \
                    AND intensityreports.eventid = eventinfo.eventid \
                    AND eventnotif.eventid = eventinfo.eventid \
                    AND eventnotif.userid = intensityreports.userid \
                    AND intensityreports.intensity >= %d \
                    AND intensityreports.intensity <= %d \
                    AND intensityreports.lat > -99 \
                    AND eventnotif.swavearrival > %s \
                    AND eventnotif.swavearrival > -99 \
                    AND eventnotif.swavearrival < 99 \
                    AND eventnotif.updateno = 0 '

NOTIFICATION_REQUEST = 'SELECT eventnotif.userid, eventnotif.estmdintensity, eventinfo.magnitude, \
                                    eventnotif.swavearrival, eventnotif.alertsite, \
                                    eventnotif.userlon, eventnotif.userlat, \
                                    eventnotif.userlon, eventnotif.userlat,  \
                                    eventinfo.longitude, eventinfo.latitude, eventinfo.depth, \
                                    eventnotif.estmdintensity, eventinfo.origintime, \
                                    eventnotif.delay \
                                    FROM eventnotif, eventinfo \
                                    WHERE eventinfo.eventid = "%s" \
                                    AND eventnotif.eventid = eventinfo.eventid \
                                    AND eventnotif.alertsite = 1 \
                                    AND eventnotif.estmdintensity >= %d \
                                    AND eventnotif.estmdintensity <= %d \
                                    AND eventnotif.swavearrival >  -99 \
                                    AND eventnotif.swavearrival <  99 \
                                    AND eventnotif.swavearrival >  %s \
                                    AND ( eventnotif.userlat > -99 OR \
                                        eventnotif.userlatpoi > -99 ) '


def myroman(n,t):
    return roman.toRoman(int(n))

# BASIC DB INTERFACE
def exec(req = "SELECT eventid, magnitude FROM eventinfo WHERE magnitude>6",
         db = None):

    with closing(sqlite3.connect(db)) as connection:
        with closing(connection.cursor()) as cursor:
            out = cursor.execute(req).fetchall()
    return out

def get_intensity(r, m,
                    a = 2.085,
                    b = 1.428,#.913,#1.06,
                    c = -1.402,#-1.107,#-0.0010,
                    d = 0.078,#.813,#-3.37,
                    s = 1,
                    m1=-0.209,
                    m2=2.042):

    rm = m1+m2*exp(m-5)

    if r <= 50:
        return a + b*m + c*log(sqrt(r**2 + rm**2))+s

    return a + b*m + c*log(sqrt(r**2 + rm**2))+d*log(r/50)+s

class ewreport():
    def __init__(self, userid, estmdintensity, magnitude, 
                 swavearrival, alertsite,
                 userlon, userlat, 
                 userlonpoi, userlatpoi,
                 longitude, latitude, depth , 
                 reportintensity=None,
                 origintime=None,
                 delay=None):
        self.userid = userid
        self.delay = delay
        self.intensity = estmdintensity
        self.swavearrival = swavearrival
        self.userlocation = (userlat, userlon)
        self.poilocation = (userlatpoi, userlonpoi)
        self.reportlocation = self.userlocation
        if alertsite == 1:
            self.reportlocation = self.userlocation
        elif alertsite == 2:
            self.reportlocation = self.poilocation
        self.eventlocation = (latitude, longitude)
        self.eventdepth = depth
        self.magnitude = magnitude
        try:
            self.epidistance = geopy.distance.geodesic(self.reportlocation, 
                                                    self.eventlocation).km
            self.eventdistance = (self.epidistance**2 + self.eventdepth**2)**.5
            #v=d/t
            #tp=d/vp
            #ts=d/vs
            ts = self.eventdistance/3.3
            tp = self.eventdistance/5.5
            tew = ts - swavearrival
            self.pwavearrival =  tp - tew
            
            self.ipeintensity = get_intensity(self.eventdistance, magnitude)
        except Exception as e:
            print(e)
            print(ts,swavearrival)
            pass#srint('No location')

        self.reportintensity = reportintensity
        self.origintime = origintime
        if origintime is not None:
            self.origintime = datetime.datetime.fromisoformat(origintime[:19])

    def __repr__(self):
        try:
            return "EMS %d @ %4g km" % (self.intensity,self.eventdistance)
        except:

            return "EMS %d (unknown event location)" % (self.intensity)
    
class intensityreport():
    def __init__(self,i,emag,ilon,ilat,elon,elat,edepth):
        self.intensity = i
        self.reportlocation = (ilon,ilat)
        self.eventlocation = (elon,elat)
        self.eventdepth = edepth
        self.magnitude = emag
        self.epidistance = float(geopy.distance.geodesic((ilat,ilon), (elat,elon)).km)
        self.eventdistance = (self.epidistance**2 + self.eventdepth**2)**.5
    def __repr__(self):
        return "EMS %d @ %4g km" % (self.intensity,self.eventdistance)
    

def get_app_data(db, event,
                   maxdistance=500,
                   maxipeerr=2, 
                    minmax=(1,8,-20)):
    
    tmp = (event, minmax[0], minmax[1], minmax[2])
    app_data = []
    for req in [REPORT_REQUEST, NOTIFICATION_REQUEST]:
        data = exec(req % tmp, db=db) 
        data = [ewreport(*d) for d in data]
        data = [ d for index,d in enumerate(data) if d.eventdistance<maxdistance and ( d.reportintensity is None or d.reportintensity<maxipeerr+get_intensity(d.eventdistance,d.magnitude) ) ]
        app_data += [data]

    return app_data


def plot_event_data(app_data, ax_map,
                    cmap='turbo',#'rainbow',#'viridis'
                   ):
    
    if app_data is None:
        return
    
    # Plot notifications
    data = app_data[1]
    longitudes = [ d.reportlocation[1] for d in data ]
    latitudes = [ d.reportlocation[0] for d in data ]
    delays = [ d.delay for d in data ]

    lolamap = ax_map(longitudes, latitudes) 
    ax_map.scatter(*lolamap,
                    facecolors='w', 
                    edgecolors='w',
                    s=10,
                    marker='o',
                    linewidths=3,
                    zorder=68,
                    #latlon=True,
                    #transform=Geodetic()
                    )
    ax_map.scatter(*lolamap,
                    facecolors='w', 
                    edgecolors='k',
                    s=10,
                    marker='o',
                    linewidths=0.5,
                    zorder=68,
                    #latlon=True,
                    #transform=Geodetic()
                    )    
    ax_map.scatter(*lolamap,
                    s=10,
                    marker='o',
                    color='0.5',
                    linewidths=0,
                    zorder=69,
                    #latlon=True,
                    #transform=Geodetic()
                    )
    
    # Plot reports
    data = app_data[0]
    longitudes = [ d.reportlocation[1] for d in data ]
    latitudes = [ d.reportlocation[0] for d in data ]
    intensities = [ d.reportintensity for d in data ]
    
    if not len(intensities):
        return
    
    indexes = np.argsort(intensities)

    longitudes = [ longitudes[index] for index in indexes ]
    latitudes = [ latitudes[index] for index in indexes ]
    intensities = [ intensities[index] for index in indexes ]

    boundaries = np.arange(int(min(intensities)), max(intensities)+2)
    norm = BoundaryNorm(boundaries=boundaries, 
                        ncolors=256, 
                        extend='min')

    lolamap = ax_map(longitudes, latitudes) 
    ax_map.scatter(*lolamap,
                    c=intensities,
                    facecolors='w', 
                    edgecolors='w',
                    s=10,
                    marker='o',
                    linewidths=3,
                    zorder=68,
                    #latlon=True,
                    #transform=Geodetic()
                    )
    ax_map.scatter(*lolamap,
                    c=intensities,
                    facecolors='w', 
                    edgecolors='k',
                    s=10,
                    marker='o',
                    linewidths=0.5,
                    zorder=68,
                    #latlon=True,
                    #transform=Geodetic()
                    )    
    ln = ax_map.scatter(*lolamap,
                        c=intensities,
                        cmap=get_cmap(cmap),
                        norm=norm,
                        s=10,
                        marker='o',
                        linewidths=0,
                        zorder=69,
                        #latlon=True,
                        #transform=Geodetic()
                        )
    
    if ax_map.ax.figure.bmap.ax.get_legend():
        legend = ax_map.ax.get_legend()

        renderer = ax_map.ax.figure.canvas.get_renderer()
        bbox_display = legend.get_window_extent(renderer)  # Get legend bbox in figure (display) coordinates
        bbox_axes = bbox_display.transformed(ax_map.ax.transAxes.inverted())  # Convert to axis coordinates

        # Extract legend position and size in axis coordinates
        legend_x0, legend_y0 = bbox_axes.x0, bbox_axes.y0  # Bottom-left corner
        legend_x1, legend_y1 = bbox_axes.x1, bbox_axes.y1  # Top-right corner
        legend_width, legend_height = bbox_axes.width, bbox_axes.height

        cax = ax_map.ax.inset_axes([legend_x1+0.01, legend_y0+0.04, 0.02, legend_height*0.9-0.04])
    else:
        cax = ax_map.ax.inset_axes([1.03, 0.5, 0.02, 0.33])

    cb = ax_map.ax.figure.colorbar(ln,cax=cax)#,orientation='horizontal')
    cb.ax.set_yticks([t+0.5 for t in cb.ax.get_yticks()][:-1])
    cb.ax.yaxis.set_major_formatter(FuncFormatter(myroman))
    cb.set_label('EMS98')
    #cb.ax.xaxis.set_label_position('top')  # Move label to top
    #cb.ax.xaxis.tick_top()  # Move ticks to the top



def plot_notification_delays(app_data, ax_map,
                   cmap='turbo',#'rainbow',#'viridis'
                   ):
    if False:
        pos = ax.get_position()
        axins = ax.figure.add_axes([pos.x1-pos.width/12+0.01,
                                    pos.y1-pos.height/6,
                                    pos.width/3,
                                    pos.height/3],
                                facecolor='None',
                                zorder=999,
                                frame_on=False)
    else:
        legend = ax_map.ax.get_legend()

        renderer = ax_map.ax.figure.canvas.get_renderer()
        bbox_display = legend.get_window_extent(renderer)  # Get legend bbox in figure (display) coordinates
        bbox_axes = bbox_display.transformed(ax_map.ax.transAxes.inverted())  # Convert to axis coordinates

        # Extract legend position and size in axis coordinates
        legend_x0, legend_y0 = bbox_axes.x0, bbox_axes.y0  # Bottom-left corner
        legend_x1, legend_y1 = bbox_axes.x1, bbox_axes.y1  # Top-right corner
        legend_width, legend_height = bbox_axes.width, bbox_axes.height

        axins = ax_map.ax.inset_axes([legend_x0-0.02, legend_y1-0.01, legend_width+0.04, legend_width+0.04])

    size = 0.3
    
    delays = []
    intensities = []
    maxint = 0
    for notification in app_data[1]:
        reported = False
        for report in app_data[0]:
            if notification.userid == report.userid:

                report_intensity = int(report.reportintensity)
                intensities += [report_intensity]
                delays += [notification.delay]

                reported = True
                if report_intensity > maxint:
                    maxint=report_intensity

                break

        if not reported:
            intensities += [0]
            delays += [notification.delay]

    if not len(delays):
        return
    
    # Define bins for the histogram
    delay_bins = np.array([np.percentile(delays,p) for p in [1,16.1,50.1,84.1,99]])    # Bins for delay
    intensity_bins = np.arange(1, maxint+1)           # Bins for intensity

    # Compute 2D histogram
    vals, delay_edges, intensity_edges = np.histogram2d(delays, intensities, bins=[delay_bins, intensity_bins])
    
    # Convert the result into a list
    vals = np.asarray(vals.tolist())

    
    boundaries = np.arange(1, maxint+2)
    norm = BoundaryNorm(boundaries=boundaries, 
                        ncolors=256, 
                        extend='min')
    
    axins.pie(vals.sum(axis=1), 
              radius=1-size, 
              startangle=90,
              counterclock=False,
              colors=['%.1f'%(1-(dt+1)/(1.1*len(vals))) for dt,_ in enumerate(vals)],
              labels=['%d-%ds'%(delay_edges[dt],delay_edges[dt+1]) for dt,_ in enumerate(vals)],
              rotatelabels=True,
              textprops={'size': 'x-small'},
              wedgeprops=dict(width=size, 
                              edgecolor='w'),
              labeldistance=(1-size)/3)

    axins.pie(vals.flatten(), 
              radius=1, 
              startangle=90,
              counterclock=False,
              colors=[get_cmap(cmap)(norm(intensity_edges[i])) for dt,intensities in enumerate(vals) for i,_ in enumerate(intensities) ],
              wedgeprops=dict(width=size, 
                              edgecolor='w'))

    #axins.set(aspect="equal", title='Pie plot with `ax.pie`')

    axins.margins(0)
    axins.set_aspect("equal")
    axins.legend([],[],
                 frameon=False,
                 labelspacing=-.2,
                 loc="lower left", bbox_to_anchor=(0.03, .9),
                 title='Delays &\nIntensities',
                 title_fontproperties={'size':'xx-small',
                                        'weight':'bold'})
    #axins.set_axis_off()


def plot_report_hist(app_data,ax_map):

    data = app_data[0]

    distances = [ d.eventdistance for d in data ]
    intensities = [ d.reportintensity for d in data ]
    magnitudes = [ d.magnitude for d in data ]
    if not len(distances):
        return
    
    pos = ax_map.ax.get_position()
    axins = ax_map.ax.figure.add_axes([pos.x0,#+0.01,
                                        pos.y0,#+0.01,
                                        pos.width/3,
                                        pos.height/3],
                                    facecolor=(1, 1, 1, 0.5),#'None',
                                    zorder=999)

    xbins = np.linspace(np.min(distances)*.8,
                        np.percentile(distances,95)*1.2,
                        64)

    if False:
        xbins = np.logspace(np.log10(np.min(distances)*.8),
                            np.log10(np.percentile(distances,99)*1.2),
                            64)
    
    h = axins.hist2d(distances,intensities,
                    cmap=pqlx,
                    bins=(xbins,
                          np.arange(min(intensities)-0.5,max(intensities)+1.5,1)),
                    label='Reports',
                    #zorder=9998,
                    cmin = 1)
    
    h[3].set_clim(vmin=-1)

    axins.plot(xbins,
                [get_intensity(r,magnitudes[0]) for r in xbins],
                alpha=0.5,
                color='k',
                #    zorder=9999,
                linewidth=3,
                label='IPE')
    
    if False:
        legopts = {'title':'Reports &\nprediction',
                    'title_fontsize':'x-small',
                    'fontsize':'x-small'}
        
        l = axins.legend(**legopts)

    cax = axins.inset_axes([0.02, 0.02, 0.4, 0.03])
    cb = axins.figure.colorbar(h[3],
                               cax=cax,
                               orientation='horizontal',
                               ticklocation='top',
                               location='top')
    
    cb.ax.tick_params(axis='both', labelsize='small')

    formatter1 = EngFormatter(unit='',
                            places=0, 
                            sep="\N{THIN SPACE}",
                            usetex=True)  # U+2009
    cb.ax.xaxis.set_major_formatter(formatter1) 
    cb.ax.set_xlabel('Report\nnumber', 
                     fontsize='small',
                     loc='left')     

    axins.tick_params(axis='both', labelsize='small')

    axins.yaxis.set_major_formatter(FuncFormatter(myroman))

    axins.set_ylabel('Reported EMS98',fontsize='small')  
    axins.set_xlabel('Distance (km, hyp.)',fontsize='small')  

    formatter1 = EngFormatter(unit='',
                              places=0, 
                              sep="\N{THIN SPACE}",
                              usetex=True)  # U+2009
    axins.xaxis.set_major_formatter(formatter1)
    #axins.xaxis.set_minor_formatter(formatter1)


    # Remove top and right spines
    axins.spines['top'].set_visible(False)
    axins.spines['right'].set_visible(False)

    # Add a white patch below x-axis label
    axins.figure.canvas.draw()  # Ensure labels are properly positioned before adding patch

    bbox = axins.get_position()  # Get axis position
    x0, y0, x1, y1 = bbox.x0, bbox.y0, bbox.x1, bbox.y1

    # Add white rectangle spanning x-axis length (in axis coordinates)
    axins.figure.patches.append(Rectangle(
        (x0, y0 - 0.07),  # Slightly below x-axis label
        x1 - x0, 0.07,  # Span the entire x-axis width
        color='white',
        transform=axins.figure.transFigure,  # Use figure coordinates
        zorder=3
    ))

    # Add white rectangle spanning y-axis length (in axis coordinates)
    axins.figure.patches.append(Rectangle(
        (x0 - 0.05, y0),  # Slightly below x-axis label
        0.05, y1 - y0,  # Span the entire x-axis width
        color='white',
        transform=axins.figure.transFigure,  # Use figure coordinates
        zorder=3
    ))


    
    