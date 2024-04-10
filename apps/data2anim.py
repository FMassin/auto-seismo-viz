#!/usr/bin/env python
from glob import glob
from addons.bmap import bmap
from cartopy.crs import Geodetic
from matplotlib.animation import FuncAnimation, writers
from matplotlib.colors import LogNorm
import numpy as np
from obspy import UTCDateTime
from os.path import splitext
from obspy import read_events

def main(folder = '/Users/fred/Documents/Projects/Collaborations/Song/04032358_Hualien_Eathquake/1712102298',
         begin = '2024-04-02T23:58:18.34',
         output = "data/Hualien.doesnotmatter",
         catalog = "/Users/fred/Documents/Projects/NaiNo-Kami/auto-seismo-viz/exports/CWA2024gnxs.xml",
         quickndirty = False,
         **args):
    # read peak ground motion data from files and write to animation.mp4

    quickndirty = bool(quickndirty)
    catalog = read_events(catalog)
    print(catalog)

    preferred_org = catalog[0].preferred_origin()

    preferred_mag = catalog[0].preferred_magnitude()
    if preferred_mag is None:
        mcts = [m.creation_info.creation_time for m in catalog[0].magnitudes]
        preferred_mag = catalog[0].magnitudes[np.argmax(mcts)]

    mags = [ m for m in catalog[0].magnitudes if m.magnitude_type == 'Mfd'  ]


    output = splitext(output)[0]
    dataset = {}
    datatimes = {}
    datafd = {}
    for f in glob('%s/data_*'%folder):
        with open(f,'r') as file:
            data = [line.split() for line in file.read().split('\n') if len(line)>0]
        index = int(f.split('_')[-1])
        dataset[index] = data[1:]
        datatimes[index] = UTCDateTime(begin) + index
        
        magindex = np.argmin([ abs(m.creation_info.creation_time - datatimes[index]-0.1) for m in mags])
        if abs( mags[magindex].creation_info.creation_time - datatimes[index]-0.1) > 1:
            datafd[index] = None
        else:
            datafd[index] = mags[magindex]
            mct = mags[magindex].creation_info.creation_time
            t = datatimes[index]
            print(f'{magindex} {mct} versus {t}')
            mags.pop(magindex)

    # cartopy map
    lon = [float(l[1]) for k in dataset for l in dataset[k]]
    lat = [float(l[0]) for k in dataset for l in dataset[k]]
    c = [float(l[4])/100 for k in dataset for l in dataset[k]]

    ax = bmap(bounds=[np.percentile(lon,2), 
                      np.percentile(lon,98), 
                      np.percentile(lat,2), 
                      np.percentile(lat,98)],
              right=False,
              left=True,
              top=False,
              bottom=True,
              quickndirty=quickndirty)
    
    ref, = ax.plot(preferred_org.longitude,
            preferred_org.latitude,
            'k*',
            alpha=0.9,
            markersize=25,
            transform=Geodetic()) 
    
    fd, = ax.plot(np.nan,
                    np.nan,
                    '*',
                    color='C1',
                    alpha=0.4,
                    markersize=25,
                    transform=Geodetic())   

    lm = ax.scatter([None,None],[None,None],s=[5,20], c='k',
                    linewidths=2,
                    transform=Geodetic(),
                    norm=LogNorm(vmin=np.percentile(c,5), vmax=max(c)))   
    
    ref.set_zorder(lm.get_zorder())
    fd.set_zorder(lm.get_zorder())
    
    ln = ax.scatter([None,None],[None,None],s=[5,20], c=[None,None],
                    transform=Geodetic(),
                    norm=LogNorm(vmin=np.percentile(c,5), vmax=max(c)))    
    
    cax = ax.inset_axes([1.03, 0.3, 0.05, 0.68])
    cb = ax.figure.colorbar(ln,cax=cax)
    cb.set_label(r'Peak ground acceleration ($m/s^2$)')

    def dt2s(dt,st):
        if dt < -1 :
            print(f'Warning: {st}: {dt} s delay')
        dt = min([dt, 15])
        dt = max([0, dt])
        s = 20 - dt
        return s
    
    def s2dt(s):
        dt = 20 - s 
        return dt
    
    kw = dict(prop="sizes", num=4, color=ln.cmap(1), fmt="{x:.0f} s",
              func=lambda s: s2dt(s))
    ax.legend(*ln.legend_elements(**kw),
              loc='lower left', bbox_to_anchor=(1, 0),#loc="lower right", 
              title="Delay", fontsize="small")


    def update(frame,ln,lm,fd):
        reftime = UTCDateTime('1970-01-01T00:00:00')
        lon = [float(l[1]) for l in dataset[frame]]
        lat = [float(l[0]) for l in dataset[frame]]
        c = [float(l[4])/100 for l in dataset[frame]]
        s = [dt2s(datatimes[frame]-(reftime+int(l[3])),l[2]) for l in dataset[frame]]
        
        order = np.argsort(c)
        lon = [lon[i] for i in order]
        lat = [lat[i]for i in order]
        s = [s[i]for i in order]
        c = [c[i] for i in order]
        
        ln.set_array(c)  
        for l in [ln, lm]:
            l.set_offsets(np.c_[lon, lat])  
            l.set_sizes(s) 

        if datafd[frame] is not None:
            o = datafd[frame].origin_id.get_referred_object()
            fd.set_data(o.longitude, o.latitude)            
            fd.set_markersize(int(25*datafd[frame].mag/preferred_mag.mag))

        time = str(datatimes[frame])[:19]
        ax.set_title(f'FinDer #{frame}\n{time}')

        #ax.relim()
        #ax.autoscale_view()
        return ln,

    ani = FuncAnimation(ax.figure, update, 
                        frames=range(len(list(dataset.keys()))), 
                        fargs=[ln,lm,fd],
                        blit=True)
    
    # Save the animation as a movie
    writer = writers['ffmpeg'](fps=1)  # Adjust fps to control the speed of the video
    ani.save(f'{output}.mp4', writer=writer)