#!/usr/bin/env python

from glob import glob
from os.path import splitext
import numpy as np

from matplotlib.animation import FuncAnimation, writers
from matplotlib.colors import LogNorm
from matplotlib.pyplot import get_cmap
from matplotlib.ticker import EngFormatter

from cartopy.crs import Geodetic

from obspy import UTCDateTime
from obspy import read_events

from addons.bmap import bmap




def main(folder = '/Users/fred/Documents/Projects/Collaborations/Song/04032358_Hualien_Eathquake/1712102298',
         begin = '2024-04-02T23:58:18.34',
         output = "data/CWA2024gnxs.doesnotmatter",
         catalog = "/Users/fred/Documents/Projects/NaiNo-Kami/auto-seismo-viz/data/CWA2024gnxs.quakeml",
         quickndirty = False,
         cmap='nipy_spectral',
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

    mags = [ m for m in catalog[0].magnitudes if m.magnitude_type == 'Mfd' and m.mag > -1.0 ] # ignore SM origins


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

    ax = bmap(bounds=[np.percentile(lon,5), 
                      np.percentile(lon,90), 
                      np.percentile(lat,5), 
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
                    transform=Geodetic())   
    
    ref.set_zorder(lm.get_zorder())
    fd.set_zorder(lm.get_zorder())
    
    ln = ax.scatter([None,None],[None,None],s=[5,20], c=[None,None],
                    transform=Geodetic(),
                    cmap=get_cmap(cmap),
                    norm=LogNorm(vmin=np.percentile(c,5), vmax=max(c)))    
    
    cax = ax.inset_axes([1.03, 0.3, 0.05, 0.68])
    
    cb = ax.figure.colorbar(ln,cax=cax)

    formatter1 = EngFormatter(unit='m/s$^2$',
                              places=0, 
                              sep="\N{THIN SPACE}",
                              usetex=True)  # U+2009
    cb.ax.yaxis.set_major_formatter(formatter1)

    cb.ax.set_title('Peak\nground\nacceleration',loc='left',fontsize="small")

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
              title="Delay PGA", fontsize="small")

    ax.set_title(f'FinDer #none\nNone',loc='left')

    def update(frame,ln,lm,fd):
        reftime = UTCDateTime('1970-01-01T00:00:00')
        lon = [float(l[1]) for l in dataset[frame]]
        lat = [float(l[0]) for l in dataset[frame]]
        c = [float(l[4])/100 for l in dataset[frame]]
        s = [dt2s(datatimes[frame]-(reftime+int(float(l[3]))),l[2]) for l in dataset[frame]]
        
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
            delay = datatimes[frame] - preferred_mag.origin_id.get_referred_object().time
            magval = datafd[frame].mag
            
            ax.set_title(f'FinDer #{frame}: Mfd {magval:.1f} @ Ot+{delay:.1f} s\n{time}',loc='left')
        else:
            print(f'Frame {frame} has no magnitude?!')

        #ax.relim()
        #ax.autoscale_view()
        return ln,

    saveopt = {'dpi':512,
               'facecolor':'none',
               'transparent':False}
    
    ax.figure.tight_layout()
    ax.figure.savefig(f'{output}_datanim.png',bbox_inches='tight',**saveopt)

    if True:
        ani = FuncAnimation(ax.figure, update, 
                            frames=range(len(list(dataset.keys()))), 
                            fargs=[ln,lm,fd],
                            blit=True)
        # Save the animation as a movie
        writer = writers['ffmpeg'](fps=1)  # Adjust fps to control the speed of the video
        ani.save(f'{output}_datanim.mp4', writer=writer)

    ax.figure.savefig(f'{output}_datanim.png',bbox_inches='tight',**saveopt)
    
    