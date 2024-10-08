#!/usr/bin/env python
from addons.fdsn import event_data
from addons.station import plotstationdata, plotallstationsdata
from os import makedirs
from matplotlib.patheffects import withStroke
from matplotlib.text import Text
path_effects=[withStroke(linewidth=2,foreground="w")]

saveopt = {}
    
def event_plots(catalog,
                eventstreams,
                eventinventories, 
                dpi=256,
                facecolor='none',
                transparent=False,
                max_num_station=12,
                max_num_single_station=4,
                **args):

    if isinstance(max_num_station,str):
        max_num_station = int(max_num_station)
        
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):
        shorteventid = str(event.resource_id).split('/')[-1]

        figs = plotallstationsdata(streams,
                                   event,
                                   inventory,
                                   max_num_station=max_num_station,
                                   **args)
        print('Saving figs')
        for fig in figs:
            figname = 'data/%s_seisgram_%s.png'%(shorteventid,fig.basename)
            fig.savefig(figname,
                        bbox_inches=None,
                        dpi=dpi,
                        facecolor=facecolor,
                        transparent=transparent)
            print(figname)

        if float(max_num_single_station) == 0:
             continue
        
        figs = plotstationdata(streams,event,inventory,max_num_station=max_num_single_station)
        for fig in figs:
            for ax in fig.axes:
                for t in ax.findobj(Text):
                    if not t.get_path_effects():
                        t.set(path_effects=path_effects) 
            figname = 'data/%s/%s.png'%(shorteventid,fig.basename)
            makedirs('data/%s'%shorteventid,exist_ok=True)
            fig.savefig(figname,
                        bbox_inches=None,
                        dpi=dpi,
                        facecolor=facecolor,
                        transparent=transparent)
            print(figname)


def main(**args):

        # Getting data
        catalog,eventstreams,eventinventories = event_data(**args)
        
        # Creating plots
        event_plots(catalog,eventstreams,eventinventories,**args)