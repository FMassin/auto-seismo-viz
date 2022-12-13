#!/usr/bin/env python
from addons.fdsn import event_data
from addons.station import plotstationdata
from os import makedirs
from matplotlib.patheffects import withStroke
from matplotlib.text import Text
path_effects=[withStroke(linewidth=2,foreground="w")]

saveopt = {}
    
def event_plots(catalog,
                eventstreams,
                eventinventories, 
                dpi=512,
                facecolor='none',
                transparent=False,
                **args):
    
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):
        shorteventid = str(event.resource_id).split('/')[-1]
        figs = plotstationdata(streams,event,inventory)
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