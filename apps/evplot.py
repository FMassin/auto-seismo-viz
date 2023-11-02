#!/usr/bin/env python
from addons.fdsn import event_data
from addons.stream import ploteqdata
from addons.core import eewmap
from addons.catalog import performance_timelines
from matplotlib.patheffects import withStroke
from matplotlib.text import Text
path_effects=[withStroke(linewidth=2,foreground="w")]

def event_plots(catalog,
                eventstreams,
                eventinventories,
                **args):

    saveopt = {'dpi':512,
               'facecolor':'none',
               'transparent':False}
    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):

        shorteventid = str(event.resource_id).split('/')[-1]

        ## Map results
        fig = eewmap({'event':event,
                    'inventory':inventory},
                    reference=False,
                    stationgroups={})

        for ax in fig.axes:
            for t in ax.findobj(Text):
                if not t.get_path_effects():
                    t.set(path_effects=path_effects)
        fig.savefig('data/%s_map.png'%shorteventid,bbox_inches='tight',**saveopt)
        print('data/%s_map.png'%shorteventid)

        ## Plot data
        fig = ploteqdata(streams['acc'].select(channel='*b')+streams['acc'].select(channel='*X'),event,inventory,lim=160)
        fig.tight_layout()
        for ax in fig.axes:
            for t in ax.findobj(Text):
                if not t.get_path_effects():
                    t.set(path_effects=path_effects)
        fig.savefig('data/%s_data.png'%shorteventid,bbox_inches=None,**saveopt)
        print('data/%s_data.png'%shorteventid)

        ## Plot results timeline
        fig = performance_timelines(event)
        for ax in fig.axes:
            for t in ax.findobj(Text):
                if not t.get_path_effects():
                    t.set(path_effects=path_effects)
        fig.savefig('data/%s_timeline.png'%shorteventid,bbox_inches='tight',**saveopt)
        #fig.savefig('data/%s_timeline.pdf'%shorteventid,bbox_inches='tight',**saveopt)
        print('data/%s_timeline.png'%shorteventid)

def main(**args):

        # Getting data
        catalog,eventstreams,eventinventories = event_data(**args)

        # Creating plots
        event_plots(catalog,eventstreams,eventinventories,**args)
