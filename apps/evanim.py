#!/usr/bin/env python
from addons.fdsn import event_data
from addons.event import animate

def event_plots(catalog,
                eventstreams,
                eventinventories, 
                **args):

    for event, streams, inventory in zip(catalog,eventstreams,eventinventories):

        shorteventid = str(event.resource_id).split('/')[-1]
        
        ## Animate data and results
        anim = animate(event,
                    streams['acc'],
                    streams['disp'],
                    inventory,
                    **args)
        anim.save('data/%s_anim.mp4'%shorteventid)#,dpi=300)
        print('data/%s_anim.mp4'%shorteventid)

def main(**args):

        # Getting data
        catalog,eventstreams,eventinventories = event_data(**args)
        
        # Creating plots
        event_plots(catalog,eventstreams,eventinventories,**args)