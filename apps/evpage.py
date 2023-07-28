#!/usr/bin/env python
from sys import argv
from string import Template
from glob import glob

ids = list(set([f.split('/')[-1].split('_')[0] for f in glob("%s/*_map.png"%argv[-1])]))


pagetemplate = '''
Event $evtid
============

${anim}

Summary map
-----------

.. figure:: ../${evtstore}/${evtid}_map.png
    :alt: About...

    Data and event parameters summary map for event ${evtid}.

Solution time evolution
-----------------------

.. figure:: ../${evtstore}/${evtid}_timeline.png
    :alt: About...

    Event parameters time evolution for event ${evtid}.

Ground motion envelopes
-----------------------

.. figure:: ../${evtstore}/${evtid}_data.png
    :alt: About...

    Ground motion envelopes for event ${evtid}.

${stations}

Downloads
---------

Metadata:
 - :download:`station metadata <../${evtstore}/${evtid}.stationxml>`
 - :download:`event metadata <../${evtstore}/${evtid}.quakeml>`

Data:
 - :download:`raw data <../${evtstore}/${evtid}.raw.mseed>`
 - :download:`acceleration data <../${evtstore}/${evtid}.acc.mseed>`
 - :download:`velocity data <../${evtstore}/${evtid}.vel.mseed>`
 - :download:`displacement data  <../${evtstore}/${evtid}.disp.mseed>`
'''

t = Template(pagetemplate)


animtemplate = Template('''
Summary Animation
-----------------

.. raw:: html

    <video controls width="100%s" src="../../../%s/${evtid}_anim.mp4"></video>

'''%('%',argv[-1]))

stationhead = '''
Stations
--------

Raw data:
'''

stationtemplate = Template('''
.. figure:: ../${evtstore}/${evtid}/${station}.${out}.png
    :alt: About...

    Seismograms and spectrograms for ${out} data at station ${station} for event ${evtid}.
''')

for id in ids:
    if argv[-2] in ids and argv[-2] != id:
        continue
    stations = list(set([f.split('/')[-1].split('.raw')[0] for f in glob("%s/%s/*.raw.png"%(argv[-1],id))]))

    with open('events/%s.rst'%id, 'w', encoding='utf-8') as file:
        anim = ''
        if glob('%s/%s_anim.mp4'%(argv[-1],id)):
            anim += animtemplate.substitute({ 'evtid': id,
                                   'evtstore': argv[-1]})
        stationsraw = ''
        subst = {'evtid': id,
                 'evtstore':argv[-1],
                 'out':'raw',
                 }
        for station in stations:
            subst['station'] = station
            stationsraw += stationtemplate.substitute(subst)
        if len(stationsraw):
            stationsraw =  stationhead+stationsraw

        eventpage = t.substitute({ 'evtid': id,
                                   'evtstore': argv[-1],
                                   'stations': stationsraw,
                                   'anim': anim
                                   })

        file.writelines(eventpage)
