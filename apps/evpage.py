#!/usr/bin/env python
from sys import argv
from string import Template
from glob import glob
from numpy import argsort
import xml.etree.ElementTree as ET
from  shutil import copyfile
from os import makedirs, system
from os.path import isdir

def qml2desc(file):
    
    # Passing the path of the
    # xml document to enable the
    # parsing process
    tree = ET.parse(file)
    
    # getting the parent tag of
    # the xml document
    root = tree.getroot()
    
    for child in root[0][0]:
        if 'preferredOriginID' in child.tag :
            prefOrId = child.text
        if 'preferredMagnitudeID' in child.tag :
            prefMagId = child.text


    for child in root[0][0]:
        if 'creationInfo' in  child.tag :
            for granchild in child:
                if 'agencyID' in  granchild.tag :
                    agency = granchild.text

        if 'origin' in  child.tag and child.attrib['publicID'] == prefOrId:
            for granchild in child:
                if 'time' in  granchild.tag :
                    time = granchild[0].text
        
        if 'magnitude' in  child.tag and child.attrib['publicID'] == prefMagId:
            for granchild in child:
                if 'mag' in  granchild.tag :
                    mag = granchild[0].text

    return mag, time, agency

ids = list(set([f.split('/')[-1].split('_')[0] for f in glob("%s/*_map.png"%argv[-1])]))


pagetemplate = '''
${agency} M${mag} event on ${date} 
======================================================================

At ${time} UTC (``$evtid``).

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
Summary animation
-----------------

.. raw:: html

    <video controls width="100%s" src="../_images/${evtid}_anim.mp4"></video>
    <center><p><small>Data and event parameters animation for event ${evtid}. Intensity and acceleration do not match accurately.</small></p></center>
                        
'''%('%'))

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

pages = []
magnitudes = []
copyfiles = []
for id in ids:
    if argv[-2] in ids and argv[-2] != id:
        continue
    stations = list(set([f.split('/')[-1].split('.raw')[0] for f in glob("%s/%s/*.raw.png"%(argv[-1],id))]))
    
    mag, time, agency = qml2desc( '%s/%s.quakeml' % ( argv[-1], id ))

    magnitudes += [float(mag)]
    pages += ['events/%s.rst'%id]

    with open('events/%s.rst'%id, 'w', encoding='utf-8') as file:
        anim = ''
        if glob('%s/%s_anim.mp4'%(argv[-1],id)):
            anim += animtemplate.substitute({ 'evtid': id,
                                   'evtstore': argv[-1]})
            
            if not isdir("_build/html/_images/"):
                makedirs('_build/html/_images/')
            copyfiles += [[ "%s/%s_anim.mp4" % ( argv[-1], id ), 
                           "_build/html/_images/%s_anim.mp4" % ( id ) 
                           ]]
            
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

        eventpage = t.substitute({ 'agency':agency,
                                   'evtid': id,
                                   'time': time[11:19],  
                                   'date': time[:10], 
                                   'mag': int(float(mag)),                 
                                   'evtstore': argv[-1],
                                   'stations': stationsraw,
                                   'anim': anim
                                   })

        file.writelines(eventpage)
        
        

eventspage = '''

.. toctree::
   :maxdepth: 3
   :caption: Events:
   
'''
for i in argsort(magnitudes)[::-1]:
    eventspage += '''
   '''+pages[i]
    
eventspage += '''
   '''
print(eventspage)
with open('events.rst', 'w', encoding='utf-8') as file:
    file.writelines(eventspage)

system('make clean html')

for args in copyfiles:
    print('copy',*args,'..')
    copyfile(*args)