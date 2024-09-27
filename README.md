# Auto Seismo-Viz
Automatic visualisation for seismic events  
![Example map](data/2022swpksz_map.pngt)

## Install 
```bash
python -m pip install --requirement requirements.txt      
```

## Usage from webservices
```bash
./seismoviz <app> <options> 
```

The available app are in `apps/`. E.g.:
```bash
./seismoviz evplot \
    catalog_uri=USGS,http://<your fdsnws event server>:8080 \
    inventory_url=http://<your fdsnws station server>/ \
    stream_url=http://<your fdsnws dataselect server>/ \
    eventid=us6000j9es \
    nseconds=9999999999 \
    debug=1 \
    limit=1
```
## Usage from file
Import event data into a local data folder:
```bash
./seismoviz evimport \
  catalog=single_event.xml  \
  stream=data.mseed  \
  inventory=inventory.xml
```
Run evplot with eventid from `single_event.xml` e.g., `g`fz2019qmev`, creates plots in local data folder:
```bash
./seismoviz evplot files="gfz2019qmev" 
```
> Export single event parameters from seiscomp with `scxmldump -f -E gfz2019qmev  -d sqlite3://event_db.sqlite -PAMFm > single_event.xml `


## To do:
- App usage helper
- More app :
  - station and event tables
- Notification: webhook
- Install via `pip`