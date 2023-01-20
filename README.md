# Auto Seismo-Viz
Automatic visualisation for seismic events  
![Example map](data/2022swpksz_map.png)

## Install 
```bash
python -m pip install --requirement requirements.txt      
```

## Usage
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

## To do:
- App usage helper
- More app 
  - station and event tables
  - html or panel report
- Notification: webhook
- Install via `pip`