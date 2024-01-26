
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

import glob, matplotlib.pyplot, matplotlib.patheffects, matplotlib.ticker

def argsort(seq):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__)

def mapcities(bmap,
              cities={'San Jose':[9.94456,-84.11525],
                      'Managua':[12.12801,-86.29536],
                      'San Salvador':[13.70178, -89.21693],
                      'Guatemala City':[14.62778, -90.51520]},
              optcities={#'weight':"bold",
                         'color':"k",
                         'fontsize':'x-small',
                         'zorder':9},
              optcitydot={'markeredgecolor':'w',
                          'markeredgewidth':.8,
                          'markerfacecolor':'none',
                          'zorder':9},
                          **opt):
    
    extent=bmap.get_extent(opt['transform'])
    
    dot_path_effects=[matplotlib.patheffects.withStroke(linewidth=2,foreground="k")]
    label_path_effects=[matplotlib.patheffects.withStroke(linewidth=2,foreground="w")]
    for city in cities:
        text = bmap.plot(*cities[city][::-1],'o',path_effects=dot_path_effects,**optcitydot,**opt)
        text = bmap.text(*[d-(extent[3]-extent[2])*.01 for d in cities[city][::-1]],city,va='top',ha='right',path_effects=label_path_effects,clip_on=True,**optcities,**opt)

import cartopy.io.shapereader as shpreader
def mapcountrynames(ax,**opt):
    extent=ax.get_extent(opt['transform'])
    shpfilename = shpreader.natural_earth(resolution='110m',
                                      category='cultural',
                                      name='admin_0_countries')

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()
    for country in countries:

        x = country.geometry.centroid.x        
        y = country.geometry.centroid.y

        if  x<extent[0] or x>extent[1]:
            continue
        if y<extent[2] or y>extent[3]:
            continue

        ax.text(x, y+(extent[3]-extent[2])*.02, country.attributes['NAME'], 
                color='k', ha='center', va='bottom', alpha=.7,clip_on=True,
                path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w", alpha=.7)], 
                **opt)

def sanitize_lonlist(lons):
    new_list = []
    oldval = 0
    # used to compare with the adjacent longitudes
    # and values exceed, disconnect linestring
    treshold = 10   
    for ix,ea in enumerate(lons):
        diff = oldval - ea
        if (ix>0):
            if (diff>treshold):
                ea = ea+360
        oldval = ea
        new_list.append(ea)
    return new_list

import geopandas,requests,glob,os
def mapfaults(ax,
              url='https://raw.githubusercontent.com/GEMScienceTools/gem-global-active-faults/master/geojson/gem_active_faults_harmonized.geojson',
              fallback_label='Active faults',
              linestyle_str = [(0, (3, 1, 1, 1, 1, 1)),':','-.','--'],
              legendfaults=True,
              **opt):
    
    

    if 'color' not in opt:
        opt['color'] = 'magenta'
    if 'alpha' not in opt:
        opt['alpha'] = .4
    if 'linewidth' not in opt:
        opt['linewidth'] = 0.8

    f=os.environ['HOME']+'/.local/share/cartopy/faults.geojson'
    if not glob.glob(f):
        print('downloading',url,'in',f)
        r = requests.get(url, allow_redirects=True)
        open(f, 'wb').write(r.content)

    faults=geopandas.read_file(f)
    extent=ax.get_extent(opt['transform'])
    labels=[]
    slip_types = list(set([slip_type for slip_type in faults.slip_type.values if None is not slip_type and  '_' in slip_type] ))

    for i,slip_type in enumerate(faults.slip_type.values):

        # grab x and y of the first geometry object
        geometry  = faults.iloc[i].geometry
        xs = geometry.xy[0].tolist()
        if not len( [x for x in xs if x>=extent[0] and x<=extent[1] ]):
            continue

        ys = geometry.xy[1].tolist()
        if not len( [x for x in ys if x>=extent[2] and x<=extent[3] ]):
            continue

        label=None

        # special cases
        tmpopt=opt.copy()
        if slip_type is not None and  '_' in slip_type:
            tmpopt['linewidth'] = opt['linewidth']*2
            tmpopt['linestyle']=linestyle_str[slip_types.index(slip_type)]
            if slip_type not in labels:
                label=slip_type.replace('_',' ')
                labels += [slip_type]

        elif fallback_label not in labels:
            label = fallback_label
            labels += [fallback_label]

        if not legendfaults:
            label=None

        # plot the geometry using sanitized values
        path_effects = [matplotlib.patheffects.withStroke(linewidth=tmpopt['linewidth']*2,#.5,
                                                      foreground="w",
                                                      alpha=.5)]
        ax.plot(xs, ys, 
                label=label,
                path_effects=path_effects, #
                **tmpopt)

import matplotlib.transforms as mtransforms
class _TransformedBoundsLocator:
    """
    Axes locator for `.Axes.inset_axes` and similarly positioned Axes.
    The locator is a callable object used in `.Axes.set_aspect` to compute the
    axes location depending on the renderer.
    """

    def __init__(self, bounds, transform):
        """
        *bounds* (a ``[l, b, w, h]`` rectangle) and *transform* together
        specify the position of the inset Axes.
        """
        self._bounds = bounds
        self._transform = transform

    def __call__(self, ax, renderer):
        # Subtracting transSubfigure will typically rely on inverted(),
        # freezing the transform; thus, this needs to be delayed until draw
        # time as transSubfigure may otherwise change after this is evaluated.
        return mtransforms.TransformedBbox(
            mtransforms.Bbox.from_bounds(*self._bounds),
            self._transform - ax.figure.transSubfigure)

def makefigax(fig=None,ax=None,axprop={},figprop={}):
    if ax is not None:
        ax2 = ax.figure.add_axes(ax.get_position(True), 
                                **axprop, 
                                axes_locator=_TransformedBoundsLocator([0, 0, 1, 1], ax.transAxes))
        fig = ax2.figure
        ax.remove()
        ax=ax2
    else:
        if fig is not None:
            pass
        else :
            fig = matplotlib.pyplot.figure(**figprop)
        ax = matplotlib.pyplot.axes(**axprop)
    return ax

def bmap(bounds=[-89, -83, 8, 14],
         pad=1/10,
         figsize=(5,5),
         right=False,
         left=False,
         top=False,
         bottom=False,
         fig=None,
         ax=None,
         aspect=None,
         #url='https://server.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}.jpg',
         #url='https://tile.osm.ch/osm-swiss-style/{z}/{x}/{y}.png',
         url='https://tile.osm.ch/switzerland/{z}/{x}/{y}.png',
         zoom=8,#10,
         legendfaults=True,
         label_style={'size':'small',
                      'path_effects':[matplotlib.patheffects.withStroke(linewidth=2,foreground="w")]}):
    padlo=(bounds[1]-bounds[0])*pad
    padla=(bounds[3]-bounds[2])*pad
    midlo=(bounds[1]+bounds[0])/2
    mila=(bounds[3]+bounds[2])/2

    projection=ccrs.Orthographic(central_longitude=midlo,central_latitude=mila)
    

    if ax is not None and aspect is None:
        aspect = ax.get_box_aspect()


    ax=makefigax(fig=fig,
                 ax=ax,
                 axprop={'projection':projection},
                 figprop={'figsize':figsize}
                 )


    ax.set_extent([bounds[0]-padlo, bounds[1]+padlo,
                   bounds[2]-padla, bounds[3]+padla])
    
    print('before map:',aspect)
    print('after map:',ax.get_box_aspect())

    gl=ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    if right or not left:
        gl.left_labels= False
        gl.right_labels= True
    else :
        gl.right_labels= False
        gl.left_labels= True
    if  top or not bottom:
        gl.bottom_labels= False
        gl.top_labels= True
    else :
        gl.top_labels= False
        gl.bottom_labels= True
    gl.ylabel_style=label_style
    gl.xlabel_style=label_style
    gl.xpadding=7
    gl.ypadding=7

    ax.yaxis.tick_right()
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.xaxis.tick_top()
    ax.tick_params(axis="y", direction="out", length=4)
    ax.tick_params(axis="x", direction="out", length=4)


    #ax.add_feature(cfeature.LAND, alpha=0.95,color='w')
    #ax.add_feature(cfeature.OCEAN, alpha=0.95,color='w')
    #ax.add_feature(cfeature.LAKES, alpha=0.95,color='w')

    ax.add_feature(cfeature.BORDERS, linewidth=2,color='w')
    ax.add_feature(cfeature.BORDERS, linewidth=1,color='.5')
    ax.add_feature(cfeature.COASTLINE, linewidth=2,color='w')    
    ax.add_feature(cfeature.COASTLINE, linewidth=1,color='.5') 
    ax.add_feature(cfeature.RIVERS, linewidth=.5)
    
    if False:
        ax.add_image(cimgt.GoogleTiles(url=url), zoom)

    if False:
        ax.add_wms(wms='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv',
                layers=['GEBCO_LATEST'],
                alpha=1/3,
                zorder=1
                )
        
    if False:
        # Create a Stamen Terrain instance.
        stamen_terrain = cimgt.Stamen('terrain-background')
        # Add the Stamen data at zoom level 8.
        ax.add_image(stamen_terrain, 7,alpha=1/2)   
    
    mapcities(ax,transform=ccrs.Geodetic())
    mapfaults(ax,transform=ccrs.Geodetic(),legendfaults=legendfaults)
    mapcountrynames(ax,transform=ccrs.Geodetic())

    return ax

   