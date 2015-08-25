"""
Module for plotting data on maps with matplotlib.

Contains the :class:`Basemap` class (which does most of the
heavy lifting), and the following functions:

:func:`interp`:  bilinear interpolation between rectilinear grids.

:func:`maskoceans`:  mask 'wet' points of an input array.

:func:`shiftgrid`:  shifts global lat/lon grids east or west.

:func:`addcyclic`: Add cyclic (wraparound) point in longitude.
"""
from matplotlib import __version__ as _matplotlib_version
from matplotlib.cbook import is_scalar, dedent
# check to make sure matplotlib is not too old.
_mpl_required_version = '0.98'
if _matplotlib_version < _mpl_required_version:
    msg = dedent("""
    your matplotlib is too old - basemap requires version %s or
    higher, you have version %s""" %
    (_mpl_required_version,_matplotlib_version))
    raise ImportError(msg)
from matplotlib import rcParams, is_interactive
from matplotlib.collections import LineCollection
from matplotlib.patches import Ellipse, Circle, Polygon, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
from mpl_toolkits.basemap import pyproj
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.image import imread
import sys, os, math
#from .proj import Proj
import numpy as np
import numpy.ma as ma
import _geoslib
import functools

# basemap data files now installed in lib/matplotlib/toolkits/basemap/data
# check to see if environment variable BASEMAPDATA set to a directory,
# and if so look for the data there.
if 'BASEMAPDATA' in os.environ:
    basemap_datadir = os.environ['BASEMAPDATA']
    if not os.path.isdir(basemap_datadir):
        raise RuntimeError('Path in environment BASEMAPDATA not a directory')
else:
    basemap_datadir = os.sep.join([os.path.dirname(__file__), 'data'])

__version__ = '1.0.7'

# module variable that sets the default value for the 'latlon' kwarg.
# can be set to True by user so plotting functions can take lons,lats
# in degrees by default, instead of x,y (map projection coords in meters).
latlon_default = False

# supported map projections.
_projnames = {'cyl'      : 'Cylindrical Equidistant',
             'merc'     : 'Mercator',
             'tmerc'    : 'Transverse Mercator',
             'omerc'    : 'Oblique Mercator',
             'mill'     : 'Miller Cylindrical',
             'gall'     : 'Gall Stereographic Cylindrical',
             'cea'      : 'Cylindrical Equal Area',
             'lcc'      : 'Lambert Conformal',
             'laea'     : 'Lambert Azimuthal Equal Area',
             'nplaea'   : 'North-Polar Lambert Azimuthal',
             'splaea'   : 'South-Polar Lambert Azimuthal',
             'eqdc'     : 'Equidistant Conic',
             'aeqd'     : 'Azimuthal Equidistant',
             'npaeqd'   : 'North-Polar Azimuthal Equidistant',
             'spaeqd'   : 'South-Polar Azimuthal Equidistant',
             'aea'      : 'Albers Equal Area',
             'stere'    : 'Stereographic',
             'npstere'  : 'North-Polar Stereographic',
             'spstere'  : 'South-Polar Stereographic',
             'cass'     : 'Cassini-Soldner',
             'poly'     : 'Polyconic',
             'ortho'    : 'Orthographic',
             'geos'     : 'Geostationary',
             'nsper'    : 'Near-Sided Perspective',
             'sinu'     : 'Sinusoidal',
             'moll'     : 'Mollweide',
             'hammer'   : 'Hammer',
             'robin'    : 'Robinson',
             'kav7'     : 'Kavrayskiy VII',
             'eck4'     : 'Eckert IV',
             'vandg'    : 'van der Grinten',
             'mbtfpq'   : 'McBryde-Thomas Flat-Polar Quartic',
             'gnom'     : 'Gnomonic',
             'rotpole'  : 'Rotated Pole',
             }
supported_projections = []
for _items in _projnames.items():
    supported_projections.append(" %-17s%-40s\n" % (_items))
supported_projections = ''.join(supported_projections)

_cylproj = ['cyl','merc','mill','gall','cea']
_pseudocyl = ['moll','robin','eck4','kav7','sinu','mbtfpq','vandg','hammer']
_dg2rad = math.radians(1.)
_rad2dg = math.degrees(1.)

# projection specific parameters.
projection_params = {'cyl'      : 'corners only (no width/height)',
             'merc'     : 'corners plus lat_ts (no width/height)',
             'tmerc'    : 'lon_0,lat_0,k_0',
             'omerc'    : 'lon_0,lat_0,lat_1,lat_2,lon_1,lon_2,no_rot,k_0',
             'mill'     : 'corners only (no width/height)',
             'gall'     : 'corners only (no width/height)',
             'cea'      : 'corners only plus lat_ts (no width/height)',
             'lcc'      : 'lon_0,lat_0,lat_1,lat_2,k_0',
             'laea'     : 'lon_0,lat_0',
             'nplaea'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'splaea'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'eqdc'     : 'lon_0,lat_0,lat_1,lat_2',
             'aeqd'     : 'lon_0,lat_0',
             'npaeqd'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'spaeqd'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'aea'      : 'lon_0,lat_0,lat_1',
             'stere'    : 'lon_0,lat_0,lat_ts,k_0',
             'npstere'  : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'spstere'  : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'cass'     : 'lon_0,lat_0',
             'poly'     : 'lon_0,lat_0',
             'ortho'    : 'lon_0,lat_0,llcrnrx,llcrnry,urcrnrx,urcrnry,no width/height',
             'geos'     : 'lon_0,satellite_height,llcrnrx,llcrnry,urcrnrx,urcrnry,no width/height',
             'nsper'    : 'lon_0,satellite_height,llcrnrx,llcrnry,urcrnrx,urcrnry,no width/height',
             'sinu'     : 'lon_0,lat_0,no corners or width/height',
             'moll'     : 'lon_0,lat_0,no corners or width/height',
             'hammer'   : 'lon_0,lat_0,no corners or width/height',
             'robin'    : 'lon_0,lat_0,no corners or width/height',
             'eck4'    : 'lon_0,lat_0,no corners or width/height',
             'kav7'    : 'lon_0,lat_0,no corners or width/height',
             'vandg'    : 'lon_0,lat_0,no corners or width/height',
             'mbtfpq'   : 'lon_0,lat_0,no corners or width/height',
             'gnom'     : 'lon_0,lat_0',
             'rotpole'  : 'lon_0,o_lat_p,o_lon_p,corner lat/lon or corner x,y (no width/height)'
             }

# The __init__ docstring is pulled out here because it is so long;
# Having it in the usual place makes it hard to get from the
# __init__ argument list to the code that uses the arguments.
_Basemap_init_doc = """

 Sets up a basemap with specified map projection.
 and creates the coastline data structures in map projection
 coordinates.

 Calling a Basemap class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y map projection coordinates
 (in meters). The inverse transformation is done if the optional keyword
 ``inverse`` is set to True.

 The desired projection is set with the projection keyword. Default is ``cyl``.
 Supported values for the projection keyword are:

 ==============   ====================================================
 Value            Description
 ==============   ====================================================
%(supported_projections)s
 ==============   ====================================================

 For most map projections, the map projection region can either be
 specified by setting these keywords:

 .. tabularcolumns:: |l|L|

 ==============   ====================================================
 Keyword          Description
 ==============   ====================================================
 llcrnrlon        longitude of lower left hand corner of the desired map
                  domain (degrees).
 llcrnrlat        latitude of lower left hand corner of the desired map
                  domain (degrees).
 urcrnrlon        longitude of upper right hand corner of the desired map
                  domain (degrees).
 urcrnrlat        latitude of upper right hand corner of the desired map
                  domain (degrees).
 ==============   ====================================================

 or these

 .. tabularcolumns:: |l|L|

 ==============   ====================================================
 Keyword          Description
 ==============   ====================================================
 width            width of desired map domain in projection coordinates
                  (meters).
 height           height of desired map domain in projection coordinates
                  (meters).
 lon_0            center of desired map domain (in degrees).
 lat_0            center of desired map domain (in degrees).
 ==============   ====================================================

 For ``sinu``, ``moll``, ``hammer``, ``npstere``, ``spstere``, ``nplaea``, ``splaea``,
 ``npaeqd``, ``spaeqd``, ``robin``, ``eck4``, ``kav7``, or ``mbtfpq``, the values of
 llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, width and height are ignored
 (because either they are computed internally, or entire globe is
 always plotted).

 For the cylindrical projections (``cyl``, ``merc``, ``mill``, ``cea``  and ``gall``),
 the default is to use
 llcrnrlon=-180,llcrnrlat=-90, urcrnrlon=180 and urcrnrlat=90). For all other
 projections except ``ortho``, ``geos`` and ``nsper``, either the lat/lon values of the
 corners or width and height must be specified by the user.

 For ``ortho``, ``geos`` and ``nsper``, the lat/lon values of the corners may be specified,
 or the x/y values of the corners (llcrnrx,llcrnry,urcrnrx,urcrnry) in the
 coordinate system of the global projection (with x=0,y=0 at the center
 of the global projection).  If the corners are not specified,
 the entire globe is plotted.

 For ``rotpole``, the lat/lon values of the corners on the unrotated sphere
 may be provided as llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat, or the lat/lon
 values of the corners on the rotated sphere can be given as
 llcrnrx,llcrnry,urcrnrx,urcrnry.

 Other keyword arguments:

 .. tabularcolumns:: |l|L|

 ==============   ====================================================
 Keyword          Description
 ==============   ====================================================
 resolution       resolution of boundary database to use. Can be ``c``
                  (crude), ``l`` (low), ``i`` (intermediate), ``h``
                  (high), ``f`` (full) or None.
                  If None, no boundary data will be read in (and
                  class methods such as drawcoastlines will raise an
                  if invoked).
                  Resolution drops off by roughly 80%% between datasets.
                  Higher res datasets are much slower to draw.
                  Default ``c``. Coastline data is from the GSHHS
                  (http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html).
                  State, country and river datasets from the Generic
                  Mapping Tools (http://gmt.soest.hawaii.edu).
 area_thresh      coastline or lake with an area smaller than
                  area_thresh in km^2 will not be plotted.
                  Default 10000,1000,100,10,1 for resolution
                  ``c``, ``l``, ``i``, ``h``, ``f``.
 rsphere          radius of the sphere used to define map projection
                  (default 6370997 meters, close to the arithmetic mean
                  radius of the earth). If given as a sequence, the
                  first two elements are interpreted as the radii
                  of the major and minor axes of an ellipsoid.
                  Note: sometimes an ellipsoid is specified by the
                  major axis and an inverse flattening parameter (if).
                  The minor axis (b) can be computed from the major
                  axis (a) and the inverse flattening parameter using
                  the formula if = a/(a-b).
 ellps            string describing ellipsoid ('GRS80' or 'WGS84',
                  for example). If both rsphere and ellps are given,
                  rsphere is ignored. Default None. See pyproj.pj_ellps
                  for allowed values.
 suppress_ticks   suppress automatic drawing of axis ticks and labels
                  in map projection coordinates.  Default False,
                  so parallels and meridians can be labelled instead.
                  If parallel or meridian labelling is requested
                  (using drawparallels and drawmeridians methods),
                  automatic tick labelling will be supressed even if
                  suppress_ticks=False.  suppress_ticks=False
                  is useful if you want to use your own custom tick
                  formatter, or  if you want to let matplotlib label
                  the axes in meters using map projection
                  coordinates.
 fix_aspect       fix aspect ratio of plot to match aspect ratio
                  of map projection region (default True).
 anchor           determines how map is placed in axes rectangle
                  (passed to axes.set_aspect). Default is ``C``,
                  which means map is centered.
                  Allowed values are
                  ``C``, ``SW``, ``S``, ``SE``, ``E``, ``NE``,
                  ``N``, ``NW``, and ``W``.
 celestial        use astronomical conventions for longitude (i.e.
                  negative longitudes to the east of 0). Default False.
                  Implies resolution=None.
 ax               set default axes instance
                  (default None - matplotlib.pyplot.gca() may be used
                  to get the current axes instance).
                  If you don not want matplotlib.pyplot to be imported,
                  you can either set this to a pre-defined axes
                  instance, or use the ``ax`` keyword in each Basemap
                  method call that does drawing. In the first case,
                  all Basemap method calls will draw to the same axes
                  instance.  In the second case, you can draw to
                  different axes with the same Basemap instance.
                  You can also use the ``ax`` keyword in individual
                  method calls to selectively override the default
                  axes instance.
 ==============   ====================================================

 The following keywords are map projection parameters which all default to
 None.  Not all parameters are used by all projections, some are ignored.
 The module variable ``projection_params`` is a dictionary which
 lists which parameters apply to which projections.

 .. tabularcolumns:: |l|L|

 ================ ====================================================
 Keyword          Description
 ================ ====================================================
 lat_ts           latitude of true scale. Optional for stereographic,
                  cylindrical equal area and mercator projections.
                  default is lat_0 for stereographic projection.
                  default is 0 for mercator and cylindrical equal area
                  projections.
 lat_1            first standard parallel for lambert conformal,
                  albers equal area and equidistant conic.
                  Latitude of one of the two points on the projection
                  centerline for oblique mercator. If lat_1 is not given, but
                  lat_0 is, lat_1 is set to lat_0 for lambert
                  conformal, albers equal area and equidistant conic.
 lat_2            second standard parallel for lambert conformal,
                  albers equal area and equidistant conic.
                  Latitude of one of the two points on the projection
                  centerline for oblique mercator. If lat_2 is not
                  given it is set to lat_1 for lambert conformal,
                  albers equal area and equidistant conic.
 lon_1            Longitude of one of the two points on the projection
                  centerline for oblique mercator.
 lon_2            Longitude of one of the two points on the projection
                  centerline for oblique mercator.
 k_0              Scale factor at natural origin (used
                  by 'tmerc', 'omerc', 'stere' and 'lcc').
 no_rot           only used by oblique mercator.
                  If set to True, the map projection coordinates will
                  not be rotated to true North.  Default is False
                  (projection coordinates are automatically rotated).
 lat_0            central latitude (y-axis origin) - used by all
                  projections.
 lon_0            central meridian (x-axis origin) - used by all
                  projections.
 boundinglat      bounding latitude for pole-centered projections
                  (npstere,spstere,nplaea,splaea,npaeqd,spaeqd).
                  These projections are square regions centered
                  on the north or south pole.
                  The longitude lon_0 is at 6-o'clock, and the
                  latitude circle boundinglat is tangent to the edge
                  of the map at lon_0.
 round            cut off pole-centered projection at boundinglat
                  (so plot is a circle instead of a square). Only
                  relevant for npstere,spstere,nplaea,splaea,npaeqd
                  or spaeqd projections. Default False.
 satellite_height height of satellite (in m) above equator -
                  only relevant for geostationary
                  and near-sided perspective (``geos`` or ``nsper``)
                  projections. Default 35,786 km.
 ================ ====================================================

 Useful instance variables:

 .. tabularcolumns:: |l|L|

 ================ ====================================================
 Variable Name    Description
 ================ ====================================================
 projection       map projection. Print the module variable
                  ``supported_projections`` to see a list of allowed
                  values.
 epsg             EPSG code defining projection (see
                  http://spatialreference.org for a list of
                  EPSG codes and their definitions).
 aspect           map aspect ratio
                  (size of y dimension / size of x dimension).
 llcrnrlon        longitude of lower left hand corner of the
                  selected map domain.
 llcrnrlat        latitude of lower left hand corner of the
                  selected map domain.
 urcrnrlon        longitude of upper right hand corner of the
                  selected map domain.
 urcrnrlat        latitude of upper right hand corner of the
                  selected map domain.
 llcrnrx          x value of lower left hand corner of the
                  selected map domain in map projection coordinates.
 llcrnry          y value of lower left hand corner of the
                  selected map domain in map projection coordinates.
 urcrnrx          x value of upper right hand corner of the
                  selected map domain in map projection coordinates.
 urcrnry          y value of upper right hand corner of the
                  selected map domain in map projection coordinates.
 rmajor           equatorial radius of ellipsoid used (in meters).
 rminor           polar radius of ellipsoid used (in meters).
 resolution       resolution of boundary dataset being used (``c``
                  for crude, ``l`` for low, etc.).
                  If None, no boundary dataset is associated with the
                  Basemap instance.
 proj4string      the string describing the map projection that is
                  used by PROJ.4.
 ================ ====================================================

 **Converting from Geographic (lon/lat) to Map Projection (x/y) Coordinates**

 Calling a Basemap class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y map projection
 coordinates (in meters).  If optional keyword ``inverse`` is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed.

 For cylindrical equidistant projection (``cyl``), this
 does nothing (i.e. x,y == lon,lat).

 For non-cylindrical projections, the inverse transformation
 always returns longitudes between -180 and 180 degrees. For
 cylindrical projections (self.projection == ``cyl``, ``mill``,
 ``cea``, ``gall`` or ``merc``)
 the inverse transformation will return longitudes between
 self.llcrnrlon and self.llcrnrlat.

 Input arguments lon, lat can be either scalar floats, sequences
 or numpy arrays.

 **Example Usage:**

 >>> from mpl_toolkits.basemap import Basemap
 >>> import numpy as np
 >>> import matplotlib.pyplot as plt
 >>> # read in topo data (on a regular lat/lon grid)
 >>> etopo = np.loadtxt('etopo20data.gz')
 >>> lons  = np.loadtxt('etopo20lons.gz')
 >>> lats  = np.loadtxt('etopo20lats.gz')
 >>> # create Basemap instance for Robinson projection.
 >>> m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
 >>> # compute map projection coordinates for lat/lon grid.
 >>> x, y = m(*np.meshgrid(lons,lats))
 >>> # make filled contour plot.
 >>> cs = m.contourf(x,y,etopo,30,cmap=plt.cm.jet)
 >>> m.drawcoastlines() # draw coastlines
 >>> m.drawmapboundary() # draw a line around the map region
 >>> m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
 >>> m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
 >>> plt.title('Robinson Projection') # add a title
 >>> plt.show()

 [this example (simpletest.py) plus many others can be found in the
 examples directory of source distribution.  The "OO" version of this
 example (which does not use matplotlib.pyplot) is called "simpletest_oo.py".]
""" % locals()

# unsupported projection error message.
_unsupported_projection = ["'%s' is an unsupported projection.\n"]
_unsupported_projection.append("The supported projections are:\n")
_unsupported_projection.append(supported_projections)
_unsupported_projection = ''.join(_unsupported_projection)

def drawparallels(bm_obj,circles,color='k',linewidth=1.,zorder=None, \
                  dashes=[1,1],labels=[0,0,0,0],labelstyle=None, \
                  fmt='%g',xoffset=None,yoffset=None,ax=None,latmax=None,
                  **kwargs):
    """
    Draw and label parallels (latitude lines) for values (in degrees)
    given in the sequence ``circles``.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keyword          Description
    ==============   ====================================================
    color            color to draw parallels (default black).
    linewidth        line width for parallels (default 1.)
    zorder           sets the zorder for parallels (if not specified,
                     uses default zorder for matplotlib.lines.Line2D
                     objects).
    dashes           dash pattern for parallels (default [1,1], i.e.
                     1 pixel on, 1 pixel off).
    labels           list of 4 values (default [0,0,0,0]) that control
                     whether parallels are labelled where they intersect
                     the left, right, top or bottom of the plot. For
                     example labels=[1,0,0,1] will cause parallels
                     to be labelled where they intersect the left and
                     and bottom of the plot, but not the right and top.
    labelstyle       if set to "+/-", north and south latitudes are
                     labelled with "+" and "-", otherwise they are
                     labelled with "N" and "S".
    fmt              a format string to format the parallel labels
                     (default '%g') **or** a function that takes a
                     latitude value in degrees as it's only argument
                     and returns a formatted string.
    xoffset          label offset from edge of map in x-direction
                     (default is 0.01 times width of map in map
                     projection coordinates).
    yoffset          label offset from edge of map in y-direction
                     (default is 0.01 times height of map in map
                     projection coordinates).
    ax               axes instance (overrides default axes instance)
    latmax           absolute value of latitude to which meridians are drawn
                     (default is 80).
    \**kwargs        additional keyword arguments controlling text
                     for labels that are passed on to
                     the text method of the axes instance (see
                     matplotlib.pyplot.text documentation).
    ==============   ====================================================

    returns a dictionary whose keys are the parallel values, and
    whose values are tuples containing lists of the
    matplotlib.lines.Line2D and matplotlib.text.Text instances
    associated with each parallel. Deleting an item from the
    dictionary removes the corresponding parallel from the plot.
    """
    # if celestial=True, don't use "N" and "S" labels.
    if labelstyle is None and bm_obj.celestial:
        labelstyle="+/-"
    # get current axes instance (if none specified).
    ax = ax or bm_obj._check_ax()
    # don't draw meridians past latmax, always draw parallel at latmax.
    if latmax is None: latmax = 80.
    # offset for labels.
    if yoffset is None:
        yoffset = (bm_obj.urcrnry-bm_obj.llcrnry)/100.
        if bm_obj.aspect > 1:
            yoffset = bm_obj.aspect*yoffset
        else:
            yoffset = yoffset/bm_obj.aspect
    if xoffset is None:
        xoffset = (bm_obj.urcrnrx-bm_obj.llcrnrx)/100.

    if bm_obj.projection in _cylproj + _pseudocyl:
        lons = np.linspace(bm_obj.llcrnrlon, bm_obj.urcrnrlon, 10001)
    elif bm_obj.projection in ['tmerc']:
        lon_0 = bm_obj.projparams['lon_0']
        # tmerc only defined within +/- 90 degrees of lon_0
        lons = np.linspace(lon_0-90,lon_0+90,100001)
    else:
        lonmin = bm_obj.boundarylonmin; lonmax = bm_obj.boundarylonmax
        lons = np.linspace(lonmin, lonmax, 10001)
    # make sure latmax degree parallel is drawn if projection not merc or cyl or miller
    try:
        circlesl = list(circles)
    except:
        circlesl = circles
    if bm_obj.projection not in _cylproj + _pseudocyl:
        if max(circlesl) > 0 and latmax not in circlesl:
            circlesl.append(latmax)
        if min(circlesl) < 0 and -latmax not in circlesl:
            circlesl.append(-latmax)
    xdelta = 0.01*(bm_obj.xmax-bm_obj.xmin)
    ydelta = 0.01*(bm_obj.ymax-bm_obj.ymin)
    linecolls = {}
    for circ in circlesl:
        lats = circ*np.ones(len(lons),np.float32)
        x,y = bm_obj(lons,lats)
        # remove points outside domain.
        # leave a little slop around edges (3*xdelta)
        # don't really know why, but this appears to be needed to
        # or lines sometimes don't reach edge of plot.
        testx = np.logical_and(x>=bm_obj.xmin-3*xdelta,x<=bm_obj.xmax+3*xdelta)
        x = np.compress(testx, x)
        y = np.compress(testx, y)
        testy = np.logical_and(y>=bm_obj.ymin-3*ydelta,y<=bm_obj.ymax+3*ydelta)
        x = np.compress(testy, x)
        y = np.compress(testy, y)
        lines = []
        if len(x) > 1 and len(y) > 1:
            # split into separate line segments if necessary.
            # (not necessary for cylindrical or pseudocylindricl projections)
            xd = (x[1:]-x[0:-1])**2
            yd = (y[1:]-y[0:-1])**2
            dist = np.sqrt(xd+yd)
            if bm_obj.projection not in ['cyl','rotpole']:
                split = dist > bm_obj.rmajor/10.
            else:
                split = dist > 1.
            if np.sum(split) and bm_obj.projection not in _cylproj:
                ind = (np.compress(split,np.squeeze(split*np.indices(xd.shape)))+1).tolist()
                xl = []
                yl = []
                iprev = 0
                ind.append(len(xd))
                for i in ind:
                    xl.append(x[iprev:i])
                    yl.append(y[iprev:i])
                    iprev = i
            else:
                xl = [x]
                yl = [y]
            # draw each line segment.
            for x,y in zip(xl,yl):
                # skip if only a point.
                if len(x) > 1 and len(y) > 1:
                    l = Line2D(x,y,linewidth=linewidth)
                    l.set_color(color)
                    l.set_dashes(dashes)
                    l.set_label('_nolabel_')
                    if zorder is not None:
                        l.set_zorder(zorder)
                    ax.add_line(l)
                    lines.append(l)
        linecolls[circ] = (lines,[])
    # draw labels for parallels
    # parallels not labelled for fulldisk orthographic or geostationary
    if bm_obj.projection in ['ortho','geos','nsper','vandg','aeqd'] and max(labels):
        if bm_obj.projection == 'vandg' or bm_obj._fulldisk:
            sys.stdout.write('Warning: Cannot label parallels on %s basemap' % _projnames[bm_obj.projection])
            labels = [0,0,0,0]
    # search along edges of map to see if parallels intersect.
    # if so, find x,y location of intersection and draw a label there.
    dx = (bm_obj.xmax-bm_obj.xmin)/1000.
    dy = (bm_obj.ymax-bm_obj.ymin)/1000.
    if bm_obj.projection in _pseudocyl:
        lon_0 = bm_obj.projparams['lon_0']
    for dolab,side in zip(labels,['l','r','t','b']):
        if not dolab: continue
        # for cylindrical projections, don't draw parallels on top or bottom.
        if bm_obj.projection in _cylproj + _pseudocyl and side in ['t','b']: continue
        if side in ['l','r']:
            nmax = int((bm_obj.ymax-bm_obj.ymin)/dy+1)
            yy = np.linspace(bm_obj.llcrnry,bm_obj.urcrnry,nmax)
            if side == 'l':
                if bm_obj.projection in _pseudocyl:
                    lats = np.linspace(-89.99,89,99,nmax)
                    if bm_obj.celestial:
                        lons = (bm_obj.projparams['lon_0']+180.)*np.ones(len(lats),lats.dtype)
                    else:
                        lons = (bm_obj.projparams['lon_0']-180.)*np.ones(len(lats),lats.dtype)
                    xx, yy = bm_obj(lons, lats)
                else:
                    xx = bm_obj.llcrnrx*np.ones(yy.shape,yy.dtype)
                    lons,lats = bm_obj(xx,yy,inverse=True)
                    # Mask the missing values, and fill with zero
                    lons = ma.masked_equal(lons,1e+30)
                    lats = ma.masked_array(lats,mask=lons.mask)
                    lons = ma.filled(lons,fill_value=0.)
                    lats = ma.filled(lats,fill_value=0.)
                    lons = lons.tolist(); lats = lats.tolist()
            else:
                if bm_obj.projection in _pseudocyl:
                    lats = np.linspace(-89.99,89,99,nmax)
                    if bm_obj.celestial:
                       lons = (bm_obj.projparams['lon_0']-180.)*np.ones(len(lats),lats.dtype)
                    else:
                       lons = (bm_obj.projparams['lon_0']+180.)*np.ones(len(lats),lats.dtype)
                    xx, yy = bm_obj(lons, lats)
                else:
                    xx = bm_obj.urcrnrx*np.ones(yy.shape,yy.dtype)
                    lons,lats = bm_obj(xx,yy,inverse=True)
                    # Mask the missing values, and fill with zero
                    lons = ma.masked_equal(lons,1e+30)
                    lats = ma.masked_array(lats,mask=lons.mask)
                    lons = ma.filled(lons,fill_value=0.)
                    lats = ma.filled(lats,fill_value=0.)
                    lons = lons.tolist(); lats = lats.tolist()
            if max(lons) > 1.e20 or max(lats) > 1.e20:
                raise ValueError('inverse transformation undefined - please adjust the map projection region')
            # adjust so 0 <= lons < 360
            lons = [(lon+360) % 360 for lon in lons]
        else:
            nmax = int((bm_obj.xmax-bm_obj.xmin)/dx+1)
            xx = np.linspace(bm_obj.llcrnrx,bm_obj.urcrnrx,nmax)
            if side == 'b':
                lons,lats = bm_obj(xx,bm_obj.llcrnry*np.ones(xx.shape,np.float32),inverse=True)
                # Mask the missing values, and fill with zero
                lons = ma.masked_equal(lons,1e+30)
                lats = ma.masked_array(lats,mask=lons.mask)
                lons = ma.filled(lons,fill_value=0.)
                lats = ma.filled(lats,fill_value=0.)
                lons = lons.tolist(); lats = lats.tolist()
            else:
                lons,lats = bm_obj(xx,bm_obj.urcrnry*np.ones(xx.shape,np.float32),inverse=True)
                # Mask the missing values, and fill with zero
                lons = ma.masked_equal(lons,1e+30)
                lats = ma.masked_array(lats,mask=lons.mask)
                lons = ma.filled(lons,fill_value=0.)
                lats = ma.filled(lats,fill_value=0.)
                lons = lons.tolist(); lats = lats.tolist()
            if max(lons) > 1.e20 or max(lats) > 1.e20:
                raise ValueError('inverse transformation undefined - please adjust the map projection region')
            # adjust so 0 <= lons < 360
            lons = [(lon+360) % 360 for lon in lons]
        for lat in circles:
            # don't label parallels for round polar plots
            if bm_obj.round: continue
            # find index of parallel (there may be two, so
            # search from left and right).
            nl = _searchlist(lats,lat)
            nr = _searchlist(lats[::-1],lat)
            if nr != -1: nr = len(lons)-nr-1
            latlab = _setlatlab(fmt,lat,labelstyle)
            # parallels can intersect each map edge twice.
            for i,n in enumerate([nl,nr]):
                # don't bother if close to the first label.
                if i and abs(nr-nl) < 100: continue
                if n >= 0:
                    t = None
                    if side == 'l':
                        if bm_obj.projection in _pseudocyl:
                            if bm_obj.celestial:
                                xlab,ylab = bm_obj(lon_0+179.9,lat)
                            else:
                                xlab,ylab = bm_obj(lon_0-179.9,lat)
                        else:
                            xlab = bm_obj.llcrnrx
                        xlab = xlab-xoffset
                        if bm_obj.projection in _pseudocyl:
                            if lat>0:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='bottom',**kwargs)
                            elif lat<0:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='top',**kwargs)
                            else:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='center',**kwargs)
                        else:
                           t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='center',**kwargs)
                        #print "Parallel side {}, xlab = {}, yy[{}] = {}".format(side,xlab,n,yy[n])
                    elif side == 'r':
                        if bm_obj.projection in _pseudocyl:
                            if bm_obj.celestial:
                               xlab,ylab = bm_obj(lon_0-179.9,lat)
                            else:
                               xlab,ylab = bm_obj(lon_0+179.9,lat)
                        else:
                            xlab = bm_obj.urcrnrx
                        xlab = xlab+xoffset
                        if bm_obj.projection in _pseudocyl:
                            if lat>0:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='bottom',**kwargs)
                            elif lat<0:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='top',**kwargs)
                            else:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='center',**kwargs)
                        else:
                           t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='center',**kwargs)
                        #print "Parallel side {}, xlab = {}, yy[{}] = {}".format(side,xlab,n,yy[n])
                    elif side == 'b':
                        t = ax.text(xx[n],bm_obj.llcrnry-yoffset,latlab,horizontalalignment='center',verticalalignment='top',**kwargs)
                        #print "Parallel side {}, xlab = {}, yy[{}] = {}".format(side,xlab,n,yy[n])
                    else:
                        t = ax.text(xx[n],bm_obj.urcrnry+yoffset,latlab,horizontalalignment='center',verticalalignment='bottom',**kwargs)
                        #print "Parallel side {}, xlab = {}, yy[{}] = {}".format(side,xlab,n,yy[n])
                    if t is not None: linecolls[lat][1].append(t)

    # set axes limits to fit map region.
    bm_obj.set_axes_limits(ax=ax)
    keys = list(linecolls.keys()); vals = list(linecolls.values())
    for k,v in zip(keys,vals):
        if v == ([], []):
            del linecolls[k]
        # add a remove method to each tuple.
        else:
            linecolls[k] = _tup(linecolls[k])
    # override __delitem__ in dict to call remove() on values.
    pardict = _dict(linecolls)
    # clip parallels for round polar plots (and delete labels).
    if bm_obj.round:
        c = Circle((0.5*(bm_obj.xmax+bm_obj.xmin),0.5*(bm_obj.ymax+bm_obj.ymin)),
            radius=0.5*(bm_obj.xmax-bm_obj.xmin),fc='none')
        if c not in ax.patches:
            p = ax.add_patch(c)
            p.set_clip_on(False)
        for par in pardict:
            lines,labs = pardict[par]
            for l in lines:
                l.set_clip_path(c)
    return pardict

def drawmeridians(self,meridians,color='k',linewidth=1., zorder=None,\
                  dashes=[1,1],labels=[0,0,0,0],labelstyle=None,\
                  fmt='%g',xoffset=None,yoffset=None,ax=None,latmax=None,
                  **kwargs):
    """
    Draw and label meridians (longitude lines) for values (in degrees)
    given in the sequence ``meridians``.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keyword          Description
    ==============   ====================================================
    color            color to draw meridians (default black).
    linewidth        line width for meridians (default 1.)
    zorder           sets the zorder for meridians (if not specified,
                     uses default zorder for matplotlib.lines.Line2D
                     objects).
    dashes           dash pattern for meridians (default [1,1], i.e.
                     1 pixel on, 1 pixel off).
    labels           list of 4 values (default [0,0,0,0]) that control
                     whether meridians are labelled where they intersect
                     the left, right, top or bottom of the plot. For
                     example labels=[1,0,0,1] will cause meridians
                     to be labelled where they intersect the left and
                     and bottom of the plot, but not the right and top.
    labelstyle       if set to "+/-", east and west longitudes are
                     labelled with "+" and "-", otherwise they are
                     labelled with "E" and "W".
    fmt              a format string to format the meridian labels
                     (default '%g') **or** a function that takes a
                     longitude value in degrees as it's only argument
                     and returns a formatted string.
    xoffset          label offset from edge of map in x-direction
                     (default is 0.01 times width of map in map
                     projection coordinates).
    yoffset          label offset from edge of map in y-direction
                     (default is 0.01 times height of map in map
                     projection coordinates).
    ax               axes instance (overrides default axes instance)
    latmax           absolute value of latitude to which meridians are drawn
                     (default is 80).
    \**kwargs        additional keyword arguments controlling text
                     for labels that are passed on to
                     the text method of the axes instance (see
                     matplotlib.pyplot.text documentation).
    ==============   ====================================================

    returns a dictionary whose keys are the meridian values, and
    whose values are tuples containing lists of the
    matplotlib.lines.Line2D and matplotlib.text.Text instances
    associated with each meridian. Deleting an item from the
    dictionary removes the correpsonding meridian from the plot.
    """
    # for cylindrical projections, try to handle wraparound (i.e. if
    # projection is defined in -180 to 0 and user asks for meridians from
    # 180 to 360 to be drawn, it should work)
    if self.projection in _cylproj or self.projection in _pseudocyl:
        def addlon(meridians,madd):
            minside = (madd >= self.llcrnrlon and madd <= self.urcrnrlon)
            if minside and madd not in meridians: meridians.append(madd)
            return meridians
        merids = list(meridians)
        meridians = []
        for m in merids:
            meridians = addlon(meridians,m)
            meridians = addlon(meridians,m+360)
            meridians = addlon(meridians,m-360)
        meridians.sort()
    # if celestial=True, don't use "E" and "W" labels.
    if labelstyle is None and self.celestial:
        labelstyle="+/-"
    # get current axes instance (if none specified).
    ax = ax or self._check_ax()
    # don't draw meridians past latmax, always draw parallel at latmax.
    if latmax is None: latmax = 80. # unused w/ cyl, merc or miller proj.
    # offset for labels.
    if yoffset is None:
        yoffset = (self.urcrnry-self.llcrnry)/100.
        if self.aspect > 1:
            yoffset = self.aspect*yoffset
        else:
            yoffset = yoffset/self.aspect
    if xoffset is None:
        xoffset = (self.urcrnrx-self.llcrnrx)/100.

    lats = np.linspace(self.latmin,self.latmax,10001)
    if self.projection not in _cylproj + _pseudocyl:
        testlat = np.logical_and(lats>-latmax,lats<latmax)
        lats = np.compress(testlat,lats)

    xdelta = 0.01*(self.xmax-self.xmin)
    ydelta = 0.01*(self.ymax-self.ymin)
    linecolls = {}
    for merid in meridians:
        lons = merid*np.ones(len(lats),np.float32)
        x,y = self(lons,lats)
        # remove points outside domain.
        # leave a little slop around edges (3*xdelta)
        # don't really know why, but this appears to be needed to
        # or lines sometimes don't reach edge of plot.
        testx = np.logical_and(x>=self.xmin-3*xdelta,x<=self.xmax+3*xdelta)
        x = np.compress(testx, x)
        y = np.compress(testx, y)
        testy = np.logical_and(y>=self.ymin-3*ydelta,y<=self.ymax+3*ydelta)
        x = np.compress(testy, x)
        y = np.compress(testy, y)
        lines = []
        if len(x) > 1 and len(y) > 1:
            # split into separate line segments if necessary.
            # (not necessary for mercator or cylindrical or miller).
            xd = (x[1:]-x[0:-1])**2
            yd = (y[1:]-y[0:-1])**2
            dist = np.sqrt(xd+yd)
            if self.projection not in ['cyl','rotpole']:
                split = dist > self.rmajor/10.
            else:
                split = dist > 1.
            if np.sum(split) and self.projection not in _cylproj:
                ind = (np.compress(split,np.squeeze(split*np.indices(xd.shape)))+1).tolist()
                xl = []
                yl = []
                iprev = 0
                ind.append(len(xd))
                for i in ind:
                    xl.append(x[iprev:i])
                    yl.append(y[iprev:i])
                    iprev = i
            else:
                xl = [x]
                yl = [y]
            # draw each line segment.
            for x,y in zip(xl,yl):
                # skip if only a point.
                if len(x) > 1 and len(y) > 1:
                    l = Line2D(x,y,linewidth=linewidth)
                    l.set_color(color)
                    l.set_dashes(dashes)
                    l.set_label('_nolabel_')
                    if zorder is not None:
                        l.set_zorder(zorder)
                    ax.add_line(l)
                    lines.append(l)
        linecolls[merid] = (lines,[])
    # draw labels for meridians.
    # meridians not labelled for sinusoidal, hammer, mollweide,
    # VanDerGrinten or full-disk orthographic/geostationary.
    if self.projection in ['sinu','moll','hammer','vandg'] and max(labels):
        sys.stdout.write('Warning: Cannot label meridians on %s basemap' % _projnames[self.projection])
        labels = [0,0,0,0]
    if self.projection in ['ortho','geos','nsper','aeqd'] and max(labels):
        if self._fulldisk and self.boundinglat is None:
            sys.stdout.write(dedent(
            """'Warning: Cannot label meridians on full-disk
            Geostationary, Orthographic or Azimuthal equidistant basemap
            """))
            labels = [0,0,0,0]
    # search along edges of map to see if parallels intersect.
    # if so, find x,y location of intersection and draw a label there.
    dx = (self.xmax-self.xmin)/1000.
    dy = (self.ymax-self.ymin)/1000.
    if self.projection in _pseudocyl:
        lon_0 = self.projparams['lon_0']
        xmin,ymin = self(lon_0-179.9,-90)
        xmax,ymax = self(lon_0+179.9,90)
    for dolab,side in zip(labels,['l','r','t','b']):
        #print "Meridian (side {}):".format(side)
        if not dolab or self.round: continue
        # for cylindrical projections, don't draw meridians on left or right.
        if self.projection in _cylproj + _pseudocyl and side in ['l','r']: continue
        if side in ['l','r']:
            nmax = int((self.ymax-self.ymin)/dy+1)
            yy = np.linspace(self.llcrnry,self.urcrnry,nmax)
            if side == 'l':
                lons,lats = self(self.llcrnrx*np.ones(yy.shape,np.float32),yy,inverse=True)
                # Mask the missing values, and fill with zero
                lons = ma.masked_equal(lons,1e+30)
                lats = ma.masked_array(lats,mask=lons.mask)
                lons = ma.filled(lons,fill_value=0.)
                lats = ma.filled(lats,fill_value=0.)
                lons = lons.tolist(); lats = lats.tolist()
            else:
                lons,lats = self(self.urcrnrx*np.ones(yy.shape,np.float32),yy,inverse=True)
                # Mask the missing values, and fill with zero
                lons = ma.masked_equal(lons,1e+30)
                lats = ma.masked_array(lats,mask=lons.mask)
                lons = ma.filled(lons,fill_value=0.)
                lats = ma.filled(lats,fill_value=0.)
                lons = lons.tolist(); lats = lats.tolist()
            if max(lons) > 1.e20 or max(lats) > 1.e20:
                raise ValueError('inverse transformation undefined - please adjust the map projection region')
            # adjust so 0 <= lons < 360
            lons = [(lon+360) % 360 for lon in lons]
        else:
            nmax = int((self.xmax-self.xmin)/dx+1)
            if self.projection in _pseudocyl:
                xx = np.linspace(xmin,xmax,nmax)
            else:
                xx = np.linspace(self.llcrnrx,self.urcrnrx,nmax)
            if side == 'b':
                lons,lats = self(xx,self.llcrnry*np.ones(xx.shape,np.float32),inverse=True)
                # Mask the missing values, and fill with zero
                lons = ma.masked_equal(lons,1e+30)
                lats = ma.masked_array(lats,mask=lons.mask)
                lons = ma.filled(lons,fill_value=0.)
                lats = ma.filled(lats,fill_value=0.)
                lons = lons.tolist(); lats = lats.tolist()
            else:
                lons,lats = self(xx,self.urcrnry*np.ones(xx.shape,np.float32),inverse=True)
                # Mask the missing values, and fill with zero
                lons = ma.masked_equal(lons,1e+30)
                lats = ma.masked_array(lats,mask=lons.mask)
                lons = ma.filled(lons,fill_value=0.)
                lats = ma.filled(lats,fill_value=0.)
                lons = lons.tolist(); lats = lats.tolist()
            if max(lons) > 1.e20 or max(lats) > 1.e20:
                raise ValueError('inverse transformation undefined - please adjust the map projection region')
            # adjust so 0 <= lons < 360
            #print " ".join(["{:7.2f}".format(x) for x in lons])
            lons = [(lon+360) % 360 for lon in lons]
            #print " ".join(["{:7.2f}".format(x) for x in lons])
        for lon in meridians:
            # adjust so 0 <= lon < 360
            lon2 = (lon+360) % 360
            #print "\tlon,lon2 = {},{}".format(lon,lon2)
            # find index of meridian (there may be two, so
            # search from left and right).
            nl = _searchlist(lons,lon2)
            nr = _searchlist(lons[::-1],lon2)
            #print "\t\tnl, nr = {},{}".format(nl,nr)
            if nr != -1: nr = len(lons)-nr-1
            lonlab = _setlonlab(fmt,lon2,labelstyle)
            # meridians can intersect each map edge twice.
            for i,n in enumerate([nl,nr]):
                lat = lats[n]/100.
                # no meridians > latmax for projections other than merc,cyl,miller.
                if self.projection not in _cylproj and lat > latmax: continue
                # don't bother if close to the first label.
                if i and abs(nr-nl) < 100: continue
                if n >= 0:
                    t = None
                    if side == 'l':
                        t = ax.text(self.llcrnrx-xoffset,yy[n],lonlab,horizontalalignment='right',verticalalignment='center',**kwargs)
                        #print "\t\tMeridian side {}, xlab = {}, yy[{}] = {}".format(side,self.llcrnrx-xoffset,yy[n])
                    elif side == 'r':
                        t = ax.text(self.urcrnrx+xoffset,yy[n],lonlab,horizontalalignment='left',verticalalignment='center',**kwargs)
                        #print "\t\tMeridian side {}, xlab = {}, yy[{}] = {}".format(side,self.urcrnrx+xoffset,yy[n])
                    elif side == 'b':
                        t = ax.text(xx[n],self.llcrnry-yoffset,lonlab,horizontalalignment='center',verticalalignment='top',**kwargs)
                        #print "\t\tMeridian side {}, xx[{}] = {}, ylab = {}".format(side,n,xx[n],self.llcrnry-yoffset)
                    else:
                        t = ax.text(xx[n],self.urcrnry+yoffset,lonlab,horizontalalignment='center',verticalalignment='bottom',**kwargs)
                        #print "\t\tMeridian side {}, xx[{}] = {}, ylab = {}".format(side,n,xx[n],self.urcrnry+yoffset)

                    if t is not None: linecolls[lon][1].append(t)
    # set axes limits to fit map region.
    self.set_axes_limits(ax=ax)
    # remove empty values from linecolls dictionary
    keys = list(linecolls.keys()); vals = list(linecolls.values())
    for k,v in zip(keys,vals):
        if v == ([], []):
            del linecolls[k]
        else:
        # add a remove method to each tuple.
            linecolls[k] = _tup(linecolls[k])
    # override __delitem__ in dict to call remove() on values.
    meridict = _dict(linecolls)
    # for round polar plots, clip meridian lines and label them.
    if self.round:
        c = Circle((0.5*(self.xmax+self.xmin),0.5*(self.ymax+self.ymin)),
            radius=0.5*(self.xmax-self.xmin),fc='none')
        if c not in ax.patches:
            p = ax.add_patch(c)
            p.set_clip_on(False)
        # label desired?
        label = False
        for lab in labels:
            if lab: label = True
        for merid in meridict:
            lines,labs = meridict[merid]
            # clip lines.
            for l in lines:
                l.set_clip_path(c)
            if not label: continue
            # label
            lonlab = _setlonlab(fmt,merid,labelstyle)
            x,y = self(merid,self.boundinglat)
            r = np.sqrt((x-0.5*(self.xmin+self.xmax))**2+
                        (y-0.5*(self.ymin+self.ymax))**2)
            r = r + np.sqrt(xoffset**2+yoffset**2)
            if self.projection.startswith('np'):
                pole = 1
            elif self.projection.startswith('sp'):
                pole = -1
            elif self.projection == 'ortho' and self.round:
                pole = 1
            if pole == 1:
                theta = (np.pi/180.)*(merid-self.projparams['lon_0']-90)
                if self.projection == 'ortho' and\
                   self.projparams['lat_0'] == -90:
                    theta = (np.pi/180.)*(-merid+self.projparams['lon_0']+90)
                x = r*np.cos(theta)+0.5*(self.xmin+self.xmax)
                y = r*np.sin(theta)+0.5*(self.ymin+self.ymax)
                if x > 0.5*(self.xmin+self.xmax)+xoffset:
                    horizalign = 'left'
                elif x < 0.5*(self.xmin+self.xmax)-xoffset:
                    horizalign = 'right'
                else:
                    horizalign = 'center'
                if y > 0.5*(self.ymin+self.ymax)+yoffset:
                    vertalign = 'bottom'
                elif y < 0.5*(self.ymin+self.ymax)-yoffset:
                    vertalign = 'top'
                else:
                    vertalign = 'center'
                # labels [l,r,t,b]
                if labels[0] and not labels[1] and x >= 0.5*(self.xmin+self.xmax)+xoffset: continue
                if labels[1] and not labels[0] and x <= 0.5*(self.xmin+self.xmax)-xoffset: continue
                if labels[2] and not labels[3] and y <= 0.5*(self.ymin+self.ymax)-yoffset: continue
                if labels[3] and not labels[2]and y >= 0.5*(self.ymin+self.ymax)+yoffset: continue
            elif pole == -1:
                theta = (np.pi/180.)*(-merid+self.projparams['lon_0']+90)
                x = r*np.cos(theta)+0.5*(self.xmin+self.xmax)
                y = r*np.sin(theta)+0.5*(self.ymin+self.ymax)
                if x > 0.5*(self.xmin+self.xmax)-xoffset:
                    horizalign = 'right'
                elif x < 0.5*(self.xmin+self.xmax)+xoffset:
                    horizalign = 'left'
                else:
                    horizalign = 'center'
                if y > 0.5*(self.ymin+self.ymax)-yoffset:
                    vertalign = 'top'
                elif y < 0.5*(self.ymin+self.ymax)+yoffset:
                    vertalign = 'bottom'
                else:
                    vertalign = 'center'
                # labels [l,r,t,b]
                if labels[0] and not labels[1] and x <=  0.5*(self.xmin+self.xmax)+xoffset: continue
                if labels[1] and not labels[0] and x >=  0.5*(self.xmin+self.xmax)-xoffset: continue
                if labels[2] and not labels[3] and y >=  0.5*(self.ymin+self.ymax)-yoffset: continue
                if labels[3] and not labels[2] and y <=  0.5*(self.ymin+self.ymax)+yoffset: continue
            t =\
            ax.text(x,y,lonlab,horizontalalignment=horizalign,verticalalignment=vertalign,**kwargs)
            meridict[merid][1].append(t)
    return meridict

def _check_ax(self):
    """
    Returns the axis on which to draw.
    Returns self.ax, or if self.ax=None returns plt.gca().
    """
    if self.ax is None:
        try:
            ax = plt.gca()
        except:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        # associate an axes instance with this Basemap instance
        # the first time this method is called.
        #self.ax = ax
    else:
        ax = self.ax
    return ax


def _searchlist(a,x):
    """
    like bisect, but works for lists that are not sorted,
    and are not in increasing order.
    returns -1 if x does not fall between any two elements"""
    # make sure x is a float (and not an array scalar)
    x = float(x)
    itemprev = a[0]
    nslot = -1
    eps = 180.
    for n,item in enumerate(a[1:]):
        if item < itemprev:
            if itemprev-item>eps:
                if ((x>itemprev and x<=360.) or (x<item and x>=0.)):
                    nslot = n+1
                    break
            elif x <= itemprev and x > item and itemprev:
                nslot = n+1
                break
        else:
            if item-itemprev>eps:
                if ((x<itemprev and x>=0.) or (x>item and x<=360.)):
                    nslot = n+1
                    break
            elif x >= itemprev and x < item:
                nslot = n+1
                break
        itemprev = item
    return nslot

def interp(datain,xin,yin,xout,yout,checkbounds=False,masked=False,order=1):
    """
    Interpolate data (``datain``) on a rectilinear grid (with x = ``xin``
    y = ``yin``) to a grid with x = ``xout``, y= ``yout``.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Arguments        Description
    ==============   ====================================================
    datain           a rank-2 array with 1st dimension corresponding to
                     y, 2nd dimension x.
    xin, yin         rank-1 arrays containing x and y of
                     datain grid in increasing order.
    xout, yout       rank-2 arrays containing x and y of desired output grid.
    ==============   ====================================================

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    checkbounds      If True, values of xout and yout are checked to see
                     that they lie within the range specified by xin
                     and xin.
                     If False, and xout,yout are outside xin,yin,
                     interpolated values will be clipped to values on
                     boundary of input grid (xin,yin)
                     Default is False.
    masked           If True, points outside the range of xin and yin
                     are masked (in a masked array).
                     If masked is set to a number, then
                     points outside the range of xin and yin will be
                     set to that number. Default False.
    order            0 for nearest-neighbor interpolation, 1 for
                     bilinear interpolation, 3 for cublic spline
                     (default 1). order=3 requires scipy.ndimage.
    ==============   ====================================================

    .. note::
     If datain is a masked array and order=1 (bilinear interpolation) is
     used, elements of dataout will be masked if any of the four surrounding
     points in datain are masked.  To avoid this, do the interpolation in two
     passes, first with order=1 (producing dataout1), then with order=0
     (producing dataout2).  Then replace all the masked values in dataout1
     with the corresponding elements in dataout2 (using numpy.where).
     This effectively uses nearest neighbor interpolation if any of the
     four surrounding points in datain are masked, and bilinear interpolation
     otherwise.

    Returns ``dataout``, the interpolated data on the grid ``xout, yout``.
    """
    # xin and yin must be monotonically increasing.
    if xin[-1]-xin[0] < 0 or yin[-1]-yin[0] < 0:
        raise ValueError('xin and yin must be increasing!')
    if xout.shape != yout.shape:
        raise ValueError('xout and yout must have same shape!')
    # check that xout,yout are
    # within region defined by xin,yin.
    if checkbounds:
        if xout.min() < xin.min() or \
           xout.max() > xin.max() or \
           yout.min() < yin.min() or \
           yout.max() > yin.max():
            raise ValueError('yout or xout outside range of yin or xin')
    # compute grid coordinates of output grid.
    delx = xin[1:]-xin[0:-1]
    dely = yin[1:]-yin[0:-1]
    if max(delx)-min(delx) < 1.e-4 and max(dely)-min(dely) < 1.e-4:
        # regular input grid.
        xcoords = (len(xin)-1)*(xout-xin[0])/(xin[-1]-xin[0])
        ycoords = (len(yin)-1)*(yout-yin[0])/(yin[-1]-yin[0])
    else:
        # irregular (but still rectilinear) input grid.
        xoutflat = xout.flatten(); youtflat = yout.flatten()
        ix = (np.searchsorted(xin,xoutflat)-1).tolist()
        iy = (np.searchsorted(yin,youtflat)-1).tolist()
        xoutflat = xoutflat.tolist(); xin = xin.tolist()
        youtflat = youtflat.tolist(); yin = yin.tolist()
        xcoords = []; ycoords = []
        for n,i in enumerate(ix):
            if i < 0:
                xcoords.append(-1) # outside of range on xin (lower end)
            elif i >= len(xin)-1:
                xcoords.append(len(xin)) # outside range on upper end.
            else:
                xcoords.append(float(i)+(xoutflat[n]-xin[i])/(xin[i+1]-xin[i]))
        for m,j in enumerate(iy):
            if j < 0:
                ycoords.append(-1) # outside of range of yin (on lower end)
            elif j >= len(yin)-1:
                ycoords.append(len(yin)) # outside range on upper end
            else:
                ycoords.append(float(j)+(youtflat[m]-yin[j])/(yin[j+1]-yin[j]))
        xcoords = np.reshape(xcoords,xout.shape)
        ycoords = np.reshape(ycoords,yout.shape)
    # data outside range xin,yin will be clipped to
    # values on boundary.
    if masked:
        xmask = np.logical_or(np.less(xcoords,0),np.greater(xcoords,len(xin)-1))
        ymask = np.logical_or(np.less(ycoords,0),np.greater(ycoords,len(yin)-1))
        xymask = np.logical_or(xmask,ymask)
    xcoords = np.clip(xcoords,0,len(xin)-1)
    ycoords = np.clip(ycoords,0,len(yin)-1)
    # interpolate to output grid using bilinear interpolation.
    if order == 1:
        xi = xcoords.astype(np.int32)
        yi = ycoords.astype(np.int32)
        xip1 = xi+1
        yip1 = yi+1
        xip1 = np.clip(xip1,0,len(xin)-1)
        yip1 = np.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype(np.float32)
        dely = ycoords-yi.astype(np.float32)
        dataout = (1.-delx)*(1.-dely)*datain[yi,xi] + \
                  delx*dely*datain[yip1,xip1] + \
                  (1.-delx)*dely*datain[yip1,xi] + \
                  delx*(1.-dely)*datain[yi,xip1]
    elif order == 0:
        xcoordsi = np.around(xcoords).astype(np.int32)
        ycoordsi = np.around(ycoords).astype(np.int32)
        dataout = datain[ycoordsi,xcoordsi]
    elif order == 3:
        try:
            from scipy.ndimage import map_coordinates
        except ImportError:
            raise ValueError('scipy.ndimage must be installed if order=3')
        coords = [ycoords,xcoords]
        dataout = map_coordinates(datain,coords,order=3,mode='nearest')
    else:
        raise ValueError('order keyword must be 0, 1 or 3')
    if masked and isinstance(masked,bool):
        dataout = ma.masked_array(dataout)
        newmask = ma.mask_or(ma.getmask(dataout), xymask)
        dataout = ma.masked_array(dataout,mask=newmask)
    elif masked and is_scalar(masked):
        dataout = np.where(xymask,masked,dataout)
    return dataout

class _tup(tuple):
    # tuple with an added remove method.
    # used for objects returned by drawparallels and drawmeridians.
    def remove(self):
        for item in self:
            for x in item:
                x.remove()
class _dict(dict):
    # override __delitem__ to first call remove method on values.
    def __delitem__(self,key):
        self[key].remove()
        super(_dict, self).__delitem__(key)

def _setlonlab(fmt,lon,labelstyle):
    # set lon label string (called by Basemap.drawmeridians)
    try: # fmt is a function that returns a formatted string
        lonlab = fmt(lon)
    except: # fmt is a format string.
        if lon>180:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    lonlabstr = r'${\/-%s\/^{\circ}}$'%fmt
                else:
                    lonlabstr = r'${%s\/^{\circ}\/W}$'%fmt
            else:
                if labelstyle=='+/-':
                    lonlabstr = u'-%s\N{DEGREE SIGN}'%fmt
                else:
                    lonlabstr = u'%s\N{DEGREE SIGN}W'%fmt
            lonlab = lonlabstr%np.fabs(lon-360)
        elif lon<180 and lon != 0:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    lonlabstr = r'${\/+%s\/^{\circ}}$'%fmt
                else:
                    lonlabstr = r'${%s\/^{\circ}\/E}$'%fmt
            else:
                if labelstyle=='+/-':
                    lonlabstr = u'+%s\N{DEGREE SIGN}'%fmt
                else:
                    lonlabstr = u'%s\N{DEGREE SIGN}E'%fmt
            lonlab = lonlabstr%lon
        else:
            if rcParams['text.usetex']:
                lonlabstr = r'${%s\/^{\circ}}$'%fmt
            else:
                lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
            lonlab = lonlabstr%lon
    return lonlab

def _setlatlab(fmt,lat,labelstyle):
    # set lat label string (called by Basemap.drawparallels)
    try: # fmt is a function that returns a formatted string
           latlab = fmt(lat)
    except: # fmt is a format string.
        if lat<0:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    latlabstr = r'${\/-%s\/^{\circ}}$'%fmt
                else:
                    latlabstr = r'${%s\/^{\circ}\/S}$'%fmt
            else:
                if labelstyle=='+/-':
                    latlabstr = u'-%s\N{DEGREE SIGN}'%fmt
                else:
                    latlabstr = u'%s\N{DEGREE SIGN}S'%fmt
            latlab = latlabstr%np.fabs(lat)
        elif lat>0:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    latlabstr = r'${\/+%s\/^{\circ}}$'%fmt
                else:
                    latlabstr = r'${%s\/^{\circ}\/N}$'%fmt
            else:
                if labelstyle=='+/-':
                    latlabstr = u'+%s\N{DEGREE SIGN}'%fmt
                else:
                    latlabstr = u'%s\N{DEGREE SIGN}N'%fmt
            latlab = latlabstr%lat
        else:
            if rcParams['text.usetex']:
                latlabstr = r'${%s\/^{\circ}}$'%fmt
            else:
                latlabstr = u'%s\N{DEGREE SIGN}'%fmt
            latlab = latlabstr%lat
    return latlab
