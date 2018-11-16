module crs;
import std.string;
import std.algorithm;
import std.conv;
import std.array;
import std.traits;
import std.range;

import proj;

/** Keys for Projection info */
enum PK {
    a = "+a",   /// Semimajor radius of the ellipsoid axis
    alpha = "+alpha",   /// Angle of rotation pole
    aperture = "+aperture", /// Aperture
    axis = "+axis", /// Axis orientation
    azi = "+azi",   /// Azimuth
    b = "+b",   /// Semiminor radius of the ellipsoid axis
    czech = "+czeck",   /// Reverse sign of output coordinates
    ellps = "+ellps",   /// Ellipsoid
    gamma = "+gamma",   /// Azimuth of centerline clockwise from north of the bearing
    guam = "+guam", /// Use Guam elliptical formulas
    h = "+h",   /// Height above earth
    k_0 = "+k_0",   /// Scaling factor
    lat_0 = "+lat_0",   /// Latitude of origin
    lat_1 = "+lat_1",   /// First standard parallel
    lat_2 = "+lat_2",   /// Second standard parallel
    lat_3 = "+lat_3",   /// Latitude of 3rd control point
    lat_b = "+lat_b",   /// Angular distance from tangency point
    lat_ts = "+lat_ts", /// Latitude of true scale
    lon_0 = "+lon_0",   /// Longitude of projection center
    lon_1 = "+lon_1",   /// Longitude of first control point
    lon_2 = "+lon_2",   /// Longitude of second control point
    lon_3 = "+lon_3",   /// Longitude of third control point
    lon_wrap = "+lon_wrap", /// Center longitude to use for wrapping
    lonc = "+lonc", /// Longitude of rotational pole point
    lsat = "+lsat", /// Landsat satelite used
    M = "+M",   /// Larger than 0
    m = "+m",
    mode = "+mode", /// Can be either plane, di, dd, or hex
    n = "+n",
    ns = "+ns", /// Use non-skewed cartesian coordinates
    no_cut = "+no_cut", /// Do not cut at hemisphere limit
    no_rot = "+no_rot", /// Do not rotate axis
    no_off = "+no_off", /// Do not offset origin
    north_square = "+north_square", /// Position of north square
    o_alpha = "+o_alpha",   /// Angle to rotate projection
    o_lat_c = "+o_lat_c",   /// Latitude of rotation point
    o_lat_p = "+o_lat_p",   /// Latitude of new pole
    o_lon_c = "+o_lon_c",   /// Longitude of rotation point
    o_lon_p = "+o_lon_p",   /// Longitude of new pole
    o_proj = "+o_proj", /// Oblique projection
    orient = "+orient", /// Can be either isea or pole
    over = "+over", /// Allow longitude output outside -180 to 180 range
    path = "+path", /// Path of satelite
    pm = "+pm", /// Alternate prime meridian
    proj = "+proj", /// Projection name
    q = "+q",
    R = "+R",   /// Radius of sphere in meters
    resolution = "+resolution", /// Resolution
    south = "+south",   /// Flag for using UTM in southern hemisphere
    south_square = "+south_square", /// Position of south square
    sweep = "+sweep",   /// Sweep angle axis (x or y)
    theta = "+theta",
    units = "+units",
    vunits = "+vunits",
    W = "+W",
    x_0 = "+x_0",   /// False easting in meters
    y_0 = "+y_0",   /// False northing in meters
    zone = "+zone"  /// UTM zone
}

/** Class representing an error from the proj library */
class PJException : Exception {
    /** PJException contructor */
    this(string msg, string file = __FILE__, size_t line = __LINE__) {
        super(msg, file, line);
    }
}

private extern (C) double identity(double x) {return x;}

/** Coordinate struct that functions identically to PJ_COORD union but with D-ish
    syntax. It can be used anywhere that accepts a PJ_COORD, and it supports indexing and
    slicing. All union members from PJ_COORD are forwarded as well.
*/
struct Coord {
    PJ_COORD rep; /// PJ_COORD underlying the Coord

    /** Coord constructor that takes 2 doubles (x, y; lon, lat) */
    this(double a, double b) @nogc nothrow {
        this.rep = proj_coord(a, b, 0, 0);
    }
    /** Coord constructor that takes 3 doubles (x, y, z; lon, lat, z) */
    this(double a, double b, double c) @nogc nothrow {
        this.rep = proj_coord(a, b, c, 0);
    }
    /** Coord constructor that takes 4 doubles (x, y, z, t; lon, lat, z, t) */
    this(double a, double b, double c, double d) @nogc nothrow {
        this.rep = proj_coord(a, b, c, d);
    }

    /** Coord constructor that takes a PJ_COORD */
    this(PJ_COORD input) @nogc nothrow {
        this.rep = input;
    }

    /** Index the Coord struct */
    double opIndex(int i) @nogc nothrow {return rep.v[i];}
    /** Slice the Coord struct */
    double[] opIndex() @nogc nothrow {return rep.v[];}
    /** Size of Coord struct (4) */
    @property int opDollar() @nogc nothrow {return 4;}
    /** Slice the Coord struct */
    double[] opSlice(int start, int end) @nogc nothrow {return rep.v[start..end];}

    alias rep this;
}

/** Cartographic projection type that wraps a PJ* and handles construction and destruction.
    Base class for all named projections. See https://proj4.org/operations/projections/index.html
    for options for each projection
*/
class Projection {

    protected PJ* rep;

    /** Construct a Projection using an associated array where keys are members of the PK enum
        and values are strings. Second optional argument is PJ_CONTEXT* that is used for
        thread safety and can be ignored in single-threaded cases. Throws a PJException
        if the proj library sets an error.
    */
    this(string[PK] init, PJ_CONTEXT* ctx=PJ_DEFAULT_CTX) {
        if (PK.proj !in init) {
            throw new Exception("Failed to specify projection");
        }
        string proj = init[PK.proj];
        init.remove(PK.proj);
        string base = "+proj=" ~ proj;
        foreach (key, item; init) {
            base ~= " " ~ key ~ "=" ~ item;
        }

        this.rep = proj_create(ctx, base.toStringz);
        if (proj_errno(this.rep)) {
            throw new PJException(to!string(proj_errno_string(proj_errno(this.rep)).fromStringz));
        }
    }

    /** Construct a Projection using a string that is sent to proj_create. Second optional
        argument is PJ_CONTEXT* that is used for thread safety and can be ignored in
        single-threaded cases. Throws a PJException if the proj library sets an error.
    */
    this(string init, PJ_CONTEXT* ctx=PJ_DEFAULT_CTX) {
        this.rep = proj_create(ctx, init.toStringz);
        if (proj_errno(this.rep)) {
            throw new PJException(to!string(proj_errno_string(proj_errno(this.rep)).fromStringz));
        }
    }

    /** Transform a Coord or PJ_COORD from geodetic coordinates to projection coordinates.
        If the second optional argument is set to true, transforms from projection coordinates
        to geodetic coordinates. Returns a Coord. Expects longitude and latitude to be in radians.
    */
    Coord transform(PJ_COORD val, bool inverse=false) nothrow {
        return Coord(proj_trans(this.rep, inverse ? PJ_INV : PJ_FWD, val));
    }

    /** Transform a lon and a lat in degrees (or x and y) from geodetic coordinates to projection
        coordinates. If the third optional argument is set to true, transforms from
        projection coordinates to geodetic coordinates. Returns [x, y]
        (or [lon, lat] in degrees if inverse is true) */
    double[2] transform(double lon, double lat, bool inverse=false) nothrow {
        immutable PJ_COORD inter = inverse ? proj_coord(lon, lat, 0, 0)
                                           : proj_coord(proj_torad(lon), proj_torad(lat), 0, 0);
        auto ret = proj_trans(this.rep, inverse ? PJ_INV : PJ_FWD, inter);
        if (inverse) {
            return [proj_todeg(ret.v[0]), proj_todeg(ret.v[1])];
        }
        else {
            return [ret.v[0], ret.v[1]];
        }
    }

    /** Transform a [lon, lat] in degrees (or [x, y]) from geodetic coordinates to projection
        coordinates. If the second optional argument is set to true, transforms from
        projection coordinates to geodetic coordinates. Returns [x, y]
        (or [lon, lat] in degrees if inverse is true) */
    double[2] transform(double[2] lonlat, bool inverse=false) nothrow {
        return transform(lonlat[0], lonlat[1], inverse);
    }

    /** Transform a range Coord or PJ_COORD from geodetic coordinates to projection coordinates.
        If the second optional argument is set to true, transforms from projection coordinates
        to geodetic coordinates. Returns a range of Coord.
        Expects longitude and latitude to be in radians. Is probably slower than proj_trans_array.
    */
    auto transform(T)(T range, bool inverse=false) nothrow
        if (isInputRange!T && is(ElementType!T == PJ_COORD) || is(ElementType!T == Coord)) {
        immutable direction = inverse ? PJ_INV : PJ_FWD;
        return range.map!(a => Coord(proj_trans(a.rep, direction, a.rep)));
    }

    /** Transform a range of [lon, lat] in degrees (or [x, y]) from geodetic coordinates
        to projection coordinates. If the second optional argument is set to true, transforms
        from projection coordinates to geodetic coordinates. Returns a range of [x, y]
        (or [lon, lat] in degrees if inverse is true) */
    auto transform(T)(T range, bool inverse=false) nothrow
        if (isInputRange!T && isRandomAccessRange!(ElementType!T)
            && __traits(isFloating, ElementType!(ElementType!T))) {
        immutable conv1 = inverse ? &identity : &proj_torad;
        immutable conv2 = inverse ? &proj_todeg : &identity;
        immutable direction = inverse ? PJ_INV : PJ_FWD;
        return range.map!(x => proj_coord(conv1(x[0]), conv1(x[1]), 0, 0))
                    .map!(x => proj_trans(this.rep, direction, x))
                    .map!(x => [conv2(x.v[0]), conv2(x.v[1])]);
    }

    /** Transform a range of lon and a range of lat in degrees (or x and y) from geodetic
        coordinates to projection coordinates. If the third optional argument is set to true,
        transforms from projection coordinates to geodetic coordinates. Returns a range of [x, y]
        (or [lon, lat] in degrees if inverse is true) */
    auto transform(T, R)(T lons, R lats, bool inverse=false) nothrow
        if (isInputRange!T && __traits(isFloating, ElementType!T)
            && isInputRange!R && __traits(isFloating, ElementType!R)) {
        immutable conv1 = inverse ? &identity : &proj_torad;
        immutable conv2 = inverse ? &proj_todeg : &identity;
        immutable direction = inverse ? PJ_INV : PJ_FWD;
        return zip(lons, lats)
                  .map!(x => proj_coord(conv1(x[0]), conv1(x[1]), 0, 0))
                  .map!(x => proj_trans(this.rep, direction, x))
                  .map!(x => [conv2(x.v[0]), conv2(x.v[1])]);
    }

    ~this() {
        proj_destroy(this.rep);
    }

    /** Calculate geodesic distance between two points in geodetic coordinates (lon, lat).
        Returns the distance in meters
    */
    double lp_dist(PJ_COORD a, PJ_COORD b) nothrow {
        return proj_lp_dist(this.rep, a, b);
    }

    /** Calculate geodesic distance between two points in geodetic coordinates (lon, lat, z).
        Returns the distance in meters
    */
    double lpz_dist(PJ_COORD a, PJ_COORD b) nothrow {
        return proj_lpz_dist(this.rep, a, b);
    }
}

private enum projs = [
    "AlbersEqualArea":"aea",
    "AzimuthalEquidistant":"aeqd",
    "Airy":"airy",
    "Aitoff":"aitoff",
    "ModifiedSteroegraphicsOfAlaska":"alsk",
    "ApianGlobular1":"apian",
    "AugustEpicycloidal":"august",
    "BaconGlobular":"bacon",
    "BipolarConicOfWesternHemisphere":"bipc",
    "BoggsEumorphic":"boggs",
    "CalCoop":"calcofi",
    "Cassini":"cass",
    "CentralConic":"cc",
    "EqualAreaCylindrical":"cea",
    "Collignon":"collg",
    "CrasterParabolic":"crast",
    "DenoyerSemiEliptical":"denoy",
    "Eckert1":"eck1",
    "Eckert2":"eck2",
    "Eckert3":"eck3",
    "Eckert4":"eck4",
    "Eckert5":"eck5",
    "Eckert6":"eck6",
    "EquidistantCylindrical":"eqc",
    "PlateCarree":"eqc",
    "EqualEarth":"eqearth",
    "ExtendedTransverseMercator":"etmerc",
    "Fahey":"fehey",
    "Foucaut":"fouc",
    "FoucautSinusoidal":"fouc_s",
    "Gall":"gall",
    "Ginsburg":"gins8",
    "GeneralSinusoidalSeries":"gn_sinu",
    "Gnomonic":"gnom",
    "GoodeHomolosine":"goode",
    "ModStereographics48US":"gs48",
    "ModStereographics50US":"gs50",
    "HammerEckertGreifendorff":"hammer",
    "HatanoAsymmetricalEqualArea":"hatano",
    "HEALPix":"healpix",
    "rHEALPix":"rhealpix",
    "InterruptedGoodeHomolosine":"igh",
    "IcosahedralSnyderEqualArea":"isea",
    "Kavraisky5":"kav5",
    "Kavraisky7":"kav7",
    "Krovak":"krovak",
    "LambertAzimuthalEqualArea":"laea",
    "Lagrange":"lagrng",
    "Larrivee":"larr",
    "Laskowski":"lask",
    "LambertConformalConic":"lcc",
    "LambertConformalConicAlt":"lcca",
    "LambertEqualAreaConic":"leac",
    "LeeOblatedStereographic":"lee_os",
    "Loximuthal":"loxim",
    "McBrydeThomasFlatPolarSine1":"mbt_s",
    "McBrydeThomasFlatPoleSine2":"mbt_fps",
    "McBrideThomasFlatPolarParabolic":"mbtfpp",
    "McBrydeThomasFlatPolarQuartic":"mbtfpq",
    "McBrydeThomasFlatPolarSinusoidal":"mbtfps",
    "Mercator":"merc",
    "MillerOblatedStereographic":"mil_os",
    "MillerCylindrical":"mill",
    "Mollweide":"moll",
    "NaturalEarth":"natearth",
    "Nell":"nell",
    "NellHammer":"nell_h",
    "NicolosiGlobular":"nicol",
    "NewZealandMapGrid":"nzmg",
    "UniversalTransverseMercator":"utm",
    "Bonne":"bonne",
    "CentralConic":"ccon",
    "ChamberlinTrimetric":"chamb",
    "EquidistantConic":"eqdc",
    "Euler":"euler",
    "Geostationary":"geos",
    "IMWPolyconic":"imw_p",
    "Lamborde":"labrd",
    "Landsat":"lsat",
    "Murdoch1":"murd1",
    "Murdoch2":"murd2",
    "Murdoch3":"murd3",
    "NearsidedPerspective":"nsper",
    "GeneralObliqueTransformation":"ob_tran",
    "ObliqueCylindricalEqualArea":"ocea",
    "OblatedEqualArea":"oea",
    "ObliqueMercator":"omerc",
    "OrteliusOval":"ortel",
    "Orthographic":"ortho",
    "PerspectiveConic":"pconic",
    "Polyconic":"poly",
    "PutninsP1":"putp1",
    "PutninsP2":"putp2",
    "PutninsP3":"putp3",
    "PutninsP3Prime":"putp3p",
    "PutninsP4Prime":"putp4p",
    "PutninsP5":"putp5",
    "PutninsP5Prime":"putp5p",
    "PutninsP6":"putp6",
    "PutninsP6Prime":"putp6p",
    "QuarticAuthalic":"qua_aut",
    "QuadrilateralizedSphericalCube":"qsc",
    "Robinson":"robin",
    "RoussilheStereographic":"rouss",
    "RectangularPolyconic":"rpoly",
    "Sinusoidal":"sinu",
    "SwissObliqueMercator":"somerc",
    "Stereographic":"stere",
    "ObliqueStereographicAlternative":"sterea",
    "GaussSchreiberTransverseMercator":"gstmerc",
    "TransverseCentralCylindrical":"tcc",
    "TransverseCylindricalEqualArea":"tcea",
    "Tissot":"tissot",
    "TransverseMercator":"tmerc",
    "TwoPointEquidistant":"tpeqd",
    "TiltedPerspective":"tpers",
    "UniversialPolarStereographic":"ups",
    "Urmaev":"urm5",
    "UrmaevFlatPolarSinusoidal":"urmfps",
    "VanDerGrinten1":"vandg",
    "VanDerGrinten2":"vandg2",
    "VanDerGrinten3":"vandg3",
    "VanDerGrinten4":"vandg4",
    "Vitkovsky":"vitk1",
    "Wagner1":"wag1",
    "Wagner2":"wag2",
    "Wagner3":"wag3",
    "Wagner4":"wag4",
    "Wagner5":"wag5",
    "Wagner6":"wag6",
    "Wagner7":"wag7",
    "WebMercator":"webmerc",
    "Werenskoild":"weren",
    "Winkel1":"wink1",
    "Winkel2":"wink2",
    "WinkelTripel":"wintri"
];

static foreach(key, item; projs) {
    mixin(
    "final class " ~ key ~ " : Projection {
        this(string[PK] init, PJ_CONTEXT* ctx=PJ_DEFAULT_CTX) {
            init[PK.proj] = \"" ~ item ~ "\";
            super(init);
        }

        this (PJ_CONTEXT* ctx=PJ_DEFAULT_CTX) {
            super([PK.proj:\"" ~ item ~ "\"]);
        }
    }\n");
}
