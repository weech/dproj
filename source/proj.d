module proj;
/** Proj bindings for D */

/** The version numbers should be updated with every release! **/
enum PROJ_VERSION_MAJOR = 6;
enum PROJ_VERSION_MINOR = 0;
enum PROJ_VERSION_PATCH = 0;

enum PJ_LOG_LEVEL {
    PJ_LOG_NONE  = 0,
    PJ_LOG_ERROR = 1,
    PJ_LOG_DEBUG = 2,
    PJ_LOG_TRACE = 3,
    PJ_LOG_TELL  = 4,
    PJ_LOG_DEBUG_MAJOR = 2, /* for proj_api.h compatibility */
    PJ_LOG_DEBUG_MINOR = 3  /* for proj_api.h compatibility */
}

enum PJ_DEFAULT_CTX = null; /// Is 0 in C header but needs to be null in D

/** Apply transformation to observation - in forward or inverse direction */
alias PJ_DIRECTION = int;
enum PJ_FWD   =  1;   /** Forward    */
enum PJ_IDENT =  0;   /** Do nothing */
enum PJ_INV   = -1;    /** Inverse    */

extern (C)
{

extern __gshared const byte[] pj_release; /** global release id string */

/* first forward declare everything needed */

struct PJ_AREA;

/** Common designation */
struct P5_FACTORS {
    double meridional_scale;               /** h */
    double parallel_scale;                 /** k */
    double areal_scale;                    /** s */

    double angular_distortion;             /** omega */
    double meridian_parallel_angle;        /** theta-prime */
    double meridian_convergence;           /** alpha */

    double tissot_semimajor;               /** a */
    double tissot_semiminor;               /** b */

    double dx_dlam, dx_dphi;
    double dy_dlam, dy_dphi;
}
alias PJ_FACTORS = P5_FACTORS;

/** Data type for projection/transformation information */
struct PJconsts;
alias PJ =  PJconsts;         /** the PJ object herself */

/** Data type for library level information */
struct PJ_INFO
{
    int         major;              /** Major release number                 */
    int         minor;              /** Minor release number                 */
    int         patch;              /** Patch level                          */
    const char* release;           /** Release info. Version + date         */
    const char* _version;           /** Full version number                  */
    const char* searchpath;        /** Paths where init and grid files are  */
                                    /* looked for. Paths are separated by   */
                                    /* semi-colons.                         */
    const char** paths;
    size_t path_count;
}

struct PJ_PROJ_INFO {
    const char*  id;                /** Name of the projection in question                       */
    const char*  description;       /** Description of the projection                            */
    const char*  definition;        /** Projection definition                                    */
    int         has_inverse;        /** 1 if an inverse mapping exists, 0 otherwise              */
    double      accuracy;           /** Expected accuracy of the transformation. -1 if unknown.  */
}

struct PJ_GRID_INFO {
    byte[32]        gridname;       /** name of grid                         */
    byte[260]       filename;      /** full path to grid                    */
    byte[8]        format;          /** file format of grid                  */
    PJ_LP       lowerleft;          /** Coordinates of lower left corner     */
    PJ_LP       upperright;         /** Coordinates of upper right corner    */
    int         n_lon, n_lat;       /** Grid size                            */
    double      cs_lon, cs_lat;     /** Cell size of grid                    */
}

struct PJ_INIT_INFO {
    byte[32]        name;           /* name of init file                        */
    byte[260]       filename;      /* full path to the init file.              */
    byte[32]        _version;        /* version of the init file                 */
    byte[32]        origin;         /* origin of the file, e.g. EPSG            */
    byte[16]        lastupdate;     /* Date of last update in YYYY-MM-DD format */
}

/** Data types for list of operations, ellipsoids, datums and units used in PROJ.4 */
struct PJ_LIST {
    const char*  id;                /** projection keyword */
    PJ** proj;     /** projection entry point */
    const char* descr;     /** description text */
}

alias PJ_OPERATIONS = PJ_LIST;

struct PJ_ELLPS {
    const char*  id;    /** ellipse keyword name */
    const char*  major; /** a= value */
    const char*  ell;   /** elliptical parameter */
    const char*  name;  /** comments */
}

struct PJ_UNITS {
    const char*  id;        /** units keyword */
    const char*  to_meter;  /** multiply by value to get meters */
    const char*  name;      /** comments */
    double factor;     /** to_meter factor in actual numbers */
}

struct PJ_PRIME_MERIDIANS {
    const char* id;        /** prime meridian keyword */
    const char* defn;      /** offset from greenwich in DMS format. */
}

/* Geodetic, mostly spatiotemporal coordinate types */
struct PJ_XYZT {
    double x, y, z, t;
}
struct PJ_UVWT {
    double u, v, w, t;
}
struct PJ_LPZT {
    double lam, phi, z, t;
}
/** Rotations: omega, phi, kappa */
struct PJ_OPK {
    double o, p, k;
}
/** East, North, Up */
struct PJ_ENU {
    double e, n, u;
}
/** Geodesic length, fwd azi, rev azi */
struct PJ_GEOD {
    double s, a1, a2;
}

/** Classic proj.4 pair/triplet types - moved into the PJ_ name space */
struct PJ_UV { double   u,   v; }
struct PJ_XY { double   x,   y; }
struct PJ_LP { double lam, phi; }

struct PJ_XYZ { double   x,   y,  z; }
struct PJ_UVW { double   u,   v,  w; }
struct PJ_LPZ { double lam, phi,  z; }

/** Avoid preprocessor renaming and implicit type-punning: Use a union to make it explicit.
    Data type for generic geodetic 3D data plus epoch information */
union PJ_COORD {
    double[4] v;   /** First and foremost, it really is "just 4 numbers in a vector" */
    PJ_XYZT xyzt;
    PJ_UVWT uvwt;
    PJ_LPZT lpzt;
    PJ_GEOD geod;
    PJ_OPK opk;
    PJ_ENU enu;
    PJ_XYZ xyz;
    PJ_UVW uvw;
    PJ_LPZ lpz;
    PJ_XY xy;
    PJ_UV uv;
    PJ_LP lp;
}


/** The context type - properly namespaced synonym for projCtx */
struct projCtx_t;
alias PJ_CONTEXT = projCtx_t;

/* A P I */

/* Functionality for handling thread contexts */
PJ_CONTEXT* proj_context_create() @nogc nothrow ;
PJ_CONTEXT* proj_context_destroy(PJ_CONTEXT* ctx) @nogc nothrow ;

/* Manage the transformation definition object PJ */
PJ* proj_create(PJ_CONTEXT* ctx, const char* definition) @nogc nothrow ;
PJ* proj_create_argv(PJ_CONTEXT* ctx, int argc, char** argv)@nogc nothrow ;
PJ* proj_create_crs_to_crs(PJ_CONTEXT* ctx, const char* srid_from, const char* srid_to, PJ_AREA* area) @nogc nothrow ;
PJ* proj_destroy(PJ* P) @nogc nothrow ;

int proj_angular_input(PJ* P, PJ_DIRECTION dir) @nogc nothrow ;
int proj_angular_output(PJ* P, PJ_DIRECTION dir) @nogc nothrow ;

PJ_COORD proj_trans(PJ* P, PJ_DIRECTION direction, PJ_COORD coord) @nogc nothrow ;
int proj_trans_array(PJ* P, PJ_DIRECTION direction, size_t n, PJ_COORD* coord)@nogc nothrow ;

size_t proj_trans_generic (
    PJ* P,
    PJ_DIRECTION direction,
    double* x, size_t sx, size_t nx,
    double* y, size_t sy, size_t ny,
    double* z, size_t sz, size_t nz,
    double* t, size_t st, size_t nt
) @nogc nothrow ;

/** Initializers */
PJ_COORD proj_coord(double x, double y, double z, double t) @nogc nothrow ;

/** Measure internal consistency - in forward or inverse direction */
double proj_roundtrip(PJ* P, PJ_DIRECTION direction, int n, PJ_COORD* coord) @nogc nothrow ;

/** Geodesic distance between two points with angular 2D coordinates */
double proj_lp_dist(const PJ* P, PJ_COORD a, PJ_COORD b) @nogc nothrow ;

/** The geodesic distance AND the vertical offset */
double proj_lpz_dist(const PJ* P, PJ_COORD a, PJ_COORD b) @nogc nothrow ;

/** Euclidean distance between two points with linear 2D coordinates */
double proj_xy_dist(PJ_COORD a, PJ_COORD b) @nogc nothrow ;

/** Euclidean distance between two points with linear 3D coordinates */
double proj_xyz_dist(PJ_COORD a, PJ_COORD b) @nogc nothrow ;

/** Geodesic distance (in meter) + fwd and rev azimuth between two points on the ellipsoid */
PJ_COORD proj_geod(const PJ* P, PJ_COORD a, PJ_COORD b) @nogc nothrow ;

/** Set or read error level */
int  proj_context_errno(PJ_CONTEXT* ctx) @nogc nothrow ;
int  proj_errno(const PJ* P) @nogc nothrow ;
int  proj_errno_set(const PJ* P, int err) @nogc nothrow ;
int  proj_errno_reset(const PJ* P) @nogc nothrow ;
int  proj_errno_restore(const PJ* P, int err) @nogc nothrow ;
const(char)* proj_errno_string(int err) @nogc nothrow ;

alias PJ_LOG_FUNCTION = void* function(void*, int, const char*) @nogc nothrow ;

PJ_LOG_LEVEL proj_log_level(PJ_CONTEXT* ctx, PJ_LOG_LEVEL log_level) @nogc nothrow ;
void proj_log_func(PJ_CONTEXT* ctx, void* app_data, PJ_LOG_FUNCTION logf) @nogc nothrow ;

/** Scaling and angular distortion factors */
PJ_FACTORS proj_factors(PJ* P, PJ_COORD lp) @nogc nothrow ;

/** Info functions - get information about various PROJ.4 entities */
PJ_INFO proj_info() @nogc nothrow ;
PJ_PROJ_INFO proj_pj_info(PJ* P) @nogc nothrow ;
PJ_GRID_INFO proj_grid_info(const char* gridname) @nogc nothrow ;
PJ_INIT_INFO proj_init_info(const char* initname) @nogc nothrow ;

/* List functions: */
/** Get lists of operations, ellipsoids, units and prime meridians. */
const(PJ_OPERATIONS)* proj_list_operations() @nogc nothrow ;
const(PJ_ELLPS)* proj_list_ellps() @nogc nothrow ;
const(PJ_UNITS)* proj_list_units() @nogc nothrow ;
const(PJ_PRIME_MERIDIANS)* proj_list_prime_meridians() @nogc nothrow ;

/* These are trivial, and while occasionally useful in real code, primarily here to      */
/* simplify demo code, and in acknowledgement of the proj-internal discrepancy between   */
/* angular units expected by classical proj, and by Charles Karney's geodesics subsystem */
double proj_torad(double angle_in_degrees) @nogc nothrow ;
double proj_todeg(double angle_in_radians) @nogc nothrow ;

/** Geographical to geocentric latitude - another of the "simple, but useful" */
PJ_COORD proj_geocentric_latitude(const PJ* P, PJ_DIRECTION direction, PJ_COORD coord) @nogc nothrow ;

double proj_dmstor(const char* _is, char** rs) @nogc nothrow ;
char*  proj_rtodms(char* s, double r, int pos, int neg) @nogc nothrow ;

}
