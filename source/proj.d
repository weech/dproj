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
    const(char)* release;           /** Release info. Version + date         */
    const(char)* _version;           /** Full version number                  */
    const(char)* searchpath;        /** Paths where init and grid files are  */
                                    /* looked for. Paths are separated by   */
                                    /* semi-colons.                         */
    const(char)** paths;
    size_t path_count;
}

struct PJ_PROJ_INFO {
    const(char)*  id;                /** Name of the projection in question                       */
    const(char)*  description;       /** Description of the projection                            */
    const(char)*  definition;        /** Projection definition                                    */
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
    const(char)*  id;                /** projection keyword */
    PJ** proj;     /** projection entry point */
    const(char)* descr;     /** description text */
}

alias PJ_OPERATIONS = PJ_LIST;

struct PJ_ELLPS {
    const(char)*  id;    /** ellipse keyword name */
    const(char)*  major; /** a= value */
    const(char)*  ell;   /** elliptical parameter */
    const(char)*  name;  /** comments */
}

struct PJ_UNITS {
    const(char)*  id;        /** units keyword */
    const(char)*  to_meter;  /** multiply by value to get meters */
    const(char)*  name;      /** comments */
    double factor;     /** to_meter factor in actual numbers */
}

struct PJ_PRIME_MERIDIANS {
    const(char)* id;        /** prime meridian keyword */
    const(char)* defn;      /** offset from greenwich in DMS format. */
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

/** Callback to resolve a filename to a full path */
alias proj_file_finder = char* function(PJ_CONTEXT*, const(char)*, void*) @nogc nothrow ;

void proj_context_set_file_finder(PJ_CONTEXT* ctx, proj_file_finder finder, void* user_data);
void proj_context_set_search_paths(PJ_CONTEXT* cts, int count_paths, char** paths);

void proj_context_use_proj4_init_rules(PJ_CONTEXT* ctx, int enable) @nogc nothrow ;
int proj_context_get_use_proj4_init_rules(PJ_CONTEXT* ctx, int from_legacy_code_path) @nogc nothrow ;

/* Manage the transformation definition object PJ */
PJ* proj_create(PJ_CONTEXT* ctx, const(char)* definition) @nogc nothrow ;
PJ* proj_create_argv(PJ_CONTEXT* ctx, int argc, char** argv) @nogc nothrow ;
PJ* proj_create_crs_to_crs(PJ_CONTEXT* ctx, const(char)* srid_from, const(char)* srid_to, PJ_AREA* area) @nogc nothrow ;
PJ* proj_normalize_for_visualization(PJ_CONTEXT* ctx, const(char)* source_crs,
                                     const(char)* target_crs, PJ_AREA* area) @nogc nothrow ;
PJ* proj_destroy(PJ* P) @nogc nothrow ;

PJ_AREA* proj_area_create() @nogc nothrow ;
void proj_area_set_bbox(PJ_AREA* area, double west_lon_degree, double south_lat_degree,
                        double east_lon_degree, double north_lat_degree) @nogc nothrow ;
void proj_area_destroy(PJ_AREA* area);

int proj_angular_input(PJ* P, PJ_DIRECTION dir) @nogc nothrow ;
int proj_angular_output(PJ* P, PJ_DIRECTION dir) @nogc nothrow ;

PJ_COORD proj_trans(PJ* P, PJ_DIRECTION direction, PJ_COORD coord) @nogc nothrow ;
int proj_trans_array(PJ* P, PJ_DIRECTION direction, size_t n, PJ_COORD* coord) @nogc nothrow ;

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

alias PJ_LOG_FUNCTION = void* function(void*, int, const(char)*) @nogc nothrow ;

PJ_LOG_LEVEL proj_log_level(PJ_CONTEXT* ctx, PJ_LOG_LEVEL log_level) @nogc nothrow ;
void proj_log_func(PJ_CONTEXT* ctx, void* app_data, PJ_LOG_FUNCTION logf) @nogc nothrow ;

/** Scaling and angular distortion factors */
PJ_FACTORS proj_factors(PJ* P, PJ_COORD lp) @nogc nothrow ;

/** Info functions - get information about various PROJ.4 entities */
PJ_INFO proj_info() @nogc nothrow ;
PJ_PROJ_INFO proj_pj_info(PJ* P) @nogc nothrow ;
PJ_GRID_INFO proj_grid_info(const(char)* gridname) @nogc nothrow ;
PJ_INIT_INFO proj_init_info(const(char)* initname) @nogc nothrow ;

/* List functions: */
/** Get lists of operations, ellipsoids, units and prime meridians. */
const(PJ_OPERATIONS)* proj_list_operations() @nogc nothrow ;
const(PJ_ELLPS)* proj_list_ellps() @nogc nothrow ;
const(PJ_UNITS)* proj_list_units() @nogc nothrow ;
const(PJ_UNITS)* proj_list_angular_units() @nogc nothrow ;
const(PJ_PRIME_MERIDIANS)* proj_list_prime_meridians() @nogc nothrow ;

/* These are trivial, and while occasionally useful in real code, primarily here to      */
/* simplify demo code, and in acknowledgement of the proj-internal discrepancy between   */
/* angular units expected by classical proj, and by Charles Karney's geodesics subsystem */
double proj_torad(double angle_in_degrees) @nogc nothrow ;
double proj_todeg(double angle_in_radians) @nogc nothrow ;

/** Geographical to geocentric latitude - another of the "simple, but useful" */
PJ_COORD proj_geocentric_latitude(const PJ* P, PJ_DIRECTION direction, PJ_COORD coord) @nogc nothrow ;

double proj_dmstor(const(char)* _is, char** rs) @nogc nothrow ;
char*  proj_rtodms(char* s, double r, int pos, int neg) @nogc nothrow ;

} // End extern (C)

// End of normal C API. Beginning of C++ API bindings (still found in proj.h)

/** Type representing a NULL terminated list of NULL-terminate strings. */
alias PROJ_STRING_LIST = char**;

/** Guessed WKT "dialect". */
enum PJ_GUESSED_WKT_DIALECT {
    /** WKT2_2018 */
    PJ_GUESSED_WKT2_2018,

    /** WKT2_2015 */
    PJ_GUESSED_WKT2_2015,

    /** WKT1 */
    PJ_GUESSED_WKT1_GDAL,

    /** ESRI variant of WKT1 */
    PJ_GUESSED_WKT1_ESRI,

    /** Not WKT / unrecognized */
    PJ_GUESSED_NOT_WKT
}

/** Object category. */
enum PJ_CATEGORY {
    PJ_CATEGORY_ELLIPSOID,
    PJ_CATEGORY_PRIME_MERIDIAN,
    PJ_CATEGORY_DATUM,
    PJ_CATEGORY_CRS,
    PJ_CATEGORY_COORDINATE_OPERATION
}

/**  Object type. */
enum PJ_TYPE
{
    PJ_TYPE_UNKNOWN,

    PJ_TYPE_ELLIPSOID,

    PJ_TYPE_PRIME_MERIDIAN,

    PJ_TYPE_GEODETIC_REFERENCE_FRAME,
    PJ_TYPE_DYNAMIC_GEODETIC_REFERENCE_FRAME,
    PJ_TYPE_VERTICAL_REFERENCE_FRAME,
    PJ_TYPE_DYNAMIC_VERTICAL_REFERENCE_FRAME,
    PJ_TYPE_DATUM_ENSEMBLE,

    /** Abstract type, not returned by proj_get_type() */
    PJ_TYPE_CRS,

    PJ_TYPE_GEODETIC_CRS,
    PJ_TYPE_GEOCENTRIC_CRS,

    /** proj_get_type() will never return that type, but
      PJ_TYPE_GEOGRAPHIC_2D_CRS or PJ_TYPE_GEOGRAPHIC_3D_CRS. */
    PJ_TYPE_GEOGRAPHIC_CRS,

    PJ_TYPE_GEOGRAPHIC_2D_CRS,
    PJ_TYPE_GEOGRAPHIC_3D_CRS,
    PJ_TYPE_VERTICAL_CRS,
    PJ_TYPE_PROJECTED_CRS,
    PJ_TYPE_COMPOUND_CRS,
    PJ_TYPE_TEMPORAL_CRS,
    PJ_TYPE_ENGINEERING_CRS,
    PJ_TYPE_BOUND_CRS,
    PJ_TYPE_OTHER_CRS,

    PJ_TYPE_CONVERSION,
    PJ_TYPE_TRANSFORMATION,
    PJ_TYPE_CONCATENATED_OPERATION,
    PJ_TYPE_OTHER_COORDINATE_OPERATION,
}

/** Comparison criterion. */
enum PJ_COMPARISON_CRITERION
{
    /** All properties are identical. */
    PJ_COMP_STRICT,

    /** The objects are equivalent for the purpose of coordinate
     operations. They can differ by the name of their objects,
     identifiers, other metadata.
     Parameters may be expressed in different units, provided that the
     value is (with some tolerance) the same once expressed in a
     common unit.
    */
    PJ_COMP_EQUIVALENT,

    /** Same as EQUIVALENT, relaxed with an exception that the axis order
     of the base CRS of a DerivedCRS/ProjectedCRS or the axis order of
     a GeographicCRS is ignored. Only to be used
     with DerivedCRS/ProjectedCRS/GeographicCRS */
    PJ_COMP_EQUIVALENT_EXCEPT_AXIS_ORDER_GEOGCRS,
}

/** WKT version. */
enum PJ_WKT_TYPE
{
    /** cf osgeo::proj::io::WKTFormatter::Convention::WKT2 */
    PJ_WKT2_2015,
    /** cf osgeo::proj::io::WKTFormatter::Convention::WKT2_SIMPLIFIED */
    PJ_WKT2_2015_SIMPLIFIED,
    /** cf osgeo::proj::io::WKTFormatter::Convention::WKT2_2018 */
    PJ_WKT2_2018,
    /** cf osgeo::proj::io::WKTFormatter::Convention::WKT2_2018_SIMPLIFIED */
    PJ_WKT2_2018_SIMPLIFIED,
    /** cf osgeo::proj::io::WKTFormatter::Convention::WKT1_GDAL */
    PJ_WKT1_GDAL,
    /** cf osgeo::proj::io::WKTFormatter::Convention::WKT1_ESRI */
    PJ_WKT1_ESRI
}

/** Specify how source and target CRS extent should be used to restrict
  * candidate operations (only taken into account if no explicit area of
  * interest is specified.) */
enum PROJ_CRS_EXTENT_USE
{
    /** Ignore CRS extent */
    PJ_CRS_EXTENT_NONE,

    /** Test coordinate operation extent against both CRS extent. */
    PJ_CRS_EXTENT_BOTH,

    /** Test coordinate operation extent against the intersection of both
        CRS extent. */
    PJ_CRS_EXTENT_INTERSECTION,

    /** Test coordinate operation against the smallest of both CRS extent. */
    PJ_CRS_EXTENT_SMALLEST
}

/** Describe how grid availability is used. */
enum PROJ_GRID_AVAILABILITY_USE {
    /** Grid availability is only used for sorting results. Operations
         where some grids are missing will be sorted last. */
    PROJ_GRID_AVAILABILITY_USED_FOR_SORTING,

    /** Completely discard an operation if a required grid is missing. */
    PROJ_GRID_AVAILABILITY_DISCARD_OPERATION_IF_MISSING_GRID,

    /** Ignore grid availability at all. Results will be presented as if
         all grids were available. */
    PROJ_GRID_AVAILABILITY_IGNORED,
}

/** PROJ string version. */
enum PJ_PROJ_STRING_TYPE
{
    /** cf osgeo::proj::io::PROJStringFormatter::Convention::PROJ_5 */
    PJ_PROJ_5,
    /** cf osgeo::proj::io::PROJStringFormatter::Convention::PROJ_4 */
    PJ_PROJ_4
}

/** Spatial criterion to restrict candidate operations. */
enum PROJ_SPATIAL_CRITERION {
    /** The area of validity of transforms should strictly contain the
         are of interest. */
    PROJ_SPATIAL_CRITERION_STRICT_CONTAINMENT,

    /** The area of validity of transforms should at least intersect the
         area of interest. */
    PROJ_SPATIAL_CRITERION_PARTIAL_INTERSECTION
}

/** Describe if and how intermediate CRS should be used */
enum PROJ_INTERMEDIATE_CRS_USE {
    /** Always search for intermediate CRS. */
    PROJ_INTERMEDIATE_CRS_USE_ALWAYS,

    /** Only attempt looking for intermediate CRS if there is no direct
         transformation available. */
    PROJ_INTERMEDIATE_CRS_USE_IF_NO_DIRECT_TRANSFORMATION,

    /* Do not attempt looking for intermediate CRS. */
    PROJ_INTERMEDIATE_CRS_USE_NEVER,
}

/** Type of coordinate system. */
enum PJ_COORDINATE_SYSTEM_TYPE
{
    PJ_CS_TYPE_UNKNOWN,

    PJ_CS_TYPE_CARTESIAN,
    PJ_CS_TYPE_ELLIPSOIDAL,
    PJ_CS_TYPE_VERTICAL,
    PJ_CS_TYPE_SPHERICAL,
    PJ_CS_TYPE_ORDINAL,
    PJ_CS_TYPE_PARAMETRIC,
    PJ_CS_TYPE_DATETIMETEMPORAL,
    PJ_CS_TYPE_TEMPORALCOUNT,
    PJ_CS_TYPE_TEMPORALMEASURE
}

/** Structure given overall description of a CRS.

  This structure may grow over time, and should not be directly allocated by
  client code.
*/
struct PROJ_CRS_INFO
{
    /** Authority name. */
    char* auth_name;
    /** Object code. */
    char* code;
    /** Object name. */
    char* name;
    /** Object type. */
    PJ_TYPE type;
    /** Whether the object is deprecated */
    int isdeprecated;
    /** Whereas the west_lon_degree, south_lat_degree, east_lon_degree and
      north_lat_degree fields are valid. */
    int bbox_valid;
    /** Western-most longitude of the area of use, in degrees. */
    double west_lon_degree;
    /** Southern-most latitude of the area of use, in degrees. */
    double south_lat_degree;
    /** Eastern-most longitude of the area of use, in degrees. */
    double east_lon_degree;
    /** Northern-most latitude of the area of use, in degrees. */
    double north_lat_degree;
    /** Name of the area of use. */
    char* area_name;
    /** Name of the projection method for a projected CRS. Might be NULL even
     for projected CRS in some cases. */
    char* projection_method_name;
}

/** Structure describing optional parameters for proj_get_crs_list();

  This structure may grow over time, and should not be directly allocated by
  client code.
*/
struct PROJ_CRS_LIST_PARAMETERS
{
    /** Array of allowed object types. Should be NULL if all types are allowed*/
    const PJ_TYPE* types;
    /** Size of types. Should be 0 if all types are allowed*/
    size_t typesCount;

    /** If TRUE and bbox_valid == TRUE, then only CRS whose area of use
     * entirely contains the specified bounding box will be returned.
     * If FALSE and bbox_valid == TRUE, then only CRS whose area of use
     * intersects the specified bounding box will be returned.
     */
    int crs_area_of_use_contains_bbox;
    /** To set to TRUE so that west_lon_degree, south_lat_degree,
     * east_lon_degree and north_lat_degree fields are taken into account. */
    int bbox_valid;
    /** Western-most longitude of the area of use, in degrees. */
    double west_lon_degree;
    /** Southern-most latitude of the area of use, in degrees. */
    double south_lat_degree;
    /** Eastern-most longitude of the area of use, in degrees. */
    double east_lon_degree;
    /** Northern-most latitude of the area of use, in degrees. */
    double north_lat_degree;

    /** Whether deprecated objects are allowed. Default to FALSE. */
    int allow_deprecated;
}

struct PJ_OBJ_LIST;

extern (C) {

void proj_string_list_destroy(PROJ_STRING_LIST list) @nogc nothrow ;

int proj_context_set_database_path(PJ_CONTEXT* ctx, const(char)* dbPath,
                                   char** auxDbPaths, char** options) @nogc nothrow ;

const(char)* proj_context_get_database_path(PJ_CONTEXT* ctx) @nogc nothrow ;

const(char)* proj_context_get_database_metadata(PJ_CONTEXT* ctx, const(char)* key) @nogc nothrow ;

PJ_GUESSED_WKT_DIALECT proj_context_guess_wkt_dialect(PJ_CONTEXT* ctx, const(char)* wkt) @nogc nothrow ;

PJ* proj_create_from_wkt(PJ_CONTEXT* ctx, const(char)* wkt, char** options,
                         PROJ_STRING_LIST* out_warnings, PROJ_STRING_LIST* out_grammar_errors) @nogc nothrow ;

PJ* proj_create_from_database(PJ_CONTEXT* ctx, const(char)* auth_name, const(char)* code,
                              PJ_CATEGORY category, int usePRJAlternativeGridNames, char** options) @nogc nothrow ;

int proj_uom_get_info_from_database(PJ_CONTEXT* ctx, const(char)* auth_name, const(char)* code,
                                    char** out_name, double* out_conv_factor, char** out_category) @nogc nothrow ;

PJ* proj_clone(PJ_CONTEXT* ctx, PJ* obj) @nogc nothrow ;

PJ_OBJ_LIST* proj_create_from_name(PJ_CONTEXT* ctx, const(char)* auth_name, const(char)* searchedName,
                              PJ_TYPE* types, size_t typeCount, int approximateMatch,
                              size_t limitResultCount, char** options) @nogc nothrow ;

PJ_TYPE proj_get_type(const(PJ)* obj) @nogc nothrow ;

int proj_is_deprecated(const(PJ)* obj) @nogc nothrow ;

PJ_OBJ_LIST* proj_get_non_deprecated(PJ_CONTEXT* ctx, const(PJ)* obj) @nogc nothrow ;

int proj_is_equivalent_to(const(PJ)* obj, const(PJ)* other, PJ_COMPARISON_CRITERION criterion) @nogc nothrow ;

int proj_is_crs(const(PJ)* obj) @nogc nothrow ;

const(char)* proj_get_name(const(PJ)* obj) @nogc nothrow ;

const(char)* proj_get_id_auth_name(const(PJ)* obj, int index) @nogc nothrow ;

const(char)* proj_get_id_code(const(PJ)* obj, int index) @nogc nothrow ;

int proj_get_area_of_use(PJ_CONTEXT* ctx, const(PJ)* obj, double* out_west_lon_degree,
                         double* out_south_lat_degree, double* out_east_lon_degree,
                         double* out_north_lat_degree, char** out_area_name) @nogc nothrow ;

const(char)* proj_as_wkt(PJ_CONTEXT* ctx, const(PJ)* obj, PJ_WKT_TYPE type, char** options) @nogc nothrow ;

const(char)* proj_as_proj_string(PJ_CONTEXT* ctx, const(PJ)* obj,
                                 PJ_PROJ_STRING_TYPE type, char** options) @nogc nothrow ;

PJ* proj_get_source_crs(PJ_CONTEXT* ctx, const(PJ)* obj) @nogc nothrow ;

PJ* proj_get_target_crs(PJ_CONTEXT* ctx, const(PJ)* obj) @nogc nothrow ;

PJ_OBJ_LIST* proj_identify(PJ_CONTEXT* ctx, const(PJ)* obj, const(char)* auth_name, char** options,
                           int** out_confidence) @nogc nothrow ;

void proj_int_list_destroy(int* list);

PROJ_STRING_LIST proj_get_authorities_from_database(PJ_CONTEXT* ctx) @nogc nothrow ;

PROJ_STRING_LIST proj_get_codes_from_database(PJ_CONTEXT* ctx, const(char)* auth_name, PJ_TYPE type,
                                              int allow_deprecated) @nogc nothrow ;

PROJ_CRS_LIST_PARAMETERS* proj_get_crs_list_parameters_create() @nogc nothrow ;

void proj_get_crs_list_parameters_destroy(PROJ_CRS_LIST_PARAMETERS* params) @nogc nothrow ;

PROJ_CRS_INFO** proj_get_crs_info_list_from_database(PJ_CONTEXT* ctx, const(char)* auth_name,
                                                     const(PROJ_CRS_LIST_PARAMETERS)* params,
                                                     int* out_result_count) @nogc nothrow ;

void proj_crs_info_list_destroy(PROJ_CRS_INFO** list);

struct PJ_OPERATION_FACTORY_CONTEXT;

PJ_OPERATION_FACTORY_CONTEXT* proj_create_operation_factory_context(PJ_CONTEXT* ctx,
                                                                    const(char)* authority)  @nogc nothrow ;

void proj_operation_factory_context_destroy(PJ_OPERATION_FACTORY_CONTEXT* ctx) @nogc nothrow ;

void proj_operation_factory_context_set_desired_accuracy(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    double accuracy) @nogc nothrow ;

void proj_operation_factory_context_set_area_of_interest(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    double west_lon_degree, double south_lat_degree,
                                    double east_lon_degree, double north_lat_degree) @nogc nothrow ;

void proj_operation_factory_context_set_crs_extent_use(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    PROJ_CRS_EXTENT_USE use) @nogc nothrow ;

void proj_operation_factory_context_set_spatial_criterion(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    PROJ_SPATIAL_CRITERION criterion) @nogc nothrow ;

void proj_operation_factory_context_set_grid_availability_use(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    PROJ_GRID_AVAILABILITY_USE use) @nogc nothrow ;

void proj_operation_factory_context_set_use_proj_alternative_grid_names(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    int usePROJNames) @nogc nothrow ;

void proj_operation_factory_context_set_allow_use_intermediate_crs(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    PROJ_INTERMEDIATE_CRS_USE use) @nogc nothrow ;

void proj_operation_factory_context_set_allowed_intermediate_crs(PJ_CONTEXT* ctx,
                                    PJ_OPERATION_FACTORY_CONTEXT* factory_ctx,
                                    char** list_of_auth_name_codes) @nogc nothrow ;

PJ_OBJ_LIST* proj_create_operations(PJ_CONTEXT* ctx, const(PJ)* source_crs, const(PJ)* target_crs,
                                    const(PJ_OPERATION_FACTORY_CONTEXT)* operationContext) @nogc nothrow ;

int proj_list_get_count(const(PJ_OBJ_LIST)* result) @nogc nothrow ;

PJ* proj_list_get(PJ_CONTEXT* ctx, const(PJ_OBJ_LIST)* result, int index) @nogc nothrow ;

void proj_list_destroy(PJ_OBJ_LIST* result) @nogc nothrow ;

PJ* proj_crs_get_geodetic_crs(PJ_CONTEXT* ctx, const(PJ)* crs) @nogc nothrow ;

PJ* proj_crs_get_horizontal_datum(PJ_CONTEXT* ctx, const(PJ)* crs) @nogc nothrow ;

PJ* proj_crs_get_sub_crs(PJ_CONTEXT* ctx, const(PJ)* crs, int index) @nogc nothrow ;

PJ* proj_crs_get_datum(PJ_CONTEXT* ctx, const(PJ)* crs) @nogc nothrow ;

PJ* proj_crs_get_coordinate_system(PJ_CONTEXT* ctx, const(PJ)* crs) @nogc nothrow ;

PJ_COORDINATE_SYSTEM_TYPE proj_cs_get_type(PJ_CONTEXT* ctx, const(PJ)* cs) @nogc nothrow ;

int proj_cs_get_axis_count(PJ_CONTEXT* ctx, const(PJ)* cs) @nogc nothrow ;

int proj_cs_get_axis_info(PJ_CONTEXT* ctx, const(PJ)* cs, int index, char** out_name,
                          char** out_abbrev, char** out_direction, double* out_unit_conv_factor,
                          char** out_unit_name, char** out_unit_auth_name, char** out_unit_code) @nogc nothrow ;

PJ* proj_get_ellipsoid(PJ_CONTEXT* ctx, const(PJ)* obj) @nogc nothrow ;

int proj_ellipsoid_get_parameters(PJ_CONTEXT* ctx, const(PJ)* ellipsoid,
                                  double* out_semi_major_metre, double* out_semi_minor_metre,
                                  int* out_is_semi_minor_computed, double* out_inv_flattening) @nogc nothrow ;

PJ* proj_get_prime_meridian(PJ_CONTEXT* ctx, const(PJ)* obj) @nogc nothrow ;

int proj_prime_meridian_get_parameters(PJ_CONTEXT* ctx, const(PJ)* prime_meridian,
                                       double* out_longitude, double* out_unit_conv_factor,
                                       char** out_unit_name) @nogc nothrow ;

PJ* proj_crs_get_coordoperation(PJ_CONTEXT* ctx, const(PJ)* crs) @nogc nothrow ;

int proj_coordoperation_get_method_info(PJ_CONTEXT* ctx, const(PJ)* coordoperation,
                                        char** out_method_name, char** out_method_auth_name,
                                        char** out_method_code) @nogc nothrow ;

int proj_coordoperation_is_instantiable(PJ_CONTEXT* ctx, const(PJ)* coordoperation) @nogc nothrow ;

int proj_coordoperation_has_ballpark_transformation(PJ_CONTEXT* ctx, const(PJ)* coordoperation) @nogc nothrow ;

int proj_coordoperation_get_param_count(PJ_CONTEXT* ctx, const(PJ)* coordoperaton) @nogc nothrow ;

int proj_coordoperation_get_param_index(PJ_CONTEXT* ctx, const(PJ)* coordoperation,
                                        const(char)* name) @nogc nothrow ;

int proj_coordoperation_get_param(PJ_CONTEXT* ctx, const(PJ)* coordoperation, int index,
                                 char** out_name, char** out_auth_name, char** out_code,
                                 double* out_value, char** out_value_string,
                                 double* out_unit_conv_factor, char** out_unit_name,
                                 char** out_unit_auth_name, char** out_unit_code,
                                 char** out_unit_category) @nogc nothrow ;

int proj_coordoperation_get_grid_used_count(PJ_CONTEXT* ctx, const(PJ)* coordoperation) @nogc nothrow ;

int proj_coordoperation_get_grid_used(PJ_CONTEXT* ctx, const(PJ)* coordoperation, int index,
                                 char** out_short_name, char** out_full_name,
                                 char** out_package_name, char** out_url, int* out_direct_download,
                                 int* out_open_license, int* out_available) @nogc nothrow ;

double proj_coordoperation_get_accuracy(PJ_CONTEXT* ctx, const(PJ)* obj) @nogc nothrow ;

int proj_coordoperation_get_towgs84_values(PJ_CONTEXT* ctx, const(PJ)* coordoperation,
                                           double* out_values, int value_count,
                                           int emit_error_if_incompatible) @nogc nothrow ;

}
