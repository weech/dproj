import std.stdio;

import proj;

void main() {
    /* Create context and projection. Single-threaded programs don't need to create the
        context and can use PJ_DEFAULT_CTX */
    PJ_CONTEXT* ctx = proj_context_create();
    // Note that this string will have to be converted to C-style if not literal
    PJ* projection = proj_create(ctx, "+proj=utm +zone=32 +ellps=GRS80");
    if (proj_errno(projection)) {
        writeln(proj_errno_string(proj_errno(projection)));
    }

    // Create a coordinate object. Note that this API works in radians (crs uses degrees)
    PJ_COORD a = proj_coord(proj_torad(12), proj_torad(55), 0, 0);

    // Transform forwards and backwards
    PJ_COORD b = proj_trans(projection, PJ_FWD, a);
    writefln("easting: %g, northing: %g", b.enu.e, b.enu.n);
    b = proj_trans(projection, PJ_INV, b);
    writefln("longitude: %g, latitude: %g", proj_todeg(b.lp.lam), proj_todeg(b.lp.phi));

    // Clean up
    proj_destroy(projection);
    proj_context_destroy(ctx);
}