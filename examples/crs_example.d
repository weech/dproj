import std.stdio;

import crs;

void main() {

    // Create the projection
    auto projection = new UniversalTransverseMercator([PK.zone:"32", PK.ellps:"GRS80"]);

    // Transform forwards and backwards
    auto point = projection.transform(12, 55);
    writefln("easting: %g, northing: %g", point[0], point[1]);
    point = projection.transform(point, true);
    writefln("longitude: %g, latitude: %g", point[0], point[1]);
}