## A D binding for PROJ
This library offers D bindings for the PROJ library. PROJ provides easy-to-use
cartographic transformations. Details about the PROJ project can be found at
[PROJ4.org](https://PROJ4.org/). This binding is for PROJ version 6.0.

This packages offers two APIs. The first is a translation of the PROJ.h header file, and
documentation can be found [here](https://PROJ4.org/development/reference/functions.html).
This API can be imported with `import proj;`. An example program using this API can be found
in [examples/proj_example.d](examples/proj_example.d). It is ill-tested, especially the newer API
functions, so use with caution and submit bug reports!

The second API uses some of the features of D to make using the library easier. It is imported with
`import crs;`. An example of using this API is at [examples/crs_example.d](examples/crs_example.d).
The central object is the `Projection` class. It is constructed either with a
string like with the `proj_create` function, or with an associative array of strings indexed
with members of the `PK` enum. The second style of construction will feel familiar to those who have
used Python interfaces to PROJ like cartopy or pyproj. There are many subclasses of `Projection`
--one for each projection listed on the
[PROJ website](https://PROJ4.org/operations/projections/index.html). The full list is in
available_classes.txt. Using these promotes more readable code.

### Install
These instructions are for Linux. If someone can get it to work on Windows/OS X, please let
me know so we can add those operating systems to this README!

First you need PROJ 5 at minimum. Linux distro package managers tend to have older versions, so
I recommend using conda to install the latest version (and so does PROJ itself). Instructions
are [here](https://proj4.org/install.html#conda). If you do use conda, be sure to add your conda
environment's `lib` to your `LD_LIBRARY_PATH` environment variable, so your linker can find it,
and your conda environment's `lib/pkgconfig` to your `PKG_CONFIG_PATH` environment variable,
so dub can find it.

For now you can install these bindings with a `git clone https://github.com/weech/dproj.git` and
`dub add-local dproj`.

