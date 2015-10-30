# Introduction #

The PD5 software package is an open-source collection of library classes and utilities for primer design. PD5 is written in object-oriented C++ with an emphasis on classes suitable for the design of bespoke primer design programs.


# Details #

General methods in the current release include primer candidate generation; primer candidate analysis, which incorporates a novel scoring method; and primer candidate selection using a multi objective optimisation approach. Specific methods for sequence searching; sequence location and the analysis of hairpin formation, dimer formation, non-specific binding and primary structure are also included.


PD5 is published under the terms of the GNU General Public License and the GNU Library General Public License which guarantee access to [source-code](https://code.google.com/p/primer-design/source/browse/) and allow redistribution and modification.

# Getting Started #

For a quick introduction to PD5 primer design download and unzip `PD5_x.x+apps+manual.zip` (where x.x is the version number) from the [downloads](http://code.google.com/p/primer-design/downloads/list) page and change directory to `/primer-design/basic_app/`. There you will find a basic PD5 application named `basic_app.cpp`.

## PD5 Basic Application ##

`basic_app.cpp` is a simple application for primer design that can be used as a starting point for more specific apps.

To build in directory `/primer-design/basic_app/` use command: -
```
   ./build_basic_app
```
To test: -
```
   ./example_command
```

This passes the sequence for the Arabidopsis speechless gene (At5g53210) in the command line as an example template.

The resulting output shows lists of forward and reverse candidate primers before and after sorting/selection, and then details of the best 6 candidate pairs.

_This information is also in `/primer-design/basic_app/Readme.txt`._

## PD5 Full Build ##

For the full build of the PD5 library and all its associated apps, use the Makefiles provided.

# PD5 Galaxy Plugin #

PD5's command line interface has a Galaxy plugin. Go into Galaxy's tools folder and checkout PD5 (so you should have "galaxy/galaxy-dist/tools/primer-design"). We have included the config file galaxy\_pd5\_cli\_tool.xml, so all you need to do is add the pd5\_cli executable to your PATH and add the following to your Galaxy tool-conf.xml:

```
  <section name="PD5 primer design" id="primerdesign">
    <tool file="primer-design/galaxy_pd5_cli_tool.xml" />
  </section>
```

# PD5 Reference #
A more detailed reference manual is available from the [downloads](http://code.google.com/p/primer-design/downloads/list) page. Download and unzip _Reference manual html.zip_ and click on /html/index.html for the main reference manual page, which should open in your default browser.