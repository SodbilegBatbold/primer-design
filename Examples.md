# Introduction #

The source code for three example applications using the PD5 software library are included in the project source directories. All three applications can be built on a Mac using: -

> make -f Makefile.OSX

or on Linux using: -

> make -f Makefile.linux


# Basic\_app #

basic\_app.cpp is a simple application for primer design that can be used as a starting point for more specific apps.

To build in directory /primer-design/basic\_app/ use command: -

> ./build\_basic\_app
To test: -

> ./example\_command
This passes the sequence for the Arabidopsis speechless gene (At5g53210) in the command line as an example template.

The resulting output shows lists of forward and reverse candidate primers before and after sorting/selection, and then details of the best 6 candidate pairs.

This information is also in /primer-design/basic\_app/Readme.txt.

# PD5\_cli #

This application provides a command line interface to the PD5 modules. Use: -

> ./example\_command

to see an example or

> ./pd5\_cli

without arguments to see a brief manual giving details of the required arguments.

# PD5\_ssr #

This application is a specific use of PD5 and is provided with README files giving more details.