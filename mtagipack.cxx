/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 */
/* This file is JUST for generating pages in Doxygen's output, dont bother to compile it. */
/** \file
 * \brief       Doxygen Documentation Mainpage. */

/**
 * \mainpage    mtagipack
 * \brief       Collection of utilities to read scanner data into matlab & octave.
 *
 * \section     index_contents          The contents of mtagipack:
 *
 * - mtagipack_setup.m            Build the mex files before using them in matlab
 * - mexLoadAgilentParams.c       Load agilent 'procpar' variables into a struct
 * - mexLoadAgilentTraces.c       Load agilent(/varian) fid file data
 * - mexSaveAgilentTraces.c       Save agilent(/varian) fid file data
 *
 * \section     index_obtain            Where to find mtagipack
 *
 * The package can be found here: 
 *      <a href="http://www.cmrr.umn.edu/~tesch/downloads"> http://www.cmrr.umn.edu/~tesch/downloads</a>
 * Maybe it has moved, in which case look around or ask me: 
 *
 * \section     index_setup             Setup
 * 
 * It maybe necessary to use 'mex -setup' in matlab to adjust your mex compilation flags
 * to REMOVE the '-ansi' flag.  Practically nobody has written "ansi C" since the mid-80s.
 *
\verbatim
tesch$ tar xzvf mtagipack-1.6.7.tgz 
tesch$ matlab -nojvm
...
>> path(path,'.../.../path/to/mtagipack-1.6.7');
>> which mtagipack_setup
.../mtagipack-1.6.7/mtagipack_setup.m
>> mtagipack_setup
setup mex in .../mtagipack-1.6.7
*** if you get errors, make sure to remove "-ansi" from mexopts.sh (see mex -setup) ***
compiling mexLoadAgilentTraces...
compiling mexLoadAgilentParams...
return to /wherever you were/
>> mexLoadAgilentTraces
Error using mexLoadAgilentTraces
usage: [traces mainh blockh] = mexLoadAgilentTraces(filename)

>> [traces mainh blockh] = mexLoadAgilentTraces('.../data2/djswiDD_16kview_test_dataset.fid/');  
read blocks:8 traces:24576 samples:6291456 channel:0/1
>> size(traces)

ans =

         256        3072           8

>> mainh

mainh = 

      nblocks: 8
      ntraces: 3072
           np: 256
       ebytes: 4
       tbytes: 2048
       bbytes: 6291484
       status: 73
      vers_id: 0
    nbheaders: 1

>> 

\endverbatim
 *
 * \section  index_changelog         Version ChangeLog.
 *      Release \ref changelog
 */

/**
 * \page        version                 AGIPACK-VERSION
 * git ident: $Id: c797996b9601dcd5ec6d887202fc49f2c1db53db $
 */
