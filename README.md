# mtagipack
c and matlab routines for reading varian/agilent nmr/mri data files

# The contents of mtagipack:
- mtagipack_setup.m            Build the mex files before using them in matlab
- mexLoadAgilentParams.c       Load agilent 'procpar' variables into a struct
- mexLoadAgilentTraces.c       Load agilent(/varian) fid file data
- mexSaveAgilentTraces.c       Save agilent(/varian) fid file data

# Where to find mtagipack
- https://github.com/tesch1/mtagipack

# Setup

It maybe necessary to use 'mex -setup' in matlab to adjust your mex compilation flags
to REMOVE the '-ansi' flag.  Practically nobody has written "ansi C" since the mid-80s.

Get the tarball/zip, and do the following:

```bash
tesch$ tar xzvf mtagipack-1.6.7.tgz 
tesch$ matlab -nojvm
```
```matlab
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
```
