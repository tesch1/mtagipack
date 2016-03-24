%! @file mtagipack_setup.m
% usage:
% put the containing directory into your matlab path:
% \code
% >> path(path, '...path/to/mtagipack/');
% \endcode
% compile the mex functions:
% \code
% >> which mtagipack_setup
% ...path/to/mtagipack/mtagipack_setup.m
% >> mtagipack_setup
% setup mex in .../mtagipack
% *** if you get errors, make sure to remove "-ansi" from mexopts.sh (see mex -setup) ***
% compiling mexLoadAgilentTraces...
% compiling mexLoadAgilentParams...
% return to /wherever you were/
% >> mexLoadAgilentTraces
% Error using mexLoadAgilentTraces
% usage: [traces mainh blockh] = mexLoadAgilentTraces(filename)
% 
% \endcode
%
debugflags = '-v -g';
debugflags = '';

oldpath=pwd;
[fpath, name, exte] = fileparts(which('mtagipack_setup'));

cd (fpath);

disp(['setup mex in ' fpath]);
disp(['*** if you get errors, make sure to remove "-ansi" from mexopts.sh (see mex -setup) ***']);

disp('compiling mexLoadAgilentTraces...')
mex -largeArrayDims mexLoadAgilentTraces.c mexCommons.c varian.c

disp('compiling mexSaveAgilentTraces...')
mex -largeArrayDims mexSaveAgilentTraces.c mexCommons.c varian.c

disp('compiling mexLoadAgilentParams...')
mex -largeArrayDims mexLoadAgilentParams.c cprocpar.c

disp(['return to ' oldpath]);
cd (oldpath);
