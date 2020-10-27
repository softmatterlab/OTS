% Series of examples to demonstrate the use of LGF.
%
% See also LGF.

%   Author: S. Masoumeh Mousavi
%   Revision: 1.0.0  
%   Date: 2015/01/01

example('Use of LGF')

%% LOAD LGF
exampletitle('LOAD LGF')

examplecode('lgf = LGF');

examplewait()

%% LOGARITM OF FACTORIAL 
exampletitle('LOGARITM OF FACTORIAL')

examplecode('lgf.getfact((2))');
examplecode('lgf.getfact((10))');
examplecode('lgf.getfact(([2 3 4 5]))');
examplecode('lgf.getfact(([2; 3; 4; 5]))');

examplecode('lgf.getfact(([2 3 4 5 ; 8 11 3 2; 0 4 6 9]))');