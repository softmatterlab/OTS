% Dipole Approximation
% Version 1.0.0
%
% Object to define a static dipole
%   StaticDipole        - Static dipole
% 
% Objects to define electromagnetic fields
%   EField              - (Abstract) : Electromagnetic field
%   EFieldPlaneWave     - < <a href="matlab:help EField">EField</a> : Plane wave
%   EFieldFocus         - < <a href="matlab:help EField">EField</a> : Focal field
%
% Object to define an induced dipole
%   InducedDipole       - < <a href="matlab:help EField">EField</a> : Induced dipole
%
% Examples
%   example_staticdipole    - Example to demonstrate the use of <a href="matlab:help StaticDipole">StaticDipole</a>
%   example_efieldplanewave - Example to demonstrate the use of <a href="matlab:help EFieldPlaneWave">EFieldPlaneWave</a>
%   example_efieldfocus     - Example to demonstrate the use of <a href="matlab:help EFieldFocus">EFieldFocus</a>
%   example_induceddipole   - Example to demonstrate the use of <a href="matlab:help InducedDipole">InducedDipole</a>
%
% See also OTS, Beam.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clc
help da