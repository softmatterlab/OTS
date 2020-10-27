% Electro-Magnetic Theory 
% Version 1.0.0
%
% Object to define coefficients
%   Coefficients	- Coefficients 
%
% Objects to define Clebsch-Gordan coefficients
%   CG      - < handle : Clebsch-Gordan coefficients
%   LGF     - < handle : logarithm of factorial
% 
% Objects to define incident field
%   IncidentField               - (Abstract) < handle : Incident field
%   IncidentFieldPlaneWave      - < <a href="matlab:help IncidentField">IncidentField</a> < handle : Incident field plane wave
%   IncidentFieldFocusedBeam    - < <a href="matlab:help IncidentField">IncidentField</a> < handle : Incident field focused beam
%
% Object to define a T-matrix
%   TMatrix             - (Abstract) < handle : T-matrix
%   TMatrixSphere       - < <a href="matlab:help TMatrix">TMatrix</a> < handle : T-Matrix of spherical particle
%   TmatrixCluster      - < <a href="matlab:help TMatrix">TMatrix</a> < handle : T-Matrix of cluster
%   TmatrixInclusions   - < <a href="matlab:help TMatrix">TMatrix</a> < handle : TMatrix inclusion
%
% Examples
%   example_coefficients                - Example to demonstrate the use of <a href="matlab:help Coefficients">Coefficients</a>
%   example_cg                          - Example to demonstrate the use of <a href="matlab:help CG">CG</a>
%   example_lgf                         - Example to demonstrate the use of <a href="matlab:help LGF">LGF</a>
%   example_incidentfieldplanewave      - Example to demonstrate the use of <a href="matlab:help IncidentFieldPlaneWave ">IncidentFieldPlaneWave </a>
%   example_incidentfieldfocusedbeam    - Example to demonstrate the use of <a href="matlab:help IncidentFieldFocusedBeam">IncidentFieldFocusedBeam</a>
%   example_tMatrixsphere               - Example to demonstrate the use of <a href="matlab:help TMatrixSphere">TMatrixSphere</a>
%   example_tmatrixcluster              - Example to demonstrate the use of <a href="matlab:help TmatrixCluster">TmatrixCluster</a>
%   example_tmatrixinclusions           - Example to demonstrate the use of <a href="matlab:help TmatrixInclusions">TmatrixInclusions</a>
%
% See also OTS, mie.

%   Author: S. Masoumeh Mousavi,Agnese Callegari, Giovanni Volpe
%   Revision: 1.0.0
%   Date: 2015/01/01

clc
help emt