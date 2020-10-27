% Brownian Motion
% Version 1.0.0
%
% Auxiliary functions
%   brenner     - diffusion coefficient in front of a surface
%
% Functions for the analysis of Brownian motion
%   blocking    - blocking of data
%   pdf         - potential density function
%   potential   - potential
%   msd         - mean square displacement
%   acf         - autocorrelation function
%   ccf         - cross-correlation function
%   psd         - power spectral density
%   driffusion  - drift and diffusion calcualtion in 1D
%
% Object to define a Brownian motion
%   BrownianMotion                      - (Abstract) : Brownian motion
%   BrownianMotion1DFree                - < <a href="matlab:help BrownianMotion">BrownianMotion</a> : 1D free Brownian motion
%   BrownianMotion1DFreeInertial        - < <a href="matlab:help BrownianMotion1DFree">BrownianMotion1DFree</a> : 1D free Brownian motion with inertia
%   BrownianMotion1DOT                  - < <a href="matlab:help BrownianMotion1DFree">BrownianMotion1DFree</a> : 1D free Brownian motion in an optical tweezers
%   BrownianMotion1DDoubleWell          - < <a href="matlab:help BrownianMotion1DFree">BrownianMotion1DFree</a> : 1D Brownian motion in a double-well potential
%   BrownianMotion1DDiffusionGradient   - < <a href="matlab:help BrownianMotion1DFree">BrownianMotion1DFree</a> : 1D free Brownian motion in a diffusion gradient
%   BrownianMotion2DFree                - < <a href="matlab:help BrownianMotion">BrownianMotion</a> : 2D free Brownian motion
%   BrownianMotion2DOT                  - < <a href="matlab:help BrownianMotion2DFree">BrownianMotion2DFree</a> : 2D free Brownian motion in an optical tweezers
%   BrownianMotion2DRotational          - < <a href="matlab:help BrownianMotion2DFree">BrownianMotion2DFree</a> : 2D free Brownian motion in a rotational force field
%   BrownianMotion3DFree                - < <a href="matlab:help BrownianMotion">BrownianMotion</a> : 3D free Brownian motion
%   BrownianMotion3DOT                  - < <a href="matlab:help BrownianMotion3DFree">BrownianMotion1DFree</a> : 3D free Brownian motion in an optical tweezers
% 
% Examples
%   example_brownianmotion1dfree                - Example to demonstrate the use of <a href="matlab:help BrownianMotion1DFree">BrownianMotion1DFree</a>
%   example_brownianmotion1dfreeinertial        - Example to demonstrate the use of <a href="matlab:help BrownianMotion1DFreeInertial">BrownianMotion1DFreeInertial</a>
%   example_brownianmotion1dot                  - Example to demonstrate the use of <a href="matlab:help BrownianMotion1DOT">BrownianMotion1DOT</a>
%   example_brownianmotion1ddoublewell          - Example to demonstrate the use of <a href="matlab:help BrownianMotion1DDoubleWell">BrownianMotion1DDoubleWell</a>
%   example_brownianmotion1ddiffusiongradient   - Example to demonstrate the use of <a href="matlab:help BrownianMotion1DDiffusionGradient">BrownianMotion1DDiffusionGradient</a>
%   example_brownianmotion1dfree                - Example to demonstrate the use of <a href="matlab:help BrownianMotion1DFree">BrownianMotion1DFree</a>
%   example_brownianmotion2dfree                - Example to demonstrate the use of <a href="matlab:help BrownianMotion2DFree">BrownianMotion2DFree</a>
%   example_brownianmotion2dot                  - Example to demonstrate the use of <a href="matlab:help BrownianMotion2DOT">BrownianMotion2DOT</a>
%   example_brownianmotion3dfree                - Example to demonstrate the use of <a href="matlab:help BrownianMotion3DFree">BrownianMotion3DFree</a>
%   example_brownianmotion3dot                  - Example to demonstrate the use of <a href="matlab:help BrownianMotion3DOT">BrownianMotion3DOT</a>
%
% See also OTS.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

clc
help bm