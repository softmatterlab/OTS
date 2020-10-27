function [U,binscenters,U_vec,p,p_err,p_vec] = potential(Vx,T,varargin)
% POTENTIAL   Potential
%
% [U,BC,Uvec,P,Perr,Pvec] = POTENTIAL(VX,T) calculates the potential U
%   associated to signal VX, assuming absolute temperature T.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   BC are the bins centers for the potential.
%   Uvec is the set of all calculated potentials.
%   P is the probability density function (PDF).
%   Perr is the PDF error.
%   Pvec is the set of all PDF calculated.
%
% [U,BC,Uvec,P,Perr,Pvec] = POTENTIAL(VX,T,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       BinsNumber      -   number of bins (default = 50)
%       BinsCenters     -   centers of bins

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

[p,binscenters,p_err,p_vec] = pdf(Vx,varargin{:});

% Analysis
U_vec = -PhysConst.kB*T*log(p_vec);
U = -PhysConst.kB*T*log(p);