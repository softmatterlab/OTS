function [msd,tau,msd_err,msd_vec] = msd(dt,Vx,varargin)
% MSD   Mean Square Displacement
%
% [MSD,TAU,MSDerr,MSDvec] = MSD(DT,VX) calculates the mean square displacement MSD
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   MSDerr is the error calculated as standard deviation of the MSD value.
%   MSDvec is the set of all MSD calculated
%
% [MSD,TAU,MSDerr,MSDvec] = MSD(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       TauMax      -   Maximum delay to be calcualted (default = +Inf)

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

% Maximum delay
taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    end
end
if isfinite(taumax)
    maxlags = ceil(taumax/dt);
else
    maxlags = size(Vx,1);
end

% Analysis
msd_vec = zeros(maxlags,size(Vx,2));
for m = 1:1:maxlags-1
    msd_vec(m+1,:) = mean( (Vx(m+1:end,:)-Vx(1:end-m,:)).^2 , 1);
end

tau = dt*[0:1:maxlags-1]';

msd = mean(msd_vec,2);
msd_err = std(msd_vec,0,2);