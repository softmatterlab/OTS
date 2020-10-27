function [acf,tau,acf_err,acf_vec] = acf(dt,Vx,varargin)
% ACF   Autocorrelation function
%
% [ACF,TAU,ACFerr,ACFvec] = ACF(DT,VX) calculates the autocorrelation function ACF
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   ACFerr is the error calculated as standard deviation of the ACF value.
%   ACFvec is the set of all ACF calculated
%
% [ACF,TAU,ACFerr,ACFvec] = ACF(DT,VX,'PropertyName',PropertyValue) permits
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
acf_vec = zeros(maxlags,size(Vx,2));
for n = 1:1:size(Vx,2)
    tmp = xcov(Vx(:,n),maxlags,'unbiased');
    tmp = fftshift(tmp);
    acf_vec(1:1:maxlags,n) = tmp(1:1:maxlags);
end
    
tau = dt*[0:1:maxlags-1]';

acf = mean(acf_vec,2);
acf_err = std(acf_vec,0,2);
