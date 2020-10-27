function [ccf,tau,ccf_err,ccf_vec] = ccf(dt,Vx,Vy,varargin)
% CCF   Cross-correlation function
%
% [CCF,TAU,CCFerr,CCFvec] = CCF(DT,VX,VY) calculates the cross-correlation function CCF
%   between signals VX and VY sampled at time intervals DT.
%   VX and VY can be matrices containing several signals (one per row, all the same length).
%   TAU are the values of the time delays.
%   CCFerr is the error calculated as standard deviation of the CCF value.
%   CCFvec is the set of all CCF calculated
%
% [CCF,TAU,CCFerr,CCFvec] = CCF(DT,VX,VY,'PropertyName',PropertyValue) permits
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
ccf_vec = zeros(maxlags,size(Vx,2));
for n = 1:1:size(Vx,2)
    tmp = xcov(Vx(:,n),Vy(:,n),maxlags,'unbiased');
    tmp = fftshift(tmp);
    ccf_vec(1:1:maxlags,n) = tmp(1:1:maxlags);
end
    
tau = dt*[0:1:maxlags-1]';

ccf = mean(ccf_vec,2);
ccf_err = std(ccf_vec,0,2);
