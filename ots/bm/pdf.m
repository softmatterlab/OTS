function [p,binscenters,p_err,p_vec] = pdf(Vx,varargin)
% PDF   Probability density function
%
% [P,BC,Perr,Pvec] = PSD(VX) calculates the probability density function P
%   associated to signal VX.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   BC are the bins centers for the potential.
%   Perr is the PDF error.
%   Pvec is the set of all PDF calculated.
%
% [P,BC,Perr,Pvec] = PSD(VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       BinsNumber      -   number of bins (default = 50)
%       BinsCenters     -   centers of bins

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%   Author: Giovanni Volpe
%   Revision 1.0.1
%   Date: 2016/08/01
%   lines 40-44 - corrected bug (identified by Mark Cronin-Golomb)

% Number of Bins
binsnumber = 50;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'binsnumber')
        binsnumber = varargin{n+1};
    end
end

% Centers of Bins
binscenters = [min(Vx):(max(Vx)-min(Vx))/binsnumber:max(Vx)];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'binscenters')
        binscenters = varargin{n+1};
        binsnumber = length(binscenters);
    end
end

% Analysis
p_vec = hist(Vx,binscenters);
% START CHANGE v1.0.1
if size(p_vec,1)==1
    p_vec = p_vec';
end
% END CHANGE
p = mean(p_vec,2)
p_err = std(p_vec,0,2)
    