function [psd,f,psd_err,psd_vec,f_vec] = psd(dt,Vx,varargin)
% PSD   Power spectral density
%
% [PSD,F,PSDerr,PSDvec,Fvec] = PSD(DT,VX) calculates the power spectral density PSD
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per row, all the same length).
%   F are the values of the frequencies.
%   PSDerr is the PSD error.
%   PSDvec is the set of all PSD calculated
%
% [PSD,F,PSDerr,PSDvec,Fvec] = PSD(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       Fmin            -   minimum frequency (default = 1/total acquisition time)
%       Fmax            -   maximum frequency (defauls = 1/DT)
%       Blocking        -   blocking (default = 'none' | 'lin' | 'log')
%       BinsNumber      -   number of bins for blocking (default = 100)

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

% Minimum frequency
fmin = 1/(dt*size(Vx,1));
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'fmin')
        fmin = varargin{n+1};
    end
end

% Maximum frequency
fmax = 1/dt;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'fmax')
        fmax = varargin{n+1};
    end
end

% Blocking: none (default) | lin | log
blocking = 'none';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'blocking')
        blocking = varargin{n+1};
    end
end

% Number of bins for blocking
binsnumber = 100;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'binsnumber')
        binsnumber = varargin{n+1};
    end
end

% Analysis
fs = 1/dt;
T = dt*size(Vx,1); % measurement time

ft_vec = dt*fft(Vx);
psd_vec = ft_vec.*conj(ft_vec)/T;

fNyq = fs/2; % Nyquist frequency
f = [1:1:size(Vx)]'/T;

% Only data up to the Nyquist frequency
% and between fmin and fmax
ind = find(f<=fNyq & f>=fmin & f<=fmax);
f = f(ind);
f_vec = f;
psd_vec = psd_vec(ind,:);

psd = mean(psd_vec,2);
psd_err = std(psd_vec,0,2);

% Blocking
if strcmpi(blocking,'lin')
    
    % Centers of bins
    binscenters = [fmin:(fmax-fmin)/binsnumber:fmax];
    
    % Widths of Bins
    binswidths(1) = (binscenters(2)-binscenters(1))/2;
    binswidths(2:length(binscenters)-1) = (binscenters(3:end)-binscenters(1:end-2))/4;
    binswidths(length(binscenters)) = (binscenters(end)-binscenters(end-1))/2;

    f_b = zeros(size(binscenters));
    psd_b = zeros(size(binscenters,size(psd_vec,2)));
    psd_b_err = zeros(size(binscenters,size(psd_vec,2)));
    for bin = 1:1:length(binscenters)

        ind = find( f>=binscenters(bin)-binswidths(bin) & f<=binscenters(bin)+binswidths(bin) );
        
        f_b(bin) = mean(f(ind));
        psd_b(bin) = mean(psd(ind));
        psd_err_b(bin) = mean(psd(ind))/sqrt(length(ind));
    end
    
    f = f_b;
    psd = psd_b;
    psd_err = psd_err_b;

elseif strcmpi(blocking,'log')

    % Centers of bins
    binscenters = [log(fmin):(log(fmax)-log(fmin))/binsnumber:log(fmax)];
    binscenters = exp(binscenters);
    
    % Widths of Bins
    binswidths(1) = (binscenters(2)-binscenters(1))/2;
    binswidths(2:length(binscenters)-1) = (binscenters(3:end)-binscenters(1:end-2))/4;
    binswidths(length(binscenters)) = (binscenters(end)-binscenters(end-1))/2;

    f_b = zeros(size(binscenters));
    psd_b = zeros(size(binscenters,size(psd_vec,2)));
    psd_b_err = zeros(size(binscenters,size(psd_vec,2)));
    for bin = 1:1:length(binscenters)

        ind = find( f>=binscenters(bin)-binswidths(bin) & f<=binscenters(bin)+binswidths(bin) );

        f_b(bin) = mean(f(ind));
        psd_b(bin) = mean(psd(ind));
        psd_err_b(bin) = mean(psd(ind))/sqrt(length(ind));
    end
    
    f = f_b;
    psd = psd_b;
    psd_err = psd_err_b;
    
end