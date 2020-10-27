function [c,d,nc,nd,binscenter,binswidth] = driffusion(dt,x,varargin)
% DRIFFUSION    calculates drift and diffusion
%
% [C,D,NC,ND,BC,BW] = DRIFFUSION(DT,VX) calculates the drift C and diffusion D 
%   of signal VX sampled at time intervals DT.
%   VX can be a matrix containing several signals (one per COLUMN, all the same length).
%   NC are the number of sampels per bin used to calculate the drift.
%   ND are the number of sampels per bin used to calculate the diffusion.
%   BC are the bins centers.
%   BW are the bins width.
%
% [C,D,NC,ND,BC,BW] = DRIFFUSION(DT,VX,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       Delay           -   delay (default = 1)
%       DelayDrift      -   delay drift (default = delay)
%       DelayDiffusion	-   delay diffusion (default = delay)
%       BinsCenter      -   centers of bins
%       BinsWidth       -   widths of bins
%       Unbiased        -   whether to use unbiased (default) or biased drift in diffusion calculation

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%   Author: Giovanni Volpe
%   Revision 1.0.1
%   Date: 2016/08/01
%   lines 82-86 & 107-112 - corrected bug (addition by Agnese Callegari)

% Delay
delay = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'delay')
        delay = varargin{n+1};
    end
end

% Delay Drift
delayc = delay;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'delaydrift')
        delayc = varargin{n+1};
    end
end

% Delay Diffusion
delayd = delay;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'delaydiffusion')
        delayd = varargin{n+1};
    end
end

% Center of Bins
binscenter = [min(x):(max(x)-min(x))/50:max(x)];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'binscenter')
        binscenter = varargin{n+1};
    end
end

% Width of Bins
binswidth(1) = (binscenter(2)-binscenter(1))/2;
binswidth(2:length(binscenter)-1) = (binscenter(3:end)-binscenter(1:end-2))/4;
binswidth(length(binscenter)) = (binscenter(end)-binscenter(end-1))/2;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'binswidth')
        if length(varargin{n+1})==1
            binswidth = varargin{n+1}*ones(size(binscenter));
        else
            binswidth = varargin{n+1};
        end
    end
end

% Unbiased (subtract drift in diffusion calculation)
unbiased = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'unbiased')
        unbiased = varargin{n+1};
    end
end

% Analysis
if (delayc==delayd)
    % START CHANGE v1.0.1
    x0 = x(1:end-delay,:); % initial states
    dx = x(delay+1:end,:) - x(1:end-delay,:); % displacements
    % END CHANGE

    c = zeros(size(binscenter));
    d = zeros(size(binscenter));
    nc = zeros(size(binscenter));
    nd = zeros(size(binscenter));
    for bin = 1:1:length(binscenter)
        bin_start = binscenter(bin)-binswidth(bin);
        bin_end = binscenter(bin)+binswidth(bin);

        cdisplacements = dx( (x0>bin_start) & (x0<bin_end) ); % conditional displacements

        c(bin) = mean(cdisplacements)/(delay*dt); % Drift
        if unbiased
            d(bin) = mean( (cdisplacements-c(bin)*delay*dt).^2 )/(2*delay*dt); % Diffusion
        else
            d(bin) = mean( cdisplacements.^2 )/(2*delay*dt); % Diffusion
        end
        nc(bin) = length(cdisplacements);
        nd(bin) = length(cdisplacements);
    end
else
    % START CHANGE v1.0.1
    x0c = x(1:end-delayc,:); % initial states drift
    dxc = x(delayc+1:end,:) - x(1:end-delayc,:); % displacements drift
    x0d = x(1:end-delayd,:); % initial states diffusion
    dxd = x(delayd+1:end,:) - x(1:end-delayd,:); % displacements diffusion
    % END CHANGE

    c = zeros(size(binscenter));
    d = zeros(size(binscenter));
    nc = zeros(size(binscenter));
    nd = zeros(size(binscenter));
    for bin = 1:1:length(binscenter)
        bin_start = binscenter(bin)-binswidth(bin);
        bin_end = binscenter(bin)+binswidth(bin);

        cdisplacementsc = dxc( (x0c>bin_start) & (x0c<bin_end) ); % conditional displacements drift
        cdisplacementsd = dxd( (x0d>bin_start) & (x0d<bin_end) ); % conditional displacements diffusion

        c(bin) = mean(cdisplacementsc)/(delayc*dt); % Drift
        if unbiased
            d(bin) = mean( (cdisplacementsd-c(bin)*delayd*dt).^2 )/(2*delayd*dt); % Diffusion
        else
            d(bin) = mean( cdisplacementsd.^2 )/(2*delayd*dt); % Diffusion            
        end
        nc(bin) = length(cdisplacementsc);
        nd(bin) = length(cdisplacementsd);
    end    
end