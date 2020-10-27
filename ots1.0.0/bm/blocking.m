function [y_b,x_b,y_b_err] = blocking(x,y,err,varargin)
% BLOCKING    Blocking of signal
%
% [Yb,Xb,Yberr] = BLOCKING(X,Y,ERR) blocks the series of data (X,Y) 
%   with errors ERR (standard deviation).
%   Yb are the values of the blocked signals.
%   Xb are the positions of the blocked signals.
%   Tberr are the errors f the blocked signals.
%
% [Yb,Xb,Yberr] = BLOCKING(X,Y,ERR,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       XMin        -   initial value for the blocking
%       XMax        -   final value for the blocking
%       Blocking    -   kind of blocking (default = 'lin' | 'log')
%       BinsNumber  -   number of bins (default = 100)

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

% Minimum
xmin = min(x);
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xmin')
        xmin = varargin{n+1};
    end
end

% Maximum
xmax = max(x);
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'xmax')
        xmax = varargin{n+1};
    end
end

% Blocking: lin (default) | log
blocking = 'lin';
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

% Blocking
if strcmpi(blocking,'log')

    % Centers of bins
    binscenters = [log(xmin):(log(xmax)-log(xmin))/binsnumber:log(xmax)];
    binscenters = exp(binscenters);
    
    % Widths of Bins
    binswidths(1) = (binscenters(2)-binscenters(1))/2;
    binswidths(2:length(binscenters)-1) = (binscenters(3:end)-binscenters(1:end-2))/4;
    binswidths(length(binscenters)) = (binscenters(end)-binscenters(end-1))/2;

    x_b = zeros(size(binscenters));
    y_b = zeros(size(binscenters,size(y,2)));
    y_b_err = zeros(size(binscenters,size(y,2)));
    for bin = 1:1:length(binscenters)

        ind = find( x>=binscenters(bin)-binswidths(bin) & x<=binscenters(bin)+binswidths(bin) );

        x_b(bin) = mean(x(ind));
        y_b(bin) = mean(y(ind));
        y_b_err(bin) = mean(err(ind))/sqrt(length(ind));
    end
        
else % strcmpi(blocking,'lin')
    
    % Centers of bins
    binscenters = [xmin:(xmax-xmin)/binsnumber:xmax];
    
    % Widths of Bins
    binswidths(1) = (binscenters(2)-binscenters(1))/2;
    binswidths(2:length(binscenters)-1) = (binscenters(3:end)-binscenters(1:end-2))/4;
    binswidths(length(binscenters)) = (binscenters(end)-binscenters(end-1))/2;

    x_b = zeros(size(binscenters));
    y_b = zeros(size(binscenters,size(y,2)));
    y_b_err = zeros(size(binscenters,size(y,2)));
    for bin = 1:1:length(binscenters)

        ind = find( x>=binscenters(bin)-binswidths(bin) & x<=binscenters(bin)+binswidths(bin) );
        
        x_b(bin) = mean(x(ind));
        y_b(bin) = mean(y(ind));
        y_b_err(bin) = mean(err(ind))/sqrt(length(ind));
    end
    
end