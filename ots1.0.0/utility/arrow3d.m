function h = arrow3d(x1,y1,z1,x2,y2,z2,varargin)
% ARROW3D   Draw a line with an arrowhead in 3D
%
% H = ARROW3D(X1,Y1,Z1,X2,Y2,Z2) draws arrow(s) from (X1,Y1,Z1) to (X2,Y2,Z2)
%   and returns the graphics handle of the arrow(s).
%   X1,Y1,Z1,X2,Y2,Z2 should be real matrices with the same dimensions.
%
% H = ARROW3D(X1,Y1,Z1,X2,Y2,Z2,'PropertyName',PropertyValue) sets the property
%   PropertyName to PropertyValue. All standard surf properties
%   can be used and also the ARROW2D properties listed below.
%
% ARROW3D properties:
%       Color       Arrow color both edges and faces [default = 'k']
%       StemWidth   Arrow stem width [default = .1]
%       HeadLength  Arrow head length [default = 1]
%       HeadWidth   Arrow head width [default = 1]
%       HeadNode    Arrow head intersection with the arrow stem [default = .5]
%       N           Number of radial sections [default = 32]
%
% See also ARROW2D, SURF.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

% Permission is granted to distribute ARROW3D with the toolboxes for the book
% "Optical Tweezers", by P. H. Jones, O. M. Marago & G. Volpe 
% (Cambridge University Press, 2015).

% Arrow color
color = 'k';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'color')
        color = varargin{n+1};
    end
end

% Arrow stem width
swidth = .1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'stemwidth')
        swidth = varargin{n+1};
    end
end
swidth = swidth.*ones(size(x1));

% Arrow head length
hlength = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headlength')
        hlength = varargin{n+1};
    end
end
hlength = hlength.*ones(size(x1));

% Arrow head width
hwidth = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headwidth')
        hwidth = varargin{n+1};
    end
end
hwidth = hwidth.*ones(size(x1));

% Arrow head intersection with the arrow stem
hnode = .5;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'headnode')
        hnode = varargin{n+1};
    end
end
hnode = hnode.*ones(size(x1));

% Number of equally spaced points around the arrow
N = 32;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'n')
        N = varargin{n+1};
    end
end
N = N.*ones(size(x1));

% Prepares coordiantes in a standard format
slength = sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2); % length
theta = acos((z2-z1)./slength) ; % polar angle 
phi = atan2(y2-y1,x2-x1); % azimuthal angle

h = [];
for m = 1:1:size(slength,1)
    for n = 1:1:size(slength,2)
        % Calculates the coordinates of the arrow in a standard reference frame
        r = [ 0; swidth(m,n); swidth(m,n); hwidth(m,n); 0 ];
        [X,Y] = cylinder(r,N(m,n));
        Z = [ 0; 0; slength(m,n)-hnode(m,n); slength(m,n)-hlength(m,n); slength(m,n) ] * ones(1,N(m,n)+1);
                
        % y-rotation
        Xt = X.*cos(theta(m,n)) + Z.*sin(theta(m,n));
        Zt = -X.*sin(theta(m,n)) + Z.*cos(theta(m,n));
        X = Xt;
        Z = Zt;
        
        % z-rotation
        Xt = X.*cos(phi(m,n)) - Y.*sin(phi(m,n));
        Yt = X.*sin(phi(m,n)) + Y.*cos(phi(m,n));
        X = Xt;
        Y = Yt;
        
        % translation
        X = X + x1(m,n);
        Y = Y + y1(m,n);
        Z = Z + z1(m,n);
        
        % Plots cylinder
        if ishold()
            h = [ h; surf(X,Y,Z,'edgecolor',color,'facecolor',color) ];
        else
            hold on
            h = [ h; surf(X,Y,Z,'edgecolor',color,'facecolor',color) ];
            hold off
        end
    end
end

% Sets other properties
for n = 1:2:length(varargin)
    if ~strcmpi(varargin{n},'color') ...
            & ~strcmpi(varargin{n},'stemwidth') ...
            & ~strcmpi(varargin{n},'headlength') ...
            & ~strcmpi(varargin{n},'headwidth') ...
            & ~strcmpi(varargin{n},'headnode') ...
            & ~strcmpi(varargin{n},'n')
        set(h,varargin{n},varargin{n+1});
    end
end