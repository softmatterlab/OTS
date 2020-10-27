classdef VecSpHarm < handle
    % VecSpHarm < handle : Vector spherical harmonics
    %   Vector spherical harmonics are defined on the coordinates theta and phi.
    %   This handle class calculates and stores the values of the vector spherical
    %   harmonics in the structures Y_vec, Z1_vec, Z2_vec.
    %
    % VecSpHarm properties:
    %   theta   -   polar coordinates [rad]
    %   phi     -   azimuthal coordinates [rad]
    %   ur      -   radial unit vector (ComplexVector)
    %   utheta  -   polar unit vector (ComplexVector)
    %   uphi    -   azimuthal unit vector (ComplexVector)
    %   Y_vec	-   Y spherical harmonics values [struct]
    %   Z1_vec  -   Z1 spherical harmonics values [struct]
    %   Z2_vec  -   Z2 spherical harmonics values [struct]
    %
    % VecSpHarm methods:
    %   VecSpHarm   -   constructor 
    %   Y           -   Y-vector spherical harmonics
    %   Z1          -   Z1-vector spherical harmonics
    %   Z2          -   Z2-vector spherical harmonics
    %   create      -   creates vector spehrical harmonics with index L
    %   plot        -   plots vector spherical harmonics
    %
    % See also SpBessel, SpHarm, Multipole, ComplexVector.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        theta   % polar coordinates [rad]
        phi     % azimuthal coordinates [rad]
        ur      % radial unit vector (ComplexVector)
        utheta  % polar unit vector (ComplexVector)
        uphi    % azimuthal unit vector (ComplexVector)
        Y_vec	% Y spherical harmonics values (ComplexVector)
        Z1_vec	% Z1 spherical harmonics values (ComplexVector)
        Z2_vec	% Z2 spherical harmonics values (ComplexVector)
    end
    methods
        function obj = VecSpHarm(theta,phi)
            % VECSPHARM(THETA,PHI) constructs a vector spherical harmonics 
            %   on the coordinates THETA and PHI.
            %
            % See also VecSpHarm, ComplexVector.

            obj.theta = theta;
            obj.phi = phi;

            [x,y,z] = Transform.Sph2Car(theta,phi,ones(size(theta)));
            obj.ur = ComplexVector(x,y,z,cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta));
            obj.utheta = ComplexVector(x,y,z,cos(phi).*cos(theta),sin(phi).*cos(theta),-sin(theta));
            obj.uphi = ComplexVector(x,y,z,-sin(phi),cos(phi),zeros(size(theta)));
        end
        function Ylm = Y(vsh,l,m)
            % Y Y-vector spherical harmonics
            %
            % Ylm = Y(VSH,L,M) returns the Y-vector spherical harmonics Ylm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also VecSpHarm, ComplexVector.
            
            if size(vsh.Y_vec,1)>=l+1 && ~isempty(vsh.Y_vec{l+1,m+l+1})
                Ylm = vsh.Y_vec{l+1,m+l+1};
            else
                vsh.create(l)
                Ylm = vsh.Y(l,m);
                
            end            
        end
        function Z1lm = Z1(vsh,l,m)
            % Z1 Z1-vector spherical harmonics
            %
            % Z1lm = Z1(VSH,L,M) returns the Z1-vector spherical harmonics Ylm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also VecSpHarm, ComplexVector.
            
            if size(vsh.Z1_vec,1)>=l+1 && ~isempty(vsh.Z1_vec{l+1,m+l+1})
                Z1lm = vsh.Z1_vec{l+1,m+l+1};
            else
                vsh.create(l)
                Z1lm = vsh.Z1(l,m);
            end
        end
        function Z2lm = Z2(vsh,l,m)
            % Z2 Z2-vector spherical harmonics
            %
            % Z2lm = Z2(VSH,L,M) returns the Z2-vector spherical harmonics Ylm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also VecSpHarm, ComplexVector.
            
            if size(vsh.Z2_vec,1)>=l+1 && ~isempty(vsh.Z2_vec{l+1,m+l+1})
                Z2lm = vsh.Z2_vec{l+1,m+l+1};
            else
                vsh.create(l);
                Z2lm = vsh.Z2(l,m);
            end
        end
        function create(vsh,l)
            % CREATE Crete vector spherical harmonics with index L
            %
            % CREATE(VSH,L) creates the vector spherical harmonics with
            %   index L (for all values of M = -L ... 0 ... L).
            %   L	-   polar index (positive integer)
            %
            % See also VecSpHarm, ComplexVector.
            
            dtheta = pi/1e+6;
            theta = vsh.theta;
            theta(vsh.theta==0) = 3*dtheta;

            sh = SpHarm(theta,vsh.phi);
            
            for m = -l:1:l
                
                vsh.Y_vec{l+1,m+l+1} = sh.Y(l,m)*vsh.ur;
                
                if l>0
                    DY = sh.dYtheta(l,m).*vsh.utheta + sin(theta).^-1.*sh.dYphi(l,m).*vsh.uphi;
                    vsh.Z1_vec{l+1,m+l+1} = -(1i/sqrt(l*(l+1)))*vsh.ur*DY;
                    vsh.Z2_vec{l+1,m+l+1} = vsh.Z1(l,m)*vsh.ur;
                else
                    vsh.Z1_vec{l+1,m+l+1} = 0*vsh.Y(l,m);
                    vsh.Z2_vec{l+1,m+l+1} = 0*vsh.Y(l,m);
                end
                
            end
            
        end
        function plot(vsh,l,m,varargin)
            % PLOT Plots vector spherical harmonics
            %
            % PLOT(VSH,L,M) plots the vector spherical harmonics VSH. 
            %   It returns a graphic handler to the plotted object.
            %
            % PLOT(SH,L,M,'Intensity',INTENSITY) selects whether the
            %   intensity plot is on (INTENSITY='on' - default) 
            %   or off (INTENSITY='off').
            %
            % PLOT(SH,L,M,'Fields',FIELDS) selects whether the
            %   field plot is on (FIELDS='on' - default) 
            %   or off (FIELDS='off').
            %
            % See also VecSpHarm, surf.
            
            % Intensity: on(default)|off
            intensity = 'on';
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'intensity')
                    intensity = varargin{n+1};
                end
            end
            
            % fields: on(default)|off
            fields = 'on';
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'fields')
                    fields = varargin{n+1};
                end
            end
            
            figure

            subplot(2,2,1)
            title(['$Y_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            [x,y,z] = Transform.Sph2Car(vsh.theta,vsh.phi,ones(size(vsh.theta)));
            surf(x,y,z,vsh.Y(l,m).norm());
            hold off
            axis equal tight
            view(3)
            
            subplot(2,2,2)
            title(['${\bf Y}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            if strcmpi(intensity,'on')
                [x,y,z] = Transform.Sph2Car(vsh.theta,vsh.phi,ones(size(vsh.theta)));
                surf(x,y,z,vsh.Y(l,m).norm());
            end
            if strcmpi(fields,'on')
                vsh.Y(l,m).plot();
            end
            hold off
            axis equal tight
            view(3)
            
            subplot(2,2,3)
            title(['${\bf Z}^{(1)}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            if strcmpi(intensity,'on')
                [x,y,z] = Transform.Sph2Car(vsh.theta,vsh.phi,ones(size(vsh.theta)));
                surf(x,y,z,vsh.Z1(l,m).norm());
            end
            if strcmpi(fields,'on')
                vsh.Z1(l,m).plot();
            end
            hold off
            axis equal tight
            view(3)
            
            subplot(2,2,4)
            title(['${\bf Z}^{(2)}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            if strcmpi(intensity,'on')
                [x,y,z] = Transform.Sph2Car(vsh.theta,vsh.phi,ones(size(vsh.theta)));
                surf(x,y,z,vsh.Z2(l,m).norm());
            end
            if strcmpi(fields,'on')
                vsh.Z2(l,m).plot();
            end
            hold off
            axis equal tight
            view(3)
        end
    end
end