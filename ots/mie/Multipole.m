classdef Multipole < handle
    % Multipole < handle :  Multipole fields
    %   Multipoles are defined on the coordinates r, theta and phi (or X, Y and Z).
    %   This handle class calculates and stores the values of the multipoles 
    %   in the structures J1_vec, J2_vec, H1_vec and H2_vec.
    %
    % Multipole properties:
    %   farfield    -   exact (default) or farfield
    %   theta       -   polar coordinates [rad]
    %   phi         -   azimuthal coordinates [rad]
    %   r           -   radial coordinate [m]
    %   X           -   x-coordinates [m]
    %   Y           -   y-coordinates [m]
    %   Z           -   z-coordinates [m]
    %   k           -   wave number in the medium [m^-1]
    %   kr          -   k*r
    %   J1_vec      -   J1 multipole values [struct]
    %   J2_vec      -   J2 multipole values [struct]
    %   H1_vec      -   H1 multipole values [struct]
    %   H2_vec      -   H2 multipole values [struct]
    %
    % Multipole methods:
    %   Multipole   -   constructor
    %   J1          -   J1 multipoles
    %   J2          -   J2 multipoles
    %   H1          -   H1 multipoles
    %   H2          -   H2 multipoles
    %   createJ     -   creates J-multipoles with index L
    %   createH     -   creates H-multipoles with index L
    %   plot        -   plots multipoles
    %
    % See also SpBessel, SpHarm, VecSpHarm, ComplexVector.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        farfield    % exact (default) or farfield
        theta       % polar coordinates [rad]
        phi         % azimuthal coordinates [rad]
        r           % radial coordinate [m]
        X           % x-coordinates [m]
        Y           % y-coordinates [m]
        Z           % z-coordinates [m]
        k           % wave number in the medium [m^-1]
        kr          % k*r
        J1_vec      % J1 multipole values (ComplexVector)
        J2_vec      % J2 multipole values (ComplexVector)
        H1_vec      % H1 multipole values (ComplexVector)
        H2_vec      % H2 multipole values (ComplexVector)
    end
    methods
        function obj = Multipole(theta,phi,r,k,varargin)
            % MULIPOLE(THETA,PHI,R,K) constructs the mulitpoles
            %   with polar index L and azimuatal index M
            %   on the coordinates THETA, PHI and R
            %   for light with wave number K.
            %
            % MULIPOLE(L,M,THETA,PHI,R,K,'FarField',true) uses the far-field approximation.
            %
            % See also Multipole, ComplexVector.
            
            % Whether the multipole should be calculated in the far-field [default = false]
            obj.farfield = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'farfield')
                    obj.farfield = varargin{n+1};
                end
            end
            
            obj.theta = theta;
            obj.phi = phi;
            obj.r = r;
            [obj.X,obj.Y,obj.Z] = Transform.Sph2Car(theta,phi,r);
            
            obj.k = k;
            obj.kr = obj.k*obj.r;
            
        end
        function J1lm = J1(multi,l,m)
            % J1 J1 multipole
            %
            % J1lm = J1(MULTI,L,M) returns the multipole J1lm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also Multipole, ComplexVector.
            
            if size(multi.J1_vec,1)>=l+1 && ~isempty(multi.J1_vec{l+1,m+l+1})
                J1lm = multi.J1_vec{l+1,m+l+1};
            else
                multi.createJ(l);
                J1lm = multi.J1(l,m);
            end
        end
        function J2lm = J2(multi,l,m)
            % J2 J2 multipole
            %
            % J2lm = J2(MULTI,L,M) returns the multipole J2lm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also Multipole, ComplexVector.
            
            if size(multi.J2_vec,1)>=l+1 && ~isempty(multi.J2_vec{l+1,m+l+1})
                J2lm = multi.J2_vec{l+1,m+l+1};
            else
                multi.createJ(l);
                J2lm = multi.J2(l,m);
            end
        end
        function H1lm = H1(multi,l,m)
            % H1 H1 multipole
            %
            % H1lm = H1(MULTI,L,M) returns the multipole H1lm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also Multipole, ComplexVector.
            
            if size(multi.H1_vec,1)>=l+1 && ~isempty(multi.H1_vec{l+1,m+l+1})
                H1lm = multi.H1_vec{l+1,m+l+1};
            else
                multi.createH(l);
                H1lm = multi.H1(l,m);
            end
        end
        function H2lm = H2(multi,l,m)
            % H2 H2 multipole
            %
            % H2lm = H2(MULTI,L,M) returns the multipole H2lm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also Multipole, ComplexVector.
            
            if size(multi.H2_vec,1)>=l+1 && ~isempty(multi.H2_vec{l+1,m+l+1})
                H2lm = multi.H2_vec{l+1,m+l+1};
            else
                multi.createH(l);
                H2lm = multi.H2(l,m);
            end
        end
        function createJ(multi,l)
            % CREATEJ Crete J multipoles with index L
            %
            % CREATEJ(MULTI,L) creates the J multipoles with
            %   index L (for all values of M = -L ... 0 ... L).
            %   L	-   polar index (positive integer)
            %
            % See also Multipole, ComplexVector.
                
            vsh = VecSpHarm(multi.theta,multi.phi);
                
            j = SpBessel.j(l,multi.kr);
            dj = SpBessel.dj(l,multi.kr);

            for m = -l:1:l

                multi.J1_vec{l+1,m+l+1} = j.*vsh.Z1(l,m);
                multi.J1_vec{l+1,m+l+1}.X = multi.X;
                multi.J1_vec{l+1,m+l+1}.Y = multi.Y;
                multi.J1_vec{l+1,m+l+1}.Z = multi.Z;

                if multi.farfield
                    multi.J2_vec{l+1,m+l+1} = -dj.*vsh.Z2(l,m);
                    multi.J2_vec{l+1,m+l+1}.X = multi.X;
                    multi.J2_vec{l+1,m+l+1}.Y = multi.Y;
                    multi.J2_vec{l+1,m+l+1}.Z = multi.Z;
                else
                    multi.J2_vec{l+1,m+l+1} = (1i*sqrt(l*(l+1))*j./multi.kr).*vsh.Y(l,m) - (j./multi.kr + dj).*vsh.Z2(l,m);
                    multi.J2_vec{l+1,m+l+1}.X = multi.X;
                    multi.J2_vec{l+1,m+l+1}.Y = multi.Y;
                    multi.J2_vec{l+1,m+l+1}.Z = multi.Z;
                end

            end
            
        end
        function createH(multi,l)
            % CREATEH Crete H multipoles with index L
            %
            % CREATEH(MULTI,L) creates the H multipoles with
            %   index L (for all values of M = -L ... 0 ... L).
            %   L	-   polar index (positive integer)
            %
            % See also Multipole, ComplexVector.
                
            vsh = VecSpHarm(multi.theta,multi.phi);
                
            h = SpBessel.h(l,multi.kr);
            dh = SpBessel.dh(l,multi.kr);

            for m = -l:1:l

                multi.H1_vec{l+1,m+l+1} = h.*vsh.Z1(l,m);
                multi.H1_vec{l+1,m+l+1}.X = multi.X;
                multi.H1_vec{l+1,m+l+1}.Y = multi.Y;
                multi.H1_vec{l+1,m+l+1}.Z = multi.Z;

                if multi.farfield
                    multi.H2_vec{l+1,m+l+1} = -dh.*vsh.Z2(l,m);
                    multi.H2_vec{l+1,m+l+1}.X = multi.X;
                    multi.H2_vec{l+1,m+l+1}.Y = multi.Y;
                    multi.H2_vec{l+1,m+l+1}.Z = multi.Z;
                else
                    multi.H2_vec{l+1,m+l+1} = (1i*sqrt(l*(l+1))*h./multi.kr).*vsh.Y(l,m) - (h./multi.kr + dh).*vsh.Z2(l,m);
                    multi.H2_vec{l+1,m+l+1}.X = multi.X;
                    multi.H2_vec{l+1,m+l+1}.Y = multi.Y;
                    multi.H2_vec{l+1,m+l+1}.Z = multi.Z;
                end
            end
            
        end
        function [fig,h1,h2,h3,h4] = plot(multi,l,m,varargin)
            % PLOT Plots multipoles
            %
            % [FIG,H1,H2,H3,H4] = PLOT(MULTI,L,M) plots multipoles MULTI. 
            %   It returns the number of the figure FIG, 
            %   and the graphic handlers H1, H2, H3, H4 to each of the plotted multipoles.
            %
            % [FIG,H1,H2,H3,H4] = PLOT(MULTI,L,M,'Figure',FIG) used figure FIG.
            %
            % [FIG,H1,H2,H3,H4] = PLOT(MULTI,L,M,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard plot properties
            %   can be used.
            %
            % See also Multipole, ComplexVector.

            % figure number
            fig = 0;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'figure')
                    fig = varargin{n+1};
                end
            end
            
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
            
            if fig==0
                fig = figure;
            else
                figure(fig)
            end
            
            subplot(2,2,1)
            title(['${\bf J}^{(1)}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            h1 = multi.J1(l,m).plot();
            hold off
            axis equal tight
            view(3)
            
            subplot(2,2,2)
            title(['${\bf J}^{(2)}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            h2 = multi.J2(l,m).plot();
            hold off
            axis equal tight
            view(3)
            
            subplot(2,2,3)
            title(['${\bf H}^{(1)}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            h3 = multi.H1(l,m).plot();
            hold off
            axis equal tight
            view(3)
            
            subplot(2,2,4)
            title(['${\bf H}^{(2)}_{' int2str(l) ',' int2str(m) '}$'],'Interpreter','Latex','FontSize',20)
            hold on
            h4 = multi.H2(l,m).plot();
            hold off
            axis equal tight
            view(3)
            
            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'figure') && ~strcmpi(varargin{n},'intensity') && ~strcmpi(varargin{n},'fields')
                    set([h1 h2 h3 h4],varargin{n},varargin{n+1});
                end
            end
            
        end
    end
end