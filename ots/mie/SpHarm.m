classdef SpHarm < handle
    % SpHarm < handle : Spherical harmonics
    %   Spherical harmonics are defined on the coordinates theta and phi.
    %   This handle class calculates and stores the values of the spherical
    %   harmonics and their derivatives in the structures Y_vec,
    %   dYtheta_vec and dYphi_vec.
    %
    % SpHarm properties:
    %   theta           -   polar coordinates [rad]
    %   phi             -   azimuthal coordinates [rad]
    %   Y_vec           -   spherical harmonics [struct]
    %   dYtheta_vec     -   theta-derivative of spherical harmonics [struct]
    %   dYphi_vec       -   phi-derivative of spherical harmonics [struct]
    %
    % SpHarm methods:
    %   SpHarm  -   constructor
    %   Y       -   spherical harmonics
    %   dYtheta -   theta-derivative of spherical harmonics
    %   dYphi   -   phi-derivative of spherical harmonics
    %   plot    -   plots spherical harmonics
    %
    % See also SpBessel, VecSpHarm, Multipole.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        theta       % polar coordiantes [rad]
        phi         % azimuthal coordinates [rad]
        Y_vec       % spherical harmonics [struct]
        dYtheta_vec % theta-derivative of spherical harmonics [struct]
        dYphi_vec   % phi-derivative of spherical harmonics [struct]
    end
    methods
        function obj = SpHarm(theta,phi)
            % SPHARM(THETA,PHI) constructs a spherical harmonics 
            %   on the coordinates THETA and PHI.
            %
            % See also SpHarm.

            obj.theta = theta;
            obj.phi = phi;
        end
        function Ylm = Y(sh,l,m)
            % Y Spherical harmonics
            %
            % Ylm = Y(SH,L,M) returns the spherical harmonics Ylm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also SpHarm.
                        
            if size(sh.Y_vec,1)>=l+1 && ~isempty(sh.Y_vec{l+1,m+l+1})

                Ylm = sh.Y_vec{l+1,m+l+1};

            else
                
                Pl = sqrt(0.5/pi)*legendre(l,cos( ...
                    reshape(sh.theta,1,numel(sh.theta)) ...
                    ),'norm');
                
                for mi = -l:1:l
                    Plm = (-1)^mi*reshape(Pl(abs(mi)+1,:),size(sh.theta));
                    if mi>=0
                        sh.Y_vec{l+1,mi+l+1} = exp(1i*mi*sh.phi).*Plm;
                    else
                        sh.Y_vec{l+1,mi+l+1} = (-1)^mi*exp(1i*mi*sh.phi).*conj(Plm);
                    end
                end
                
                Ylm = sh.Y(l,m);
                
            end
        end
        function dYthetalm = dYtheta(sh,l,m)
            % DYTHETA theta-derivative of spherical harmonics
            %
            % dYthetalm = DYTHETA(SH,L,M) returns the theta-derivative of spherical harmonics 
            %   dYthetalm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also SpHarm.
                        
            if size(sh.dYtheta_vec,1)>=l+1 && ~isempty(sh.dYtheta_vec{l+1,m+l+1})

                dYthetalm = sh.dYtheta_vec{l+1,m+l+1};

            else
                
                dtheta = pi/1e+6;
                theta = sh.theta;
                theta(sh.theta==0) = 3*dtheta;

                SH2theta = SpHarm(theta+dtheta,sh.phi);
                SH1theta = SpHarm(theta-dtheta,sh.phi);
                
                for mi = -l:1:l
                    sh.dYtheta_vec{l+1,mi+l+1} = (SH2theta.Y(l,mi)-SH1theta.Y(l,mi))/(2*dtheta);
                end
                
                dYthetalm = sh.dYtheta(l,m);
                
            end
        end
        function dYphilm = dYphi(sh,l,m)
            % DYPHI phi-derivative of spherical harmonics
            %
            % dYphilm = DYPHI(SH,L,M) returns the phi-derivative of spherical harmonics 
            %   dYthetalm with indices L and M.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % See also SpHarm.
                        
            if size(sh.dYphi_vec,1)>=l+1 && ~isempty(sh.dYphi_vec{l+1,m+l+1})

                dYphilm = sh.dYphi_vec{l+1,m+l+1};

            else
                
                dtheta = pi/1e+6;
                theta = sh.theta;
                theta(sh.theta==0) = 3*dtheta;
                
                dphi = pi/1e+6;
                SH2phi = SpHarm(theta,sh.phi+dphi);
                SH1phi = SpHarm(theta,sh.phi-dphi);
                
                for mi = -l:1:l
                    sh.dYphi_vec{l+1,mi+l+1} = (SH2phi.Y(l,mi)-SH1phi.Y(l,mi))/(2*dphi);
                end
                
                dYphilm = sh.dYphi(l,m);
                
            end
        end
        function h = plot(sh,l,m,varargin)
            % PLOT Plots spehrical harmonics
            %
            % H = PLOT(SH,L,M) plots the spherical harmonics SH with indices L and M. 
            %   It returns a graphic handler to the plotted object.
            %   L	-   polar index (positive integer)
            %   M	-   azimuthal index (M in -l, ..., 0, ..., +l)
            %
            % H = PLOT(SH,L,M,'Kind',KIND) selects the kind of plot.
            %   KIND can be one of the following:
            %       'unit'  -   (default) absolute value of the spherical harmonics on the unit sphere
            %       'abs'   -   absolute value of the spherical harmonics on a scattergram
            %       'real'  -   real part of the spherical harmonics on a scattergram
            %       'imag'  -   imaginary part of the spherical harmonics on a scattergram
            %
            % H = PLOT(SH,L,M,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. All standard surf properties
            %   can be used.
            %
            % See also SpHarm, surf.

            % Kind of plot: unit (default)|abs|real|imag
            kind = 'unit';
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'kind')
                    kind = varargin{n+1};
                end
            end
                        
            if strcmpi(kind,'abs')
                [x,y,z] = Transform.Sph2Car(sh.theta,sh.phi,abs(sh.Y(l,m)));
                h = surf(x,y,z,abs(sh.Y(l,m)));
            elseif strcmpi(kind,'real')
                [x,y,z] = Transform.Sph2Car(sh.theta,sh.phi,abs(real(sh.Y(l,m))));
                h = surf(x,y,z,real(sh.Y(l,m)));
            elseif strcmpi(kind,'imag')
                [x,y,z] = Transform.Sph2Car(sh.theta,sh.phi,abs(imag(sh.Y(l,m))));
                h = surf(x,y,z,imag(sh.Y(l,m)));
            else
                [x,y,z] = Transform.Sph2Car(sh.theta,sh.phi,ones(size(sh.theta)));
                h = surf(x,y,z,abs(sh.Y(l,m)));
            end
            axis equal
            view(3)

            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'kind')
                    set(h,varargin{n},varargin{n+1});
                end
            end
        end
    end
end