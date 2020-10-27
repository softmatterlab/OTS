classdef Beam
    % Beam (Abstract) : Paraxial light beam
    %   A paraxial light beam is defined by its complex radial Er and
    %   azimuthal Ephi electric fields defined on the polar coordiantes (r,phi), 
    %   by the relative dielectric permittivity er and
    %   the relative magnetic permeability mr of the medium, 
    %   and by the vacuum light wavelength lambda0.
    %   Instances of this class cannot be created. Use one of the subclasses 
    %   (e.g., BeamGauss, BeamHG, BeamLG).
    % 
    % Beam properties:
    %   r       -   (read only) radial coordinate [m]
    %   phi     -   (read only) azimuthal coordinate [m]
    %   Er      -   radial electric field [V/m]
    %   Ephi    -   azimuthal electric field [V/m]
    %   er      -   relative dielectric permittivity
    %   mr      -   relative magnetic permeability
    %   lambda0 -   vacuum wavelength [m]
    %
    % Beam methods:
    %   Beam            -   constructor (accessible only by the subclasses)
    %   plot            -   plots a beam
    %   Br              -   radial magnetic field
    %   Bphi            -   azimuthal magnetic field
    %   uplus           -   +b (= b)
    %   uminus          -   -b (inverts fields)
    %   plus            -   b1+b2 (sums fields)
    %   minus           -   b1-b2 (subtracts fields)
    %   mtimes          -   a*b (fields by a scalar)
    %   intensity       -   intensity [W/m^2]
    %   power           -   power [W]
    %   expand          -   expander
    %   attenuate       -   attenuator
    %   normalize       -   normalization in power
    %   iris            -   iris / aperture stop
    %   focus           -   focal field
    %   focusinterface  -   focal field with interface
    %
    % See also BeamGauss, BeamHG, BeamLG.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties (GetAccess = public, SetAccess = private)
        r       % radial coordinate [m]
        phi     % azimuthal coordinate [rad]
    end
    properties
        Er      % radial electric field [V/m]
        Ephi    % azimuthal electric field [V/m]
        er      % relative dielectric permittivity
        mr      % relative magnetic permeability
        lambda0 % vacuum wavelength [m]
    end
    methods (Access = protected)
        function obj = Beam(R,Nphi,Nr,varargin) 
            % BEAM(R,Nphi,Nr) constructs the coordinates a Beam
            %   The radial coordinate goes up to R with Nr divisions and
            %   the azimuthal coordinate has Nphi divisions.
            %   This method is only accessible by the subclasses of Beam.
            %
            % BEAM(R,Nphi,Nr,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also Beam.

            Check.isreal('R must be a positive real number',R,'>',0)
            Check.isinteger('Nphi must be a positive integer',Nphi,'>',0)
            Check.isinteger('Nr must be a positive integer',Nr,'>',0)
            
            % Vacuum wavelength
            obj.lambda0 = 532e-9;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'lambda0')
                    obj.lambda0 = varargin{n+1};
                    Check.isreal('lambda0 must be a positive real number',obj.lambda0,'>',0)
                end
            end
            
            % Relative dielectric permittivity
            obj.er = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'er')
                    obj.er = varargin{n+1};
                    Check.isnumeric('er must be a number',obj.er)
                end
            end

            % Relative magnetic permeability
            obj.mr = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'mr')
                    obj.mr = varargin{n+1};
                    Check.isnumeric('mr must be a number',obj.mr)
                end
            end
            
            [obj.r,obj.phi] = meshgrid([0.5*R/Nr:R/Nr:R],[2*pi/Nphi:2*pi/Nphi:2*pi]);
        end
    end
    methods (Access = public)
        function plot(obj,varargin) 
            % PLOT Plot beam
            %
            % PLOT(B) plots intensity and polarization of beam B.
            %
            % PLOT(B,'Levels',L) sets the number of levels for the
            %   intensity plot to L (default L = 128).
            %
            % PLOT(B,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue.
            %
            % See also Beam.
            
            
            % Number of levels for the intensity
            L = 128;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'levels')
                    L = varargin{n+1};
                    Check.isinteger('L must be a positive integer',L,'>',0)
                end
            end
            
            % Maximum number of arrows
            V = 25;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'arrows')
                    V = varargin{n+1};
                    Check.isinteger('V must be a positive integer',V,'>',0)
                end
            end
            
            if max(max(abs(obj.r))) > 1e-3
                obj.r = obj.r*1e+3;
                dimension = '[mm]';
            elseif max(max(abs(obj.r))) > 1e-6
                obj.r = obj.r*1e+6;
                dimension = '[\mum]';
            elseif max(max(abs(obj.r))) > 1e-9
                obj.r = obj.r*1e+9;
                dimension = '[nm]';
            else
                dimension = '[m]';
            end
                                    
            hold on
            
            % Cartesian coordinates
            [x,y] = Transform.Pol2Car(obj.phi,obj.r);
            [Ex,Ey] = Transform.Pol2CarVector(obj.phi,obj.Ephi,obj.Er);
            
            % Intensity of the beam in free space
            I = obj.intensity();
            contourf([zeros(size(x,1)+1,1),[x;x(1,:)]],[zeros(size(y,1)+1,1),[y;y(1,:)]],[[I(:,1);I(1,1)],[I;I(1,:)]], L)
            colormap(ones(L,3)-[0:1/(L-1):1]'*[0 1 1])
            shading flat
            
            % Direction in the xy plane of the electric field 
            % (real part in black, imaginary part in blue)
            Etot = sqrt(Ex.*conj(Ex) + Ey.*conj(Ey));
            ind1 = [1:ceil(size(x,1)/V):size(x,1)];
            ind2 = [1:ceil(size(x,2)/V):size(x,2)];
            quiver(x(ind1,ind2), y(ind1,ind2), real(Ex(ind1,ind2))./Etot(ind1,ind2), real(Ey(ind1,ind2))./Etot(ind1,ind2), 0.5, 'k')
            quiver(x(ind1,ind2), y(ind1,ind2), imag(Ex(ind1,ind2))./Etot(ind1,ind2), imag(Ey(ind1,ind2))./Etot(ind1,ind2), 0.5, 'b')
            
            hold off
            
            box on
            axis equal tight
            
            xlabel(cat(2,'x position ',dimension))
            ylabel(cat(2,'y position ',dimension))
            
            % Sets properties
            for n = 1:2:length(varargin)
                if ~strcmpi(varargin{n},'levels') && ~strcmpi(varargin{n},'arrows')
                    set(h,varargin{n},varargin{n+1});
                end
            end
        end
        function B = Br(obj)
            % BR Radial component of magnetic field
            %
            % B = BR(OBJ) Radial component of the magnetic induction field.
            %
            % See also Beam.

            B = -sqrt(obj.er*obj.mr)/PhysConst.c0*obj.Ephi;
        end
        function B = Bphi(obj)
            % BPHI Azimuthal component of magnetic field
            %
            % B = BPHI(OBJ) Azimuthal component of the magnetic induction field.
            %
            % See also Beam.

            B = sqrt(obj.er*obj.mr)/PhysConst.c0*obj.Er;
        end
        function b_p = uplus(b)
            % UPLUS Unitary plus (fields)
            %
            % Bp = UPLUS(B) Unitary plus (Bp = +B).
            %   +B = B
            %
            % See also Beam.
            
            b_p = b;
        end
        function b_m = uminus(b)
            % UMINUS Unitary minus (fields)
            %
            % Bm = UMINUS(B) Unitary minus (Bm = -B).
            %   -B inverts the fields components of B.
            %
            % See also Beam.
            
            b_m = b;
            b_m.Er = -b_m.Er;
            b_m.Ephi = -b_m.Ephi;
        end
        function b = plus(b1,b2)
            % PLUS Binary addition (fields)
            %
            % B = PLUS(B1,B2) Binary addition (B = B1+B2).
            %   The field components Er and Ephi of B are the sum of the
            %   fields components of B1 and B2.
            %
            % See also Beam.
            
            b = b1;
            b.Er = b1.Er+b2.Er;
            b.Ephi = b1.Ephi+b2.Ephi;
        end
        function b = minus(b1,b2)
            % MINUS Binary subtraction (fields)
            %
            % B = MINUS(B1,B2) Binary subtraction (B = B1-B2).
            %   The field components Er and Ephi of B are the difference of the
            %   fields components of B1 and B2.
            %
            % See also Beam.

            b = b1;
            b.Er = b1.Er-b2.Er;
            b.Ephi = b1.Ephi-b2.Ephi;
        end
        function ab = mtimes(a,b)
            % MTIMES Product by a scalar (fields)
            %
            % aB = MTIMES(a,B) Product by a scalar (aB = a*B).
            %   The field components Er and Ephi of aB are the 
            %   fields components of B1 and B2 multiplied by a.
            %
            % See also Beam.
            
            ab = b;
            ab.Er = a*b.Er;
            ab.Ephi = a*b.Ephi;
        end
        function I = intensity(obj)
            % INTENSITY Calculates beam intensity [W/m^2]
            %
            % I = INTENSITY(B) calculates the intensity I of beam B. 
            %   I is a matrix expressed in W/m^2 corresponding to the
            %   coordiantes B.r and B.phi.
            %
            % See also Beam.
            
            I = 0.5/PhysConst.Z0*sqrt(obj.er/obj.mr)*(obj.Er.*conj(obj.Er)+obj.Ephi.*conj(obj.Ephi));
        end       
        function P = power(obj) 
            % POWER Calculates beam power [W]
            %
            % P = POWER(B) calculates the power P of beam B.
            %   P is a scalar expressed in W.
            %
            % See also Beam.

            dr = obj.r(1,2)-obj.r(1,1);
            dphi = obj.phi(2,1)-obj.phi(1,1);
            
            P = sum(sum(obj.intensity().*obj.r*dr*dphi));
        end
        function res = expand(obj,T) 
            % EXPAND Expands beam
            %
            % E = EXPAND(B,T) expands beam B by T.
            %
            % See also Beam.
            
            res = obj;
            res.r = res.r*T;
            res.Er = res.Er/T;
            res.Ephi = res.Ephi/T;
        end
        function res = attenuate(obj,A) 
            % ATTENUATE Attenuates beam
            %
            % R = ATTENUATE(B,A) attenuates the power of beam B by A.
            %
            % See also Beam.
            
            res = obj;
            res.Er = res.Er/sqrt(A);
            res.Ephi = res.Ephi/sqrt(A);
        end
        function res = normalize(obj,P0) 
            % NORMALIZE Normalizes beam power to a value
            %
            % R = NORMALIZE(B,P0) normalizes the power of beam B to P0.
            %
            % See also Beam.

            if nargin<2, P0 = 1; end
            res = obj.attenuate(obj.power()/P0);
        end
        function res = iris(obj,L)
            % IRIS Iris / Aperture stop 
            %
            % R = IRIS(B,L) eliminates the fields outside a radius L.
            %
            % See also Beam.

            res = obj;
            res.Er(res.r>L) = 0;
            res.Ephi(res.r>L) = 0;
        end
        function E = focus(obj,f,X,Y,Z,varargin)
            % FOCUS Focal fields
            % 
            % E = FOCUS(B,F,X,Y,Z) calculates the focal fields of B 
            %   at positions X, Y, Z for a lens with focal lenght F
            % 
            % E = FOCUS(B,F,X,Y,Z,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       er  -   relative electric permittivity [default: B.er]
            %       mr  -   relative magnetic permeability [default: B.mr]
            %
            % See also Beam.

            % Relative dielectric permittivity
            er1 = obj.er; % before objective
            er2 = er1; % after objective
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'er')
                    er2 = varargin{n+1};
                    Check.isnumeric('er must be a number',er2)
                end
            end

            % Relative magnetic permeability
            mr1 = obj.mr; % before objective
            mr2 = mr1; % after objective
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'mr')
                    mr2 = varargin{n+1};
                    Check.isnumeric('mr must be a number',mr2)
                end
            end
            
            % ni = sqrt(er1*mr1);
            % nt = sqrt(er2*mr2);
            % kt = 2*pi*nt/obj.lambda0;
            % 
            % phi = obj.phi;
            % dphi = phi(2,1)-phi(1,1);
            % 
            % theta = asin(obj.r/f);
            % dr = obj.r(1,2)-obj.r(1,1);
            % R = dr*size(obj.r,2);
            % dtheta = ones(size(obj.r,1),1)*(asin([dr:dr:R]/f)-asin([0:dr:R-dr]/f));
            % 
            % utheta = ComplexVector(zeros(size(phi)), zeros(size(phi)), zeros(size(phi)), ...
            %     cos(phi).*cos(theta), ...
            %     sin(phi).*cos(theta), ...
            %     sin(theta));
            % uphi = ComplexVector(zeros(size(phi)),zeros(size(phi)),zeros(size(phi)), ...
            %     -sin(phi), ...
            %     cos(phi), ...
            %     zeros(size(phi)));
            % 
            % Et = (obj.Ephi*uphi + obj.Er*utheta)*(sqrt(ni/nt)*sqrt(cos(theta)));
            % ktx = -kt*sin(theta).*cos(phi);
            % kty = -kt*sin(theta).*sin(phi);
            % ktz = kt*cos(theta);
            % 
            % E = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));
            % for n = 1:1:numel(Et)
            %     phase = exp(1i*ktx(n)*X+1i*kty(n)*Y).*exp(1i*ktz(n)*Z);
            %     E = E + ComplexVector(E.X,E.Y,E.Z, ...
            %         Et.Vx(n)*phase*sin(theta(n))*dphi*dtheta(n), ...
            %         Et.Vy(n)*phase*sin(theta(n))*dphi*dtheta(n), ...
            %         Et.Vz(n)*phase*sin(theta(n))*dphi*dtheta(n) ...
            %         );
            % end
            % E = 1i*kt*f*exp(-1i*kt*f)/(2*pi)*E;

            % Numerically more stable version
            phi = obj.phi;
            dphi = phi(2,1)-phi(1,1);
            
            theta = asin(obj.r/f);
            dr = obj.r(1,2)-obj.r(1,1);
            R = dr*size(obj.r,2);
            dtheta = ones(size(obj.r,1),1)*(asin([dr:dr:R]/f)-asin([0:dr:R-dr]/f));
            
            utheta = Point(cos(phi).*cos(theta), sin(phi).*cos(theta), -sin(theta)).tovector();
            uphi = Point(-sin(phi), cos(phi), zeros(size(phi))).tovector();
            
            k = 2*pi*sqrt(er2*mr2)/obj.lambda0;
            rho = sqrt(X.^2+Y.^2);
            varphi = atan2(Y,X);
            
            n1 = sqrt(er1);
            n2 = sqrt(er2);
            Et = (obj.Ephi*uphi + obj.Er*utheta)*(sqrt(n1/n2)*sqrt(cos(theta)));
            
            E = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));
            for n = 1:1:numel(Et)
                phase = exp(1i*(k*cos(theta(n))*Z+k*sin(theta(n))*rho.*cos(varphi-phi(n))));
                E = E + ComplexVector(E.X,E.Y,E.Z, ...
                    Et.Vx(n)*phase*sin(theta(n))*dphi*dtheta(n), ...
                    Et.Vy(n)*phase*sin(theta(n))*dphi*dtheta(n), ...
                    Et.Vz(n)*phase*sin(theta(n))*dphi*dtheta(n) ...
                    );
            end
            E = 1i*k*f*exp(-1i*k*f)/(2*pi)*E;

        end
        function E = focusinterface(obj,f,X,Y,Z,z0,varargin)
            % FOCUSINTERFACE Focal fields
            % 
            % E = FOCUSINTERFACE(B,F,X,Y,Z,Z0) calcualtes the focal fields of B 
            %   at positions X, Y, Z for a lens with focal lenght F with an
            %   interface at Z0.
            % 
            % E = FOCUSINTERFACE(B,F,X,Y,Z,Z0,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       er  -   relative electric permittivity (before interface) [default: B.er]
            %       mr	-   relative magnetic permeability (before interface) [default: B.mr]
            %       es  -   relative electric permittivity (after interface) [default: B.er]
            %       ms  -   relative magnetic permeability (after interface) [default: B.mr]
            %
            % See also Beam.

            % Relative dielectric permittivity
            er1 = obj.er; % before objective
            er2 = er1; % after objective
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'er')
                    er2 = varargin{n+1};
                    Check.isnumeric('er must be a number',er2)
                end
            end
            es = er1; % in the sample
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'es')
                    es = varargin{n+1};
                    Check.isnumeric('es must be a number',es)
                end
            end

            % Relative magnetic permeability
            mr1 = obj.mr; % before objective
            mr2 = mr1; % after objective
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'mr')
                    mr2 = varargin{n+1};
                    Check.isnumeric('mr must be a number',mr2)
                end
            end
            ms = mr1; % in the sample
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'ms')
                    ms = varargin{n+1};
                    Check.isnumeric('ms must be a number',ms)
                end
            end

            ni = sqrt(er1*mr1);
            nt = sqrt(er2*mr2);
            kt = 2*pi*nt/obj.lambda0;
            ns = sqrt(es*ms);
            ks = 2*pi*ns/obj.lambda0;
            
            phi = obj.phi;
            dphi = phi(2,1)-phi(1,1);
            
            theta = asin(obj.r/f);
            dr = obj.r(1,2)-obj.r(1,1);
            R = dr*size(obj.r,2);
            dtheta = ones(size(obj.r,1),1)*(asin([dr:dr:R]/f)-asin([0:dr:R-dr]/f));
            
            thetas = asin(nt/ns*sin(theta));

            % utheta = ComplexVector(zeros(size(phi)), zeros(size(phi)), zeros(size(phi)), ...
            %     cos(phi).*cos(theta), ...
            %     sin(phi).*cos(theta), ...
            %     sin(theta));
            us = ComplexVector(zeros(size(phi)),zeros(size(phi)),zeros(size(phi)), ...
                cos(phi).*real(cos(thetas)), ...
                sin(phi).*real(cos(thetas)), ...
                min(ones(size(thetas)),real(sin(thetas))));
            uphi = ComplexVector(zeros(size(phi)),zeros(size(phi)),zeros(size(phi)), ...
                -sin(phi), ...
                cos(phi), ...
                zeros(size(phi)));

            % Et = (obj.Ephi*uphi + obj.Er*ur)*(sqrt(ni/nt)*sqrt(cos(theta)));
            % ktx = -kt*sin(theta).*cos(phi);
            % kty = -kt*sin(theta).*sin(phi);
            ktz = kt*cos(theta);
            
            ts = 2*nt*cos(theta)./(nt*cos(theta)+ns*cos(thetas));
            tp = 2*nt*cos(theta)./(nt*cos(thetas)+ns*cos(theta));
            Es = (ts.*obj.Ephi*uphi + tp.*obj.Er*us)*(sqrt(ni/nt)*sqrt(cos(theta)));
            ksx = -kt*sin(theta).*cos(phi); % -ks*sin(thetas).*cos(phi); 
            ksy = -kt*sin(theta).*sin(phi); % -ks*sin(thetas).*sin(phi);
            ksz = kt*sqrt(ns^2/nt^2-sin(theta).^2); % ks*cos(thetas);
            
            E = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));
            for n = 1:1:numel(Es)
                phase = exp(1i*(ksx(n)*X+ksy(n)*Y+ktz(n)*z0+ksz(n)*(Z-z0)));
                E = E + ComplexVector(E.X,E.Y,E.Z, ...
                    Es.Vx(n)*phase*sin(theta(n))*dphi*dtheta(n), ...
                    Es.Vy(n)*phase*sin(theta(n))*dphi*dtheta(n), ...
                    Es.Vz(n)*phase*sin(theta(n))*dphi*dtheta(n) ...
                    );
            end
            E = 1i*kt*f*exp(-1i*kt*f)/(2*pi)*E;
        end        
    end
end