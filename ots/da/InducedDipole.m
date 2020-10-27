classdef InducedDipole < EField
    % InducedDipole < EField : Induced dipole
    %
    % InducedDipole properties:
    %   lambda0 -   vacuum wavelength [m] < EField
    %   er      -   relative dielectric permittivity < EField
    %   mr      -   relative magnetic permeability < EField
    %   alpha   -   polarizability [Cm^2/V]
    %   rd      -   position (Point) [m]
    %
    % InducedDipole methods:
    %   InducedDipole   -   constructor 
    %   n               -   refractive index < EField
    %   lambda          -   wavelenght in the medium [m]  < EField
    %   k               -   wave number in the medium [m^-1]  < EField
    %   omega           -   angular frequency [Hz]  < EField
    %   B               -   magnetic field [T] < EField
    %   S               -   Poynting vector (Vector) [W/m^2] < EField
    %   Ls            	-   spin density [kg m^2/s] < EField
    %   E               -   electric field [V/m]
    %   dipolemoment    -   dipole moment
    %   Estandard       -   electric field for a dipole with p = 1 at the origin and orineted along z [V/m]
    %   sext            -   extinction cross-section [m^-2]
    %   sscat           -   scattering cross-section [m^-2]
    %   sabs            -   absorption cross-section [m^-2]
    %   force           -   force on the dipole in a EM field [N]
    %
    % InducedDipole static method:
    %   polarizability  -   polarizability
    %
    % See also EField.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        alpha   % polarizability [Cm^2/V]
        rd      % position (Point) [m]
    end
    methods
        function id = InducedDipole(alpha,varargin)
            % INDUCEDDIPOLE(ALPHA) constructs an induced dipole of polarizability ALPHA.
            %
            % EFIELDFOCUS(ALPHA,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %       rd          -	Position (Point) [default: Point(0,0,0)]
            %
            % See also InducedDipole, EField.

            id = id@EField(varargin{:});            

            Check.isnumeric('The polarizability must be a number',alpha)

            id.alpha = alpha;

            % Position
            id.rd = Point(0,0,0);
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'rd')
                    id.rd = varargin{n+1};  % V 1.0.2
                    Check.isa('The position of a dipole must be a Point',id.rd,'Point')
                end
            end
            
        end
        function E = E(id,R,Ei)
            % E Electric field [V/m]
            %
            % E = E(ID,R,Ei) calculates the electric field at positions R (Point)
            %   for the induced dipole ID usbject to the electric field Ei.
            %   E is a ComplexVector.
            %
            % See also InducedDipole, InducedDipole.Estandard, Point, ComplexVector.
                        
            Check.isa('The set of positions where to calculate the electric field must be a Point',R,'Point')
            Check.isa('The inducing electric field must be a ComplexVector',Ei,'ComplexVector')

            p = id.dipolemoment(Ei);

            R = R-id.rd;
            E = p.Vx * id.Estandard(R.yrotation(-pi/2)).yrotation(pi/2) ...
                + p.Vy * id.Estandard(R.xrotation(pi/2)).xrotation(-pi/2) ...
                + p.Vz * id.Estandard(R);
            E.X = E.X+id.rd.X;
            E.Y = E.Y+id.rd.Y;
            E.Z = E.Z+id.rd.Z;
        end
        function p = dipolemoment(id,Ei)
            % DIPOLEMOMENT Electric field
            %
            % P = DIPOLEMOMENT(ID,Ei) calculates the dipole moment of the
            %   induced dipole ID usbject to the electric field Ei.
            %   P is a ComplexVector.
            %
            % See also InducedDipole, ComplexVector.

            Check.isa('The inducing electric field must be a ComplexVector',Ei,'ComplexVector')

            p = id.alpha*Ei;
        end
        function E = Estandard(id,R)
            % ESTANDARD Electric field for a dipole with p = 1 at the origin and orineted along z [V/m]
            %
            % E = ESTANDARD(ID,R) calculates the electric field at positions R (Point)
            %   for the induced dipole with p = 1 at the origin and orineted along z.
            %   E is a ComplexVector.
            %
            % See also InducedDipole, Point, ComplexVector.
            
            ONES = ones(size(R));
            ZEROS = zeros(size(R));

            [theta,phi,r] = Transform.Car2Sph(R.X,R.Y,R.Z);
            kr = id.k()*r;
            
            [Vx,Vy,Vz] = Transform.Sph2CarVector(theta,phi,ZEROS,ZEROS,ONES);
            ur = ComplexVector(R.X,R.Y,R.Z,Vx,Vy,Vz);

            [Vx,Vy,Vz] = Transform.Sph2CarVector(theta,phi,ONES,ZEROS,ZEROS);
            utheta = ComplexVector(R.X,R.Y,R.Z,Vx,Vy,Vz);
            
            E = id.k()^3/(4*pi*PhysConst.e0*id.er) * exp(1i*kr)./kr ...
                .* ( ...
                2*cos(theta).*(kr.^-2-1i*kr.^-1) * ur ...
                + sin(theta).*(kr.^-2-1i*kr.^-1-1) * utheta ...
                );
        end
        function s = sext(id)
            % SEXT Extinction cross-section [m^-2]
            %
            % S = SEXT(ID) calculates the extinction cross-section of ID.
            %
            % See also InducedDipole.

            s = id.k()/(PhysConst.e0*id.er)*imag(id.alpha);
        end
        function s = sscat(id)
            % SSCAT Scattering cross-section [m^-2]
            %
            % S = SSCAT(ID) calculates the extinction cross-section of ID.
            %
            % See also InducedDipole.

            s = id.k()^4/(6*pi*(PhysConst.e0*id.er)^2)*abs(id.alpha)^2;
        end        
        function s = sabs(id)
            % SABS Absorption cross-section [m^-2]
            %
            % S = SABS(ID) calculates the absorption cross-section of ID.
            %
            % See also InducedDipole.

            s = id.sext()-id.sscat();
        end
        function [F,Fgrad,Fscat,Fsc] = force(id,r,ef,varargin)
            % FORCE Force on the dipole in a EM field [N]
            %
            % [F,Fgrad,Fscat,Fsc] = FORCE(ID,R,EF) calculates the force exerted 
            %   on ID by the electric field EF (EField) at positions R (Point).
            %
            % [F,Fgrad,Fscat,Fsc] = FORCE(ID,R,EF,'dr',DR) sets the increment 
            %   in the calcualtion of the derivatives to DR [default = 1e-10 m].
            %
            % See also InducedDipole, EField, Point, Vector.
            
            % increment [m]
            dr = 1e-10;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'dr')
                    dr = varargin{n+1};
                    Check.isnumeric('dr must be a positive real number',dr,'>',0)
                end
            end
            
            Fgrad = .25*real(id.alpha)*Vector(r.X,r.Y,r.Z, ...
                ( norm(ef.E(r+Point(dr,0,0))).^2 - norm(ef.E(r+Point(-dr,0,0))).^2 )/(2*dr), ...
                ( norm(ef.E(r+Point(0,dr,0))).^2 - norm(ef.E(r+Point(0,-dr,0))).^2 )/(2*dr), ...
                ( norm(ef.E(r+Point(0,0,dr))).^2 - norm(ef.E(r+Point(0,0,-dr))).^2 )/(2*dr) ...
                );
            
            Fscat = id.sext()*id.n()/PhysConst.c0*ef.S(r,varargin{:});
            
            [Ls,DxLs] = ef.Ls(r,varargin{:});
            
            Fsc = id.sext()*PhysConst.c0/id.n()*DxLs;
            
            F = Fgrad + Fscat + Fsc;
            
        end
    end
    methods (Static)
        function alpha = polarizability(kind,a,ep,varargin)
            % POLARIZABILITY Polarizability
            %
            % ALPHA = POLARIZABILITY(KIND,A,EP) calcualtes the
            %   polarizability of spherical particle of radius A and
            %   relative refractive index EP.
            %   KIND = 'Corrected' uses the formula with radiative correction
            %   KIND = 'Clausius-Mossotti' uses the Clausius-Mossotti formula, 
            %
            % ALPHA = POLARIZABILITY(KIND,A,EP,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       em          -   relative electric permittivity [default: 1]
            %   
            % See also InducedDipole
            
            Check.isreal('The radius must be a positive real number',a,'>',0)
            Check.isnumeric('The particle polarisability must be a real number',ep)

            % medium relative dielectric constant
            em = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'em')
                    em = varargin{n+1};
                    Check.isreal('em must be a real number greater than or equal to 1',em,'>=',1)
                end
            end
            
            er = ep/em;
            
            switch lower(kind)
                
                case 'clausius-mossotti'
                    V = 4/3*pi*a.^3;
                    alpha = 3*V*PhysConst.e0*em * (er-1)/(er+2);
                    
                otherwise % polarizability with radiative correction

                    % vacuum wavelength [m]
                    lambda0 = 532e-9;
                    for n = 1:2:length(varargin)
                        if strcmpi(varargin{n},'lambda0')
                            lambda0 = varargin{n+1};
                            Check.isreal('lambda0 must be a positive real number',lambda0,'>',0)
                        end
                    end
                    
                    lambda = lambda0/sqrt(em);
                    ka = 2*pi/lambda*a;
            
                    alpha0 = InducedDipole.polarizability('Clausius-Mossotti',a,ep,varargin{:});
                    alpha = alpha0.*(1-(er-1)/(er+2)*( ka.^2 + 2*1i/3*ka.^3) ).^-1;
                    
            end
        end
    end
end