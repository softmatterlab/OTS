classdef EField
    % EField (Abstract) : Electromagnetic field
    %   An electromagnetic field associates to each point in space an
    %   electric an a magnetic field.
    %   Instances of this class cannot be created. Use one of the subclasses 
    %   (e.g., EFieldPlaneWave,EFieldFocus, InducedDipole).
    % 
    % EField properties:
    %   lambda0 -   vacuum wavelength [m]
    %   er      -   relative dielectric permittivity
    %   mr      -   relative magnetic permeability
    %
    % EField methods:
    %   EField	-   constructor (accessible only by the subclasses)
    %   n       -   refractive index
    %   lambda  -   wavelenght in the medium [m]
    %   k       -   wave number in the medium [m^-1]
    %   omega   -   angular frequency [Hz]
    %   B       -   magnetic field [T]
    %   S       -   Poynting vector (Vector) [W/m^2]
    %   Ls      -   spin density [kg m^2/s]
    %
    % EField methods (abstract):
    %   E       -   electric field [V/m]
    %
    % See also EFieldPlaneWave,EFieldFocus, InducedDipole.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        lambda0 % vacuum wavelength [m]
        er      % relative dielectric permittivity
        mr      % relative magnetic permeability
    end
    methods (Access = protected)
        function ef = EField(varargin)
            % EFIELD() constructs an EField
            %   This method is only accessible by the subclasses of Beam.
            %
            % EFIELD('PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also EField.
            
            % Vacuum wavelength
            ef.lambda0 = 532e-9;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'lambda0')
                    ef.lambda0 = varargin{n+1};
                    Check.isreal('lambda0 must be a positive real number',ef.lambda0,'>',0)
                end
            end
            
            % Relative dielectric permittivity
            ef.er = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'er')
                    ef.er = varargin{n+1};
                    Check.isnumeric('er must be a number',ef.er)
                end
            end

            % Relative magnetic permeability
            ef.mr = 1;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'mr')
                    ef.mr = varargin{n+1};
                    Check.isnumeric('mr must be a number',ef.mr)
                end
            end            
        end
    end
    methods (Abstract) 
        E(ef,r,varargin)  % Electromagnetic field at position r
    end
    methods
        function n = n(ef)
            % N Refractive index
            %
            % N = N(EF) calcualtes the refractive index.
            %
            % See also EField.

            n = sqrt(ef.er*ef.mr);
        end
        function lambda = lambda(ef)
            % LAMBDA Wavelength [m]
            %
            % LAMBDA = LAMBDA(EF) calcualtes the wavelength in the medium.
            %
            % See also EField.

            lambda = ef.lambda0/ef.n();
        end
        function k = k(ef)
            % K Wave number [m^-1]
            %
            % K = K(EF) calcualtes the Wave number in the medium.
            %
            % See also EField.
            
            k = 2*pi/ef.lambda();
        end
        function omega = omega(ef)
            % OMEGA Angular frequency [Hz]
            %
            % OMEGA = OMEGA(EF) calcualtes the angular frequency in the medium.
            %
            % See also EField.
            
            omega = 2*pi*PhysConst.c0/ef.lambda0;
        end
        function B = B(ef,r,varargin)
            % B Magnetic field [T]
            %
            % B = B(EF,R) calculates the magnetic field at positions R (Point).
            %   B is a ComplexVector.
            %
            % B = B(EF,R,'dr',DR) sets the increment in the calcualtion of
            %   the derivatives to DR [default = 1e-12 m]. 
            %
            % See also EField, Point, ComplexVector.
            
            Check.isa('The set of positions where to calculate the magnetic field must be a Point',r,'Point')

            % increment [m]
            dr = 1e-12;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'dr')
                    dr = varargin{n+1};
                    Check.isnumeric('dr must be a positive real number',dr,'>',0)
                end
            end
            
            dxE = ( ef.E(r+Point(dr,0,0),varargin{:}) - ef.E(r+Point(-dr,0,0),varargin{:}) )./(2*dr);
            dyE = ( ef.E(r+Point(0,dr,0),varargin{:}) - ef.E(r+Point(0,-dr,0),varargin{:}) )./(2*dr);
            dzE = ( ef.E(r+Point(0,0,dr),varargin{:}) - ef.E(r+Point(0,0,-dr),varargin{:}) )./(2*dr);
            
            B = -1i/ef.omega() * ComplexVector(r.X,r.Y,r.Z, ...
                dyE.Vz-dzE.Vy, ...
                -dxE.Vz+dzE.Vx, ...
                dxE.Vy-dyE.Vx ...
                );
        end
        function S = S(ef,r,varargin)
            % S Poynting vector [W/m^2]
            %
            % S = S(EF,R) calculates the Poynting vector at positions R (Point).
            %   S is a ComplexVector.
            %
            % S = S(EF,R,'dr',DR) sets the increment in the calcualtion of
            %   the derivatives to DR [default = 1e-10 m]. 
            %
            % See also EField, Point, ComplexVector.

            S = .5/(PhysConst.m0*ef.mr)*real(ef.E(r,varargin{:})*conj(ef.B(r,varargin{:})));
        end
        function [Ls,DxLs] = Ls(ef,r,varargin)
            % LS Spin density [kg m^2/s]
            %
            % LS = LS(EF,R) calculates the spin density LS at positions R (Point).
            %   S is a ComplexVector.
            %
            % [LS,DxLS = LS(EF,R) calculates the spin density LS and its rotor 
            %   at positions R (Point).
            %   S and DxLS are a ComplexVector.
            %
            % S = S(EF,R,'dr',DR) sets the increment in the calcualtion of
            %   the derivatives to DR [default = 1e-12 m]. 
            %
            % See also EField, Point, ComplexVector.

            Ls = -1i*PhysConst.e0*ef.er/(4*ef.omega)*ef.E(r,varargin{:})*conj(ef.E(r,varargin{:}));
            
            if nargout>1
                % increment [m]
                dr = 1e-12;
                for n = 1:2:length(varargin)
                    if strcmpi(varargin{n},'dr')
                        dr = varargin{n+1};
                        Check.isnumeric('dr must be a positive real number',dr,'>',0)
                    end
                end
                
                dxLs = ( ef.Ls(r+Point(dr,0,0),varargin{:}) - ef.Ls(r+Point(-dr,0,0),varargin{:}) )./(2*dr);
                dyLs = ( ef.Ls(r+Point(0,dr,0),varargin{:}) - ef.Ls(r+Point(0,-dr,0),varargin{:}) )./(2*dr);
                dzLs = ( ef.Ls(r+Point(0,0,dr),varargin{:}) - ef.Ls(r+Point(0,0,-dr),varargin{:}) )./(2*dr);

                DxLs = ComplexVector(r.X,r.Y,r.Z, ...
                    dyLs.Vz-dzLs.Vy, ...
                    -dxLs.Vz+dzLs.Vx, ...
                    dxLs.Vy-dyLs.Vx ...
                    );
            end
        end
    end
end