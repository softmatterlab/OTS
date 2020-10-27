classdef EFieldPlaneWave < EField
    % EFieldPlaneWave < EField : Electromagnetic field of a plane wave
    %
    % EFieldPlaneWave properties:
    %   lambda0 -   vacuum wavelength [m] < EField
    %   er      -   relative dielectric permittivity < EField
    %   mr      -   relative magnetic permeability < EField
	%   E0      -   complex amplitude [V/m]
	%   uk      -   propagation direction (Vector) [m]
	%   e       -   polarisation (ComplexVector)
    %
    % EFieldPlaneWave methods:
    %   EFieldPlaneWave     -   constructor 
    %   n                   -   refractive index < EField
    %   lambda              -   wavelenght in the medium [m]  < EField
    %   k                   -   wave number in the medium [m^-1]  < EField
    %   omega               -   angular frequency [Hz]  < EField
    %   B                   -   magnetic field [T] < EField
    %   S                   -   Poynting vector (Vector) [W/m^2] < EField
    %   Ls                  -   spin density [kg m^2/s] < EField
    %   E                   -   electric field [V/m]
    %
    % See also EField.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        E0  % complex amplitude [V/m]
        uk  % propagation direction (Vector) [m]
        e   % polarisation (ComplexVector)
    end
    methods
        function ef = EFieldPlaneWave(E0,uk,e,varargin)
            % EFIELDPLANEWAVE(E0,UK,E) constructs a plane wave with
            %   amplitude E0, propagation direction UK (Vector)
            %   and polarization direction E (Vector).
            %
            % EFIELDPLANEWAVE(E0,UK,E,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also EFieldPlaneWave, EField.

            ef = ef@EField(varargin{:});
            
            Check.isa('The propagation direction must be a Vector',uk,'Vector')
            Check.isa('The polarization direction must be a Vector',e,'Vector')

            ef.E0 = E0;
            ef.uk = normalize(uk);
            ef.e = normalize(e-(e.*ef.uk)*ef.uk);
        end
        function E = E(ef,r)
            % E Electric field [V/m]
            %
            % E = E(EF,R) calculates the electric field at positions R (Point).
            %   E is a ComplexVector.
            %
            % See also EFieldPlaneWave, Point, ComplexVector.
            
            Check.isa('The set of positions where to calculate the electric field must be a Point',r,'Point')

            E = ef.E0*exp(1i*(ef.k()*ef.uk.topoint().*r))*ef.e;
            E.X = r.X;
            E.Y = r.Y;
            E.Z = r.Z;
        end
    end
end