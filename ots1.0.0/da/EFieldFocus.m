classdef EFieldFocus < EField
    % EFieldFocus < EField : Electromagnetic field of a focused beam
    %
    % EFieldFocus properties:
    %   lambda0 -   vacuum wavelength [m] < EField
    %   er      -   relative dielectric permittivity < EField
    %   mr      -   relative magnetic permeability < EField
    %   beam    -   Beam
	%   f       -   focal length [m]
    %
    % EFieldFocus methods:
    %   EFieldFocus -   constructor 
    %   n           -   refractive index < EField
    %   lambda      -   wavelenght in the medium [m]  < EField
    %   k           -   wave number in the medium [m^-1]  < EField
    %   omega       -   angular frequency [Hz]  < EField
    %   B           -   magnetic field [T] < EField
    %   S           -   Poynting vector (Vector) [W/m^2] < EField
    %   Ls      	-   spin density [kg m^2/s] < EField
    %   E           -   electric field [V/m]
    %
    % See also EField, Beam.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        beam    % Beam
        f       % focal length [m]
    end
    methods
        function ef = EFieldFocus(beam,f,varargin)
            % EFIELDFOCUS(BEAM,F) constructs the focal field of BEAM
            %   obtained using an aplanatic objective with focal lenght F.
            %
            % EFIELDFOCUS(BEAM,F,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       lambda0     -   vacuum wavelength [default: 532e-9 m]
            %       er          -   relative electric permittivity [default: 1]
            %       mr          -   relative magnetic permeability [default: 1]
            %
            % See also EFieldFocus, EField.

            ef = ef@EField(varargin{:});

            Check.isa('The optical beam must be a Beam',beam,'Beam')
            Check.isreal('The focal length must be a positive real number',f,'>',0)
            
            ef.beam = beam;
            ef.f = f;
        end
        function E = E(ef,r,varargin)
            % E Electric field [V/m]
            %
            % E = E(EF,R) calculates the electric field at positions R (Point).
            %   E is a ComplexVector.
            %
            % See also EFieldFocus, Beam.focus, Point, ComplexVector.

            Check.isa('The set of positions where to calculate the electric field must be a Point',r,'Point')

            E = ef.beam.focus(ef.f,r.X,r.Y,r.Z,'er',ef.er,'mr',ef.mr);
        end
    end
end