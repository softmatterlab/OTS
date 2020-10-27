classdef IncidentField
    % IncidentField (Abstract) :  Incident field
    %   Incident field taht illuminates the particle.
    %   Instances of this class cannot be created. Use one of the subclasses 
    %   (e.g., IncidentFieldPlaneWave, IncidentFieldFocusedBeam).
    %
    % IncidentField properties:
    %   k0       -      vacuum wavenumber 
    %   nm       -      medium refractive index
    %
    % IncidentField abstract methods:
    %   Wi              -   get multipole coefficients (Coefficients)
    %   incoming_exact  -   incoming field (exact) 
    %
    % IncidentField methods:
    %   IncidentField       -   constructor
    %   incoming            -   incoming field  
    %   incoming_expansion	-   Incoming field calculated with multipole expansion  
    %
    % See also IncidentFieldPlaneWave, IncidentFieldFocusedBeam.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        k0  % vacuum wavelength
        nm  % medium refactive index
    end
    methods (Access = protected)
        function obj = IncidentField(k0,nm)
            % INCIDENTFIELD(K0,NM) constructs a incident field. 
            %   K0 is medium wave number and NM is medium refractive index.
            %
            % See also IncidentField, IncidentFieldPlaneWave, IncidentFieldFocusedBeam.

            obj.k0 = k0;
            obj.nm = nm;
        end
    end
    methods (Abstract)
        Wi(field,L)  % get multipole coefficients (Coefficients)
        incoming_exact(field,theta,phi,r)  % incoming field (exact)  
    end
    methods
        function [Ei,Bi] = incoming(field,L,theta,phi,r,varargin)
            % INCOMING Incoming field  
            % 
            % [Ei,Bi] = INCOMING(FIELD,L,THETA,PHI,R) determine the
            %   incoming field at coordinates THETA, PHI and R.
            %   L is maximum of the multipole index for which the multipole
            %   expansion is calculated.
            %
            % [Ei,Bi] = INCOMING(FIELD,L,THETA,PHI,R,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       multipoles  -  Multipole
            %
            % See also IncidentField, IncidentFieldPlaneWave, IncidentFieldFocusedBeam, Multipole.
            
            if L==Inf
                [Ei,Bi] = incoming_exact(field,theta,phi,r);
            else
                [Ei,Bi] = incoming_expansion(field,L,theta,phi,r,varargin{:});
            end
        end
        function [Ei,Bi] = incoming_expansion(field,L,theta,phi,r,varargin)
            % INCOMING_EXPANSION Incoming field calculated with multipole expansion  
            % 
            % [Ei,Bi] = INCOMING_EXPANSION(FIELD,L,THETA,PHI,R) determines the
            %   incoming field at coordinates  THETA,PHI and R
            %   using a multipole expansion up to the multipole order L.
            %
            % [Ei,Bi] = INCOMING_EXPANSION(FIELD,L,THETA,PHI,R,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       multipoles  -   Multipole
            %
            % See also IncidentField, IncidentFieldPlaneWave, IncidentFieldFocusedBeam, Multipole.
            
            % multipoles
            multi = Multipole(theta,phi,r,field.nm*field.k0,varargin{:});
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'multipoles')
                    multi = varargin{n+1};
                end
            end
            
            % [x,y,z] = Transform.Sph2Car(theta,phi,r);
            x = multi.X;
            y = multi.Y;
            z = multi.Z;
            
            Wi = field.Wi(L);
            
            Ei = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            Bi = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            for l = 1:1:L
                for m = -l:1:l
                    Ei = Ei + ...
                        Wi.C(Coefficients.index(1,l,m)) * multi.J1(l,m) + ...
                        Wi.C(Coefficients.index(2,l,m)) * multi.J2(l,m);
                    Bi = Bi + ...
                        -1i*field.nm/PhysConst.c0 * Wi.C(Coefficients.index(2,l,m)) * multi.J1(l,m) + ...
                        -1i*field.nm/PhysConst.c0 * Wi.C(Coefficients.index(1,l,m)) * multi.J2(l,m);
                end
            end
            
        end
    end
end