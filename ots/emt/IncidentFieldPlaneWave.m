classdef IncidentFieldPlaneWave < IncidentField
    % IncidentFieldPlaneWave < IncidentField : Plane wave
    %   A plane wave is defined by its propagation and polarization vectors.
    %
    % IncidentFieldPlaneWave properties:
    %   k0	-	vacuum wavenumber < IncidentField
    %   nm  -	medium refractive index < IncidentField
    %   ki	-   propagation unit vector (ComplexVector)
    %   Ei	-	polarization unit vector (ComplexVector)
    %
    % IncidentFieldPlaneWave methods:
    %   IncidentFieldPlaneWave	-	constructor
    %   angles                  -   plane wave polar and azimutal angles related to propagation direction
    %   Wi                      -   get multipole coefficients (Coefficients)
    %   incoming_exact          -   incoming field (exact) 
    %   incoming                -   incoming field < IncidentField
    %   incoming_expansion      -   incoming field calculated with multipole expansion < IncidentField
    %
    % See also IncidentField, IncidentFieldFocusedBeam.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        ki  % propagation unit vector (ComplexVector)
        Ei  % polarization unit vector (ComplexVector)
    end
    methods
        function obj = IncidentFieldPlaneWave(ki,Ei,k0,nm)
            % INCIDENTFIELDPLANEWAVE(Ki,Ei,K0,NM) construct a plane wave
            %   that propagates in Ki direction (ComplexVector) and 
            %   has polarization in Ei direction (ComplexVector).
            %   K0 is the medium wavenumber and NM is the medium refractive index.
            %
            % See also IncidentFieldPlaneWave, IncidentField, IncidentFieldFocusedBeam.
            
            obj = obj@IncidentField(k0,nm);           
            obj.ki = ki;
            obj.Ei = Ei;
        end
        function [theta,phi] = angles(field)
            % ANGLES Angles related to the plane wave propagation direction
            %
            % [THETA,PHI] = ANGLES(FIELD) determine the plane wave 
            %   polar (THETA) and azimutal (PHI) angles related to the propagation direction
            %
            % Sea also IncidentFieldPlaneWave, IncidentField, IncidentFieldFocusedBeam.

            [theta,phi] = Transform.Car2Sph(field.ki.Vx,field.ki.Vy,field.ki.Vz);
        end
        function coeff = Wi(field,L)
            % WI Incident field amplitudes
            %
            % COEFF = WI(FIELD,L) calculates the J-multipole expansion
            %   amplitude coefficients of the incident plane wave FIELD.
            %   L is maximum number of multipoles.
            %
            % See also IncidentFieldPlaneWave, IncidentField.

            [theta,phi] = field.angles();
            
            C = sparse(zeros(Coefficients.index(2,L,L),1));
            
            vsh = VecSpHarm(theta,phi);
            for l = 1:1:L
                for m = -l:1:l
                    C(Coefficients.index(1,l,m)) = 4*pi*1i^l * field.Ei.*conj(vsh.Z1(l,m));
                    C(Coefficients.index(2,l,m)) = 4*pi*1i^(l+1) * field.Ei.*conj(vsh.Z2(l,m));
                end
            end
            
            coeff = Coefficients(C);
        end
        function [Ei,Bi] = incoming_exact(field,theta,phi,r)
            % INCOMING_EXACT Incoming plane wave (exact) 
            %
            % [Ei,Bi] = INCOMING_EXACT(FIELD,THETA,PHI,R) calculates
            %   incoming electric and magnetic fields Ei and Bi at coordinates THETA, PHI, R.
            %   The incoming electric field FIELD has amplitude Ei = 1 W/m, 
            %   is propagating in the z-direciton and is lineraly polarized
            %   along the x-direction.
            %
            % See also IncidentFieldPlaneWave, IncidentField.
            
            km = field.nm*field.k0;
            
            [x,y,z] = Transform.Sph2Car(theta,phi,r);
            
            % x-polarization
            Ei = ComplexVector(x,y,z,exp(1i*km*z),zeros(size(theta)),zeros(size(theta)));
            Bi = ComplexVector(x,y,z,zeros(size(theta)),field.nm/PhysConst.c0*exp(1i*km*z),zeros(size(theta)));
        end
    end
end