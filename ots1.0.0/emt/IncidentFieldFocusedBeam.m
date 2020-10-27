classdef IncidentFieldFocusedBeam < IncidentField
    % IncidentFieldFocusedBeam < IncidentField : Focused paraxial beam
    %   The properties of the focal field are defined by the beam, the
    %   focal length of the objective and refraction index of medium.
    %
    % IncidentFieldFocusedBeam properties:
    %   k0     -   vacume wavenumber < IncidentField
    %   nm     -   medium refractive index < IncidentField
    %   b      -   Paraxial beam (Beam)
    %   f      -   Focal length [m]
    %
    % IncidentFieldFocusedBeam methods:
    %   IncidentFieldFocusedBeam       -  constructor
    %   Wi                             -  incident field amplitudes
    %   incoming_exact                 -  incoming plane wave 
    %   incoming                       -  incoming field <  IncidentField
    %   incoming_expansion             -  Incoming field calculated with multipole expansion < IncidentField
    %
    % See also IncidentField, IncidentFieldPlaneWave.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        b   % Paraxial beam (Beam)
        f   % Focal length [m]
    end
    methods
        function obj = IncidentFieldFocusedBeam(b,f,nm)
            % INCIDENTFIELDFOCUSEDBEAM(B,F,NM) constructs a focused beam.
            %   B is the beam, F is the focal length of the objective 
            %   and NM is medium refractive index.
            %
            % See also IncidentFieldFocusedBeam, Beam, IncidentField, IncidentFieldPlaneWave.
            
            k0 = 2*pi/b.lambda0;
            obj = obj@IncidentField(k0,nm);

            obj.b = b;
            obj.f = f;
        end
        function [coeff,coeff_pw] = Wi(field,L,P,coeff_pw)
            % Wi Focused beam amplitudes 
            % 
            % [COEFF,COEFF_PW] = Wi(FIELD,L,P,COEFF_PW) calculates the J-multipole expansion
            %   amplitude coefficients of the focused beam FIELD.
            %   L is number of coefficients.
            %   P (optional, default P=Point(0,0,0)) is the position of the focus (Point).
            %   COEFF_PW (optional) are the plane wave coefficients.
            %
            % See also IncidentFieldFocusedBeam, Beam, IncidentFieldPlaneWave, IncidentField.
            
            % Plane wave coefficients
            if nargin<4
                coeff_pw = [];
            end
            
            % P is the origin of the reference frame
            if nargin<3
                P = Point(0,0,0);
            end
            
            phi = field.b.phi;
            dphi = phi(2,1)-phi(1,1);
            
            theta = asin(field.b.r/field.f);
            dr = field.b.r(1,2)-field.b.r(1,1);
            R = dr*size(field.b.r,2);
            dtheta = ones(size(field.b.r,1),1)*(asin([dr:dr:R]/field.f)-asin([0:dr:R-dr]/field.f));
            
            utheta = Point(cos(phi).*cos(theta), sin(phi).*cos(theta), -sin(theta)).tovector();
            uphi = Point(-sin(phi), cos(phi), zeros(size(phi))).tovector();
            
            kt = 2*pi*field.nm/field.b.lambda0;
            ktx = kt*sin(theta).*cos(phi);
            kty = kt*sin(theta).*sin(phi);
            ktz = kt*cos(theta);
            
            nb = sqrt(field.b.er*field.b.mr);
            Et = (field.b.Ephi.*uphi + field.b.Er.*utheta)*(sqrt(nb/field.nm).*sqrt(cos(theta)));
            % Et = ComplexVector(X,Y,Z,Et.Vx,Et.Vy,Et.Vz);
            i = Coefficients.index(2,L,L);
            coeff = Coefficients(zeros(1,i));
            for n = 1:1:numel(Et)
            
                if length(coeff_pw)<n
                    
                    ki = ComplexVector(0,0,0,ktx(n),kty(n),ktz(n));
            
                    Et_n = ComplexVector(0,0,0,Et.Vx(n),Et.Vy(n),Et.Vz(n));

                    pw = IncidentFieldPlaneWave(ki,Et_n,field.k0,field.nm);
                    coeff_pw(n).Wi = pw.Wi(L);
                end
            
                % Phase factor due to the reference frame shift
                phase = exp( 1i * (ktx(n)*P.X+kty(n)*P.Y+ktz(n)*P.Z) );
            
                coeff = coeff + coeff_pw(n).Wi * phase * sin(theta(n))*dphi*dtheta(n);
            end
            coeff = 1i*kt*field.f*exp(-1i*kt*field.f)/(2*pi)*coeff;
        end
        function [Ei,Bi] = incoming_exact(field,theta,phi,r)
            % INCOMING_EXACT Incoming focused beam (exact)
            %
            % [Ei,Bi] = INCOMING_EXACT(FIELD,THETA,PHI,R) calculates the
            %   focal fields Ei and Bi at coordinates THETA, PHI, R.
            %
            % !!! Note that at the moment Bi is not calculated !!!
            %   
            % See also IncidentFieldFocusedBeam, IncidentFieldPlaneWave,MieParticle.

            [X, Y, Z] = Transform.Sph2Car(theta,phi,r);
            Ei = field.b.focus(field.f,X,Y,Z,'er',field.nm^2,'mr',1);
        end
    end
    
end