classdef TMatrix
    % TMatrix (Abstract) : Transition matrix (T-matrix)
    %   A transition matrix relates the coefficients of incident and scattered field.
    %   Instances of this class cannot be created. Use one of its subclasses 
    %   (e.g., TMatrixSphere, TMatrixCluster, TMatrixInclusion).
    %
    % TMatrix properties:
    %   T             -   T-Matrix (matrix)
    %
    % TMatrix methods:
    %   TMatrix       -   constructor (accessible only by the subclasses)
    %   Li            -   incident field amplitude coefficient number
    %   Ls            -   scattered field amplitude coefficient number
    %   As            -   scattered field coefficients
    %   force         -   optical force
    %   torque        -   optical torque
    %   incoming      -   incoming field
    %   scattering    -   scattered field
    %   total         -   total field
    %   rotate        -   T-matrix rotation
    %   translate     -   T-matrix translation
    %
    % TMatrix abstract methods:
    %   k0            -   medium wavenumber
    %   nm            -   medium refraction index
    %
    % TMatrix static methods:
    %   index         -   index coressponding to plm
    %   indices       -   indices coressponding to plm and p'l'm'
    %   plms          -   plm and p'l'm' corresponding to indices
    %   O             -   function O (for force calculation)
    %   C1            -   Clebsch-Gordan coefficients for l1 = 1
    %   K             -   function K (for force calculation)
    %   I             -   function I (for force calculation)
    %   wigner        -   Wigner rotation matrix
    %   Gaunt         -   Gaunt integral
    %   Gj            -   Nozawa's addition theorem for Bessel harmonic
    %   Gh            -   Nozawa's addition theorem for Hankel harmonic
    %   JJ            -   translation matrix coresponding to Gj
    %   HH            -   translation matrix coresponding to Gh
    %
    % See also TMatrix, TMatrixSphere, TMatrixCluster, TMatrixInclusion.
    
    %   Author: S. Masoumeh Mousavi, Agnese Callegari, Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        T   % T-Matrix
    end
    methods (Access = protected)
        function tm = TMatrix()
            % TMATRIX() constructs a TMatrix
            %   This method is only accessible by the subclasses of TMatrix
            %   (e.g., TMatrixSphere, TMatrixCluster, TMatrixInclusion).
            %
            % See also TMatrix, TMatrixSphere, TMatrixCluster, TMatrixInclusion.
        end
    end
    methods
        function L = Li(tm)
            % LI Incident field amplitude coefficient number
            %
            % L = LI(tm) determines number of coefficients of the incident
            %   field amplitude.
            %
            % See also TMatrix.

            L = floor(sqrt((size(tm.T,2)-1)/2));
        end
        function L = Ls(tm)
            % LS Scattered field amplitude coefficient number
            %
            % L = LS(tm) determines number of coefficients of the scattered
            %   field amplitude.
            %
            % See also TMatrix.

            L = floor(sqrt((size(tm.T,1)-1)/2));
        end
        function coeff = As(tm,Wi)
            % AS Scattered field coefficients
            %
            % coeff = AS(tm,Wi) determines the coefficients of the scattered field
            %   Wi is the incident field amplitude.
            %
            % See also TMatrix, IncidentField, IncidentFieldFocusedBeam.

            coeff = Coefficients(tm.T*Wi.C);
            
            % L = Wi.lmax();
            % 
            % A1 = zeros(L,2*L+1);
            % A2 = zeros(L,2*L+1);
            % 
            % for lp = 1:1:L
            %     for mp = -lp:1:lp
            %         for l = 1:1:L
            %             for m = -l:1:l
            %                 A1(lp,mp+L+1) = A1(lp,mp+L+1) ...
            %                     + tm.T(1,1,l,lp,m,mp)*Wi.c1(l,m) ...
            %                     + tm.T(2,1,l,lp,m,mp)*Wi.c2(l,m);
            %                 A2(lp,mp+L+1) = A2(lp,mp+L+1) ...
            %                     + tm.T(1,2,l,lp,m,mp)*Wi.c1(l,m) ...
            %                     + tm.T(2,2,l,lp,m,mp)*Wi.c2(l,m);
            %             end
            %         end
            %     end
            % end
            % 
            % coeff = Coefficients(A1,A2);
        end
        function F = force(tm,Wi)
            % FORCE Optical force
            %
            % F = FORCE(tm,Wi) calculates the optical force.
            %   Wi is the incident field amplitude.
            %
            % See also TMatrix, IncidentField, IncidentFieldFocusedBeam.

            Li = tm.Li();
            Ls = tm.Ls();
            
            As = tm.As(Wi);
            
            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:size(tm.T,1)],[1:1:size(tm.T,2)]);
            [PP,P] = meshgrid(p,pp);
            [LP,L] = meshgrid(l,lp);
            [MP,M] = meshgrid(m,mp);
            
            pre = (1i.^l .* (As.C+Wi.C)');
            post = (1i.^-lp' .* As.C);
            
            % These formulas wok also in a medium as
            % em = e0*nm^2
            % and
            % km = k0*nm
            F = Vector(0,0,0, ...
                -0.5 * PhysConst.e0 * tm.k0^-2 * real( pre * TMatrix.I(pi/2,0,P,PP,L,LP,M,MP) * post ) ...
                ,...
                -0.5 * PhysConst.e0 * tm.k0^-2 * real( pre * TMatrix.I(pi/2,pi/2,P,PP,L,LP,M,MP) * post ) ...
                , ...
                -0.5 * PhysConst.e0 * tm.k0^-2 * real( pre * TMatrix.I(0,0,P,PP,L,LP,M,MP) * post ) ...
                );

            % Fx = 0;
            % Fy = 0;
            % Fz = 0;
            % for lp = 1:1:Ls
            %     for mp = -lp:1:lp
            %         for l = 1:1:Li
            %             for m = -l:1:l
            %                 Fx = Fx - ( ...
            %                     conj(As.c1(l,m)+Wi.c1(l,m)) * As.c1(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,0,1,1,l,lp,m,mp) ...
            %                     + conj(As.c1(l,m)+Wi.c1(l,m)) * As.c2(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,0,1,2,l,lp,m,mp) ...
            %                     + conj(As.c2(l,m)+Wi.c2(l,m)) * As.c1(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,0,2,1,l,lp,m,mp) ...
            %                     + conj(As.c2(l,m)+Wi.c2(l,m)) * As.c2(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,0,2,2,l,lp,m,mp) ...
            %                     );
            %                 Fy = Fy - ( ...
            %                     conj(As.c1(l,m)+Wi.c1(l,m)) * As.c1(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,pi/2,1,1,l,lp,m,mp) ...
            %                     + conj(As.c1(l,m)+Wi.c1(l,m)) * As.c2(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,pi/2,1,2,l,lp,m,mp) ...
            %                     + conj(As.c2(l,m)+Wi.c2(l,m)) * As.c1(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,pi/2,2,1,l,lp,m,mp) ...
            %                     + conj(As.c2(l,m)+Wi.c2(l,m)) * As.c2(lp,mp) * 1i^(l-lp) * TMatrix.I(pi/2,pi/2,2,2,l,lp,m,mp) ...
            %                     );
            %                 Fz = Fz - ( ...
            %                     conj(As.c1(l,m)+Wi.c1(l,m)) * As.c1(lp,mp) * 1i^(l-lp) * TMatrix.I(0,0,1,1,l,lp,m,mp) ...
            %                     + conj(As.c1(l,m)+Wi.c1(l,m)) * As.c2(lp,mp) * 1i^(l-lp) * TMatrix.I(0,0,1,2,l,lp,m,mp) ...
            %                     + conj(As.c2(l,m)+Wi.c2(l,m)) * As.c1(lp,mp) * 1i^(l-lp) * TMatrix.I(0,0,2,1,l,lp,m,mp) ...
            %                     + conj(As.c2(l,m)+Wi.c2(l,m)) * As.c2(lp,mp) * 1i^(l-lp) * TMatrix.I(0,0,2,2,l,lp,m,mp) ...
            %                     );
            %             end
            %         end
            %     end
            % end
            %
            % % These formulas wok also in a medium as
            % % em = e0*nm^2
            % % and
            % % km = k0*nm
            % Fx = 0.5 * PhysConst.e0 * tm.k0^-2 * real(Fx);
            % Fy = 0.5 * PhysConst.e0 * tm.k0^-2 * real(Fy);
            % Fz = 0.5 * PhysConst.e0 * tm.k0^-2 * real(Fz);
            % 
            % F = Vector(0,0,0,Fx,Fy,Fz);
        end
        function T = torque(tm,Wi)
            % TORQUE Optical torque
            %
            % T = TORQUE(tm,Wi) Calculates the optical torque.
            %   Wi is the incident field amplitude.
            %
            % See also TMatrix, IncidentField, IncidentFieldFocusedBeam.
            
            Li = tm.Li();
            Ls = tm.Ls();
            
            As = tm.As(Wi);
            
            WCu = sparse(zeros(Coefficients.index(2,Li,Li),1));
            WCd = sparse(zeros(Coefficients.index(2,Li,Li),1));
            ACu = sparse(zeros(Coefficients.index(2,Ls,Ls),1));
            ACd = sparse(zeros(Coefficients.index(2,Ls,Ls),1));
            for p = [1 2]
                for l = 1:1:Ls
                    for m = -l:1:l
                        if m+1 <= l
                            WCu(Coefficients.index(p,l,m)) = Wi.C(Coefficients.index(p,l,m+1));
                            ACu(Coefficients.index(p,l,m)) = As.C(Coefficients.index(p,l,m+1));
                        end
                        if m-1 >= -l
                            WCd(Coefficients.index(p,l,m)) = Wi.C(Coefficients.index(p,l,m-1));
                            ACd(Coefficients.index(p,l,m)) = As.C(Coefficients.index(p,l,m-1));
                        end
                    end
                end
            end
            
            Wi_up = WCu;
            Wi_down = WCd;
            As_up = ACu;
            As_down = ACd;
            
            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:size(tm.T,1)],[1:1:size(tm.T,2)]);
            
            smin = sqrt((l-m).*(l+1+m));
            spos = sqrt((l+m).*(l+1-m));
            
            T = Vector(0,0,0, ...
                -0.25 * PhysConst.e0 * tm.k0^-3*tm.nm^-1 *( real(smin*(Wi_up.*conj(As.C)+As_up.*conj(As.C))...
                +spos*(Wi_down.*conj(As.C)+As_down.*conj(As.C)))) ...
                ,...
                -0.25 * PhysConst.e0 * tm.k0^-3*tm.nm^-1 *( imag(-smin*(Wi_up.*conj(As.C)+As_up.*conj(As.C))...
                +spos*(Wi_down.*conj(As.C)+As_down.*conj(As.C)) )) ...
                , ...
                -0.5 * PhysConst.e0 * tm.k0^-3*tm.nm^-1 *( m*(real(Wi.C.*conj(As.C))+abs(As.C).^2))) ...
                ;
            
        end
        function [Ei,Bi] = incoming(tm,Wi,theta,phi,r,varargin)
            % INCOMING Incoming field
            %
            % [Ei,Bi] = INCOMING(tm,Wi,THETA,PHI,R) Calculates incident
            %   electric and magnetic fields Ei and Bi at coordinates THETA, PHI and R.
            %   THETA is the polar angle [0 pi], PHI is the azimuthal angle [-pi pi] 
            %   and R is a radial distance.
            %   Wi is the incident field amplitude.
            %
            % [Ei,Bi] = INCOMING(tm,Wi,THETA,PHI,R,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       multipoles      -   Multipole variable (for preallocation)
            %
            % See also TMatrix, IncidentField, IncidentFieldFocusedBeam, Multipole.

            % multipoles
            multi = Multipole(theta,phi,r,tm.nm*tm.k0,varargin{:});
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'multipoles')
                    multi = varargin{n+1};
                end
            end
            
            % [x,y,z] = Transform.Sph2Car(theta,phi,r);
            x = multi.X;
            y = multi.Y;
            z = multi.Z;
                        
            Es = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            Bs = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            for l = 1:1:tm.Ls
                for m = -l:1:l
                    Es = Es + ...
                        Wi.C(Coefficients.index(1,l,m)) * multi.J1(l,m) + ...
                        Wi.C(Coefficients.index(2,l,m)) * multi.J2(l,m);
                    Bs = Bs + ...
                        -1i*tm.nm/PhysConst.c0 * Wi.C(Coefficients.index(2,l,m)) * multi.J1(l,m) + ...
                        -1i*tm.nm/PhysConst.c0 * Wi.C(Coefficients.index(1,l,m)) * multi.J2(l,m);
                end
            end
            
        end
        function [Es,Bs] = scattering(tm,Wi,theta,phi,r,varargin)
            % SCATTERING Scattered field
            %
            % [Es,Bs] = SCATTERING(tm,Wi,THETA,PHI,R) Calculates scattered
            %   electric and magnetic fields Es and Bs at coordinates THETA, PHI and R.
            %   THETA is the polar angle [0 pi], PHI is the azimuthal angle [-pi pi] 
            %   and R is a radial distance.
            %   Wi is the incident field amplitude.
            %
            % [Es,Bs] = SCATTERING(tm,Wi,theta,phi,r,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       multipoles      -   Multipole variable (for preallocation)
            %
            % See also TMatrix, IncidentField, IncidentFieldFocusedBeam, Multipole.
            
            % multipoles
            multi = Multipole(theta,phi,r,tm.nm*tm.k0,varargin{:});
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'multipoles')
                    multi = varargin{n+1};
                end
            end
            
            % [x,y,z] = Transform.Sph2Car(theta,phi,r);
            x = multi.X;
            y = multi.Y;
            z = multi.Z;
            
            As = tm.As(Wi);
            
            Es = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            Bs = ComplexVector(x,y,z,zeros(size(theta)),zeros(size(theta)),zeros(size(theta)));
            for l = 1:1:tm.Ls
                for m = -l:1:l
                    Es = Es + ...
                        As.C(Coefficients.index(1,l,m)) * multi.H1(l,m) + ...
                        As.C(Coefficients.index(2,l,m)) * multi.H2(l,m);
                    Bs = Bs + ...
                        -1i*tm.nm/PhysConst.c0 * As.C(Coefficients.index(2,l,m)) * multi.H1(l,m) + ...
                        -1i*tm.nm/PhysConst.c0 * As.C(Coefficients.index(1,l,m)) * multi.H2(l,m);
                end
            end
            
        end
        function [Et,Bt] = total(tm,Wi,theta,phi,r,varargin)
            % TOTAL Total field
            %
            % [Et,Bt] = TOTAL(tm,Wi,THETA,PHI,R) calculates total electric
            %   and magnetic fields at coordinates THETA, PHI, R. 
            %   THETA is the polar angle [0 pi], PHI is the azimuthal angle [-pi pi] 
            %   and R is a radial distance.
            %   Wi is the incident field amplitude.
            %
            % [Et,Bt] =  TOTAL(tm,Wi,THETA,PHI,R,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       multipoles      -   Multipole variable (for preallocation)
            %
            % See also TMatrix, IncidentField, IncidentFieldFocusedBeam, Multipole.
            
            [Ei,Bi] = tm.incoming(theta,phi,r);
            [Es,Bs] = tm.scattering(Wi,theta,phi,r,varargin{:});
            Et = Ei+Es;
            Bt = Bi+Bs;
        end
        function tm_r = rotate(tm,alpha,beta,gamma,lgf)
            % ROTATE T-matrix rotation
            %
            % tm_r = ROTATE(tm,alpha,beta,gamma,lgf) rotate the T-matrix coefficients 
            %   from the lab reference frame to the particle reference frame 
            %   by rotating 
            %   (1) alpha around the z axis,
            %   (2) beta around the y axis and
            %   (3) gamma around the z axis.
            %   lgf is logaritm of factorial (for preallocation).
            % 
            % see also TMatrix, LGF.
                        
            % [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:size(tm.T,1)],[1:1:size(tm.T,2)]);
            % [P,PP] = meshgrid(p,pp);
            % [L,LP] = meshgrid(l,lp);
            % [M,MP] = meshgrid(m,mp);
            % 
            % tm_r = tm;
            % D = tm.wigner(P,PP,L,LP,M,MP,alpha,beta,gamma);
            % tm_r.T = D' * tm.T * D;

            % because we need a translation from the lab reference frame to the particle reference frame
            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:size(tm.T,1)],[1:1:size(tm.T,2)]);
            [PP,P] = meshgrid(p,pp);
            [LP,L] = meshgrid(l,lp);
            [MP,M] = meshgrid(m,mp);
            
            tm_r = tm;
            D = tm.wigner(P,PP,L,LP,M,MP,lgf,alpha,beta,gamma);
            tm_r.T = D * tm.T * D';
            
        end
        function tm_t = translate(tm,R,rp,cg,varargin)
            % TRANSLATE T-matrix translation
            %
            % tm_t = TRANSLATE(tm,R,rp,cg) translate the T-matrix coefficients 
            %   from the lab reference frame to the particle reference frame.
            %   R is the tranlation vector, rp is a radial distance for
            %   computing the scattering (it shoudl enclose the cluster) 
            %   and cg are the Clebsch-Gordan coefficients (for
            %   preallocation).
            %
            % tm_t = TRANSLATE(tm,R,rp,cg,'PropertyName',PropertyValue) sets the property
            %   PropertyName to PropertyValue. The properties listed below
            %   can be used:
            %       Lmax    -   Maximum value of L for Gh and Gj.
            %
            % See also TMatrix, CG.
            
            R = -R;  % because we need a translation from the lab reference frame to the particle reference frame
            
            km = tm.nm*tm.k0;
            
            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:size(tm.T,1)],[1:1:size(tm.T,2)]);
            % [P,PP] = meshgrid(p,pp);
            [L,LP] = meshgrid(l(1:2:end),lp(1:2:end));
            [M,MP] = meshgrid(m(1:2:end),mp(1:2:end));
            
            tm_t = tm;
            if norm(rp) < norm(R);
                tm_t.T = TMatrix.HH(LP,L,MP,M,R,km,cg,varargin{:})*tm.T*TMatrix.JJ(LP,L,MP,M,-R,km,cg,varargin{:});
            else % norm(rp) >= norm(R);
                tm_t.T = TMatrix.JJ(LP,L,MP,M,R,km,cg,varargin{:})*tm.T*TMatrix.JJ(LP,L,MP,M,-R,km,cg,varargin{:});
            end
        end
    end    
    methods (Abstract)
        k0(tm)  % vacuum wavelength
        nm(tm)  % medium refractive index
    end    
    methods (Static)
        function i = index(p,l,m)
            % INDEX Scalar index corresponding to p,l,m
            %
            % I = INDEX(P,L,M) calculates the index corresponding to P,L,M. 
            %   P,L,M are the indices of multipole fields.
            %
            % See also TMatrix, Coefficients.

            i = Coefficients.index(p,l,m);
        end
        function [i,ip] = indices(p,pp,l,lp,m,mp)
            % INDICES Scalar indeces corresponding to p,l,m and p',l',m'
            %
            % [I,Ip] = INDICES(P,Pp,L,Lp,M,Mp) calculates the indices
            %   corresponding to P,L,M and Pp,Lp,Mp.
            %
            % See also TMatrix, Coefficients.

            i = TMatrix.index(p,l,m);
            ip = TMatrix.index(pp,lp,mp);
        end
        function [p,pp,l,lp,m,mp] = plms(i,ip)
            % PLMS Multipole p,l,m and p',l',m' corresponding to scalar indices
            %
            % [P,Pp,L,Lp,M,Mp] = PLMS(I,Ip) calculates multipole indices
            %   P,L,M and Pp, Lp, Mp corresponding to I and Ip.
            %
            % See also TMatrix, Coefficients.

            [p,l,m] = Coefficients.plm(i);
            [pp,lp,mp] = Coefficients.plm(ip);
        end
        function o = O(p,pp,l,lp)
            % O Function O
            %
            % R = O(P,Pp,L,Lp) calculate the O function for indices for
            %   multipole indices P,Pp,L and Lp.
            %   This function is used in computing optical force. 
            %
            % See also Coefficients, TMatrix
            
            o = sparse(zeros(size(p)));
            
            ind = lp==l-1 & p==pp;
            o(ind) = sqrt( (l(ind)-1).*(l(ind)+1) ./ (l(ind).*(2*l(ind)+1)) );
            
            ind = lp==l & p~=pp;
            o(ind) = -1./sqrt(l(ind).*(l(ind)+1));
            
            ind = lp==l+1 & p==pp;
            o(ind) = -sqrt( (l(ind).*(l(ind)+2)) ./ ((l(ind)+1).*(2*l(ind)+1)) );
            
            % if lp==l-1 && p==pp
            %     o = sqrt( (l-1)*(l+1) / (l*(2*l+1)) );
            % elseif lp==l && p~=pp
            %     o = -1/sqrt(l*(l+1));
            % elseif lp==l+1 && p==pp
            %     o = -sqrt( (l*(l+2)) / ((l+1)*(2*l+1)) );
            % else
            %     % in particular for lp < l-1 || lp > l+1
            %     o = 0;
            % end
        end
        function c = C1(l,lp,mu,mp)
            % C1 Clebsch-Gordan coefficients for l1 = 1
            %
            % R = C1(L,Lp,MU,Mp) determine the Clebsch-Gordan coefficients for
            %   l1 = 1, l2 = L, lp = Lp and m1 = MU, m2 = Mp and M = MU+Mp.
            %
            % See also  TMatrix, Coefficients.
            
            % Note that in the Book is used C1(l,lp,mu,mp-mu)
            
            c = sparse(zeros(size(l)));
            
            ind = mu==-1 & l==lp-1;
            c(ind) = sqrt( (lp(ind)-1-mp(ind)).*(lp(ind)-mp(ind)) ./ (2*lp(ind).*(2*lp(ind)-1)) );
                
            ind = mu==-1 & l==lp;
            c(ind) = -sqrt( (lp(ind)-mp(ind)).*(lp(ind)+1+mp(ind)) ./ (2*lp(ind).*(lp(ind)+1)) );
                
            ind = mu==-1 & l==lp+1;
            c(ind) = sqrt( (lp(ind)+2+mp(ind)).*(lp(ind)+1+mp(ind)) ./ (2*(lp(ind)+1).*(2*lp(ind)+3)) );
                
            ind = mu==0 & l==lp-1;
            c(ind) = sqrt( (lp(ind)-mp(ind)).*(lp(ind)+mp(ind)) ./ (lp(ind).*(2*lp(ind)-1)) );
                
            ind = mu==0 & l==lp;
            c(ind) = -mp(ind)./sqrt(lp(ind).*(lp(ind)+1));
                
            ind = mu==0 & l==lp+1;
            c(ind) = -sqrt( (lp(ind)+1-mp(ind)).*(lp(ind)+1+mp(ind)) ./ ((lp(ind)+1).*(2*lp(ind)+3)) );
                
            ind = mu==+1 & l==lp-1;
            c(ind) = sqrt( (lp(ind)-1+mp(ind)).*(lp(ind)+mp(ind)) ./ (2*lp(ind).*(2*lp(ind)-1)) );
            
            ind = mu==+1 & l==lp;
            c(ind) = sqrt( (lp(ind)+mp(ind)).*(lp(ind)+1-mp(ind)) ./ (2*lp(ind).*(lp(ind)+1)) );
            
            ind = mu==+1 & l==lp+1;
            c(ind) = sqrt( (lp(ind)+1-mp(ind)).*(lp(ind)+2-mp(ind)) ./ (2*(lp(ind)+1).*(2*lp(ind)+3)) );
            
            % if mu==-1 && l==lp-1
            %     c = sqrt( (lp-1-mp).*(lp-mp) ./ (2*lp.*(2*lp-1)) );
            % elseif mu==-1 && l==lp
            %     c = -sqrt( (lp-mp).*(lp+1+mp) ./ (2*lp.*(lp+1)) );
            % elseif mu==-1 && l==lp+1
            %     c = sqrt( (lp+2+mp).*(lp+1+mp) ./ (2*(lp+1).*(2*lp+3)) );
            % elseif mu==0 && l==lp-1
            %     c = sqrt( (lp-mp).*(lp+mp) ./ (lp.*(2*lp-1)) );
            % elseif mu==0 && l==lp
            %     c = -mp./sqrt(lp.*(lp+1));
            % elseif mu==0 && l==lp+1
            %     c = -sqrt( (lp+1-mp).*(lp+1+mp) ./ ((lp+1).*(2*lp+3)) );
            % elseif mu==+1 && l==lp-1
            %     c = sqrt( (lp-1+mp).*(lp+mp) ./ (2*lp.*(2*lp-1)) );
            % elseif mu==+1 && l==lp
            %     c = sqrt( (lp+mp).*(lp+1-mp) ./ (2*lp.*(lp+1)) );
            % elseif mu==+1 && l==lp+1
            %     c = sqrt( (lp+1-mp).*(lp+2-mp) ./ (2*(lp+1).*(2*lp+3)) );
            % else
            %     c = 0;
            % end
        end
        function k = K(mu,p,pp,l,lp,m,mp)
            % K Function K
            %
            % R = K(MU,P,Pp,L,Lp,M,Mp) calculate the function K for indices 
            %   MU, P, Pp, L, Lp, M, Mp multipoles. 
            %   This function is used in computing optical force.        
            %
            % See also Coefficients, TMatrix.

            k = sparse(zeros(size(mu)));
            
            ind = lp>=l-1 & lp<=l+1 & mu+mp==m;
            k(ind) = TMatrix.C1(lp(ind),l(ind),mu(ind),m(ind)).*TMatrix.O(p(ind),pp(ind),l(ind),lp(ind));
            % k = 16*pi^2*sqrt(3/(4*pi))*C1(l,lp,mu,mp)*1i^(l-lp)*O(p,pp,l,lp);
            
            % if lp>=l-1 && lp<=l+1 && mu+mp==m
            %     k = TMatrix.C1(lp,l,mu,m)*TMatrix.O(p,pp,l,lp);
            %     % k = 16*pi^2*sqrt(3/(4*pi))*C1(l,lp,mu,mp)*1i^(l-lp)*O(p,pp,l,lp);
            % else
            %     k = 0;
            % end
        end
        function r = I(theta,phi,p,pp,l,lp,m,mp)
            % I I Integral
            %
            % R = I(THETA,PHI,P,Pp,L,Lp,M,Mp) calculate the I Integral for specified 
            %   values of THETA [0 pi] and PHI [-pi pi] and multipole indices 
            %   P, Pp, L, Lp, M, Mp.
            %   This function is used in computing optical force.       
            %
            % See also TMatrix, Coefficients.

            sh = SpHarm(theta,phi);
            Y1m = sh.Y(1,-1);
            Y10 = sh.Y(1,0);
            Y1p = sh.Y(1,+1);

            r = sqrt(4*pi/3) * ( ...
                conj(Y1m).*TMatrix.K(-ones(size(p)),p,pp,l,lp,m,mp) ...
                + conj(Y10).*TMatrix.K(zeros(size(p)),p,pp,l,lp,m,mp) ...
                + conj(Y1p).*TMatrix.K(ones(size(p)),p,pp,l,lp,m,mp) ...
                );
            % r = 4*pi/3 * ( ...
            % 	conj(Y1m.Y)*1i^(lp-l)/(16*pi^2)*K(-1,p,pp,l,lp,m,mp) ...
            % 	+ conj(Y10.Y)*1i^(lp-l)/(16*pi^2)*K(0,p,pp,l,lp,m,mp) ...
            % 	+ conj(Y1p.Y)*1i^(lp-l)/(16*pi^2)*K(1,p,pp,l,lp,m,mp) ...
            %     );

            % sh = SpHarm(theta,phi);
            % Y1m = sh.Y(1,-1);
            % Y10 = sh.Y(1,0);
            % Y1p = sh.Y(1,+1);
            % 
            % r = sqrt(4*pi/3) * ( ...
            %     conj(Y1m)*TMatrix.K(-1,p,pp,l,lp,m,mp) ...
            %     + conj(Y10)*TMatrix.K(0,p,pp,l,lp,m,mp) ...
            %     + conj(Y1p)*TMatrix.K(1,p,pp,l,lp,m,mp) ...
            %     );
            % % r = 4*pi/3 * ( ...
            % % 	conj(Y1m.Y)*1i^(lp-l)/(16*pi^2)*K(-1,p,pp,l,lp,m,mp) ...
            % % 	+ conj(Y10.Y)*1i^(lp-l)/(16*pi^2)*K(0,p,pp,l,lp,m,mp) ...
            % % 	+ conj(Y1p.Y)*1i^(lp-l)/(16*pi^2)*K(1,p,pp,l,lp,m,mp) ...
            % %     );
        end
        function D = wigner(p,pp,l,lp,m,mp,lgf,alpha,beta,gamma)
            % WIGNER Wigner Rotation Matrix
            %
            % D = WIGNER(P,Pp,L,Lp,M,Mp,LGF,ALPHA,BETA,GAMMA) calculate the 
            %   Wigner rotation matrix for given value of indices
            %   P, Pp, L, Lp, M, Mp. 
            %   LGF is a logaritm of the factorial(for preallocation). 
            %   ALPHA is a rotation angle around the z axis,
            %   BETA is a rotation angle around the y axis and 
            %   GAMMA is a rotation angle around the z axis. 
            %
            % See also TMatrix, LGF
            
            D = exp(-1i*m*alpha).*d(p,pp,l,lp,m,mp,beta).*exp(-1i*mp*gamma);
            
            function w = d(p,pp,l,lp,m,mp,beta)
                % D Wigner matrix 
                %
                % R = D(P,Pp,L,Lp,M,Mp,BETA) calculate the D_MMp(BETA) function 
                %   inside the Wigner rotation matrix for given value 
                %   of indices P, Pp, L, Lp, M, Mp. 
                %   BETA is a rotation angle around the y axis.  
                %
                % See also TMatrix, LGF
            
                w = zeros(size(l));
                indices = find(l+m>=0 & l-m>=0 & l+mp>=0 & l-mp>=0 & l==lp & p==pp);
                w(indices) = (-1).^(m(indices)-mp(indices)) .* ...
                    sqrt( ...
                    exp(lgf.getfact(l(indices)+m(indices)) + ...
                    lgf.getfact(l(indices)-m(indices)) + ...
                    lgf.getfact(l(indices)+mp(indices)) + ...
                    lgf.getfact(l(indices)-mp(indices))) ...
                    )' .* ...
                    ssum(l(indices),m(indices),mp(indices),beta);
            end
            
            function S = ssum(l,m,mp,beta)
                % SSUM Summation for computation of Wigner rotation matrix
                %
                % S = SSUM(L,M,Mp,BETA) calculate the summation inside the wigner rotation matrix.
                %   L, M, Mp are multipole indices and BETA is a rotation angle around y axes.
                %   This summation performed for all integer values that give
                %   non-negative factorials.
                %
                % See also TMatrix.

                lpt = l;
                l = reshape(l,1,numel(l));
                m = reshape(m,1,numel(m));
                mp = reshape(mp,1,numel(mp));
                S = zeros(size(l));
                smin = max(min(mp-m),0);
                smax = min(max(l-m),max(l+mp));
                for s = smin:1:smax

                    indices = find(s>=0 & l+mp-s>=0 & m-mp+s>=0 & l-m-s>=0 );
                    S(indices) = S(indices) + ...
                        (-1)^s .* ...
                        cos(beta/2).^(2*l(indices)+mp(indices)-m(indices)-2*s) .* ...
                        sin(beta/2).^(m(indices)-mp(indices)+2*s) ./ ...
                        exp( ...
                        lgf.getfact((l(indices)+mp(indices)-s)) + ...
                        lgf.getfact(s) + ...
                        lgf.getfact((m(indices)-mp(indices)+s)) + ...
                        lgf.getfact((l(indices)-m(indices)-s) ...
                        ));
                end
                S = reshape(S,size(lpt));
            end

        end        
        function I = Gaunt(l,L,lp,m,M,mp,cg)
            % GAUNT Gaunt integrals
            %
            % I = GAUNT(L,LL,Lp,M,MM,Mp,CG) calculate the Gaunt integrals for 
            %   multipole indices L, LL, Lp, M, MM, Mp.
            %   CG is the Clebsch-Gordan coefficients (for preallocation).
            %
            % See also TMatrix, CG.

            I = sqrt((2*l+1).*(2*L+1)./(4*pi*(2*lp+1))).*cg.C(l,L,lp,m,M,mp).*cg.C(l,L,lp,zeros(size(l)),zeros(size(L)),zeros(size(lp)));
        end
        function g = Gj(lp,l,mp,m,kR,sh,cg,varargin)
            % GJ Nozawa's addition theorem for spherical Bessel harmonic
            %
            % G = GJ(Lp,L,Mp,M,kR,SH,CG) calculate the Gaunt addition theorem for 
            %   indices Lp, L, Mp, M of multipoles field.
            %   kR is a amount of translation times to medium wave number.
            %   SH is a spherical Bessel harmonics (for preallocation).
            %   CG is the Clebsch-Gordan coefficients (for preallocation).
            %
            % G = GJ(Lp,L,Mp,M,KR,SH,CG,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       Lmax      -   maximum number for multipole index LL 
            %
            % See also TMatrix, CG, SpHarm, SpBessel.

            % Lmax
            Lmax = cg.Jmax_created;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'lmax')
                    lmax = varargin{n+1};
                end
            end
            
            M = mp-m;

            g = sparse(zeros(size(l)));
            for L = 0:1:Lmax
                
                % conjshY = zeros(size(M));
                % for t = 1:1:numel(M)
                %     if abs(l(t)-L)<=lp(t) && lp(t)<=l(t)+L &&  abs(M(t))<=L 
                %         conjshY(t) = conj(sh.Y(L,M(t)));
                %     end
                % end
                % 
                % g = g + 4*pi*1i.^(lp+L-l).* ... 
                %     SpBessel.j(L,kR).* ...
                %     conjshY.* ...
                %     TMatrix.Gaunt(l,L*ones(size(l)),lp,m,M,mp,cg);

                % indices = find(abs(l-L)<=lp & lp<=l+L & abs(M)<=L & mod(l+lp+L,2)==0);
                indices = find(L>=abs(l-lp) & L<=l+lp & abs(M)<=L & mod(l+lp+L,2)==0);
                
                shY = zeros(size(indices));
                for t = 1:1:numel(indices)
                    shY(t) = sh.Y(L,M(indices(t)));
                end
                conjshY = conj(shY);
                
                g(indices) = g(indices) + 4*pi*1i.^(lp(indices)+L-l(indices)).* ... 
                    SpBessel.j(L,kR).* ...
                    conjshY.* ...
                    TMatrix.Gaunt(l(indices),L*ones(size(indices)),lp(indices),m(indices),M(indices),mp(indices),cg);
            end
            
        end
        function g = Gh(lp,l,mp,m,kR,sh,cg,varargin)
            % GH Nozawa's addition theorem for spherical Hankel harmonic
            %
            % G = GH(Lp,L,Mp,M,KR,SH,CG) calculate the Gaunt addition theorem for 
            %   indices Lp, L, Mp, M of multipoles field.
            %   KR is a amount of translation times to medium wave number.
            %   SH is a spherical Hankel harmonics (for preallocation). 
            %   CG is the Clebsch-Gordan coefficients (for preallocation).
            %
            % G = GH(Lp,L,Mp,M,KR,SH,CG,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       Lmax      -   maximum number for multipole index LL
            %
            % See also TMatrix, CG, SpHarm, SpBessel.
            
            % Lmax
            Lmax = cg.Jmax_created;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'lmax')
                    Lmax = varargin{n+1};
                end
            end
            
            M = mp-m;
            
            g = sparse(zeros(size(l)));
            for L = 0:1:Lmax
                
                % indices = find(abs(l-L)<=lp & lp<=l+L & abs(M)<=L & mod(l+lp+L,2)==0);
                indices = find(L>=abs(l-lp) & L<=l+lp & abs(M)<=L & mod(l+lp+L,2)==0);
                
                shY = zeros(size(indices));
                for t = 1:1:numel(indices)
                    shY(t) = sh.Y(L,M(indices(t)));
                end
                conjshY = conj(shY);
                
                
                g(indices) = g(indices) + 4*pi*1i.^(lp(indices)+L-l(indices)).* ...
                    SpBessel.h(L,kR).* ...
                    conjshY.* ...
                    TMatrix.Gaunt(l(indices),L*ones(size(indices)),lp(indices),m(indices),M(indices),mp(indices),cg);
            end
        end
        function J = JJ(lp,l,mp,m,R,k,cg,varargin)
            % JJ Translation matrix
            %
            % J = JJ(Lp,L,Mp,M,R,K,CG) calculate the translation matrix for 
            %   indices Lp, L, Mp, M of spherical Bessel harmonics. 
            %   R is a translated veector (vector-point),K is wave number in medium.
            %   CG is the Clebsch-Gordan coefficients (for preallocation).
            %
            % J = JJ(Lp,L,Mp,M,R,K,CG,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       Lmax      -   maximum number for multipole index LL 
            %
            % See also TMatrix, CG, Point.

            kR = k*norm(R);
            [theta,phi] = Transform.Car2Sph(R.X,R.Y,R.Z);
            sh = SpHarm(theta,phi);
            
            u = ones(size(l));
            
            J11 = sparse(zeros(size(l)));
            J12 = sparse(zeros(size(l)));
            for mu = -1:1:1
                
                % J11 = J11 + cg.C(u,lp,lp,-mu*u,mp+mu,mp).*cg.C(u,l,l,-mu*u,m+mu,m).*TMatrix.Gj(lp,l,mp+mu,m+mu,R,k,cg,varargin{:});
                % J12 = J12 + 1i*sqrt((2*lp+1)./(lp+1)).* cg.C(u,lp-1,lp,-mu*u,mp+mu,mp).*cg.C(u,l,l,-mu*u,m+mu,m).*TMatrix.Gj(lp-1,l,mp+mu,m+mu,R,k,cg,varargin{:});
                
                Ct = cg.C(u,l,l,-mu*u,m+mu,m);
                C11 = cg.C(u,lp,lp,-mu*u,mp+mu,mp) .* Ct;
                C12 = cg.C(u,lp-1,lp,-mu*u,mp+mu,mp) .* Ct;
                
                indices11 = find(C11~=0);
                indices12 = find(C12~=0);
                
                J11(indices11) = J11(indices11) + C11(indices11) .* TMatrix.Gj(lp(indices11),l(indices11),mp(indices11)+mu,m(indices11)+mu,kR,sh,cg,varargin{:});
                J12(indices12) = J12(indices12) + C12(indices12) .* TMatrix.Gj(lp(indices12)-1,l(indices12),mp(indices12)+mu,m(indices12)+mu,kR,sh,cg,varargin{:});
            end
            J12 = 1i*sqrt((2*lp+1)./(lp+1)) .* J12;
            
            % placing the J11 J12 J21 J22 elements inside J matrix
            J = sparse(zeros(2*size(l)));
            J(1:2:end,1:2:end) = J11; % J11
            J(2:2:end,2:2:end) = J11; % J22
            J(1:2:end,2:2:end) = J12; % J12 
            J(2:2:end,1:2:end) = J12; % J22 
            
        end
        function H = HH(lp,l,mp,m,R,k,cg,varargin)
            % HH Translation matrix
            %
            % H = HH(Lp,L,Mp,M,R,K,CG) calculate the translation matrix for 
            %   indices Lp, L, Mp, M of spherical Hankel harmonics. 
            %   R is a translated vector (vector-point).
            %   K is wave number in medium. 
            %   CG is the Clebsch-Gordan coefficients (for preallocation).
            %
            % H = HH(Lp,L,Mp,M,R,K,CG,'PropertyName',PropertyValue) 
            %   sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       Lmax      -   maximum number for multipole index LL 
            %
            % See also TMatrix, CG, Point.
            
            kR = k*norm(R);
            [theta,phi] = Transform.Car2Sph(R.X,R.Y,R.Z);
            sh = SpHarm(theta,phi);

            u = ones(size(l));
            
            H11 = sparse(zeros(size(l)));
            H12 = sparse(zeros(size(l)));
            if norm(R)~=0  % fi norm(R)==0, H=0
                for mu=-1:1:1
                    % H11 = H11 + cg.C(u,lp,lp,-mu*u,mp+mu,mp).*cg.C(u,l,l,-mu*u,m+mu,m).*TMatrix.Gh(lp,l,mp+mu,m+mu,R,k,cg,varargin);
                    % H12 = H12 + 1i*sqrt((2*lp+1)./(lp+1)).* cg.C(1*u,lp-1,lp,-mu*u,mp+mu,mp).*cg.C(u,l,l,-mu*u,m+mu,m).*TMatrix.Gh(lp-1,l,mp+mu,m+mu,R,k,cg,varargin);

                    Ct = cg.C(u,l,l,-mu*u,m+mu,m);
                    C11 = cg.C(u,lp,lp,-mu*u,mp+mu,mp) .* Ct;
                    C12 = cg.C(u,lp-1,lp,-mu*u,mp+mu,mp) .* Ct;

                    indices11 = find(C11~=0);
                    indices12 = find(C12~=0);

                    H11(indices11) = H11(indices11) + C11(indices11) .* TMatrix.Gh(lp(indices11),l(indices11),mp(indices11)+mu,m(indices11)+mu,kR,sh,cg,varargin{:});
                    H12(indices12) = H12(indices12) + C12(indices12) .* TMatrix.Gh(lp(indices12)-1,l(indices12),mp(indices12)+mu,m(indices12)+mu,kR,sh,cg,varargin{:});
                end
                H12 = 1i*sqrt((2*lp+1)./(lp+1)) .* H12;
            end
            
            % placing the H11 H12 H21 H22 elements inside H matrix
            H = sparse(zeros(2*size(l)));
            H(1:2:end,1:2:end) = H11; % H11
            H(2:2:end,2:2:end) = H11; % H22  
            H(1:2:end,2:2:end) = H12; % H12 
            H(2:2:end,1:2:end) = H12; % H21 
            
        end
    end
end