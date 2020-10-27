classdef TMatrixInclusions < TMatrix
    % TMatrixInclusions < TMatrix : transition matrix for inclusion
    %   A transition matrix for inclusion relates the coefficients of incident and
    %   scattered field.
    %
    % TMatrixInclusions properties:
    %   sphere         -   mie for host sphere 
    %   inclusions     -   mie for inclusions (mie, cell)
    %   R              -   inclusions positions (Point, cell)
    %   Lmax           -   number of coefficients
    %   cg             -   Clebsch-Gordan coefficients (for prealocation)
    %
    % TMatrixInclusion methods:
    %   TMatrixInclusions -   constructor
    %   Li                -   incident field amplitude coefficient number < TMatrix
    %   Ls                -   scattered field amplitude coefficient number < TMatrix
    %   As                -   scattered field coefficients < TMatrix
    %   force             -   optical force < TMatrix
    %   torque            -   optical torque < TMatrix
    %   incoming          -   incoming field < TMatrix
    %   scattering        -   scattered field < TMatrix
    %   total             -   total field < TMatrix
    %   rotate            -   T-matrix rotation < TMatrix
    %   translate         -   T-matrix translation < TMatrix
    %   k0                -   medium wavenumber
    %   nm                -   medium refraction index
    %
    % TMatrixInclusion static methods:
    %   index         -   index coressponding to plm < TMatrix
    %   indices       -   indices coressponding to plm and p'l'm' < TMatrix
    %   plms          -   plm and p'l'm' corresponding to indices < TMatrix
    %   O             -   function O (for force calculation) < TMatrix
    %   C1            -   Clebsch-Gordan coefficients for l1 = 1 < TMatrix
    %   K             -   function K (for force calculation) < TMatrix
    %   I             -   function I (for force calculation) < TMatrix
    %   wigner        -   Wigner rotation matrix < TMatrix
    %   Gaunt         -   Gaunt integral < TMatrix
    %   Gj            -   Nozawa's addition theorem for Bessel harmonic < TMatrix
    %   Gh            -   Nozawa's addition theorem for Hankel harmonic < TMatrix
    %   JJ            -   translation matrix coresponding to Gj < TMatrix
    %   HH            -   translation matrix coresponding to Gh < TMatrix
    %
    % See also TMatrix, TMatrixCluster, MieParticle.
    
    %   Author: S. Masoumeh Mousavi
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    
    properties
        sphere  % host sphere
        inclusions  % Mie particles inside host sphere (TMatrixSphere, vector)
        R  % particle positions related to host sphere (Point, vector)
        Lmax  % number of coefficients
        cg  % Clebsch-Gordan coefficients
    end
    
    methods
        function obj = TMatrixInclusions(sphere,inclusions,Lmax,R,cg)
            % TMatrixInclusions Inclusions T-matrix
            %
            % TM = TMatrixInclusions(SPHERE,INCLUSIONS,LMAX,R,CG) constructs a
            %   t-matrix for inclusions ,SPHERE is mie constructer for spherical host sphere.
            %   INCLUSIONS is a mie constructed for spherical inclusions
            %   with different radius and refraction index (cell). LMAX is maxium number of mie
            %   coefficients between INCLUSIONS. R is cell that shows the position of the center of mass for
            %   spherical inclusion related to the center of mass for host
            %   sphere.CG is prealocated Clebesh-Gordan coefficints.
            %
            % See also TMatrixInclusion, TMatrix, MieParticle, CG.
            
            obj.sphere = sphere; % mie for host sphere
            Rho0 = obj.sphere.R();
            obj.inclusions = inclusions; % mie for inclusions
            obj.R = R;
            obj.cg = cg;
            N = length(inclusions);  % number of inclusions
            obj.Lmax = Lmax;  % maximum number of mie coefficient between inclusions
            
            % distance between inclusions
            for a = 1:1:N
                for ap = 1:1:N
                    dR{a,ap} = R{a}-R{ap};
                end
            end
            
            % T-matrix for inclusions
            for a = 1:1:N
                tm_inclusions{a} = TMatrixSphere(inclusions{a},'Li',Lmax,'Ls',Lmax);  % t-matrix for sphere a
            end
            
            km = obj.k0*obj.nm;  % wavenumber in medium
            kp = obj.k0()*sphere.np(); % wave number for host sphere (as a medium for inclusions)
            nm = obj.nm;
            np = sphere.np();
            
            % T-matrix calculation
            dn = TMatrix.index(2,Lmax,Lmax)-2;  % Length of a T-matrix block
            R0 = sparse(zeros(dn));
            RW = sparse(zeros(dn));
            MW = sparse(zeros(dn));
            M0 = sparse(zeros(dn));
            
            for l = 1:1:Lmax % l==lp
                for m = -l:1:l % m==mp
                    % external medium
                    um = SpBessel.u(l,km*Rho0);
                    dum = SpBessel.du(l,km*Rho0);
                    wm = SpBessel.w(l,km*Rho0);
                    dwm = SpBessel.dw(l,km*Rho0);
                    
                    % host sphere
                    u = SpBessel.u(l,kp*Rho0);
                    du = SpBessel.du(l,kp*Rho0);
                    w = SpBessel.w(l,kp*Rho0);
                    dw = SpBessel.dw(l,kp*Rho0);
                    
                    R0(TMatrix.index(1,l,m)-2,TMatrix.index(1,l,m)-2) =  ( 1i*np )/( nm*u*dwm - np*du*wm );  % p==pp==1
                    R0(TMatrix.index(2,l,m)-2,TMatrix.index(2,l,m)-2) = ( 1i*np )/( np*u*dwm - nm*wm*du ); % p==pp==2
                    RW(TMatrix.index(1,l,m)-2,TMatrix.index(1,l,m)-2) = ( nm*w*dwm - np*wm*dw ) / ( 1i*np ); % p==pp==1
                    RW(TMatrix.index(2,l,m)-2,TMatrix.index(2,l,m)-2) = ( np*w*dwm - nm*wm*dw ) / ( 1i*np ); % p==pp==2
                    MW(TMatrix.index(1,l,m)-2,TMatrix.index(1,l,m)-2) = ( nm*w*dum - np*um*dw ) / ( -1i*np ); % p==pp==1
                    MW(TMatrix.index(2,l,m)-2,TMatrix.index(2,l,m)-2) = ( np*w*dum - nm*um*dw ) / ( -1i*np ); % p==pp==2
                    M0(TMatrix.index(1,l,m)-2,TMatrix.index(1,l,m)-2) =  ( nm*u*dum - np*du*um )/( -1i*np ) ; % p==pp==1
                    M0(TMatrix.index(2,l,m)-2,TMatrix.index(2,l,m)-2) = ( np*u*dum - nm*du*um ) /( -1i*np );% p==pp==2
                end
            end
            
            ni = [0:dn:N*dn];  % block initial coordinates
            
            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:dn+2],[1:1:dn+2]);
            [l,lp] = meshgrid(l(1:2:end),lp(1:2:end));
            [m,mp] = meshgrid(m(1:2:end),mp(1:2:end));
            
            invJ0 = sparse(zeros(dn, N*dn));
            J0 = sparse(zeros(N*dn, dn));
            
            for a = 1:1:N
                J1 = TMatrix.JJ(lp,l,mp,m,-R{a},km,cg);
                J2 = TMatrix.JJ(lp,l,mp,m,R{a},km,cg);
                invJ0(:, ni(a)+1:1:a*dn) = J1(3:end,3:end);
                J0(ni(a)+1:1:a*dn, :) = J2(3:end,3:end);
            end
            % calculation of M matrix in inclusion
            M11 = sparse(zeros(dn*N));
            dM = mat2cell(M11, dn*ones(N,1), dn*ones(N,1));
            for a = 1:1:N
                for ap = 1:1:N
                    if a==ap
                        dM{a,a} = inv(-tm_inclusions{a}.T(3:end,3:end)); % -tm_mie(a).T(3:end,3:end).^-1;
                    else
                        H1 = TMatrix.HH(lp,l,mp,m,dR{a,ap},km,cg);
                        dM{a,ap} = H1(3:end,3:end);
                    end
                end
            end
            M11 = cell2mat(dM);
            M21 = RW*invJ0;
            M12 = J0;
            M22 = inv(R0);
            M = sparse(zeros((N+1)*dn));
            M(1:1:N*dn,1:1:N*dn) = M11;
            M(N*dn+1:1:(N+1)*dn,1:1:N*dn) = M21;
            M(1:1:N*dn,N*dn+1:1:(N+1)*dn) = M12;
            M(N*dn+1:1:(N+1)*dn,N*dn+1:1:(N+1)*dn) = M22;
            % calculation inverse of M matrix
            [Lm,Um,pm] = lu(M); % LU decomposition
            invM = inv(Um)*inv(Lm)*pm;
            % elements of M neede for calculation T-matrix
            invM12 = invM(1:1:N*dn,N*dn+1:1:(N+1)*dn);
            invM22 = invM(N*dn+1:1:(N+1)*dn,N*dn+1:1:(N+1)*dn);
            % T-matrix
            T = M0*invM22+MW*invJ0*invM12;
            obj.T = sparse(zeros(dn+2));
            obj.T(3:end,3:end) = T;
            
        end
        function res = k0(tm)
            % K0 Medium wave number
            %
            % RES = K0(TM) determines the medium wave number.
            %
            % See also TMatrixInclusion, TMatrix, MieParticle.

            res = tm.sphere.k0;
        end
        function res = nm(tm)
            % NM Medium refractive index
            %
            % RES = NM(TM) determines the refractive index of medium.
            %
            % See also TMatrixInclusion, TMatrix, MieParticle.

            res = tm.sphere.nm;
        end
        
    end
end
