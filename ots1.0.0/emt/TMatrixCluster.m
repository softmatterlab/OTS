classdef TMatrixCluster < TMatrix
    % TMatrixCluster < TMatrix : Transition matrix for cluster
    %   A transition matrix for cluster relates the coefficients of incident and
    %   scattered field.
    %
    % TMatrixCluster properties:
    %   mieparticles   -   Mie particles (MieParticle, cell)
    %   R              -   spheres positions (Point, cell)
    %   N              -   number of spheres
    %   dR             -   distance between particles (Point, cell)
    %   Lmax           -   number of coefficients
    %   cg             -   Clebsch-Gordan coefficients (for prealocation)
    %   Mi             -   inverse of M matrix
    %
    % TMatrixCluster methods:
    %   TMatrixCluster    -   constructor
    %   Minverse          -   inverse of M matrix (matrix)
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
    % TMatrixCluster static methods:
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
    % See also TMatrix, TMatrixSphere, TMatrixInclusions.
    
    %   Author: S. Masoumeh Mousavi
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        mieparticles  % Mie particles (MieParticle, vector)
        R  % particle positions (Point, vector)
        N  % number of particles
        dR  % distance between particles (Point, matrix)
        Lmax  % number of coefficients
        cg  % Clebsch-Gordan coefficients
        Mi  % inverse of M-matrix
    end
    methods
        function obj = TMatrixCluster(mieparticles,Lmax,R,cg,Mi)
            % TMATRIXCLUSTER Cluster T-matrix 
            %
            % TM = TMATRIXCLUSTER(MIEPARTICLES,LMAX,R,CG,MI) constructs a t-matrix for aggregates of spheres,
            %   MIEPARTICLES is a cell constructed from mie particle for spheres with different
            %   radius and refraction index. LMAX is maxium number of mie
            %   coefficients between mieparticles.
            %   R is cell that shows the position of the center of mass for each
            %   sphere in aggregates,CG is prealocated Clebesh-Gordan coefficints.
            %   MI is inverse of M matrix (optional).
            %
            % See also TMatrixCluster, TMatrix, Mie Particle, CG.
            
            obj.mieparticles = mieparticles;
            obj.Lmax = Lmax;
            obj.R = R;
            obj.cg = cg;
            N = length(mieparticles);  % number of particles
            obj.N = N;
            
            % distance between the center of particles to each other
            for a = 1:1:N
                for ap = 1:1:N
                    dR{a,ap} = R{a}-R{ap};
                end
            end
            obj.dR = dR;
            
            % wavenumber in medium
            km = obj.k0*obj.nm;
            
            % T-matrix for spheres
            for a = 1:1:N
                tm_mie{a} = TMatrixSphere(obj.mieparticles{a},'Li',Lmax,'Ls',Lmax);
            end
            
            % T-matrix calculation
            dn = TMatrix.index(2,Lmax,Lmax)-2;  % Length of a T-matrix block
            ni = [0:dn:N*dn];  % block initial coordinates

            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:dn+2],[1:1:dn+2]);
            [l,lp] = meshgrid(l(1:2:end),lp(1:2:end));
            [m,mp] = meshgrid(m(1:2:end),mp(1:2:end));

            Jpre = sparse(zeros(dn, N*dn));
            Jpost = sparse(zeros(N*dn, dn));
            for a = 1:1:N
                J1 = TMatrix.JJ(lp,l,mp,m,-R{a},km,cg);
                J2 = TMatrix.JJ(lp,l,mp,m,R{a},km,cg);
                Jpre(:, ni(a)+1:1:a*dn) = J1(3:end,3:end);
                Jpost(ni(a)+1:1:a*dn, :) = J2(3:end,3:end);
            end
            if nargin<5
                obj.Mi = obj.Minverse(tm_mie,dn,dR);
            else
                obj.Mi = Mi;
            end
            T = -Jpre*obj.Mi*Jpost;
            obj.T = sparse(zeros(dn+2));
            obj.T(3:end,3:end) = T;
        end
        function res = k0(tm)
            % K0 Medium wavenumber
            %
            % RES = K0(TM) determines the medium wave number.
            %
            % See also TMatrixCluster, TMatrix, MieParticle.

            res = tm.mieparticles{1}.k0;
        end
        function res = nm(tm)
            % NM Medium refractive index
            %
            % RES = NM(TM) determines the refractive index of medium.
            %
            % See also TMatrixCluster, TMatrix, MieParticle.

            res = tm.mieparticles{1}.nm;
        end
        function Mi = Minverse(tm,tm_mie,dn,dR)
            % MINVERSE Inverse of function M
            %
            % Mi = MINVERSE(TM,TM_MIE,Dn,DR) is a matrix which relates the amplitudes of the incident
            %   field to those of the fields scattered by each sphere in the aggregate.
            %   TM_MIE is a cell constructed by T-matrix for spheres in cluster
            %   Dn is Length of a T-matrix block and DR is a distance between
            %   spheres in cluster(Point,cell).
            %
            % See also TMatrixCluster, TMatrix.
            
            N = tm.N;  % number of particles
            
            k = tm.k0*tm.nm;  % wavenumber in medium

            % Calculation of M
            [p,pp,l,lp,m,mp] = TMatrix.plms([1:1:dn+2],[1:1:dn+2]);
            [l,lp] = meshgrid(l(1:2:end),lp(1:2:end));
            [m,mp] = meshgrid(m(1:2:end),mp(1:2:end));

            M = sparse(zeros(dn*N));
            dM = mat2cell(M, dn*ones(N,1), dn*ones(N,1));
            for a = 1:1:N
                for ap = 1:1:N
                    if a==ap
                        dM{a,a} = inv(-tm_mie{a}.T(3:end,3:end)); % -tm_mie(a).T(3:end,3:end).^-1;
                    else
                        H1 = TMatrix.HH(lp,l,mp,m,dR{a,ap},k,tm.cg);
                        dM{a,ap} = H1(3:end,3:end);
                    end
                end
            end
            M = cell2mat(dM);
            
            % inverse calculation
            [Lm,Um,p] = lu(M);
            Mi = inv(Um)*inv(Lm)*p;
            
        end
    end
end