classdef BrownianMotion1DFreeInertial < BrownianMotion1DFree
	% BrownianMotion1DFreeInertial < BrownianMotion1DFree < BrownianMotion : 1D free Brownian motion
    %   Free Brownian motion in 1D taking into account the particle inertia.
    %   
    % BrownianMotion1DFreeInertial properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
    %   d       -   density [kg m^-3]
    %
    % BrownianMotion1DFreeInertial methods:
    %   BrownianMotion  -   constructor
    %   kBT             -   thermal energy [J] < BrownianMotion
    %   D               -   diffusion constant [m^2/s] < BrownianMotion
    %   gamma           -   friction coefficient [Kg/s] < BrownianMotion
    %   times           -   sample times [s] < BrownianMotion
    %   simulate        -   run simulaiton of Brownian motion < BrownianMotion
    %   plot            -   plot Brownian motion < BrownianMotion
    %   play            -   play Brownian motion < BrownianMotion
    %   dimensions      -   numebr of dimensions < BrownianMotion1DFree
    %   forsimulate     -   simulates the Brownian motion
    %   forplot         -   plots Brownian motion < BrownianMotion1DFree
    %   V               -   particle volume [m^3]
    %   M               -   particle mass [Kg]
    %
    % See also BrownianMotion, BrownianMotion1DFree.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        d   % density [kg m^-3]
    end
    methods
        function obj = BrownianMotion1DFreeInertial(dt,R,eta,T,d)
            % BROWNIANMOTION1DFREEINERTIAL(DT,R,ETA,T,D) consstructs a 1D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA,
            %   absolute temperature T and particle density D
            %
            % See also BrownianMotion, BrownianMotion1DFree.
            
            Check.isreal('d must be a positive real number',d,'>',0)
            Check.samesize('d must be a positive real number',d,0)

            obj = obj@BrownianMotion1DFree(dt,R,eta,T);
            
            obj.d = d;
        end
        function res = V(bm)
            % V Particle volume [m^3]
            %
            % v = V(BM) returns the particle volume [m^3].
            %
            % See also BrownianMotion1DFreeInertial.

            res = 4/3*pi*bm.R^3;
        end
        function res = m(bm)
            % M Particle mass [Kg]
            %
            % m = M(BM) returns the particle mass [Kg].
            %
            % See also BrownianMotion1DFreeInertial.

            res = bm.V*bm.d;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 1D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion1DFreeInertial.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real vector 2x1',r0)
            Check.samesize('r0 must be a real vector 2x1',r0,[0 0]')
            
            % pre-calculation coefficients
            dt = bm.dt;
            C1 = 1-bm.gamma/bm.m*dt;
            hs = sqrt(2*bm.kBT*bm.gamma*dt)/bm.m*h; % scaled noise
            
            % inizialization
            r = zeros(N,1);
            r(1) = r0(1);
            r(2) = r0(2);
            
            v = zeros(N,1);
            v(2) = r0(2)-r0(1);
            
            % simulation
            for n = 3:1:N
                r(n) = r(n-1) + v(n-1)*dt;
                v(n) = C1*v(n-1) + hs(n); 
            end
        end
    end
end