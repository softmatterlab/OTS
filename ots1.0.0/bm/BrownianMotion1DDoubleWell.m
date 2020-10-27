classdef BrownianMotion1DDoubleWell < BrownianMotion1DFree
	% BrownianMotion1DDoubleWell < BrownianMotion1DFree < BrownianMotion : 1D free Brownian motion
    %   Brownian motion in a 1D double-well potential.
    %   
    % BrownianMotion1DDoubleWell properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
    %   a       -   a parameter of double-well potential [N m^-3]
    %   b       -   b parameter of double-well potential [N m^-1]
    %
    % BrownianMotion1DDoubleWell methods:
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
    %
    % See also BrownianMotion.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        a	% [N m^-3]
        b   % [N m^-1]
    end
    methods
        function obj = BrownianMotion1DDoubleWell(dt,R,eta,T,a,b)
            % BROWNIANMOTION1DDOUBLEWELL(DT,R,ETA,T,A,B) consstructs a 1D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %   A and B characterize the doulbe-well potential.
            %
            % See also BrownianMotion1DDoubleWell, BrownianMotion, BrownianMotion1DFree.
            
            Check.isreal('a must be a positive real number',a,'>',0)
            Check.samesize('a must be a positive real number',a,0)
            Check.isreal('b must be a real number',b)
            Check.samesize('b must be a positive real',b,0)

            obj = obj@BrownianMotion1DFree(dt,R,eta,T);
            
            obj.a = a;
            obj.b = b;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 1D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion1DDoubleWell.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real number',r0)
            Check.samesize('r0 must be a real number',r0,[0])
            
            % pre-calculation coefficients
            C1 = 1+ bm.b/bm.gamma*bm.dt;
            C2 = bm.a/bm.gamma*bm.dt;
            hs = sqrt(2*bm.D*bm.dt)*h; % scaled noise

            % inizialization
            r = zeros(N,1);
            r(1) = r0;
            
            % simulation
            for n = 2:1:N
                r(n) = C1*r(n-1) - C2*r(n-1)^3 + hs(n);
            end
        end
    end
end