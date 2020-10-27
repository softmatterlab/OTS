classdef BrownianMotion1DFree < BrownianMotion
	% BrownianMotion1DFree < BrownianMotion : 1D free Brownian motion
    %   Free Brownian motion in 1D.
    %   
    % BrownianMotion1DFree properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
    %
    % BrownianMotion1DFree methods:
    %   BrownianMotion  -   constructor
    %   kBT             -   thermal energy [J] < BrownianMotion
    %   D               -   diffusion constant [m^2/s] < BrownianMotion
    %   gamma           -   friction coefficient [Kg/s] < BrownianMotion
    %   times           -   sample times [s] < BrownianMotion
    %   simulate        -   run simulaiton of Brownian motion < BrownianMotion
    %   plot            -   plot Brownian motion < BrownianMotion
    %   play            -   play Brownian motion < BrownianMotion
    %   dimensions      -   numebr of dimensions
    %   forsimulate     -   simulates the Brownian motion
    %   forplot         -   plots Brownian motion
    %
    % See also BrownianMotion.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    methods
        function obj = BrownianMotion1DFree(dt,R,eta,T)
            % BROWNIANMOTION1DFREE(DT,R,ETA,T) consstructs a 1D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %
            % See also BrownianMotion, BrownianMotion1DFree.
            
            obj = obj@BrownianMotion(dt,R,eta,T);
        end
        function M = dimensions(bm)
            % DIMENSIONS Number of dimensions (1)
            %
            % M = DIMENSIONS(BM) returns the number of dimensions of the
            %   Brownian motion (1).
            %
            % See also BrownianMotion1DFree.

            M = 1;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 1D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion1DFree.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real number',r0)
            Check.samesize('r0 must be a real number',r0,[0])
            
            % pre-calculation coefficients
            hs = sqrt(2*bm.D*bm.dt)*h; % scaled noise
            
            % inizialization
            r = zeros(N,1);
            r(1) = r0;
            
            % simulation
            for n = 2:1:N
                r(n) = r(n-1) + hs(n);
            end
        end
        function forplot(bm,t,r,dimension)
            % FORPLOT Plots 1D Brownian motion
            %
            % FORPLOT(BM,T,R,UNITS) plots the Brownian motion (T,X) in a
            %   plot time versus position. 
            %   UNITS represents the units of the Brownian motion.
            %
            % See also BrownianMotion1D.

            hold on
            plot(t,r)
            xlabel('time [s]')
            ylabel(['x ' dimension ])
            axis tight
        end
    end
end