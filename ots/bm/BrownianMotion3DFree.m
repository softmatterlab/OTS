classdef BrownianMotion3DFree < BrownianMotion
	% BrownianMotion3DFree < BrownianMotion : 3D free Brownian motion
    %   Free Brownian motion in 3D.
    %   
    % BrownianMotion3DFree properties:
    %   dt      -   time step [s] < BrownianMotion
    %   R       -   particle radium [m] < BrownianMotion
    %   eta     -   medium bulk viscosity [Pa s] < BrownianMotion
    %   T       -   temperature [K] < BrownianMotion
    %   t       -   time [s] < BrownianMotion
    %   r       -   trajectory [m] < BrownianMotion
    %   h       -   noise < BrownianMotion
    %
    % BrownianMotion3DFree methods:
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
        function obj = BrownianMotion3DFree(dt,R,eta,T)
            % BROWNIANMOTION3DFREE(DT,R,ETA,T) consstructs a 3D Brownian motion with
            %   time step DT, particle radius R, fluid viscosity ETA and
            %   absolute temperature T.
            %
            % See also BrownianMotion, BrownianMotion3DFree.
            
            obj = obj@BrownianMotion(dt,R,eta,T);
        end
        function M = dimensions(bm)
            % DIMENSIONS Number of dimensions (3)
            %
            % M = DIMENSIONS(BM) returns the number of dimensions of the
            %   Brownian motion (3).
            %
            % See also BrownianMotion3DFree.

            M = 3;
        end
        function r = forsimulate(bm,N,r0,h)
            % FORSIMULATE Simulates 3D Brownian motion
            %
            % R = FORSIMULATE(BM,N,R0,H) simulates the Brownian motion 
            %   starting at R0 for N time steps and using noise H.
            %
            % See also BrownianMotion3DFree.
            
            Check.isinteger('N must be a positive integer',N,'>',0)
            Check.isreal('r0 must be a real vector 1x2',r0)
            Check.samesize('r0 must be a real vector 1x2',r0,[0 0 0])
            
            % pre-calculation coefficients
            hs = sqrt(2*bm.D*bm.dt)*h; % scaled noise
            
            % inizialization
            r = zeros(N,3);
            r(1,:) = r0;
            
            % simulation
            for n = 2:1:N
                r(n,:) = r(n-1,:) + hs(n,:);
            end
        end
        function forplot(bm,t,r,dimension)
            % FORPLOT Plots 3D Brownian motion
            %
            % FORPLOT(BM,T,R,UNITS) plots the Brownian motion R in a xyz-plot.
            %   UNITS represents the units of the Brownian motion.
            %
            % See also BrownianMotion3D.
                        
            hold on
            plot3(r(:,1),r(:,2),r(:,3))
            xlabel(['x ' dimension ])
            ylabel(['y ' dimension ])
            zlabel(['z ' dimension ])
            axis equal tight
            view(3)

        end
    end
end