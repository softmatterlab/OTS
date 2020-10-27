classdef PhaseMaskGS < PhaseMask
	% PhaseMaskGS < PhaseMask : Constant phase mask
    %   A phase mask with random phase values.
    %
    % PhaseMaskGS properties:
	%   slm         -   spatial light modulator < PhaseMask 
	%   phase       -   phase shift values < PhaseMask 
	%   alpha_vec   -   polar deflection angles < PhaseMaskSuperposition
	%   phi_vec     -   azimuthal deflection angles < PhaseMaskSuperposition
	%   f_vec       -   focal lengths < PhaseMaskSuperposition
	%   phase0_vec  -   initial phase values [default = 0 rad] < PhaseMaskSuperposition
    %
    % PhaseMaskGS methods:
    %   PhaseMaskGS  -   constructor
    %   disp                    -   displays phase mask name
    %   title                   -   phase mask title
    %   plot                    -   plots phase mask < PhaseMask 
    %   uplus                   -   returns unchanged phase mask (+PM) < PhaseMask 
    %   uminus                  -   inverses phase matrix (-PM) < PhaseMask 
    %   plus                    -   sums two phase masks (PM1+PM2) < PhaseMask 
    %   minus                   -   substracts two phase masks (PM1-PM2) < PhaseMask 
    %   grating                 -   adds a grating to the phase mask < PhaseMask 
    %   fresnel                 -   adds a Fresnel lens to the phase mask < PhaseMask 
    %   focus                   -   calculates the focal field < PhaseMask 
    %   gratinglengths          -   grating lengths
    %
    % See also PhaseMask, PhaseMaskSuperposition, SLM, Holography.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
    end
    methods
        function pm = PhaseMaskGS(slm,alpha_vec,phi_vec,f_vec,MaxIterations,varargin)
            % PHASEMASKGRATING(SLM,ALPHAvec,PHIvec,Fvec,MaxIterations) constructs a phase
            %   mask representing the superpositon of gratings (ALPHAvec and PHIvec) 
            %   and Fresnel lenses (Fvec). 
            %   MaxIteration is the maximum number of iterations of the
            %   Gerchberg-Saxton algorithm.
            % 
            % PHASEMASKGRATING(SLM,ALPHAvec,PHIvec,Fvec,MaxIterations,'Phase0',PHASE0vec) adds a
            %   constant phase offset equal to PHASE0vec to each of the superposition elements
            %   in the intial guess.
            %
            % See also PhaseMaskGS.

            pm = pm@PhaseMaskSuperposition(slm,alpha_vec,phi_vec,f_vec,varargin{:});
                        
            
            for k = 1:1:MaxIterations
                disp([' iteration ' num2str(k) '/' num2str(MaxIterations)]);
                phi = angle(sum(exp(1i.*deltatrap).*Vtraps it ...
                    ./sqrt(Vtraps it.*conj(Vtraps it)),3));
                for j = 1:T
                    Vtraps_it(:,:,j) = 1/N/M*sum(sum(exp(1i ...
                        *(phi-deltatrap(:,:,j))),1),2)*ones(M,N);
                end
            end
            phi = mod(phi,2*pi);
            pm.phase = round(phi/(2*pi/P))*2*pi/P;
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskGS.

            disp(['<a href="matlab:help PhaseMaskGS">PhaseMaskGS</a> ' num2str(numel(pm.alpha_vec)) ' traps']);
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskGS.

            txt = ['PhaseMaskGS - ' num2str(numel(pm.alpha_vec)) ' traps'];
        end
        function gl_vec = gratinglengths(pm)
            % GRATINGLENGTHS Grating lenghts
            %
            % GLvec = GRATINGLENGTHS(PM) returns the grating lengths.
            %
            % See also PhaseMaskGS.

            gl_vec = pm.slm.lambda./tan(pm.alpha_vec);
        end
    end
end