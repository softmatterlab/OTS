classdef PhaseMaskRS < PhaseMask
	% PhaseMaskRS < PhaseMask : Constant phase mask
    %   A phase mask obtained from the random superposition of other phase masks.
    %
    % PhaseMaskRS properties:
	%   slm     -   spatial light modulator < PhaseMask 
	%   phase	-   phase shift values < PhaseMask 
	%   p       -   trap probabilities
    %
    % PhaseMaskRS methods:
    %   PhaseMaskRS         -   constructor
    %   disp                -   displays phase mask name
    %   title               -   phase mask title
    %   plot                -   plots phase mask < PhaseMask 
    %   uplus               -   returns unchanged phase mask (+PM) < PhaseMask 
    %   uminus              -   inverses phase matrix (-PM) < PhaseMask 
    %   plus                -   sums two phase masks (PM1+PM2) < PhaseMask 
    %   minus               -   substracts two phase masks (PM1-PM2) < PhaseMask 
    %   grating             -   adds a grating to the phase mask < PhaseMask 
    %   fresnel             -   adds a Fresnel lens to the phase mask < PhaseMask 
    %   focus               -   calculates the focal field < PhaseMask 
    %
    % See also PhaseMask, SLM, Holography.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        p  % trap probabilities
    end
    methods
        function pm = PhaseMaskRS(slm,pms,varargin)
            % PHASEMASKRS(SLM,PMs) constructs a phase mask by random
            %   superpositon of the phase masks PMs.
            % 
            % PHASEMASKRS(SLM,PMs,'P',P) weights the masks PMs according to
            %   probabilities P.
            %
            % See also PhaseMaskRS.
            
            pm = pm@PhaseMask(slm);
            
            T = length(pms);  % number of traps

            % probabilities
            pm.p = ones(1,T);
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'p')
                    pm.p = varargin{n+1};
                end
            end
            
            psum = cumsum(pm.p);
            r = rand(slm.N,slm.M)*psum(end);
            pm.phase = pms{1}.phase;
            for t = 2:1:T
                pm.phase(r>psum(t-1)) = pms{t}.phase(r>psum(t-1));
            end
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskRS.

            disp('<a href="matlab:help PhaseMaskRS">PhaseMaskRS</a>');
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskRS.

            txt = ['PhaseMaskRS'];
        end
    end
end