classdef PhaseMaskLG < PhaseMask
	% PhaseMaskLG < PhaseMask : Constant phase mask
    %   A phase mask with a phase singularity of order L.
    %
    % PhaseMaskLG properties:
	%   slm     -   spatial light modulator < PhaseMask 
	%   phase	-   phase shift values < PhaseMask 
    %   l       -   Laguerre-Gaussian beam order
    %   phase0  -   phase value [default = 0 rad]
    %
    % PhaseMaskLG methods:
    %   PhaseMaskLG         -   constructor
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
        l  % Laguerre-Gaussian beam order
        phase0  % phase value [default = 0 rad]
    end
    methods
        function pm = PhaseMaskLG(slm,l,varargin)
            % PHASEMASKLG(SLM,L) constructs phase mask with a phase
            %   singularity of order L.
            % 
            % PHASEMASKLG(SLM,L,'Phase0',PHASE0) adds a
            %   constant phase offset equal to PHASE0.
            %
            % See also PhaseMaskLG.

            pm = pm@PhaseMask(slm);
                        
            % grating phase offset [default = 0 rad]
            phase0 = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'phase0')
                    phase0 = varargin{n+1};
                end
            end

            pm.l = l;
            pm.phase0 = phase0;
            
            [X,Y] = pm.slm.pmeshgrid();  % SLM pixel positions [m]
            pm.phase = mod(phase0+l*atan2(Y,X),2*pi);
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskLG.

            disp(['<a href="matlab:help PhaseMaskLG">PhaseMaskLG</a> l=' int2str(pm.l) ' phase0=' num2str(pm.phase0/(2*pi)*360) 'deg']);
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskLG.

            txt = ['PhaseMaskLG - $l=' int2str(pm.l) '$ $\varphi_0=' num2str(pm.phase0/(2*pi)*360) '^{\circ}$'];
        end
    end
end