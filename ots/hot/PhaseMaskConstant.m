classdef PhaseMaskConstant < PhaseMask
	% PhaseMaskConstant < PhaseMask : Constant phase mask
    %   A phase mask were all pixels have the same value.
    %
    % PhaseMaskConstant properties:
	%   slm     -   spatial light modulator < PhaseMask 
	%   phase	-   phase shift values < PhaseMask 
    %   phase0  -   phase value [default = 0 rad]
    %
    % PhaseMaskConstant methods:
    %   PhaseMaskConstant   -   constructor
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
        phase0  % phase value [default = 0 rad]
    end
    methods
        function pm = PhaseMaskConstant(slm,varargin)
            % PHASEMASKCONSTANT(SLM) constructs an phase mask with constant phase equal 0.
            % 
            % PHASEMASKCONSTANT(SLM,'Phase0',PHASE0) constructs an PhaseMask 
            %   with constant phase equal PHASE0.
            %
            % See also PhaseMaskConstant.

            pm = pm@PhaseMask(slm);
                       
            % phase value [default = 0 rad]
            phase0 = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'phase0')
                    phase0 = mod(varargin{n+1},2*pi);
                end
            end

            pm.phase = phase0*ones(slm.N,slm.M);
            pm.phase0 = phase0;
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskConstant.

            disp(['<a href="matlab:help PhaseMaskConstant">PhaseMaskConstant</a> phase0=' num2str(pm.phase0/(2*pi)*360) 'deg']);
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskConstant.

            txt = ['PhaseMaskConstant - $\varphi_0=' num2str(pm.phase0/(2*pi)*360) '^{\circ}$'];
        end
    end
end