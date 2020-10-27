classdef PhaseMaskFresnel < PhaseMask
	% PhaseMaskFresnel < PhaseMask : Constant phase mask
    %   A phase mask representing a Fresnel lens with focal length F.
    %
    % PhaseMaskFresnel properties:
	%   slm     -   spatial light modulator < PhaseMask 
	%   phase	-   phase shift values < PhaseMask 
    %   f       -   focal length
    %   phase0  -   phase value [default = 0 rad]
    %
    % PhaseMaskFresnel methods:
    %   PhaseMaskFresnel    -   constructor
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
        f  % focal length
        phase0  % phase value [default = 0 rad]
    end
    methods
        function pm = PhaseMaskFresnel(slm,f,varargin)
            % PHASEMASKFRESNEL(SLM,F) constructs a Fresnel lens with focal length F.
            % 
            % PHASEMASKFRESNEL(SLM,F,'Phase0',PHASE0) adds a
            %   constant phase offset equal to PHASE0.
            %
            % See also PhaseMaskFresnel.

            pm = pm@PhaseMask(slm);

            % grating phase offset [default = 0 rad]
            phase0 = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'phase0')
                    phase0 = varargin{n+1};
                end
            end

            [X,Y] = pm.slm.pmeshgrid();  % SLM pixel positions [m]
            r = hypot(X,Y); % [m]
            pm.phase = mod(2*pi*r.^2/(2*slm.lambda*f),2*pi);
            pm.phase = mod(pm.phase + phase0, 2*pi);
            pm.f = f;
            pm.phase0 = phase0;
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskFresnel.

            disp(['<a href="matlab:help PhaseMaskFresnel">PhaseMaskFresnel</a> alpha=' num2str(pm.f) 'm']);
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskFresnel.

            txt = ['PhaseMaskFresnel - $f=' num2str(pm.f) '{\rm m}$'];
        end
    end
end