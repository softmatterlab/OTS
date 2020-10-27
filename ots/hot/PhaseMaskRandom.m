classdef PhaseMaskRandom < PhaseMask
	% PhaseMaskRandom < PhaseMask : Constant phase mask
    %   A phase mask with random phase values.
    %
    % PhaseMaskRandom properties:
	%   slm     -   spatial light modulator < PhaseMask 
	%   phase	-   phase shift values < PhaseMask 
    %
    % PhaseMaskRandom methods:
    %   PhaseMaskRandom     -   constructor
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

    methods
        function pm = PhaseMaskRandom(slm)
            % PHASEMASKRANDOM(SLM) constructs a random phase mask.
            %
            % See also PhaseMaskRandom.

            pm = pm@PhaseMask(slm);

            pm.phase = 2*pi*rand(slm.N,slm.M);
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskRandom.

            disp('<a href="matlab:help PhaseMaskRandom">PhaseMaskRandom</a>');
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskRandom.

            txt = ['PhaseMaskRandom'];
        end
    end
end