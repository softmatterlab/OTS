classdef PhaseMaskGrating < PhaseMask
	% PhaseMaskGrating < PhaseMask : Constant phase mask
    %   A phase mask representing a diffraciton grating with deflection
    %   angles alpha and phi.
    %
    % PhaseMaskGrating properties:
	%   slm     -   spatial light modulator < PhaseMask 
	%   phase	-   phase shift values < PhaseMask 
	%   alpha   -   polar deflection angle
	%   phi     -   azimuthal deflection angle
    %   phase0  -   phase value [default = 0 rad]
    %
    % PhaseMaskGrating methods:
    %   PhaseMaskGrating   -   constructor
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
    %   gratinglength       -   grating length
    %
    % See also PhaseMask, SLM, Holography.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        alpha  % polar deflection angle
        phi  % azimuthal deflection angle
        phase0  % phase value [default = 0 rad]
    end
    methods
        function pm = PhaseMaskGrating(slm,alpha,phi,varargin)
            % PHASEMASKGRATING(SLM,ALPHA,PHI) constructs a diffraction
            %   grating that produces a polar deflection of ALPHA and an
            %   azimuthal deflection of PHI.
            % 
            % PHASEMASKGRATING(SLM,ALPHA,PHI,'Phase0',PHASE0) adds a
            %   constant phase offset equal to PHASE0.
            %
            % See also PhaseMaskGrating.

            pm = pm@PhaseMask(slm);
                        
            % grating phase offset [default = 0 rad]
            phase0 = 0;
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'phase0')
                    phase0 = varargin{n+1};
                end
            end

            pm.alpha = alpha;
            pm.phi = phi;
            pm.phase0 = phase0;
            
            % grating length [m]
            if alpha~=0
                [X,Y] = pm.slm.pmeshgrid();  % SLM pixel positions [m]
                pm.phase = mod(cos(phi).*X+sin(phi).*Y,pm.gratinglength())/pm.gratinglength()*2*pi;
                pm.phase = mod(pm.phase + phase0, 2*pi);
            end
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskGrating.

            disp(['<a href="matlab:help PhaseMaskGrating">PhaseMaskGrating</a> alpha=' num2str(pm.alpha/(2*pi)*360) 'deg phi=' num2str(pm.phi/(2*pi)*360) 'deg phase0=' num2str(pm.phase0/(2*pi)*360) 'deg']);
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskGrating.

            txt = ['PhaseMaskGrating - $\alpha=' num2str(pm.alpha/(2*pi)*360) '^{\circ}$ $\phi=' num2str(pm.phi/(2*pi)*360) '^{\circ}$ $\varphi_0=' num2str(pm.phase0/(2*pi)*360) '^{\circ}$'];
        end
        function gl = gratinglength(pm)
            % GRATINGLENGTH Grating lenght
            %
            % GL = GRATINGLENGTH(PM) returns the grating length.
            %
            % See also PhaseMaskGrating.
            
            gl = pm.slm.lambda/tan(pm.alpha);
        end
    end
end