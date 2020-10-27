classdef PhaseMaskSuperposition < PhaseMask
	% PhaseMaskSuperposition < PhaseMask : Constant phase mask
    %   A phase mask with random phase values.
    %
    % PhaseMaskSuperposition properties:
	%   slm         -   spatial light modulator < PhaseMask 
	%   phase       -   phase shift values < PhaseMask 
	%   alpha_vec   -   polar deflection angles
	%   phi_vec     -   azimuthal deflection angles
	%   f_vec       -   focal lengths
	%   phase0_vec  -   initial phase values [default = 0 rad]
    %
    % PhaseMaskSuperposition methods:
    %   PhaseMaskSuperposition  -   constructor
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
    % See also PhaseMask, SLM, Holography.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    properties
        alpha_vec  % polar deflection angles
        phi_vec  % azimuthal deflection angles
        f_vec  % focal lengths
        phase0_vec  % initial phase values [default = 0 rad]
    end
    methods
        function pm = PhaseMaskSuperposition(slm,alpha_vec,phi_vec,f_vec,varargin)
            % PHASEMASKGRATING(SLM,ALPHAvec,PHIvec,Fvec) constructs a phase
            %   mask representing the superpositon of gratings (ALPHAvec and PHIvec) 
            %   and Fresnel lenses (Fvec).
            % 
            % PHASEMASKGRATING(SLM,ALPHAvec,PHIvec,Fvec,'Phase0',PHASE0vec) adds a
            %   constant phase offset equal to PHASE0vec to each of the
            %   superposition elements.
            %
            % See also PhaseMaskSuperposition.

            pm = pm@PhaseMask(slm);
                        
            % grating phase offset [default = 0 rad]
            phase0_vec = zeros(size(alpha_vec));
            for n = 1:1:length(varargin)
                if strcmpi(varargin{n},'phase0')
                    phase0_vec = varargin{n+1};
                end
            end

            pm.alpha_vec = alpha_vec;
            pm.phi_vec = phi_vec;
            pm.f_vec = f_vec;
            pm.phase0_vec = phase0_vec;

            [X,Y] = pm.slm.pmeshgrid();  % SLM pixel positions [m]
            gl_vec = gratinglengths(pm);
            r = hypot(X,Y); % [m]
            s = zeros(size(X));
            for i = 1:1:length(alpha_vec)
                if alpha_vec(i)~=0
                    s = s + exp( ...
                        2*pi*1i*mod(cos(phi_vec(i)).*X+sin(phi_vec(i)).*Y,gl_vec(i))/gl_vec(i) ...
                        + 1i*mod(2*pi*r.^2/(2*slm.lambda*f_vec(i)),2*pi) ...
                        + 1i*phase0_vec(i) ...
                        );
                else
                    s = s + exp( ...
                        1i*mod(2*pi*r.^2/(2*slm.lambda*f_vec(i)),2*pi) ...
                        + 1i*phase0_vec(i) ...
                        );
                end
            end
            pm.phase = angle(s);
        end
        function disp(pm)
            % DISP Displays phase mask name
            %
            % DISP(PM) displays the phase mask name.
            %
            % See also PhaseMaskSuperposition.

            disp(['<a href="matlab:help PhaseMaskSuperposition">PhaseMaskSuperposition</a> ' num2str(numel(pm.alpha_vec)) ' traps']);
        end
        function txt = title(pm)
            % TITLE Displays phase mask name
            %
            % TXT = TITLE(PM) returns the title for the phase mask plot.
            %
            % See also PhaseMaskSuperposition.

            txt = ['PhaseMaskSuperposition - ' num2str(numel(pm.alpha_vec)) ' traps'];
        end
        function gl_vec = gratinglengths(pm)
            % GRATINGLENGTHS Grating lenghts
            %
            % GLvec = GRATINGLENGTHS(PM) returns the grating lengths.
            %
            % See also PhaseMaskSuperposition.

            gl_vec = pm.slm.lambda./tan(pm.alpha_vec);
        end
    end
end