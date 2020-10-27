classdef SLM
    % SLM : Spatial light modulator
    %   A spatial light modulator is a screen with M times N pixels of
    %   pixel size psize and capable of modulating the phase of a beam of
    %   wavelength lambda with a numebr of steps equal to levels.
    %
    % SLM properties:
	%   M       -   SLM pixel rows
	%   N       -   SLM pixel columns
	%   psize   -   pixel size [m]
	%   lambda  -   light wavelength
	%   levels  -   number of grey levels for 2pi modulation
    %
    % SLM methods:
    %   SLM         - constructor
    %   pmeshgrid   - pixels positions
    %
    % See also Holography.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        M       % SLM pixel rows
        N       % SLM pixel columns
        psize   % pixel size [m]
        lambda  % light wavelength
        levels  % number of grey levels for 2pi modulation
    end
    methods
        function slm = SLM(M,N,psize,lambda,levels)
            % SLM(M,N,PSIZE,LAMBDA) constructs an SLM with MxN pixels of
            %   size PSIZE and sets the wavelength at LAMBDA.
            %
            % SLM(M,N,PSIZE,LAMBDA,LEVELS) sets the number of levels to
            %   LEVELS (defauls = 255).
            %
            % See also SLM.
            
            if nargin<5
                levels = 255;
            end
            
            slm.M = M;
            slm.N = N;
            slm.psize = psize;
            slm.lambda = lambda;
            slm.levels = levels;
        end
        function [X,Y] = pmeshgrid(slm)
            % PMESHGRID SLM pixels positions
            %
            % [X,Y] = PMESHGRID(SLM) returns the SLM pixels positions.
            %
            % See also SLM.

            [X,Y] = meshgrid( ...
                ([.5:1:slm.M-.5]-slm.M/2)*slm.psize, ...
                ([.5:1:slm.N-.5]-slm.N/2)*slm.psize ...
                );  % SLM pixel positions [m]
        end
    end
end