classdef Particle < Shape
    % Particle (Abstract) < Shape : Optically trappable particle
    %   A particle is an object that can be optically trapped, i.e., when
    %   an electromagnetic field impinges on it, scattering occurs and
    %   forces and torques are produced.
    % 
    % Particle methods (Abstract):
    %   plot            -   (Abstract) plots shape set in 3D < Shape
    %   disp            -   (Abstract) prints shape set < Shape
    %   translate       -   (Abstract) translation < Shape
    %   xrotation       -   (Abstract) rotation around x-axis < Shape
    %   yrotation       -   (Abstract) rotation around y-axis < Shape
    %   zrotation       -   (Abstract) rotation around z-axis < Shape
    %   numel           -   (Abstract) number of shapes < Shape
    %   size            -   (Abstract) size of shape set < Shape
    %   barycenter      -   (Abstract) particle center of mass
    %   scattering      -   (Abstract) scattered rays
    %   force           -   (Abstract) force due to a set of rays
    %   torque          -   (Abstract) torque due to a set of rays
    %
    % See also Shape, ParticleSpherical, ParticleEllipsoidal, ParticleCylindrical, Ray.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    methods (Abstract)
        barycenter(par)  % particle center of mass
        scattering(par,r)  % scattered rays
        force(par,r)  % force due to a set of rays
        torque(par,r)  % torque due to a set of rays
    end
end