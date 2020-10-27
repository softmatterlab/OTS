classdef PositionDVM
    % PositionDVM (Abstract) : Particle positions obtained by digital video microscopy
    %   List of the positions of all particles detected in a certain frame
    %   by digital video microscopy. Each particle position should be
    %   indexed by a progressive integer numebr starting from 1.
    %   The units of the positions should be given in pixels.
    %   
    % PositionDVM abstract methods:
    %   plot                -   plot particle positions
    %   numel               -   number of particles
    %   extract             -   position of subset of particles
    %   remove              -   remove particles
    %   distance            -   distance between two sets of particles
    %   gettrajectories     -   transform a position to a series of trajectories
    %   
    % See also TrajectoryDVM, PositionDVM2D.

    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    methods (Abstract)
        plot(pos)  % plot the particle positions
        numel(pos)  % number of particles
        extract(pos,index)  % position of subset of particles
        remove(pos,index)  % remove particles
        distance(pos1,pos2)  % distance between two sets of particles
        gettrajectories(pos,t0)  % transform a position to a series of trajectories
    end
end