classdef TrajectoryDVM
    % TrajectoryDVM (Abstract) : Particle trajectory obtained by digital video microscopy
    %   Trajectory of a particle detected by digital video microscopy.
    %   The units of the positions should be given in pixels and the unit
    %   of time should be given in frame number.
    %   
    % TrajectoryDVM abstract methods:
    %   plot            -   plot trajectory
    %   getposition     -   get the positon corresponding to index [pixels]
    %   gettime         -   get the time corresponding to index [frame number]
    %   append          -   append position to the end of the trajectory
    %   
    % See also PositionDVM, TrajectoryDVM2D.

    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    methods (Abstract)
        plot(traj,varargin)  % plot trajectory
        getposition(traj,index)  % get the positon corresponding to index [pixels]
        gettime(traj,index)  % get the time corresponding to index [frame number]
        append(traj,pos)  % append position to the end of the trajectory
    end
end