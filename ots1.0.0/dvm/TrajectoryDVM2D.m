classdef TrajectoryDVM2D < TrajectoryDVM
    % TrajectoryDVM2D < <a href="matlab:help TrajectoryDVM">TrajectoryDVM</a> : 2D particle trajectory obtained by digital video microscopy
    %   Each particle has a position X and Y measured either in pixels.
    %   Furthermore, each particle is also characterized by its area Area in pixel^2.
    %   This implementation works together with Position2D, DVM2DThreshold and DVM2DGauss.
    %
    % TrajectoryDVM2D properties:
    %   T       -   time [frame number]
    %   X       -   x-position [pixel]
    %   Y       -   y-position [pixel]
    %   Area    -   area [pixel^2]
    %
    % TrajectoryDVM2D methods:
    %   TrajectoryDVM2D -   constructor
    %   plot            -   plot trajectory
    %   getposition     -   get the positon corresponding to index
    %   gettime         -   get the time corresponding to index
    %   append          -   append position to the end of the trajectory
    %   
    % See also TrajectoryDVM, Position2D, DVM2DThreshold, DVM2DGauss.

    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties
        T       % time [frame number]
        X       % x-position [pixel]
        Y       % y-position [pixel]
        Area    % area [pixel^2]
    end
    methods
        function obj = TrajectoryDVM2D(T,X,Y,Area)
            % TRAJECTORYDVM2D(T,X,Y,AREA) constructs a trajectory with times
            %   T, coordinates X and Y, and area AREA.
            %
            % See also TrajectoryDVM2D.
            
            if ~isreal(T) || ~isreal(X) || ~isreal(Y) || numel(T)~=numel(X) || numel(T)~=numel(Y) || numel(X)~=numel(Y)
                error('T, X and Y should be real matrices with the same number of elements.')
            end
            
            obj.T = reshape(T,1,numel(T));
            obj.X = reshape(X,1,numel(X));
            obj.Y = reshape(Y,1,numel(Y));
            obj.Area = reshape(Area,1,numel(Area));
        end
        function h = plot(traj,varargin)
            % PLOT Plot trajectory
            %
            % H = PLOT(TRAJ) plots the trajectory TRAJ. 
            %   It returns a graphic handler H to the plotted set of positions.
            %
            % H = PLOT(TRAJ,'PropertyName',PropertyValue) sets the property PropertyName 
            %   to PropertyValue. All standard plot properties can be used.
            %
            % See also TrajectoryDVM2D.

            h = plot(traj.X,traj.Y,'k');
            axis equal
            for n = 1:2:length(varargin)
                set(h,varargin{n},varargin{n+1});
            end
        end
        function pos = getposition(traj,index)
            % GETPOSITION Get the positon corresponding to index
            %
            % POS = GETPOSITION(TRAJ,INDEX) returns the position POS
            %   corresponding to the INDEX-th element of the trajectory TRAJ.
            %
            % See also TrajectoryDVM2D, Position2D.

            if nargin==1
                index = length(traj.T);
            end
            
            if ~isnumeric(index) || ~isreal(index) || index<1 || index>length(traj.T)
                error('The index must be an integer number greater than 0 and smaller than the trajectory length.')
            end
            
            pos = Position2D(traj.X(index),traj.Y(index),traj.Area(index));
        end
        function t = gettime(traj,index)
            % GETTIME Get the time corresponding to index
            %
            % T = GETTIME(TRAJ,INDEX) returns the time T
            %   corresponding to the INDEX-th element of the trajectory TRAJ.
            %
            % See also TrajectoryDVM2D.

            if nargin==1
                index = length(traj.T);
            end
            
            if ~isnumeric(index) || ~isreal(index) || index<1 || index>length(traj.T)
                error('The index must be an integer number greater than 0 and smaller than the trajectory length.')
            end
            
            t = traj.T(index);
        end
        function traj_a = append(traj,pos,index,t)
            % APPEND Append position to the end of the trajectory
            %
            % TRAJ_A = APPEND(TRAJ,POS,INDEX,T) returns the trajectory TRAJ_A 
            %   corresponding to the trajectory TRAJ with POS(INDEX) appended 
            %   with time stap T.
            %
            % TRAJ_A = APPEND(TRAJ,POS,INDEX) T is set to the last time of TRAJ plus 1.
            %
            % See also TrajectoryDVM2D, Position2D.            
            
            if nargin==3
                t = traj.T(end)+1;
            end
            if ~isnumeric(t) || ~isreal(t)
                error('The time of the trajectory must be a real number.')
            end
            
            if ~isa(pos,'Position2D')
                error('The position to be appended must be a DVMPosition2D.')
            end
            
            traj_a = TrajectoryDVM2D([traj.T t],[traj.X pos.X(index)],[traj.Y pos.Y(index)],[traj.Area pos.Area(index)]);
        end
    end
end