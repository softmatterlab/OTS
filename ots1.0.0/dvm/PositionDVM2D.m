classdef PositionDVM2D < PositionDVM
    % PositionDVM2D < <a href="matlab:help PositionDVM">PositionDVM</a> : 2D particle positions obtained by digital video microscopy
    %   Each particle has a position X and Y measured either in pixels.
    %   Furthermore, each particle is also characterized by its area Area in pixel^2.
    %   This implementation works together with Trajectory2D, DVM2DThreshold and DVM2DGauss.
    %
    % PositionDVM2D properties:
    %   X       -   x-position [pixel]
    %   Y       -   y-position [pixel]
    %   Area    -   area [pixel^2]
    %
    % PositionDVM2D methods:
    %   PositionDVM2D   -   constructor
    %   plot            -   plot particle positions
    %   numel           -   number of particles
    %   extract         -   position of subset of particles
    %   remove          -   remove particles
    %   distance        -   distance between two sets of particles
    %   gettrajectories -   transform a position to a series of trajectories
    %   
    % See also PositionDVM, Trajectory2D, DVM2DThreshold, DVM2DGauss.

    %   Author: Giuseppe Pesce, Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        X       % x-position [pixel]
        Y       % y-position [pixel]
        Area    % area [pixel^2]
    end
    methods
        function obj = PositionDVM2D(X,Y,Area)
            % POSITIONDVM2D(X,Y,AREA) constructs a set of positions of particles 
            %   with coordinates X and Y and area AREA. X, Y and AREA should be real matrices with the
            %   same number of elements, otherwise this function return an error.
            %
            % See also PositionDVM2D.
            
            if ~isreal(X) || ~isreal(Y) || numel(X)~=numel(Y)
                error('X and Y should be real matrices with the same number of elements.')
            end
            
            obj.X = reshape(X,numel(X),1);
            obj.Y = reshape(Y,numel(Y),1);
            obj.Area = reshape(Area,numel(Area),1);
            
        end
        function h = plot(pos,varargin)
            % PLOT Plot particle positions
            %
            % H = PLOT(POS) plots the particle positions POS. 
            %   It returns a graphic handler H to the plotted set of positions.
            %
            % H = PLOT(POS,'PropertyName',PropertyValue) sets the property PropertyName 
            %   to PropertyValue. All standard plot properties can be used.
            %
            % See also PositionDVM2D.
            
            h = plot(pos.X,pos.Y,'.k');
            axis equal
            for n = 1:2:length(varargin)
                set(h,varargin{n},varargin{n+1});
            end
        end
        function n = numel(pos)
            % NUMEL Number of particles
            %
            % N = NUMEL(POS) returns the number N of particles of POS.
            %
            % See also PositionDVM2D.

            n = length(pos.X);
        end
        function pos_e = extract(pos,index)
            % EXTRACT Position of subset of particles
            %
            % POS_E = EXTRACT(POS,INDEX) returns the positions POS_E 
            %   of the particles with indices INDEX extracted from POS.
            %
            % See also PositionDVM2D.

            pos_e = PositionDVM2D(pos.X(index),pos.Y(index),pos.Area(index));
        end
        function pos_r = remove(pos,index)
            % REMOVE Remove particles
            %
            % POS_R = REMOVE(POS,INDEX) returns the positions POS_R 
            %   of all the particles of POS exept the ones with indices INDEX.
            %
            % See also PositionDVM2D.

            pos_r = PositionDVM2D(pos.X([1:index-1 index+1:end]),pos.Y([1:index-1 index+1:end]),pos.Area([1:index-1 index+1:end]));
        end
        function d = distance(pos1,pos2)
            % DISTANCE Distance between two sets of particles
            %
            % D = DISTANCE(POS1,POS2) calculates the matrix D of the distances 
            %   between the particles in POS1 and the particles in POS2.
            %   D(i,j) is the distance between the i-th particle of POS1 
            %   and the j-th particle of POS2.
            %
            % See also PositionDVM2D.

            d = sqrt( ...
                ( pos1.X*ones(1,pos2.numel()) - ones(pos1.numel(),1)*pos2.X' ).^2 + ...
                ( pos1.Y*ones(1,pos2.numel()) - ones(pos1.numel(),1)*pos2.Y' ).^2 ...
                );
        end
        function Trajectories = gettrajectories(pos,t0)
            % GETTRAJECTORIES Transform a position to a series of trajectories
            %
            % TRAJECTORIES = GETTRAJECTORIES(POS,T0) transforms the position POS to a series of
            % trajectories (<a href="matlab:help Trajectory2D">Trajectory2D</a>) with initial time T0.
            % This function returns an error if T0 is not a real number.
            %
            % TRAJECTORIES = GETTRAJECTORIES(POS) sets the initial time at 0.
            %
            % See also PositionDVM2D, Trajectory2D.
            
            % initial time [default = 0 frames]
            if nargin==1
                t0 = 0;
            end
            if ~isnumeric(t0) || ~isreal(t0)
                error('The initial time of the trajectory must be a real number.')
            end
            
            for i = 1:1:pos.numel()
                Trajectories(i) = Trajectory2D(t0,pos.X(i),pos.Y(i),pos.Area(i));
            end
            
        end
    end
end