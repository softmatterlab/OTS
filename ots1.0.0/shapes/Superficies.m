classdef Superficies < Shape
    % Superficies (Abstract) < Shape : 2D surface embedded in 3D space
    %   A surface is a shape for which is possible to define a normal line
    %   and a tangential plane at each point.
    % 
    % Superficies abstract methods:
    %   intersectionpoint   -   (Abstract) intersection point set with line/vector set
    %   perpline            -   (Abstract) perpendicular line at point
    %   tangentplane        -   (Abstract) tangent plane set passing by point set
    %
    % See also Shape, Point, Vector, SLine, Plane, Spherical, Ellipsoidal, Cylindrical.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    methods (Abstract)
        intersectionpoint(sup,d,n)  % intersection point set with line/vector set
        perpline(sup,p)  % perpendicular line at point
        tangentplane(sup,p)  % tangent plane set passing by point set
    end
end