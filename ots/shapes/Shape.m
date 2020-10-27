classdef Shape
    % Shape (Abstract) : 3D shape
    %   A shape is a geometrical object, or a set of geometric objects,
    %   that can be plotted, translated and rotated in 3D. 
    % 
    % Shape abstract methods:
    %   plot            -   (Abstract) plots shape set in 3D
    %   disp            -   (Abstract) prints shape set
    %   translate       -   (Abstract) translation
    %   xrotation       -   (Abstract) rotation around x-axis
    %   yrotation       -   (Abstract) rotation around y-axis
    %   zrotation       -   (Abstract) rotation around z-axis
    %   numel           -   (Abstract) number of shapes
    %   size            -   (Abstract) size of shape set
    %
    % See also Point, Vector, SLine, Superficies, Plane, Spherical, Elliprical, Cylindrical.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01

    methods (Abstract)
        plot(obj)  % plots shape set in 3D
        disp(obj)  % prints shape set
        translate(obj,V)  % translates shape set by V
        xrotation(obj,phi)  % rotates shape set around x-axis by phi
        yrotation(obj,phi)  % rotates shape set around y-axis by phi
        zrotation(obj,phi)  % rotates shape set around z-axis by phi
        numel(obj)  % number of shapes
        size(obj,dim)  % size of shape set matrix
    end
end