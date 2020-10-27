classdef PhysConst
    % PhysConst : Basic physical constants
    % 
    % PhysConst properties:
    %   h   -   (Constant) Planck constant [Js]
    %   e0  -   (Constant) electric permittivity of vacuum [F/m]
    %   m0  -   (Constant) magnetic permeability of vacuum [H/m]
    %   c0  -   (Constant) speed of light in vacuum [m/s]
    %   Z0  -   (Constant) impedance of vacuum [Ohms]
    %   e   -   (Constant) elementary charge [C]
    %   me  -   (Constant) electron mass [kg]
    %   mp  -   (Constant) proton mass [kg]
    %   mn  -   (Constant) neutron mass [kg]
    %   NA  -   (Constant) Avogadro constant [mol^-1]
    %   R   -   (Constant) ideal gas constant [J/(K mol)]
    %   kB  -   (Constant) Boltzmann constant [J/K]

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0  
    %   Date: 2015/01/01
    
    properties (Constant)
        h = 6.62606957e-34; % Planck constant [Js]
        e0 = 8.8541878176e-12; % electric permittivity of vacuum [F/m]
        m0 = 4*pi*1e-7; % magnetic permeability of vacuum [H/m]
        c0 = sqrt(1/(PhysConst.m0*PhysConst.e0)); % speed of light in vacuum [m/s]
        Z0 = sqrt(PhysConst.m0/PhysConst.e0); % impedance of vacuum [Ohms]
        e = 1.602176565e-19; % elementary charge [C]
        me = 9.109382913e-31; % electron mass [kg]
        mp = 1.672621777e-27; % proton mass [kg]
        mn = 1.674927351e-27; % neutron mass [kg]
        NA = 6.02214129e+23 % Avogadro constant [mol^-1]
        R = 8.3144621; % ideal gas constant [J/(K mol)]
        kB = PhysConst.R/PhysConst.NA; % Boltzmann constant [J/K]
    end
end