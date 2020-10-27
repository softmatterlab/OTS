% Series of examples to demonstrate the use of PhaseMask.
%
% See also PhaseMask, SLM.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

close all; 
clear all; 
clc;

slm = SLM(792,600,20e-6,1064e-9);

pm = PhaseMask(slm)

pm.title()

fig = pm.plot()