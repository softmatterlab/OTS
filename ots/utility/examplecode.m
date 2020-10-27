function examplecode(codeline,pausetime)
% EXAMPLECODE   Auxiliary function for examples - Execute code
%
% EXAMPLECODE(CODELINE,PAUSETIME) executes and displays CODELINE 
%   and then waits for PAUSETIME seconds.
%
% EXAMPLECODE(CODELINE) executes and displays CODELINE 
%   and then waits for 1 second.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

if nargin<2
    pausetime = 1;
end

fprintf(['>> ' codeline '\n'])
evalin('caller',codeline)
pause(pausetime)