function [D,lambda,lambda_app,D0] = brenner(h,R,eta,T,fig)
% BRENNER   Diffusion coefficient in front of a surface
%
% [D,L,L1,D0] = BRENNER(H,R,ETA,T,FIG) calculates the diffusion coefficient D
%   in front of a surface at distances H assuming a spherical particle of
%   radius R in afluid with viscosity ETA and temperature T.
%   FIG determines whether to show or not a figure.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

if (nargin<5)
    fig = false;
end
if (nargin==0)
    h = [1e-9:1e-8:.2e-5];
    T = 293;
    eta = 0.001;
    R = 1e-6;
    fig = true;
end

kB = 1.38e-23;
D0 = kB*T/(6*pi*eta*R);

alpha = acosh(1+h/R);

S = 0;
for n = 1:1:10
    S = S + ...
        n*(n+1)/((2*n-1)*(2*n+3)) * ...
        ( ( 2*sinh((2*n+1)*alpha) + (2*n+1) * sinh(2*alpha) ) ./...
        ( 4*sinh((n+0.5)*alpha).^2 - (2*n+1)^2 * sinh(alpha).^2 ) - 1 );
end
lambda = 4/3*sinh(alpha).*S;

lambda_app = R*h.^-1 + 0.2*log(R*h.^-1) + 0.9712;

D = D0*lambda.^-1;

if (fig==true)
    f = figure
    hold on
    plot(h/R,lambda.^-1)
    plot(h/R,lambda_app.^-1,'--')
    plot(h/R,h/R,':')
    axis([0 2 0 1])
    xlabel('z/R')
    ylabel('D/D_0')
    box on
    saveas(f,'brenner.eps','eps')
end