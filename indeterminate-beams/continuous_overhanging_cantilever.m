% Script to solve for the constants of integration
% and unknown reaction forces in a continuous
% overhanging cantilever beam subject to a tip load.
%
% Setup: p(x) = 0, cantilever support at x=0,
%        tip load -W at x=L, and roller supports at
%        x=L/3 and x=2*L/3.  

% Choose numerical values for constants:
W = 2.0;
L = 12.0;
EI = 1.0e3;

% Matrix--vector form of equation system 
% derived in lecture:
A = [1.0           , 1.0   , 1.0               , 0.0              ;
     2.0*L/3.0     , L/3.0 , L                 , 1.0              ;
     0.0           , 0.0   , (L/3.0)^3/6.0     , (L/3.0)^2/2.0    ;
     (L/3.0)^3/6.0 , 0.0   , (2.0*L/3.0)^3/6.0 , (2.0*L/3.0)^2/2.0];
b = [W  ;
     0.0;
     0.0;
     0.0];

% Use linear solver (equivalent to x = inv(A)*b
% but better practice for large systems):
x = A\b;

% Unpack vector of unknowns to symbols
% used in notes:
By = x(1)
Cy = x(2)
C1 = x(3)
C2 = x(4)

% Postprocess results:

% Heaviside function:
H = @(x)(x>0); 
% Deflection derived in lecture:
v = @(x)((1.0/EI) ...
          *(By/6.0*((x-L/3.0).^3).*H(x-L/3.0) ...
            + Cy/6.0*((x-2.0*L/3.0).^3).*H(x-2.0*L/3.0) ...
            + C1*x.^3/6.0 + C2*x.^2/2.0));
% Internal moment resultants (M = EIv''):
M = @(x)(By*(x-L/3.0).*H(x-L/3.0) ...
         + Cy*(x-2.0*L/3.0).*H(x-2.0*L/3.0) ...
         + C1*x + C2);
% Internal shear resultants (V = (EIv'')'):
V = @(x)(By*H(x-L/3.0) + Cy*H(x-2.0*L/3.0) + C1);

% Plots:
N_pts = 10000;
x = linspace(0,L,N_pts);
figure;
plot(x,v(x));
figure;
plot(x,M(x));
figure;
plot(x,V(x));
