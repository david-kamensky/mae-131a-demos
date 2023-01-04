% Script to solve for the constants of integration
% and unknown reaction force from the propped
% overhanging cantilever beam example from lecture.
%
% Setup: p(x) = P*x/L, cantilever support at x=0,
%        free end at x=L, and roller support at
%        x=L/2.  

% Choose numerical values for constants:
P = 2.0;
L = 10.0;
EI = 1.0e3;

% Matrix--vector form of equation system 
% derived in lecture:
A = [0     , (L/2.0)^3/6.0 , (L/2.0)^2/2.0 ;
     L/2.0 , L             , 1.0           ;
     1.0   , 1.0           , 0.0            ];
b = [-P*(L/2.0)^5/(120*L) ;
     -P*L^2/6.0           ;
     -P*L/2.0              ];

% Use linear solver (equivalent to x = inv(A)*b
% but better practice for large systems):
x = A\b;

% Unpack vector of unknowns to symbols
% used in notes:
By = x(1)
C1 = x(2)
C2 = x(3)

% Postprocess results:

% Heaviside function:
H = @(x)(x>0); 
% Deflection derived in lecture:
v = @(x)((1.0/EI) ...
          *(P*x.^5/(120.0*L) ...
            + (By/6.0)*((x-L/2.0).^3).*H(x-L/2.0) ...
            + C1*x.^3/6.0 + C2*x.^2/2.0));
% Internal moment resultants (M = EIv''):
M = @(x)(P*(x.^3)/(6.0*L) + By*(x-L/2.0).*H(x-L/2.0) ...
         + C1*x + C2);
% Internal shear resultants (V = (EIv'')'):
V = @(x)(P*(x.^2)/(2.0*L) + By*H(x-L/2.0) + C1);

% Plots:
N_pts = 10000;
x = linspace(0,L,N_pts);
figure;
plot(x,v(x));
figure;
plot(x,M(x));
figure;
plot(x,V(x));
