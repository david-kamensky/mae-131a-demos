% This script directly solves for the forces in 
% elements of a statically-determinate end-to-end 
% axial deformation assemblage, then post-processes
% a displacement field from those.

%%%%%%% Input data %%%%%%%

% Array of Young's moduli of elements:
E = [1,1,1,1];
% Array of cross-sectional areas of elements:
A = [1,1,1,1];
% Array of lengths of nodal positions:
x = [0,0.1,0.2,0.3,0.4];
% Array of nodal forces (no force at constrained left end):
P = [0,0,1,1];
% Displacement BC at left end:
u_0 = 0;

%%%%%%% Derived quantities %%%%%%%

% Number of elements:
N_el = numel(E);
% Number of nodes:
N_node = N_el + 1;
% Array of element lengths:
L = x(2:N_node) - x(1:N_el);
% Array of stiffnesses of elements:
k = E.*A./L;

%%%%%%% Solve for forces one-by-one %%%%%%%

F = zeros(1,N_el);
F(N_el) = P(N_node-1);
for i=N_el-1:-1:1
    F(i) = F(i+1) + P(i);
end

%%%%%%% Post-process results %%%%%%%

% Stress in each element:
sigma = F./A;
% Strain in each element:
epsilon = sigma./E;
% Elongation of each element:
elong = epsilon.*L;
%elong = F./k; % (equivalent)
% Displacement of each node:
u = zeros(1,N_node);
u(1) = u_0;
for i=2:N_node
    u(i) = u(i-1) + elong(i-1);
end

plot(x,u);

