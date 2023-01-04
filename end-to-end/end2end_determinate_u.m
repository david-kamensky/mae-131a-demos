% This script solves directly for nodal displacements in
% a statically determinate end-to-end axial deformation
% assemblage.  Compare the resulting displacements with
% those post-processed from solving directly for
% forces in the script end2end_determinate_F.m (for
% the same input data).

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

%%%%%%% Assemble the equation system %%%%%%%

K = zeros(N_node,N_node);
K(1,1) = 1;
for i=2:N_node-1
    K(i,i-1:i+1) = [-k(i-1) , k(i-1)+k(i) , -k(i)];
end
K(N_node,N_node-1:N_node) = [-k(N_el) , k(N_el)];

%%%%%%% Solve system for nodal displacements %%%%%%%

u = K\[u_0;transpose(P)];

%%%%%%% Post-process results %%%%%%%

% Strain in each element:
epsilon = transpose(u(2:N_node) - u(1:N_el))./L;
% Stress in each element:
sigma = E.*epsilon;
% Force in each element:
F = sigma.*A;

plot(x,u);
