% This script solves for nodal displacements in a
% statically indeterminate end-to-end axial
% deformation assemblage, then post-processes
% forces from those.

%%%%%%% Input data %%%%%%%

% Array of Young's moduli of elements:
E = [1,1,1,1];
% Array of cross-sectional areas of elements:
A = [1,1,1,1];
% Array of lengths of nodal positions:
x = [0,0.1,0.2,0.3,0.4];
% Array of nodal forces (at interior nodes only):
P = [0,10,0];
% Displacement BCs at both ends:
u_0 = 0;
u_L = 1;

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
K(N_node,N_node) = 1;
for i=2:N_node-1
    K(i,i-1:i+1) = [-k(i-1) , k(i-1)+k(i) , -k(i)];
end

%%%%%%% Solve system for nodal displacements %%%%%%%

u = K\[u_0;transpose(P);u_L];

%%%%%%% Post-process results %%%%%%%

% Strain in each element:
epsilon = transpose(u(2:N_node) - u(1:N_el))./L;
% Stress in each element:
sigma = E.*epsilon;
% Force in each element:
F = sigma.*A;

plot(x,u);
