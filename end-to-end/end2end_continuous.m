% This script approximates the indeterminate case of the
% tapered bar covered in lecture as an end-to-end
% assemblage of uniform axial deformation elements, 
% and compares the result with the analytical solution
% derived in lecture.

%%%%%%% Input data %%%%%%%

% Young's modulus:
E = 1;
% Geometry:
L_tot = 10; % Total length of bar
R_0 = 1; % Radius at x=0
R_L = 0.5; % Radius at x=L_tot
% Total number of nodes:
N_node = 10;
% Displacement BCs at both ends:
u_0 = 0;
u_L = 1;
% Forces at interior nodes:
P = zeros(1,N_node-2);

%%%%%%% Derived quantities %%%%%%%

% "Mesh" of nodal positions:
x = linspace(0,L_tot,N_node);
% Number of elements:
N_el = N_node - 1;
% Array of element lengths:
L = x(2:N_node) - x(1:N_el);
% Array of cross-sectional areas at nodes:
R_func = @(x) R_0 - ((R_0-R_L)/L_tot)*x;
A_nodal = pi*R_func(x).^2;
% Cross-sectional area for each element:
A = 0.5*(A_nodal(1:N_el) + A_nodal(2:N_node));
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

% Compare with ODE solution from lecture:
u_exact = @(x) (u_L*R_L/(R_0-R_L))*(R_0/(R_0-((R_0-R_L)/L_tot)*x)-1);
plot(x,arrayfun(u_exact,x));
hold on;
scatter(x,u);

