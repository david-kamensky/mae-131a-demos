% This script solves for the member forces in a statically determinate
% planar truss.  If the given truss is not determinate, the resulting 
% system of equations will be singular.  This should be review from
% statics, but we will solve a number of determinate trusses for review
% and for setting up some determinate solid mechanics problems.  (For 
% indeterminate planar trusses, see the script "truss_indeterminate.m", 
% which solves for nodal displacement vectors using solid mechanics 
% instead of just forces using statics, but that is outside the scope of 
% the course.)

%%%%%%% Specify problem data %%%%%%%

% Each row gives the x,y position of a node in the truss.
nodes = [0,0;2,0;4,0;0,-2;2,-1];

% Each row gives two node indices (corresponding to rows in `nodes`),
% indicating the two endpoints of a truss element.
elements = [1,2;2,3;1,5;2,5;1,4;4,5;5,3];

% Each row gives nonzero or 0 for each component of displacement at a node,
% where a nonzero indicates that the corresponding component is fixed.  
% Note that there must be exactly three ones for the truss to be 
% determinate, and the nonzero values entered for these should be 1, 2,
% and 3, which will be used as indices for the three unknown reaction
% forces associated with these constraints.
BCs = [1,2;0,0;0,0;3,0;0,0];

% Each row gives the x,y components of a force applied at the 
% corresponding joint.
node_forces = [0,0;0,0;0,-1;0,0;0,0];

%%%%%%% End of problem setup %%%%%%%

% Number of nodes:
N_node = size(nodes,1);

% Number of elements:
N_el = size(elements,1);

% Total degrees of freedom (element forces and reaction forces):
N_dof = N_el + 3;

% Set up matrix and vector for the system of equations, where
% each row corresponds to equilibrium of a given node in a Cartesian
% direction, with ordering
% 
% [x equilibrium for node 1;
%  x equilibrium for node 2;
%  ...
%  x equilibrium for node N_node;
%  y equilibrium for node 1;
%  ...
%  y equilibrium for node N_node] ,
%
% implemented as (direction==1 => x, direction==2 => y):
row_index = @(node,direction)(node + (direction-1)*N_node);
%
% The vector of unknowns is ordered as:
%
% [force in element 1;
%  force in element 2;
%  ...
%  force in element N_el;
%  reaction force 1;
%  reaction force 2;
%  reaction force 3] .

A = zeros(2*N_node,N_dof);
B = zeros(N_dof,1);

% Loop over elements and add contributions to corresponding nodal
% equilibrium equations.
for el=1:N_el

    % "Un-pack" data for the current element from the arrays used to 
    % define our problem.

    % Node numbers of endpoints of current element:
    n1 = elements(el,1);
    n2 = elements(el,2);
    % x and y coordinates of first endpoint:
    x1 = nodes(n1,1);
    y1 = nodes(n1,2);
    % x and y coordinates of second endpoint:
    x2 = nodes(n2,1);
    y2 = nodes(n2,2);

    % Element geometry:
    Dx = x2-x1;
    Dy = y2-y1;
    L = sqrt(Dx^2 + Dy^2);
    c = Dx/L; % cos(counterclockwise angle of vector 1->2 from x-axis)
    s = Dy/L; % sin(...)

    % Row indices in linear system:
    row1x = row_index(n1,1);
    row1y = row_index(n1,2);
    row2x = row_index(n2,1);
    row2y = row_index(n2,2);

    % Matrix contributions from unknown truss member forces; remember
    % each row of A is dotted with the vector of unknown forces in
    % the corresponding algebraic equation and we follow the convention
    % that positive member forces are in tension and negative ones are
    % in compression.
    A(row1x,el) = c;
    A(row2x,el) = -c;
    A(row1y,el) = s;
    A(row2y,el) = -s;

end % end for el

% Loop over nodes and add contributions from reaction and external forces
% in x and y directions:
for node=1:N_node
    % Loop over directions:
    for direction=1:2
        
        % Add contribution from reaction force, if present for current
        % node and direction:
        if(BCs(node,direction) > 0)
            A(row_index(node,direction),BCs(node,direction)+N_el) = 1;
        end

        % Add external force; negative sign is because it is on the 
        % opposite side of the matrix equation 
        % 
        %  A*F = B
        %
        % from the contributions of the unknown member and reaction 
        % forces in the vector F.
        B(row_index(node,direction)) = -node_forces(node,direction);

    end % end for direction
end % end for node

% Solve for and print the member and reaction forces:
F = A\B;

% Print out solution of member and reaction forces:
fprintf("Truss member forces:\n");
disp(F(1:N_el));
fprintf("Reaction forces:\n");
disp(F(N_el+1:N_el+3));