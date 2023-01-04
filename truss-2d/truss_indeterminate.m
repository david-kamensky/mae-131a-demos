% This script is outside the scope of the course material covered in
% lectures, assignments, and exams, but it may be of interest to some
% students.  It solves for the nodal displacements of a
% (possibly-indeterminate) planar truss, as a 2D extension of our analysis
% of 1D end-to-end assemblages.  In this case, extra bookkeeping is needed
% to manage the connectivity of the truss structure and the presence of
% multiple degrees of freedom (both x- and y-components of displacement)
% at each node.

%%%%%%% Specify problem data %%%%%%%

% Each row gives the x,y position of a node in the truss.
nodes = [0,0;0,1;0,2;1,0;1,1;2,0;2,1;3,0];

% Each row gives two node indices (corresponding to rows in `nodes`),
% indicating the two endpoints of a truss element.
elements = [1,2;2,3;1,4;2,4;2,5;3,7;4,6;4,7;5,7;6,8;7,8;4,5;6,7;3,5];

% Each row gives 1 or 0 for each component of displacement at a node,
% where a 1 indicates that the corresponding component is fixed.
BCs = [1,0;1,0;1,1;0,0;0,0;0,0;0,0;0,0];

% Each row gives the stiffness of the corresponding element:
k = 1e2*[1;1;1;1;1;1;1;1;1;1;1;1;1;1];

% Each row gives the x,y components of a force applied at the 
% corresponding joint.
node_forces = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,-3];

%%%%%%% End of problem setup %%%%%%%

% Number of nodes:
N_node = size(nodes,1);

% Number of elements:
N_el = size(elements,1);

% Number of degrees of freedom (DoFs):
N_dof = 2*N_node;

% Mapping from a node index and a Cartesian direction (where 
% direction=1 => x and direction=2 => y) to a degree of freedom (DoF) 
% index:
dof_index = @(node,direction)(2*(node-1)+(direction-1)+1);

% Set up stiffness matrix and load vector:
K = zeros(N_dof,N_dof);
F = zeros(N_dof,1);

% Loop over elements and add stiffness contributions to corresponding 
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
    ke = k(el);
    % Indices of x and y displacements in the vector of all
    % unknown scalar displacement components:
    dof1x = dof_index(n1,1);
    dof1y = dof_index(n1,2);
    dof2x = dof_index(n2,1);
    dof2y = dof_index(n2,2);

    % Element geometry:
    Dx = x2-x1;
    Dy = y2-y1;
    L = sqrt(Dx^2 + Dy^2);

    % Add element's contributions to linear system of equilibrium
    % equations at nodes.

    % Linearized change in length (valid for small rotations only):
    %
    %  DL = (u2x-u1x)*Dx/L + (u2y-u2y)*Dy/L
    %
    % This is the dot product of the displacement difference between
    % nodes and a unit vector in the direction of the member. See Section 
    % 3.10 of the Craig textbook. Figures 3.23--3.24 and 
    % Equation (3.29) cover the special case where node 1 is fixed.
    %
    % Element force:
    %
    %  fe = ke*DL = ke*((u2x-u1x)*Dx/L + (u2y-u1y)*Dy/L)
    %
    % x,y components of element force acting on nodes 1 and 2:
    %
    %  fe1x = fe*Dx/L
    %  fe1y = fe*Dy/L
    %  fe2x = -fe*Dx/L
    %  fe2y = -fe*Dy/L

    K(dof1x,dof1x) = K(dof1x,dof1x) + ke*(Dx/L)^2;
    K(dof1x,dof2x) = K(dof1x,dof2x) - ke*(Dx/L)^2;
    K(dof1x,dof1y) = K(dof1x,dof1y) + ke*(Dx/L)*(Dy/L);
    K(dof1x,dof2y) = K(dof1x,dof2y) - ke*(Dx/L)*(Dy/L);

    K(dof1y,dof1x) = K(dof1y,dof1x) + ke*(Dx/L)*(Dy/L);
    K(dof1y,dof2x) = K(dof1y,dof2x) - ke*(Dx/L)*(Dy/L);
    K(dof1y,dof1y) = K(dof1y,dof1y) + ke*(Dy/L)^2;
    K(dof1y,dof2y) = K(dof1y,dof2y) - ke*(Dy/L)^2;

    K(dof2x,dof1x) = K(dof2x,dof1x) - ke*(Dx/L)^2;
    K(dof2x,dof2x) = K(dof2x,dof2x) + ke*(Dx/L)^2;
    K(dof2x,dof1y) = K(dof2x,dof1y) - ke*(Dx/L)*(Dy/L);
    K(dof2x,dof2y) = K(dof2x,dof2y) + ke*(Dx/L)*(Dy/L);

    K(dof2y,dof1x) = K(dof2y,dof1x) - ke*(Dx/L)*(Dy/L);
    K(dof2y,dof2x) = K(dof2y,dof2x) + ke*(Dx/L)*(Dy/L);
    K(dof2y,dof1y) = K(dof2y,dof1y) - ke*(Dy/L)^2;
    K(dof2y,dof2y) = K(dof2y,dof2y) + ke*(Dy/L)^2;
end

% Loop over nodes and fill in entries of load vector:
for node=1:N_node
    F(dof_index(node,1)) = node_forces(node,1);
    F(dof_index(node,2)) = node_forces(node,2);
end

% Apply BCs to the system:
for node=1:N_node
    for direction=1:2
        if(BCs(node,direction))
            dof = dof_index(node,direction);
            for j=1:N_dof
                K(dof,j) = (dof==j);
            end
            F(dof) = 0;
        end
    end
end

% Solve for the displacements:
U = K\F;

% Plot the original and deformed trusses:
figure(1)
hold on;
for el=1:N_el
    n1 = elements(el,1);
    n2 = elements(el,2);
    X1 = nodes(n1,1);
    X2 = nodes(n2,1);
    Y1 = nodes(n1,2);
    Y2 = nodes(n2,2);
    dof1x = dof_index(n1,1);
    dof2x = dof_index(n2,1);
    dof1y = dof_index(n1,2);
    dof2y = dof_index(n2,2);
    u1x = U(dof1x);
    u2x = U(dof2x);
    u1y = U(dof1y);
    u2y = U(dof2y);
    x1 = X1 + u1x;
    x2 = X2 + u2x;
    y1 = Y1 + u1y;
    y2 = Y2 + u2y;
    plot([X1,X2],[Y1,Y2],':b');
    plot([x1,x2],[y1,y2],'-b');
end
