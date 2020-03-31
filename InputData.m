%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       INPUT   DATA       %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Element Properties
% Saved into a vector for multiple materials
% A(mat) = Gives desired A
% If (mat) = 1  ---> Mat 1
% If (mat) = 2  ---> Mat 2

m1.E = 210e3 ;         % [MPa] Young's modul
m1.d = 5;              % [mm]  Circular beam radius
m1.A = pi*(m1.d/2)^2;  % [mm^2] Beam section area
m1.rho = 2300;         % [kg/m^3] material density
m1.l = 1;              % [m] Hor bar longitude

% material matrix deifinition
input.mat = [% Young M.   Section A.     Density
                m1.E,      m1.A,         m1.rho; % Material (1)
                m1.E,      m1.A,         m1.rho; % Material (2)
            ]';
% 

input.Mind = [% Mat. index
               1; % Element (1)
               1; % Element (2)
        ]';



%% Nodal Points
input.position = [%  X      Y
                    0,      0;    % Node (1)
                 2*m1.l,   0;     % Node (2)
                    0,     -m1.l; % Node (3)
                  ]';

%% Connectivity Matrix T(i_elem,i_node)
% Node 1 = Initial beam node
% Node 2 = Final beam node
input.T = [%  Node1   Node2
              1,     2; % 1-2   - Element 1
              2,     3; % 2-3   - Element 2
           ]';
    
%% Boundary Conditions
% Fix nodes matrix creation
input.fixNod = [% Node      DOF  Magnitude
                  1,         1,         0;
                  1,         2,         0;
                  3,         1,         0;
                  3,         2,         0;
]';
                
% [N] External force
force.F1 = 1;
force.F2 = 0;
% DOF = 1 - X Axis
% DOF = 2 - Y Axis
input.Fext = [%   Node   DOF     Magnitude
                    2,     2,     -force.F1;        % Force 1
                    2,     1,     -force.F2;        % Force 2
             ]';
         