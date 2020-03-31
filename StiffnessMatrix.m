%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  STIFFNESS MATRIX  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function in charge of generating the Global Stiffness Matrix


function [output] = StiffnessMatrix(dim,input,T2)

% Stiffness matrix for each bar element
Kel = zeros(dim.NnodesXelement*dim.Ndim,dim.NnodesXelement*dim.Ndim,dim.Nelements);

for e=1:dim.Nelements
    
    % Element Position
    x1e=input.position(1,input.T(1,e));
    x2e=input.position(1,input.T(2,e));
    y1e=input.position(2,input.T(1,e));
    y2e=input.position(2,input.T(2,e));
    
    % Element Lenght and orientation
    le=sqrt((x2e-x1e)^2 +(y2e-y1e)^2);
    se=((y2e-y1e)/le);                  % sine   of Element
    ce=((x2e-x1e)/le);                  % cosine of Element
    
    
    %Element Matrix 
    %       Ae * Ee   /   le
    Ee = input.mat(1,input.Mind(1,e));
    Ae = input.mat(2,input.Mind(1,e));
    ke=((Ae*Ee)/le)...
        *([ce^2 ce*se -ce^2 -ce*se;
           ce*se se^2 -ce*se -se^2;
           -ce^2 -ce*se ce^2 ce*se;
           -ce*se -se^2 ce*se se^2;
        ]);
    
    % Store Element Matrix
    for r=1:dim.Ndim*dim.NnodesXelement
        for s=1:dim.Ndim*dim.NnodesXelement
            Kel(r,s,e)=ke(r,s);
        end
    end
end

%  Global stiffnes Storage
output.KG = zeros(dim.Ndofs,dim.Ndofs);

for e=1:dim.Nelements                             % For each Beam
    for i=1:dim.Ndim*dim.NnodesXelement           % Local degree of freedom
        I=T2(i,e);                                % Global dof
        for j= 1:dim.Ndim*dim.NnodesXelement      % Local dof
            J=T2(j,e);                            % Global dog
            output.KG(I,J)=output.KG(I,J)+ Kel(i,j,e);
        end
    end
end

end