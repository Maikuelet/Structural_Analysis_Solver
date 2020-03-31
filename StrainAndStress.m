%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%       STRAIN & STRESS      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = StrainAndStress(dim,output,input,T2)

output.strain = zeros (dim.Nelements);
output.stress = zeros (dim.Nelements);

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
    
    % Obtain element displacements
    for r=1:dim.Ndim*dim.NnodesXelement        
        p    = T2(r,e);        % Global degree of fredom
        dg(r) = output.U(p);          % Element displacements at global coordinates
    end
    
    R = [ ce se 0 0;
         -se ce 0 0;
          0  0 ce se;
          0  0 -se ce;
        ];
    % Calculate local element displacements
    dl =  R*dg';                     % Element displacements at local coordinates
    
    Ee = input.mat(1,input.Mind(1,e));

    % Calculate strain and stress
    output.strain(e) = (1 / le) * [-1 0 1 0]*dl;
    output.stress(e) = Ee*output.strain(e);
end

end