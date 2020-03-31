%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  REACTION AND DISPALCEMENT  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = ReaccionsAndDisplacements(input,output,dim)

%% Global force vector assembly
% Force 1 Applyed

f = zeros(dim.Ndofs,1);

for i = 1:size(input.Fext,2)
% Vector Assembly
        % Node*2dim -2 + DOF
fposi = input.Fext(1,i)*2 -2 + input.Fext(2,i);
f(fposi,1) = input.Fext(3,i); 
end


%% Restricted degrees of freedom

% Vector conaining all the fixed degrees of freedom
NfixNod = size(input.fixNod,2);              
vR = zeros(1,NfixNod);
uR = zeros(NfixNod,1);      % Imposed zero displacement value
     
fixnodeposi = input.fixNod(1,:)*2-2 + input.fixNod(2,:);   % Node*2dim -2 + DOF
vR(1,fixnodeposi) = input.fixNod(1,:)*2-2 + input.fixNod(2,:);

vR = vR(vR~=0);     %Reducing the 0's of previously stored arrat

% Vector conaining all the free degrees of freedom
vL = setdiff(1:dim.Ndofs,vR);



%% Reduced system matrices
KLL = output.KG(vL,vL);
KLR = output.KG(vL,vR);
KRL = output.KG(vR,vL);
KRR = output.KG(vR,vR);
fextL = f(vL,1);
fextR = f(vR,1);

%% Solve linear system
uL = KLL\(fextL - KLR*uR);
Rd = KRR*uR + KRL*uL - fextR;

% Assembly of global displacements
output.U(vL) = uL;
output.U(vR) = uR;

% Assembly of global reactions
output.rf(vL) = 0;              % This is 0 because of the own value deffintion 
output.rf(vR) = Rd;

end 