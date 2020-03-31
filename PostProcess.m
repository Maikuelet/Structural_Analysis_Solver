%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  POSTPROCESS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot deformed and undeformed shape of 2D structure

function PostProcess (input ,output,m1)


% Reorganizinf position and displacements vector  U in order to be able 
% to fit the conditions
input.position = input.position';
U_mod = zeros(size(input.position,1),size(input.position,2));
for i = 1:size(input.position,1) 
    j = 2*i-1;
    U_mod(i,:) = output.U(1,j:j+1);
end

%% Plot - 1
h1 = figure(1);
mFactor = 10000;  % Multiplying factor to see clearly the deformed structure
x = 'X - Axis';
y = 'Y - Axis';
titulo = 'Structure Representation';

output.U = U_mod*mFactor; 
d = input.position + output.U;

nodesX1 = input.position(input.T(: ,1) ,1)';
nodesX2 = input.position(input.T(: ,2) ,1)';
nodesY1 = input.position(input.T(: ,1) ,2)';
nodesY2 = input.position(input.T(: ,2), 2)';
plot ([ nodesX1  ; nodesX2]  ,....
      [ nodesY1  ; nodesY2 ] ,'k',...
      [ d(input.T(: ,1),1)' ; d(input.T(: ,2) ,1)' ] ,[ d( input.T(: ,1) ,2)' ;...
    d(input.T(: ,2) ,2)'] ,'r','LineWidth',2) ;

 
xlabel(x,'interpreter','latex' ), ylabel(y,'interpreter','latex' );
title ( titulo,'interpreter','latex' ) ;

legend({'Original','','Deformed'},'Location','northeastoutside')
%legend('Original','Deformed','Location','northeastoutside')

dimm = [0.25 0.1 0.3 0.1];
str = sprintf('Displacements multiplied by a factor of %0.2e',mFactor);
annotation('textbox',dimm,'String',str,'FitBoxToText','on');
set(gcf, 'Position',  [100, 100, 750, 500]);



%% Displays

disp('Stiffness Matrix:');
disp(output.KG);
disp('-------------------------------------------------------------');
disp('Reactions:');
disp(output.rf);
disp('-------------------------------------------------------------');
disp('Displacements:');
disp(output.U);
disp('-------------------------------------------------------------');
disp('Strain:');
disp(output.strain);
disp('-------------------------------------------------------------');
disp('Stress:');
disp(output.stress);

%% Plot Storage
pathh     = pwd;
myfolder = 'Plots';
f1 = fullfile(pathh , myfolder);
mkdir(f1);       % Creates Folder

f = fullfile(f1 , sprintf('Plots_Rho%d_E.png', m1.rho));
saveas(h1,f);



end
