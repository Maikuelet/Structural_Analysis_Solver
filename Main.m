%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  STATIC STRUCTURE 2D SOLVER   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Miquel Altadill Llasat  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%     Julia Blanch Sánchez        %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot and save deformed and undeformed shape of 2D structure

clc
clear all
close all

%% =========================== INPUT DATA ============================== %%
InputData

%% ==========================  PREPROCESS =============================  %%
[dim, T2 ] = PreProcess(input);

%% ============================= SOLVER ===============================  %%
[output]  = StiffnessMatrix(dim,input,T2);
[output] = ReaccionsAndDisplacements(input,output,dim);
[output] = StrainAndStress(dim,output,input,T2);

%% ========================== POSTPROCESS =============================  %%
PostProcess (input ,output,m1);