% ########################
% ### 2D deforestation ###
% ########################

% Last version: Aug 12 2025

% ######################################################
% ### Deforestation Model with Agribusiness Dynamics ###
% ######################################################

%clear all; clc; close all;

addpath ..\00_functions\ %..\ is used to get back in the folders
addpath ..\02_inputs\ %..\ is used to get back in the folders

format long

i = 2; % MT
       % Legal Amazon
       % Legal Amazon, different params
if i ==1
    addpath ..\03_outputs\exercise_mt\
    % LOAD SPATIAL DATA
    protected_areas = readgeoraster("protected001.tif");
    productivity = readgeoraster("productivity001.tif");
    deforestation = readgeoraster("deforestation001.tif");

    [m, n] = size(productivity);

    params = struct('N',max([m,n]),...
                'beta',0.10,...
                'alpha',0.75,...
                'pi',0.5, ...
                'Tlag',3,...
                'Tmax',36,...
                'delta',1,...
                'tau',1.2,...
                'wage',1, ...
                'elasticity',1.05,...
                'Sbar',400,...
                'pI',1,...
                'rho',1.5); 

    mask = zeros(params.N);
    mask(1:m,1:n) = 1;
    
    A = zeros(params.N, params.N);                  % padded output (zeros elsewhere)
    A(1:m, 1:n) = productivity/1000;
    
    F = 3.45*ones(params.N);               % Homogeneous deforestation cost
                                          % target size
    L = zeros(params.N, params.N, 'uint8');                  % padded output (zeros elsewhere)
    L(1:m, 1:n) = protected_areas;

    D = zeros(params.N, params.N, params.Tmax, 'uint8');                  % padded output (zeros elsewhere)
    D(1:m, 1:n,1) = deforestation;
    
    w = params.wage * ones(params.N,params.N);            % Constant wage
        
    O = zeros(params.N, params.N, params.Tmax); O_aux = rand(params.N,params.N); 
    aux = O_aux; aux(~D(:,:,1)) = 0; O(:,:,1) = aux;
    Agribusiness = zeros(params.N, params.N, params.Tmax); Agribusiness(:,:,1) = D(:,:,1); % Agribusiness starts in one of the vertices
        
    ret = SolveDeforestationAMZ2(params,A,F,O,D,L,Agribusiness,O_aux,1,false,mask);
end

if i ==2
    addpath ..\03_outputs\exercise_amabiome\
    % LOAD SPATIAL DATA
    protected_areas = readgeoraster("conservation001.tif") + readgeoraster("indiginous001.tif");
    productivity = readgeoraster("productivity001.tif");
    deforestation = readgeoraster("deforestation001.tif");
    agri = readgeoraster("agribusiness001.tif");

    [m, n] = size(productivity);
    params = struct( 'm',m,...
                'n',n,...
                'N',max([m,n]),...
                'beta',0.10,...
                'alpha',0.75,...
                'pi',0.5, ...
                'Tlag',2,...
                'Tmax',3,...
                'delta',1,...
                'tau',1.2,...
                'wage',1, ...
                'elasticity',1.05,...
                'Sbar',400,...
                'pI',1,...
                'rho',1.5); 
    
    aux = readgeoraster("mask001.tif");
    mask = zeros(params.N, params.N);
    aux(isnan(aux)) = 0;
    mask(1:m,1:n) = aux;

    A = zeros(params.N, params.N, 'like', productivity);                  % padded output (zeros elsewhere)
    A(1:m, 1:n) = productivity/1000;
    
    F = 3.4*ones(params.N);               % Homogeneous deforestation cost
                                          % target size
    L = zeros(params.N, params.N, 'like', protected_areas);                  % padded output (zeros elsewhere)
    L(1:m, 1:n) = protected_areas;

    D = zeros(params.N, params.N, params.Tmax);                  % padded output (zeros elsewhere)
    D(1:m, 1:n,1) = deforestation;

    Agribusiness = zeros(params.N, params.N, params.Tmax);                  % padded output (zeros elsewhere)
    Agribusiness(1:m, 1:n,1) = agri;
    
    w = params.wage * ones(params.N,params.N);            % Constant wage
        
    O = zeros(params.N, params.N, params.Tmax); O_aux = rand(params.N,params.N); O(:,:,1) = D(:,:,1).*O_aux;
        
    ret = SolveDeforestationAMZ2(params,A,F,O,D,L,Agribusiness,O_aux,1,false,mask);
end

if i ==3
    addpath ..\03_outputs\exercise_amabiome\
    % LOAD SPATIAL DATA
    protected_areas = readgeoraster("conservation001.tif") + readgeoraster("indiginous001.tif");
    productivity = readgeoraster("productivity001.tif");
    deforestation = readgeoraster("deforestation001.tif");

    [m, n] = size(productivity);
    params = struct('N',max([m,n]),...
                'beta',0.10,...
                'alpha',0.75,...
                'pi',0.75, ...
                'Tlag',2,...
                'Tmax',36,...
                'delta',1,...
                'tau',1.2,...
                'wage',1, ...
                'elasticity',1.05,...
                'Sbar',400,...
                'pI',1,...
                'rho',1.5); 
    
    aux = readgeoraster("mask001.tif");
    mask = zeros(params.N, params.N);
    mask(1:m,1:n) = aux;

    A = zeros(params.N, params.N, 'like', productivity);                  % padded output (zeros elsewhere)
    A(1:m, 1:n) = productivity/1000;
    
    F = 3.8*ones(params.N);               % Homogeneous deforestation cost
                                          % target size
    L = zeros(params.N, params.N, 'like', protected_areas);                  % padded output (zeros elsewhere)
    L(1:m, 1:n) = protected_areas;

    D = zeros(params.N, params.N, params.Tmax);                  % padded output (zeros elsewhere)
    D(1:m, 1:n,1) = deforestation;
    
    w = params.wage * ones(params.N,params.N);            % Constant wage
        
    O = zeros(params.N, params.N, params.Tmax); O_aux = rand(params.N,params.N); O(:,:,1) = D(:,:,1).*O_aux;
    Agribusiness = zeros(params.N, params.N, params.Tmax); Agribusiness(:,:,1) = D(:,:,1); % Agribusiness starts in one of the vertices
        
    ret = SolveDeforestationAMZ2(params,A,F,O,D,L,Agribusiness,O_aux,1,true,mask);
end