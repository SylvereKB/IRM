close all
clear

fpath = '/Users/kezetasylvere/Documents/CODES AOUT 2025';     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;  % Line width
LETTERS = 'A':'Z';
% Colors
color = [0.89 0.26 0.20;...      % vermillion (red)
         0.255 0.411 0.882;...   % royal blue 
         0.291 0 0.517;...       % indigo
         0.133 0.698 0.133];     % forest green

%% Set of fixed parameters
VarM0 = 0.01/5;             % M0 variance
A = 30;                     % mosquito age limit
AvrMosquiLifeSpan = 21;     % average mosquito lifespan in days
MinMosquiAgeLayEggs = 3;    % min age for AFM to lay eggs in days
da = 1;                     % age step
xmin = 0;                   % min value of IR level
xmax = 1.1;                 % max value of IR level
dx = VarM0/2;               % IR level step (dx<VarM0)
age = 0:da:A;               % discretisation of AFM age
IRLevel = xmin:dx:xmax;     % discretisation of insect. resis. level
Na = length(age);           % number of age steps
Nx = length(IRLevel);       % number of insect. resist. level steps

% Natural death rate of eggs laid by unexposed/exposed AFM
mu0 = 0.2*ones(1,Nx);       mu1 = mu0;    
mu2 = mu0;                        

% Hatching rate of eggs laid by unexposed/exposed AFM
gamma0 = 0.75*ones(1,Nx);   gamma1 = gamma0;    
gamma2 = gamma0;                   

% Proportion of hatched eggs laid by AFM that reach adulthood
tau = 0.67*ones(1,Nx);

% Key parameters that determine array sizes
Ny = 12;                    % Total years
DurPeriod = 3;              % Duration of each period    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EmergenceThreshold = 10;    % in percentage
RatioEpsi = 2;

% Calculate dependent constants
nb_cycles = ceil(Ny/DurPeriod);  % Number of cycles based on Ny and DurPeriod

%% SECTION: Mos vs Rot vs MosOpti - OPTIMIZED VERSION
if(1)
    format long
    tic
    
    % Parameter vectors - modify these as needed
    VectStrategy = 3;                    % lStra=length(VectStrategy); %Mos=1, Rot=2, MosOpti=3
    VectVar0 = [0.01, 0.01/5, 0.01/10];  % lVar0=length(VectVar0);
    VectScenario = 2:4;                  % lScena=length(VectScenario);
    VectDurCycle = 3;                    % lDurCycle=length(VectDurCycle);          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VectP0 = 0.05:0.1:0.15;                       % lP0=length(VectP0);
    
    % Calculate lengths
    lStra = length(VectStrategy);
    lVar0 = length(VectVar0);
    lScena = length(VectScenario);
    lDurCycle = length(VectDurCycle);
    lP0 = length(VectP0);
    
    % Total number of rows in experiment plan
    n_rows = lStra * lVar0 * lScena * lDurCycle * lP0;
    
    % OPTIMIZED: Use ndgrid to generate all combinations efficiently
    fprintf('Generating experiment plan using ndgrid...\n');
    [G1, G2, G3, G4, G5] = ndgrid(VectStrategy, VectVar0, VectScenario, VectDurCycle, VectP0);
    
    % Reshape to create experiment plan matrix - much faster than nested loops
    plan_exp = [G1(:), G2(:), G3(:), G4(:), G5(:)];
    
    % Verify we have the correct number of rows
    assert(size(plan_exp, 1) == n_rows, 'Mismatch in experiment plan size');
    
    % Pre-allocate arrays for INPUT storage
    Strategy = plan_exp(:, 1);
    VM0 = plan_exp(:, 2); 
    REFF = plan_exp(:, 3);
    DurCycle = plan_exp(:, 4);
    P0 = plan_exp(:, 5);

    % Pre-allocate arrays for OUTPUT storage (using calculated nb_cycles)
    PropI1 = zeros(n_rows, nb_cycles);
    PropI2 = zeros(n_rows, nb_cycles);
    Rgain = zeros(n_rows, 1);
    Temerg = zeros(n_rows, 1); 
    
    % Parallel processing setup
    numWorkers = 6;                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('Using parallel processing with %d workers...\n', numWorkers);
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        if numWorkers > 1
            parpool(numWorkers);
            fprintf('Parallel pool created with %d workers\n', numWorkers);
        else
            fprintf('Running in serial mode (1 worker)\n');
        end
    else
        fprintf('Parallel pool already exists with %d workers\n', poolobj.NumWorkers);
    end

    % NOUVELLE SECTION DE PROGRESSION
    fprintf('\n=== DEBUT DU CALCUL ===\n');
    fprintf('Total de combinaisons de parametres: %d\n', n_rows);
    fprintf('Format d''affichage: run actuel/nombre total (pourcentage)\n');
    fprintf('Debut du traitement...\n\n');
    
    % Pre-allocate temporary cell arrays for parfor compatibility
    temp_PropI1 = cell(n_rows, 1);
    temp_PropI2 = cell(n_rows, 1);
    
    % Main computation loop avec nouveau systeme de progression
    parfor nl = 1:n_rows 
        param = plan_exp(nl, :);
        strategy = param(1);
        VarM0_local = param(2);
        Scenario = param(3);
        DurPeriod_local = param(4);
        p0 = param(5);
        
        % Calculate local nb_cycles for this iteration
        nb_cycles_local = ceil(Ny/DurPeriod_local);
        
        [p1, p2, RelatGain, Temergence] = FunctionMos_Rot_MosOpti(strategy, VarM0_local, Scenario, p0, ...
            Ny, DurPeriod_local, tau, age, IRLevel, dx, Nx, da, Na, gamma0, gamma1, gamma2, mu0, mu1, mu2, ...
            MinMosquiAgeLayEggs, AvrMosquiLifeSpan, RatioEpsi, EmergenceThreshold, nb_cycles_local);
        
        % Store scalar results
        Rgain(nl) = RelatGain;
        Temerg(nl) = Temergence;
        
        % Store variable-length results in cell arrays
        temp_PropI1{nl} = p1;
        temp_PropI2{nl} = p2;
        
        % AFFICHAGE DE LA PROGRESSION - Format: run actuel/total (%)
        fprintf('Progress: %d/%d (%.1f%%) - Strategy:%d, Var0:%.5f, Scenario:%d\n', ...
                nl, n_rows, 100*nl/n_rows, strategy, VarM0_local, Scenario);
    end
    
    fprintf('\n=== TRAITEMENT TERMINE ===\n');
    fprintf('Post-traitement des resultats en cours...\n');
    
    % After parfor loop, convert cell arrays to matrices
    for nl = 1:n_rows
        p1_temp = temp_PropI1{nl};
        p2_temp = temp_PropI2{nl};
        
        % Handle variable-length outputs safely
        if length(p1_temp) <= nb_cycles
            PropI1(nl, 1:length(p1_temp)) = p1_temp;
            if length(p1_temp) < nb_cycles
                PropI1(nl, length(p1_temp)+1:end) = NaN;  % Fill remaining with NaN
            end
        else
            PropI1(nl, :) = p1_temp(1:nb_cycles);
        end
        
        if length(p2_temp) <= nb_cycles
            PropI2(nl, 1:length(p2_temp)) = p2_temp;
            if length(p2_temp) < nb_cycles
                PropI2(nl, length(p2_temp)+1:end) = NaN;  % Fill remaining with NaN
            end
        else
            PropI2(nl, :) = p2_temp(1:nb_cycles);
        end
    end
    
    % Output table name
    tabName = sprintf('TabForMosRotMosOpti_Strategy3_DurPeriod%d05.xlsx', DurPeriod); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate output table
    fprintf('Generation du tableau de resultats...\n');
    Tab = table(Strategy, VM0, REFF, DurCycle, P0, PropI1, PropI2, Rgain, Temerg);  
    
    fprintf('Ecriture des resultats dans: %s...\n', tabName);
    writetable(Tab, fullfile(fpath, tabName));
    
    fprintf('\n=== CALCUL TERMINE AVEC SUCCES ===\n');
    fprintf('Fichier de sortie: %s\n', tabName);
    fprintf('Nombre total de runs traites: %d/%d\n', n_rows, n_rows);
    toc
end

%% Intermediate functions

function [RelatMosqui] = FunctionMainOpti(p1_and_p2, Ny, DurPeriod, VarM0, Scenario, ...
    tau, age, IRLevel, dx, Nx, da, Na, gamma0, gamma1, gamma2, mu0, mu1, mu2, MinMosquiAgeLayEggs, ...
    AvrMosquiLifeSpan, RatioEpsi, InitialE0, InitialE1, InitialE2, ...
    InitialA0, InitialA1, InitialA2, TotalMosquiEqui)

    SizePeriod = ceil(Ny/DurPeriod);
    nn0 = SizePeriod;  % More explicit sizing
    
    p1 = p1_and_p2(1:nn0); 
    p2 = p1_and_p2(nn0+1:2*nn0);

    [T, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
        ~, ~, ~, ~, ~, ~, ~, TotalMosqui] = MainFunctionIRM2(Ny, DurPeriod, ...
        VarM0, p2, p1, Scenario, tau, age, IRLevel, dx, Nx, da, Na, ...
        gamma0, gamma1, gamma2, mu0, mu1, mu2, MinMosquiAgeLayEggs, AvrMosquiLifeSpan, ...
        RatioEpsi, InitialE0, InitialE1, InitialE2, InitialA0, InitialA1, InitialA2);

    RelatMosqui = TotalMosqui / (T * TotalMosquiEqui);
end

function [p1, p2, RelatGain, Temergence] = FunctionMos_Rot_MosOpti(strategy, VarM0, Scenario, p0, ...
    Ny, DurPeriod, tau, age, IRLevel, dx, Nx, da, Na, gamma0, gamma1, gamma2, mu0, mu1, mu2, ...
    MinMosquiAgeLayEggs, AvrMosquiLifeSpan, RatioEpsi, EmergenceThreshold, nb_cycles_local)

    SizePeriod = nb_cycles_local;  % Use passed parameter for consistency
 
    if strategy == 1  % Mosaic
        p1 = 0.5 * ones(1, SizePeriod);
        p2 = p1;
    elseif strategy == 2  % Rotation
        p2 = mod(1:SizePeriod, 2);     % [1,0,1,0,...] for SizePeriod=4: [1,0,1,0]
        p1 = 1 - p2;                    % [0,1,0,1,...] for SizePeriod=4: [0,1,0,1]
    elseif strategy == 3  % MosOpti
        [Theta0, Theta1, Theta2, ~, ~, ~, ~, ~] = ...
            FuncTheta(tau, age, IRLevel, Na, Nx, Scenario, gamma0, gamma1, gamma2, ...
            mu0, mu1, mu2, MinMosquiAgeLayEggs, AvrMosquiLifeSpan);
    
        % For p0=p1 fixed, equilibrate system without p2 
        NyEqui = 5;
        p1equi = p0 * ones(1, SizePeriod); 
        p2equi = 0 * p1equi;
        
        % Find the maximum value of Theta1/Theta2 and its position
        Theta = (1 - p1equi(1) - p2equi(1)) * Theta0 + p1equi(1) * Theta1 + p2equi(1) * Theta2;
        [~, maxId] = max(Theta1); max1 = IRLevel(maxId); 
        [~, maxId] = max(Theta2); max2 = IRLevel(maxId); IdMaxTheta2 = maxId;
        [~, maxId] = max(Theta); maxGlobal = IRLevel(maxId); IdMaxGlobal = maxId;

        % Start close enough to equilibrium to reach it in Ny=5 years
        InitialE0 = 1 * normpdf(IRLevel, maxGlobal, VarM0);  
        InitialE1 = 1 * normpdf(IRLevel, maxGlobal, VarM0); 
        InitialE2 = zeros(1, Nx); 
        InitialA0 = zeros(Na, Nx); 
        InitialA1 = zeros(Na, Nx); 
        InitialA2 = zeros(Na, Nx);

        [~, ~, ~, ~, ~, ~, ~, AA0T, AA1T, AA2T, ~, ~, ...
            ~, E0, E1, E2, A0, A1, A2, ~] = MainFunctionIRM2(NyEqui, DurPeriod, ...
            VarM0, p2equi, p1equi, Scenario, tau, age, IRLevel, dx, Nx, da, Na, gamma0, ...
            gamma1, gamma2, mu0, mu1, mu2, MinMosquiAgeLayEggs, AvrMosquiLifeSpan, ...
            RatioEpsi, InitialE0, InitialE1, InitialE2, InitialA0, InitialA1, InitialA2);

        % Calculate total mosquito count at equilibrium   
        TotalMosquiEqui = AA0T(end) + AA1T(end) + AA2T(end);

        % Retrieve system equilibrium without I2
        InitialE0 = E0(end, :); InitialA0 = A0(:, :, end); 
        InitialE1 = E1(end, :); InitialA1 = A1(:, :, end);
        InitialE2 = E2(end, :); InitialA2 = A2(:, :, end);
        
        % Objective function for p1_and_p2 variables only
        objFun = @(p1_and_p2) FunctionMainOpti(p1_and_p2, Ny, DurPeriod, VarM0, Scenario, ...
        tau, age, IRLevel, dx, Nx, da, Na, gamma0, gamma1, gamma2, mu0, mu1, mu2, MinMosquiAgeLayEggs, ...
        AvrMosquiLifeSpan, RatioEpsi, InitialE0, InitialE1, InitialE2, ...
        InitialA0, InitialA1, InitialA2, TotalMosquiEqui);
        
        % Solving the system with Optimal strategy
        Ops = optimset('Display', 'iter', 'Algorithm', 'interior-point', 'MaxFunEvals', 10^5, ...
                      'MaxIter', 2000, 'TolX', 0.0001, 'TolFun', 0.0001);                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Constraint matrix for optimization (p1 + p2 <= 1 for each period)
        AOpti = [eye(SizePeriod),eye(SizePeriod)];
        bOpti = ones(SizePeriod, 1);
        
        % Initial guess
        InitP1 = 0.8*ones(1, SizePeriod);
        InitP2 = 1 - InitP1;
        p1_and_p2_Init = [InitP1, InitP2];
        
        % Bounds
        lb = zeros(1, 2*SizePeriod);
        ub = ones(1, 2*SizePeriod);
        
        tic
        Prop = fmincon(objFun, p1_and_p2_Init, AOpti, bOpti, [], [], lb, ub, [], Ops);
        toc
        
        p1 = Prop(1:SizePeriod);
        p2 = Prop(SizePeriod+1:2*SizePeriod);
    end
    
    [~, ~, ~, ~, ~, ~, ...
        ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, RelatGain, Temergence] = FunctionMainRun(Ny, DurPeriod, ...
        VarM0, EmergenceThreshold, p2, p1, p0, Scenario, tau, age, IRLevel, dx, Nx, da, Na, ...
        gamma0, gamma1, gamma2, mu0, mu1, mu2, MinMosquiAgeLayEggs, AvrMosquiLifeSpan, ...
        RatioEpsi);
end