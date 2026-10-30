function[Theta0,Theta1,Theta2,AA0T,AA1T,AA2T,...
    AA0,AA1,AA2,E0,E1,E2,A0,A1,A2,time,T,RelatGain,Temergence,PopAdat2] = FunctionMainRun(Ny,DurPeriod, ...
    VarM0,EmergenceThreshold,p2,p1,p0,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
    gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
    RatioEpsi)
   
    %% computation of fitness functions
    [Theta0,Theta1,Theta2,~,~,~,~,~]=...
            FuncTheta(tau,age,IRLevel,Na,Nx,Scenario,gamma0,gamma1,gamma2,...
            mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan);
        
        %% we set the system at equilibrium without I2
        NyEqui=5;
        p1equi =[p0,p0,p0,p0]; p2equi =0*p1equi;
        
        % find the maximum value of Theta1/Theta2 and its position
        Theta=(1-p1equi(1)-p2equi(1))*Theta0+ p1equi(1)*Theta1+p2equi(1)*Theta2;
        
        [maxTheta1,maxId]= max(Theta1); max1= IRLevel(maxId); 
        [maxTheta2,maxId]= max(Theta2); max2= IRLevel(maxId); IdMaxTheta2=maxId;
        [maxTheta,maxId]= max(Theta); maxGlobal= IRLevel(maxId); IdMaxGlobal=maxId;
        
        % we start sufficiently near the previous equilibrium to reach it
        % in Ny=5 ans
        InitialE0=1*normpdf(IRLevel,maxGlobal,VarM0);  
        InitialE1= 1*normpdf(IRLevel,maxGlobal,VarM0); 
        InitialE2= zeros(1,Nx); 
        InitialA0= zeros(Na,Nx); InitialA1= zeros(Na,Nx); InitialA2= zeros(Na,Nx);
        
        [~,time,~,~,~,~,~,AA0T,AA1T,AA2T,~,~,...
            ~,E0,E1,E2,A0,A1,A2,~]=MainFunctionIRM2(NyEqui,DurPeriod, ...
            VarM0,p2equi,p1equi,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
            gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
            RatioEpsi,InitialE0,InitialE1,InitialE2,InitialA0,InitialA1,InitialA2);
        
        % we compute the total number of mosquitoes at equlibrium   
        TotalMosquiEqui= AA0T(length(time)) + AA1T(length(time)) + ...
            AA2T(length(time));
        
        % we set initial conditions at equilibrium without I2
        InitialE0=E0(length(time),:); InitialA0= A0(:,:,length(time)); 
        InitialE1= E1(length(time),:); InitialA1= A1(:,:,length(time));
        InitialE2= E2(length(time),:); InitialA2=A2(:,:,length(time));
        
        %% we intriduce p2
        
        [T,time,~,~,~,~,~,AA0T,AA1T,AA2T,AA0,AA1,...
            AA2,E0,E1,E2,A0,A1,A2,TotalMosqui]=MainFunctionIRM2(Ny,DurPeriod, ...
            VarM0,p2,p1,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
            gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
            RatioEpsi,InitialE0,InitialE1,InitialE2,InitialA0,InitialA1,InitialA2);
        
        
        %% Computation of Temergence and retaive gain 
        RelatGain=1-TotalMosqui/(T*TotalMosquiEqui);
        
        Group1=AA0(:,IdMaxGlobal)+AA1(:,IdMaxGlobal)+AA2(:,IdMaxGlobal);
        Group2=AA0(:,IdMaxTheta2)+AA1(:,IdMaxTheta2)+AA2(:,IdMaxTheta2);
        PopAdat2=100*Group2./(Group1+Group2);
        
        Row_Temg= find(PopAdat2>EmergenceThreshold);
        if isempty(Row_Temg)  
           id_Temg= length(time);
        else 
            GRAD= gradient(Row_Temg);
            id_Temg= Row_Temg(1);
            for np=2:length(Row_Temg)
                if GRAD(np)>1, id_Temg= Row_Temg(np); 
                end
            end
        end
        Temergence=time(id_Temg);
end