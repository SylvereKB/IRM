function [T,time,Theta0,Theta1,Theta2,r,rX,AA0T,AA1T,AA2T,AA0,AA1,...
    AA2,E0,E1,E2,A0,A1,A2,TotalMosqui]=MainFunctionIRM2(Ny,DurPeriod, ...
    VarM0,p2,p1,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
    gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
    RatioEpsi,InitialE0,InitialE1,InitialE2,InitialA0,InitialA1,InitialA2)

    %% Description of inputs
   % Ny= Number of years of insecticide deployment (integer)
   % DurPeriod= insecticide usage duration (integer)
   % VarM0= M0 variance (real)
   % EmergenceThreshold= Resistance emergence threshold (in percentage)
   % p= constant value of I1 exposure rate (real)
   % p2= overall exposure rate to insecticide (real)
   % p1= vector containing values of I1 exposure rates (vector)
   % p1+p2 <=1
   % Scenario= deployment strategy (vector)
   % tau= proportion of hatched eggs that reach adulthood (real)
   % age= discretisation of AFM age (vector)
   % IRLevel= discretisation of insect. resis. level (vector)
   % dx= IR level steps (real)  ; Nx= number of IR level steps (integer)
   % da= age step (real)        ; Na= number of age steps (integer)
   % RatioEpsi= the ratio of epsi1 compared to epsi0 (real)
   % gamma_i= hatching rate of eggs laid by A_i mosqui, i=0,1,2(vector)
   % mu_i= natural death rate of E_i eggs, i=0,1,2 (vector)
   % MinMosquiAgeLayEggs= min age for AFM to lay eggs in days (integer)
   % AvrMosquiLifeSpan= average mosquito lifespan in days (integer)
   % InitialE_i= initial values for E_i, i=0,1,2 (vector)
   % InitialA_i= initial values for A_i, i=0,1,2 (vector)
   % RatioEpsi= ratio of VarM1 relatively at VarM0 (real)
   % InitialE_i and InitialA_i= initial value (vector and matrix)
    %% Fixed parameters
    T= 300*Ny;                % time limit
    dt= 1;                    % time step
    VarM1= RatioEpsi*VarM0;   % M1 variance
    VarM2= RatioEpsi*VarM0;   % M2 variance
    k=0.5;                    % Logistic function exponent
    

         
    %% Discretisation of time
    time= 0:dt:T;             % discretisation of total time
    
    Nt= length(time);        % number of time steps
    
    %% insecticide exposure rate
    cp1= zeros(1,Nt); cp2= zeros(1,Nt);
    for j=1:(ceil(Ny/DurPeriod))
       Tdebut=(j-1)*DurPeriod*300; IdDebut=1+floor(Tdebut/dt);
       Tfin=min(Tdebut + DurPeriod * 300, T); IdFin=1+floor(Tfin/dt); 
       if j <= length(p1)
       cp1(IdDebut:IdFin) = p1(j); cp2(IdDebut:IdFin) = p2(j);
       else
       cp1(IdDebut:IdFin) = p1(end);  % Répète la dernière valeur de p1
       cp2(IdDebut:IdFin) = p2(end);
       end
     end
    
    %% Variable Parameters
    %mutation kernel of unexposed/exposed AFM
    FunM0= @(x)normpdf(x,0,VarM0);  %mutation function from IR level x to y 
    FunM1= @(x)normpdf(x,0,VarM1);  %mutation function from IR level x to y 
    FunM2= @(x)normpdf(x,0,VarM2);  %mutation function from IR level x to y 
    M0= FunM0(dx*(1-Nx:Nx-1));
    M1= FunM1(dx*(1-Nx:Nx-1));
    M2= FunM2(dx*(1-Nx:Nx-1));
    
    %% computation of r and Theta/0/1/2/12 (fitness function of unex./exposed AFM)
    [Theta0,Theta1,Theta2,r,rX,d0,d1,d2]=...
    FuncTheta(tau,age,IRLevel,Na,Nx,...
    Scenario,gamma0,gamma1,gamma2,...
    mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan);

    %% Computation of E0, E1, E2, A0, A1 and A2.
    %Initialization 
    E0= zeros(Nt,Nx);           A0= zeros(Na,Nx,Nt);
    E1= zeros(Nt,Nx);           A1= zeros(Na,Nx,Nt);
    E2= zeros(Nt,Nx);           A2= zeros(Na,Nx,Nt);
   
    E0(1,:)=InitialE0;     E1(1,:)=InitialE1;     E2(1,:)=InitialE2; 
    A0(:,:,1)=InitialA0;   A1(:,:,1)=InitialA1;   A2(:,:,1)=InitialA2; 

    for t=1:(Nt-1)
        %computation of H(E^n(t))
        A_Int=zeros(1,Nx);
        for x=1:Nx
            A_Int(x)=trapz(age,(A0(:,x,t)+A1(:,x,t)+A2(:,x,t))/tau(x));
        end
        % Et=trapz(IRLevel,E0(t,:)+E1(t,:)+E2(t,:)+A_Int*0);
        Et=trapz(IRLevel,E0(t,:)+E1(t,:)+E2(t,:));
        HEt=(1+Et)^(-k);
 
        %\int_Omega mi(x-y)[\int_0^\infty ri(a,x)Ai(a,x,t)da]dy (i=0,1,2,12)
        rA0 = trapz(age, r .* A0(:,:,t));
        rA1 = trapz(age, r .* A1(:,:,t));
        rA2 = trapz(age, r .* A2(:,:,t));      

        ConvoKernel0=dx*conv(rA0,M0,'same');  
        ConvoKernel1=dx*conv(rA1,M1,'same');  
        ConvoKernel2=dx*conv(rA2,M2,'same');
    
        %computation of E0, E1 and E2
        E0(t+1,:)=(E0(t,:)/dt + HEt*ConvoKernel0)./(1/dt + mu0 + gamma0);
        E1(t+1,:)=(E1(t,:)/dt + HEt*ConvoKernel1)./(1/dt + mu1 + gamma1);
        E2(t+1,:)=(E2(t,:)/dt + HEt*ConvoKernel2)./(1/dt + mu2 + gamma2);
    
        %computation of A0, A1 and A2
        A0(1,:,t+1)=(1-cp1(t)-cp2(t))*(tau.*gamma0.*E0(t+1,:)+...
            tau.*gamma1.*E1(t+1,:)+tau.*gamma2.*E2(t+1,:));
        A1(1,:,t+1)=cp1(t)*(tau.*gamma0.*E0(t+1,:) +...
        tau.*gamma1.*E1(t+1,:)+tau.*gamma2.*E2(t+1,:));
        A2(1,:,t+1)=cp2(t)*(tau.*gamma0.*E0(t+1,:) +...
            tau.*gamma2.*E2(t+1,:)+tau.*gamma1.*E1(t+1,:));
        
        for a=2:Na
            A0(a,:,t+1)=(A0(a,:,t)/dt + A0(a-1,:,t+1)/da)./...
                (1/dt + 1/da + d0(a,:));
            A1(a,:,t+1)=(A1(a,:,t)/dt + A1(a-1,:,t+1)/da)./...
                (1/dt + 1/da + d1(a,:));
            A2(a,:,t+1)=(A2(a,:,t)/dt + A2(a-1,:,t+1)/da)./...
                (1/dt + 1/da + d2(a,:));
        end
    end
    
    %% Computation of total pop of eggs and AFM at time t
    %size A0=(a,x,t). trapz(age,A0,1) return a matrix of size (1,x,t).
    %squeeze(trapz(age, A0, 1)) return a matrix of size (x,t) ie supression of first dim
    %squeeze(trapz(age, A0, 1))' return the transpose matrix of size (t,x)
    AA0= (squeeze(trapz(age, A0, 1)))';
    AA1= (squeeze(trapz(age, A1, 1)))';
    AA2= (squeeze(trapz(age, A2, 1)))';
       
    AA0T= squeeze(trapz(IRLevel, AA0, 2));
    AA1T= squeeze(trapz(IRLevel, AA1, 2));
    AA2T= squeeze(trapz(IRLevel, AA2, 2));
    
    TotalMosqui= trapz(time, AA0T+AA1T+AA2T);

end
