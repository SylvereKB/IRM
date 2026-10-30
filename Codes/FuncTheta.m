function[Theta0,Theta1,Theta2,r,rX,d0,d1,d2]=...
    FuncTheta(tau,age,IRLevel,Na,Nx,...
    Scenario,gamma0,gamma1,gamma2,...
    mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan)

%% computation of Theta0/1/2/12 (fitness function of unex./exposed AFM) 

    % AFM laying eggs rate
    rrm= 200;          % maximum number of eggs laid by sAFM
    rr0= 0.875*rrm;    % average number of eggs laid by sAFM
    rr1= 0.750*rrm;    % average number of eggs laid by rAFM

% mortality rate of unexposed/exposed AFM
        if Scenario==1    %Theta2(maxTheta1) 25% more efficient than Theta1(maxTheta1)
            ProbaSurviFullySensi1= 10^-10;
            RatioPs2_sur_Ps1=10^4.575;
        end
        if Scenario==2    %Theta2(maxTheta1) 50% more efficient than Theta1(maxTheta1)
            ProbaSurviFullySensi1= 10^-10;
            RatioPs2_sur_Ps1=10^5.67;
        end
        if Scenario==3    %Theta2(maxTheta1) 75% more efficient than Theta1(maxTheta1)
            ProbaSurviFullySensi1= 10^-10;
            RatioPs2_sur_Ps1=10^6.637;
        end
        if Scenario==4    %Theta2(maxTheta1) 90% more efficient than Theta1(maxTheta1)
            ProbaSurviFullySensi1= 10^-10;
            RatioPs2_sur_Ps1=10^7.48;
        end
    
    ProbaSurviFullySensi2=RatioPs2_sur_Ps1*ProbaSurviFullySensi1; 
    dd01=-log(ProbaSurviFullySensi1);  %death rate of sAFM exposed to I1 
    dd02=-log(ProbaSurviFullySensi2);  
    dd1=-log(1-ProbaSurviFullySensi1); %death rate of rAFM exposed to I1
    dd2=-log(1-ProbaSurviFullySensi2); 
    NaturalDeath= 10*(age>AvrMosquiLifeSpan);  
    d0= zeros(Na,Nx);                 d1= zeros(Na,Nx);   
    d2= zeros(Na,Nx);                 
    DeathInsecticide1= zeros(1,Nx);
    DeathInsecticide2= zeros(1,Nx);
    for s=1:Nx
        DeathInsecticide1(s)=dd01*(dd1/dd01)^(IRLevel(s)); 
        DeathInsecticide2(s)=dd02*(dd2/dd02)^(IRLevel(s));
        d0(:,s)= NaturalDeath(:);
        d1(:,s)= NaturalDeath(:)+DeathInsecticide1(s);
        d2(:,s)= NaturalDeath(:)+DeathInsecticide2(s);
    end 

    %computation of survival probability of unexposed/exposed AFM
    % Calcul des intégrales cumulées de mortalité
    CumD0 = cumtrapz(age, d0);      % Na x Nx
    CumD1 = cumtrapz(age, d1);
    CumD2 = cumtrapz(age, d2);

    % Probabilité de survie = exp(-∫ mortalité)
    SurvProb0 = exp(-CumD0);        % Na x Nx
    SurvProb1 = exp(-CumD1);
    SurvProb2 = exp(-CumD2);

    
    % laying rate of unexposed/exposed AFM
    Barr= age>MinMosquiAgeLayEggs;
    r= zeros(Na,Nx);
    for s=1:Nx
        r(:,s)= rrm/(1+(rrm/rr0-1)*...
            (rr0*(rrm-rr1)/(rr1*(rrm-rr0)))^IRLevel(s))*Barr(:);
    end  
    
    % egg-laying rate as function of insecticide resistance level
    rX=squeeze(trapz(age,r,1));
    
    %computation of Theta0/Theta1/Theta2/Theta12 
    Theta0= tau.*gamma0./(mu0+gamma0).*...
            squeeze(trapz(age,r.*SurvProb0,1)); 
    Theta1= tau.*gamma1./(mu1+gamma1).*...
            squeeze(trapz(age,r.*SurvProb1,1)); 
    Theta2= tau.*gamma2./(mu2+gamma2).*...
            squeeze(trapz(age,r.*SurvProb2,1)); 
    
end
