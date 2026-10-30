close all
clear

fpath = '/Users/kezetasylvere/Documents/CODES AOUT 2025';


LW= 2;  % Line width
LETTERS='A':'Z';
%Couleurs
color= [ 0.133 0.698 0.133;...  %forest green
         0.890 0.260 0.200;...  %vermillon (rouge)
         0.255 0.411 0.882;...  %royal blue 
         0.291 0.000 0.517];    %indigo

%% Set of fixed parameters
VarM0= 0.01/5;             % M0 variance
A= 30;                   % mosquito age limit
AvrMosquiLifeSpan=21;    % average mosquito lifespan in days
MinMosquiAgeLayEggs=3;   % min age for AFM to lay eggs in days
da= 1;                  % age step
xmin= 0;                % min value of IR level
xmax= 1.1;               % max value of IR level
dx= VarM0/2;             % IR level step (dx<VarM0)
age= 0:da:A;              % discretisation of AFM age
IRLevel= xmin:dx:xmax;    % discretisation of insect. resis. level
Na= length(age);         % number of age steps
Nx= length(IRLevel);     % number of insect. resist. level steps

% natural death rate of eggs laid by unexposed/exposed AFM
mu0= 0.2*ones(1,Nx);             mu1= mu0;    
mu2= mu0;                        

% hatching rate of eggs laid by unexposed/exposed AFM
gamma0= 0.75*ones(1,Nx);          gamma1= gamma0;    
gamma2= gamma0;                   

% proportion of hatched eggs laid by AFM that reach adulthood
tau= 0.67*ones(1,Nx);

Ny=12;
DurPeriod=3;
nb_cycles=ceil(Ny/DurPeriod);
EmergenceThreshold=10; % in percentage
RatioEpsi=2;

%% SECTION 1: TAB FOR AS
if(0)
    
    format long
    tic

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vect_VarM0=VarM0; l_VarM0=length(Vect_VarM0); % on est un peu obligé de le faire pour chaque var0 dut fait de dx=var0/2
    Vect_RatioEpsi=[2,3]; l_RatioEpsi=length(Vect_RatioEpsi);
    Vect_p0=[0.1,0.25,0.33,0.5,0.67,0.8];    l_p0=length(Vect_p0);
    Vect_p1val=[0.1,0.25,0.5,0.75,0.9];  l_p1val=length(Vect_p1val);
    Vect_Scenario=1:4;        l_Scenario=length(Vect_Scenario);

    n_rows=l_RatioEpsi*l_p0*l_p1val*l_Scenario*l_VarM0;


    %STRORAGE OF INPUT
    RatEpsi=zeros(n_rows,1);
    P0=zeros(n_rows,1); 
    P1val=zeros(n_rows,1); 
    Scena=zeros(n_rows,1);
    VM0=zeros(n_rows,1);

    %STRORAGE OF OUPUT
    Rgain=zeros(n_rows,1);
    Temerg=zeros(n_rows,1); 
    Temerg2=zeros(n_rows,1); 


    %BEGIN EXPERIENCE_PLAN
    row_number=0;
    plan_exp=-111*ones(n_rows,5);

    for l0=1:l_RatioEpsi
        for l1=1:l_p0
            for l2=1:l_p1val
                for l3=1:l_Scenario
                    for l4=1:l_VarM0
                        row_number=row_number+1;
                        plan_exp(row_number,1)=Vect_RatioEpsi(l0);
                        plan_exp(row_number,2)=Vect_p0(l1);
                        plan_exp(row_number,3)=Vect_p1val(l2);     
                        plan_exp(row_number,4)=Vect_Scenario(l3);
                        plan_exp(row_number,5)=Vect_VarM0(l4);
                    end
                end
            end
        end
    end
    
    %BEGIN computation on experience plan
    fprintf('Progress:\n');
    fprintf(['\n' repmat('.',1,n_rows ) '\n\n']);
    
    parfor nl=1:n_rows 
       fprintf('\b|\n'); 

       param=plan_exp(nl,:);
       RatioEpsi=param(1);
       p0=param(2);
       p1val=param(3);
       Scenario=param(4);
       
       p1=ones(1,ceil(Ny/DurPeriod))*p1val*p0;
       p2=ones(1,ceil(Ny/DurPeriod))*(1-p1val)*p0;
       
       [Theta0,Theta1,Theta2,AA0T,AA1T,AA2T,...
        AA0,AA1,AA2,E0,E1,E2,A0,A1,A2,time,T,RelatGain,Temergence] = FunctionMainRun(Ny,DurPeriod, ...
        VarM0,EmergenceThreshold,p2,p1,p0,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
        gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
        RatioEpsi);
    
        RatEpsi(nl)=RatioEpsi;
        P0(nl)=p0; 
        P1val(nl)=p1val; 
        Scena(nl)=Scenario;
        VM0(nl)=VarM0;
    
        Rgain(nl)=RelatGain;
        Temerg(nl)=Temergence/300; 

    end
    
    %nom du tableau des sorties
    tabName='TabForASEpsiTresPetit.xlsx';
    
    %Ici on gènère le contenu du tableau des simus
     Tab=table(RatEpsi,P0,P1val,Scena,VM0,Rgain,Temerg);  
     
     writetable(Tab,fullfile(fpath,tabName))
     
     toc

end


%% SECTION 2: Initial conditions prior to I2 deployment=f(Scenario,RatioEpsi,p0)
if(0)      
    %parameters setting
    Scenario=3;                     RatioEpsi=2;                  p0=0.5;
    NyEqui= 5;     p1equi= p0*ones(1,nb_cycles);        p2equi= 0*p1equi;
    %fitness functions calculation
    [Theta0,Theta1,Theta2,~,~,~,~,~]=FuncTheta(tau,age,IRLevel,Na,Nx,Scenario, ...
        gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan);
         
    %find the maximum value of Gamma1/Gamma2 and its position
    ThetaEqui= (1-p1equi(1)-p2equi(1))*Theta0 + p1equi(1)*Theta1 + p2equi(1)*Theta2;
    
    [~,maxId]= max(Theta1); max1= IRLevel(maxId); IdMaxGamma1=maxId;
    [~,maxId]= max(Theta2); max2= IRLevel(maxId); IdMaxGamma2=maxId;
    [~,maxId]= max(ThetaEqui); maxEqui= IRLevel(maxId); IdMaxEqui=maxId;
    
    %start close enough to this equilibrium to reach it in NyEqui=5 years
    InitialE0= 1*normpdf(IRLevel,maxEqui,VarM0); InitialA0= zeros(Na,Nx);
    InitialE1= 1*normpdf(IRLevel,maxEqui,VarM0); InitialA1= zeros(Na,Nx);
    InitialE2= zeros(1,Nx);                      InitialA2= zeros(Na,Nx);
    
    %run MainFunctionIRM2
    [~,time,~,~,~,~,~,AA0T,AA1T,AA2T,AA0,AA1,~,E0,E1,E2,A0,A1,A2,~]= MainFunctionIRM2(Ny,DurPeriod, ...
    VarM0,p2equi,p1equi,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
    gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
    RatioEpsi,InitialE0,InitialE1,InitialE2,InitialA0,InitialA1,InitialA2);

    %calculate the total number of mosquitoes at equilibrium  
    TotalMosquiEqui= AA0T(end) + AA1T(end) + AA2T(end);
    
    %we get the equilibrium of the system in the absence of I2
    InitialE0= E0(end,:); InitialA0= A0(:,:,end); 
    InitialE1= E1(end,:); InitialA1= A1(:,:,end);
    InitialE2= E2(end,:); InitialA2= A2(:,:,end);

    %create graphics: Evolutionary dynamics under exclusive use of I1 at a proportion p_0=0.5
    figure;
    set(gcf, 'Position', [100, 100, 900, 360]);
    %fitness functions
    GrNum= 1;
    subplot('Position', [0.08, 0.25, 0.25, 0.68]);  
    plot(IRLevel,Theta0,'LineWidth',LW,'Color',color(1,:),'LineStyle','--');
    xlim([0 1]);
    hold on
    plot(IRLevel,Theta1,'LineWidth',LW,'Color',color(2,:),'LineStyle','-');
    plot(IRLevel,ThetaEqui,'LineWidth',LW,'Color',color(4,:),'LineStyle','-.');
    plot([max1 max1], [0 max(Theta1)], 'Color',color(2,:), 'LineStyle', '--');
    plot([maxEqui maxEqui], [0 max(ThetaEqui)], 'Color',color(4,:), 'LineStyle', '--');
    
    % Marques sur l'axe x pour maxTheta1 et maxTheta
    ylim_current = ylim;
    plot(max1, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
    color(2,:),'MarkerEdgeColor',color(2,:));
    plot(maxEqui, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
    color(4,:),'MarkerEdgeColor',color(4,:));
    
    % Texte sur les marques
    text(max1+0.02,ylim_current(1)+0.02*diff(ylim_current),'$x_1$','Interpreter','latex', ...
    'HorizontalAlignment','left','FontSize',15,'Color',color(2,:));
    text(maxEqui-0.02,ylim_current(1)+0.02*diff(ylim_current),'$\bar{x}_1$','Interpreter','latex', ...
    'HorizontalAlignment','right','FontSize',15,'Color',color(4,:));
    hold off
    
    % add labels and formatting
    xlabel('Phenotype ($x$)', 'Interpreter', 'latex', 'FontSize',18);
    ylabel('Fitness', 'Interpreter', 'latex', 'FontSize',18);
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);
    legend('$\Gamma_0$','$\Gamma_1$','$\Theta_{[p_0,0]}$','Interpreter', ...
        'latex','location','best','Orientation','vertical','fontsize',15);
    legend boxoff
    
    GrNum= GrNum+1;
    Tot_PopAtEqui=AA0(end,:)+AA1(end,:);
    
    subplot('Position', [0.38, 0.25, 0.25, 0.68]);  % Deuxième subplot
    plot(IRLevel,Tot_PopAtEqui,'LineWidth',LW,'Color',color(4,:),'LineStyle','-');
    xlim([0 1]);
    hold on
    plot([maxEqui maxEqui], [0 max(Tot_PopAtEqui)],'Color',color(4,:),'LineStyle','--');
    
    % Marque sur l'axe x pour maxTheta
    ylim_current = ylim;
    plot(maxEqui,ylim_current(1),'o','MarkerSize', 6,'MarkerFaceColor', ...
    color(4,:),'MarkerEdgeColor',color(4,:));
    
    % Texte sous la marque
    text(maxEqui+0.02, ylim_current(1)+0.02*diff(ylim_current),'$\bar{x}_1$', ...
    'Interpreter','latex','HorizontalAlignment','left','FontSize',16, ...
    'Color',color(4,:));
    hold off
    
    % add labels and formatting
    xlabel('Phenotype ($x$)','Interpreter','latex','FontSize',18);
    ylabel('Tot. mosqui. equi.','Interpreter','latex','FontSize',18);
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);
    
    GrNum= GrNum+1;
    Tot_AA=AA0+AA1;
    
    subplot('Position', [0.68, 0.25, 0.25, 0.68]);  % Troisième subplot
    surf(time,IRLevel,Tot_AA','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:600:300*Ny);
    xticklabels(0:2:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',18,'rotation',-40);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',18,'rotation',40,'HorizontalAlignment','center');
    zlabel('Total mosqui.','Interpreter','latex','fontsize',16);
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);
    
    %save the figure
    % set(gcf,'Units','centimeters','PaperUnits','centimeters','PaperPositionMode','manual');
    set(gca,'LooseInset',get(gca,'TightInset'));

    fileName1 = 'fig2_JMB.pdf';
    exportgraphics(gcf,fullfile(fpath,fileName1),'ContentType','image', ...
        'BackgroundColor','white','Resolution',600);
    
end 
%% SECTION 3: Fitness functions in a mixed insecticide landscape=f(Scenario,RatioEpsi,p0,p1,p2)
if(0)
    %setting parameters values      
    Scenario=3;                   RatioEpsi=2;                  p0=0.5;
    p1=0.2*ones(1,nb_cycles);     p2=0.4*ones(1,nb_cycles);

    %run FonctionMainRun
    [Theta0,Theta1,Theta2,AA0T,AA1T,AA2T,AA0,AA1,AA2,E0,E1,E2,A0,A1,A2,...
        time,T,RelatGain,Temergence,PopAdapt2] = FunctionMainRun(Ny,DurPeriod, ...
        VarM0,EmergenceThreshold,p2,p1,p0,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
        gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
        RatioEpsi);

    %calculate combined fitness function for new landscape with I1 and I2
    Theta=(1-p1(1)-p2(1))*Theta0+ p1(1)*Theta1 + p2(1)*Theta2;
    [~,maxId]= max(Theta1);    max1= IRLevel(maxId);    IdMaxGamma1=maxId;
    [~,maxId]= max(Theta2);    max2= IRLevel(maxId);    IdMaxGamma2=maxId;
    [maxTheta,maxId]= max(Theta); maxGlobal= IRLevel(maxId);
    
    %create graphic
    figure
    set(gcf, 'Position', [100, 100, 240, 240]);
    %Fitness functions
    plot(IRLevel,Theta0,'LineWidth',LW,'Color',color(1,:),'LineStyle',':');
    xlim([0 1]);
    hold on 
    plot(IRLevel,Theta1,'LineWidth',LW,'Color',color(2,:),'LineStyle','-');
    plot(IRLevel,Theta2,'LineWidth',LW,'Color',color(3,:),'LineStyle','-.');
    plot(IRLevel,Theta,'LineWidth',LW,'Color',color(4,:),'LineStyle','--');
    plot([max1 max1],[0 max(Theta1)],'Color',color(2,:),'LineStyle','--');
    plot([max2 max2],[0 max(Theta2)],'Color',color(3,:),'LineStyle','--');
    plot([maxGlobal maxGlobal],[0 max(Theta)],'Color',color(4,:),'LineStyle','--');
    % Markers on x-axis for maxTheta1, maxTheta2 and maxTheta
    ylim_current = ylim;
    plot(max1, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(2,:),'MarkerEdgeColor',color(2,:));
    plot(max2, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(3,:),'MarkerEdgeColor',color(3,:));
    plot(maxGlobal, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(4,:),'MarkerEdgeColor',color(4,:));
    % Text under markers
    text(max1+0.02, ylim_current(1)+0.02*diff(ylim_current), '$x_1$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'left', 'FontSize', 13, 'Color', color(2,:));
    text(max2+0.02, ylim_current(1)+0.02*diff(ylim_current), '$x_2$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'left', 'FontSize', 13, 'Color', color(3,:));
    text(maxGlobal-0.02, ylim_current(1)+0.02*diff(ylim_current), '$\bar{x}_2$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'right', 'FontSize', 13, 'Color', color(4,:));
    
    % Fixed: Calculate values for horizontal lines using correct variable names
    Gamma1_at_maxGamma1 = Theta1(IdMaxGamma1);
    Gamma2_at_maxGamma1 = Theta2(IdMaxGamma1);
    
    % First horizontal line: from y-axis to Gamma2(x1)
    horizontal_line1_y = Gamma2_at_maxGamma1;
    plot([0 max1], [horizontal_line1_y horizontal_line1_y], ...
         'Color','k','LineStyle',':');
    
    % Second horizontal line: from y-axis to Gamma1(x1)
    horizontal_line2_y = Gamma1_at_maxGamma1;
    plot([0 max1], [horizontal_line2_y horizontal_line2_y], ...
         'Color','k','LineStyle',':');
    
    % Add legends for horizontal lines
    text_x_pos = 0.02;
    
    % Text for first horizontal line
    text(text_x_pos, horizontal_line1_y+0.02, '$(1-{\rm r}_{\rm eff})\Gamma_1(x_1)$', ...
         'Interpreter','latex','LineWidth',1,'FontSize',11,'Color','k', ...
         'VerticalAlignment', 'top');
    
    % Text for second horizontal line
    text(text_x_pos, horizontal_line2_y+0.02, '$\Gamma_1(x_1)$', ...
         'Interpreter','latex','LineWidth',1,'FontSize',12,'Color','k', ...
         'VerticalAlignment', 'top');
    hold off
    
    % add labels and formatting
    xlabel('Phenotype ($x$)','Interpreter','latex','FontSize',15);
    ylabel('Fitness','Interpreter','latex','FontSize',15);
    legend('$\Gamma_0$','$\Gamma_1$','$\Gamma_2$','$\Theta_{[p_1,p_2]}$','Interpreter','latex', ...
        'location','east','Orientation','vertical','fontsize',13)
    legend boxoff
    % title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], 'interpreter','latex','FontSize',18);    
    
    %save the figure
    fileName2='Fitness functions_JMB.pdf';
    set(gcf,'Units','centimeters','PaperUnits','centimeters','PaperPositionMode','manual');
    set(gca,'LooseInset',get(gca,'TightInset'));
    exportgraphics(gcf,fullfile(fpath,fileName2),'ContentType','image', ...
        'BackgroundColor','white','Resolution',900);
end

%% SECTION 4: figure typical dynamics=f(Scenario,p0,p1,p2)
if(1)

    % setting parameters values
    % Scenario=2; 
    % p0=0.5;       
    % p1=0.2*ones(1,ceil(Ny/DurPeriod));   p2=0.5*ones(1,ceil(Ny/DurPeriod));

    % Scenario=3; 
    % p0=0.6;  
    % p1= [0.24, 0.20, 0.16, 0.10];    %example of deployment strategy for Scenario=3 and p0=0.6
    % p2= [0.48, 0.60, 0.72, 0.80];    %example of deployment strategy for Scenario=3 and p0=0.6
    
    Scenario=3; 
    p0=0.6; 
    p1= [0.0023, 0.0025, 0.0027, 0.0031];  %optimal p1 for Scenario=3 and p0=0.6
    p2= [0.4898, 0.9931, 0.9943, 0.9938];  %optimal p2 for Scenario=3 and p0=0.6

    % Scenario=4; 
    % p0=0.6; 
    % p1= [0.0017, 0.0027, 0.0020, 0.0022];  %optimal p1 for Scenario=4 and p0=0.6 long format
    % p2= [0.4651, 0.9838, 0.9958, 0.9957];  %optimal p2 for Scenario=4 and p0=0.6 long format
   
    % Run mosquito population model
    [Theta0,Theta1,Theta2,AA0T,AA1T,AA2T,...
    AA0,AA1,AA2,E0,E1,E2,A0,A1,A2,time,T,RelatGain,Temergence,PopAdat2] = FunctionMainRun(Ny,DurPeriod, ...
    VarM0,EmergenceThreshold,p2,p1,p0,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
    gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
    RatioEpsi);

    % Calculate combined fitness function and find optimal phenotype
    Theta=(1-p1(1)-p2(1))*Theta0+ p1(1)*Theta1+p2(1)*Theta2;
        
    [maxTheta1,maxId]= max(Theta1); max1= IRLevel(maxId); IdMaxTheta1=maxId;
    [maxTheta2,maxId]= max(Theta2); max2= IRLevel(maxId); IdMaxTheta2=maxId;
    [maxTheta,maxId]= max(Theta); maxGlobal= IRLevel(maxId); IdMaxGlobal=maxId;

    %figure 3
    figure;
    set(gcf,'position',[100,100,1200,900]);
    axes('fontsize',15)

    %Fitness functions
    GrNum= 1;
    subplot(2,3,1.5);
    plot(IRLevel,Theta0,'LineWidth',LW,'Color',color(1,:),'LineStyle',':');
    xlim([0 1.1]);
    hold on 
    plot(IRLevel,Theta1,'LineWidth',LW,'Color',color(2,:),'LineStyle','-');
    plot(IRLevel,Theta2,'LineWidth',LW,'Color',color(3,:),'LineStyle','-.');
    plot(IRLevel,Theta,'LineWidth',LW,'Color',color(4,:),'LineStyle','--');
    plot([max1 max1],[0 max(Theta1)],'Color',color(2,:),'LineStyle','--');
    plot([max2 max2],[0 max(Theta2)],'Color',color(3,:),'LineStyle','--');
    plot([maxGlobal maxGlobal],[0 max(Theta)],'Color',color(4,:),'LineStyle','--');

    % Marques sur l'axe x pour maxTheta1, maxTheta2 et maxTheta
    ylim_current = ylim;
    plot(max1, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(1,:),'MarkerEdgeColor',color(2,:));
    plot(max2, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(2,:),'MarkerEdgeColor',color(3,:));
    plot(maxGlobal, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(3,:),'MarkerEdgeColor',color(4,:));

    % Texte sous les marques
    text(max1+0.02, ylim_current(1)+0.02*diff(ylim_current), '$x_1$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'left', 'FontSize', 15, 'Color', color(2,:));
    text(max2+0.02, ylim_current(1)+0.02*diff(ylim_current), '$x_2$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'left', 'FontSize', 15, 'Color', color(3,:));
    text(maxGlobal-0.02, ylim_current(1)+0.02*diff(ylim_current), '$\bar{x}_2$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'right', 'FontSize', 15, 'Color', color(4,:));
    
    % Calcul des valeurs pour les lignes horizontales
    Theta1_at_maxTheta1 = Theta1(IdMaxTheta1);
    Theta2_at_maxTheta1 = Theta2(IdMaxTheta1);
    
    % Première ligne horizontale : de l'axe des ordonnées à Theta2(x1)
    horizontal_line1_y = Theta2_at_maxTheta1;
    plot([0 max1], [horizontal_line1_y horizontal_line1_y], ...
         'Color','k','LineWidth',1,'LineStyle',':');
    
    % Deuxième ligne horizontale : de l'axe des ordonnées à Theta1(x1)
    horizontal_line2_y = Theta1_at_maxTheta1;
    plot([0 max1], [horizontal_line2_y horizontal_line2_y], ...
         'Color','k','LineWidth',1,'LineStyle',':');
    
    % Ajout des légendes pour les lignes horizontales
    % Position des textes (à droite, sur l'axe des ordonnées)
    text_x_pos = 0.02;
    
    % Texte pour la première ligne horizontale
    text(text_x_pos, horizontal_line1_y+0.02, '$(1-r_{\rm eff})\Gamma_1(x_1)$', ...
         'Interpreter','latex','FontSize',13,'Color','k', 'LineWidth', 1, ...
         'VerticalAlignment', 'bottom');
    
    % Texte pour la deuxième ligne horizontale
    text(text_x_pos, horizontal_line2_y+0.02, '$\Gamma_1(x_1)$', ...
         'Interpreter','latex','FontSize',13,'Color','k', 'LineWidth', 1, ...
         'VerticalAlignment', 'cap');
    hold off
    
    % add labels and formatting
    xlabel('Phenotype ($x$)','Interpreter','latex','FontSize',18);
    ylabel('Fitness','Interpreter','latex','FontSize',18);
    %title('\textbf{Fitness functions}','interpreter','latex','FontSize',15);
    legend('$\Gamma_0$','$\Gamma_1$','$\Gamma_2$','$\Theta_{[p_1,p_2]}$','Interpreter','latex', ...
        'location','east','Orientation','vertical','fontsize',15)
    legend boxoff
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);

    %Rgain and Temg  
    GrNum= GrNum+1;
    TotPop=AA0T+AA1T+AA2T;
    
    subplot(2,3,2.5);
    yyaxis left
    %colororder(color(2,:))
    plot(time,TotPop,'LineWidth',LW,'linestyle','-');
    xlabel('Time $t$ (year)','Interpreter','latex','fontsize',18);
    xlim([0 T]);
    xticks(0:300:300*Ny);
    xticklabels(0:1:Ny);
    ylabel('Total mosquitoes','Interpreter','latex','fontsize',18);
    
    PopAdat2=PopAdat2/100;
    hold on
    yyaxis right 
    %colororder(color(3,:))
    plot(time,PopAdat2,'LineWidth',LW,'linestyle','-.')
    ylim([0 1]);
    plot([Temergence Temergence], [0 1], 'LineWidth', 1, 'color', 'k', 'linestyle', ':');
    ylabel('$Q_{x_2}$','Interpreter','latex','fontsize',18)
    hold off
    title(['\fontsize{15}{0}\selectfont' '\qquad \qquad \textbf{(' LETTERS(GrNum) ')}' ...
    '\fontsize{15}{0}\selectfont$\:{\rm R}_{\rm gain}=$' num2str(RelatGain,3) ...
    '\ ; \fontsize{15}{0}\selectfont$\:{\rm T}_{\rm emerg}=$' num2str(Temergence/300,3)], ...
    'interpreter','latex', 'FontSize',20);
    legend('Tot. mosqui.','$Q_{x_2}$','Interpreter','latex', ...
        'location','southeast','Orientation','vertical','fontsize',15)
    legend boxoff


    %unexposed AM
    GrNum= GrNum+1;
    subplot(2,3,4);
    surf(time,IRLevel,AA0','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:600:300*Ny);
    xticklabels(0:2:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',16, 'rotation',-50);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15, 'rotation',50);
    zlabel('Unexposed mosqui.','Interpreter','latex','fontsize',15);
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);


    % I1 exposed AM
    GrNum= GrNum+1;
    subplot(2,3,5);
    surf(time,IRLevel,AA1','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:600:300*Ny);
    xticklabels(0:2:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',16, 'rotation',-50);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15, 'rotation',50);
    zlabel('$I_1$ exposed mosqui.','Interpreter','latex','fontsize',15);
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);


    % I2 exposed AM
    GrNum= GrNum+1;
    subplot(2,3,6);
    surf(time,IRLevel,AA2','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:600:300*Ny);
    xticklabels(0:2:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',16, 'rotation',-40);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15, 'rotation',40);
    zlabel('$I_2$ exposed mosqui.','Interpreter','latex','fontsize',15);
    title(['\fontsize{15}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);

    
    fileName3 = 'fig4Opti_JMB.pdf';
    set(gcf,'Units','centimeters','PaperUnits','centimeters','PaperPositionMode','manual');
    set(gca,'LooseInset',get(gca,'TightInset'));
    exportgraphics(gcf,fullfile(fpath,fileName3),'ContentType','image', ...
        'BackgroundColor','white','Resolution',400);
end

%% SECTION 5: optimisation= fction(Scenario,p0)
if(0)
    format long
    Scenario=4; 
    [Theta0,Theta1,Theta2,~,~,~,~,~]=...
        FuncTheta(tau,age,IRLevel,Na,Nx,Scenario,gamma0,gamma1,gamma2,...
        mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan);
    
   %pour p0=p1 fixé, on met le système à l'équilibre en l'absence de p2
   %set initial parameters (p0, NyEqui, p1equi and p2equi)
   p0=0.6;
   NyEqui=5;
   p1equi= [p0,p0,p0,p0]; p2equi= 0*p1equi;

   %find the maximum value of Theta1/Theta2 and its position
   Theta=(1-p1equi(1)-p2equi(1))*Theta0+ p1equi(1)*Theta1+p2equi(1)*Theta2;
   [~,maxId]= max(Theta1); max1= IRLevel(maxId);      IdMaxTheta1=maxId;
   [~,maxId]= max(Theta2); max2= IRLevel(maxId);      IdMaxTheta2=maxId;
   [~,maxId]= max(Theta);  maxGlobal= IRLevel(maxId); IdMaxGlobal=maxId;

    %démarre suffisament proche de cet équilibre pour l'atteindre en Ny=5 ans
    InitialE0= 1*normpdf(IRLevel,maxGlobal,VarM0);  
    InitialE1= 1*normpdf(IRLevel,maxGlobal,VarM0); 
    InitialE2= zeros(1,Nx); 
    InitialA0= zeros(Na,Nx); InitialA1= zeros(Na,Nx); InitialA2= zeros(Na,Nx);

    [~,~,~,~,~,~,~,AA0T,AA1T,AA2T,~,~,...
        ~,E0,E1,E2,A0,A1,A2,~]=MainFunctionIRM2(Ny,DurPeriod,VarM0,p2equi,p1equi,...
        Scenario,tau,age,IRLevel,dx,Nx,da,Na,gamma0,gamma1,gamma2,mu0,mu1,...
        mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,RatioEpsi,InitialE0,...
        InitialE1,InitialE2,InitialA0,InitialA1,InitialA2);

    %on calcul le nbre de moustique total à l'équilibre   
    TotalMosquiEqui= AA0T(end) + AA1T(end) + AA2T(end);

    % on récupére l'équilibre du système en l'absence de I2
    InitialE0= E0(end,:);       InitialA0= A0(:,:,end); 
    InitialE1= E1(end,:);       InitialA1= A1(:,:,end);
    InitialE2= E2(end,:);       InitialA2= A2(:,:,end);
   
    % Fonction objective pour seule variable p1_and_p2
    objFun = @(p1_and_p2) FunctionMainOpti(p1_and_p2,Ny,DurPeriod,VarM0,Scenario, ...
    tau,age,IRLevel,dx,Nx,da,Na,gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,...
    AvrMosquiLifeSpan,RatioEpsi,InitialE0,InitialE1,InitialE2,...
    InitialA0,InitialA1,InitialA2,TotalMosquiEqui);
    
    %Solving the system with Optimal strategy
    Ops=optimset('Display','iter','Algorithm','interior-point','MaxFunEvals',10^5,'MaxIter',2000,'TolX',0.0001,'TolFun',0.0001);
    %fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

    SizePeriod=ceil(Ny/DurPeriod);
    AOpti= [eye(SizePeriod),eye(SizePeriod)];
    bOpti= ones(SizePeriod,1);

    InitP1=rand(1,SizePeriod);
    InitP2=1-InitP1;
    p1_and_p2_Init=[InitP1,InitP2];

    lb=zeros(1,2*SizePeriod);
    ub=ones(1,2*SizePeriod);

    tic
    Prop=fmincon(objFun,p1_and_p2_Init,AOpti,bOpti,[],[],lb,ub,[],Ops);
    toc

    p1=Prop(1:SizePeriod)                 
    %p1 =[0.0023, 0.0025, 0.0027, 0.0031]  optimal p1 for Scenario=3 and p0=0.6 
    %p1 =[0.002299023682400, 0.002451932849448, 0.002700805550697, 0.003090283559530]  optimal p1 for Scenario=3 and p0=0.6 long format
    %p1 =[0.001720708438592, 0.002719652684623, 0.001951743278483, 0.002157233434967]  optimal p1 for Scenario=4 and p0=0.6 long format

    p2=Prop(SizePeriod+1:2*SizePeriod)    
    %p2 =[0.4898, 0.9931, 0.9943, 0.9938]  optimal p2 for Scenario=3 and p0=0.6
    %p2 =[0.489811182413074, 0.993098031865287, 0.994301057591897, 0.993810319875747]  optimal p2 for Scenario=3 and p0=0.6 long format
    %p2 =[0.465055750419891, 0.983798221781905, 0.995799984591818, 0.995657969753122]  optimal p2 for Scenario=4 and p0=0.6 long format

    % on tourne le model avec l'optimal
    [~,~,~,AA0T,AA1T,AA2T,...
    AA0,AA1,AA2,E0,E1,E2,A0,A1,A2,time,T,RelatGain,Temergence,PopAdat2] = ...
        FunctionMainRun(Ny,DurPeriod,VarM0,EmergenceThreshold,p2,p1,p0,...
        Scenario,tau,age,IRLevel,dx,Nx,da,Na,gamma0,gamma1,gamma2,mu0,mu1,...
        mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,RatioEpsi);

    % Calculate combined fitness function and find optimal phenotype
    Theta=(1-p1(1)-p2(1))*Theta0+ p1(1)*Theta1+p2(1)*Theta2;
        
    [maxTheta1,maxId]= max(Theta1); max1= IRLevel(maxId);      IdMaxTheta1=maxId;
    [maxTheta2,maxId]= max(Theta2); max2= IRLevel(maxId);      IdMaxTheta2=maxId;
    [maxTheta,maxId]= max(Theta);   maxGlobal= IRLevel(maxId); IdMaxGlobal=maxId;

    %figure 4
    figure;
    set(gcf,'position',[100,100,1200,900]);
    axes('fontsize',15)

    %Fitness functions
    GrNum= 1;
    subplot(2,3,1.5);
    plot(IRLevel,Theta0,'LineWidth',LW,'Color',color(4,:),'LineStyle',':');
    xlim([0 1]);
    hold on 
    plot(IRLevel,Theta1,'LineWidth',LW,'Color',color(1,:),'LineStyle','-');
    plot(IRLevel,Theta2,'LineWidth',LW,'Color',color(2,:),'LineStyle','-.');
    plot(IRLevel,Theta,'LineWidth',LW,'Color',color(3,:),'LineStyle','--');
    plot([max1 max1],[0 max(Theta1)],'Color',color(1,:),'LineStyle','--');
    plot([max2 max2],[0 max(Theta2)],'Color',color(2,:),'LineStyle','--');
    plot([maxGlobal maxGlobal],[0 max(Theta)],'Color',color(3,:),'LineStyle','--');

    % Marques sur l'axe x pour maxTheta1, maxTheta2 et maxTheta
    ylim_current = ylim;
    plot(max1, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(1,:),'MarkerEdgeColor',color(1,:));
    plot(max2, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(2,:),'MarkerEdgeColor',color(2,:));
    plot(maxGlobal, ylim_current(1),'o','MarkerSize',6,'MarkerFaceColor', ...
        color(3,:),'MarkerEdgeColor',color(3,:));

    % Texte sous les marques
    text(max1+0.02, ylim_current(1)+0.02*diff(ylim_current), '$x_1$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'left', 'FontSize', 15, 'Color', color(1,:));
    text(max2+0.02, ylim_current(1)+0.02*diff(ylim_current), '$x_2$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'left', 'FontSize', 15, 'Color', color(2,:));
    text(maxGlobal-0.02, ylim_current(1)+0.02*diff(ylim_current), '$\bar{x}_2$', 'Interpreter', 'latex', ...
         'HorizontalAlignment', 'right', 'FontSize', 15, 'Color', color(3,:));
    
    % Calcul des valeurs pour les lignes horizontales
    Theta1_at_maxTheta1 = Theta1(IdMaxTheta1);
    Theta2_at_maxTheta1 = Theta2(IdMaxTheta1);
    
    % Première ligne horizontale : de l'axe des ordonnées à Theta2(x1)
    horizontal_line1_y = Theta2_at_maxTheta1;
    plot([0 max1], [horizontal_line1_y horizontal_line1_y], ...
         'Color','k','LineStyle',':');
    
    % Deuxième ligne horizontale : de l'axe des ordonnées à Theta1(x1)
    horizontal_line2_y = Theta1_at_maxTheta1;
    plot([0 max1], [horizontal_line2_y horizontal_line2_y], ...
         'Color','k','LineStyle',':');
    
    % Ajout des légendes pour les lignes horizontales
    % Position des textes (à droite, sur l'axe des ordonnées)
    text_x_pos = 0.02;
    
    % Texte pour la première ligne horizontale
    text(text_x_pos, horizontal_line1_y+0.02, '$(1-r_{\rm eff})\Gamma_1(x_1)$', ...
         'Interpreter','latex','FontSize',13,'Color','k', 'LineWidth', 1, ...
         'VerticalAlignment', 'cap');
    
    % Texte pour la deuxième ligne horizontale
    text(text_x_pos, horizontal_line2_y+0.02, '$\Gamma_1(x_1)$', ...
         'Interpreter','latex','FontSize',13,'Color','k', 'LineWidth', 1, ...
         'VerticalAlignment', 'cap');
    hold off
    
    % add labels and formatting
    xlabel('Phenotype ($x$)','Interpreter','latex','FontSize',18);
    ylabel('Fitness','Interpreter','latex','FontSize',18);
    %title('\textbf{Fitness functions}','interpreter','latex','FontSize',15);
    legend('$\Gamma_0$','$\Gamma_1$','$\Gamma_2$','$\Theta_{[p_1,p_2]}$','Interpreter','latex', ...
        'location','east','Orientation','vertical','fontsize',15)
    legend boxoff
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);

    %Rgain and Temg  
    GrNum= GrNum+1;
    TotPop=AA0T+AA1T+AA2T;
    
    subplot(2,3,2.5);
    yyaxis left
    %colororder(color(2,:))
    plot(time,TotPop,'LineWidth',LW,'linestyle','-');
    xlabel('Time $t$ (year)','Interpreter','latex','fontsize',18);
    xlim([0 T]);
    xticks(0:300:300*Ny);
    xticklabels(0:1:Ny);
    ylabel('Total mosquitoes','Interpreter','latex','fontsize',18);
    
    PopAdat2=PopAdat2/100;
    hold on
    yyaxis right 
    %color order(color(3,:))
    plot(time,PopAdat2,'LineWidth',LW,'linestyle','-.')
    ylim([0 1]);
    plot([Temergence Temergence], [0 1], 'LineWidth', 1, 'color', 'k', 'linestyle', ':');
    ylabel('$Q_{x_2}$','Interpreter','latex','fontsize',18)
    hold off
    title(['\fontsize{10}{0}\selectfont' '\qquad \qquad \textbf{(' LETTERS(GrNum) ')}' ...
    '\fontsize{10}{0}\selectfont$\:{\rm R}_{\rm gain}=$' num2str(RelatGain,3) ...
    '\ ; \fontsize{10}{0}\selectfont$\:{\rm T}_{\rm emerg}=$' num2str(Temergence/300,3)], ...
    'interpreter','latex', 'FontSize',20);
    legend('Tot. mosqui.','$Q_{x_2}$','Interpreter','latex', ...
        'location','southeast','Orientation','vertical','fontsize',15)
    legend boxoff


    %unexposed AM
    GrNum= GrNum+1;
    subplot(2,3,4);
    surf(time,IRLevel,AA0','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:300:300*Ny);
    xticklabels(0:1:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',16, 'rotation',-40);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15, 'rotation',40);
    zlabel('Unexposed mosqui.','Interpreter','latex','fontsize',15);
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);


    % I1 exposed AM
    GrNum= GrNum+1;
    subplot(2,3,5);
    surf(time,IRLevel,AA1','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:300:300*Ny);
    xticklabels(0:1:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',16, 'rotation',-40);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15, 'rotation',40);
    zlabel('$I_1$ exposed mosqui.','Interpreter','latex','fontsize',15);
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);


    % I2 exposed AM
    GrNum= GrNum+1;
    subplot(2,3,6);
    surf(time,IRLevel,AA2','FaceColor','interp', 'EdgeColor','none', ...
    'FaceLighting','gouraud')
    axis tight
    view(-136,60)
    camlight left
    colormap(gca, 'winter') % Colormap à dominante bleue
    xticks(0:300:300*Ny);
    xticklabels(0:1:Ny);
    
    xlabel('time ($t$)','Interpreter','latex','fontsize',16, 'rotation',-40);
    ylabel('Phenotype ($x$)','Interpreter','latex','fontsize',15, 'rotation',40);
    zlabel('$I_2$ exposed mosqui.','Interpreter','latex','fontsize',15);
    title(['\fontsize{10}{0}\selectfont' '\textbf{(' LETTERS(GrNum) ')}'], ...
    'interpreter','latex','FontSize',20);
 
    %save figure    
    fileName4 = 'fig3OptiB.pdf';
    set(gcf,'Units','centimeters','PaperUnits','centimeters','PaperPositionMode','manual');
    set(gca,'LooseInset',get(gca,'TightInset'));
    exportgraphics(gcf,fullfile(fpath,fileName4),'ContentType','image', ...
        'BackgroundColor','white','Resolution',600);
end

%% SECTION 6: Mos vs Rot vs MosOpti
if(0)
    format long
    tic
    
    % Parameter vectors - modify these as needed
    VectStrategy = 3;                    % lStra=length(VectStrategy); %Mos=1, Rot=2, MosOpti=3
    VectVar0 = [0.01, 0.01/5, 0.01/10];  % lVar0=length(VectVar0);
    VectScenario = 2:4;                  % lScena=length(VectScenario);
    VectDurCycle = 3;                    % lDurCycle=length(VectDurCycle);          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VectP0 = 0.05:0.1:0.45;                       % lP0=length(VectP0);
    
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
    tabName = sprintf('TabForMosRotMosOpti_Strategy3_DurPeriod%d0545.xlsx', DurPeriod); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

%% des functions intermédiaires

function[RelatMosqui]=FunctionMainOpti(p1_and_p2,Ny,DurPeriod,VarM0,Scenario, ...
    tau,age,IRLevel,dx,Nx,da,Na,gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,...
    AvrMosquiLifeSpan,RatioEpsi,InitialE0,InitialE1,InitialE2,...
    InitialA0,InitialA1,InitialA2,TotalMosquiEqui)


nn0=length(p1_and_p2)/2;
p1=p1_and_p2(1:nn0); 
p2=p1_and_p2(nn0+1:2*nn0);


[T,~,~,~,~,~,~,~,~,~,~,~,...
    ~,~,~,~,~,~,~,TotalMosqui]=MainFunctionIRM2(Ny,DurPeriod, ...
    VarM0,p2,p1,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
    gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
    RatioEpsi,InitialE0,InitialE1,InitialE2,InitialA0,InitialA1,InitialA2);

RelatMosqui=TotalMosqui/(T*TotalMosquiEqui);
        
end

function[p1,p2,RelatGain,Temergence]=FunctionMos_Rot_MosOpti(strategy,VarM0,Scenario,p0,...
    Ny,DurPeriod,tau,age,IRLevel,dx,Nx,da,Na,gamma0,gamma1,gamma2,mu0,mu1,mu2,...
    MinMosquiAgeLayEggs,AvrMosquiLifeSpan,RatioEpsi,EmergenceThreshold)

 
    if strategy==1%Mosaic
        %p1=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]; % DurCycle=1
        %p1=[0.5,0.5,0.5,0.5,0.5,0.5]; % DurCycle=2
        %p1=[0.5,0.5,0.5,0.5]; % DurCycle=3
        p1=[0.5,0.5,0.5];      % DurCycle=4
        p2=p1;
    elseif strategy==2%Rot
        %p1=[0,1,0,1,0,1,0,1,0,1,0,1];   %Durcycle=1
        %p1=[0,1,0,1,0,1];   %Durcycle=2
        %p1=[0,1,0,1];   %Durcycle=3
        p1=[0,1,0];   %Durcycle=4
        p2=1-p1;
    elseif strategy==3%MosOpti
        [Theta0,Theta1,Theta2,~,~,~,~,~]=...
            FuncTheta(tau,age,IRLevel,Na,Nx,Scenario,gamma0,gamma1,gamma2,...
            mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan);
    
      
       %pour p0=p1 fixé, on met le système à l'équilibre en l'absence de p2 
       NyEqui=5;
       p1equi =[p0,p0,p0,p0]; p2equi =0*p1equi;
       % find the maximum value of Theta1/Theta2 and its position
       Theta=(1-p1equi(1)-p2equi(1))*Theta0+ p1equi(1)*Theta1+p2equi(1)*Theta2;
       [maxTheta1,maxId]= max(Theta1); max1= IRLevel(maxId); 
       [maxTheta2,maxId]= max(Theta2); max2= IRLevel(maxId); IdMaxTheta2=maxId;
       [maxTheta,maxId]= max(Theta); maxGlobal= IRLevel(maxId); IdMaxGlobal=maxId;

        % démarre suffisament proche de cet équilibre pour l'atteindre en
        % Ny=5 ans
        InitialE0=1*normpdf(IRLevel,maxGlobal,VarM0);  
        InitialE1= 1*normpdf(IRLevel,maxGlobal,VarM0); 
        InitialE2= zeros(1,Nx); 
        InitialA0= zeros(Na,Nx); InitialA1= zeros(Na,Nx); InitialA2= zeros(Na,Nx);

        [~,time,~,~,~,~,~,AA0T,AA1T,AA2T,~,~,...
            ~,E0,E1,E2,A0,A1,A2,~]=MainFunctionIRM2(NyEqui,DurPeriod, ...
            VarM0,p2equi,p1equi,Scenario,tau,age,IRLevel,dx,Nx,da,Na,gamma0,...
            gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
            RatioEpsi,InitialE0,InitialE1,InitialE2,InitialA0,InitialA1,InitialA2);

        %on calcul le nbre de moustique total à l'équilibre   
        TotalMosquiEqui= AA0T(length(time)) + AA1T(length(time)) + ...
            AA2T(length(time));

        % on récupére l'équilibre du système en l'absence de I2
        InitialE0=E0(length(time),:); InitialA0= A0(:,:,length(time)); 
        InitialE1= E1(length(time),:); InitialA1= A1(:,:,length(time));
        InitialE2= E2(length(time),:); InitialA2=A2(:,:,length(time));
        
        % Fonction objective pour seule variable p1_and_p2
        objFun = @(p1_and_p2) FunctionMainOpti(p1_and_p2,Ny,DurPeriod,VarM0,Scenario, ...
        tau,age,IRLevel,dx,Nx,da,Na,gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,...
        AvrMosquiLifeSpan,RatioEpsi,InitialE0,InitialE1,InitialE2,...
        InitialA0,InitialA1,InitialA2,TotalMosquiEqui);
        
        %Solving the system with Optimal strategy
        Ops=optimset('Display','iter','Algorithm','interior-point','MaxFunEvals',10^5,'MaxIter',2000,'TolX',0.0001,'TolFun',0.0001);
        %fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
        SizePeriod=ceil(Ny/DurPeriod);
        AOpti=zeros(SizePeriod,2*SizePeriod);
        for jj=1:SizePeriod
            AOpti(jj,jj)=1;
            AOpti(jj,jj+SizePeriod)=1;
        end
        bOpti=ones(SizePeriod,1);
        InitP1=0.5*rand(1,SizePeriod);
        InitP2=0.5*rand(1,SizePeriod);
        p1_and_p2_Init=[InitP1,InitP2];
        lb=zeros(1,2*SizePeriod);
        ub=ones(1,2*SizePeriod);
        tic
        Prop=fmincon(objFun,p1_and_p2_Init,AOpti,bOpti,[],[],lb,ub,[],Ops);
        toc
        p1=Prop(1:SizePeriod);
        p2=Prop(SizePeriod+1:2*SizePeriod);
    end
    
    [~,~,~,~,~,~,...
        ~,~,~,~,~,~,~,~,~,~,~,RelatGain,Temergence] = FunctionMainRun(Ny,DurPeriod, ...
        VarM0,EmergenceThreshold,p2,p1,p0,Scenario,tau,age,IRLevel,dx,Nx,da,Na,...
        gamma0,gamma1,gamma2,mu0,mu1,mu2,MinMosquiAgeLayEggs,AvrMosquiLifeSpan,...
        RatioEpsi);

end

