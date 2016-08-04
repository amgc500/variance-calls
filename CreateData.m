%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file contains code to accompany the paper 
% "Optimal robust bounds for variance options" 
% by Alexander M G Cox and Jiajie Wang
%
% Copyright (c) 2013 Alexander Cox and Jiajie Wang
%
% Please send comments, corrections etc. to
% a.m.g.cox@bath.ac.uk
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to generate output for our paper on Variance Calls.

try
    clear all;
    
    %% Set relevant parameters
    
    % Parameters for a Heston model
    T = 1;
    rho = -0.65;
    V0 = 0.04;
    theta = 0.035;
    kappa = 1.2;
    r =  0.0;
    xi = 0.5;
    s0 = 2;
    
    % Parameters for rho-comparison
    T_rho = 4;
    RhoVarCallK = 0.035; % Strike for variance options as a function of Rho
    NumRho = 11; % 21 - Number of different values of rho to compute to identify effect of rho on price bounds in Heston

    % Parameters for rho-time comparison
    T_rho_max = 10; % Largest value of T to compute
    T_rho_min = 1; % Smallest value of T to compute
    RhoTVarCallK = 0.035; % Strike for variance options will scale with T (so multiply by T)
    NumTRho = 11; % 21 - Number of different values of T to compute for large-T behaviour
    
    % Parameters for a 'model-risk' variant:
    T2 = T;
    rho2 = 0.5;
    V02 = V0;
    theta2 = 0.07;
    kappa2 = 2.4;
    r2 =  r;
    xi2 = xi;
    s02 = s0;
    
    % Parameters for plots
    u_plot = 8;
    xmin = 1.25; % Minimum asset value for plots
    xmax = 2.75; % Maximum asset value for plots

    % Parameters for sizes - Commented numbers are large(!) defaults
    space_size_large = 8000; % 8000 - Number of grid points for step one barriers
    space_size = 4000; % 4000 - Number of grid points for barriers
    time_steps_Root_large = 16000; % 16000 - Number of time steps for Root
    time_steps_Rost_large = 16000; % 16000 - Number of time steps for Rost
    time_steps_Root = 2000; % 2000 - Number of time steps for Root
    time_steps_Rost = 2000; % 2000 - Number of time steps for Rost
    time_step_changes_Rost = 10; % 10 - Number of different time steps to use in Rost
    
    NumBarriers = 6; % Number of barriers to compute for different maturities    
    NumVarStrikes = 17; % Number of strikes for computing variance calls
    MaxVarStrike = 0.25; % Largest strike for a variance call. Will be equally distributed
    
    
    % Paramters for simulations
    NoSims = 50000; % 50000- Number of simulations to run to check barriers embed correctly
    NoHedgeSims = 2000; % 2000 - Number of simulations to run for hedge distributions
    hist_boxes = 20; % Number of boxes for histogram plots
    pct_low_hist = 0.01; % Percentage of the tail to censor for histogram
    time_steps_sims = 20000; % 20000 - Number of (equally spaced) time steps to use in the simulations.
    
    % Color vector
    X = jet(16);
    colRoot = 1; % Root Colour
    colRost = 16; % Rost Colour
    colOther = 8; % Heston Colour
    
    disp('Warning: This code may take a long time to complete! Make the default parameters smaller for faster compile times.');
    
    %% Set up directory to hold data
        
    current_date = sprintf('%04u%02u%02u%02u%02u%02u',fix(clock));
    
    dirpath = strcat('data/',current_date,'/');
    
    if ~mkdir(dirpath)
        error('Barriers:Output:NoDirCreated','Unable to Create Directory')
    end
    
    diary(strcat(dirpath,'output.txt'));
    
    diary on;
        
    sourcefiles = {'BarrierClass.m','HestonModel.m',...
        'CreateData.m','PriceModel.m','RootBarrierClass.m',...
        'RostBarrierClass.m','CPA.m'};
    
    for i = 1:length(sourcefiles)
        if ~copyfile(sourcefiles{i},strcat(dirpath,sourcefiles{i}))
            error('Barriers:Output:NoSave','Unable to Save Files')
        end
    end

    %% Use Heston to generate Call prices
    
    HP = HestonModel(s0,V0,xi,r,T,kappa,theta,rho);
    
    display('Computing barriers')
    Rost = RostBarrierClass(HP,space_size_large,time_steps_Rost_large,time_step_changes_Rost,1,u_plot);
    Root = RootBarrierClass(HP,space_size_large,time_steps_Root_large,1,u_plot);
    
    hold off
    Rost.SetPlotType(1);
    Rost.PlotBarrier(colRost);
    Rost.SetPlotType(0);
    
    %title('Rost Barrier for Heston Call Price')
    ylabel('Integrated Variance')
    xlabel('Asset price')
    xlim([xmin xmax]);
    
    saveas(gcf,strcat(dirpath,'Heston_Rost_Barrier.eps'),'psc2')
    
    hold off
    
    Rost.u_plot()
    %title('Rost Viscosity Solution for Heston Call Price')
    ylabel('Value of u')
    xlabel('Asset price')
    xlim([xmin xmax]);
    saveas(gcf,strcat(dirpath,'Heston_Rost_u.eps'),'psc2')
    
    
    Root.SetPlotType(1);
    Root.PlotBarrier(colRoot);
    Root.SetPlotType(0);
    
    %title('Root Barrier for Heston Call Price')
    ylabel('Integrated Variance')
    xlabel('Asset price')
    xlim([xmin xmax]);
    y = ylim;
    ylim([0,y(2)]);
    
    saveas(gcf,strcat(dirpath,'Heston_Root_Barrier.eps'),'psc2')

    Root.u_plot()
    %title('Root Viscosity Solution for Heston Call Price')
    ylabel('Value of u')
    xlabel('Asset price')
    xlim([xmin xmax]);

    saveas(gcf,strcat(dirpath,'Heston_Root_u.eps'),'psc2')
    
    display('Running simulations')
    Rost.simulate(NoSims);
    Root.simulate(NoSims);
    
    Rost.PlotEmpCall(colRost);
    hold on
    Root.PlotEmpCall(colRoot);
    
    %title('Heston Call Price and Empirical Call Prices')
    xlabel('Strike Price')
    ylabel('Call Price')
    legend('Rost','Actual','Root')
    xlim([xmin xmax]);

    saveas(gcf,strcat(dirpath,'Heston_Empirical_Calls.eps'),'psc2')

    hold off
    
    Rost.CompareImpliedVol(colRost);
    hold on
    Root.CompareImpliedVol(colRoot);
    
    %title('Heston Implied Volatility and Empirical Implied Volatilities')
    xlabel('Strike Price')
    ylabel('Implied Volatility')
    legend('Rost','Actual','Root')
    xlim([xmin xmax]);

    saveas(gcf,strcat(dirpath,'Heston_Empirical_ImpVol.eps'),'psc2')
    
    hold off

    T_vec = linspace(T/NumBarriers,T,NumBarriers);
    
    HPCell = cell(NumBarriers,1);
    RostCell = cell(NumBarriers,1);
    RootCell = cell(NumBarriers,1);
    
    HPCell{NumBarriers} = HP;
    RostCell{NumBarriers} = Rost;
    RootCell{NumBarriers} = Root;
    
    display('Computing multiple barriers')
    
    for i = 1:(NumBarriers-1)
        HPCell{i} = HestonModel(s0,V0,xi,r,T_vec(i),kappa,theta,rho);
        RostCell{i} = RostBarrierClass(HPCell{i},space_size,time_steps_Rost,time_step_changes_Rost,1);
        RootCell{i} = RootBarrierClass(HPCell{i},space_size,time_steps_Root,1);
    end
    
    RostCell{1}.PlotBarrier(colRost);
    hold on
    for i = 2:NumBarriers
        RostCell{i}.PlotBarrier(colRost);
    end
    
    %title('Reversed Barriers for Heston Calls as Maturity Increases')
    xlabel('Integrated Variance')
    ylabel('Asset Price')
    ylim([xmin xmax]);
    
    saveas(gcf,strcat(dirpath,'Heston_Rost_Barriers_in_Maturity.eps'),'psc2')
    
    hold off
    
    clf
    
    RootCell{1}.PlotBarrier(colRoot);
    hold on
    for i = 2:NumBarriers
        RootCell{i}.PlotBarrier(colRoot);
    end
    
    %title('Barriers for Heston Calls as Maturity Increases')
    xlabel('Integrated Variance')
    ylabel('Asset Price')
    ylim([xmin xmax]);
    
    saveas(gcf,strcat(dirpath,'Heston_Root_Barriers_in_Maturity.eps'),'psc2')
    
    hold off

    VarK = linspace(0,MaxVarStrike,NumVarStrikes);
    RostVarCall = zeros(NumVarStrikes,NumBarriers);
    HestonVarCall = zeros(NumVarStrikes,NumBarriers);
    
    for i = 1:NumBarriers
        for j = 1:NumVarStrikes
            f = str2func(strcat('@(t) 1*(t>',num2str(VarK(j),10),')'));
            F = str2func(strcat('@(t) max(t-',num2str(VarK(j),10),',0)'));
            RostVarCall(j,i) = RostCell{i}.price(f);
            HestonVarCall(j,i) = HPCell{i}.VaroptionPrice(F);
        end
    end

    plot(VarK,HestonVarCall(:,1),'Color',[0 0 0.5],'LineStyle','-.');

    hold on;
    
    plot(VarK,HestonVarCall(:,NumBarriers),'Color',[0 0.5 0],'LineStyle','-.');
    
    for i = 1:NumBarriers
        plot(VarK,RostVarCall(:,i),'Color',X(colRost,:));
    end
    
    
    %title('Heston Variance Call Upper Bounds as Maturity Increases')
    xlabel('Strike Price')
    ylabel('Option Price')
    legend(strcat('Actual, T = ',num2str(T_vec(1),4)),strcat('Actual, T = ',num2str(T_vec(NumBarriers),4)));
        
    saveas(gcf,strcat(dirpath,'Heston_Var_Call_Upper_Bound_in_Maturity.eps'),'psc2')

    hold off
    

    %% Simulate Hedge, upper and lower bounds

    
    t_vec = (0:1:time_steps_sims)*(1/time_steps_sims);

    [p,iv] = HP.simulate(t_vec);

    VarK = T*(theta+V0)/4;
    
    % f = str2func(strcat('@(t) 1*(t>',num2str(VarK,10),')'));
    % F = str2func(strcat('@(t) max(t-',num2str(VarK,10),',0)'));

    f = str2func(strcat('@(t) min(2*t,2*',num2str(VarK,10),')'));
    F = str2func(strcat('@(t) ((t<',num2str(VarK,10),').*t.^2+(t>=',num2str(VarK,10),').*(2*t*',num2str(VarK,10),'-',num2str(VarK,10),'*',num2str(VarK,10),'))'));
    
    [Rost_a,Rost_c] = Rost.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
    [Root_a,Root_c] = Root.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
    
    plot(t_vec,F(iv),'Color',X(colOther,:))
    
    hold on
    
    plot(t_vec,Rost_a+Rost_c,'Color',X(colRost,:))
    plot(t_vec,Root_a+Root_c,'Color',X(colRoot,:))
    
    hold off

    
    %title('Upper and Lower hedges of Variance Option payoff')
    xlabel('Time')
    ylabel('Value')
    legend('Intrinsic value','Superhedge (Rost)','Subhedge (Root)')
    legend('Location','NorthWest')
   
    saveas(gcf,strcat(dirpath,'Heston_Hedge_True_Model.eps'),'psc2')
    
    % Try under the wrong model.
    
    HP2 = HestonModel(s02,V02,xi2,r2,T2,kappa2,theta2,rho2);
        
    [p,iv] = HP2.simulate(t_vec);
    
    f = str2func(strcat('@(t) min(2*t,2*',num2str(VarK,10),')'));
    F = str2func(strcat('@(t) ((t<',num2str(VarK,10),').*t.^2+(t>=',num2str(VarK,10),').*(2*t*',num2str(VarK,10),'-',num2str(VarK,10),'*',num2str(VarK,10),'))'));

    
    [Rost_a,Rost_c] = Rost.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
    [Root_a,Root_c] = Root.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
    
    plot(t_vec,F(iv),'Color',X(colOther,:))
    
    hold on
    
    plot(t_vec,Rost_a+Rost_c,'Color',X(colRost,:))
    plot(t_vec,Root_a+Root_c,'Color',X(colRoot,:))
    
    hold off

    
    %title('Upper and Lower hedges of Variance Option payoff in the wrong model')
    xlabel('Time')
    ylabel('Value')
    legend('Intrinsic value','Superhedge (Rost)','Subhedge (Root)')
    legend('Location','NorthWest')

    saveas(gcf,strcat(dirpath,'Heston_Hedge_Wrong_Model.eps'),'psc2')
    
    %% Do repeated simulations, under correct and incorrect model, plot histograms
    
    Pyff = zeros(NoHedgeSims,1);
    Pyff_Wrong = zeros(NoHedgeSims,1);
    Final_Correct_Root = zeros(NoHedgeSims,1);
    Final_Wrong_Root = zeros(NoHedgeSims,1);
    Final_Correct_Rost = zeros(NoHedgeSims,1);
    Final_Wrong_Rost = zeros(NoHedgeSims,1);
    
    for i = 1:NoHedgeSims
        [p,iv] = HP.simulate(t_vec);
        [Rost_a,Rost_c] = Rost.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
        [Root_a,Root_c] = Root.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
        Pyff(i) = F(iv(time_steps_sims));
        Final_Correct_Root(i) = Root_a(time_steps_sims) + Root_c(time_steps_sims);
        Final_Correct_Rost(i) = Rost_a(time_steps_sims) + Rost_c(time_steps_sims);
        
        [p,iv] = HP2.simulate(t_vec);
        [Rost_a,Rost_c] = Rost.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
        [Root_a,Root_c] = Root.PathwiseSimulateHedgeLong(f,p,t_vec,iv);
        Pyff_Wrong(i) = F(iv(time_steps_sims));
        Final_Wrong_Root(i) = Root_a(time_steps_sims) + Root_c(time_steps_sims);
        Final_Wrong_Rost(i) = Rost_a(time_steps_sims) + Rost_c(time_steps_sims);
    end
    
    
    x1 = sort([Pyff-Final_Correct_Root;Pyff_Wrong-Final_Wrong_Root]);
    xout1 = linspace(x1(floor(pct_low_hist*length(x1))),x1(ceil((1-pct_low_hist)*length(x1))),hist_boxes);
    [n2,~] = hist(Pyff-Final_Correct_Root,xout1);
    [n3,~] = hist(Pyff_Wrong-Final_Wrong_Root,xout1);

    bar(xout1',[n2' n3'])
    colormap winter
    

    xlabel('Payoff - Hedge Value')
    ylabel('Frequency')
    title(strcat('Distribution of hedging errors, censored at ',num2str(xout1(1)),' and ',num2str(xout1(length(xout1)))));
    legend('Correct model','Wrong Model')
    
    saveas(gcf,strcat(dirpath,'Heston_Subhedge_Error_Compare.eps'),'psc2')   

    x1 = sort([Final_Correct_Rost-Pyff;Final_Wrong_Rost-Pyff_Wrong]);
    xout1 = linspace(x1(floor(pct_low_hist*length(x1))),x1(ceil((1-pct_low_hist)*length(x1))),hist_boxes);
    [n2,~] = hist(Final_Correct_Rost-Pyff,xout1);
    [n3,~] = hist(Final_Wrong_Rost-Pyff_Wrong,xout1);
    
    bar(xout1',[n2' n3'])
    colormap winter
    
    xlabel('Hedge Value - Payoff')
    ylabel('Frequency')
    legend('Correct model','Wrong Model')
    title(strcat('Distribution of hedging errors, censored at ',num2str(xout1(1)),' and ',num2str(xout1(length(xout1)))));
    
    saveas(gcf,strcat(dirpath,'Heston_Superhedge_Error_Compare.eps'),'psc2')    
    
    
    
    
    %% Compare effect of rho, Root & Rost:
    
    rhoVec = linspace(-1,1,NumRho);
    
    HPCellRho = cell(NumRho,1);
    RostCellRho = cell(NumRho,1);
    RootCellRho = cell(NumRho,1);
    
    HPVarCallRho = zeros(1,NumRho);
    RostVarCallRho = zeros(1,NumRho);
    RootVarCallRho = zeros(1,NumRho);

    FRho = @(t) max(t-RhoVarCallK,0);
    fRho = @(t) 1.0*(t>RhoVarCallK);
    
    for i = 1:NumRho
        HPCellRho{i} = HestonModel(s0,V0,xi,r,T_rho,kappa,theta,rhoVec(i));
        RostCellRho{i} = RostBarrierClass(HPCellRho{i},space_size,time_steps_Rost,time_step_changes_Rost,1);
        RootCellRho{i} = RootBarrierClass(HPCellRho{i},space_size,time_steps_Root,1);
        
        RostVarCallRho(i) = RostCellRho{i}.price(fRho);
        RootVarCallRho(i) = RootCellRho{i}.price(fRho);       
    end

    HPVarCallRho(1) = HPCellRho{1}.VaroptionPrice(FRho);
    HPVarCallRho(:) = HPVarCallRho(1);
    
    plot(rhoVec,HPVarCallRho,'Color',X(colOther,:))
    hold on;
    plot(rhoVec,RostVarCallRho,'Color',X(colRost,:))
    plot(rhoVec,RootVarCallRho,'Color',X(colRoot,:))
    
    %title('Heston Variance Call: Actual, Lower and Upper Bounds against rho')
    xlabel('rho')
    ylabel('Option Price')
    legend('Heston','Rost','Root')
        
    saveas(gcf,strcat(dirpath,'Heston_Var_Call_Bounds_in_rho.eps'),'psc2')

    hold off;
    
    RootCellRho{1}.PlotBarrier(1)
    hold on
    for i = 2:NumRho
        RootCellRho{i}.PlotBarrier(rem(i,length(X))+1)
    end

    xlabel('Integrated Variance')
    ylabel('Asset Price')
    ylim([xmin xmax]);
    saveas(gcf,strcat(dirpath,'Heston_Root_Barriers_in_rho.eps'),'psc2')

    %% Compare effect of time on Upper and Lower bounds
    
    TVec = linspace(T_rho_min,T_rho_max,NumTRho);
    
    HPCellTRho = cell(NumTRho,1);
    RostCellTRho = cell(NumTRho,1);
    RootCellTRho = cell(NumTRho,1);
    
    HPVarCallTRho = zeros(1,NumTRho);
    RostVarCallTRho = zeros(1,NumTRho);
    RootVarCallTRho = zeros(1,NumTRho);

    fCellTRho = cell(NumTRho,1);
    FCellTRho = cell(NumTRho,1);
    
    for i = 1:NumTRho
        FCellTRho{i} = @(t) max(t-RhoVarCallK*TVec(i),0);
        fCellTRho{i} = @(t) 1.0*(t>RhoVarCallK*TVec(i));
    
        HPCellTRho{i} = HestonModel(s0,V0,xi,r,TVec(i),kappa,theta,-1);
        RostCellTRho{i} = RostBarrierClass(HPCellTRho{i},space_size,time_steps_Rost,time_step_changes_Rost,1);
        RootCellTRho{i} = RootBarrierClass(HPCellTRho{i},space_size,time_steps_Root,1);
        
        RostVarCallTRho(i) = RostCellTRho{i}.price(fCellTRho{i});
        RootVarCallTRho(i) = RootCellTRho{i}.price(fCellTRho{i});       

        HPVarCallTRho(i) = HPCellTRho{i}.VaroptionPrice(FCellTRho{i});

    end

    hold off;
    
    plot(TVec,HPVarCallTRho,'Color',X(colOther,:))
    hold on;
    plot(TVec,RostVarCallTRho,'Color',X(colRost,:))
    plot(TVec,RootVarCallTRho,'Color',X(colRoot,:))
    
    %title('Heston Variance Call: Actual, Lower and Upper Bounds against rho')
    xlabel('Maturity')
    ylabel('Option Price')
    legend('Heston','Rost','Root')
        
    saveas(gcf,strcat(dirpath,'Heston_Var_Call_Bounds_in_time_fixed_rho.eps'),'psc2')

    hold off;
    
    RootCellTRho{1}.PlotBarrier(1)
    hold on
    for i = 2:NumTRho
        RootCellTRho{i}.PlotBarrier(rem(i,length(X))+1)
    end

    xlabel('Integrated Variance')
    ylabel('Asset Price')
    ylim([xmin xmax]);
    saveas(gcf,strcat(dirpath,'Heston_Root_Barriers_in_T_rho_is_minus_one.eps'),'psc2')

    
    %% Save data and close diary
    save(strcat(dirpath,'data.mat'));
    
    diary off
    
catch err
    display('Something went wrong!')
    
    diary off
    
    rethrow(err);
end

