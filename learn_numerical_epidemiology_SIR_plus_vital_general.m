function learn_numerical_epidemiology_SIR_plus_vital_general(problem,solver,filename)
% The SIR model with vital dynamics
% S is the stock of susceptible population.
% I is the stock of infected.
% R is the stock of recovered population.
% NC = S + I + R conservation of individuals, only when birth rate equals death rate, otherwise 
% Sum of these equations leads to dN/dt = (alpha - mu)*N, where alpha is birth rate and mu is death rate

% beta is average number of contacts between susceptible and infective, which leads to new infection, see initialization      
% 1/gamma is average time of infection, see initialization

% Compartmental models in epidemiology (ep•i•de•mi•ol•o•gy)
% https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
% http://home.iitk.ac.in/~peeyush/mth426/Epidemiology.pdf
% https://www.hindawi.com/journals/cmmm/2017/8613878/
% https://institutefordiseasemodeling.github.io/Documentation/general/model-si.html

% run the code with a number of options
% learn_numerical_epidemiology_SIR_plus_vital_loop('SIR_1','absolute','World_data.txt')          % single case
% learn_numerical_epidemiology_SIR_plus_vital_loop('SIR_loop_1','absolute_delay','Korea_data.txt')     % sweep of N
% learn_numerical_epidemiology_SIR_plus_vital_loop('SIR_multiple_1','absolute','Italy_data.txt') % multiple sites of infection, each starting at different times

 % Tracking number of cases in China, Italy, Korea
 % https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6
 % https://en.wikipedia.org/wiki/Timeline_of_the_2019%E2%80%9320_Wuhan_coronavirus_outbreak

    close 'all'
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)]; 
     
% ======================================================================================================================    
    second_conservation(scrsz)
  
% ====================================================================================================================== 
    formatSpec  = 'reading: %s\n';
    fprintf(formatSpec,filename);

    fileID      = fopen(filename);
    data_array  = dlmread(filename,',');  % reads the entire array, delimiter is ','
    TD          = data_array(1:end,1); % double, read first column of data
    Infections  = data_array(1:end,2); % double, read second column of data
    Recovered   = data_array(1:end,3); % double, read third column of data
    Deaths      = data_array(1:end,4); % double, read fourth column of data 
    fclose(fileID);
       
    [pI,pR,pD,aI,aR,aD] = align_data(filename);    
    [I_egr,R_egr,D_egr] = plot_data(TD,Infections,Recovered,Deaths,pI,pR,pD,scrsz);        
      
% ======================================================================================================================    
    points   = 1000;
    loop_max = 1;
    if strcmp(problem,'SIR_loop_1'); loop_max = 10; end
    if strcmp(problem,'SIR_multiple_1'); loop_max = 4; end
    
    S  = zeros(points,loop_max);
    I  = zeros(points,loop_max);
    R  = zeros(points,loop_max);
    D  = zeros(points,loop_max);
    NC = zeros(points,loop_max);
    SM = zeros(points,loop_max);
    IM = zeros(points,loop_max);
    RM = zeros(points,loop_max);
    DM = zeros(points,loop_max);
    
    for loop = 1:loop_max
        
        [N,R0,alpha,beta,gamma,mu,e,eta,options,tspan,xi,label] = initialize(problem,loop,I_egr,R_egr,D_egr); %#ok<ASGLU>
    
        % x(:,1) = S, x(:,2) = I, x(:,3) = R
        if strcmp(solver,'absolute')
            [t,x]      = ode15s(@ode_SIR_absolute,tspan,xi,options);
        end
        if strcmp(solver,'absolute_delay')
            [t,x]      = ode15s(@ode_SIR_absolute_delay,tspan,xi,options);
        end
        nx             = x(:,1)+ x(:,2) + x(:,3);
        [T,S(:,loop)]  = equal_spacing(t,x(:,1),tspan,points); %#ok<ASGLU>
        [T,I(:,loop)]  = equal_spacing(t,x(:,2),tspan,points); %#ok<ASGLU>
        [T,R(:,loop)]  = equal_spacing(t,x(:,3),tspan,points); %#ok<ASGLU>
        [T,D(:,loop)]  = equal_spacing(t,nx(1,1) - nx(:,1)./exp(alpha*t(:,1)),tspan,points); %#ok<ASGLU>
        [T,NC(:,loop)] = equal_spacing(t,x(:,1)+ x(:,2) + x(:,3),tspan,points); % number conservation
        %DM             = NC(1,:) - NC(:,:)./exp(alpha*T);
        
    end

% ======================================================================================================================   
    if strcmp(problem,'SIR_1'),plot_SIR_1(T,TD,N,t,x,nx,alpha,beta,gamma,mu,tspan,scrsz),end 
    if strcmp(problem,'SIR_1') || strcmp(problem,'SIR_loop_1')
        plot_SIR_1_SIR_loop_1(T,TD,S,I,Infections,R,Recovered,D,Deaths,NC,alpha,beta,gamma,mu,tspan,aI,aR,aD,scrsz)
    end
% ======================================================================================================================        
    if strcmp(problem,'SIR_2')
        
        [~,itoday] = min(abs(T - TD(end)));
        
        figure('Name','I Solution Multiple','NumberTitle','off','position',scrsz)
            semilogy(T,I(:,1),'r-')
            hold off
            ylim([1e0,1e4])
            legend('Infectious','Location','northeast')
            title(['Model Infectious: ',num2str(round(I(itoday,1)))])
            xlabel('days'),ylabel('Number Infectious')
            set(gca,'FontSize',20) % refers to current axis
            
    end
    
% ======================================================================================================================         
     if strcmp(problem,'SIR_multiple_1')
    
        days_index = round(10/(T(2)-T(1)));
    
        figure('Name','S Solution Multiple','NumberTitle','off','position',scrsz)
            SM(:,1)  = S(:,1);
            for j=2:loop_max
                SM(:,j)              = S(:,j);
                SM(:,j)              = circshift(SM(:,j),(j-1)*days_index);
                SM(1:(j-1)*days_index,j) = 0;
            end      
            SMS  = sum(SM,2);
            semilogy(T,SM(:,1),'r-')
            hold on
            semilogy(T,SM(:,2),'b-')
            semilogy(T,SM(:,3),'g-')
            semilogy(T,SM(:,4),'c-')
            semilogy(T,SMS,'k--')
            hold off
            ylim([1e0,1e6])
            legend('Susceptible(1)','Susceptible(2)','Susceptible(3)','Susceptible(4)','Susceptible','Location','southeast')
            title('Total Susceptible - Multiple Model')
            xlabel('days'),ylabel('Number Susceptible')
            set(gca,'FontSize',20) % refers to current axis
         
        figure('Name','I Solution Multiple','NumberTitle','off','position',scrsz)
           IM(:,1)  = I(:,1);
            for j=2:loop_max
                IM(:,j)              = I(:,j);
                IM(:,j)              = circshift(IM(:,j),(j-1)*days_index);
                IM(1:(j-1)*days_index,j) = 0;
            end      
            IMS  = sum(IM,2);
            semilogy(T,IM(:,1),'r-')
            hold on
            semilogy(T,IM(:,2),'b-')
            semilogy(T,IM(:,3),'g-')
            semilogy(T,IM(:,4),'c-')
            semilogy(T,IMS,'k--')
            semilogy(TD+aI,Infections,'ko') 
            hold off
            ylim([1e0,1e6])
            legend('Infectious(1)','Infectious(2)','Infectious(3)','Infectious(4)','Infectious','Location','northwest')
            title('Total Infected - Multiple Model')
            xlabel('days'),ylabel('Number Infectious')
            set(gca,'FontSize',20) % refers to current axis
    
        figure('Name','R Solution Multiple','NumberTitle','off','position',scrsz)
            RM(:,1)  = R(:,1);
            for j=2:loop_max
                RM(:,j)              = R(:,j);
                RM(:,j)              = circshift(RM(:,j),(j-1)*days_index);
                RM(1:(j-1)*days_index,j) = 0;
            end      
            RMS  = sum(RM,2);
            semilogy(T,RM(:,1),'r-')
            hold on
            semilogy(T,RM(:,2),'b-');
            semilogy(T,RM(:,3),'g-');
            semilogy(T,RM(:,4),'c-');
            semilogy(T,RMS,'k--'),
            semilogy(TD+aR,Recovered,'ko')
            semilogy(TD+aI,Infections,'ro')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e6])
            legend('Recovered(1)','Recovered(2)','Recovered(3)','Recovered(4)','Recovered Model','Recovered Data','Infections Data','Location','southeast')
            title('Total Recovered - Multiple Model')
            xlabel('days'),ylabel('Number Recovered')
            set(gca,'FontSize',20) % refers to current axis
      
    end 
   
% ======================================================================================================================    
    %  Generate test file
 
    %fileID = fopen('test_data.txt','w');
    
    model      = zeros(length(T),4);
    model(:,1) = T(:,1);
    model(:,2) = IM(:,1);
    model(:,3) = RM(:,1);
    model(:,4) = DM(:,1);
    dlmwrite('test_data.txt',model)
    
    %for n = 1:length(T)
        %formatSpec = '%24.15e %24.15e %24.15e %24.15e\n';
        %fprintf(fileID,formatSpec,T(n,1),IM(n,1),RM(n,1),DM(n,1)); 
    %end
    %fclose(fileID);
    
% ======================================================================================================================                      
    function dxdt=ode_SIR_absolute(t,x) % ~ means ignore argument for this function
        
        % sum of these equations lead to dN/dt = (alpha - mu)*N when e = 0
        
        dxdt      =  zeros(3,1);
        N         = x(1)+x(2)+x(3);
        dxdt(1,1) =  alpha*N + e*N/(eta + t) - mu*x(1) - beta*x(2)*x(1)/N; % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - (gamma + mu)*x(2); % rate of change in infected
        dxdt(3,1) =  gamma*x(2) - mu*x(3); % rate of change of recovered
      
    end

    function dxdt=ode_SIR_absolute_delay(t,x) % ~ means ignore argument for this function
        
        % sum of these equations lead to dN/dt = (alpha - mu)*N when e = 0
        
        dxdt      =  zeros(3,1);
        N         = x(1)+x(2)+x(3);
        H         = heaviside_general_rise(t-10,0.5,2);
        dxdt(1,1) =  alpha*N + e*N/(eta + t) - H*mu*x(1) - beta*x(2)*x(1)/N; % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - H*(gamma + mu)*x(2);
        dxdt(3,1) =  H*(gamma*x(2) - mu*x(3)); % rate of change of recovered
    end

end

function [N,R0,alpha,beta,gamma,mu,e,eta,options,tspan,xi,label] = initialize(problem,loop,I_egr,R_egr,D_egr) %#ok<INUSD>
      
      % numerical solver control
      tolerance = 1.0e-12; 
      options   = odeset('RelTol',tolerance,'AbsTol',[tolerance,tolerance,tolerance]);

      % parameters and initial conditions
      switch problem
        
          case 'SIR_1'
            e      = 0;                  % switch for increase N due to area increase, crude diffusion concept of spreading
            eta    = 2000/30;           % initial radius divided by expansion velocity
            alpha  = 0.0115/365;         % 11.5 per 1000 per year, birth rate per year, if too large(travel increase population size) infection can be cyclic
            % Initial exponential growth rate is beta - (gamma + mu) = beta*(R0 - 1)/R0
            % R0   = beta/(gamma+mu)
            R0     = 3; % best guess at this point in time, cannot really determine until epidemic has run its course
            beta   = I_egr*R0/(R0-1);
            mu     = 15*0.008/365;       % natural death rate per year 8 per 1000 per year is base rate, 15 times natural death rate
            gamma  = beta/R0 - mu;
            tspan  = [0,100];             % problem time interval in days
            xi     = [100000-1,1,0];     % x = [S,I,R] initial conditions
            N      = sum(xi);
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
            case 'SIR_2' % ship
            % approximaetly 621 infected in about 20 days
            e      = 0;                  % switch for increase N due to area increase, crude diffusion concept of spreading
            eta    = 2000/30;           % initial radius divided by expansion velocity
            alpha  = 0.0115/365;         % 11.5 per 1000 per year, birth rate per year, if too large(travel increase population size) infection can be cyclic
            % Initial exponential growth rate is beta - (gamma + mu) = beta*(R0 - 1)/R0
            % R0   = beta/(gamma+mu)
            R0     = 3; % best guess at this point in time, cannot really determine until epidemic has run its course
            beta   = I_egr*R0/(R0-1);
            mu     = 15*0.008/365;       % natural death rate per year 8 per 1000 per year is base rate, 15 times natural death rate
            gamma  = beta/R0 - mu;
            tspan  = [0,20];             % problem time interval in days
            xi     = [3700-1,1,0];     % x = [S,I,R] initial conditions
            N      = sum(xi);
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
          case 'SIR_loop_1'
       
            N      = 100000 + (loop -1)*100000; 
            e      = 0;                  % switch for increase due to area increase
            eta    = 1000/100;           % initial radius divided by expansion velocity
            alpha  = 0.0115/365;         % 11.5 per 1000 per year, birth rate per year, if too large(travel increase population size) infection can be cyclic
            % Initial exponential growth rate is beta - (gamma + mu) = beta*(R0 - 1)/R0
            % R0     = beta/(gamma+mu)
            R0     = 3;
            beta   = I_egr*R0/(R0-1);
            mu     = 15*0.008/365;       % death rate per year 8 per 1000 per year
            gamma  = beta/R0 - mu;
            tspan  = [0,60];             % problem time interval in days
            xi     = [N-1,1,0];      % x = [S,I,R] initial conditions
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
          case 'SIR_multiple_1'
              
            %N0     = [100000,150000,200000,250000];
            N0     = [100000,120000,180000,200000];
            N      = N0(loop); 
            e      = 0;                  % switch for increase due to area increase
            eta    = 1000/100;           % initial radius divided by expansion velocity
            alpha  = 0.0115/365;         % 11.5 per 1000 per year, birth rate per year, if too large(travel increase population size) infection can be cyclic
            % Initial exponential growth rate is beta - (gamma + mu) = beta*(R0 - 1)/R0
            % R0     = beta/(gamma+mu)
            R0     = 3;
            beta   = I_egr*R0/(R0-1);
            mu     = 15*0.008/365;       % death rate per year 8 per 1000 per year
            gamma  = beta/R0 - mu;
            tspan  = [0,80];             % problem time interval in days
            xi     = [N-1,1,0];      % x = [S,I,R] initial conditions
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
         otherwise
            warning('Unexpected Proplem Type. No Option Possible.');
       end   
     
end

function [pI,pR,pD,aI,aR,aD] = align_data(filename)
    
    if strcmp(filename,'China_data.txt')
        pI    = [1,10];
        pR    = [5,12];
        pD    = [14,28];
        aI = 12;
        aR = 12;
        aD = 0; 
    end
    
    if strcmp(filename,'Korea_data.txt')
        pI    = [1,10];
        pR    = [30,35];
        pD    = [30,35];
        aI = 0;
        aR = 0;
        aD = 0; 
    end
    
    if strcmp(filename,'Italy_data.txt')
        pI = [3,11];
        pR = [3,11];
        pD = [3,11];
        aI = 12;
        aR = 0;
        aD = 0; 
    end
    
    if strcmp(filename,'World_data.txt')
        pI    = [1,10];
        pR    = [5,12];
        pD    = [14,28];
        aI = 12;
        aR = 0;
        aD = 0; 
    end
    end
   
function [xi,yi] = equal_spacing(x,y,tspan,n)
      
      % subfunction, interpolate nonuniform mesh onto a uniform mesh
      
      dx = (tspan(2)-tspan(1))/(n-1);
      xi = (tspan(1):dx:tspan(2))';
      yi = interp1(x,y,xi); 
      
end

function [H] = heaviside_general_rise(arg,start_amplitude,t_rise)
   % finite heaviside function for turn on of inelastic collision
   shift_o = -tanh((start_amplitude-0.5D0)/0.5D0)*t_rise;
   H       = 0.5D0*(1.0D0 + tanh((arg-shift_o)/t_rise));
end

function second_conservation(scrsz)
% second conservation equation, valid only for birth and death rates of zero

    S_infinity = linspace(0.001,0.999,1000);
    R0_infinity = log(S_infinity)./(S_infinity-1);
    
    figure('Name','S_{infinity} vs R_{0}','NumberTitle','off','position',scrsz)
        plot(S_infinity,R0_infinity,'k-')
        xlabel('S_{infinity}'),ylabel('R_{0}')
        title('Second Conservation Law, Valid for Birth and Death Rate of Zero')
        set(gca,'FontSize',20) % refers to current axis
end

function [I_egr,R_egr,D_egr] = plot_data(TD,Infections,Recovered,Deaths,pI,pR,pD,scrsz)
% https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6
    
    I_egr = log(Infections(pI(2))/Infections(pI(1)))/(length(Infections(pI(1):pI(2)))-1);  
    figure('Name','Data Infectied','NumberTitle','off','position',scrsz)
        semilogy(TD,Infections,'bo')
        hold on
        semilogy(TD,Infections(pI(1))*exp(-pI(1)*I_egr)*exp(I_egr*TD),'k--')
        hold off
        xlim([0,TD(end)])
        xlabel('Days'),ylabel('Number Infected')
        title(['Infection Rate : ',num2str(I_egr)])
        legend('Infected','Exponential Growth','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis
        
     R_egr = log(Recovered(pR(2))/Recovered(pR(1)))/(length(Recovered(pR(1):pR(2)))-1);   
     figure('Name','Data Recovered','NumberTitle','off','position',scrsz)
        semilogy(TD,Recovered,'bo')
        hold on
        semilogy(TD,Recovered(pR(1))*exp(-pR(1)*R_egr)*exp(R_egr*TD),'k--')
        hold off
        xlim([0,TD(end)])
        xlabel('Days'),ylabel('Number Recovered')
        title(['Recovery Rate : ',num2str(R_egr)])
        legend('Recovered','Exponential Growth','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis 
        
     D_egr = log(Deaths(pD(2))/Deaths(pD(1)))/(length(Deaths(pD(1):pD(2)))-1);
       
     figure('Name','Data Deaths','NumberTitle','off','position',scrsz)
        semilogy(TD,Deaths,'bo')
        hold on
        semilogy(TD,Deaths(pD(1))*exp(-pD(1)*D_egr)*exp(D_egr*TD),'k--')
        hold off
        xlim([0,TD(end)])
        xlabel('Days'),ylabel('Number Deaths') 
        title(['Death Rate : ',num2str(D_egr)])
        legend('Deaths','Exponential Growth','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis  
        
end
        
function plot_SIR_1(T,TD,N,t,x,nx,alpha,beta,gamma,mu,tspan,scrsz)

        [~,itoday] = min(abs(T - TD(end)));
        
        figure('Name','SIR Solution','NumberTitle','off','position',scrsz)
   
            subplot(2,1,1)
            plot(t,x(:,1)/N,'k-') % This is the last loop value of x and N
            hold on
            plot(t,x(:,2)/N,'r-')
            plot(t,x(:,3)/N,'b-')
            plot(t,(nx(1,1)-nx(:,1)./exp(alpha*t))/nx(1,1),'g-')
            plot([T(itoday),T(itoday)],[0,1.2],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([0,1.2])
            legend('Susceptible','Infectious','Recovered','Deaths')
            xlabel('days'),ylabel('Normalized Number')
            title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/(gamma+mu))])
            set(gca,'FontSize',10) % refers to current axis
            
            subplot(2,1,2)
            semilogy(t,x(:,1)/N,'k-')
            hold on
            semilogy(t,x(:,2)/N,'r-')
            semilogy(t,x(:,3)/N,'b-')
            semilogy(t,(nx(1,1) - nx(:,1))/nx(1,1),'g-')
            %semilogy(t,(nx(1,1) - nx(:,1)./exp(alpha*t(:,1)))/nx(1,1),'g-')
            semilogy([T(itoday),T(itoday)],[1e-6,1e1],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e-6,1e1])
            legend('Susceptible','Infectious','Recovered','Deaths','Location','southeast')
            xlabel('days'),ylabel('Normalized Number')
            title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/(gamma + mu))])
            set(gca,'FontSize',10) % refers to current axis
            
end
    
function plot_SIR_1_SIR_loop_1(T,TD,S,I,Infections,R,Recovered,D,NC,Deaths,alpha,beta,gamma,mu,tspan,aI,aR,aD,scrsz) 

    [~,itoday] = min(abs(T - TD(end)));
    
        figure('Name','S Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,S,'k-')
            hold on
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e6])
            legend('Susceptible','Location','southeast')
            xlabel('days'),ylabel('Number Susceptible')
            title(['Model Susceptible : ',num2str(round(S(itoday)))])
            set(gca,'FontSize',20) % refers to current axis
         
        figure('Name','I Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,I,'k-')
            hold on
            semilogy(TD+aI,Infections,'bo')
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            legend('Model','Data','Location','northeast')
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e6])
            xlabel('days'),ylabel('Number Infectious')
            title(['Model Infectious: ',num2str(round(I(itoday)))])
            set(gca,'FontSize',20) % refers to current axis
            
        figure('Name','R Solution Daily','NumberTitle','off','position',scrsz) % shifted because appears re
            semilogy(T,I,'k-')
            hold on
            semilogy(T,R,'k--')
            semilogy(TD,Recovered,'ko')
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e6])
            %legend('Model','Data','Time Shifted Data')
            legend('Model I','Model R','Data R','Location','southeast')
            xlabel('days'),ylabel('Number Recovered')
            title(['Model Recovered: ',num2str(round(R(itoday)))]) 
            set(gca,'FontSize',20) % refers to current axis
            text(5,2e-2,'Data shifted backwards 8 days to match model time dependence - slower reporting?','FontSize',14)
    
        figure('Name','D Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,I,'k-')
            hold on
            semilogy(T,D,'r--')
            semilogy(T,R,'b-')
            semilogy(TD,Deaths,'ko')  
            semilogy([T(itoday),T(itoday)],[1e-1,1e7],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e6])
            legend('Model I','Model D','Model D','Data D')
            xlabel('days'),ylabel('Number Deaths')
            title(['Model Deaths : ',num2str(round(NC(1,1) - NC(itoday,1)/exp(alpha*T(itoday,1))))])  
            set(gca,'FontSize',20) % refers to current axis
            %text(5,2e-1,'Data shifted backwards 8 days to match model time dependence - slower reporting?','FontSize',14)
        
        figure('Name','R/I Solution Daily','NumberTitle','off','position',scrsz) 
            semilogy(T,R./I,'k-')
            hold on
            semilogy(TD,Recovered./Infections,'bo')
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e-2,1e1])
            legend('Model','Data')
            xlabel('days'),ylabel('Number Recovered/Number Infectious')
            title(['Model Recovered/Infectious : ',num2str(R(itoday)/I(itoday))]) 
            set(gca,'FontSize',20) % refers to current axis          
           
        figure('Name','Conservation Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,NC,'k-')
            xlabel('days'),ylabel('Number Total')
            xlim([tspan(1),tspan(2)]),ylim([1e4,1e7])
            title(['Total Number of Individuals, Including Birth and Death :',' R0 = ',num2str(beta/(gamma+mu))])
            set(gca,'FontSize',20) % refers to current axis 
    
            end