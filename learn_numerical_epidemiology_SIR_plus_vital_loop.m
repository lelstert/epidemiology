function learn_numerical_epidemiology_SIR_plus_vital_loop(problem)
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

% run the code, three options
% learn_numerical_epidemiology_SIR_plus_vital_loop('SIR_1')          % single case
% learn_numerical_epidemiology_SIR_plus_vital_loop('SIR_loop_1')     % sweep of N
% learn_numerical_epidemiology_SIR_plus_vital_loop('SIR_multiple_1') % multiple sites of infection, each starting at different times

    close 'all'
    
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)]; 
    
% ====================================================================================================================== 
    % Tracking number of cases in China
    % https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6
    % https://en.wikipedia.org/wiki/Timeline_of_the_2019%E2%80%9320_Wuhan_coronavirus_outbreak
    
    % [ ... ] January Dates up to 31, 32 is Feb 1
    Days      = [19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,...
                 43,44,45,46,47,48,49,50,51,52,53,54,55,56,57] - 4;

    Cases     = [198,291,440,571,830,1287,1975,2744,4515,5974,7711,9692,11791,14380,17205,20440,24234,28018,31161,34954,... 
                 37592,40645,43141,45206,59804,63851,66492,68500,70548,72436,74185,75778,76787,77932,78971,80350,81245,...
                 82594,83774]; %johns hopkins
    Recovered = [25,25,25,25,34,38,49,50,60,103,124,171,243,328,475,632,892,1153,1530,2391,... 
                 2922,3578,4228,5123,5911,6723,8096,9419,10844,12522,14376,16882,18864,21259,23384,27878,30311,...
                 33253,36654];% johns hopkins        
    Deaths    = [3,6,9,17,25,41,56,80,106,130,170,213,259,304,361,425,490,563,637,725,... 
                 814,910,1018,1118,1259,1380,1523,1669,1770,1868,2004,2130,2240,2362,2468,2705,2770,2810,2867]; 
         
    figure('Name','China Data','NumberTitle','off','position',scrsz)
        T    = linspace(0,Days(end),Days(end)+1)';
        % average value of exponential growth rates for first 10 points, small signal growth rate
        % this is used in initialize to set some coefficients so initial growth rate matches the data
        I_egr  = sum(log(Cases(1:10))./Days(1:10))/length(Days(1:10));
        % average over all data
        R_egr  = sum(log(Recovered)./Days)/length(Days);
        %average over all data
        D_egr  = sum(log(Deaths)./Days)/length(Days); 
        semilogy(Days,Cases,'bo')
        hold on
        semilogy(T,exp(I_egr*T),'b-')
        semilogy(Days,Recovered,'ro')
        semilogy(T,exp(R_egr*T),'r-')
        semilogy(Days,Deaths,'go')
        semilogy(T,exp(D_egr*T),'g-')
        hold off
        xlabel('days'),ylabel('Number Confirmed Cases, Recovered, and Deaths')
        title(['Case Rate : ',num2str(I_egr),' Recover Rate : ',num2str(R_egr),' Death Rate : ',num2str(D_egr)])
        legend('Case','Exponential Fit','Recover','Exponential Fit','Death','Exponential Fit','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis
      
% ======================================================================================================================    
    points   = 1000;
    loop_max = 1;
    if strcmp(problem,'SIR_loop_1'); loop_max = 10; end
    if strcmp(problem,'SIR_multiple_1'); loop_max = 4; end
    
    S  = zeros(points,loop_max);
    I  = zeros(points,loop_max);
    R  = zeros(points,loop_max);
    NC = zeros(points,loop_max);
    SM = zeros(points,loop_max);
    IM = zeros(points,loop_max);
    RM = zeros(points,loop_max);
    
    for loop = 1:loop_max
        
        [N,R0,alpha,beta,gamma,mu,e,eta,options,tspan,xi,label] = initialize(problem,loop,I_egr,D_egr); %#ok<ASGLU>
    
        % x(:,1) = S, x(:,2) = I, x(:,3) = R
        [t,x]          = ode15s(@ode_SIR_absolute,tspan,xi,options);
        nx             = x(:,1)+ x(:,2) + x(:,3);
        [T,S(:,loop)]  = equal_spacing(t,x(:,1),tspan,points); %#ok<ASGLU>
        [T,I(:,loop)]  = equal_spacing(t,x(:,2),tspan,points); %#ok<ASGLU>
        [T,R(:,loop)]  = equal_spacing(t,x(:,3),tspan,points); %#ok<ASGLU>
        [T,NC(:,loop)] = equal_spacing(t,x(:,1)+ x(:,2) + x(:,3),tspan,points); % number conservation
        
    end
    
% ======================================================================================================================    
    % check second conservation equation, valid only for birth and death rates of zero, and s = s_infinity 
    %R0_check = log(S(end,1)/S(1,1))/(S(end,1)/S(1,1)-1);
    
    %disp('check initial R0 to calculated R0 from model data')
    %disp('input R0')
    %disp(R0)
    %disp('R0 from model data')
    %disp(R0_check)
    

% ======================================================================================================================   
    if strcmp(problem,'SIR_1')
        
        [~,itoday] = min(abs(T - T(end)));
        
        figure('Name','SIR Solution','NumberTitle','off','position',scrsz)
   
            subplot(2,1,1)
            plot(t,x(:,1)/N,'k-') % This is the last loop value of x and N
            hold on
            plot(t,x(:,2)/N,'r-')
            plot(t,x(:,3)/N,'b-')
            plot(t,(nx(1,1)-nx(:,1)./exp(alpha*t))/nx(1,1),'g-')
            plot([T(itoday),T(itoday)],[0,1.2],'k--')
            hold off
            ylim([0,1.2])
            legend('Susceptible','Infectious','Recovered','Deaths')
            xlabel('days'),ylabel('Normalized Number')
            title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/(gamma+mu))])
            set(gca,'FontSize',10) % refers to current axis
            
            subplot(2,1,2)
            semilogy(t,x(:,1)/N,'k-')
            hold on
            semilogy(t,abs(x(:,2))/N,'r-')
            semilogy(t,x(:,3)/N,'b-')
            semilogy(t,(nx(1,1) - nx(:,1))/nx(1,1),'g-')
            %semilogy(t,(nx(1,1) - nx(:,1)./exp(alpha*t(:,1)))/nx(1,1),'g-')
            semilogy([T(itoday),T(itoday)],[1e-6,1e1],'k--')
            hold off
            ylim([1e-6,1e1])
            legend('Susceptible','Infectious','Recovered','Deaths','Location','southeast')
            xlabel('days'),ylabel('Normalized Number')
            title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/(gamma + mu))])
            set(gca,'FontSize',10) % refers to current axis
            
    end
    
% ======================================================================================================================        
    if strcmp(problem,'SIR_2')
        
        [~,itoday] = min(abs(T - T(end)));
        
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
            semilogy(Days,Cases,'ko')
            hold off
            ylim([1e0,1e6])
            legend('Infectious(1)','Infectious(2)','Infectious(3)','Infectious(4)','Infectious','Location','northeast')
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
            semilogy(Days,Recovered,'ko')
            semilogy(Days-8,3*Recovered,'ro')
            hold off
            ylim([1e0,1e6])
            legend('Recovered(1)','Recovered(2)','Recovered(3)','Recovered(4)','Recovered','Location','southeast')
            xlabel('days'),ylabel('Number Recovered')
            set(gca,'FontSize',20) % refers to current axis
      
    end 
     
% ======================================================================================================================    
    if strcmp(problem,'SIR_1') || strcmp(problem,'SIR_loop_1')
        
        [~,itoday] = min(abs(T - Days(end)));
    
        figure('Name','S Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,S,'k-')
            hold on
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            ylim([1e0,1e7])
            legend('Susceptible','Location','southeast')
            xlabel('days'),ylabel('Number Susceptible')
            title(['Model Susceptible : ',num2str(round(S(itoday)))])
            set(gca,'FontSize',20) % refers to current axis
         
        figure('Name','I Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,abs(I),'k-')
            hold on
            semilogy(Days,Cases,'bo')
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            legend('Model','Data','Location','northeast')
            ylim([1e0,1e7])
            xlabel('days'),ylabel('Number Infectious')
            title(['Model Infectious: ',num2str(round(I(itoday)))])
            set(gca,'FontSize',20) % refers to current axis
            
            for nn = 1:loop_max
                
                [v,i] = max(I(:,nn));
                v     = v/NC(1,nn)
            end
            
        figure('Name','R Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,R,'k-') % shifted because appears recovered not reported as quickly ????
            hold on
            semilogy(Days,Recovered,'bo')
            semilogy(Days-8,Recovered,'ro')
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            ylim([1e0,1e7])
            legend('Model','Data','Time Shifted Data')
            xlabel('days'),ylabel('Number Recovered')
            title(['Model Recovered: ',num2str(round(R(itoday)))]) 
            set(gca,'FontSize',20) % refers to current axis
            text(5,2e-2,'Data shifted backwards 8 days to match model time dependence - slower reporting?','FontSize',14)
    
        figure('Name','D Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,NC(1,:) - NC(:,:)./exp(alpha*T),'k-')
            %semilogy(T,mu*R,'k-')
            hold on
            semilogy(Days,Deaths,'bo')
            semilogy(Days-8,Deaths,'ro')
            semilogy([T(itoday),T(itoday)],[1e-1,1e7],'k--')
            hold off
            ylim([1e0,1e7])
            legend('Model','Data','Time Shifted Data')
            xlabel('days'),ylabel('Number Deaths')
            title(['Model Deaths : ',num2str(round(NC(1,1) - NC(itoday,1)/exp(alpha*T(itoday,1))))])  
            set(gca,'FontSize',20) % refers to current axis
            %text(5,2e-1,'Data shifted backwards 8 days to match model time dependence - slower reporting?','FontSize',14)
        
        figure('Name','R/I Solution Daily','NumberTitle','off','position',scrsz) 
            semilogy(T,R./I,'k-')
            hold on
            semilogy(Days,Recovered./Cases,'bo')
            semilogy([T(itoday),T(itoday)],[1e-10,1e10],'k--')
            hold off
            ylim([1e-2,1e1])
            legend('Model','Data')
            xlabel('days'),ylabel('Number Recovered/Number Infectious')
            title(['Model Recovered/Infectious : ',num2str(R(itoday)/I(itoday))]) 
            set(gca,'FontSize',20) % refers to current axis          
           
        figure('Name','Conservation Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T,NC,'k-')
            xlabel('days'),ylabel('Number Total')
            ylim([1e4,1e7])
            title(['Total Number of Individuals, Including Birth and Death :',' R0 = ',num2str(beta/(gamma+mu))])
            set(gca,'FontSize',20) % refers to current axis 
    end

% ======================================================================================================================                      
    function dxdt=ode_SIR_absolute(t,x) % ~ means ignore argument for this function
        
        % sum of these equations lead to dN/dt = (alpha - mu)*N when e = 0
        
        dxdt      =  zeros(3,1);
        N         = x(1)+x(2)+x(3);
        dxdt(1,1) =  alpha*N + e*N/(eta + t) - mu*x(1) - beta*x(2)*x(1)/N; % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - (gamma + mu)*x(2);                 % rate of change in infected
        dxdt(3,1) =  gamma*x(2) - mu*x(3);                                 % rate of change of recovered
      
    end

end

function [N,R0,alpha,beta,gamma,mu,e,eta,options,tspan,xi,label] = initialize(problem,loop,I_egr,D_egr) %#ok<INUSD>
      
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
            tspan  = [0,60];             % problem time interval in days
            xi     = [100000-1,1,0];     % x = [S,I,R] initial conditions
            N      = sum(xi);
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
            case 'SIR_2'
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
            tspan  = [0,150];             % problem time interval in days
            xi     = [N-1,1,0];      % x = [S,I,R] initial conditions
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
         otherwise
            warning('Unexpected Proplem Type. No Option Possible.');
       end   
     
end
   
function [xi,yi] = equal_spacing(x,y,tspan,n)
      
      % subfunction, interpolate nonuniform mesh onto a uniform mesh
      
      dx = (tspan(2)-tspan(1))/(n-1);
      xi = (tspan(1):dx:tspan(2))';
      yi = interp1(x,y,xi); 
      
end
   