function learn_numerical_epidemiology_SIR_plus_vital(problem)
% Compartmental models in epidemiology ( ep•i•de•mi•ol•o•gy )
% https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
% The SIR model with vital dynamics where death rate equals birth rate
% Tracking number of cases in China
% https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6

% S is the stock of susceptible population.
% I is the stock of infected.
% R is the stock of recovered population.
% N = S + I + R conservation of individuals.
% Sum of these equations lead to dN/dt = (alpha - mu)*N, where alpha is birth rate and mu is death rate
        

% typical time between contacts is 1/beta, see initialization
% typical time of recovery is 1/gamma, see initialization

    close 'all'
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)]; 
    
    % Number of Confirmed Cases
    % https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6
    Cases = [278,326,547,639,916,2000,2700,4400,5500];
    Days  = 5 + [10,11,12,13,14,15,16,17,18];
       
    [N,alpha,beta,gamma,mu,options,tspan,xi,label] = initialize(problem);
    
     % x(:,1) = S, x(:,2) = I, x(:,3) = R
     [t,x] = ode15s(@ode_SIR,tspan,xi,options);
     [T,S]  = equal_spacing(t,x(:,1),tspan); %#ok<ASGLU>
     [T,I]  = equal_spacing(t,x(:,2),tspan); %#ok<ASGLU>
     [T,R]  = equal_spacing(t,x(:,3),tspan); %#ok<ASGLU>
     [T,NC] = equal_spacing(t,x(:,1)+ x(:,2) + x(:,3),tspan); % number conservation
            
    figure('Name','S Solution Daily','NumberTitle','off','position',scrsz)
        semilogy(T,S,'k-')
        legend('Susceptible','Location','southwest')
        xlabel('days'),ylabel('Number Susceptible')
        title(label)
        set(gca,'FontSize',20) % refers to current axis
    
    figure('Name','I Solution Daily','NumberTitle','off','position',scrsz)
        semilogy(T,abs(I),'k-')
        hold on
        semilogy(Days,Cases,'k--')
        hold off
        legend('Infectious','Data')
        xlabel('days'),ylabel('Number Infectious')
        title(label)
        set(gca,'FontSize',20) % refers to current axis
        
    figure('Name','R Solution Daily','NumberTitle','off','position',scrsz)
        semilogy(T,R,'k-')
        legend('Recovered')
        xlabel('days'),ylabel('Number Recovered')
        title(label) 
        set(gca,'FontSize',20) % refers to current axis
        
    figure('Name','Conservation Solution Daily','NumberTitle','off','position',scrsz)
        plot(T,NC/N,'k-')
        ylim([0.9,1.1])
        xlabel('days'),ylabel('Number Total')
        title(label)
        set(gca,'FontSize',20) % refers to current axis 
       
    figure('Name','SIR Solution','NumberTitle','off','position',scrsz)
    
        subplot(2,1,1)
        plot(t,x(:,1)/N,'k-')
        hold on
        plot(t,x(:,2)/N,'r-')
        plot(t,x(:,3)/N,'b-')
        hold off
        legend('Susceptible','Infectious','Recovered')
        xlabel('days'),ylabel('Number')
        title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/gamma)])
        set(gca,'FontSize',10) % refers to current axis
    
        subplot(2,1,2)
        semilogy(t,x(:,1)/N,'k-')
        hold on
        semilogy(t,abs(x(:,2))/N,'r-')
        semilogy(t,x(:,3)/N,'b-')
        hold off
        ylim([1e-6,1e1])
        legend('Susceptible','Infectious','Recovered')
        xlabel('days'),ylabel('Number')
        title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/gamma)]) 
        set(gca,'FontSize',10) % refers to current axis
        
    figure('Name','Confirm Cases in China','NumberTitle','off','position',scrsz)
        semilogy(Days,Cases)
        xlabel('days'),ylabel('Number')
        title(['The Number of Confirmed Cases :',' R0 = ',num2str(beta/gamma)]) 
        set(gca,'FontSize',20) % refers to current axis
          
    function dxdt=ode_SIR(~,x) % ~ means ignore argument for this function
        
        % sum of these equations lead to dN/dt = (alpha - mu)*N
        
        dxdt      =  zeros(3,1);
        dxdt(1,1) = alpha*N - mu*x(1) - beta*x(2)*x(1)/N;          % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - (gamma + mu)*x(2);         % rate of change in infected
        dxdt(3,1) =  gamma*x(2) - mu*x(3);                         % rate of change of recovered
      
    end

end

function [N,alpha,beta,gamma,mu,options,tspan,xi,label] = initialize(problem)
      
      % numerical solver control
      tolerance = 1.0e-12; 
      options   = odeset('RelTol',tolerance,'AbsTol',[tolerance,tolerance,tolerance]);

      % parameters and initial conditions
      switch problem
        
          case 'test' 
            alpha = 0.04/365;           % birth rate per year
            beta  = 0;
            gamma = 0.050*15;      
            mu    = 0.24/365;           % death rate per year
            tspan = [0,365];            % problem time interval in days
            xi    = [10000000,1,0];     % x = [S,I,R] initial conditions
            N     = sum(xi);
            label = ['The SI Model With Vital Dynamics :',' R0 = ',num2str(beta/gamma)];
         
          case 'SIR'
            alpha = 0.04/365;           % birth rate per year
            beta  = 0.075*15;
            gamma = 0.050*15;      
            mu    = 0.02/365;           % death rate per year
            tspan = [0,30];            % problem time interval in days
            xi    = [10000000,1,0];     % x = [S,I,R] initial conditions
            N     = sum(xi);
            label = ['The SI Model With Vital Dynamics :',' R0 = ',num2str(beta/gamma)];
            
          case 'SIR_1'
            alpha = 0.04/365;           % birth rate per year
            beta  = 0.075*15;
            gamma = 0.050*15;      
            mu    = 0.02/365;           % death rate per year
            tspan = [0,22];            % problem time interval in days
            xi    = [10000000,1,0];     % x = [S,I,R] initial conditions
            N     = sum(xi);
            label = ['The SI Model With Vital Dynamics :',' R0 = ',num2str(beta/gamma)];
            
         otherwise
            warning('Unexpected Proplem Type. No Option Possible.');
       end   
     
end
   
function [xi,yi] = equal_spacing(x,y,tspan)
      
      % subfunction, interpolate nonuniform data into uniform data set
      
      dx = (tspan(2)-tspan(1))/(length(x)-1);
      xi = 0:dx:tspan(2);
      yi = interp1(x,y,xi); 
      
   end
