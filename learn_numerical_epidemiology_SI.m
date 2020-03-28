
function learn_numerical_epidemiology_SI(problem)
% Compartmental models in epidemiology 
% https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
% The SI model without vital dynamics
% Tracking number of cases in China
% https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6

% S is the stock of susceptible population.
% I is the stock of infected.
% N = S + I conservation of individuals.

% typical time between contacts is 1/beta, see initialization
% typical time of recovery is 1/gamma, see initialization

    close 'all'
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)];
    
    % Number of Confirmed Cases
    % https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6
    Cases = [278,326,547,639,916,2000,2700,4400,4600];
    Days  = 5 + [10,11,12,13,14,15,16,17,18];
    
    [N,beta,gamma,options,tspan,xi] = initialize(problem);
    
     % x(:,1) = S, x(:,2) = I
     [t,x] = ode15s(@ode_SI,tspan,xi,options);
     [T,S]  = equal_spacing(t,x(:,1),tspan); %#ok<ASGLU>
     [T,I]  = equal_spacing(t,x(:,2),tspan);  %#ok<ASGLU>
     [T,NC] = equal_spacing(t,x(:,1)+ x(:,2),tspan); % number conservation
            
    figure('Name','S Solution Daily','NumberTitle','off','position',scrsz)
        semilogy(T,S,'k-')
        legend('Susceptible')
        xlabel('days'),ylabel('Number')
        title(['The SI Model Without Vital Dynamics :',' R0 = ',num2str(beta/gamma)])
        set(gca,'FontSize',20) % refers to current axis
    
    figure('Name','I Solution Daily','NumberTitle','off','position',scrsz)
        semilogy(T,abs(I),'k-')
        hold on 
        semilogy(Days,Cases,'k--')
        hold off
        legend('Infectious','Data')
        xlabel('days'),ylabel('Number')
        title(['The SI Model Without Vital Dynamics :',' R0 = ',num2str(beta/gamma)])
        xlim([1,30])
        set(gca,'FontSize',20) % refers to current axis
        
    figure('Name','Conservation Solution Daily','NumberTitle','off','position',scrsz)
        plot(T,NC/N,'k-')
        ylim([0.5,1.5])
        xlabel('days'),ylabel('Number Conservation')
        title(['The SI Model Without Vital Dynamics :',' R0 = ',num2str(beta/gamma)])
        set(gca,'FontSize',20) % refers to current axis 
    
    figure('Name','SI Solution','NumberTitle','off','position',scrsz)
    
        subplot(2,1,1)
        plot(t,x(:,1)/N,'k-')
        hold on
        plot(t,x(:,2)/N,'r-')
        hold off
        legend('Susceptible','Infectious')
        xlabel('days'),ylabel('Number')
        title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/gamma)])
        set(gca,'FontSize',10) % refers to current axis
    
        subplot(2,1,2)
        semilogy(t,x(:,1)/N,'k-')
        hold on
        semilogy(t,abs(x(:,2))/N,'r-')
        hold off
        legend('Susceptible','Infectious')
        ylim([1e-6,1e1])
        xlabel('days'),ylabel('Number')
        title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/gamma)])
        set(gca,'FontSize',10) % refers to current axis
        
    figure('Name','Confirm Cases in China','NumberTitle','off','position',scrsz)
        semilogy(Days,Cases)
        xlabel('days'),ylabel('Number')
        title(['The Number of Confirmed Cases :',' R0 = ',num2str(beta/gamma)]) 
        set(gca,'FontSize',20) % refers to current axis
            
    function dxdt=ode_SI(~,x) % ~ means ignore argument for this function
        
        dxdt      =  zeros(2,1);
        dxdt(1,1) = -beta*x(2)*x(1)/N + gamma*x(2);              % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - gamma*x(2);              % rate of change in infected
                             
    end

end

function [N,beta,gamma,options,tspan,xi] = initialize(problem)
      
      % numerical solver control
      tolerance = 1.0e-12; 
      options   = odeset('RelTol',tolerance,'AbsTol',[tolerance,tolerance]);

      % parameters and initial conditions
      switch problem
       
         case 'SI'
            beta  = 0.075*15;
            gamma = 0.050*15;   
            tspan = [0,200];       % problem time interval in days
            xi    = [10000000,1]; % x = [S,I,R] initial conditions
            N     = sum(xi);
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
