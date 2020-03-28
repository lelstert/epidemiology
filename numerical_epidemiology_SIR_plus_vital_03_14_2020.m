
function numerical_epidemiology_SIR_plus_vital_03_14_2020(problem,solver,data_file_name)

    close 'all'
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)]; 
    
    [data,I_egr,R_egr,D_egr,align] = read_data(data_file_name,scrsz); 
    
    % ======================================================================================================================    
    points   = 1000;
    loop_max = 1;
    if strcmp(problem,'SIR_loop_1'); loop_max = 10; end
    if strcmp(problem,'SIR_multiple_1'); loop_max = 4; end
  
    S_model  = zeros(points,loop_max);
    I_model  = zeros(points,loop_max);
    R_model  = zeros(points,loop_max);
    D_model  = zeros(points,loop_max);
    N_model  = zeros(points,loop_max);
    
    for loop = 1:loop_max
        
        [N,R0,alpha,beta,gamma,mu,e,eta,options,tspan,xi,label] = initialize(problem,loop,I_egr,R_egr,D_egr,data_file_name);  %#ok<ASGLU>
    
        % x(:,1) = S, x(:,2) = I, x(:,3) = R
        if strcmp(solver,'absolute')
            [t,x]      = ode15s(@ode_SIR_absolute,tspan,xi,options);
            nx         = x(:,1)+ x(:,2) + x(:,3);
        end
        if strcmp(solver,'absolute_delay')
            [t,x]      = ode15s(@ode_SIR_absolute_delay,tspan,xi,options);
            nx         = x(:,1)+ x(:,2) + x(:,3);
        end
        
        [T_model,S_model(:,loop)]  = equal_spacing(t,x(:,1),tspan,points); %#ok<ASGLU>
        [T_model,I_model(:,loop)]  = equal_spacing(t,x(:,2),tspan,points); %#ok<ASGLU>
        [T_model,R_model(:,loop)]  = equal_spacing(t,x(:,3),tspan,points); %#ok<ASGLU>
        [T_model,D_model(:,loop)]  = equal_spacing(t,nx(1,1) - nx(:,1)./exp(alpha*t(:,1)),tspan,points); %#ok<ASGLU>
        [T_model,N_model(:,loop)]  = equal_spacing(t,nx,tspan,points); % number conservation
              
    end
    
    plot_model(data,t,x,nx,T_model,S_model,I_model,R_model,D_model,N_model,alpha,beta,gamma,mu,tspan,scrsz)
    
    if strcmp(problem,'SIR_1') || strcmp(problem,'SIR_multiple_1') 
        plot_compare(data,T_model,S_model,I_model,R_model,D_model,N_model,alpha,beta,gamma,mu,tspan,loop_max,align,scrsz)
    end
    
    if strcmp(problem,'SIR_loop_1') 
        plot_loop(data,T_model,S_model,I_model,R_model,D_model,N_model,alpha,beta,gamma,mu,tspan,loop_max,align,scrsz) 
    end
    
    function dxdt=ode_SIR_absolute(t,x) 
        
        % sum of these equations lead to dN/dt = (alpha - mu)*N when e = 0
        
        dxdt      =  zeros(3,1);
        N         = x(1)+x(2)+x(3);
        dxdt(1,1) =  alpha*N + e*N/(eta + t) - mu*x(1) - beta*x(2)*x(1)/N; % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - (gamma + mu)*x(2); % rate of change in infected
        dxdt(3,1) =  gamma*x(2) - mu*x(3); % rate of change of recovered
      
    end

    function dxdt=ode_SIR_absolute_delay(t,x) 
        
        % sum of these equations lead to dN/dt = (alpha - mu)*N when e = 0
        
        dxdt      =  zeros(3,1);
        N         = x(1)+x(2)+x(3);
        H         = heaviside_general_rise(t-10,0.5,2);
        dxdt(1,1) =  alpha*N + e*N/(eta + t) - H*mu*x(1) - beta*x(2)*x(1)/N; % rate of change in susceptible
        dxdt(2,1) =  beta*x(2)*x(1)/N - H*(gamma + mu)*x(2);
        dxdt(3,1) =  H*(gamma*x(2) - mu*x(3)); % rate of change of recovered
    end    
    
end

function [p,align] = align_data(filename)
    
    if strcmp(filename,'China_data.txt') 
        p     = [1,10; 5,12; 14,28];
        align = [19; 0; 0]; 
    end
    
     if strcmp(filename,'Italy_data.txt')
        p     = [3,8; 10,15; 6,20];
        align = [0; 0; 0];  
    end
    
    if strcmp(filename,'Korea_data.txt')
        p     = [1,10; 5,12; 14,28];
        align = [0; 0; 0];  
    end
    
    if strcmp(filename,'World_data.txt')
       p     = [1,10; 5,12; 14,28];
       align = [0; 0; 0];  
    end
    
 end
      
function [data,I_egr,R_egr,D_egr,align] = read_data(filename,scrsz)
        
    formatSpec  = 'reading: %s\n';
    fprintf(formatSpec,filename)
    
    directory = 'C:\Users\LesWork\Documents\Mathematics_Physics_Science\MATLAB\epidemiology\learn_epidemiology\data\';   
    
    fileID       = fopen([directory,filename]);
    data         = dlmread([directory,filename],',');  % reads the entire array, delimiter is ','
    fclose(fileID); 
    
    [p,align] = align_data(filename);
    
    I_egr = log(data(p(1,2),2)/data(p(1,1),2))/(length( data( p(1,1):p(1,2), 2) ) - 1);
    R_egr = log(data(p(2,2),3)/data(p(2,1),3))/(length( data( p(2,1):p(2,2), 3) ) - 1);    
    D_egr = log(data(p(3,2),4)/data(p(3,1),4))/(length( data( p(3,1):p(3,2), 4) ) - 1);
    
    figure('Name','Data','NumberTitle','off','position',scrsz)
        semilogy(data(:,1),data(:,2),'ro')
        hold on
        semilogy(data(:,1),data(:,3),'go')
        semilogy(data(:,1),data(:,4),'ko')
        semilogy(data(:,1),data(p(1,1),2)*exp(-p(1,1)*I_egr)*exp(I_egr*data(:,1)),'k--')
        semilogy(data(:,1),data(p(2,1),3)*exp(-p(2,1)*R_egr)*exp(R_egr*data(:,1)),'k--')
        semilogy(data(:,1),data(p(3,1),4)*exp(-p(3,1)*D_egr)*exp(D_egr*data(:,1)),'k--')
        hold off
        xlim([0,data(end,1)])
        xlabel('Days'),ylabel('Number')
        title(['I: ',num2str(I_egr),', R: ',num2str(R_egr),', D: ',num2str(D_egr)])
        legend('Infections','Recoveries','Deaths','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis
        
end

function [N,R0,alpha,beta,gamma,mu,e,eta,options,tspan,xi,label] = initialize(problem,loop,I_egr,~,~,data_file_name) 
      
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
            tspan  = [0,100];             % problem time interval in days
            xi     = [N-1,1,0];      % x = [S,I,R] initial conditions
            label = ['The SIR Model With Vital Dynamics :',' R0 = ',num2str(beta/(gamma+mu))];
            
          case 'SIR_multiple_1'
              
            if strcmp(data_file_name,'Italy_data.txt')  
                N0     = [100000,150000,200000,250000];
            end
            if strcmp(data_file_name,'China_data.txt')
                N0     = [100000,120000,180000,200000];
            end
            if strcmp(data_file_name,'Korea_data.txt')
                N0     = [100,50000,180000,200000];
            end
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

function plot_model(data,t,x,nx,T_model,S_model,I_model,R_model,D_model,N_model,alpha,beta,gamma,mu,tspan,scrsz)

        [~,itoday] = min(abs(T_model - data(end,1) ));
        
        figure('Name','SIR Solution','NumberTitle','off','position',scrsz)
   
            subplot(2,1,1)
            plot(T_model,S_model(:,1),'k-') 
            hold on
            plot(T_model,I_model(:,1),'r-')
            plot(T_model,R_model(:,1),'b-')
            plot(T_model,D_model(:,1),'c-')
            plot(T_model,S_model,'k-') 
            plot(T_model,I_model,'r-')
            plot(T_model,R_model,'b-')
            plot(T_model,D_model,'c-')
            plot([T_model(itoday),T_model(itoday)],[1.0e-0,N_model(1,end)],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1.0e0,N_model(1,end)])
            legend('Susceptible','Infectious','Recovered','Deaths')
            xlabel('days'),ylabel('Normalized Number')
            title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/(gamma+mu))])
            set(gca,'FontSize',10) % refers to current axis
            
            subplot(2,1,2)
            semilogy(t,x(:,1)/nx(1,1),'k-')
            hold on
            semilogy(t,x(:,2)/nx(1,1),'r-')
            semilogy(t,x(:,3)/nx(1,1),'b-')
            semilogy(t,(nx(1,1) - nx(:,1)./exp(alpha*t(:,1)))/nx(1,1),'c-')
            semilogy([T_model(itoday),T_model(itoday)],[1e-6,1e0],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e-6,1e0])
            legend('Susceptible','Infectious','Recovered','Deaths','Location','southeast')
            xlabel('days'),ylabel('Normalized Number')
            title(['Average Number of Contacts Before Infected Individual has Recovered :',' R0 = ',num2str(beta/(gamma + mu))])
            set(gca,'FontSize',10) % refers to current axis
            
end

function plot_compare(data,T_model,S_model,I_model,R_model,D_model,N_model,alpha,beta,gamma,mu,tspan,loop_max,align,scrsz) %#ok<INUSL>
    
        days_index = round(10/(T_model(2)-T_model(1)));
    
        figure('Name','S Solution Multiple','NumberTitle','off','position',scrsz)
            SM(:,1)  = S_model(:,1);
            if loop_max > 1
                for j=2:loop_max
                    SM(:,j)              = S_model(:,j);
                    SM(:,j)              = circshift(SM(:,j),(j-1)*days_index);
                    SM(1:(j-1)*days_index,j) = 0;
                end      
            end
            SMS  = sum(SM,2);
            semilogy(T_model,SMS,'k-')
            hold on
            if loop_max > 1
                semilogy(T_model,SM(:,2),'b-')
                semilogy(T_model,SM(:,3),'g-')
                semilogy(T_model,SM(:,4),'c-')
            end
            semilogy(T_model,SM(:,1),'r-')
            semilogy(T_model,SMS,'k-')
            %semilogy(data(:,1)+align(1),N_model(1,1)-data(:,2),'ko')
            hold off
            ylim([1e0,1e6])
            legend('Susceptible','Location','southeast')
            title('Total Susceptible - Model')
            xlabel('days'),ylabel('Number Susceptible')
            set(gca,'FontSize',20) % refers to current axis
         
        figure('Name','I Solution Multiple','NumberTitle','off','position',scrsz)
           IM(:,1)  = I_model(:,1);
            if loop_max > 1
                for j=2:loop_max
                    IM(:,j)              = I_model(:,j);
                    IM(:,j)              = circshift(IM(:,j),(j-1)*days_index);
                    IM(1:(j-1)*days_index,j) = 0;
                end      
            end
            IMS  = sum(IM,2);
            semilogy(T_model,IMS,'k--')
            hold on
            if loop_max > 1
                semilogy(T_model,IM(:,2),'b-')
                semilogy(T_model,IM(:,3),'g-')
                semilogy(T_model,IM(:,4),'c-')
                %semilogy(data(:,1)+align(1),data(:,2),'ko')
                semilogy(data(:,1)+align(1),data(:,2)-data(:,3),'ko')
            end
            semilogy(T_model,IM(:,1),'k--')
            hold off
            ylim([1e0,1e6])
            legend('Infected','Location','southeast')
            title('Total Infected - Model')
            xlabel('days'),ylabel('Number Infected')
            set(gca,'FontSize',20) % refers to current axis
    
        figure('Name','R Solution Multiple','NumberTitle','off','position',scrsz)
            RM(:,1)  = R_model(:,1);
            if loop_max > 1
                for j=2:loop_max
                    RM(:,j)              = R_model(:,j);
                    RM(:,j)              = circshift(RM(:,j),(j-1)*days_index);
                    RM(1:(j-1)*days_index,j) = 0;
                end      
            end
            RMS  = sum(RM,2);
            semilogy(T_model,RMS,'k--')
            hold on
            if loop_max > 1
                semilogy(T_model,RM(:,2),'b-')
                semilogy(T_model,RM(:,3),'g-')
                semilogy(T_model,RM(:,4),'c-')
                semilogy(data(:,1)+align(2),data(:,3),'ko')
            end
            semilogy(T_model,RM(:,1),'k--')
            hold off
            ylim([1e0,1e6])
            legend('Recovered','Location','southeast')
            title('Total Recovered - Model')
            xlabel('days'),ylabel('Number Recovered')
            set(gca,'FontSize',20) % refers to current axis

end

function plot_loop(data,T_model,S_model,I_model,R_model,D_model,N_model,alpha,beta,gamma,mu,tspan,loop_max,align,scrsz) %#ok<INUSL>

    [~,itoday] = min(abs(T_model - data(end,1) ));
    
    figure('Name','S Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T_model,S_model,'k-')
            hold on
            semilogy([T_model(itoday),T_model(itoday)],[1e-10,1e10],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e7])
            legend('Susceptible','Location','southeast')
            xlabel('Days'),ylabel('Number Susceptible')
            title(['Model Susceptible (Today) : ',num2str(round(S_model(itoday)))])
            set(gca,'FontSize',15) % refers to current axis
         
        figure('Name','I Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T_model,I_model,'k-')
            hold on
            semilogy(data(:,1)+align(1),data(:,2),'bo')
            semilogy([T_model(itoday),T_model(itoday)],[1e-10,1e10],'k--')
            hold off
            legend('Model','Data','Location','northeast')
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e7])
            xlabel('Days'),ylabel('Number Infectious')
            title(['Model Infectious (Today) : ',num2str(round(I_model(itoday)))])
            set(gca,'FontSize',15) % refers to current axis
            
        figure('Name','R Solution Daily','NumberTitle','off','position',scrsz) % shifted because appears re
            semilogy(T_model,R_model(:,1),'k-')
            hold on
            semilogy(data(:,1),data(:,3),'ko')
            semilogy(T_model,R_model,'k-')
            semilogy([T_model(itoday),T_model(itoday)],[1e-10,1e10],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e7])
            legend('Model R','Data R','Location','southeast')
            xlabel('Days'),ylabel('Number Recovered')
            title(['Model Recovered (Today) : ',num2str(round(R_model(itoday)))]) 
            set(gca,'FontSize',15) % refers to current axis
            %text(5,2e-2,'Data shifted backwards 8 days to match model time dependence - slower reporting?','FontSize',14)
    
        figure('Name','D Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T_model,D_model(:,1),'k-')
            hold on
            semilogy(data(:,1),data(:,4),'ko') 
            semilogy(T_model,D_model,'k-')  
            semilogy([T_model(itoday),T_model(itoday)],[1e-1,1e7],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e0,1e6])
            legend('Model D','Data D','Location','southeast')
            xlabel('Days'),ylabel('Number Deaths')
            title(['Model Deaths (Today) : ',num2str(round(N_model(1,1) - N_model(itoday,1)/exp(alpha*T_model(itoday,1))))])  
            set(gca,'FontSize',15) % refers to current axis
            %text(5,2e-1,'Data shifted backwards 8 days to match model time dependence - slower reporting?','FontSize',14)
        
        figure('Name','R/I Solution Daily','NumberTitle','off','position',scrsz) 
            semilogy(T_model,R_model(:,1)./I_model(:,1),'k-')
            hold on
            semilogy(data(:,1),data(:,3)./data(:,2),'ko')
            semilogy(T_model,R_model./I_model,'k-')
            semilogy([T_model(itoday),T_model(itoday)],[1e-10,1e10],'k--')
            hold off
            xlim([tspan(1),tspan(2)]),ylim([1e-10,1e10])
            legend('Model','Data','Location','southeast')
            xlabel('Days'),ylabel('Number Recovered/Number Infectious')
            title(['Recovered/Infectious (Today) : ',num2str(R_model(itoday)/I_model(itoday))]) 
            set(gca,'FontSize',15) % refers to current axis          
           
        figure('Name','Conservation Solution Daily','NumberTitle','off','position',scrsz)
            semilogy(T_model,N_model,'k-')
            xlabel('Days'),ylabel('Number Total')
            xlim([tspan(1),tspan(2)]),ylim([1e4,1e7])
            title(['Total Number of Individuals, Including Birth and Death :',' R0 = ',num2str(beta/(gamma+mu))])
            set(gca,'FontSize',15) % refers to current axis 
    
            end