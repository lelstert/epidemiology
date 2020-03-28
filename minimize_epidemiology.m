function minimize_epidemiology(filename)

    close 'all'
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)]; 
    
    data = read_data(filename,scrsz); 
     
    p = [1.5e2,1.0e-3,4.0e1,4.0e-3,2.0e-3,20.0];
   
    options   = optimset('Display','iter','PlotFcns',@optimplotfval, 'MaxFunEvals', 300,'TolFun',1.0e-5,'TolX',1.0e-5);   
    %options   = optimset('PlotFcns',@optimplotfval, 'MaxFunEvals', 1000,'TolFun',1.0e-5,'TolX',1.0e-5);
    %edit optimplotfval.m
    %optimset('fminsearch')
    handle  = @(p)parameterfun(p,data);
    [p_min,f_min] = fminsearch(handle,p,options);
    
    [cost,model] = epidemiology(p_min,data);
    
    figure('Name','Infected Data and Model','NumberTitle','off','position',scrsz)
        semilogy(data(:,1),data(:,2),'ko')
        hold on
        semilogy(model(:,1),model(:,2),'k-')
        hold off
        xlim([0,datal(end,1)])
        xlabel('Days'),ylabel('Number Infected')
        legend('Model','Data','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis
        
end
      
function data = read_data(filename,scrsz)
        
    formatSpec  = 'reading: %s\n';
    fprintf(formatSpec,filename)
    
    fileID       = fopen(filename);
    data         = dlmread(filename,',');  % reads the entire array, delimiter is ','
  
    fclose(fileID);  
    
    figure('Name','Data','NumberTitle','off','position',scrsz)
        semilogy(data(:,1),data(:,2),'ro')
        hold on
        semilogy(data(:,1),data(:,3),'go')
        semilogy(data(:,1),data(:,4),'ko')
        hold off
        xlim([0,data(end,1)])
        xlabel('Days'),ylabel('Number')
        legend('Infections','Recoveries','Deaths','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis
        
end

function [cost,varargout] = epidemiology(p,data)

    [x0,xi,m,k,r,a,t_center,start_amplitude,plimit] = initialize(p,tspan);
   
    options        = odeset('RelTol',1.0e-8,'AbsTol',[1.0e-8,1.0e-8,1.0e-8,1.0e-8]);  
    [ts,xs]        = ode15s(@(t,x) ode_release(t,x,p,m,k,r,t_center,start_amplitude),tspan,xi,options);
   
    cost= sum(abs(model(2,:) - data(2,:)))/numel(data(1,:));
    
    punish_factor = 1.0e3;
   	punish        = sum(p < plimit(1,:)) + sum(p >= plimit(2,:));		            			          
    if punish 
         cost = punish_factor*cost;
    end

    if nargout > 1
        varargout{1} = model;
    end

end

function y = parameterfun(p,data)
    y =epidemiology(p,data);
end

