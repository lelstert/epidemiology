function learn_Italy_epidemiology_data
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

    close 'all'
    
    format long e;
    scrsz = get(0,'ScreenSize');  % determine the screen size
    f     = 0.15;
    scrsz = [f*scrsz(3),f*scrsz(4),(1-2*f)*scrsz(3),(1-2*f)*scrsz(4)];

% ====================================================================================================================== 
    % starting data is Feb 21,2019  
    filename    = 'Italy_data.txt';
    formatSpec  = 'reading: %s\n';
    fprintf(formatSpec,filename);

    fileID      = fopen(filename);
    data_array  = dlmread(filename,',');  % reads the entire array, delimiter is ','
    T           = data_array(1:end,1); % double, read first column of data
    Infections  = data_array(1:end,2); % double, read second column of data
    Recovered   = data_array(1:end,3); % double, read third column of data
    Deaths      = data_array(1:end,4); % double, read fourth column of data 
    fclose(fileID);
    
% ====================================================================================================================== 
    % Tracking number of cases in Italy
    % https://gisanddata.maps.arcgis.com/apps/opsdashboard/index.html#/bda7594740fd40299423467b48e9ecf6
    
    lambda = log(Infections(end)/Infections(3))/(length(Infections)-3);
       
    figure('Name','Italy Data Infectied','NumberTitle','off','position',scrsz)
        semilogy(Infections,'bo')
        hold on
        semilogy(Infections(3)*exp(-3*lambda)*exp(lambda*T),'k--')
        hold off
        xlim([1,T(end)])
        xlabel('Days'),ylabel('Number Infected')
        title(['Infection Rate : ',num2str(lambda)])
        legend('Infected','Exponential Growth','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis
        
     figure('Name','Italy Data Recovered','NumberTitle','off','position',scrsz)
        semilogy(Recovered,'bo')
        xlim([1,T(end)])
        xlabel('Days'),ylabel('Number Recovered')
        title('Recovery Rate : ')
        legend('Recovered','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis 
        
     lambda = log(Deaths(end)/Deaths(4))/(length(Deaths)-4);
        
     figure('Name','Italy Data Deaths','NumberTitle','off','position',scrsz)
        semilogy(Deaths,'bo')
        hold on
        semilogy(Deaths(4)*exp(-4*lambda)*exp(lambda*T),'k--')
        hold off
        xlim([1,T(end)])
        xlabel('Days'),ylabel('Number Deaths')
        title(['Death Rate : ',num2str(lambda)])
        legend('Deaths','Exponential Growth','Location','northwest')
        set(gca,'FontSize',20) % refers to current axis      
        
end
