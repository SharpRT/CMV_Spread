
clearvars
clc

%% PARAM SETUP
path='C:\';
disp(['Path = ' path])

description = '';
disp(['description = ' description])

model_isMutantSwitch = 1;
tau = 0;
rhoRep = 0;
model_isFDDT = 1;
lambda_i = 0.00768;

% control type (parameter to vary)
%{
    ParamControlType        Value  
    NoControl               0       
    beta/lambda             1
    mu/rho                  2
    Dv (vectorDisp)         3       bi br
    plantingRate(ResHost)   4       sigma
    P                       5
    alpha(insecticide)      6
    theta(CSS)              7
    p (selectivePlant)      8
    zeta (trade)            9
    DT (tradeDisp)          10      bzeta
    gamma                   11
    r (vectorBirth)         12      b
    K (vectorCapacity)      13
    m                       14
    omega                   15
    eta_lambda              16
    eta_gamma               17
%}   
paramControlType = 1;

% control method (how is control applied?)
%{
    Control Type            Value
    No Control              0
    Always/Everywhere       1
    Time-dependent          2
    Periodic                3
    Density-Depedent        4    
    Spatial                 5
%}
controlType = 0;

adjustSpatial = 1; % halve spatial domain in mutant case.

% domain setup
width = 20000*1000;
nx = 2^(10); %2^10 = 1024
dx = width/nx;

if model_isMutantSwitch && adjustSpatial
    width = width / 2;
    nx = nx / 2;
end

x = linspace(-width/2,width/2-dx,nx);
ndays = 365.25*100; %365.25*50;

controlPressure = 0.5;  
controlStartTime = 0; %365.25*20;
controlRunTime = inf;

% controlType=3 parameters
if controlType == 3
    if ~model_isMutantSwitch
       ndays = ndays / 2; 
    end
    controlStartTime = 365.25*20;
end
maxApplications = (ndays-controlStartTime)/365.25;
periodTime = 365.25; %time off and time on - full period = 2*periodTime

% controlType=4 parameters
lower_control_threshold = 5e-2;%5e-2; % 0.05;
upper_control_threshold = inf;   

% controlType=5 parameters
spatialThresholdKM = 1000;
if model_isMutantSwitch && adjustSpatial
    spatialThresholdKM = 1000 / 2;
end
spatialControlLoc = zeros(1,nx);
spatialControlLoc(abs(x)>spatialThresholdKM*1000) = 1;

% model parameters
muLambda = 1.25;

sigma_min = 0;  %standard: 0.002
sigma_max = 0.003;
omega_min = 0.002;
omega_max = 0.004;
lambda_min = 0; %standard: 0.002
lambda_max = 0.032;
rho_min = 0;
rho_max = 0.033;
r_min = 0.1;
r_max = 0.3;
alpha_min = 0.06;
alpha_max = r_max;   %standard: 0.18
gamma_min = 0;  %standard: 0.002
gamma_max = 0.016;
K_min = 0;
K_max = 250;
theta_min = 0;
theta_max = 1;
p_min = 0;
p_max = 1;
zeta_min = 0;
zeta_max = 1;
m_min = 0;
m_max = 1;
D_m_min = 0;
D_m_max = 5000;
D_zeta_min = 0;
D_zeta_max = 60000;

sigma = 0.003;
omega = 0.003;

lambda_r = lambda_i/muLambda;
rho_i = 0.003;
rho_r = 0.003;
gamma_i = 0.004;
gamma_r = 0.004;

theta = 0.05;
zeta = 0.5;
p = 0.4;

P=50;
alpha_var = 0.12;
r=0.2;
K=P*r/(r-alpha_var);
P = @(alpha,r,K) (r-alpha).*K./r;

D_m = 1000;
D_zeta = 30000;

m = 0.025;

%% HOUSEKEEPING

if(model_isMutantSwitch)
    filenameInvaderType = '_Mutant';
else
    filenameInvaderType = '_Migrant';
end

switch paramControlType %init control_string & eps_C
    case 0 %NO_Control
        stringParamControlType = '_NO_Control';
    case 1 %beta_Controlled / lambda_Controlled
        stringParamControlType = '_lambda_Controlled';
        stringParamControlName = 'lambda';
        eps_beta_C = controlPressure;
    case 2 %mu_Controlled / rho_Controlled
        stringParamControlType = '_rho_Controlled';
        stringParamControlName = 'rho';
        eps_mu_C = controlPressure;
    case 3 %D_Controlled
        stringParamControlType = '_D_Controlled';
        stringParamControlName = 'D';            
        eps_D_C = controlPressure;
    case 4 %plantRate_Controlled
        stringParamControlType = '_plantRate_Controlled';
        stringParamControlName = 'sigma';
        eps_sigma_C = controlPressure;
    case 5 %P_Controlled
        stringParamControlType = '_P_Controlled';
        stringParamControlName = 'P';
        eps_P_C = controlPressure;
    case 6 %alpha_Controlled
        stringParamControlType = '_alpha_Controlled';
        stringParamControlName = 'alpha';
        eps_alpha_C = controlPressure;
    case 7 %theta(CSS)_Controlled
        stringParamControlType = '_theta_Controlled';
        stringParamControlName = 'theta';
        eps_theta_C = controlPressure;
    case 8 %p_Controlled
        stringParamControlType = '_p_Controlled';
        stringParamControlName = 'p';
        eps_p_C = controlPressure;
    case 9 %zeta_Controlled
        stringParamControlType = '_zeta_Controlled';
        stringParamControlName = 'zeta';
        eps_zeta_C = controlPressure;
    case 10 %Dzeta_Controlled
        stringParamControlType = '_Dzeta_Controlled';
        stringParamControlName = 'Dzeta';
        eps_Dzeta_C = controlPressure;
    case 11 %gamma_Controlled
        stringParamControlType = '_gamma_Controlled';
        stringParamControlName = 'gamma';
        eps_gamma_C = controlPressure;
    case 12 %r_Controlled
        stringParamControlType = '_r_Controlled';
        stringParamControlName = 'r';
        eps_r_C = controlPressure;
    case 13 %K_Controlled
        stringParamControlType = '_K_Controlled';
        stringParamControlName = 'K';
        eps_K_C = controlPressure;
    case 14 %Dmigrate_Controlled
        stringParamControlType = '_m_Controlled';
        stringParamControlName = 'm';
        eps_Dmigrate_C = controlPressure;
    case 15 %omega_Controlled
        stringParamControlType = '_omega_Controlled';
        stringParamControlName = 'omega';
        eps_omega_C = controlPressure;
    case 16 %_Controlled
        stringParamControlType = '_etaLambda_Controlled';
        stringParamControlName = 'etaLambda';
        eps_etaLambda_C = controlPressure;
    case 17 %omega_Controlled
        stringParamControlType = '_etaGamma_Controlled';
        stringParamControlName = 'etaGamma';
        eps_etaGamma_C = controlPressure;
    otherwise
        error('Invalid paramControlType - terminating');
end
disp(['stringParamControlType = ' stringParamControlType])

switch controlType
    case 0 %No Control
        filenameControlType = '_-_NoControl';
        min_on = -1;
        max_on = -1;        
    case 1 %Always-Everywhere 
        filenameControlType = '_-_Always-Everywhere';
        min_on = 0;
        max_on = inf;        
    case 2 %Time Dependent Control  
        if ndays < controlStartTime+controlRunTime
            endTime = ndays;
        else
            endTime = controlStartTime+controlRunTime;
        end        
        filenameControlType = strcat('_-_TimeDependent(t=',num2str(controlStartTime),...
            '_endTime=',num2str(endTime),')');
        min_on = controlStartTime;
        max_on = min_on + controlRunTime;
    case 3 %Periodic Control
        filenameControlType = strcat('_-_Periodic(every_t=',num2str(periodTime*2),...
            '_after_t=', num2str(controlStartTime),...
            '_maxApplications=',num2str(maxApplications),')');
        min_on = (zeros(1,maxApplications)+controlStartTime)+(0:2:2*(maxApplications-1))*periodTime;
        max_on = min_on+periodTime;        
    case 4 %Density Dependent Control
        filenameControlType = strcat('_-_DensityDependent(t=',num2str(controlStartTime)...
            ,', between_Ii=',num2str(lower_control_threshold),'-',num2str(upper_control_threshold),')');
        min_on = controlStartTime;
        max_on = min_on + controlRunTime;     
    case 5 %Spatial
        filenameControlType = '_-_SpatialControl';
        min_on = 0;
        max_on = inf;          
    otherwise %Implement No Control
        error('Invalid paramControlType - terminating');        
end
disp(['filenameControlType = ' filenameControlType])
  
%% INITIAL-CONDITIONS/EQUILIBRIA
if model_isMutantSwitch            
    [HstarSoln,IrstarSoln,ZrstarSoln]=sim_initSolve(...
        theta,p,sigma,omega,rho_r,lambda_r,alpha_var,tau,gamma_r,P,r,K,rhoRep...
    );
else
    HstarSoln = sigma/omega;
    YstarSoln=P(alpha_var,r,K);
end

if model_isMutantSwitch
    Ii = zeros(1,nx);
    Ir = IrstarSoln * ones(1,nx);
    H = HstarSoln * ones(1,nx);

    Zr = [ZrstarSoln * ones(1,nx/2), 0, ZrstarSoln * ones(1,-1+nx/2)];
    Zi = ZrstarSoln - Zr;      
else
    Ii = zeros(1,nx);
    H = HstarSoln * ones(1,nx);
    Zi = [zeros(1,nx/2), P(alpha_var,r,K) * 0.01, zeros(1,-1+nx/2)];

    Ir = zeros(1,nx);
    Zr = zeros(1,nx);
end
  
H0 = H';
Ir0 = Ir';
Ii0 = Ii';
Zr0 = Zr';
Zi0 = Zi';

init = [H0; Ir0; Ii0; Zr0; Zi0];

%% SOLVE

nt = 20000;
tRange = 0:(ndays/(nt-1)):ndays;

f = @(t,N) myode(...
    t,N,nx,dx,x,model_isMutantSwitch,paramControlType,rhoRep,model_isFDDT,...
    sigma,omega,lambda_i,lambda_r,rho_i,rho_r,gamma_i,gamma_r,theta,zeta,p,alpha_var,r,K,m,D_m,D_zeta,...
    sigma_min,lambda_min,gamma_min,zeta_min,r_min,K_min,m_min,D_m_min,D_zeta_min,...
    omega_max,rho_max,theta_max,p_max,alpha_max,...
    controlType,controlPressure,min_on,max_on,lower_control_threshold,upper_control_threshold,spatialControlLoc...
);

disp(' ')
disp('Running ode45. Measuring runtime...')
tic

[t,N] = ode45(f,tRange,init);    

toc
disp(' ')

%% CALCULATE SPEED DATA

threshold = 1e-5;
x_centerIndex = find(x==0);
Ii_all = N(:,(2*nx+x_centerIndex):(3*nx));
nx_pos = (3*nx)-(2*nx+x_centerIndex)+1;
timeX = NaN(nx_pos,1);
rateX = NaN(nx_pos,1);

%get the time that threshold was exceeded at x(i)
for i=1:nx_pos    
    Ii_xi = Ii_all(:,i);
    
    tmp = t(find(Ii_xi>threshold,1));
    if isempty(tmp)
        break;
    else
        timeX(i) = tmp;
    end
end

%code to handle case where threshold exceeds more than one cell in a time step
numCons = diff([0 find(diff(timeX')) nx_pos]);
numRem = NaN(1,nx_pos);
numConsCount = NaN(1,nx_pos);
numCount = 1;
for i = 1:numel(numCons)
    for j = 1:numCons(i)
       numRem(numCount) = numCons(i) - j;
       numConsCount(numCount) = j - 1;
       numCount = numCount+1;
    end            
end

start_tx = [NaN (1:(nx_pos)-1)-numConsCount(2:nx_pos)];
end_tx = (1:nx_pos)+numRem(1:nx_pos);

%calculate wave speed at x(i)
for i=1:nx_pos
    if i==1 || start_tx(i)<=0
        rateX(i) = NaN;
    else
        rateX(i) = ((end_tx(i)-start_tx(i))*dx / (timeX(end_tx(i)) - timeX(start_tx(i)))); %convert to km/year
    end
end

timeX = timeX / 365.25; % convert to years
rateX = rateX * 365.25 / 1000; % convert to km/year

%% CALCULATE CONTROL MATRIX

nt_sample = 200;

control_matrix = ones(nx+1,nt_sample+1)*-1;
control_matrix(1,1) = NaN;
control_matrix(2:end,1) = x;
control_matrix(1,2:end) = t(1:(nt/nt_sample):nt)/365.25;
sampleCount = 1;
%Constructing Control Matrix
for i = 1:(nt/nt_sample):nt
    Nti = N(i,:);
    Ii_ti = Nti((2*nx+1):(3*nx))';
    eps_control_vector = control_vector(...
        t(i), nx, controlType, controlPressure, min_on, max_on,...
        lower_control_threshold, upper_control_threshold, Ii_ti...
    );    
    control_matrix(2:end,sampleCount+1) = eps_control_vector;
    sampleCount = sampleCount + 1;
end

%% SAVE DATA

filename_data = strcat(stringParamControlType,filenameControlType,filenameInvaderType,description,'.csv');
disp(['filename_data = ' filename_data])
tstr1 = {strcat(stringParamControlName,'Time') strcat(stringParamControlName,'Rate')};
tstr2 = strrep(cellstr(num2str(threshold','%.0e')),' ','');
timeStr = strcat(tstr1(1),{'_'},tstr2);
rateStr = strcat(tstr1(2),{'_'},tstr2);
outStr = [timeStr'; rateStr'];
finalStr = outStr(:)';
header={'x',finalStr{:},stringParamControlName};

fid = fopen(strcat(path,filename_data), 'w');
fprintf(fid, strcat(sprintf('%s',header{1}),sprintf(',%s',header{2:end}),'\n'));
fclose(fid);

Ntn = N(end,:);
%H_tn = Ntn((0*nx+1):(1*nx))';
%Ir_tn = Ntn((1*nx+1):(2*nx))';
Ii_tn = Ntn((2*nx+1):(3*nx))';
%Zr_tn = Ntn((3*nx+1):(4*nx))';
%Zi_tn = Ntn((4*nx+1):(5*nx))';

mainDataOut = [x',[NaN(x_centerIndex-1,1); timeX],[NaN(x_centerIndex-1,1); rateX],Ii_tn];
dlmwrite(strcat(path,filename_data), mainDataOut, '-append');

filename_controlMatrix = strcat(...
    stringParamControlType,filenameControlType,filenameInvaderType,'ControlMatrix',description,'.csv'...
);
disp(['filename_controlMatrix = ' filename_controlMatrix])
csvwrite(strcat(path,filename_controlMatrix),control_matrix)
