%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%        The Growing Trend of Colorectal Cancer During Pregnancy:        %
%        A Study of Prolactin Signaling and Cancer Aggressiveness        %
%                                                                        %
%    Maria Lopez-Cavestany*, Olivia A. Wright*, Alexandria T. Carter*,   %
%           Brittany O'Brien^, Cathy Eng^, Michael R. King*              %
%                                                                        %
% *Department of Biomedical Engineering, Vanderbilt University,          %
%   Nashville, TN 37212                                                  %
%                                                                        %
% ^Division of Hematology Oncology, Vanderbilt University Medical Center,%
%   Nashville, TN 37212                                                  %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is th main file that initiates the Monte Carlo Simulation. It runs 
% the computational model n number of times with randomly calculated 
% initial protein values and kinetic constants.

% To sucessfully run this model you need the following scripts open:
%    > cell_1_SAmechanism.m


%DEFINING THE TIME SPAN

TimeSpan = [0:60:48*3600];

%DEFINING INITIAL CONCENTRATIONS

%Set non-zero values
%A lot of variables are in structures to make the workspace cleaner
protein.PRL = 120; %ng/mL
protein.PRL = protein.PRL*0.04545; %Calculate the concentration in nM
protein.RJ = 5.0802; 
protein.S3c = 117.1798;
protein.SHP2 = 78.23736788;
protein.PPX = 70.3604;
protein.PPN = 143.870996;
protein.pS3c = 0;
protein.S3n = 101.8617;
protein.JAG = 0.1528;

%Set zero values and fill in non-random. non-zero values
y0 = zeros(1,35);
y0(1) = protein.PRL;
y0(2) = protein.RJ;
y0(3) = protein.S3c;
y0(4) = protein.SHP2;
y0(5) = protein.PPX;
y0(6) = protein.PPN;
y0(20) = protein.S3n;
y0(34) = protein.JAG;

%Set up the kinetic constants
kin.k1 = 3.08E-03;               kin.k2f = 7.7971E-04; 
kin.k2r = 0.005;                 kin.k3f = 0.063415984; 
kin.k3r = 0.2;                   kin.k4 = 0.0052;
kin.k5f = 0.0124;                kin.k5r = 8.00E-01;
kin.k6 = 0.4417;                 kin.kdeg = 2.57E-04; 
kin.deg_ratio = 7.1442;          kin.k8f = 0.1814; 
kin.k8r = 1.00E-01;              kin.k9f = 0.000879762; 
kin.k9r = 2.00E-01;              kin.k10 = 3.00E-03;
kin.k11f = 0.0028;               kin.k11r = 2.00E-01;
kin.k12 =  0.0031;               kin.k13f = 2.00E-07; 
kin.k13r = 0.225430866;          kin.k14 = 0.005886461;
kin.k15f = 0.001460885;          kin.k15r = 2.00E-01;
kin.k16 = 0.0145287;             kin.k17f = 3.55E-02; 
kin.k17r = 0.062303131;          kin.k18_1 = 1.00E-02; 
kin.k18_2 = 4.00E+02;            kin.k19 = 0.001603988;
kin.k20 = 1.00E-02;              kin.k21f = 0.052168498; 
kin.k21r = 1.00E-01;             kin.k22 = 0.000506751;
kin.k23 = 6.8485E-04;            kin.k24 = 0.00124969;
kin.k25_1 = 0.0056;              kin.k25_2 = 4.00E+02;
kin.k26 = 1.00E-03;              kin.k27 = 6.12582E-05;
kin.k28 = 0.0132;                kin.k29 = 1.00E-02;
kin.k30_1 = 0.0021;              kin.k30_2 = 4.00E+02;
kin.k31 = 0.0011;                kin.k32 = 5.00E-04;
kin.k33 = 0.0099;                kin.k34 = 1.4518E-04;
kin.Vratio = 5.00E-01;
n = 50; %number of kinetic constants + 1 for the base calculations

%SET UP MATRICES FOR DATA STORAGE AND CALCULATIONS THROUGHOUT THE SIMULTAION
y0_stored = zeros(numel(y0), n);
pSTAT_ALL = zeros((numel(TimeSpan)), n); %The range is the TimeSpan - the lag time
JAG_ALL = zeros((numel(TimeSpan)), n);
pSTAT_RMSE = zeros(1,n);
JAG_RMSE = zeros(1, n);

%This loop runs the mechanisn script for each of the 
%substitute kinetic constant values to zero 1 by 1

for eye=1:n
    
    %LOAD THE KINETC CONSTANTS FRESH
    kinetics = [kin.k1 kin.k2f kin.k2r kin.k3f kin.k3r kin.k4 kin.k5f ...
    kin.k5r kin.k6 kin.kdeg kin.deg_ratio kin.k8f kin.k8r kin.k9f kin.k9r ...
    kin.k10 kin.k11f kin.k11r kin.k12 kin.k13f kin.k13r kin.k14 kin.k15f ...
    kin.k15r kin.k16 kin.k17f kin.k17r kin.k18_1 kin.k18_2 kin.k19 ...
    kin.k20 kin.k21f kin.k21r kin.k22 kin.k23 kin.k24 kin.k25_1 kin.k25_2 ...
    kin.k26 kin.k27 kin.k28 kin.k29 kin.k30_1 kin.k30_2 kin.k31 kin.k32 ...
    kin.k33 kin.k34 kin.Vratio];

    % For the baseline iteration (eye==1), leave all rates unchanged
    if eye == 1
        %do nothing

    % For denominators, use a small non-zero value; otherwise zero out
    else
        if ismember(eye, [29,38,44,49])
            kinetics(eye) = 1E-06;
        else
            kinetics(eye) = 0;
        end
    end

    
    %This here helps keep track of where you are in the loop
    disp(eye);

    %RUN ODE SOLVER
    options = odeset('RelTol',1e-5,'AbsTol',1e-6,'NonNegative',[1:length(y0)]);
    [~,y] = ode15s(@cell_1_SAMechanism, TimeSpan, y0, options, kinetics);
    
    %SAVE THE pSTAT3 AND JAG1 OUTPUTS
    pSTAT_ALL(:,eye) = y(:,11);
    JAG_ALL(:,eye) = y(:,34);

    %CALCULATE THE RMSE
    pSTAT_RMSE(eye) = sqrt(mean((pSTAT_ALL(:,eye) - pSTAT_ALL(:,1)).^2));
    JAG_RMSE(eye)   = sqrt(mean((JAG_ALL(:,eye)   - JAG_ALL(:,1)).^2));

    end
