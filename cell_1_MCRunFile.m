%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%        Elevated Prolactin Signaling During Pregnancy Increases         %
%                   Colorectal Cancer Aggressiveness                     %
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
%    > cell_1_MCmechanism.m
%    > cell_1_MCFigureSTAT;
%    > cell_1_MCFigureJAG;


tic

%DEFINING THE TIME SPAN
LagTime = 0.05; %The lag time is 3 mins for model start up
LagTimeMins = LagTime*60; %What is the lag time in minutes

%Set up the time points we are testing to incorporate full run time 
% plus a lag/start up time
TimeSpan = [0:60:((48+LagTime)*3600)];


%DEFINE THE NUMBER OF ITERATIONS
 n=30000;


%DEFINING INITIAL CONCENTRATIONS

%Set non-zero values
%A lot of variables are in structures to make the workspace cleaner
protein.PRL = 120; %ng/mL
protein.PRL = protein.PRL*0.04545; %Calculate the concentration in nM
protein.RJ = 5.261347242; 
protein.S3c = 109.8519799;
protein.SHP2 = 78.23736788;
protein.PPX = 71.47583834;
protein.PPN = 143.870996;
protein.pS3c = 0;
protein.S3n = 101.6786441;
protein.JAG = 0.201537129;

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

%Set up the kinetics constants
kin.k1 = 3.08E-03;               kin.k2f = 0.000739219; 
kin.k2r = 0.004788793;           kin.k3f = 0.063415984; 
kin.k3r = 0.2;                   kin.k4 = 0.00589277;
kin.k5f = 0.012191872;           kin.k5r = 8.00E-01;
kin.k6 = 0.418065333;            kin.kdeg = 2.57E-04; 
kin.deg_ratio = 7.188464301;     kin.k8f = 0.146852819; 
kin.k8r = 1.00E-01;              kin.k9f = 0.000879762; 
kin.k9r = 2.00E-01;              kin.k10 = 3.00E-03;
kin.k11f = 0.002808389;          kin.k11r = 2.00E-01;
kin.k12 =  0.002883157;          kin.k13f = 2.00E-07; 
kin.k13r = 0.225430866;          kin.k14 = 0.005886461;
kin.k15f = 0.001460885;          kin.k15r = 2.00E-01;
kin.k16 = 0.0145287;             kin.k17f = 3.55E-02; 
kin.k17r = 0.062303131;          kin.k18_1 = 1.00E-02; 
kin.k18_2 = 4.00E+02;            kin.k19 = 0.001603988;
kin.k20 = 1.00E-02;              kin.k21f = 0.052168498; 
kin.k21r = 1.00E-01;             kin.k22 = 0.000506751;
kin.k23 = 0.000724421;           kin.k24 = 0.00124969;
kin.k25_1 = 0.006177181;         kin.k25_2 = 4.00E+02;
kin.k26 = 1.00E-03;              kin.k27 = 6.12582E-05;
kin.k28 = 0.013334245;           kin.k29 = 1.00E-02;
kin.k30_1 = 0.002090295;         kin.k30_2 = 4.00E+02;
kin.k31 = 0.001118064;           kin.k32 = 5.00E-04;
kin.k33 = 0.008438882;           kin.k34 = 0.000144097;
kin.Vratio = 5.00E-01;

kinetics = [kin.k1 kin.k2f kin.k2r kin.k3f kin.k3r kin.k4 kin.k5f ...
    kin.k5r kin.k6 kin.kdeg kin.deg_ratio kin.k8f kin.k8r kin.k9f kin.k9r ...
    kin.k10 kin.k11f kin.k11r kin.k12 kin.k13f kin.k13r kin.k14 kin.k15f ...
    kin.k15r kin.k16 kin.k17f kin.k17r kin.k18_1 kin.k18_2 kin.k19 ...
    kin.k20 kin.k21f kin.k21r kin.k22 kin.k23 kin.k24 kin.k25_1 kin.k25_2 ...
    kin.k26 kin.k27 kin.k28 kin.k29 kin.k30_1 kin.k30_2 kin.k31 kin.k32 ...
    kin.k33 kin.k34 kin.Vratio];

%Set up random variation in protein expression and binding kinetics
Cmax = 0.2; %Maximum percent change, was 1 for MC1 and 0.2 for MC2
random.randRJ = abs(protein.RJ - (Cmax*protein.RJ*(randn(1,n)/3)));
random.randPPX = abs(protein.PPX - (Cmax*protein.PPX*(randn(1,n)/3)));
random.randS3c = abs(protein.S3c - (Cmax*protein.S3c*(randn(1,n)/3)));
random.randS3n = abs(protein.S3n - (Cmax*protein.S3n*(randn(1,n)/3)));
random.randJAG = abs(protein.JAG - (Cmax*protein.JAG*randn(1,n)));
random.randk2f = abs(kin.k2f - (Cmax*kin.k2f*(randn(1,n)/3)));
random.randk2r = abs(kin.k2r - (Cmax*kin.k2r*(randn(1,n)/3)));
random.randk4 = abs(kin.k4 - (Cmax*kin.k4*(randn(1,n)/3)));
random.randk5f = abs(kin.k5f - (Cmax*kin.k5f*(randn(1,n)/3)));
random.randk6 = abs(kin.k6 - (Cmax*kin.k6*(randn(1,n)/3)));
random.randk8f = abs(kin.k8f - (Cmax*kin.k8f*(randn(1,n)/3)));
random.randk11f = abs(kin.k11f - (Cmax*kin.k11f*(randn(1,n)/3)));
random.randk12 = abs(kin.k12 - (Cmax*kin.k12*(randn(1,n)/3)));
random.randk23 = abs(kin.k23 - (Cmax*kin.k23*(randn(1,n)/3)));
random.randk25_1 = abs(kin.k25_1 - (Cmax*kin.k25_1*(randn(1,n)/3)));
random.randk28 = abs(kin.k28 - (Cmax*kin.k28*(randn(1,n)/3)));
random.randk30_1 = abs(kin.k30_1 - (Cmax*kin.k30_1*(randn(1,n)/3)));
random.randk31 = abs(kin.k31 - (Cmax*kin.k31*(randn(1,n)/3)));
random.randk33 = abs(kin.k33 - (Cmax*kin.k33*(randn(1,n)/3)));
random.randk34 = abs(kin.k34 - (Cmax*kin.k34*(randn(1,n)/3)));
random.randdeg_ratio = abs(kin.deg_ratio - (Cmax*kin.deg_ratio*(randn(1,n)/3)));


%SET UP MATRICES FOR DATA STORAGE AND CALCULATIONS THROUGHOUT THE SIMULTAION
y0_stored = zeros(numel(y0), n);
pSTAT_ALL = zeros((numel(TimeSpan)-LagTimeMins), n); %The range is the TimeSpan - the lag time
JAG_ALL = zeros((numel(TimeSpan)-LagTimeMins), n);

pSTAT_RMSE = zeros(1,n);
JAG_RMSE = zeros(1, n);

plots.keepingtrack = [];
kinetics_stored = zeros(49,n);

plots.pSTAT_InVitro = [1.47 1.14 0.99 1.06];
plots.pSTAT_InVitro_t = [5 15 30 60];
plots.JAG_InVitro = [1.93 2.01 1.8];
plots.JAG_InVitro_t = [12 24 48];


%This loop runs the mechanisn script for each of the 
%randomly generated protein expression and kinetic constant values

for eye=1:n
    
    %Fill in random values
    y0(2) = random.randRJ(1,eye); 
    y0(5) = random.randPPX(1,eye);
    y0(3) = random.randS3c(1,eye); 
    y0(20) = random.randS3n(1,eye);
    y0(34) = random.randJAG(1,eye);

    kinetics(2) = random.randk2f(1,eye);
    kinetics(3) = random.randk2r(1,eye);
    kinetics(6) = random.randk4(1,eye);
    kinetics(7) = random.randk5f(1,eye);
    kinetics(9) = random.randk6(1,eye);
    kinetics(12) = random.randk8f(1,eye);
    kinetics(17) = random.randk11f(1,eye);
    kinetics(19) = random.randk12(1,eye);
    kinetics(35) = random.randk23(1,eye);
    kinetics(37) = random.randk25_1(1,eye);
    kinetics(41) = random.randk28(1,eye);
    kinetics(43) = random.randk30_1(1,eye);
    kinetics(45) = random.randk31(1,eye);
    kinetics (47) = random.randk33(1,eye);
    kinetics (48) = random.randk34(1,eye);
    kinetics(11) = random.randdeg_ratio(1,eye);

    kinetics_stored (:,eye) = kinetics.'; %Stores the random kinetic values for each of the model iterations
    y0_stored (:,eye) = y0.'; %Stores the model outputs for each of the model iterations

    
    %This helps keep track of where you are in the loop
    disp(eye);
    plots.keepingtrack = [plots.keepingtrack eye];


    %RUN ODE SOLVER

    options = odeset('RelTol',1e-7,'AbsTol',1e-7,'NonNegative',[1:length(y0)]);
    [~,y] = ode15s(@cell_1_MCMechanism, TimeSpan, y0, options, kinetics);

    
    % REMOVE THE LAG TIME

    y(1:LagTimeMins,:)=[]; %remove the lag time

    
    %SAVE THE pSTAT3 AND JAG1 OUTPUTS

    pSTAT_ALL(:,eye) = y(:,11); %Stores the pSTAT3 values for each iteration in a separate table
    JAG_ALL(:,eye) = y(:,34); %Stores the JAG1 values for each iteration in a separate table

    
    %CALCULATE THE MODEL FITS
    
    calculate.pSTAT_Model = [y(6,11) y(16,11) y(31,11) y(61,11)]; %This identifies the model outputs for pSTAT3 to compare to our in vitro results
    calculate.JAG_Model = [y(721,34) y(1441,34) y(2881,34)]; %This identifies the model outputs for total JAG1 (intracellular and extracellular) to compare to our in vitro results

    % If you are interested you can also calculate the RMSE
    calculate.pSTAT_RMSE = sqrt(mean((calculate.pSTAT_Model - plots.pSTAT_InVitro).^2)); %RMSE calculation for pSTAT3
    pSTAT_RMSE(1,eye) = calculate.pSTAT_RMSE; %Stores the pSTAT3 RMSE values for each iteration in a separate table
    calculate.JAG_RMSE = sqrt(mean((calculate.JAG_Model - plots.JAG_InVitro).^2)); %RMSE calculation for JAG1
    JAG_RMSE(1,eye) = calculate.JAG_RMSE; %Stores the JAG1 RMSE values for each iteration in a separate table

end


%FINDING THE BEST FIT

%First we combined the Rsquared values for pSTAT3 and JAG into a single 
% matrix. Then we find the iteration that has the best fit. This is the one
% where the Rsquared point of pSTAT3 and JAG1 is closest to the point
% (1,1). The variable BestFit tells you the iteration number that has the 
% best fit to the in vitro data for pSTAT3 and JAG1.

RMSE_ALL = [pSTAT_RMSE ;JAG_RMSE]';
BestFit = dsearchn(RMSE_ALL,[0 0]);
BestFit_Data = [pSTAT_ALL(:,BestFit) JAG_ALL(:,BestFit)];


%NOW LET'S GRAPH OUR MONTE CARLO SIMULATION

TimeSpan(:,1:LagTimeMins)=[];
TimeSpan = TimeSpan - (LagTimeMins*60);

cell_1_MCFigureSTAT;
cell_1_MCFigureJAG;
%cell_1_MCFigureRMSEscatter;
%cell_1_MCFigureRSQUAREscatter;

toc