
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

close all
clear

% This script sets up the main model computational elements to run a single
% time with the input initial protein values and kinetic constants. This
% code does NOT run the monte carlo simulation or the sensitivity analysis.

% The initial conditions and kinetic constants reflect our fitted values
% from the Monte Carlo Simulation

% To sucessfully run this model you need the following scripts open:
%    > cell_1_mechanism.m


%DEFINING INITIAL CONCENTRATIONS

%Set non-zero initial values
PRL = 123; %ng/mL
PRL = PRL*0.04545; %Calculate the concentration in nM
RJ = 5.0802; 
S3c = 117.1798;
SHP2 = 78.23736788;
PPX = 70.3604;
PPN = 143.870996;
pS3c = 0;
S3n = 101.8617; 
JAG = 0.1528;

%Set zero values and fill in non-zero initial values
y0 = zeros(1,35);
y0(1) = PRL;
y0(2) = RJ;
y0(3) = S3c;
y0(4) = SHP2;
y0(5) = PPX;
y0(6) = PPN;
y0(11) = pS3c;
y0(20) = S3n;
y0(34) = JAG;


%DEFINING THE TIME SPAN

TimeSpan = [0:60:48*3600];


%RUN ODE SOLVER

options = odeset('RelTol',1e-7,'AbsTol',1e-7,'NonNegative',[1:length(y0)]);
[t,y] = ode15s(@cell_1_mechanism, TimeSpan, y0, options);


%CALCULATE QUANTATIES

total_c_pSTAT = y(:,11);  % Total pSTAT in the cytoplasm
total_c_pSTAT_norm = total_c_pSTAT ./ total_c_pSTAT(16); % Normalized to 15 minutes
JAG1_norm = y(:,34)./y(573,34); %JAG1c normalized to 30 minutes

Time_mins = (t(:,:)/60);
Time_hrs = (t(:,:)/3600);


%CALCULATE THE RMSE
plots.pSTAT_InVitro = [1.47 1.14 0.99 1.06];
plots.pSTAT_InVitro_t = [5 15 30 60];
plots.JAG_InVitro = [1.93 2.01 1.8];
plots.JAG_InVitro_t = [12 24 48];

calculate.pSTAT_Model = [y(6,11) y(16,11) y(31,11) y(61,11)];
calculate.JAG_Model = [y(721,34) y(1441,34) y(2881,34)];

calculate.pSTAT_RMSE = sqrt(mean((calculate.pSTAT_Model - plots.pSTAT_InVitro).^2));
calculate.JAG_RMSE = sqrt(mean((calculate.JAG_Model - plots.JAG_InVitro).^2));
