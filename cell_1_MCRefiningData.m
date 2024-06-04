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


%In this script we refine the data based on the RMSE score results to
%remove everything where JAG_RMSE and pSTAT_RMSE are in the bottom 10%.


%REMOVE TOP 99.5% OF RMSE VALUES IN BOTH PROTEINS

%Make binary matrices where data in the bottom 0.5% = 1
Refining.pSTAT_RMSE = pSTAT_RMSE;
Refining.JAG_RMSE = JAG_RMSE;

Refining.p_pSTAT = min(pSTAT_RMSE) + 0.005 * (max(pSTAT_RMSE) - min(pSTAT_RMSE));
Refining.pSTAT_RMSE( Refining.pSTAT_RMSE > Refining.p_pSTAT ) = 0;
Refining.pSTAT_RMSE(Refining.pSTAT_RMSE ~= 0) = 1;

Refining.p_JAG = min(JAG_RMSE) + 0.005 * (max(JAG_RMSE) - min(JAG_RMSE));
Refining.JAG_RMSE( Refining.JAG_RMSE > Refining.p_JAG ) = 0;
Refining.JAG_RMSE(Refining.JAG_RMSE ~= 0) = 1;

%Now we combine the binary matrices and remove the data in each group over
%in the top 99.5%
Refining.binary_RMSE = Refining.pSTAT_RMSE .* Refining.JAG_RMSE;
Refining.RMSEleftover = sum(Refining.binary_RMSE);

Refining.pSTAT_bottom = pSTAT_RMSE .* Refining.binary_RMSE;
Refining.JAG_bottom = JAG_RMSE .* Refining.binary_RMSE;
Refining.iterations_bottom = plots.keepingtrack .* Refining.binary_RMSE;


%PUT THE REFINED DATA TOGETHER

Refining.RMSE_all = [Refining.iterations_bottom; Refining.pSTAT_bottom; Refining.JAG_bottom];
Refining.RMSE_all( :, ~any(Refining.RMSE_all,1) ) = [];

Refining.pSTAT_ALL = ones(2881, n);
Refining.pSTAT_ALL = Refining.pSTAT_ALL(:,:) .* Refining.binary_RMSE;
Refining.pSTAT_ALL = Refining.pSTAT_ALL .* pSTAT_ALL;
Refining.pSTAT_ALL( :, ~any(Refining.pSTAT_ALL,1) ) = [];

Refining.JAG_ALL = ones(2881, n);
Refining.JAG_ALL = Refining.JAG_ALL(:,:) .* Refining.binary_RMSE;
Refining.JAG_ALL = Refining.JAG_ALL .* JAG_ALL;
Refining.JAG_ALL( :, ~any(Refining.JAG_ALL,1) ) = [];



