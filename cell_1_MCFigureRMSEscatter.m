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


%This script was developed based on an automatically generated MATLAB
%script for a customized figure

% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create scatter for ALL RMSE
scatter(pSTAT_RMSE, JAG_RMSE,...
    'MarkerEdgeColor',[0.584313725490196 0.717647058823529 0.811764705882353],...
    'LineWidth',1);

% Create scatter for CLEANED RMSE
%scatter(pSTAT_RMSE, JAG_RMSE,...
    %'MarkerEdgeColor',[0.584313725490196 0.717647058823529 0.811764705882353],...
    %'LineWidth',1);

% Create x and y axis labels
ylabel('JAG1 RMSE','FontWeight','bold','FontSize',14,'FontName','Arial');
xlabel('pSTAT3 RMSE','FontWeight','bold','FontSize',14,'FontName','Arial');

%Set the axis properties
xlim(axes1,'off');
hold(axes1,'off');
%set(axes1,'FontName','Arial','FontSize',12,'LineWidth',1.5,'TickDir','out',...
    %'XTick',[0 1 2 3 4],'YTick',[0 20 40 60 80 100 120 140]);
%set(axes1,'FontName','Arial','FontSize',12,'LineWidth',1.5,'TickDir','out');

%OR! If you need a log y-axis comment the properties above and uncomment
%the ones below
set(axes1,'FontName','Arial','FontSize',12,'LineWidth',1.5,'TickDir','out',...
    'XTick',[0 1 2 3 4],'YMinorTick','on','YScale','log','YTick',...
    [0.1 1 10 100 1000]);