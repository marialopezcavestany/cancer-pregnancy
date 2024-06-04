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

% Create figure and calculate time in hours
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
thours = TimeSpan/3600;

% Create axes
axes1 = axes('Parent', figure1);
hold(axes1,'on');

%Choose the colors that all of the lines are going to be shown as
colororder([0 0.447058823529412 0.741176470588235;0.466666666666667 0.674509803921569 0.188235294117647;0.929411764705882 0.694117647058824 0.125490196078431;0.650980392156863 0.650980392156863 0.650980392156863]);

%Use plot for all the lines from the computaitonal model, but use scatter
%for all the individual points from the In Vitro data
plot(thours,JAG_ten,'Parent',axes1);
%plot(thours,Refining.JAG_ALL,'Parent',axes1);
scatter(plots.JAG_InVitro_t, plots.JAG_InVitro,'MarkerEdgeColor',[0 0 0],'LineWidth',2.5);

% Create ylabel and xlabel
ylabel('JAG1 Expression (nM)','HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontName','Arial');
xlabel('Time (min)','HorizontalAlignment','center','FontWeight','bold',...
    'FontName','Arial');

%Set up the axes limits
ylim(axes1,[0 4]);
xlim(axes1,[0 48]);
hold(axes1,'off');

% Set the remaining axes properties
set(axes1,'FontName','Arial','FontSize',12,'LineWidth',1.5,'TickDir','out','XTick',...
    [0 12 24 36 48]);