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

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
histogram(y0_stored(20,:),'Parent',axes1,'LineWidth',1, 'FaceColor',[0.498039215686275 0.6 0.788235294117647],'BinMethod','auto');
ylabel('Count','FontWeight','bold','FontSize',14);
xlabel(['S3n Variation'],'FontWeight','bold','FontSize',14);
hold(axes1,'off');
set(axes1,'FontName','Arial','FontSize',12,'LineWidth',1,'TickDir','out');