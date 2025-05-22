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

function dydt = cell_1_mechanism(t,y)

%DEFINING VARIABLE NAMES
PRL = y(1);                RJ = y(2); 
S3c = y(3);                SHP2 = y(4); 
PPX = y(5);                PPN = y(6); 
PRLRJ = y(7);              PRLRJ2 = y(8); 
PRLRJ2a = y(9);            PRLRJ2aS3c = y(10); 
pS3c = y(11);              pS3cpS3c = y(12);
PRLRJ2aSHP2 = y(13);       PPXpS3c = y(14); 
PPXpS3cpS3c = y(15);       pS3cS3c = y(16);
pS3npS3n = y(17);          pS3n = y(18);
PPNpS3n = y(19);           S3n = y(20); 
PPNpS3npS3n = y(21);       pS3nS3n = y(22); 
mRNASOCS1n = y(23);        mRNASOCS1c = y(24);
SOCS1 = y(25);             SOCS1PRLRJ2a = y(26); 
PRLRJ2aS3cSHP2 = y(27);    SOCS1PRLRJ2aSHP2 = y(28); 
mRNARJn = y(29);           mRNARJc = y(30);
RJc = y(31);               mRNAJAG1n = y(32); 
mRNAJAG1c = y(33);         JAG1c = y(34);
JAG1 = y(35);


%DEFINING THE RATE CONSTANTS

%Reaction 1: -> RJ
k1 = 3.08E-03;

%Reaction 2: PRL + RJ <-> PRLRJ
k2f = 7.7971E-04; k2r = 0.005;

%Reaction 3: 2 PRLRJ <-> PRLRJ2
k3f = 0.063415984; k3r = 0.2;

%Reaction 4: PRLRJ2 -> PRLRJ2a
k4 = 0.0052;

%Reaction 5: PRLRJ2a + S3c <-> PRLRJ2aS3c
k5f = 0.0124; k5r = 8.00E-01;

%Reaction 6: PRLRJ2aS3c -> pS3c + PRLRJ2a
k6 = 0.4417;

%Reaction 7: RJ ->, PRLRJ ->, PRLR2 ->, PRLJ2a ->
kdeg = 2.57E-04; deg_ratio = 7.1442;

%Reaction 8: pS3c + pS3c <-> pS3pS3c
k8f = 0.1814; k8r = 1.00E-01;

%Reaction 9: PRLRJ2a + SHP <-> PRLRJ2aSHP
k9f = 0.000879762; k9r = 2.00E-01;

%Reaction 10: PRLRJ2aSHP -> PRLRJ2 + SHP
k10 = 3.00E-03;

%Reaction 11: PPX + pS3cpS3c <-> PPXpS3cpS3c, PPX + pS3c <-> PPXpS3c
k11f = 0.0028; k11r = 2.00E-01;

%Reaction 12: PPXpS3cpS3c -> PPX + pS3c + S3c, PPXpS3c -> PPX + S3c
k12 = 0.0031;

%Reaction 13: pS3c + S3c <-> pS3cS3c, pS3n + S3n <-> pS3nS3n
k13f = 2.00E-07; k13r = 0.225430866;

%Reaction 14: pS3pS3c -> pS3pS3n
k14 = 0.005886461;

%Reaction 15: PPN + pS3npS3n <-> PPNpS3npS3n, PPN + pS3n <-> PPNS3n
k15f = 0.001460885; k15r = 2.00E-01;

%Reaction 16: PPNpS3npS3n -> PPN + pS3n + S3n, PPNS3n -> PPN + S3n
k16 = 0.0145287;

%Reaction 17: S3c <-> S3n
k17f = 3.55E-02; k17r = 0.062303131;

%Reaction 18: pS33n -> mRNASOCS1n + pS33n (2-step reaction)
k18_1 = 1.00E-02; k18_2 = 4.00E+02;

%Reaction 19: mRNASOCS1n -> mRNASOCSc
k19 = 0.001603988;

%Reaction 20: mRNASOCS1c -> SOCS1 + mRNASOCS1c
k20 = 1.00E-02;

%Reaction 21: SOCS1 + PRLRJ2a <-> SOCS1PRLRJ2a
k21f = 0.052168498; k21r = 1.00E-01;

%Reaction 22: mRNASOCS1c ->
k22 = 0.000506751;

%Reaction 23: SOCS1 ->
k23 = 6.8485E-04;

%Reaction 24: SOCS1PRLRJ2a ->
k24 = 0.00124969;

%Reaction 25: pS33n -> mRNARJn + pS33n (2-step reaction)
k25_1 = 0.0056; k25_2 = 4.00E+02;

%Reaction 26: mRNARJn -> mRNARJc
k26 = 1.00E-03;

%Reaction 27: mRNARJc ->
k27 = 6.12582E-05;

%Reaction 28: mRNARJc -> RJc + mRNARJc
k28 = 0.0132;

%Reaction 29: RJc -> RJ
k29 = 1.00E-02;

%Reaction 30: pS3npS3n -> mRNAJAG1n + pS3npS3n (2-step reaction)
k30_1 = 0.0021; k30_2 = 4.00E+02;

%Reaction 31: mRNAJAG1n -> mRNAJAG1c
k31 = 0.0011;

%Reaction 32: mRNAJAG1c ->
k32 = 5.00E-04;

%Reaction 33: mRNAJAG1c -> JAG1c + mRNAJAG1c
k33 = 0.0099;

%Reaction 34: JAG1c -> JAG1
k34 = 1.4518E-04;

%Other parameters
Vratio = 5.00E-01;


%DEFINIING THE MASS BALANCES
 %PRL = -R2
 dydt(1,1) = (k2r*PRLRJ - k2f*PRL*RJ) * (1.39*10^ - 4);

 %RJ = R1 - R2 - R3 - R7 + R29
 dydt(2,1) = k1 + (k2r*PRLRJ - k2f*PRL*RJ) + (k3r*PRLRJ2 - k3f*PRLRJ*RJ) - kdeg*RJ + k29*RJc;

 %S3c = -R5 + (R7 x 2) + R10 + R12 - R13 - R17
 dydt(3,1) = (k5r*PRLRJ2aS3c - k5f*PRLRJ2a*S3c) + kdeg*deg_ratio*(PRLRJ2aS3c + PRLRJ2aS3cSHP2) + k10*PRLRJ2aS3cSHP2 + k12*PPXpS3c + (k13r*pS3cS3c - k13f*pS3c*S3c) + (k17r*S3n.*Vratio - k17f*S3c);
 
 %SHP2 = (R7 x3) + (-R9 x 3) + (R10 x 3) + R24
 dydt(4,1) =  kdeg*deg_ratio*(PRLRJ2aSHP2 + PRLRJ2aS3cSHP2 + SOCS1PRLRJ2aSHP2) + (k9r*PRLRJ2aSHP2 - k9f*PRLRJ2a*SHP2) + (k9r*SOCS1PRLRJ2aSHP2 - k9f*SOCS1PRLRJ2a*SHP2) + (k9r*PRLRJ2aS3cSHP2 - k9f*PRLRJ2aS3c*SHP2) + k10*PRLRJ2aS3cSHP2 + k10*PRLRJ2aSHP2 + k10*SOCS1PRLRJ2aSHP2 + k24*SOCS1PRLRJ2aSHP2;
 
 %PPX = (-R11 x 2) + (R12 x 2)
 dydt(5,1) = (k11r*PPXpS3c - k11f*PPX*pS3c) + (k11r*PPXpS3cpS3c - k11f*PPX*pS3cpS3c) + k12*PPXpS3c + k12*PPXpS3cpS3c;
 
 %PPN = (-R15 x 2) + (R16 x 2)
 dydt(6,1) = (k15r*PPNpS3n - k15f*PPN*pS3n) + (k15r*PPNpS3npS3n - k15f*PPN*pS3npS3n) + k16*PPNpS3n + k16*PPNpS3npS3n;
 
 %PRLRJ = R2 - R3 - R7
 dydt(7,1) = (k2f*PRL*RJ - k2r*PRLRJ) + (k3r*PRLRJ2 - k3f*PRLRJ*RJ) - kdeg*deg_ratio*PRLRJ;
 
 %PRLRJ2 = R3 - R4 - R7 + (R10 x 3)
 dydt(8,1) = (k3f*PRLRJ*RJ - k3r*PRLRJ2) - k4*PRLRJ2 - kdeg*deg_ratio*PRLRJ2 + k10*PRLRJ2aSHP2 + k10*SOCS1PRLRJ2aSHP2 + k10*PRLRJ2aS3cSHP2;
 
 %PRLRJ2a = R4 - R5 + R6 - R7 - R9 - R21 + R24
 dydt(9,1) = k4*PRLRJ2 + (k5r*PRLRJ2aS3c - k5f*PRLRJ2a*S3c) + k6*PRLRJ2aS3c - kdeg*deg_ratio*PRLRJ2a + (k9r*PRLRJ2aS3cSHP2 - k9f*PRLRJ2a*SHP2) + (k21r*SOCS1PRLRJ2a - k21f*SOCS1*PRLRJ2a) + k24*SOCS1PRLRJ2a;

 %PRLRJ2aS3c = R5 - R6 -R7 - R9
 dydt(10,1) = (k5f*PRLRJ2a*S3c - k5r*PRLRJ2aS3c) - k6*PRLRJ2aS3c - kdeg*deg_ratio*PRLRJ2aS3c + (k9r*PRLRJ2aS3cSHP2 - k9f*PRLRJ2aS3c*SHP2);

 %pS3c = R6 - 2*R8 - R11 - R13
 dydt(11,1) = k6*PRLRJ2aS3c + 2*(k8r*pS3cpS3c - k8f*pS3c*pS3c) + (k11r*PPXpS3c - k11f*PPX*pS3c) + (k13r*pS3cS3c - k13f*pS3c*S3c);

 %pS3cpS3c = R8 - R11 - R14
 dydt(12,1) = (k8f*pS3c*pS3c - k8r*pS3cpS3c) + (k11r*PPXpS3cpS3c - k11f*PPX*pS3cpS3c) - k14*pS3cpS3c;

 %PRLRJ2aSHP2 = -R7 + R9 - R10 + R23
 dydt(13,1) =  -kdeg*deg_ratio*PRLRJ2aSHP2 + (k9f*PRLRJ2a*SHP2 - k9r*PRLRJ2aSHP2) - k10*PRLRJ2aSHP2 + k23*SOCS1PRLRJ2aSHP2;

 %PPXpS3c = R11 - R12
 dydt(14,1) = (k11f*PPX*pS3c - k11r*PPXpS3c) - k12*PPXpS3c;

 %PPXpS3cpS3c = R11 - R12
 dydt(15,1) = (k11f*PPX*pS3cpS3c - k11r*PPXpS3cpS3c) - k12*PPXpS3cpS3c;

 %pS3cS3c = R12 + R13
 dydt(16,1) = k12*PPXpS3cpS3c + (k13f*pS3c*S3c - k13r*pS3cS3c);
 
 %pS3npS3n = R8 + R14 - R15 - R30 - R18 - R25
 dydt(17,1) = (k8f*pS3n*pS3n - k8r*pS3npS3n) + k14*pS3cpS3c./Vratio + (k15r*PPNpS3npS3n - k15f*PPN*pS3npS3n);

 %pS3n = -2*R8 - R13 - R15
 dydt(18,1) = 2*(k8r*pS3npS3n - k8f*pS3n*pS3n) + (k13r*pS3nS3n - k13f*pS3n*S3n) + (k15r*PPNpS3n - k15f*PPN*pS3n);

 %PPNpS3n = R15 - R16
 dydt(19,1) = (k15f*PPN*pS3n - k15r*PPNpS3n) - k16*PPNpS3n;

 %S3n = -R13 + R16 + R17
 dydt(20,1) = (k13r*pS3nS3n - k13f*pS3n*S3n) + k16*PPNpS3n + (k17f*S3c./Vratio - k17r*S3n);

 %PPNpS3npS3n = R15 - R16
 dydt(21,1) = (k15f*PPN*pS3npS3n - k15r*PPNpS3npS3n) - k16*PPNpS3npS3n;

 %pS3nS3n = R13 + R16
 dydt(22,1) = (k13f*pS3n*S3n - k13r*pS3nS3n) + k16*PPNpS3npS3n;

 %mRNASOCS1n = R18 - R19
 dydt(23,1) = (k18_1*pS3npS3n)./(k18_2+pS3npS3n) - k19*mRNASOCS1n;

 %mRNASOCS1c = R19 - R22
 dydt(24,1) = k19*mRNASOCS1n.*Vratio - k22*mRNASOCS1c;

 %SOCS1 = (R7 x 2) + R10 + R20 - R21 - R23
 dydt(25,1) = kdeg*deg_ratio*(SOCS1PRLRJ2a + SOCS1PRLRJ2aSHP2) + k10*SOCS1PRLRJ2aSHP2 + k20*mRNASOCS1c + (k21r*SOCS1PRLRJ2a - k21f*SOCS1*PRLRJ2a) - k23*SOCS1;

 %SOCS1PRLRJ2a = -R7 + R21 - R9 - R23 - R24
 dydt(26,1) =  -kdeg*deg_ratio*SOCS1PRLRJ2a + (k9r*SOCS1PRLRJ2aSHP2 - k9f*SOCS1PRLRJ2a*SHP2) + (k21f*SOCS1*PRLRJ2a - k21r*SOCS1PRLRJ2a) - k23*SOCS1PRLRJ2a - k24*SOCS1PRLRJ2a;

 %PRLRJ2aS3cSHP2 = -R7 + R9 - R10
 dydt(27,1) = -kdeg*deg_ratio*PRLRJ2aS3cSHP2 + (k9f*PRLRJ2aS3c*SHP2 - k9r*PRLRJ2aS3cSHP2) - k10*PRLRJ2aS3cSHP2;

 %SOCS1PRLRJ2aSHP2 =  - R7 + R9 - R10 - R23 - R24
 dydt(28,1) = -kdeg*deg_ratio*SOCS1PRLRJ2aSHP2 + (k9f*SOCS1PRLRJ2a*SHP2 - k9r*SOCS1PRLRJ2aSHP2) - k10*SOCS1PRLRJ2aSHP2 - k23*SOCS1PRLRJ2aSHP2 - k24*SOCS1PRLRJ2aSHP2;

 %mRNARJn = R25 - R26
 dydt(29,1) = (k25_1*pS3npS3n)./(k25_2+pS3npS3n) - k26*mRNARJn;

 %mRNARJc = R26 - R27
 dydt(30,1) = k26*mRNARJn.*Vratio - k27*mRNARJc;

 %RJc = R28 - R29
 dydt(31,1) = k28*mRNARJc - k29*RJc;

 %mRNAJAG1n = R30 - R31
 dydt(32,1) = (k30_1*pS3npS3n)./(k30_2+pS3npS3n) - k31*mRNAJAG1n;

 %mRNAJAG1c = R31 - R32
 dydt(33,1) = k31*mRNAJAG1n.*Vratio - k32*mRNAJAG1c;

 %JAG1c = R33 - R34
 dydt(34,1) = k33*mRNAJAG1c - k34*JAG1c;
 
 %JAG1 = R34
 dydt(35) = k34*JAG1c;

end
