% ALL METABOLIC PARAMETERS: Main text Fig. 4, Fig. 5, & Supplementary material Fig. S5
% 2016 Ashley E. Maloney (ashjames@uw.edu)
% Stationary 6-box model of hydrogen and deuterium cycling in phytoplankton 
% For manuscript "Exploring lipid ²H/¹H fractionation mechanisms in 
% response to salinity with continuous cultures of the diatom Thalassiosira
% pseudonana" published in Organic Geochemistry by A.E. Maloney, A.L.C. Shinneman, 
% K. Hemeon, and J.P. Sachs of the University of Washington. http://dx.doi.org/10.1016/j.orggeochem.2016.08.015
% Tests sensitivity of dDlipid & other pools to metabolic parameters 
% (D = ²H = deuterium)
clear all; clc; close all;
%% FRACTIONATION FACTORS
A1=1;       %media-H to cell-water-H
A2=1;       %cell-water-H to media-H 
A3=1;       %cell-water-H to other-organics-H (equilibrium exchange)
A4=1;       %other-organic-H to cell-water-H (equilibrium exchange)
A5=1;       %cell-water-H to other-organics-H
A6=1;       %other-organics-H to cell-water-H
A7=0.4;     %cell-water-H to PSI-NADPH (ferredoxin reduces NADP+ to NADPH)
A8=0.75;    %other-orgnaic-H to MET-NAD(P)H (G6PDH reduces NADP+ to NADPH)
A9=1;       %cell-water-H to Lipid-H
A10=1;      %PS1-NADPH to other-organic-H
A11=1;      %PS1-NADPH to Lipid-H
A12=1;      %MET-NAD(P)H to other-organic-H 
A13=1;      %MET-NAD(P)H to Lipid-H
A14=1;      %Lipid-H to other-organic-H
A15=1;      %other-organic-H to Lipid-H
A16=1;      %other-organic-H to media-H
A17=1;      %Lipid-H to cell-water-H
A18=1;      %Lipid-H to media-H
A19=1;      %Lipid-H to MET-NAD(P)H
%% Solution for Fig.4 & Fig.S5 panel a) PSI vs MET PARAMETER
%<<<< TESTABLE METABOLIC PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sachsX=0.001;  %"PS1 vs MET NAD(P)H Parameter" ****testing this parameter in loop below****            
M=1;           %"Metabolic Reductant Parameter"
B=1;           %"Biosynthetic Cell Water Parameter"
E=1;           %"Equilibrium Exchange Parameter"
HT=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PER=0.03;      %"Percent Extracellular Release"
%additional parameters not tested in manuscript:
L_PER=0.2;     %"Lipid Exudate Parameter"
R=0.2;         %"B-oxidation/Recyle Parameter" 
               %*R can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7=1;                            %JWAT_PSI: cell-water-H to PSI-NADPH
j10=j7*(1-sachsX);               %JPSI_Org: PS1-NADPH to other-organic-H
j11=j7*sachsX;                   %JPS1_Lip: PS1-NADPH to Lipid-H
j8=j7*M;                         %JORG_MET: other-orgnaic-H to MET-NAD(P)H 
j13=j8*(1-sachsX);               %JMET_Lip: MET-NAD(P)H to Lipid-H
j9=(j13+j11)/2;                  %JBIO_Wat_Lip: cell-water-H to Lipid-H
j15=(j13+j11)/2;                 %JBIO_Org_Lip: other-organic-H to Lipid-H
j1=j7*HT;                        %Jin: media-H to cell-water-H
j16=j7*PER*(1-L_PER);            %Jrelease_Org: other-organic-H to media-H
j18=j7*PER*L_PER;                %Jrelease_Lip: Lipid-H to cell-water-H
j3=j7*E;                         %Jeq_Org_Wat: cell-water-H to other-organics-H
j4=j3;                           %Jeq_Wat_Org: other-organic-H to cell-water-H
j14=(j15+j13+j11+j9-j18)*(1-2*R);%JCAT_Lip_Org: Lipid-H to other-organic-H (Mass Balance)
j17=(j15+j13+j11+j9-j18)*R;      %JCAT_Lip_Wat: Lipid-H to cell-water-H (Mass Balance)
j19=(j15+j13+j11+j9-j18)*R;      %JCAT_Lip_Met: Lipid-H to MET-NAD(P)H (Mass Balance)
j12=(j8*sachsX)+j19;             %JMET_Org: MET-NAD(P)H to other-organic-H 
j5=(j10+j12)*B;                  %JBIO_Wat_Org: cell-water-H to other-organics-H
j6=j12+j14+j5+j10-j8-j15-j16;    %JCAT_Org_Wat: other-organics-H to cell-water-H (Mass Balance)
j2=j6+j17+j1-j7-j5-j9;           %Jout: cell-water-H to media-H (Mass Balance)
%<<<< SOLVE FOR R (D/H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MAT=[(j2*A2),(j16*A16),(j18*A18),0,0;
    -(j3*A3+j5*A5+j7*A7+j2*A2+j9*A9),(j4*A4+j6*A6),(j17*A17),0,0;
     (j3*A3+j5*A5),-(j4*A4+j8*A8+j6*A6+j16*A16+j15*A15),(j14*A14),(j10*A10),(j12*A12);
     (j7*A7),0,0,-(j11*A11+j10*A10),0;
     0,(j8*A8),(j19*A19),0,-(j12*A12+j13*A13);
     (j9*A9),(j15*A15),-(j17*A17+j14*A14+j18*A18+(j19*A19)),(j11*A11),(j13*A13)];
VEC=[((j1*A1*RM)); -(j1*A1*RM); 0; 0; 0; 0];
X=MAT\VEC;
RW=X(1);RO=X(2);RL=X(3);RNP=X(4);RNm=X(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDM=(RM-1)*1000;dDW=(RW-1)*1000;dDO=(RO-1)*1000;dDNP=(RNP-1)*1000;dDNm=(RNm-1)*1000;dDL=(RL-1)*1000;
%<<<< CALCULATE APPARENT FRACTIONATION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaL=(dDL+1000)/(dDM+1000);
alphaO=(dDO+1000)/(dDM+1000);
%<<<< LOOP THROUGH; GRAPH TESTABLE VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.1;
for i = 2:11
sachsX(i)=sachsX(i-1)+add;
j10(i)=j7*(1-sachsX(i));         
j11(i)=j7*sachsX(i);
j13(i)=j8*(1-sachsX(i));
j15(i)=(j13(i)+j11(i))/2;
j9(i)=(j13(i)+j11(i))/2;
j14(i)=(j15(i)+j13(i)+j11(i)+j9(i)-j18)*(1-R-R);
j17(i)=(j15(i)+j13(i)+j11(i)+j9(i)-j18)*R;
j19(i)=(j15(i)+j13(i)+j11(i)+j9(i)-j18)*R;
j12(i)=(j8*sachsX(i))+j19(i);
j5(i)=B*(j10(i)+j12(i));
j6(i)=j12(i)+j14(i)+j5(i)+j10(i)-j8-j15(i)-j16;
j2(i)=j6(i)+j17(i)+j1-j7-j5(i)-j9(i);                  
%D/H ratio
MAT2=[(j2(i)*A2),(j16*A16),(j18*A18),0,0;
     -(j3*A3+j5(i)*A5+j7*A7+j2(i)*A2+j9(i)*A9),(j4*A4+j6(i)*A6),(j17(i)*A17),0,0;
      (j3*A3+j5(i)*A5),-(j4*A4+j8*A8+j6(i)*A6+j16*A16+j15(i)*A15),(j14(i)*A14),(j10(i)*A10),(j12(i)*A12);
      (j7*A7),0,0,-(j11(i)*A11+j10(i)*A10),0;
      0,(j8*A8),(j19(i)*A19),0,-(j12(i)*A12+j13(i)*A13);
      (j9(i)*A9),(j15(i)*A15),-(j17(i)*A17+j14(i)*A14+j18*A18+(j19(i)*A19)),(j11(i)*A11),(j13(i)*A13)];
VEC2=[((j1*A1*RM)); -(j1*A1*RM); 0; 0; 0; 0];
X2=MAT2\VEC2;
RW(i)=X2(1);RO(i)=X2(2);RL(i)=X2(3);RNP(i)=X2(4);RNm(i)=X2(5);
%dD
dDM(i)=(RM-1)*1000;dDW(i)=(RW(i)-1)*1000;dDO(i)=(RO(i)-1)*1000;dDNP(i)=(RNP(i)-1)*1000;dDNm(i)=(RNm(i)-1)*1000;dDL(i)=(RL(i)-1)*1000;
%alpha
alphaL(i)=(dDL(i)+1000)/(dDM(i)+1000);
alphaO(i)=(dDO(i)+1000)/(dDM(i)+1000);
end
%% Solution for Fig.4 & Fig.S5 panel b) METABOLIC REDUCTANT PARAMETER
%<<<< TESTABLE METABOLIC PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sachsXM=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
MM=.0001;       %"Metabolic Reductant Parameter" ****testing this paramter in loop below****
BM=1;           %"Biosynthetic Cell Water Parameter"
EM=1;           %"Equilibrium Exchange Parameter"
HTM=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PERM=0.03;      %"Percent Extracellular Release"
%additional parameters not tested in manuscript:
L_PERM=0.2;     %"Lipid Exudate Parameter"
ReM=0.2;        %"B-oxidation/Recyle Parameter" 
                %*Re can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7M=1;                                  %JWAT_PSI: cell-water-H to PSI-NADPH
j10M=j7M*(1-sachsXM);                   %JPSI_Org: PS1-NADPH to other-organic-H
j11M=j7M*sachsXM;                       %JPS1_Lip: PS1-NADPH to Lipid-H
j8M=j7M*MM;                             %JORG_MET: other-orgnaic-H to MET-NAD(P)H 
j13M=j8M*(1-sachsXM);                   %JMET_Lip: MET-NAD(P)H to Lipid-H
j9M=(j13M+j11M)/2;                      %JBIO_Wat_Lip: cell-water-H to Lipid-H
j15M=(j13M+j11M)/2;                     %JBIO_Org_Lip: other-organic-H to Lipid-H
j1M=j7M*HTM;                            %Jin: media-H to cell-water-H
j16M=(j7M*PERM)*(1-L_PERM);             %Jrelease_Org: other-organic-H to media-H        
j18M=(j7M*PERM)*L_PERM;                 %Jrelease_Lip: Lipid-H to cell-water-H
j3M=j7M*EM;                             %Jeq_Org_Wat: cell-water-H to other-organics-H
j4M=j3M;                                %Jeq_Wat_Org: other-organic-H to cell-water-H
j14M=(j15M+j13M+j11M+j9M-j18M)*(1-(2*ReM));%JCAT_Lip_Org: Lipid-H to other-organic-H (Mass Balance)     
j17M=(j15M+j13M+j11M+j9M-j18M)*ReM;     %JCAT_Lip_Wat: Lipid-H to cell-water-H (Mass Balance)    
j19M=(j15M+j13M+j11M+j9M-j18M)*ReM;     %JCAT_Lip_Met: Lipid-H to MET-NAD(P)H (Mass Balance)     
j12M=(j8M*sachsXM)+j19M;                %JMET_Org: MET-NAD(P)H to other-organic-H     
j5M=(j10M+j12M)*BM;                     %JBIO_Wat_Org: cell-water-H to other-organics-H
j6M=j12M+j14M+j5M+j10M-j8M-j15M-j16M;   %JCAT_Org_Wat: other-organics-H to cell-water-H (Mass Balance)
j2M=j6M+j17M+j1M-j7M-j5M-j9M;           %Jout: cell-water-H to media-H (Mass Balance)
%<<<< SOLVE FOR R (D/H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MATM=[(j2M*A2),(j16M*A16),(j18M*A18),0,0;
    -(j3M*A3+j5M*A5+j7M*A7+j2M*A2+j9M*A9),(j4M*A4+j6M*A6),(j17M*A17),0,0;
     (j3M*A3+j5M*A5),-(j4M*A4+j8M*A8+j6M*A6+j16M*A16+j15M*A15),(j14M*A14),(j10M*A10),(j12M*A12);
     (j7M*A7),0,0,-(j11M*A11+j10M*A10),0;
     0,(j8M*A8),(j19M*A19),0,-(j12M*A12+j13M*A13);
     (j9M*A9),(j15M*A15),-(j17M*A17+j14M*A14+j18M*A18+j19M*A19),(j11M*A11),(j13M*A13)];
VECM=[((j1M*A1*RM)); -(j1M*A1*RM); 0; 0; 0; 0];
XM=MATM\VECM;
RWM=XM(1);ROM=XM(2);RLM=XM(3);RNPM=XM(4);RNmM=XM(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDMM=(RM-1)*1000;dDWM=(RWM-1)*1000;dDOM=(ROM-1)*1000;dDNPM=(RNPM-1)*1000;dDNmM=(RNmM-1)*1000;dDLM=(RLM-1)*1000;
%<<<< CALCULATE APPARENT FRACTIONATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaLM=(dDLM+1000)/(dDMM+1000);
alphaOM=(dDOM+1000)/(dDMM+1000);
%<<<< LOOP THROUGH; GRAPH TESTABLE VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.1;
for i = 2:110
MM(i)=MM(i-1)+add;
j8M(i)=j7M*MM(i);
j13M(i)=j8M(i)*(1-sachsXM);
j15M(i)=(j13M(i)+j11M)/2;
j9M(i)=(j13M(i)+j11M)/2;
j14M(i)=(j15M(i)+j13M(i)+j11M+j9M(i)-j18M)*(1-ReM-ReM);
j17M(i)=(j15M(i)+j13M(i)+j11M+j9M(i)-j18M)*ReM;
j19M(i)=(j15M(i)+j13M(i)+j11M+j9M(i)-j18M)*ReM;         
j12M(i)=(j8M(i)*sachsXM)+j19M(i);
j5M(i)=BM*(j10M+j12M(i));
j6M(i)=j12M(i)+j14M(i)+j5M(i)+j10M-j8M(i)-j15M(i)-j16M;
j2M(i)=j6M(i)+j17M(i)+j1M-j7M-j5M(i)-j9M(i);                  
%D/H ratio
MAT2M=[(j2M(i)*A2),(j16M*A16),(j18M*A18),0,0;
     -(j3M*A3+j5M(i)*A5+j7M*A7+j2M(i)*A2+j9M(i)*A9),(j4M*A4+j6M(i)*A6),(j17M(i)*A17),0,0;
      (j3M*A3+j5M(i)*A5),-(j4M*A4+j8M(i)*A8+j6M(i)*A6+j16M*A16+j15M(i)*A15),(j14M(i)*A14),(j10M*A10),(j12M(i)*A12);
      (j7M*A7),0,0,-(j11M*A11+j10M*A10),0;
      0,(j8M(i)*A8),(j19M(i)*A19),0,-(j12M(i)*A12+j13M(i)*A13);
      (j9M(i)*A9),(j15M(i)*A15),-(j17M(i)*A17+j14M(i)*A14+(j18M*A18)+j19M(i)*A19),(j11M*A11),(j13M(i)*A13)];
VEC2M=[(j1M*A1*RM); -(j1M*A1*RM); 0; 0; 0; 0];
X2M=MAT2M\VEC2M;
RWM(i)=X2M(1);ROM(i)=X2M(2);RLM(i)=X2M(3);RNPM(i)=X2M(4);RNmM(i)=X2M(5);
%dD
dDMM(i)=(RM-1)*1000;dDWM(i)=(RWM(i)-1)*1000;dDOM(i)=(ROM(i)-1)*1000;dDNPM(i)=(RNPM(i)-1)*1000;dDNmM(i)=(RNmM(i)-1)*1000;dDLM(i)=(RLM(i)-1)*1000;
%alpha
alphaLM(i)=(dDLM(i)+1000)/(dDMM(i)+1000);
alphaOM(i)=(dDOM(i)+1000)/(dDMM(i)+1000);
end
%% Solution for Fig.4 & Fig.S5 panel c) BIOSYNTHETIC CELL WATER PARAMTER - same result as EQUILIBRIUM EXCHANGE PARAMATER
%<<<< TESTABLE METABOLIC PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sachsXb=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
Mb=1;           %"Metabolic Reductant Parameter"    
Bb=.00001;      %"Biosynthetic Cell Water Parameter" ****testing this paramter in loop below****
Eb=1;           %"Equilibrium Exchange Parameter" ****testung this gives same result - see individual files for this****
HTb=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PERb=0.03;      %"Percent Extracellular Release"
%additional parameters not tested in manuscript:
L_PERb=0.2;     %"Lipid Exudate Parameter"
Reb=0.2;        %"B-oxidation/Recyle Parameter" 
                %*Re can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7b=1;                          
j10b=j7b*(1-sachsXb);                 
j11b=j7b*sachsXb;                     
j8b=j7b*Mb;                       
j13b=j8b*(1-sachsXb);              
j9b=(j13b+j11b)/2;               
j15b=(j13b+j11b)/2;              
j1b=j7b*HTb;                       
j16b=j7b*PERb*(1-L_PERb);                    
j18b=j7b*PERb*L_PERb;               
j3b=j7b*Eb;                           
j4b=j3b;                          
j14b=(j15b+j13b+j11b+j9b-j18b)*(1-2*Reb);     
j17b=(j15b+j13b+j11b+j9b-j18b)*Reb;                 
j21b=(j15b+j13b+j11b+j9b-j18b)*Reb; 
j12b=(j8b*sachsXb)+j21b;                     
j5b=(j10b+j12b)*Bb;                 
j6b=j12b+j14b+j5b+j10b-j8b-j15b-j16b;   
j2b=j6b+j17b+j1b-j7b-j5b-j9b;                  
%<<<< SOLVE FOR R (D/H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MATb=[(j2b*A2),(j16b*A16),(j18b*A18),0,0;
    -(j3b*A3+j5b*A5+j7b*A7+j2b*A2+j9b*A9),(j4b*A4+j6b*A6),(j17b*A17),0,0;
     (j3b*A3+j5b*A5),-(j4b*A4+j8b*A8+j6b*A6+j16b*A16+j15b*A15),(j14b*A14),(j10b*A10),(j12b*A12);
     (j7b*A7),0,0,-(j11b*A11+j10b*A10),0;
     0,(j8b*A8),(j21b*A19),0,-(j12b*A12+j13b*A13);
     (j9b*A9),(j15b*A15),-(j17b*A17+j14b*A14+j18b*A18+j21b*A19),(j11b*A11),(j13b*A13)];
VECb=[(j1b*A1*RM); -(j1b*A1*RM); 0; 0; 0; 0];
Xb=MATb\VECb;
RWb=Xb(1);ROb=Xb(2);RLb=Xb(3);RNPb=Xb(4);RNmb=Xb(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDMb=(RM-1)*1000;dDWb=(RWb-1)*1000;dDOb=(ROb-1)*1000;dDNPb=(RNPb-1)*1000;dDNmb=(RNmb-1)*1000;dDLb=(RLb-1)*1000;
%<<<< CALCULATE APPARENT FRACTIONATIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaLb=(dDLb+1000)/(dDMb+1000);
alphaOb=(dDOb+1000)/(dDMb+1000);
%<<<< LOOP THROUGH; GRAPH TESTABLE VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.1;
for i = 2:110
Bb(i)=Bb(i-1)+add;
j5b(i)=Bb(i)*(j10b+j12b); 
j6b(i)=j12b+j14b+j5b(i)+j10b-j8b-j15b-j16b;
j2b(i)=j6b(i)+j17b+j1b-j7b-j5b(i)-j9b; 
j13b(i)=j8b*(1-sachsXb);
j11b(i)=j7b*sachsXb;
%D/H ratio
MAT2b=[(j2b(i)*A2),(j16b*A16),(j18b*A18),0,0;
     -(j3b*A3+j5b(i)*A5+j7b*A7+j2b(i)*A2+j9b*A9),(j4b*A4+j6b(i)*A6),(j17b*A17),0,0;
      (j3b*A3+j5b(i)*A5),-(j4b*A4+j8b*A8+j6b(i)*A6+j16b*A16+j15b*A15),(j14b*A14),(j10b*A10),(j12b*A12);
      (j7b*A7),0,0,-(j11b(i)*A11+j10b*A10),0;
      0,(j8b*A8),j21b*A19,0,-(j12b*A12+j13b(i)*A13);
      (j9b*A9),(j15b*A15),-(j17b*A17+j14b*A14+j18b*A18+j21b*A19),(j11b(i)*A11),(j13b(i)*A13)];
VEC2b=[(j1b*A1*RM); -(j1b*A1*RM); 0; 0; 0; 0];
X2b=MAT2b\VEC2b;
RWb(i)=X2b(1);ROb(i)=X2b(2);RLb(i)=X2b(3);RNPb(i)=X2b(4);RNmb(i)=X2b(5);
%dD
dDMb(i)=(RM-1)*1000;dDWb(i)=(RWb(i)-1)*1000;dDOb(i)=(ROb(i)-1)*1000;dDNPb(i)=(RNPb(i)-1)*1000;dDNmb(i)=(RNmb(i)-1)*1000;dDLb(i)=(RLb(i)-1)*1000;
%alpha
alphaLb(i)=(dDLb(i)+1000)/(dDMb(i)+1000);
alphaOb(i)=(dDOb(i)+1000)/(dDMb(i)+1000);
end
%% Solution for Fig.4 & Fig.S5 panel d) and Fig. 5 panel a) HYDROGEN TRANSPORT PARAMETER
%<<<< TESTABLE METABOLIC PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sachsXt=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
Mt=1;           %"Metabolic Reductant Parameter"    
Bt=1;           %"Biosynthetic Cell Water Parameter" 
Et=1;           %"Equilibrium Exchange Parameter"
HTt=.03;        %"Hydrogen Transfer Parameter" *can't be smaller than "PER"! ****testing this paramter in loop below****
PERt=0.03;      %"Percent Extracellular Release"
%additional parameters not tested in manuscript:
L_PERt=0.2;     %"Lipid Exudate Parameter"
Ret=0.2;        %"B-oxidation/Recyle Parameter" 
                %*Re can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7t=1;                           
j10t=j7t*(1-sachsXt);                 
j11t=j7t*sachsXt;                     
j8t=j7t*Mt;                       
j13t=j8t*(1-sachsXt);                 
j9t=(j13t+j11t)/2;                 
j15t=(j13t+j11t)/2;                
j1t=j7t*HTt;                       
j16t=j7t*PERt*(1-L_PERt);                     
j18t=j7t*PERt*L_PERt;                
j3t=j7t*Et;                           
j4t=j3t;                          
j14t=(j15t+j13t+j11t+j9t-j18t)*(1-Ret-Ret);     
j17t=(j15t+j13t+j11t+j9t-j18t)*Ret;                 
j19t=(j15t+j13t+j11t+j9t-j18t)*Ret;
j12t=(j8t*sachsXt)+j19t;                     
j5t=(j10t+j12t)*Bt;                 
j6t=j12t+j14t+j5t+j10t-j8t-j15t-j16t;   
j2t=j6t+j17t+j1t-j7t-j5t-j9t;                  
%<<<< SOLVE FOR R (D/H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MATt=[(j2t*A2),(j16t*A16),(j18t*A18),0,0;
    -(j3t*A3+j5t*A5+j7t*A7+j2t*A2+j9t*A9),(j4t*A4+j6t*A6),(j17t*A17),0,0;
     (j3t*A3+j5t*A5),-(j4t*A4+j8t*A8+j6t*A6+j16t*A16+j15t*A15),(j14t*A14),(j10t*A10),(j12t*A12);
     (j7t*A7),0,0,-(j11t*A11+j10t*A10),0;
     0,(j8t*A8),(j19t*A19),0,-(j12t*A12+j13t*A13);
     (j9t*A9),(j15t*A15),-(j17t*A17+j14t*A14+j18t*A18+(j19t*A19)),(j11t*A11),(j13t*A13)];
VECt=[((j1t*A1*RM)); -(j1t*A1*RM); 0; 0; 0; 0];
Xt=MATt\VECt;
RWt=Xt(1);ROt=Xt(2);RLt=Xt(3);RNPt=Xt(4);RNmt=Xt(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDMt=(RM-1)*1000;dDWt=(RWt-1)*1000;dDOt=(ROt-1)*1000;dDNPt=(RNPt-1)*1000;dDNmt=(RNmt-1)*1000;dDLt=(RLt-1)*1000;
%<<<< CALCULATE APPARENT FRACTIONATIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaLt=(dDLt+1000)/(dDMt+1000);
alphaOt=(dDOt+1000)/(dDMt+1000);
%<<<< LOOP THROUGH; GRAPH TESTABLE VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.01;
for i = 2:210
HTt(i)=HTt(i-1)+add;
j1t(i)=j7t*HTt(i);                   %Hydrogen into cells
j2t(i)=j6t+j17t+j1t(i)-j7t-j5t-j9t;  %Hydrogen out of cells
%D/H ratio
MAT2t=[(j2t(i)*A2),(j16t*A16),(j18t*A18),0,0;
     -(j3t*A3+j5t*A5+j7t*A7+j2t(i)*A2+j9t*A9),(j4t*A4+j6t*A6),(j17t*A17),0,0;
      (j3t*A3+j5t*A5),-(j4t*A4+j8t*A8+j6t*A6+j16t*A16+j15t*A15),(j14t*A14),(j10t*A10),(j12t*A12);
      (j7t*A7),0,0,-(j11t*A11+j10t*A10),0;
      0,(j8t*A8),(j19t*A19),0,-(j12t*A12+j13t*A13);
      (j9t*A9),(j15t*A15),-(j17t*A17+j14t*A14+j18t*A18+(j19t*A19)),(j11t*A11),(j13t*A13)];
VEC2t=[((j1t(i)*A1*RM)); -(j1t(i)*A1*RM); 0; 0; 0; 0];
X2t=MAT2t\VEC2t;
RWt(i)=X2t(1);ROt(i)=X2t(2);RLt(i)=X2t(3);RNPt(i)=X2t(4);RNmt(i)=X2t(5);
%dD
dDMt(i)=(RM-1)*1000;dDWt(i)=(RWt(i)-1)*1000;dDOt(i)=(ROt(i)-1)*1000;dDNPt(i)=(RNPt(i)-1)*1000;dDNmt(i)=(RNmt(i)-1)*1000;dDLt(i)=(RLt(i)-1)*1000;
%alpha
alphaLt(i)=(dDLt(i)+1000)/(dDMt(i)+1000);
alphaOt(i)=(dDOt(i)+1000)/(dDMt(i)+1000);  
end
%% Solution for Fig.4 & Fig.S5 panel e) and Fig. 5 panel b) PERCENT EXUDATE RELEASE
%<<<< TESTABLE METABOLIC PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sachsXp=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
Mp=1;           %"Metabolic Reductant Parameter"    
Bp=1;           %"Biosynthetic Cell Water Parameter" 
Ep=1;           %"Equilibrium Exchange Parameter"
HTp=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PERp=0.0000001; %"Percent Extracellular Release" ****testing this paramter in loop below****
%additional parameters not tested in manuscript:
L_PERp=0.2;     %"Lipid Exudate Parameter"
Rep=0.2;        %"B-oxidation/Recyle Parameter" 
                %*Re can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7p=1;  
j10p=j7p*(1-sachsXp);  
j11p=j7p*sachsXp;   
j8p=j7p*Mp;           
j13p=j8p*(1-sachsXp);   
j9p=(j13p+j11p)*.5;    
j15p=(j13p+j11p)*.5;     
j1p=j7p*HTp;               
j16p=j7p*PERp*(1-L_PERp);            
j18p=j7p*PERp*L_PERp;               
j3p=j7p*Ep;                           
j4p=j3p;                          
j14p=(j15p+j13p+j11p+j9p-j18p)*(1-Rep-Rep);     
j17p=(j15p+j13p+j11p+j9p-j18p)*Rep;                 
j19p=(j15p+j13p+j11p+j9p-j18p)*Rep; 
j12p=(j8p*sachsXp)+j19p;                     
j5p=(j10p+j12p)*Bp;                 
j6p=j12p+j14p+j5p+j10p-j8p-j15p-j16p;   
j2p=j6p+j17p+j1p-j7p-j5p-j9p;                  
%<<<< SOLVE FOR R (D/H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MATp=[(j2p*A2),(j16p*A16),(j18p*A18),0,0;
    -(j3p*A3+j5p*A5+j7p*A7+j2p*A2+j9p*A9),(j4p*A4+j6p*A6),(j17p*A17),0,0;
     (j3p*A3+j5p*A5),-(j4p*A4+j8p*A8+j6p*A6+j16p*A16+j15p*A15),(j14p*A14),(j10p*A10),(j12p*A12);
     (j7p*A7),0,0,-(j11p*A11+j10p*A10),0;
     0,(j8p*A8),(j19p*A19),0,-(j12p*A12+j13p*A13);
     (j9p*A9),(j15p*A15),-(j17p*A17+j14p*A14+j18p*A18+(j19p*A19)),(j11p*A11),(j13p*A13)];
VECp=[((j1p*A1*RM)); -(j1p*A1*RM); 0; 0; 0; 0];
Xp=MATp\VECp;
RWp=Xp(1);ROp=Xp(2);RLp=Xp(3);RNPp=Xp(4);RNmp=Xp(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDMp=(RM-1)*1000;dDWp=(RWp-1)*1000;dDOp=(ROp-1)*1000;dDNPp=(RNPp-1)*1000;dDNmp=(RNmp-1)*1000;dDLp=(RLp-1)*1000;
%<<<< CALCULATE APPARENT FRACTIONATIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaLp=(dDLp+1000)/(dDMp+1000);
alphaOp=(dDOp+1000)/(dDMp+1000);
%<<<< LOOP THROUGH; GRAPH TESTABLE VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.01;
for i = 2:110
PERp(i)=PERp(i-1)+add;
j16p(i)=j7p*PERp(i)*(1-L_PERp);                     
j18p(i)=j7p*PERp(i)*L_PERp; 
j14p(i)=(j15p+j13p+j11p+j9p-j18p(i))*(1-Rep-Rep); 
j17p(i)=(j15p+j13p+j11p+j9p-j18p(i))*Rep;
j19p(i)=(j15p+j13p+j11p+j9p-j18p(i))*Rep;
j12p(i)=(j8p*sachsXp)+j19p(i);          
j5p(i)=(j10p+j12p(i))*Bp; 
j6p(i)=j12p(i)+j14p(i)+j5p(i)+j10p-j8p-j15p-j16p(i);
j2p(i)=j6p(i)+j17p(i)+j1p-j7p-j5p(i)-j9p; 
%D/H ratio
MAT2p=[(j2p(i)*A2),(j16p(i)*A16),(j18p(i)*A18),0,0;
     -(j3p*A3+j5p(i)*A5+j7p*A7+j2p(i)*A2+j9p*A9),(j4p*A4+j6p(i)*A6),(j17p(i)*A17),0,0;
      (j3p*A3+j5p(i)*A5),-(j4p*A4+j8p*A8+j6p(i)*A6+j16p(i)*A16+j15p*A15),(j14p(i)*A14),(j10p*A10),(j12p(i)*A12);
      (j7p*A7),0,0,-(j11p*A11+j10p*A10),0;
      0,(j8p*A8),(j19p(i)*A19),0,-(j12p(i)*A12+j13p*A13);
      (j9p*A9),(j15p*A15),-(j17p(i)*A17+j14p(i)*A14+j18p(i)*A18+(j19p(i)*A19)),(j11p*A11),(j13p*A13)];
VEC2p=[((j1p*A1*RM)); -(j1p*A1*RM); 0; 0; 0; 0];
X2p=MAT2p\VEC2p;
RWp(i)=X2p(1);ROp(i)=X2p(2);RLp(i)=X2p(3);RNPp(i)=X2p(4);RNmp(i)=X2p(5);
%dD
dDMp(i)=(RM-1)*1000;dDWp(i)=(RWp(i)-1)*1000;dDOp(i)=(ROp(i)-1)*1000;dDNPp(i)=(RNPp(i)-1)*1000;dDNmp(i)=(RNmp(i)-1)*1000;dDLp(i)=(RLp(i)-1)*1000;
%alpha
alphaLp(i)=(dDLp(i)+1000)/(dDMp(i)+1000);
alphaOp(i)=(dDOp(i)+1000)/(dDMp(i)+1000);
end
%% FIGURES
%Main text Fig. 4
% Fig. 4: Sensitivity of aLipids-Media to changes in metabolic parameters using default values where d2HMedia = 0 permil. 
% Lipid-Media 2H/1H fractionation decreases (i.e., aLipids-Media approaches unity) as 
% (a) (1-x) (the amount of hydrogen from metabolic NAD(P)H vs photosynthetic NADPH) increases, 
% (b) M (size of the metabolic NAD(P)H flux relative to PS1) increases, 
% (c) B and E increase (where B is the amount of organic-hydrogen source from cell-water-hydrogen relative to total NAD(P)H, and E is the equilibrium exchange between organic and cell-water-hydrogen), 
% (d) HT (Hydrogen Transport) is less than photosynthesis and decreases, and 
% (e) PER (Percent Extracellular Release) increases.
figure
%alpha X
subplot(1,5,1);
plot ((1-sachsX),alphaL,'-b','LineWidth',3);
axis([0 1 0.6 .9]);
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
title('(a)','FontWeight','Bold');
xlabel('1-x','FontWeight','Bold','FontAngle','Italic');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%alpha M
subplot(1,5,2);
plot (MM,alphaLM,'-b','LineWidth',3);
axis([0 10 0.6 .9]);
title('(b)','FontWeight','Bold');
xlabel('M','FontWeight','Bold','FontAngle','Italic');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%alpha B and E
subplot(1,5,3);
plot (Bb,alphaLb,'-b','LineWidth',3);
axis([0 10 0.6 .9]);
title('(c)','FontWeight','Bold');
xlabel('B & E','FontWeight','Bold','FontAngle','Italic');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%alpha HT
subplot(1,5,4);
plot (HTt,alphaLt,'-b','LineWidth',3);
axis([0 1 0.6 .9]);
xlabel('HT','FontWeight','Bold','FontAngle','Italic');
title('(d)','FontWeight','Bold');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%alpha PER
subplot(1,5,5);
plot (PERp,alphaLp,'-b','LineWidth',3);
axis([0 .7 0.6 .9]);
title('(e)','FontWeight','Bold');
xlabel('PER','FontWeight','Bold','FontAngle','Italic');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%%
%Main text Fig. 5
% Fig. 4: Sensitivity of d2HResevoir to changes in metabolic parameters using default values and d2HMedia = 0 permil. 
% All non-media boxes become 2H-enriched as (a) HT (Hydrogen Transport) decreases and 
% (b) PER (Percent Extracellular Release) increases.
figure
%HT
subplot(2,3,1);
plot (HTt,dDMt,'-k','LineWidth',5);
hold on
plot (HTt,dDWt,'--c','LineWidth',2);
plot (HTt,dDOt,'LineWidth',3,'Color',[1,0,0]);
plot (HTt,dDLt,'LineWidth',3,'Color',[0,0,1]);
plot (HTt,dDNPt,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (HTt,dDNmt,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 1 -700 200]);
xlabel('HT','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%PER
subplot(2,3,2);
plot (PERp,dDMp,'-k','LineWidth',5);
hold on
plot (PERp,dDWp,'--c','LineWidth',2);
plot (PERp,dDOp,'LineWidth',3,'Color',[1,0,0]);
plot (PERp,dDLp,'LineWidth',3,'Color',[0,0,1]);
plot (PERp,dDNPp,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (PERp,dDNmp,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 .7 -700 200]);
xlabel('PER','FontWeight','Bold','FontAngle','Italic'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%trick for the legend
subplot(2,3,3)
plot (1,1,'--c','LineWidth',4);
hold on
plot (1,1,'-k','LineWidth',5);
plot (1,1,'LineWidth',3,'Color',[1,0,0]);
plot (1,1,'LineWidth',3,'Color',[0,0,1]);
plot (1,1,'--','LineWidth',2,'Color',[1,0.5,.1]);
plot (1,1,'--','LineWidth',2,'Color',[0,1,.1]); 
axis off; set(gca,'LineWidth',2);
legend('  Cell Water','  Media','  Other Organics','  Lipids','  NAD(P)H MET','  NADPH PS1');
h = findobj('type', 'axes'); % Find all sets of axes
bcg = get(gcf,'Color'); % Setting bcg as the background color of the figure
set(h(1), 'xcolor', bcg) % Making the vertical lines blend in with the background
set(h(1), 'ycolor', bcg) % Making the horizontal lines blend in with the background
%%
%Supplementary Material Fig. S5
% Fig. S5. Changes in d2HResevoir, aLipid-Media, and aOrganics-Media as metabolic parameters vary: 
% (a) 1-x = PSI vs MET NAD(P)H Parameter, 
% (b) M = Metabolic Reductant Parameter, 
% (c) B = Biosynthetic Cell Water Parameter and E = Equilibrium Exchange Parameter, 
% (d) HT = Hydrogen Transport Parameter, 
% (e) PER = Percent Exudate Release.
figure
% X ______________________________________
%alpha
subplot(2,6,1);
plot ((1-sachsX),alphaO,'-r','LineWidth',3);
hold on
plot ((1-sachsX),alphaL,'-b','LineWidth',3);
axis([0 1 0.6 1]);
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
title('(a)','FontWeight','Bold');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx),'XTickLabel',[]);
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,6,7);
plot ((1-sachsX),dDM,'-k','LineWidth',5);
hold on
plot ((1-sachsX),dDW,'--c','LineWidth',2);
plot ((1-sachsX),dDO,'LineWidth',3,'Color',[1,0,0]);
plot ((1-sachsX),dDL,'LineWidth',3,'Color',[0,0,1]);
plot ((1-sachsX),dDNP,'--','LineWidth',2,'Color',[0,1,.1]); 
plot ((1-sachsX),dDNm,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 1 -700 200]);
xlabel('1-x','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
% M ______________________________________
%alpha
subplot(2,6,2);
plot (MM,alphaOM,'-r','LineWidth',3);
hold on
plot (MM,alphaLM,'-b','LineWidth',3);
axis([0 10 0.6 1]);
title('(b)','FontWeight','Bold');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx),'XTickLabel',[]);
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%dD
subplot(2,6,8);
plot (MM,dDMM,'-k','LineWidth',5);
hold on
plot (MM,dDWM,'--c','LineWidth',2); 
plot (MM,dDOM,'LineWidth',3,'Color',[1,0,0]);
plot (MM,dDLM,'LineWidth',3,'Color',[0,0,1]);
plot (MM,dDNPM,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (MM,dDNmM,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 10 -700 200]);
xlabel('M','FontWeight','Bold','FontAngle','Italic'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
% B & E ____________________________________
%alpha
subplot(2,6,3);
plot (Bb,alphaOb,'-r','LineWidth',3);
hold on
plot (Bb,alphaLb,'-b','LineWidth',3);
axis([0 10 0.6 1]);
title('(c)','FontWeight','Bold');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx),'XTickLabel',[]);
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%dD
subplot(2,6,9);
plot (Bb,dDMb,'-k','LineWidth',5);
hold on
plot (Bb,dDWb,'--c','LineWidth',2);
plot (Bb,dDOb,'LineWidth',3,'Color',[1,0,0]);
plot (Bb,dDLb,'LineWidth',3,'Color',[0,0,1]);
plot (Bb,dDNPb,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (Bb,dDNmb,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 10 -700 200]);
xlabel('B & E','FontWeight','Bold','FontAngle','Italic'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
% HT ______________________________________
%alpha
subplot(2,6,4);
plot (HTt,alphaOt,'-r','LineWidth',3);
hold on
plot (HTt,alphaLt,'-b','LineWidth',3);
axis([0 1 0.6 1]);
title('(d)','FontWeight','Bold');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx),'XTickLabel',[]);
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%dD
subplot(2,6,10);
plot (HTt,dDMt,'-k','LineWidth',5);
hold on
plot (HTt,dDWt,'--c','LineWidth',2);
plot (HTt,dDOt,'LineWidth',3,'Color',[1,0,0]);
plot (HTt,dDLt,'LineWidth',3,'Color',[0,0,1]);
plot (HTt,dDNPt,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (HTt,dDNmt,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 1 -700 200]);
xlabel('HT','FontWeight','Bold','FontAngle','Italic'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
% PER _________________________________________
%alpha
subplot(2,6,5);
plot (PERp,alphaOp,'-r','LineWidth',3);
hold on
plot (PERp,alphaLp,'-b','LineWidth',3);
axis([0 .7 0.6 1]);
title('(e)','FontWeight','Bold');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx),'XTickLabel',[]);
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%dD
subplot(2,6,11);
plot (PERp,dDMp,'-k','LineWidth',5);
hold on
plot (PERp,dDWp,'--c','LineWidth',2);
plot (PERp,dDOp,'LineWidth',3,'Color',[1,0,0]);
plot (PERp,dDLp,'LineWidth',3,'Color',[0,0,1]);
plot (PERp,dDNPp,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (PERp,dDNmp,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 .7 -700 200]);
xlabel('PER','FontWeight','Bold','FontAngle','Italic'); 
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy),'YTickLabel',[]);
%trick for the legend to fit
subplot(2,6,6)
plot (1,1,'-r','LineWidth',3);
hold on
plot (1,1,'-b','LineWidth',3);
axis off; set(gca,'LineWidth',2);
legend(sprintf( '%s\n%s', '   alpha', '   Organics-Media' ), sprintf( '%s\n%s', '   alpha', '   Lipids-Media' ));
h = findobj('type', 'axes'); % Find all sets of axes
bcg = get(gcf,'Color'); % Setting bcg as the background color of the figure
set(h(1), 'xcolor', bcg) % Making the vertical lines blend in with the background
set(h(1), 'ycolor', bcg) % Making the horizontal lines blend in with the background
subplot(2,6,12);
plot (1,1,'--c','LineWidth',4);
hold on
plot (1,1,'-k','LineWidth',5);
plot (1,1,'LineWidth',3,'Color',[1,0,0]);
plot (1,1,'LineWidth',3,'Color',[0,0,1]);
plot (1,1,'--','LineWidth',2,'Color',[1,0.5,.1]);
plot (1,1,'--','LineWidth',2,'Color',[0,1,.1]); 
axis off; set(gca,'LineWidth',2);
legend('  Cell Water','  Media','  Other Organics','  Lipids','  NAD(P)H MET','  NADPH PS1');
h = findobj('type', 'axes'); 
bcg = get(gcf,'Color'); 
set(h(1), 'xcolor', bcg) 
set(h(1), 'ycolor', bcg) 