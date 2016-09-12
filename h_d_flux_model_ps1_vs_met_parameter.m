% METABOLIC PARAMETER "x"
% 2016 Ashley E. Maloney (ashjames@uw.edu)
% Stationary 6-box model of hydrogen and deuterium cycling in phytoplankton 
% For manuscript "Exploring lipid ²H/¹H fractionation mechanisms in 
% response to salinity with continuous cultures of the diatom Thalassiosira
% pseudonana" published in Organic Geochemistry by A.E. Maloney, A.L.C. Shinneman, 
% K. Hemeon, and J.P. Sachs of the University of Washington. http://dx.doi.org/10.1016/j.orggeochem.2016.08.015
% Tests sensitivity of dDlipid & other pools to metabolic parameters 
% (D = ²H = deuterium)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
%%
%<<<< FRACTIONATION FACTORS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
%<<<< TESTABLE METABOLIC PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sachsX=0.001;  %"PS1 vs MET NAD(P)H Parameter" 
               % as in Sachs and Kawka 2015 (Plos ONE) 
               % see also Sachs et al. 2016 (GCA)
               % *** testing this parameter in loop below ***
               % when testing other parameters default sachsX=0.5 
               % bigger "sachsX" means "Lipids" gets both
               %   -> larger D-depleted PSI-NADPH flux (J11)
               %   -> smaller D-enriched MET-NAD(P)H flux (J13)
               % and the opposit for "Other Organics"
                
M=1;           %"Metabolic Reductant Parameter"
B=1;           %"Biosynthetic Cell Water Parameter"
E=1;           %"Equilibrium Exchange Parameter"
HT=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PER=0.03;      %"Percent Extracellular Release"
%additional parameters not tested in manuscript:
L_PER=0.2;     %"Lipid Exudate Parameter"
Re=0.2;        %"B-oxidation/Recyle Parameter" 
               %*Re can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7=1;                             %JWAT_PSI: cell-water-H to PSI-NADPH
j10=j7*(1-sachsX);                %JPSI_Org: PS1-NADPH to other-organic-H
j11=j7*sachsX;                    %JPS1_Lip: PS1-NADPH to Lipid-H
j8=j7*M;                          %JORG_MET: other-orgnaic-H to MET-NAD(P)H 
j13=j8*(1-sachsX);                %JMET_Lip: MET-NAD(P)H to Lipid-H
j9=(j13+j11)/2;                   %JBIO_Wat_Lip: cell-water-H to Lipid-H
j15=(j13+j11)/2;                  %JBIO_Org_Lip: other-organic-H to Lipid-H
j1=j7*HT;                         %Jin: media-H to cell-water-H
j16=j7*PER*(1-L_PER);             %Jrelease_Org: other-organic-H to media-H
j18=j7*PER*L_PER;                 %Jrelease_Lip: Lipid-H to cell-water-H
j3=j7*E;                          %Jeq_Org_Wat: cell-water-H to other-organics-H
j4=j3;                            %Jeq_Wat_Org: other-organic-H to cell-water-H
j14=(j15+j13+j11+j9-j18)*(1-2*Re);%JCAT_Lip_Org: Lipid-H to other-organic-H (Mass Balance)
j17=(j15+j13+j11+j9-j18)*Re;      %JCAT_Lip_Wat: Lipid-H to cell-water-H (Mass Balance)
j19=(j15+j13+j11+j9-j18)*Re;      %JCAT_Lip_Met: Lipid-H to MET-NAD(P)H (Mass Balance)
j12=(j8*sachsX)+j19;              %JMET_Org: MET-NAD(P)H to other-organic-H 
j5=(j10+j12)*B;                   %JBIO_Wat_Org: cell-water-H to other-organics-H
j6=j12+j14+j5+j10-j8-j15-j16;     %JCAT_Org_Wat: other-organics-H to cell-water-H (Mass Balance)
j2=j6+j17+j1-j7-j5-j9;            %Jout: cell-water-H to media-H (Mass Balance)
%<<<< SOLVE FOR R (D/H = ²H/¹H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MAT=[(j2*A2),(j16*A16),(j18*A18),0,0;
    -(j3*A3+j5*A5+j7*A7+j2*A2+j9*A9),(j4*A4+j6*A6),(j17*A17),0,0;
     (j3*A3+j5*A5),-(j4*A4+j8*A8+j6*A6+j16*A16+j15*A15),(j14*A14),(j10*A10),(j12*A12);
     (j7*A7),0,0,-(j11*A11+j10*A10),0;
     0,(j8*A8),(j19*A19),0,-(j12*A12+j13*A13);
     (j9*A9),(j15*A15),-(j17*A17+j14*A14+j18*A18+(j19*A19)),(j11*A11),(j13*A13)];
VEC=[((j1*A1*RM)); -(j1*A1*RM); 0; 0; 0; 0];
SOLVE=MAT\VEC;
RW=SOLVE(1);RO=SOLVE(2);RL=SOLVE(3);RNP=SOLVE(4);RNm=SOLVE(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDM=(RM-1)*1000;    %Media
dDW=(RW-1)*1000;    %Cell Water
dDO=(RO-1)*1000;    %Other Organics
dDNP=(RNP-1)*1000;  %Photosynthetic NADPH
dDNm=(RNm-1)*1000;  %Metabolic NAD(P)H
dDL=(RL-1)*1000;    %Lipids
%<<<< CALCULATE APPARENT FRACTIONATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaL=(dDL+1000)/(dDM+1000); %fractionation b/t Media & Lipids
alphaO=(dDO+1000)/(dDM+1000); %fractionation b/t Media & Other Organics
%<<<< LOOP THROUGH: INCREASE "sachsX" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.1;
for i = 2:11
sachsX(i)=sachsX(i-1)+add;
j10(i)=j7*(1-sachsX(i));         
j11(i)=j7*sachsX(i);
j13(i)=j8*(1-sachsX(i));
j15(i)=(j13(i)+j11(i))/2;
j9(i)=(j13(i)+j11(i))/2;
j14(i)=(j15(i)+j13(i)+j11(i)+j9(i)-j18)*(1-Re-Re);
j17(i)=(j15(i)+j13(i)+j11(i)+j9(i)-j18)*Re;
j19(i)=(j15(i)+j13(i)+j11(i)+j9(i)-j18)*Re;
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
SOLVE2=MAT2\VEC2;
RW(i)=SOLVE2(1);RO(i)=SOLVE2(2);RL(i)=SOLVE2(3);RNP(i)=SOLVE2(4);RNm(i)=SOLVE2(5);
%dD
dDM(i)=(RM-1)*1000;
dDW(i)=(RW(i)-1)*1000;
dDO(i)=(RO(i)-1)*1000;
dDNP(i)=(RNP(i)-1)*1000;
dDNm(i)=(RNm(i)-1)*1000;
dDL(i)=(RL(i)-1)*1000;
%alpha
alphaL(i)=(dDL(i)+1000)/(dDM(i)+1000);
alphaO(i)=(dDO(i)+1000)/(dDM(i)+1000);
end
%<<<< GRAPH ALPHA & dD & TEST MASS BALANCE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
figure
%alpha
subplot(2,2,1);
plot ((1-sachsX),alphaO,'-r','LineWidth',3);
hold on
plot ((1-sachsX),alphaL,'-b','LineWidth',3);
axis([0 1 0.6 1.2]);
xlabel('1-x','FontWeight','Bold','FontAngle','Italic');
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
legend(' alpha Organics-Media', ' alpha Lipids-Media');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,2,2);
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
legend(' Media',' Cell Water',' Other Organics',' Lipids',' NADPH PS1',' Metabolic NAD(P)H');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%check mass balance: what goes in must come out for stationary assumption
subplot(2,2,3) 
plot((1-sachsX),(j1+j17+j6+j4-j2-j7-j9-j5-j3),'c','LineWidth',5); %Cell Water
hold on
plot((1-sachsX),(j18+j2+j16-j1),'k','LineWidth',4); %Media
plot((1-sachsX),(j9+j15+j11+j13-(j17+j18+j14+j19)),'b','LineWidth',1); %Lipids
plot((1-sachsX),(j5+j14+j12+j10+j3-j6-j15-j8-j16-j4),'r','LineWidth',2);%Other Organics
plot((1-sachsX),(j7-j10-j11),'g','LineWidth',1); %Photosynthetic NADPH
plot((1-sachsX),(j8+j19-j12-j13),'--m','LineWidth',1); %Metabolic NAD(P)H
%legend('cWater','media','lipids','organics','NADPH PS1','NAD(P)H MET');
title('steady state tests');
ylabel('flux in - flux out');
axis([0 1 -1 1]);