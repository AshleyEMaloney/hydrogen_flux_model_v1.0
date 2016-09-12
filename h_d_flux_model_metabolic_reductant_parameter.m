% METABOLIC PARAMETER "M"
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
%%
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
figure
%alpha
subplot(2,2,1);
plot (MM,alphaOM,'-r','LineWidth',3);
hold on
plot (MM,alphaLM,'-b','LineWidth',3);
axis([0 10 0.6 1]);
xlabel('M','FontWeight','Bold','FontAngle','Italic');
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
legend('   alpha Organics-Media', '   alpha Lipids-Media');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,2,2);
plot (MM,dDMM,'-k','LineWidth',5); 
hold on
plot (MM,dDWM,'--c','LineWidth',2); 
plot (MM,dDOM,'LineWidth',3,'Color',[1,0,0]);
plot (MM,dDLM,'LineWidth',3,'Color',[0,0,1]);
plot (MM,dDNPM,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (MM,dDNmM,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 10 -700 200]);
xlabel('M','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
legend('  Media','  Cell Water','  Other Organics','   Lipids','  NADPH PS1','  Metabolic NAD(P)H');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
subplot(2,2,3) %check mass balance
plot(MM,(j1M+j17M+j6M+j4M-j2M-j7M-j9M-j5M-j3M),'c','LineWidth',5); %cell water
hold on
plot(MM,(j18M+j2M+j16M-j1M),'k','LineWidth',4); %Media
plot(MM,(j9M+j15M+j11M+j13M-(j17M+j18M+j14M+j19M)),'b','LineWidth',1); %Lipids
plot(MM,(j5M+j14M+j12M+j10M+j3M-j6M-j15M-j8M-j16M-j4M),'r','LineWidth',2);%organics
plot(MM,(j7M-j10M-j11M),'g','LineWidth',1); %PS1
plot(MM,(j8M+j19M-j12M-j13M),'--m','LineWidth',1); %MET
%legend('cWater','media','lipids','organics','PS1','MET');
title('steady state tests');
ylabel('flux in - flux out');
axis([0 1 -1 1]);