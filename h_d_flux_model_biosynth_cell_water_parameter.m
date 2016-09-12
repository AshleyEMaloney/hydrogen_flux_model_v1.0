% METABOLIC PARAMETER "B"
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
sachsXb=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
Mb=1;           %"Metabolic Reductant Parameter"    
Bb=.00001;      %"Biosynthetic Cell Water Parameter" ****testing this paramter in loop below****
Eb=1;           %"Equilibrium Exchange Parameter"
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
figure
%alpha
subplot(2,2,1);
plot (Bb,alphaOb,'-r','LineWidth',3);
hold on
plot (Bb,alphaLb,'-b','LineWidth',3);
axis([0 10 0.6 1]);
xlabel('B','FontWeight','Bold','FontAngle','Italic');
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
legend('   alpha Organics-Media', '   alpha Lipids-Media');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,2,2);
plot (Bb,dDMb,'-k','LineWidth',5);
hold on
plot (Bb,dDWb,'--c','LineWidth',2);
plot (Bb,dDOb,'LineWidth',3,'Color',[1,0,0]);
plot (Bb,dDLb,'LineWidth',3,'Color',[0,0,1]);
plot (Bb,dDNPb,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (Bb,dDNmb,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 10 -700 200]);
xlabel('B','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
legend('  Media','  Cell Water','  Other Organics','   Lipids','  NADPH PS1','  Metabolic NAD(P)H');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%check mass balance
subplot(2,2,3) 
plot(Bb,(j1b+j17b+j6b+j4b-j2b-j7b-j9b-j5b-j3b),'c','LineWidth',5); %cell water
hold on
plot(Bb,(j18b+j2b+j16b-j1b),'k','LineWidth',4); %Media
plot(Bb,(j9b+j15b+j11b+j13b-(j17b+j18b+j14b+j21b)),'b','LineWidth',4); %Lipids
plot(Bb,(j5b+j14b+j12b+j10b+j3b-j6b-j15b-j8b-j16b-j4b),'r','LineWidth',2);%organics
plot(Bb,(j7b-j10b-j11b),'g','LineWidth',1); %PS1
plot(Bb,(j8b+j21b-j12b-j13b),'--m','LineWidth',1); %MET
title('steady state tests');
ylabel('flux in - flux out');
axis([0 1 -1 1]);