% METABOLIC PARAMETER "HT"
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
figure
%alpha
subplot(2,2,1);
plot (HTt,alphaOt,'-r','LineWidth',3);
hold on
plot (HTt,alphaLt,'-b','LineWidth',3);
axis([0 2 0.6 1]);
xlabel('HT','FontWeight','Bold','FontAngle','Italic');
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
legend('   alpha Organics-Media', '   alpha Lipids-Media');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,2,2);
plot (HTt,dDMt,'-k','LineWidth',5);
hold on
plot (HTt,dDWt,'--c','LineWidth',2);
plot (HTt,dDOt,'LineWidth',3,'Color',[1,0,0]);
plot (HTt,dDLt,'LineWidth',3,'Color',[0,0,1]);
plot (HTt,dDNPt,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (HTt,dDNmt,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 2 -700 200]);
xlabel('HT','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
legend('  Media','  Cell Water','  Other Organics','   Lipids','  NADPH PS1','  Metabolic NAD(P)H');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%SAFETY
if j1t(1)<0; j1t(i)=NaN; display('NEGATIVE FLUX j1'); end;
if j2t(1)<0; j2t(i)=NaN; display('WARNING: NEGATIVE FLUX j2'); end;
if j3t(1)<0; j3t(i)=NaN; display('NEGATIVE FLUX j3'); end;
if j5t(1)<0; j5t(i)=NaN; display('NEGATIVE FLUX j5'); end;
if j6t(1)<0; j6t(i)=NaN; display('NEGATIVE FLUX j6'); end;
if j8t(1)<0; j8t(i)=NaN; display('NEGATIVE FLUX j8'); end;
if j9t(1)<0; j9t(i)=NaN; display('NEGATIVE FLUX j9'); end;
if j10t(1)<0; j10t(i)=NaN; display('NEGATIVE FLUX j10'); end;
if j12t(1)<0; j12t(i)=NaN; display('NEGATIVE FLUX j12'); end;
if j14t(1)<0; j14t(i)=NaN; display('NEGATIVE FLUX j14'); end;
if j17t(1)<0; j17t(i)=NaN; display('NEGATIVE FLUX j17'); end;
if j18t(1)<0; j18t(i)=NaN; display('NEGATIVE FLUX j18'); end;
if j19t(1)<0; j19t(i)=NaN; display('NEGATIVE FLUX j19'); end;
%check mass balance
subplot(2,2,3) 
plot(HTt,(j1t+j17t+j6t+j4t-j2t-j7t-j9t-j5t-j3t),'c','LineWidth',5); %cell water
hold on
plot(HTt,(j18t+j2t+j16t-j1t),'k','LineWidth',4); %Media
plot(HTt,(j18t+j2t+j16t-j1t),'k','LineWidth',5); %Media_no_20
plot(HTt,(j9t+j15t+j11t+j13t-(j17t+j18t+j14t+j19t)),'b','LineWidth',1); %Lipids
plot(HTt,(j5t+j14t+j12t+j10t+j3t-j6t-j15t-j8t-j16t-j4t),'r','LineWidth',2);%organics
plot(HTt,(j7t-j10t-j11t),'g','LineWidth',1); %PS1
plot(HTt,(j8t+j19t-j12t-j13t),'--m','LineWidth',1); %OPP
%legend('cWater','media','no 20','lipids','organics','PS1','OPP');
title('steady state tests');
ylabel('flux in - flux out');
axis([0 1 -1 1]);

% Checking Fluxes
% figure
% plot (HTt,j1t, HTt,j2t)
% figure
% plot(HTt,j3t,HTt,j4t,HTt,j5t,HTt,j6t,HTt,j7t,HTt,j8t,HTt,j9t)
% figure
% plot(HTt,j10t,HTt,j11t,HTt,j12t,HTt,j13t,HTt,j14t,HTt,j15t,HTt,j16t,HTt,j17t,HTt,j18t)