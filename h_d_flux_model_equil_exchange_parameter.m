% METABOLIC PARAMETER "E"
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
sachsXe=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
Me=1;           %"Metabolic Reductant Parameter"    
Be=1;           %"Biosynthetic Cell Water Parameter" 
Ee=.01;         %"Equilibrium Exchange Parameter" ****testing this paramter in loop below****
HTe=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PERe=0.03;      %"Percent Extracellular Release"
%additional parameters not tested in manuscript:
L_PERe=0.2;     %"Lipid Exudate Parameter"
Ree=0.2;        %"B-oxidation/Recyle Parameter" 
                %*Re can't result in a J19 flux bigger than (j13+j11+j9+j15)
%<<<<< FLUXES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
j7e=1;                             
j10e=j7e*(1-sachsXe);                 
j11e=j7e*sachsXe;                     
j8e=j7e*Me;                        
j13e=j8e*(1-sachsXe);                 
j9e=(j13e+j11e)/2;                 
j15e=(j13e+j11e)/2;                
j1e=j7e*HTe;                       
j16e=j7e*PERe*(1-L_PERe);        
j18e=j7e*PERe*L_PERe;            
j3e=j7e*Ee;                        
j4e=j3e;                           
j14e=(j15e+j13e+j11e+j9e-j18e)*(1-Ree-Ree);     
j17e=(j15e+j13e+j11e+j9e-j18e)*Ree;                 
j19e=(j15e+j13e+j11e+j9e-j18e)*Ree;
j12e=(j8e*sachsXe)+j19e;              
j5e=(j10e+j12e)*Be;                
j6e=j12e+j14e+j5e+j10e-j8e-j15e-j16e;   
j2e=j6e+j17e+j1e-j7e-j5e-j9e;      
%<<<< SOLVE FOR R (D/H) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDexperiment=0;
RM=((dDexperiment/1000)+1); 
MATe=[(j2e*A2),(j16e*A16),(j18e*A18),0,0;
    -(j3e*A3+j5e*A5+j7e*A7+j2e*A2+j9e*A9),(j4e*A4+j6e*A6),(j17e*A17),0,0;
     (j3e*A3+j5e*A5),-(j4e*A4+j8e*A8+j6e*A6+j16e*A16+j15e*A15),(j14e*A14),(j10e*A10),(j12e*A12);
     (j7e*A7),0,0,-(j11e*A11+j10e*A10),0;
     0,(j8e*A8),(j19e*A19),0,-(j12e*A12+j13e*A13);
     (j9e*A9),(j15e*A15),-(j17e*A17+j14e*A14+j18e*A18+(j19e*A19)),(j11e*A11),(j13e*A13)];
VECe=[((j1e*A1*RM)); -(j1e*A1*RM); 0; 0; 0; 0];
Xe=MATe\VECe;
RWe=Xe(1);ROe=Xe(2);RLe=Xe(3);RNPe=Xe(4);RNme=Xe(5);
%<<<< CONVERT TO dD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dDMe=(RM-1)*1000;dDWe=(RWe-1)*1000;dDOe=(ROe-1)*1000;dDNPe=(RNPe-1)*1000;dDNme=(RNme-1)*1000;dDLe=(RLe-1)*1000;
%<<<< CALCULATE APPARENT FRACTIONATIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
alphaLe=(dDLe+1000)/(dDMe+1000);
alphaOe=(dDOe+1000)/(dDMe+1000);
%<<<< LOOP THROUGH; GRAPH TESTABLE VARIABLE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add=.1;
for i = 2:110
Ee(i)=Ee(i-1)+add;
j3e(i)=Ee(i);                           
j4e(i)=j3e(i);
%D/H ratio
MAT2e=[(j2e*A2),(j16e*A16),(j18e*A18),0,0;
     -(j3e(i)*A3+j5e*A5+j7e*A7+j2e*A2+j9e*A9),(j4e(i)*A4+j6e*A6),(j17e*A17),0,0;
      (j3e(i)*A3+j5e*A5),-(j4e(i)*A4+j8e*A8+j6e*A6+j16e*A16+j15e*A15),(j14e*A14),(j10e*A10),(j12e*A12);
      (j7e*A7),0,0,-(j11e*A11+j10e*A10),0;
      0,(j8e*A8),(j19e*A19),0,-(j12e*A12+j13e*A13);
      (j9e*A9),(j15e*A15),-(j17e*A17+j14e*A14+(j18e*A18)+(j19e*A19)),(j11e*A11),(j13e*A13)];
VEC2e=[((j1e*A1*RM)); -(j1e*A1*RM); 0; 0; 0; 0];
X2e=MAT2e\VEC2e;
RWe(i)=X2e(1);ROe(i)=X2e(2);RLe(i)=X2e(3);RNPe(i)=X2e(4);RNme(i)=X2e(5);
%dD
dDMe(i)=(RM-1)*1000;dDWe(i)=(RWe(i)-1)*1000;dDOe(i)=(ROe(i)-1)*1000;dDNPe(i)=(RNPe(i)-1)*1000;dDNme(i)=(RNme(i)-1)*1000;dDLe(i)=(RLe(i)-1)*1000;
%alpha
alphaLe(i)=(dDLe(i)+1000)/(dDMe(i)+1000);
alphaOe(i)=(dDOe(i)+1000)/(dDMe(i)+1000);
end
figure
%alpha
subplot(2,2,1);
plot (Ee,alphaOe,'-r','LineWidth',3);
hold on
plot (Ee,alphaLe,'-b','LineWidth',3);
axis([0 10 0.6 1]);
xlabel('E','FontWeight','Bold','FontAngle','Italic');
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
legend('   alpha Organics-Media', '   alpha Lipids-Media');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,2,2);
plot (Ee,dDMe,'-k','LineWidth',5);
hold on
plot (Ee,dDWe,'--c','LineWidth',2);
plot (Ee,dDOe,'LineWidth',3,'Color',[1,0,0]);
plot (Ee,dDLe,'LineWidth',3,'Color',[0,0,1]);
plot (Ee,dDNPe,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (Ee,dDNme,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 10 -700 200]);
xlabel('E','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
legend('  Media','  Cell Water','  Other Organics','   Lipids','  NADPH PS1','  Metabolic NAD(P)H');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%check mass balance
subplot(2,2,3)
plot(Ee,(j1e+j17e+j6e+j4e-j2e-j7e-j9e-j5e-j3e),'c','LineWidth',5); %cell water
hold on
plot(Ee,(j18e+j2e+j16e-j1e),'k','LineWidth',4); %Media
plot(Ee,(j9e+j15e+j11e+j13e-(j17e+j18e+j14e+j19e)),'b','LineWidth',1); %Lipids
plot(Ee,(j5e+j14e+j12e+j10e+j3e-j6e-j15e-j8e-j16e-j4e),'r','LineWidth',2);%organics
plot(Ee,(j7e-j10e-j11e),'g','LineWidth',1); %PS1
plot(Ee,(j8e+j19e-j12e-j13e),'--m','LineWidth',1); %MET
%legend('cWater','media','mediano20','lipids','organics','PS1','MET');
title('steady state tests');
ylabel('flux in - flux out');
axis([0 1 -1 1]);
