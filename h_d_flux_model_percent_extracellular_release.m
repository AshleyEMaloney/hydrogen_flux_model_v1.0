% METABOLIC PARAMETER "PER"
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
sachsXp=.5;     %"PS1 vs MET NAD(P)H Parameter"                 
Mp=1;           %"Metabolic Reductant Parameter"    
Bp=1;           %"Biosynthetic Cell Water Parameter"
Ep=1;           %"Equilibrium Exchange Parameter"
HTp=1;          %"Hydrogen Transfer Parameter" *can't be smaller than "PER"!
PERp=0.0000001; %"Percent Extracellular Release"  *** testing this parameter in loop below ***
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
figure
%alpha
subplot(2,2,1);
plot (PERp,alphaOp,'-r','LineWidth',3);
hold on
plot (PERp,alphaLp,'-b','LineWidth',3);
axis([0 1 0.6 1]);
xlabel('PER','FontWeight','Bold','FontAngle','Italic');
ylabel('a','FontName','Symbol','FontWeight','Bold'); 
legend('   alpha Organics-Media', '   alpha Lipids-Media');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%dD
subplot(2,2,2);
plot (PERp,dDMp,'-k','LineWidth',5);
hold on
plot (PERp,dDWp,'--c','LineWidth',2);
plot (PERp,dDOp,'LineWidth',3,'Color',[1,0,0]);
plot (PERp,dDLp,'LineWidth',3,'Color',[0,0,1]);
plot (PERp,dDNPp,'--','LineWidth',2,'Color',[0,1,.1]); 
plot (PERp,dDNmp,'--','LineWidth',2,'Color',[1,0.5,.1]); 
axis([0 1 -700 200]);
xlabel('PER','FontWeight','Bold','FontAngle','Italic'); 
ylabel('d','FontName','Symbol'); 
legend('  Cell Water','  Media','  Other Organics','   Lipids','  NADPH PS1','  Metabolic NAD(P)H');
NumTicksx=3;L=get(gca,'XLim');set(gca,'XTick',linspace(L(1),L(2),NumTicksx));
NumTicksy=3;K=get(gca,'YLim');set(gca,'YTick',linspace(K(1),K(2),NumTicksy));
%check mass balance
subplot(2,2,3) 
plot(PERp,(j1p+j17p+j6p+j4p-j2p-j7p-j9p-j5p-j3p),'c','LineWidth',5); %cell water
hold on
plot(PERp,(j18p+j2p+j16p-j1p),'k','LineWidth',4); %Media
plot(PERp,(j9p+j15p+j11p+j13p-(j17p+j18p+j14p+j19p)),'b','LineWidth',1); %Lipids
plot(PERp,(j5p+j14p+j12p+j10p+j3p-j6p-j15p-j8p-j16p-j4p),'r','LineWidth',2);%organics
plot(PERp,(j7p-j10p-j11p),'g','LineWidth',1); %PS1
plot(PERp,(j8p+j19p-j12p-j13p),'--m','LineWidth',1); %MET
title('steady state tests');
ylabel('flux in - flux out'); 
axis([0 1 -1 1]);