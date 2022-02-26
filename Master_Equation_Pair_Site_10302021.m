%Name:Xin (Jason) Zhang, Creation Date:6/26/2021, Updated Since:6/26/2021 Version:R2019B

%% Preliminary Information
%Context:
%{
System consists of a single-site, in which each site has 3 possible states:
->Carbon Monoxide (CO)
->Oxygen (O)
->Empty (CO)
There are 9 possible site/state combinations.

The reaction is CO oxidation on an RuO2 facet (Catalytic Reaction)
Two types of reaction mechanisms:
->Langmuir-Hinshelwood (LH)
    -Reactants are adsorbed into neighboring sites and react
    -CO and O on neighboring sites form CO2
->Eley-Rideal (ER)
    -One of the molecules adsorb to a site and react with a gas-phase one
    -O on a site reacting with gas-phase CO to form CO2(Can't do other way)
%}

%Table of Contents: 
%{
Code sections are as followed:
-Preliminary Information
-Parameters (All the Rates)
-Master Equations(Differential Equation Solver)
-Graphing 
-Functions (Rate Functions)
%}
%% Parameters
Temperature = 300; %Kelvin
M_O2 = 32;         %Molar mass of O2
M_CO = 28;         %Molar mass of CO
Area = 5^2;        %Area per unit cell [Angstrom^2]
PP_CO = 506.25;    %Partial Pressure of CO [mbar]
PP_O2 = 5060.25;    %Partial Pressure of O2 [mbar]

P_CO = adsorption(M_CO, PP_CO, Area, Temperature); %CO Adsorption
d_CO = Arrhenius(1.097407, Temperature);           %CO Desorption
h_CO = Arrhenius(0.699167, Temperature);           %CO Hopping

P_O2 = adsorption(M_O2, PP_O2, Area, Temperature); %O Hopping
d_O2 = Arrhenius(2.8683, Temperature);             %O Desorption
h_O = Arrhenius(0.6075,  Temperature);             %O Hopping

k_OCO = Arrhenius(1.08,Temperature);               %LH Reaction original is 1.08
e_OCO = 0.00005*Arrhenius(0.45833,Temperature);    %ER Reaction 0.54 0.4175
%% Master Equations

%{
x(1)=(O1,O2)  x(2)=(O1,CO2)  x(3)=(O1,E2)  x(4)=(O2,O3)   x(5)=(O2,CO3)   x(6)=(O2,E3)        
x(7)=(CO1,O2) x(8)=(CO1,CO2) x(9)=(CO1,E2) x(10)=(CO2,O3) x(11)=(CO2,CO3) x(12)=(CO2,E3)
x(13)=(E1,O2) x(14)=(E1,CO2) x(15)=(E2,O3)  x(16)=(E2,CO3) 
   
%}
%(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))
%(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))

%16 turned into 15
%17 turned into 16
%18 turned into 17
fx = @(t,x) [ P_O2*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))+h_O*(x(3)*x(15)-x(1)*x(6))-d_O2*x(1)-k_OCO*x(1)*x(5)-2*e_OCO*x(1);                                  %(O1,O2) [1]
              P_CO*x(3)+h_CO*(x(3)*x(16)-x(2)*x(12))-d_CO*x(2)-k_OCO*2*x(2)*x(10)-e_OCO*x(2);                                %(O1,CO2)[2]
              d_CO*x(2)+d_O2*x(1)*x(4)+k_OCO*(x(1)*x(5)+x(2)*x(10))+h_CO*(x(2)*x(12)-x(3)*x(16))+h_O*(x(1)*x(6)+x(13)-x(3)*x(15))+e_OCO*(x(1)-x(3));%(O1,E2)[3]
              P_O2*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))+h_O*(x(3)*x(15)-x(13)*x(4))-d_O2*x(4)-k_OCO*x(7)*x(4)-2*e_OCO*x(4);                                 %(O2,O3)[4]
              P_CO*x(6)+P_O2*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))*x(16)+h_O*(x(3)*x(16)-x(13)*x(5))-d_CO*x(5)-k_OCO*((1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))+x(7)*x(5))-e_OCO*x(5)-d_O2*x(1)*x(2)*x(5); %(O2,CO3)[5]
              d_CO*x(5)+P_O2*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))+h_O*(x(3)*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))-2*x(13)*x(6))-d_O2*x(1)*x(6)+e_OCO*(x(4)-x(6))-P_CO*x(6)-k_OCO*x(7)*x(12);    %(O2,E3) [6]
              P_CO*x(13)+P_O2*x(9)*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))+h_O*(x(9)*x(15)-x(7)*x(6))-d_O2*x(7)*x(4)-d_CO*x(7)-k_OCO*2*x(7)*x(5)-e_OCO*x(7);   %(CO1,O2) [7]
              P_CO*(x(14)+x(9))-h_CO*x(8)*x(12)-d_CO*2*x(8)-k_OCO*x(7)*x(10);                                                %(CO1,CO2)[8]
              P_CO*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))+d_O2*x(7)*x(4)+k_OCO*x(8)*x(10)+e_OCO*x(7)+h_CO*(x(14)-2*x(9)*x(16))+h_O*(x(7)*x(6)-x(9)*x(15))-d_CO*x(9)-P_O2*x(9)*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)));%(CO1,E2) [9]
              P_CO*x(15)+h_CO*(x(9)*x(15)-x(14)*x(10))-d_CO*x(7)-k_OCO*(2*x(2)*x(10))-e_OCO*x(10);                           %(CO2,O3) [10]
              P_CO*(x(16)+x(12))+h_CO*(x(9)*x(16)-x(14)*x(11))-d_CO*2*x(11)-k_OCO*x(2)*x(11);                                %(CO2,CO3)[11]
              P_CO*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))+e_OCO*x(10)+h_CO*(x(9)*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))-x(14)*x(12))-d_CO*x(12)-k_OCO*x(2)*x(12);                              %(CO2,E3) [12]
              d_CO*x(7)+P_O2*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))+h_O*(x(3)*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))-2*x(13)*x(4))-e_OCO*x(13)-k_OCO*x(13)*x(5);                         %(E1,O2)  [13]
              P_CO*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))+e_OCO*x(2)+h_CO*(2*x(9)*x(16)-2*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))*x(12))+d_CO*(x(8)-x(14))-k_OCO*x(14)*x(10);                   %(E1,CO2) [14]
              d_CO*x(10)+d_O2*x(1)*x(4)+h_O*(x(6)-x(13)*x(6))+e_OCO*(x(4)-x(15))-P_CO*x(15)-P_O2*(1-(x(1)+x(2)+x(3))-(x(7)+x(8)+x(9))-(x(13)+x(14)))*x(15)+k_OCO*(x(7)*x(4)+x(2)*x(10))+h_CO*x(14)*x(10);%(E2,O3) [15]
              P_CO*(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))+e_OCO*x(5)+h_CO*(x(12)-x(16)-x(9)*x(16))+d_CO*(x(11)-x(16))+k_OCO*(x(7)*x(5)+x(2)*x(11))+e_OCO*x(5)+d_O2*x(1)*x(5); %(E2,CO3) [16] 
                                  ];

                             
t_int = [0,0.005];                      %Time-Interval [Seconds]
initial_cond = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];   %Initial Conditions from x(1)-(1-(x(4)+x(5)+x(6))-(x(10)+x(11)+x(12))-(x(15)+x(16)))
[t,xa] = ode15s(fx,t_int,initial_cond); %Integrates Equations [1]-[18]
x1 = xa(:,1);     %x(1)
x2 = xa(:,2);     %x(2)
x3 = xa(:,3);     %x(3)
x4 = xa(:,4);     %x(4)
x5 = xa(:,5);     %x(5)
x6 = xa(:,6);     %x(6)
x7 = xa(:,7);     %x(7)
x8 = xa(:,8);     %x(8)
x9 = xa(:,9);     %x(9)
x10 = xa(:,10);   %x(10)
x11 = xa(:,11);   %x(11)
x12 = xa(:,12);   %x(12)
x13 = xa(:,13);   %x(13)
x14 = xa(:,14);   %x(14)
x15 = xa(:,15);   %x(15)
x16 = xa(:,16);   %x(16)

O1 = x1+x2+x3;              %O1 Total 
O2 = x4+x5+x6;              %O2 Total 
O3 = x4+x10+x15;            %O3 Total 

CO1 = x7+x8+x9;             %CO1 Total
CO2 = x10+x11+x12;          %CO2 Total
CO3 = x5+x11+x16;           %CO3 Total

E1 = 1-O1-CO1;              %E1 Total
E2 = 1-O2-CO2;              %E2 Total
E3 = 1-O3-CO3;              %E3 Total

Site_1 = [O1 CO1 E1];
Site_2 = [O2 CO2 E2];
Site_3 = [O3 CO3 E3];
%% Graphing
figure(1);
tile = tiledlayout(3,1);

ax1 = nexttile; %Site 1
plot(t,Site_1,'--','Linewidth',2);
xlim([-0.00001,0.003]);
ylim([-0.01,1.05]);
lgd1 = legend('$P(O_1)$','$P(CO_1)$','$P(E_1)$','Interpreter','latex','Location','northeast','Fontsize',20);
title("Probability Vs Time",'Fontsize',20,'Fontname','times');
grid(ax1,'on')

ax2 = nexttile; %Site 2
plot(t,Site_2,'--','Linewidth',2);
xlim([-0.00001,0.003]);
ylim([-0.01,1.05]);
lgd2 = legend('$P(O_2)$','$P(CO_2)$','$P(E_2)$','Interpreter','latex','Location','northeast','Fontsize',20);
grid('on')

ax3 = nexttile; %Site 3
plot(t,Site_3,'--','Linewidth',2);
xlim([-0.00001,0.003]);
ylim([-0.01,1.05]);
lgd3 = legend('$P(O_3)$','$P(CO_3)$','$P(E_3)$','Interpreter','latex','Location','northeast','Fontsize',20);


%Graph Formatting
linkaxes([ax1,ax2,ax3],'xy');
ylabel(tile,"Probability [P(x)]",'Fontsize',25,'fontname','times');
xlabel(tile,"Time [s]",'Fontsize',25,'fontname','times');
ax1.FontSize = 20;
ax2.FontSize = 20;
ax3.FontSize = 20;
ax1.FontName = 'times';
ax2.FontName = 'times';
ax3.FontName = 'times';
htitle1 = get(lgd1,'Title');
htitle2 = get(lgd2,'Title');
htitle3 = get(lgd3,'Title');
set(htitle1,'String','Site-1 Probability');
set(htitle2,'String','Site-2 Probability');
set(htitle3,'String','Site-3 Probability');
grid('on');

%% Functions

function a = Arrhenius(Ea,T)     %(Activation Energy, Temp.)
    Nu = 10^13;                  %Prefactor range is 10^12-10^13 1/second
    kB = 8.617333262145*(10^-5); %Boltzman Constant in eV/T
    a = Nu*exp(-Ea/(kB*T));      %Arrhenius Equation
end

function f = adsorption(MM, PP, A, T) %(Molar Mass, Partial-P, Area, Temp)
    kB = 8.617333262145*(10^-5);      %Boltzman Constant in eV/T
    f = (A*PP)/sqrt(2*pi*MM*kB*T);    %Adsorption Rate Equation
end

function e = roughrate(P)
    e = P*10^5;
end