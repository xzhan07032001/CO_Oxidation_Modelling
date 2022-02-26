%Name:Xin (Jason) Zhang, Creation Date:6/26/2021, Updated Since:9/30/2021 Version:R2019B

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
Temperature = 300;  %Kelvin
M_O2 = 32;          %Molar mass of O2
M_CO = 28;          %Molar mass of CO
Area = 5^2;         %Area per unit cell [Angstrom^2]
PP_CO = 1012.5;     %Partial Pressure of CO [mbar]
PP_O2 = 1012.5;    %Partial Pressure of O2 [mbar]

P_CO = adsorption(M_CO, PP_CO, Area, Temperature); %CO Adsorption
d_CO = Arrhenius(1.097407, Temperature);           %CO Desorption
h_CO = Arrhenius(0.699167, Temperature);           %CO Hopping

P_O2 = adsorption(M_O2, PP_O2, Area, Temperature); %O Hopping
d_O2 = Arrhenius(2.8683, Temperature);             %O Desorption
h_O = Arrhenius(0.6075,  Temperature);             %O Hopping

k_OCO = Arrhenius(1.08,Temperature);               %LH Reaction original is 1.08
e_OCO = 0*Arrhenius(0.45833,Temperature);       %ER Reaction 0.54 0.45833 is original 0.00005
e_OCOPREF = Arrhenius(0.45833,Temperature);

PP_CORR = roughrate(PP_CO);                        %CO Adsorption (Approx.)
PP_O2RR = roughrate(PP_O2);                        %O2 Adsorption (Approx.) 
%% Master Equations

%{
x(1)=(O1)=(O3) x(2)=(O2) x(3)=(CO1)=(CO3)         
x(4)=(CO2) x(5)=(E1)
   
Simplifications:
E2 = 1-x(2)-x(4)
E1 = (1-x(1)-x(3))
O1 = O3
CO1=CO3
%}

fx = @(t,x) [h_O*(1-x(1)-x(3))*x(2)+P_O2*(1-x(1)-x(3))*(1-x(2)-x(4))-d_O2*x(1)*x(2)-h_O*x(1)*(1-x(2)-x(4))-k_OCO*x(1)*x(4);          %d/dt=O1 [1] COMPLETE
             h_O*(1-x(2)-x(4))*(x(1)+x(5))+P_O2*2*(1-x(2)-x(4))*((1-x(1)-x(3))+(1-x(5)-x(6)))-d_O2*x(1)*(x(2)+x(6))-h_O*x(2)*(1-x(5)-x(6))-k_OCO*x(2)*(x(3)+x(6));%d/dt=O2 [2] COMPLETE 
             P_CO*(1-x(1)-x(3))+h_CO*(1-x(1)-x(3))*x(4)-d_CO*x(3)-k_OCO*x(3)*x(2)-h_CO*x(3)*(1-x(2)-x(4));                                           %d/dt=CO1[3] COMPLETE
             P_CO*(1-x(2)-x(4))+h_CO*(1-x(2)-x(4))*(x(3)+x(4))-d_CO*x(4)-k_OCO*2*x(4)*x(1)-h_CO*x(3)*((1-x(1)-x(3))+(1-x(5)-x(6)));                                     %d/dt=CO2[4]
             h_O*(1-x(5)-x(6))*x(2)+P_O2*(1-x(2)-x(4))*(1-x(5)-x(6))-d_O2*x(2)*x(5)-h_O*x(5)*(1-x(2)-x(4))-k_OCO*x(5)*x(4);  %d/dt=O3 [5] CHECK
             P_CO*(1-x(5)-x(6))+h_CO*((1-x(5)-x(6))*x(4)+(1-x(2)-x(4))*x(6))-d_CO*x(6)-k_OCO*x(6)*x(2);%d/dt=CO3 [6]
                                  ];
%h_O*((1-x(5)-x(6))*x(2)+((1-x(2)-x(4))*x(5))+P_O2*((1-x(5)-x(6))*(1-x(2)-x(4))))-d_O2*x(2)*x(5)-e_OCO*(x(5)*x(4)+x(5)*x(2))-k_OCO*x(5)*x(4);
t_int = [0,1];                      %Time-Interval [Seconds] original is .005
initial_cond = [0,0,0,0,0,0];           %Initial Conditions from x(1)-x(18)
[t,xa] = ode15s(fx,t_int,initial_cond); %Integrates Equations [1]-[18]
x1 = xa(:,1);   %x(1) P(O1)
x2 = xa(:,2);   %x(2) P(O2)
x3 = xa(:,3);   %x(3) P(CO1)
x4 = xa(:,4);   %x(4) P(CO2)
x5 = xa(:,5);   %x(5) P(O3)
x6 = xa(:,6);   %x(6) P(CO3)
x7 = 1-x1-x3;   %P(E1)
x8 = 1-x2-x4;   %P(E2)
x9 = 1-x5-x6;   %P(E3)
Site_1 = [x1 x3 x7];
Site_2 = [x2 x4 x8];
Site_3 = [x5 x6 x9];
%38 is .0005 44 is .001 50 is .002 seconds
%Table_Value = [ x1(15) x1(20) x1(25) x1(32)
%                x2(15) x2(20) x2(25) x2(32)
%                x3(15) x3(20) x3(25) x3(32)
%                x4(15) x4(20) x4(25) x4(32)
%                x5(15) x5(20) x5(25) x5(32)
%                x6(15) x6(20) x6(25) x6(32) 
%                x7(15) x7(20) x7(25) x7(32) 
%                x8(15) x8(20) x8(25) x8(32) 
%                x9(15) x9(20) x9(25) x9(32)];
%disp(Table_Value);
%% Graphing
figure(1);
tile = tiledlayout(3,1);

ax1 = nexttile; %Site 1
plot(t,Site_1,'--','Linewidth',2);
xlim([-0.00001,1.002/1]);
ylim([-0.01,1.05]);
lgd1 = legend('$P(O_1)$','$P(CO_1)$','$P(E_1)$','Interpreter','latex','Location','northeast','Fontsize',20);
title("Probability Vs Time",'Fontsize',20,'Fontname','times');
grid(ax1,'on')

ax2 = nexttile;
plot(t,Site_2,'--','Linewidth',2);
xlim([-0.00001,1.002/1]);
ylim([-0.01,1.05]);
lgd2 = legend('$P(O_2)$','$P(CO_2)$','$P(E_2)$','Interpreter','latex','Location','northeast','Fontsize',20);
grid(ax2,'on')

ax3 = nexttile;
plot(t,Site_3,'--','Linewidth',2);
xlim([-0.00001,1.002/1]);
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
