%Name:Xin (Jason) Zhang, Creation Date:6/17/2021, Updated Since:9/17/2021
%Version:R2021a

%% Preliminary Information
%Context:
%{
System consists a triple-site, in which each site has 3 possible states:
->Carbon Monoxide (CO)
->Oxygen (O)
->Empty (CO)
There are 3^3 = 27 possible site/state combinations.

The reaction is CO oxidation on an RuO2 (1 1 1) facet (Catalytic Reaction)
The two types of reaction mechanisms:
->Langmuir-Hinshelwood (LH)
    -Reactants are adsorbed into neighboring sites and react
    -CO and O on neighboring sites form CO2
->Eley-Rideal (ER) [Reaction type is currently assumed to not occur]
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
Area = sqrt(25);   %Area per unit cell
PP_CO = 506.25*2;       %Partial Pressure of CO [mbar] 
PP_O2 = 506.25*2;       %Partial Pressure of O2 [mbar]

P_CO = adsorption(M_CO, PP_CO, Area, Temperature); %CO Adsorption
d_CO = Arrhenius(1.097407, Temperature);     %CO Desorption
h_CO = Arrhenius(0.699167, Temperature);     %CO Hopping

P_O2 = adsorption(M_O2, PP_O2, Area, Temperature); %O Hopping
d_O2 = Arrhenius(2.8683, Temperature);     %O Desorption
h_O = Arrhenius(0.6075,  Temperature);     %O Hopping

k_OCO = Arrhenius(1.08,Temperature);    %LH Reaction
%e_OCO = Arrhenius(0.4175,Temperature); % ER-RXN [Don't use until checked]
%% Master Equations

%{
x(1)=(CO,CO,CO) x(2)=(CO,CO,O) x(3)=(CO,CO,E)   
x(4)=(CO,O,CO)  x(5)=(CO,O,O)  x(6)=(CO,O,E)        
x(7)=(CO,E,CO)  x(8)=(CO,E,O)  x(9)=(CO,E,E)
x(10)=(O,CO,CO) x(11)=(O,CO,O) x(12)=(O,CO,E)
x(13)=(O,O,CO)  x(14)=(O,O,O)  x(15)=(O,O,E)  
x(16)=(O,E,CO)  x(17)=(O,E,O)  x(18)=(O,E,E)
x(19)=(E,CO,CO) x(20)=(E,CO,O) x(21)=(E,CO,E)
x(22)=(E,O,CO)  x(23)=(E,O,O)  x(24)=(E,O,E)
x(25)=(E,E,CO)  x(26)=(E,E,O)  x(27)=(E,E,E)
%}
fx = @(t,x) [-3*d_CO*x(1)+P_CO*(x(3)+x(7)+x(19));                                                           %d/dt P(CO,CO,CO)[1]
             -((2*d_CO)+k_OCO)*x(2)+P_CO*(x(8)+x(20));                                                      %d/dt P(CO,CO,O) [2]
             d_CO*x(1)-(P_CO+2*d_CO+h_CO)*x(3)+h_CO*x(7)+P_CO*(x(9)+x(21));                                 %d/dt P(CO,CO,E) [3]
             -((2*d_CO)+(2*k_OCO))*x(4)+P_CO*(x(6)+x(22));                                                  %d/dt P(CO,O,CO) [4]
             -(d_O2+d_CO+2*k_OCO)*x(5)+P_O2*x(9)+P_CO*x(23);                                                %d/dt P(CO,O,O)  [5]
             d_CO*x(4)-(h_O+P_CO+d_CO+k_OCO)*x(6)+h_O*x(8)+P_CO*x(24);                                      %d/dt P(CO,O,E)  [6]
             d_CO*x(1)+h_CO*x(3)-(P_CO+2*d_CO+2*h_CO)*x(7)+P_CO*(x(9)+x(25))+h_CO*x(19);                    %d/dt P(CO,E,CO) [7]
             d_CO*x(2)+h_O*x(6)-(h_O+P_CO+d_CO+h_CO)*x(8)+h_CO*x(20)+P_CO*x(26);                            %d/dt P(CO,E,O)  [8]
             k_OCO*(x(2)+x(4))+d_CO*(x(3)+x(7))+d_O2*x(5)-(d_O2+2*P_CO+d_CO+h_CO)*x(9)+h_CO*x(21)+P_CO*x(27);%d/dt P(CO,E,E) [9]
             -((2*d_CO)+k_OCO)*x(10)+P_CO*x(12)+P_CO*x(16);                                                 %d/dt P(O,CO,CO)[10]
             -(d_CO+2*k_OCO)*x(11)+P_CO*x(17);                                                              %d/dt P(O,CO,O) [11]
             d_CO*x(10)-(P_CO+d_CO+h_CO+k_OCO)*x(12)+h_CO*x(16)+P_CO*x(18);                                 %d/dt P(O,CO,E) [12]
             -(d_O2+d_CO+k_OCO)*x(13)+P_CO*x(15)+P_O2*x(25);                                                %d/dt P(O,O,CO) [13]
             -(2*d_O2)*x(14)+P_O2*(x(18)+x(26));                                                            %d/dt P(O,O,O)  [14]
             d_CO*x(13)-(d_O2+h_O+P_CO)*x(15)+h_O*x(17)+P_O2*x(27);                                         %d/dt P(O,O,E)  [15] 
             d_CO*x(10)+h_CO*x(12)-(h_O+P_CO+d_CO+h_CO)*x(16)+P_CO*x(18)+h_O*x(22);                         %d/dt P(O,E,CO) [16]
             d_CO*x(11)+h_O*x(15)-((2*h_O)+P_CO)*x(17)+h_O*x(23);                                           %d/dt P(O,E,O)  [17]
             k_OCO*x(11)+d_CO*x(12)+k_OCO*x(13)+d_O2*x(14)+d_CO*x(16)-(P_O2+h_O+(2*P_CO))*x(18)+h_O*x(24);  %d/dt P(O,E,E)  [18]
             d_CO*x(1)+h_CO*x(7)-(P_CO+(2*d_CO)+h_CO)*x(19)+P_CO*x(21)+P_CO*x(25);                          %d/dt P(E,CO,CO)[19]
             d_CO*x(2)+h_CO*x(8)-(P_CO+d_CO+h_CO+k_OCO)*x(20)+P_CO*x(26);                                   %d/dt P(E,CO,O) [20]
             d_CO*x(3)+h_CO*x(9)+d_CO*x(19)-((2*P_CO)+d_CO+(2*h_CO))*x(21)+h_CO*x(25)+P_CO*x(27);           %d/dt P(E,CO,E) [21]
             d_CO*x(4)+h_O*x(16)-(h_O+P_CO+d_CO+k_OCO)*x(22)+P_CO*x(24);                                           %d/dt P(E,O,CO) [22]
             d_CO*x(5)+h_CO*x(17)-(d_O2+h_O+P_CO)*x(23)+P_O2*x(27);                                                %d/dt P(E,O,O)  [23]
             d_CO*x(6)+h_O*x(18)+d_CO*x(22)-((2*h_CO)+(2*P_CO))*x(24)+h_O*x(26);                                   %d/dt P(E,O,E)  [24]
             k_OCO*(x(4)+x(10))+d_CO*(x(7)+x(19))+d_O2*(x(13))+h_CO*x(21)-(P_O2+2*P_CO+d_CO+h_CO)*x(25)+P_CO*x(27);%d/dt P(E,E,CO) [25]
             k_OCO*(x(5)+x(11))+d_CO*(x(8)+x(20))+d_O2*x(14)+h_O*x(24)-x(26)*(P_O2+h_O+2*P_CO);                    %d/dt P(E,E,O)  [26]
             k_OCO*(x(6)+x(12)+x(20)+x(22))+d_CO*(x(9)+x(21)+x(25))+d_O2*(x(15)+x(23))- x(27)*(P_O2*2+P_CO*3);     %d/dt P(E,E,E)  [27]
                                  ];

t_int = [0,0.1];                                                        %Time-Interval [Seconds]
initial_cond = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0];%Initial Conditions from x(1)-x(18) (Full Empty Site Coverage)
[t,xa] = ode45(fx,t_int,initial_cond);                                   %Integrates Equations [1]-[18]
x1 = xa(:,1);   %P(CO,CO,CO)
x2 = xa(:,2);   %P(CO,CO,O)
x3 = xa(:,3);   %P(CO,CO,E)
x4 = xa(:,4);   %P(CO,O,CO)
x5 = xa(:,5);   %P(CO,O,O)
x6 = xa(:,6);   %P(CO,O,E)
x7 = xa(:,7);   %P(CO,E,CO)
x8 = xa(:,8);   %P(CO,E,O)
x9 = xa(:,9);   %P(CO,E,E)
x10 = xa(:,10); %P(O,CO,CO)
x11 = xa(:,11); %P(O,CO,O)
x12 = xa(:,12); %P(O,CO,E)
x13 = xa(:,13); %P(O,O,CO)
x14 = xa(:,14); %P(O,O,O)
x15 = xa(:,15); %P(O,O,E)
x16 = xa(:,16); %P(O,E,CO)
x17 = xa(:,17); %P(O,E,O)
x18 = xa(:,18); %P(O,E,E)
x19 = xa(:,19); %P(E,CO,CO)
x20 = xa(:,20); %P(E,CO,O)
x21 = xa(:,21); %P(E,CO,E)
x22 = xa(:,22); %P(E,O,CO)
x23 = xa(:,23); %P(E,O,O)
x24 = xa(:,24); %P(E,O,E)
x25 = xa(:,25); %P(E,E,CO)
x26 = xa(:,26); %P(E,E,O)
x27 = xa(:,27); %P(E,E,E)

O1 = x10+x11+x12+x13+x14+x15+x16+x17+x18;%O1 Total 
O2 = x4+x5+x6+x13+x14+x15+x22+x23+x24;   %O2 Total 
O3 = x2+x5+x8+x11+x14+x17+x20+x23+x26;   %O3 Total 

CO1 = x1+x2+x3+x4+x5+x6+x7+x8+x9;        %CO1 Total
CO2 = x1+x2+x3++x10+x11+x12+x19+x20+x21; %CO2 Total
CO3 = x1+x4+x7+x10+x13+x16+x19+x22+x25;  %CO3 Total

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
xlim([-0.00001,0.015]);
ylim([-0.1,1.05]);
lgd1 = legend('$P(O_1)$','$P(CO_1)$','$P(E_1)$','Interpreter','latex','Location','northeast','Fontsize',20);
title("Probability Vs Time",'Fontsize',20,'Fontname','times');
grid(ax1,'on')

ax2 = nexttile; %Site 2
plot(t,Site_2,'--','Linewidth',2);
xlim([-0.00001,0.015]);
ylim([-0.1,1.05]);
lgd2 = legend('$P(O_2)$','$P(CO_2)$','$P(E_2)$','Interpreter','latex','Location','northeast','Fontsize',20);
grid('on')

ax3 = nexttile; %Site 3
plot(t,Site_3,'--','Linewidth',2);
xlim([-0.00001,0.015]);
ylim([-0.1,1.05]);
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

function a = Arrhenius(Ea,T)   %(Activation Energy, Temp.)
    Nu = 10^12.5;              %Prefactor range is 10^12-10^13 1/second
    kB = 8.617333262145*10^-5; %Boltzman Constant in eV/T
    a = Nu*exp(-Ea/(kB*T));    %Arrhenius Equation
end

function f = adsorption(MM, PP, A, T) %(Molar Mass, Partial-P, Area, Temp)
    kB = 8.617333262145*10^-5;        %Boltzman Constant in eV/T
    f = (A*PP)/(2*pi*MM*kB*T);        %Adsorption Rate Equation
end
