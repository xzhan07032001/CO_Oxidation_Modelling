%Name:Xin (Jason) Zhang, Creation Date:9/28/2021, Updated Since:9/30/2021 Version:R2021a

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
-Functions (Rate and Phase Functions)
%}
%% Parameters
Temperature = 300;  %Kelvin
M_O2 = 32;          %Molar mass of O2
M_CO = 28;          %Molar mass of CO
Area = 5^2;         %Area per unit cell [Angstrom^2]

d_CO = Arrhenius(1.097407, Temperature);           %CO Desorption
h_CO = Arrhenius(0.699167, Temperature);           %CO Hopping

d_O2 = Arrhenius(2.8683, Temperature);             %O Desorption
h_O = Arrhenius(0.6075,  Temperature);             %O Hopping

k_OCO = Arrhenius(1.08,Temperature);               %LH Reaction original is 1.08
%% Master Equation Initialization and TOF Calculations
PP_CO_MATRIX = logspace(1,3,200); %Partial Pressure Range of CO [mbar]
PP_O2_MATRIX= logspace(1,3,200); %Partial Pressure Range of O2 [mbar]



t_int = [0,0.2];                       %Time-Interval [Seconds] original is .005
initial_cond = [0,0,0,0,0,0];           %Initial Conditions from x(1)-x(18)
Phase_Matrix1 = zeros(200);                 %Initialize Phase Zero Matrices [Make sure Pressure matrix dimensions match!!!)
Phase_Matrix2 = zeros(200);
Phase_Matrix3 = zeros(200);

for i = 1:length(PP_CO_MATRIX)
    PP_CO = PP_CO_MATRIX(i);
    P_CO = adsorption(M_CO, PP_CO, Area, Temperature); %CO Adsorption
    disp(i);
    for j = 1:length(PP_O2_MATRIX)
        PP_O2 = PP_O2_MATRIX(j);
        P_O2 = adsorption(M_O2, PP_O2, Area, Temperature); %O Adsorption
        
        fx = @(t,x) [h_O*(1-x(1)-x(3))*x(2)+P_O2*(1-x(1)-x(3))*(1-x(2)-x(4))-d_O2*x(1)*x(2)-h_O*x(1)*(1-x(2)-x(4))-k_OCO*x(1)*x(4);                                       %d/dt=O1 [1] COMPLETE
             h_O*(1-x(2)-x(4))*(x(1)+x(5))+P_O2*2*(1-x(2)-x(4))*((1-x(1)-x(3))+(1-x(5)-x(6)))-d_O2*x(1)*(x(2)+x(6))-h_O*x(2)*(1-x(5)-x(6))-k_OCO*x(2)*(x(3)+x(6));%d/dt=O2 [2] COMPLETE 
             P_CO*(1-x(1)-x(3))+h_CO*(1-x(1)-x(3))*x(4)-d_CO*x(3)-k_OCO*x(3)*x(2)-h_CO*x(3)*(1-x(2)-x(4));                                                        %d/dt=CO1[3] COMPLETE
             P_CO*(1-x(2)-x(4))+h_CO*(1-x(2)-x(4))*(x(3)+x(4))-d_CO*x(4)-k_OCO*2*x(4)*x(1)-h_CO*x(3)*((1-x(1)-x(3))+(1-x(5)-x(6)));                               %d/dt=CO2[4]
             h_O*(1-x(5)-x(6))*x(2)+P_O2*(1-x(2)-x(4))*(1-x(5)-x(6))-d_O2*x(2)*x(5)-h_O*x(5)*(1-x(2)-x(4))-k_OCO*x(5)*x(4);                                       %d/dt=O3 [5] CHECK
             P_CO*(1-x(5)-x(6))+h_CO*((1-x(5)-x(6))*x(4)+(1-x(2)-x(4))*x(6))-d_CO*x(6)-k_OCO*x(6)*x(2);                                                           %d/dt=CO3 [6]
                                  ];
        [t,xa] = ode15s(fx,t_int,initial_cond); %Integrates Equations [1]-[18]
        x1 = xa(end,1);   %x(1) P(O1)
        x2 = xa(end,2);   %x(2) P(O2)
        x3 = xa(end,3);   %x(3) P(CO1)
        x4 = xa(end,4);   %x(4) P(CO2)
        x5 = xa(end,5);   %x(5) P(O3)
        x6 = xa(end,6);   %x(6) P(CO3)
        % "1" is O, "2" is CO, and "3" is Empty (E), "0" is the equivalent
        % to ? in Wang et al
        
        if x1 > 0.8 
            Phase_Matrix1(i,j) = 1;
        elseif x3 > 0.8
            Phase_Matrix1(i,j) = 2;
        elseif (1-x1-x3) > 0.8
            Phase_Matrix1(i,j) = 3;
        end
        
        if x2 > 0.8 
            Phase_Matrix2(i,j) = 1;
        elseif x4 > 0.8
            Phase_Matrix2(i,j) = 2;
        elseif (1-x2-x4) > 0.8
            Phase_Matrix2(i,j) = 3;
        end
        
        if x5 > 0.8 
            Phase_Matrix3(i,j) = 1;
        elseif x6 > 0.8
            Phase_Matrix3(i,j) = 2;
        elseif (1-x5-x6) > 0.8
            Phase_Matrix3(i,j) = 3;
            
        end
    end
end


%% Graphing
% Implement Later Implement above to below in particular 1 = O2, 2 = CO, 3
% = Empty, Separate 3 states into individual sites
figure(1);
tile = tiledlayout(3,1);

ax1 = nexttile; %Site 1
[X,Y] = meshgrid(PP_CO_MATRIX,PP_O2_MATRIX);
loglog(X,Y);
surface(X,Y,Phase_Matrix1);
shading flat
colorbar
colormap winter

ax2 = nexttile; %Site 2
[X,Y] = meshgrid(PP_CO_MATRIX,PP_O2_MATRIX);
loglog(X,Y);
surface(X,Y,Phase_Matrix2);
shading flat
colorbar
colormap winter

ax3 = nexttile; %Site 3
[X,Y] = meshgrid(PP_CO_MATRIX,PP_O2_MATRIX);
loglog(X,Y);
surface(X,Y,Phase_Matrix3);
shading flat
colorbar
colormap winter

%Graph Formatting
linkaxes([ax1,ax2,ax3],'xy');
ylabel(tile,"Pressure O2 [mbar]",'Fontsize',25,'fontname','times');
xlabel(tile,"Pressure CO [mbar]",'Fontsize',25,'fontname','times');
ax1.FontSize = 20;
ax2.FontSize = 20;
ax3.FontSize = 20;
ax1.FontName = 'times';
ax2.FontName = 'times';
ax3.FontName = 'times';

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
