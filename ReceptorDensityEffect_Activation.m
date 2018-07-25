clear

%% Define scenario
% promoter P1 produces a activator A that activates promoter P2
% binding of inducer I to A with cooperativity n will render A active
% P2 drives expression of the output gene GFP
%
% since only steady state solution is considered here, I assume degradation
% term is incoporated as a non-identifiable parameter that is already
% incoporated into the parameters of basal and max tx rates

%parameters
Vmax_P1 = 100;              %tx rate for repressor
alpha_P2 = 1;       %basal tx rate for GFP
Vmax_P2 = 10000;          %max tx rate for GFP
Ka = 1000;                %K_D between P2 and A
Ki = 100;               %K_D between A and I
n = 2;                  %hill coefficient for binding between A and I
d_1 = 0.5;              %degradation rate for A
d_2 = 0.5;              %degradation rate for GFP

%passing parameters into parameter vector for ODE fun
para = zeros(8,1);
para(1) = Vmax_P1;
para(2) = alpha_P2;
para(3) = Vmax_P2;
para(4) = Ka;
para(5) = Ki;
para(6) = n;
para(7) = d_1;
para(8) = d_2;

%input conditions
%I = 200;                %inducer concentration to be tested

inducerRange = [linspace(0,100,101) linspace(0,1e4,101)];    %inducer range (fixed and unchanged for one expt)
inducerRange = inducerRange';

VmaxP1_range = [200 400 800 1600 6400 12800];

rc = true;  %decides if resource competition happens, note that rc only controls which OED to execute and graph appearance
%% Solution Solving

%set simulation time
t_len=10;

%set initial condition
species_t0 = zeros(2,1);
species_t0(1) = 0;  %initial R concentration
species_t0(2) = 0;  %initial GFP concentration

GFP_final = zeros(length(inducerRange),1);

for j = 1:length(VmaxP1_range)
    para(1) = VmaxP1_range(j);
    
    for i = 1:length(inducerRange)
        I = inducerRange(i);
        if rc == false
            [t_output,species_output] = ode45(@(t,y) activationOED(t,y,para,I), [0,t_len], species_t0, []);
        else
            [t_output,species_output] = ode45(@(t,y) activationOED_RC(t,y,para,I), [0,t_len], species_t0, []);
        end
        
        R_output = species_output(:,1);
        GFP_output = species_output(:,2);
        
        GFP_final(i,j) = GFP_output(end);
    end
end

%% Plot Graphs

displayTimeProfile = false;

if displayTimeProfile == true
fig1 = figure(1);
set(fig1,'Name','Time Profile of Species');
graph1 = plot(t_output, R_output, 'b', t_output, GFP_output, 'r');
title('Time Profile of Species');
xlabel('time (au)');
ylabel('species concentration (au)');
set(graph1,'LineWidth',3);
legend('Repressor', 'GFP');
end

displayResponseCurve = false;

if displayResponseCurve == true
fig2 = figure(2);
set(fig2,'Name','Response Curve');
graph2 = semilogx(inducerRange, GFP_final);
title('Response Curve');
xlabel('inducer concentration (au)');
ylabel('GFP concentration (au)');
set(graph2,'LineWidth',3);
%axis([10,1e4,200,2000]);
end

displayMultipleResponseCurve = true;
if displayMultipleResponseCurve == true;
%Figure to illustrate how response curve changes with Vmax_P1
fig3 = figure(3);
set(fig3,'Name','Response Curves under different P1');

for j = 1:length(VmaxP1_range);

    graph3 = semilogx(inducerRange, GFP_final);

end
title('Activation: Response curves under different P_c');
xlabel('inducer concentration (au)');
ylabel('GFP concentration (au)');
set(graph3,'LineWidth',3);
lgn = legend( strcat('VmaxP1 = ', num2str( VmaxP1_range' ) ) ); 
lgn.Location = 'northwest';
refline(0,200);
if rc == false 
    axis([1,1e4,0,20000]);
else
    axis([1,1e4,0,8000]);
end
end