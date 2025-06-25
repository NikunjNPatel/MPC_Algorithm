clear all;
close all;
clc;

% Generate simulation data used in model
sim_data = init_constants();

% Load the values needed
Ts = sim_data{'Ts'};
no_of_outputs = sim_data{'outputs'};
hz = sim_data{'hz'};
no_of_inputs = sim_data{'inputs'};
trajectory = sim_data{'trajectory'};

% Time vector for the simulation
t = 0:Ts:sim_data{'time_length'};
sim_length = length(t);

% Generate reference signals using generate trajectory function
plot_trajectory = true;
[xDot_ref, yDot_ref, psi_ref, X_ref, Y_ref] = generate_trajectory(t, plot_trajectory);

% Concatenate reference signal into single vector
refSignals = zeros(1, sim_length*no_of_outputs);
idx = 1;
for i = 1:no_of_outputs:length(refSignals)
    refSignals(i) = xDot_ref(idx);
    refSignals(i+1) = psi_ref(idx);
    refSignals(i+2) = X_ref(idx);
    refSignals(i+3) = Y_ref(idx);
    idx = idx + 1;
end

% Initial state and storage arrays
x_dot = xDot_ref(1);
y_dot = yDot_ref(1);
psi = psi_ref(1);
psi_dot = 0;
X = X_ref(1);
Y = Y_ref(1);

states = [x_dot, y_dot, psi, psi_dot, X, Y];
storeStates = zeros(sim_length, length(states));
storeStates(1,:) = states;

% Initial acceleration and storage arrays
x_dot_dot = 0;
y_dot_dot = 0;
psi_dot_dot = 0;

acceleration = [x_dot_dot, y_dot_dot, psi_dot_dot];
storeAcceleration = zeros(sim_length, length(acceleration));

% Initial inputs and storage arrays
delta = 0;  % Steering angle
net_acceleration = 0;
storeInputs = zeros(sim_length, no_of_inputs);
storeInputs(1,1) = delta;
storeInputs(1,2) = net_acceleration;

deltaInputs = zeros(no_of_inputs*hz, 1); % To store delta inputs obtained as result from MPC

% Main loop
k = 1;
for i=1:sim_length-1
    
    % Discrete state matrics of LPV system
    [Ad, Bd, Cd, Dd] = state_space(states, delta);

    % Create augmented state vector including inputs
    x_aug_t = [states'; delta; net_acceleration];

    k=k+no_of_outputs;
    if k+no_of_outputs*hz-1 <= length(refSignals)
        r=refSignals(k:k+no_of_outputs*hz-1);
    else
        r=refSignals(k:length(refSignals));
        hz=hz-1;
    end
    
    % Generate simplification matrices for the cost function
    [Hdb,Fdbt,Cdb,Adc,G,ht] = mpc_simplification(Ad,Bd,Cd,Dd,hz,x_aug_t,deltaInputs);
    ft=[x_aug_t',r]*Fdbt;

    %% Calling the optimizer (quadprog)
    
    % Cost function in quadprog: min(du)*1/2*du'Hdb*du+f'du
    % Hdb must be positive definite for the problem to have finite minimum.
    options = optimoptions('quadprog','Display', 'off','LinearSolver','dense');

    [deltaInputs,fval]=quadprog(Hdb,ft,G,ht,[],[],[],[],[],options);
    
    if isempty(deltaInputs)
        disp('The solver could not find the solution');
    end
    
    delta=delta+deltaInputs(1);
    net_acceleration=net_acceleration+deltaInputs(2);
    
    storeInputs(i+1,1)=delta;
    storeInputs(i+1,2)=net_acceleration;
    
    % Simulate the new states
    time_interval=(Ts)/30;
    T = (Ts)*(i-1):time_interval:Ts*(i-1)+(Ts);
    [T,x]=ode45(@(t,x) open_loop_new_states(t,x,[delta,net_acceleration]),T,states);
    
    states=x(end,:);
    storeStates(i+1,:)=states;
    
    % Accelerations
    x_dot_dot=(x(end,1)-x(end-1,1))/time_interval;
    y_dot_dot=(x(end,2)-x(end-1,2))/time_interval;
    psi_dot_dot=(x(end,4)-x(end-1,4))/time_interval;
    
    acceleration=[x_dot_dot,y_dot_dot,psi_dot_dot];
    storeAcceleration(i+1,:)=acceleration;
    
    if mod(i,500)==0
        fprintf('Progress (%f) \n', i/sim_length*100)
    end
end

figure(1);
plot(X_ref, Y_ref, '--b', 'LineWidth', 2)
hold on
plot(storeStates(:,5), storeStates(:,6), 'r', 'LineWidth', 1)
grid on
xlabel('x - position [m]')
ylabel('y - position [m]')
