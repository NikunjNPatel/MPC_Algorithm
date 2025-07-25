function [Ad, Bd, Cd, Dd] = state_space(states, delta)

    % Simulation constants
    sim_data = init_constants();
    g = sim_data{'g'};
    m = sim_data{'m'};
    Iz = sim_data{'Iz'};
    Cf = sim_data{'Cf'};
    Cr = sim_data{'Cr'};
    lf = sim_data{'lf'};
    lr = sim_data{'lr'};
    Ts = sim_data{'Ts'};
    mu = sim_data{'mu'};

    % states
    x_dot = states(1);
    y_dot = states(2);
    psi = states(3);

    % Calculate continuous state matrix A, B, C, D elements
    A11=-mu*g/x_dot;
    A12=Cf*sin(delta)/(m*x_dot);
    A14=Cf*lf*sin(delta)/(m*x_dot)+y_dot;
    A22=-(Cr+Cf*cos(delta))/(m*x_dot);
    A24=-(Cf*lf*cos(delta)-Cr*lr)/(m*x_dot)-x_dot;
    A34=1;
    A42=-(Cf*lf*cos(delta)-lr*Cr)/(Iz*x_dot);
    A44=-(Cf*lf.^2*cos(delta)+lr.^2*Cr)/(Iz*x_dot);
    A51=cos(psi);
    A52=-sin(psi);
    A61=sin(psi);
    A62=cos(psi);
    
    B11=-1/m*sin(delta)*Cf;
    B12=1;
    B21=1/m*cos(delta)*Cf;
    B41=1/Iz*cos(delta)*Cf*lf;
    
    A=[A11, A12, 0, A14, 0, 0; ...
       0, A22, 0, A24, 0, 0; ...
       0, 0, 0, A34, 0, 0; ...
       0, A42, 0, A44, 0, 0; ...
       A51, A52, 0, 0, 0, 0; ...
       A61, A62, 0, 0, 0, 0];
    B=[B11, B12; B21, 0; 0, 0; B41, 0; 0, 0; 0, 0];
    C=[1, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 1];
    D=[0, 0; 0, 0; 0, 0; 0, 0];
    
    % Discretise the system (forward Euler)
    Ad=eye(size(A, 2)) + Ts*A;
    Bd=Ts*B;
    Cd=C;
    Dd=D;

end