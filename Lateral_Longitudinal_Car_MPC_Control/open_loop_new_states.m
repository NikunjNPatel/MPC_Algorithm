function dStates = open_loop_new_states(t, old_states, U)
    
    sim_data = init_constants();
    g = sim_data{'g'};
    m = sim_data{'m'};
    Iz = sim_data{'Iz'};
    Cf = sim_data{'Cf'};
    Cr = sim_data{'Cr'};
    lf = sim_data{'lf'};
    lr = sim_data{'lr'};
    mju = sim_data{'mu'};

    x_dot = old_states(1);
    y_dot = old_states(2);
    psi = old_states(3);
    psi_dot = old_states(4);

    delta = U(1);
    a = U(2);

    Fyf=Cf*(delta-y_dot/x_dot-lf*psi_dot/x_dot);
    Fyr=Cr*(-y_dot/x_dot+lr*psi_dot/x_dot);
   
    % The nonlinear equation describing the dynamics of the drone
    dStates(1,1)=a+(-Fyf*sin(delta)-mju*m*g)/m+psi_dot*y_dot;
    dStates(2,1)=(Fyf*cos(delta)+Fyr)/m-psi_dot*x_dot;
    dStates(3,1)=psi_dot;
    dStates(4,1)=(Fyf*lf*cos(delta)-Fyr*lr)/Iz;
    dStates(5,1)=x_dot*cos(psi)-y_dot*sin(psi);
    dStates(6,1)=x_dot*sin(psi)+y_dot*cos(psi);
end