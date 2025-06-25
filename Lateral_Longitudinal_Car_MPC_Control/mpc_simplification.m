function [Hdb, Fdbt, Cdb, Adc, G, ht] = mpc_simplification(Ad, Bd, Cd, Dd, hz, x_aug_t, du)
    
    [A_aug, B_aug, C_aug, ~] = generate_augmented_matrices(Ad, Bd, Cd, Dd);

    % Get simulation constants
    sim_data = init_constants();

    Q = sim_data{'Q'};
    S = sim_data{'S'};
    R = sim_data{'R'};
    Cf = sim_data{'Cf'};
    g = sim_data{'g'};
    m = sim_data{'m'};
    mju = sim_data{'mu'};
    lf = sim_data{'lf'};
    inputs = sim_data{'inputs'};

    % Inputs Constraints
    d_delta_max = pi/300;
    d_a_max = 0.1;
    d_delta_min = -pi/300;
    d_a_min = -0.1;

    ub_global = zeros(1, inputs*hz);
    lb_global = zeros(1, inputs*hz);

    for i=1:inputs*hz
        if mod(i-1,2)==0
            ub_global(i) = d_delta_max;
            lb_global(i) = -d_delta_min;
        else
            ub_global(i) = d_a_max;
            lb_global(i) = -d_a_min;
        end
    end

    ub_global = ub_global(1:inputs*hz);
    lb_global = lb_global(1:inputs*hz);
    ublb_global = [ub_global, lb_global];

    I_global = eye(inputs*hz);
    I_global_negative = -I_global;
    I_mega_global = [I_global; I_global_negative];

    y_asterisk_max_global = [];
    y_asterisk_min_global = [];

    C_asterisk = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1];
    C_asterisk_global = zeros(size(C_asterisk, 1)*hz, size(C_asterisk, 2)*hz);

    CQC = C_aug'*Q*C_aug;
    CSC = C_aug'*S*C_aug;
    QC = Q*C_aug;
    SC = S*C_aug;

    Qdb = zeros(size(CQC, 1)*hz, size(CQC, 2)*hz);
    Tdb = zeros(size(QC, 1)*hz, size(QC, 2)*hz);
    Rdb = zeros(size(R, 1)*hz, size(R, 2)*hz);
    Cdb = zeros(size(B_aug, 1)*hz, size(B_aug, 2)*hz);
    Adc = zeros(size(A_aug, 1)*hz, size(A_aug, 2));

    % Filling above matrices for LPV system
    A_product = A_aug;
    state_predicted_aug = x_aug_t;
    A_aug_size = size(A_aug);
    B_aug_size = size(B_aug);
    A_aug_collection = zeros(A_aug_size(1), A_aug_size(2), hz);
    B_aug_collection = zeros(B_aug_size(1), B_aug_size(2), hz);

    for i= 1:hz
        if i==hz
            Qdb(1+size(CSC,1)*(i-1):size(CSC, 1)*i, 1+size(CSC, 2)*(i-1):size(CSC, 2)*i) = CSC;
            Tdb(1+size(SC,1)*(i-1):size(SC, 1)*i, 1+size(SC, 2)*(i-1):size(SC, 2)*i) = SC;
        else
            Qdb(1+size(CQC,1)*(i-1):size(CQC, 1)*i, 1+size(CQC, 2)*(i-1):size(CQC, 2)*i) = CQC;
            Tdb(1+size(QC,1)*(i-1):size(QC, 1)*i, 1+size(QC, 2)*(i-1):size(QC, 2)*i) = QC;
        end

        Rdb(1+size(R, 1)*(i-1):size(R, 1)*i, 1+size(R, 2)*(i-1):size(R, 2)*i) = R;

        Adc(1+A_aug_size(1)*(i-1):A_aug_size(1)*i, :) = A_product;
        A_aug_collection(:,:,i) = A_aug;
        B_aug_collection(:,:,i) = B_aug;

        % Forming state contraints
        x_dot_max = 30;
        if 0.17*state_predicted_aug(1) < 3
            y_dot_max = 0.17*state_predicted_aug(1);
        else
            y_dot_max = 3;
        end
        delta_max = pi/6;
        Fyf=Cf*(state_predicted_aug(7)-state_predicted_aug(2)/state_predicted_aug(1)-lf*state_predicted_aug(4)/state_predicted_aug(1));
        a_max=1+(Fyf*sin(state_predicted_aug(7))+mju*m*g)/m-state_predicted_aug(4)*state_predicted_aug(2);

        x_dot_min = 1;
        if -0.17*state_predicted_aug(1) > -3
            y_dot_min = -0.17*state_predicted_aug(1);
        else
            y_dot_min = -3;
        end
        delta_min = -pi/6;
        a_min=-4+(Fyf*sin(state_predicted_aug(7))+mju*m*g)/m-state_predicted_aug(4)*state_predicted_aug(2);

        y_asterisk_max = [x_dot_max, y_dot_max, delta_max, a_max];
        y_asterisk_min = [x_dot_min, y_dot_min, delta_min, a_min];

        y_asterisk_max_global = [y_asterisk_max_global, y_asterisk_max];
        y_asterisk_min_global = [y_asterisk_min_global, y_asterisk_min];
        
        C_asterisk_global(1+size(C_asterisk, 1)*(i-1):size(C_asterisk, 1)*i, 1+size(C_asterisk, 2)*(i-1):size(C_asterisk, 2)*i) = C_asterisk;

        if i<hz
            du1 = du(inputs*i);
            du2 = du(inputs*(i+1)-1);
            state_predicted_aug = A_aug*state_predicted_aug + B_aug*[du1, du2]';
            state_predicted = state_predicted_aug(1:6)';
            delta_predicted = state_predicted_aug(7);
            a_predicted = state_predicted_aug(8);
            [Ad, Bd, Cd, Dd] = state_space(state_predicted, delta_predicted);
            [A_aug, B_aug, C_aug, ~] = generate_augmented_matrices(Ad, Bd, Cd, Dd);
            A_product = A_aug*A_product;
        end
    end
    
    for i=1:hz
        for j=1:hz
            if j<=i
                AB_product = eye(A_aug_size(1));
                for ii=i:-1:j
                    if ii>j
                        AB_product = AB_product*A_aug_collection(:,:,ii);
                    else
                        AB_product = AB_product*B_aug_collection(:,:,ii);
                    end
                end
                Cdb(1+B_aug_size(1)*(i-1):B_aug_size(1)*i, 1+B_aug_size(2)*(j-1):B_aug_size(2)*j) = AB_product;
            end
        end
    end

    Cdb_constraints = C_asterisk_global*Cdb;
    Cdb_constraints_negative = -Cdb_constraints;
    Cdb_constraints_global = [Cdb_constraints; Cdb_constraints_negative];

    Adc_constraints = C_asterisk_global*Adc;
    Adc_constraints_x0 = (Adc_constraints*x_aug_t)';

    y_max_Adc_difference=y_asterisk_max_global-Adc_constraints_x0;
    y_min_Adc_difference=-y_asterisk_min_global+Adc_constraints_x0;
    y_Adc_difference_global=[y_max_Adc_difference,y_min_Adc_difference];
    
    G=[I_mega_global;Cdb_constraints_global];
    ht=[ublb_global,y_Adc_difference_global];
    
    Hdb=Cdb'*Qdb*Cdb+Rdb;
    Fdbt=[Adc'*Qdb*Cdb;-Tdb*Cdb];

end