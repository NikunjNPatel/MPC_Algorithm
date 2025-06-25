function sim_data = init_constants()
    
    % Create empty dictionary to store all constants 
    sim_data = dictionary();
    sim_data{'g'} = 9.81;
    sim_data{'m'} = 1500;
    sim_data{'Iz'} = 3000;
    sim_data{'Cf'} = 38000;
    sim_data{'Cr'} = 66000;
    sim_data{'rho'} = 1.225;
    sim_data{'lf'} = 2;
    sim_data{'lr'} =  3;
    sim_data{'Ts'} =  0.02;
    sim_data{'mu'} = 0.02;
    sim_data{'outputs'} = 4;
    sim_data{'inputs'} = 2;
    sim_data{'trajectory'} = 3; % Choose trajectory from 1, 2, 3

    if sim_data{'trajectory'} == 1
        sim_data{'hz'} = 10;
        sim_data{'time_length'} = 60;
        sim_data{'Q'} = [1 0 0 0;0 200 0 0;0 0 50 0;0 0 0 50];
        sim_data{'S'} = [1 0 0 0;0 200 0 0;0 0 50 0;0 0 0 50];
        sim_data{'R'} = [100 0;0 1];
    elseif sim_data{'trajectory'} == 2
        sim_data{'hz'} = 10;
        sim_data{'time_length'} = 140;
        sim_data{'Q'} = [10000 0 0 0;0 10000 0 0;0 0 100 0;0 0 0 100];
        sim_data{'S'} = [10000 0 0 0;0 10000 0 0;0 0 100 0;0 0 0 100];
        sim_data{'R'} = [100 0;0 1];
    else
        sim_data{'hz'} = 10;
        sim_data{'Q'} = [100 0 0 0;0 100000 0 0;0 0 5000 0;0 0 0 5000];
        sim_data{'S'} = [100 0 0 0;0 100000 0 0;0 0 5000 0;0 0 0 5000];
        sim_data{'R'} = [10000 0;0 100];
        first_section = 14;
        other_section = 14;
        sim_data{'time_length'} = first_section + other_section*10;
        delay=zeros(1,12);
        for x = 2:length(delay)
            delay(x)=first_section+(x-2)*other_section;
        end
        sim_data{'delay'} = delay;
    end
end

