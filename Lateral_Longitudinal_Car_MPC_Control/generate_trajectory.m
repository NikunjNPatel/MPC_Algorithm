function [x_dot_body, y_dot_body, psiInt, X, Y] = generate_trajectory(t_vector, plot_flag)

    if nargin < 2
        plot_flag = false;
    end
    
    sim_data = init_constants();
    Ts = sim_data{'Ts'};
    trajectory = sim_data{'trajectory'};

    if trajectory == 1
        X = 15*t_vector;
        Y = 750/900.^2*X.^2+250;
    elseif trajectory == 2
        X1=15*t_vector(1:40/Ts);
        Y1=50*sin(2*pi*0.75/40*t_vector(1:40/Ts))+250;
    
        X2=300*cos(2*pi*0.5/60*(t_vector(40/Ts+1:100/Ts)-40)-pi/2)+600;
        Y2=300*sin(2*pi*0.5/60*(t_vector(40/Ts+1:100/Ts)-40)-pi/2)+500;
    
        X3=600-15*(t_vector(100/Ts+1:140/Ts+1)-100);
        Y3=50*cos(2*pi*0.75/40*(t_vector(100/Ts+1:140/Ts+1)-100))+750;
    
        X=[X1,X2,X3];
        Y=[Y1,Y2,Y3];
    else
        delay = sim_data{'delay'};

        % X & Y levels
        f_x=[0,60,110,140,160,110,40,10,40,70,110,150];
        f_y=[40,20,20,60,100,140,140,80,60,60,90,90];
        
        % X & Y derivatives
        f_x_dot=[1,2,1,1,0,-1,-1,0,1,1,1,1]*3;
        f_y_dot=[0,0,0,1,1,0,0,-1,0,0,0,0]*3;
        
        X=[];
        Y=[];
        for i = 1:length(delay)-1
           % Extract the time elements for each section separately
           if not(i==length(delay)-1)
               t_temp=t_vector(delay(i)/Ts+1:delay(i+1)/Ts);
           else
               t_temp=t_vector(delay(i)/Ts+1:delay(i+1)/Ts+1);
           end
    
           % Generate data for a subtrajectory
           M=[1,t_temp(1),t_temp(1).^2,t_temp(1).^3; ...
               1,t_temp(end),t_temp(end).^2,t_temp(end).^3; ...
               0,1,2*t_temp(1),3*t_temp(1).^2; ...
               0,1,2*t_temp(end),3*t_temp(end).^2];
           
           c_x=[f_x(i);f_x(i+1)-f_x_dot(i+1)*Ts;f_x_dot(i);f_x_dot(i+1)];
           c_y=[f_y(i);f_y(i+1)-f_y_dot(i+1)*Ts;f_y_dot(i);f_y_dot(i+1)];
           
           a_x=inv(M)*c_x;
           a_y=inv(M)*c_y;
           
           % Compute X and Y values
           X_temp=a_x(1)+a_x(2)*t_temp+a_x(3)*t_temp.^2+a_x(4)*t_temp.^3;
           Y_temp=a_y(1)+a_y(2)*t_temp+a_y(3)*t_temp.^2+a_y(4)*t_temp.^3;
           
           X=[X,X_temp];
           Y=[Y,Y_temp];
        end
    end
    
    % Generate yaw angle and body frame xDot, yDot
    X = round(X, 8);
    Y = round(Y, 8);
    
    dX=X(2:end)-X(1:end-1);
    dY=Y(2:end)-Y(1:end-1);

    XDot = dX./Ts; YDot =  dY./Ts;
    XDot = [XDot(1), XDot];
    YDot = [YDot(1), YDot];
    
    XDot = round(XDot, 8); YDot = round(YDot, 8);
    
    
    psi=zeros(1,length(X));
    psiInt=psi;
    psi(1)=atan2(dY(1),dX(1));
    psi(2:end)=atan2(dY(:),dX(:));
    dpsi=psi(2:end)-psi(1:end-1);
    
    psiInt(1)=psi(1);
    for i = 2:length(psiInt)
        if dpsi(i-1)<-pi
            psiInt(i)=psiInt(i-1)+(dpsi(i-1)+2*pi);
        elseif dpsi(i-1)>pi
            psiInt(i)=psiInt(i-1)+(dpsi(i-1)-2*pi);
        else
            psiInt(i)=psiInt(i-1)+dpsi(i-1);
        end
    end

    x_dot_body=cos(psiInt).*XDot+sin(psiInt).*YDot;
    y_dot_body=-sin(psiInt).*XDot+cos(psiInt).*YDot;
    
    psiInt=round(psiInt,8);
    x_dot_body=round(x_dot_body,8);
    y_dot_body=round(y_dot_body,8);

    if trajectory==2
        % Modify x_dot_body (only for trajectory 2)
        x_dot_body_temp=x_dot_body(40/Ts+1:60/Ts);
        x_dot_body_temp=linspace(x_dot_body_temp(1),16.065,1000);
        x_dot_body(40/Ts+1:60/Ts)=x_dot_body_temp;
    
        x_dot_body_temp=x_dot_body(80/Ts+3:100/Ts+2);
        x_dot_body_temp=linspace(16.065,x_dot_body_temp(end),1000);
        x_dot_body(80/Ts+3:100/Ts+2)=x_dot_body_temp;
    
        x_dot_body(60/Ts+1:80/Ts+2)=16.065;

    end

    if plot_flag
        plot(X, Y)
        xlabel('X - Position [m]')
        ylabel('Y - Position [m]')
        title(sprintf('Trajectory - %d', trajectory))
    end
end