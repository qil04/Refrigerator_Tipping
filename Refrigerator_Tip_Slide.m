close all;
clear all;
clc;


% Inputs: Inputs include the Mass of the block, angle of pull, force of pull, 
%         coefficients of kinetic and static friction.
% 
% Output: Output includes the acceleration (ax and ay) of the block.
% This file generates a series of condition check problems such as tip/slide 
%       that determine the refrigerator by a given value. The user have to
%       put their assumptions into the output dialog box. If the
%       assumption wroks, it will generate the displacement and velocity
%       with Euler's method.
%File created by Chris Li 2021/3/24

%Part1 Assign values to the system (can be changed)
m = 10;          %kg
angle_pull = pi/3;  %rad
P = 100;            %N
uk = 0.5;           %dimensionless
us = 0.8;           %dimensionless
B = 0.6;            %width m
H = 1.8;            %height of refrigerator m
h = 1.5;            %height of pulling force m
I = (1/12)*m*(B^2+H^2); %Mass Moment of inertia
g = 9.81;             %gravity acceleration

syms Ff Na Nb ax ay alpha
%Define Newtonia equation (Sum Fx, Sum Fy, Sum moment)
equ_x = P*cos(angle_pull) - Ff == m*ax;
equ_y = P*sin(angle_pull) + Na + Nb - m*g == m*ay;
equ_m = -P*sin(angle_pull)*(B/2) - P*cos(angle_pull)*(h-H/2) - Na*(B/2) + Nb*(B/2) - (H/2)*Ff == I*alpha;
sol_nm = solve([equ_x, equ_y, equ_m], [Ff, Na, Nb]);
Ff_general = sol_nm.Ff;
Na_general = sol_nm.Na;
Nb_general = sol_nm.Nb;

%For users provides information
% motion_type = input("Enter a case to check ");
motion_type = 'PointA Tip';
%Switch method consist the whole 6 situations:      Checking Condition
%   1. No motion in statis equalibrium          Na>0, Nb>0, Ff<=us(Na+Nb)
%   2. Point A Tips                             Nb>0, Fs<=us(Na+Nb)
%   3. Point B Tips                             Na>0, Fs<=us(Na+Nb)
%   4. Slip                                     Na>0, Nb>0
%   5. Point A Tip and System Slip              Nb>0,
%   6. Point B Tip and System Slip              Na>0,
switch motion_type
    case "No motion"
        %Na>0 Nb>0 Ff<= us(Na+Nb)
        ax = 0; ay = 0; alpha = 0;
        Na_nm = subs(Na_general);
        Nb_nm = subs(Nb_general);
        Ff_nm = subs(Ff_general);
        if Na_nm > 0 && Nb_nm >0 && Ff_nm <= us*(Na+Nb)
            %The answer using for plotting later
            f1 = @(x,y) ax;
            f2 = @(x,y) ay;
            f3 = @(x,y) alpha;
            [Vx1,Vy1] = euler(f1, 0, 1, 10, 500); %Velocity in x direction
            [Vx2,Vy2] = euler(f2, 0, 1, 10, 500); %Velocity in y direction
            [Vx3,Vy3] = euler(f3, 0, 1, 10, 500); %Angular Velocity
%             figure(1)
%             plot(Vx1,Vy1,Vx2,Vy2,Vx3,Vy3)
%             legend('Velocity in x direction','Velocity in y direction', 'Angular Velocity');
%             title('Simulate Velocity with Time');
%             ylabel('Velocity [meter/second]');
%             xlabel('Time [seconds]');
            %Plotting Displacement
            f4 = @(x,y) ax.*x;
            f5 = @(x,y) ay.*x;
            f6 = @(x,y) alpha.*x;
            [Yx1,Yy1] = euler(f4, 0, 1, 10, 500);   %Position in x direction
            [Yx2,Yy2] = euler(f5, 0, 1, 10, 500);   %Position in y direction
            [Yx3,Yy3] = euler(f6, 0, 1, 10, 500);   %Angular Displacement
%             figure(2)
%             plot(Yx1,Yy1,Yx2,Yy2,Yx3,Yy3)
%             legend('Location in x direction','Location in y direction', 'Angular Displacement');
%             title('Simulate Location with Time');
%             ylabel('Location [meter]');
%             xlabel('Time [seconds]');
        else
            disp("Your checking situation won't happen given information provided!")
        end
    case "PointA Tip"
        %Applied General Motion Plane
        equ_ax = ax == (-H/2)*alpha;
        equ_ay = ay == (-B/2)*alpha;
        sol_AT = solve([equ_x, equ_y, equ_m, equ_ax, equ_ay], [Ff, Nb, ax, ay, alpha]);
        Ff_AT = sol_AT.Ff; Nb_AT = sol_AT.Nb; ax_AT = sol_AT.ax; ay_AT = sol_AT.ay; alpha_AT = sol_AT.alpha;
        %Because A tips, No normal force at point A
        Na = 0;
        Ff_AT = subs(Ff_AT); Nb_AT = subs(Nb_AT); ax_AT = subs(ax_AT); ay_AT = subs(ay_AT); alpha_AT = subs(alpha_AT);
        if Nb_AT > 0 && Ff_AT <= us*Nb_AT
            %Euler method Plotting Velocity
            f1 = @(x,y) ax_AT;
            f2 = @(x,y) ay_AT;
            f3 = @(x,y) alpha_AT;
            [Vx1,Vy1] = euler(f1, 0, 1, 10, 500); %Velocity in x direction
            [Vx2,Vy2] = euler(f2, 0, 1, 10, 500); %Velocity in y direction
            [Vx3,Vy3] = euler(f3, 0, 1, 10, 500); %Angular Velocity
            figure(1)
            plot(Vx1,Vy1,Vx2,Vy2,Vx3,Vy3)
            legend('Velocity in x direction','Velocity in y direction', 'Angular Velocity');
            title('Simulate Velocity with Time');
            ylabel('Velocity [meter/second]');
            xlabel('Time [seconds]');
            Normal_Force = zeros(length(Vy3), 1);
            index = 1;
            for ii = Vy3
                Normal_Force(index) = m*(ay_AT-(B/2)*alpha_AT-(H/2)*ii^2) + m*9.81 - P*sin(angle_pull);
                index = index + 1;
            end
            %Plotting Displacement
            f4 = @(x,y) ax_AT.*x;
            f5 = @(x,y) ay_AT.*x;
            f6 = @(x,y) alpha_AT.*x;
            [Yx1,Yy1] = euler(f4, 0, 1, 10, 500);   %Position in x direction
            [Yx2,Yy2] = euler(f5, 0, 1, 10, 500);   %Position in y direction
            [Yx3,Yy3] = euler(f6, 0, 1, 10, 500);   %Angular Displacement
            figure(2)
            plot(Yx1,Yy1,Yx2,Yy2,Yx3,Yy3)
            legend('Location in x direction','Location in y direction', 'Angular Displacement');
            title('Simulate Location with Time');
            ylabel('Location [meter]');
            xlabel('Time [seconds]');
            figure(3)
            plot(Yy3, Normal_Force);
            title('Normal Force at Point B with Time');
            ylabel('Normal Force [Newtons]');
            xlabel('Time [seconds]');
        else
            disp("Your checking situation won't happen given information provided!")
        end
    case "PointB Tip"
        %Applied General Motion Plane
        equ_ax = ax == (-H/2)*alpha;
        equ_ay = ay == (B/2)*alpha;
        sol_BT = solve([equ_x, equ_y, equ_m, equ_ax, equ_ay], [Ff, Na, ax, ay, alpha]);
        Ff_BT = sol_BT.Ff; Na_BT = sol_BT.Na; ax_BT = sol_BT.ax; ay_BT = sol_BT.ay; alpha_BT = sol_BT.alpha;
        %Because B tips, No normal force at point B
        Nb = 0;
        Ff_BT = subs(Ff_BT); Na_BT = subs(Na_BT); ax_BT = subs(ax_BT); ay_BT = subs(ay_BT); alpha_BT = subs(alpha_BT);
        if Na_BT > 0 && Ff_BT <= us*Na_BT
            %Euler method
            f1 = @(x,y) ax_BT;
            f2 = @(x,y) ay_BT;
            f3 = @(x,y) alpha_BT;
            [Vx1,Vy1] = euler(f1, 0, 1, 10, 500); %Velocity in x direction
            [Vx2,Vy2] = euler(f2, 0, 1, 10, 500); %Velocity in y direction
            [Vx3,Vy3] = euler(f3, 0, 1, 10, 500); %Angular Velocity
            figure(1)
            plot(Vx1,Vy1,Vx2,Vy2,Vx3,Vy3)
            legend('Velocity in x direction','Velocity in y direction', 'Angular Velocity');
            title('Simulate Velocity with Time');
            ylabel('Velocity [meter/second]');
            xlabel('Time [seconds]');
            %Plotting Displacement
            f4 = @(x,y) ax_BT.*x;
            f5 = @(x,y) ay_BT.*x;
            f6 = @(x,y) alpha_BT.*x;
            [Yx1,Yy1] = euler(f4, 0, 1, 10, 500);   %Position in x direction
            [Yx2,Yy2] = euler(f5, 0, 1, 10, 500);   %Position in y direction
            [Yx3,Yy3] = euler(f6, 0, 1, 10, 500);   %Angular Displacement
            figure(2)
            plot(Yx1,Yy1,Yx2,Yy2,Yx3,Yy3)
            legend('Location in x direction','Location in y direction', 'Angular Displacement');
            title('Simulate Location with Time');
            ylabel('Location [meter]');
            xlabel('Time [seconds]');
        else
            disp("Your checking situation won't happen given information provided!")
        end
    case "Slip"
        equ_fa = uk*(Na+Nb) == Ff; 
        equ_alpha = -P*sin(angle_pull)*(B/2) - P*cos(angle_pull)*(h-H/2) - Na*(B/2) + Nb*(B/2) - (H/2)*Ff == 0;
        equ_slip_ay = P*sin(angle_pull) + Na + Nb - m*g == 0;
        sol_Sp = solve([equ_x, equ_fa, equ_alpha, equ_slip_ay], [ax, Na, Nb, Ff]);
        ax_Sp = sol_Sp.ax; Na_Sp = sol_Sp.Na; Nb_Sp = sol_Sp.Nb;
        %Because Slip, no rotation, alpha = 0, ay=0
        if Na_Sp > 0 && Nb_Sp > 0 
            f1 = @(x,y) ax_Sp;
            [Vx1,Vy1] = euler(f1, 0, 1, 10, 500); %Velocity in x direction
            figure(1)
            plot(Vx1,Vy1)
            legend('Velocity in x direction');
            title('Simulate Velocity with Time');
            ylabel('Velocity [meter/second]');
            xlabel('Time [seconds]');
            %Plotting Displacement
            f4 = @(x,y) ax_Sp.*x;
            [Yx1,Yy1] = euler(f4, 0, 1, 10, 500);   %Position in x direction
            figure(2)
            plot(Yx1,Yy1,Yx2,Yy2,Yx3,Yy3)
            legend('Location in x direction','Location in y direction', 'Angular Displacement');
            title('Simulate Location with Time');
            ylabel('Location [meter]');
            xlabel('Time [seconds]');
        else
            disp("Your checking situation won't happen given information provided!")
        end
    
end

%% 画图并生成动画；
figure;
ah = gca;       % figure句柄；
% inStruct： 输入结构体参数；
    % inStruct.refW: 冰箱宽度;
    % inStruct.refH: 冰箱高度：
    % inStruct.refX: 冰箱横坐标（冰箱中线的下面那个点）
    % inStruct.refY: 冰箱纵坐标（冰箱中线的下面那个点）
    % inStruct.ang：冰箱角度（rad）（冰箱底线向右方向和水平地面右方向的夹角，冰箱前倾对应角度为负值，冰箱后仰对应角度为正值）；
inStruct.refW = B;
inStruct.refH = H;
tmp1 = ceil(sqrt(B^2+H^2));
frameData(length(Vx3)) = struct('cdata',[],'colormap',[]);
for iii=1:length(Vx3)
    drawnow;
    inStruct.refX = Yy1 - B/2;
    inStruct.refY = Yy2 - H/2;
    sp = [inStruct.refX, inStruct.refY];
    inStruct.ang = deg2rad(Vy3(iii));
    plot([sp(1)-tmp1, sp(1)+tmp1], [0, 0], 'k-', 'LineWidth', 2);
    hold on;
    axis([sp(1)-tmp1, sp(1)+tmp1, floor(0.125*tmp1), ceil(1.25*tmp1)]);
    axis equal;
    frameData(iii) = getP(inStruct, ah);
    hold off;
end

movie(frameData);
save([mfilename, '.mat']);

function [x,y]=euler(f,xinit,yinit,xfinal,n)
%Input: function, x initial value (time), y initial value
%       (velocity/displacement, x final value (time), number of increments
%Generate the numerical solution for any first order ODE with inital value
h=(xfinal-xinit)/n;
x=[xinit zeros(1,n)]; y=[yinit zeros(1,n)];
%Execute arithmetic euler's method
for i=1:n
    x(i+1)=x(i)+h;
    ynew=y(i)+h*f(x(i),y(i)); y(i+1)=y(i)+(h/2)*(f(x(i),y(i))+f(x(i+1),ynew));
end

end