function [T,X]=Main_Ref()
clear

m = 10;             %kg 
P = 100;            %N
uk = 0.5;           %dimensionless
us = 0.8;           %dimensionless
B = 0.6;            %width m
H = 1.8;            %height of refrigerator m
d = 1.5;            %height of pulling force m
I = (1/12)*m*(B^2+H^2); %Mass Moment of inertia
g = 9.81;             %gravity acceleration

% Opt    = odeset('Events', @myEvent);
myode = @myode;
[T,y]=ode45(myode,[0,10],[0,0,0,0,0,0]);

for i=1:length(y)
    [~, NA, NB, Ff] = myode(T(i), y(i,:));
    NAA(i) = NA;
    NBB(i) = NB;
    FFf(i) = Ff;
end

%Plot Normal Force
figure(1)
plot(T,[NAA; NBB; FFf], 'LineWidth',2);
legend('Normal_A','Normal_B', 'Friction_B');
title('Reaction Force');
ylabel('Force [Newtons]');
xlabel('Time [Seconds]');

figure(2)
plot(T, [y(:,4), y(:,5), y(:,6)], 'LineWidth',2)
legend('a_X','a_Y', 'alpha');
title('Linear and angular Acceleration With Time');
ylabel('Acceleration [Meter/Seconds^2]');
xlabel('Time [Seconds]');

Yy1 = y(:,4);
Yy2 = y(:,5);
Vy3 = y(:,6);
figure;
ah = gca;       % figure handler；
inStruct.refW = B;
inStruct.refH = H;
tmp1 = ceil(sqrt(B^2+H^2));
frameData(length(T)) = struct('cdata',[],'colormap',[]);
for iii=1:length(T)
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

end