function [T,X]=Main_Ref()
clear

myode = @myode;
[T,y]=ode45(myode,[0,0.5],[0,0,0,0,0,0]);

for i=1:length(y)
    [~, NA, NB, Ff] = myode(T(i), y(i,:));
    NAA(i) = NA;
    NBB(i) = NB;
    FFf(i) = Ff;
end

%Plot Normal Force
figure(1)
plot(T,[NAA; NBB; FFf], 'LineWidth',2);
legend('Normal_A','Normal_B', 'Friction');
title('Reaction Force');
ylabel('Force [newtons]');
xlabel('Time [seconds]');

figure(2)
plot(T, [y(:,4), y(:,5), y(:,6)], 'LineWidth',2)
legend('Accel_X','Accel_Y', 'Alpha');
title('Acceleration XY_alpha');
ylabel('Acceleration [meter/seconds^2]');
xlabel('Time [seconds]');

% figure(3)
% New_Array = diff(y(:,6), 2);
% plot(T(1:end-2), New_Array);
% title('Angular Displacement with time');
% ylabel('Angular Displacement [radian]');
% xlabel('Time [seconds]');
end
