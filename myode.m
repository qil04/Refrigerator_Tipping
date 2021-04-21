function [dydt, NA, NB, Ff] = myode(t, y)
%%Part1 Assign values to the system (can be changed)
m = 10;             %kg 
P = 100;            %N
uk = 0.5;           %dimensionless
us = 0.8;           %dimensionless
B = 0.6;            %width m
H = 1.8;            %height of refrigerator m
d = 1.5;            %height of pulling force m
I = (1/12)*m*(B^2+H^2); %Mass Moment of inertia
g = 9.81;             %gravity acceleration

theta = y(5);
Sth = sin(theta);
Cth = cos(theta);
n = 6;
switch n
    case 3  %No slip; No tip
        %[NA, NB, Ff, ax, ay, alpha]
        A = [0, 0, 1,  0,  0, 0;...
             1, 1, 0,  0,  0, 0;...
             -B/2, B/2, -(H/2), 0, 0, 0];
        B = [P; m*g; (d-H/2)*P];
        X = A\B;
        dydt = [y(2); X(4); y(4); X(5); y(6); X(6)];
        NA = X(1); NB = X(2); Ff = X(3);
    case 4  %No tip; Slip
        %[NA, NB, Ff, ax, ay, alpha]
        A = [0, 0, 1,  m,  0, 0;...
             1, 1, 0,  0,  0, 0;...
             -B/2, B/2, -(H/2), 0, 0, 0;...
             uk, uk, -1, 0, 0, 0];
        B = [P; m*g; (d-H/2)*P; 0];
        X = A\B;
        dydt = [y(2); X(4); y(4); X(5); y(6); X(6)];
        NA = X(1); NB = X(2); Ff = X(3);
    case 5  %No slip; B Tip
        %[NA, NB, Ff, ax, ay, alpha]
        A = [0, 0, 1/2,  m, 0, 0;...
             0, 1, 0, 0, -m, 0;...
             0, B/2, -(H/4), 0, 0, I;...
             0, 0, 0, -1, 0, (B/2)*Sth-(H/2)*Cth;...
             0, 0, 0, 0, -1, -(B/2)*Cth-(H/2)*Sth];
        B = [P; m*g; (d-H/2)*P; ... 
             y(6)^2*(-B*Cth/2-H*Sth/2); y(6)^2*(-B*Sth/2+H*Cth/2)];
        X = A\B;
        dydt = [y(2); X(4); y(4); X(5); y(6); X(6)];
        NA = X(1); NB = X(2); Ff = X(3);
    case 6  %Slip; Tip
        %[NA, NB, Ff, ax, ay, alpha]
        A = [0, 0, 1/2,  m, 0, 0;...
             0, 1, 0, 0, -m, 0;...
             0, B/2, -(H/4), 0, 0, I;...
             0, 0, 0, -1, 0, (B/2)*Sth-(H/2)*Cth;...
             0, 0, 0, 0, -1, -(B/2)*Cth-(H/2)*Sth;...
             0, -uk, 1, 0, 0, 0];
        B = [P; m*g; (d-H/2)*P; ... 
             y(6)^2*(-B*Cth/2-H*Sth/2); y(6)^2*(-B*Sth/2+H*Cth/2); 0];
        X = A\B;
        NA = 0; NB = X(2); Ff = X(3);
        dydt = [y(2); X(4); y(4); X(5); y(6); X(6)];
end
