% Get the state of the refrigerator at a certain time;

function frameData = getP(inStruct, ah)

% inStruct£º Enter the structure parameters£»
    % inStruct.refW: Width;
    % inStruct.refH: Height£º
    % inStruct.refX: Refrigerator x loc (the point below the center line of the refrigerator
    % inStruct.refY: Refrigerator y loc£¨the point below the center line of the refrigerator£©
    % inStruct.ang£ºRefrigerator angle£¨The Angle between the bottom line of the refrigerator to the right \
                                        %and the right direction of the horizontal ground, the Angle corresponding \
                                        %to the forward tilt of the refrigerator is negative, and the Angle corresponding \
                                        %to the backward tilt of the refrigerator is positive£©£»
% ah: figure handle£»

% frameData£ºFigure image data, as a frame of animation data£»


%% Draw a normal outline of the refrigerator;
refW = inStruct.refW;       % W£»
refH = inStruct.refH;       % H£»
ang = inStruct.ang;         % Angle;
LW = 2;                     % Line width£»
sp = [inStruct.refX, inStruct.refY];    % The coordinates of the lines in the refrigerator;

rotateMat = [cos(ang), -sin(ang); sin(ang), cos(ang)];  

coe = 1/40;     % The ratio of the distance between the door seam and the center line to the width of the refrigerator;
% Draw the outer line first;
refX1 = [sp(1), sp(1)+refW*coe, sp(1)+refW/2, sp(1)+refW/2, sp(1)+refW*coe, sp(1)-refW*coe, sp(1)-refW/2, sp(1)-refW/2, sp(1)-refW*coe, sp(1)];
refY1 = [sp(2), sp(2), sp(2), sp(2)+refH, sp(2)+refH, sp(2)+refH, sp(2)+refH, sp(2), sp(2), sp(2)];
% Rotate
tmp1 = rotateMat*[refX1; refY1];
refX1 = tmp1(1, :);
refY1 = tmp1(2, :);
% To prevent the refrigerator from falling below the horizon, so raise it a little bit;
% To prevent the refrigerator from hanging, so lower it a little bit;
tmp1 = min(refY1);
refY1 = refY1-tmp1;


% Draw a crack between a door and its frame£»
refX2 = [refX1(2), refX1(5)];
refY2 = [refY1(2), refY1(5)];
refX3 = [refX1(9), refX1(6)];
refY3 = [refY1(9), refY1(6)];

%generate the image
axes(ah)
plot(refX1, refY1, 'k-', 'LineWidth', LW);
hold on;
plot(refX2, refY2, 'b-', 'LineWidth', LW);
plot(refX3, refY3, 'b-', 'LineWidth', LW);
hold off;
frameData = getframe;

end

