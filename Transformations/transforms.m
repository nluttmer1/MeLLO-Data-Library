clear
clc
close all

%%

filename = 'Walker Static.csv';  % The name of the CSV file
data = readtable(filename);

IMUR_T = get_pt(data,5,12:14);
IMUR_M = get_pt(data,5,15:17);
IMUR_B = get_pt(data,5,18:20);

marker_x = IMUR_B;
marker_m = IMUR_M;
marker_b = IMUR_T;

close all
figure

[vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b);
plot_markers_and_transforms(vx,vy,vz,marker_x,marker_m,marker_b,'R')
T_R = [vx', vy', vz', marker_m'; 0 0 0 1];


IMUL_T = get_pt(data,5,3:5);
IMUL_M = get_pt(data,5,6:8);
IMUL_B = get_pt(data,5,9:11);

marker_x = IMUL_T;
marker_m = IMUL_M;
marker_b = IMUL_B;

[vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b);
plot_markers_and_transforms(vx,vy,vz,marker_x,marker_m,marker_b,'L')
T_L = [vx', vy', vz', marker_m'; 0 0 0 1];

T_RL = inv(T_R)*T_L

%%

filename = 'Wheelbarrow Static.csv';  % The name of the CSV file
data = readtable(filename);

IMUR_T = get_pt(data,5,18:20);
IMUR_M = get_pt(data,5,21:23);
IMUR_B = get_pt(data,5,24:26);

marker_x = IMUR_B;
marker_m = IMUR_M;
marker_b = IMUR_T;

close all
figure

[vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b);
plot_markers_and_transforms(vx,vy,vz,marker_x,marker_m,marker_b,'R')
T_R = [vx', vy', vz', marker_m'; 0 0 0 1];

IMUL_T = get_pt(data,5,27:29); % called R
IMUL_M = get_pt(data,5,30:32);
IMUL_B = get_pt(data,5,33:35); % called L

marker_x = IMUL_T;
marker_m = IMUL_M;
marker_b = IMUL_B;

[vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b);
plot_markers_and_transforms(vx,vy,vz,marker_x,marker_m,marker_b,'L')
T_L = [vx', vy', vz', marker_m'; 0 0 0 1];

T_RL = inv(T_R)*T_L


%%

filename = 'Shopping Cart Static.csv';  % The name of the CSV file
data = readtable(filename);

IMUR_T = get_pt(data,5,15:17);
IMUR_M = get_pt(data,5,18:20);
IMUR_B = get_pt(data,5,21:23);

marker_x = IMUR_B;
marker_m = IMUR_M;
marker_b = IMUR_T;

close all
figure

[vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b);
plot_markers_and_transforms(vx,vy,vz,marker_x,marker_m,marker_b,'R')
T_R = [vx', vy', vz', marker_m'; 0 0 0 1];

IMUL_T = get_pt(data,5,24:26);
IMUL_M = get_pt(data,5,27:29);
IMUL_B = get_pt(data,5,30:32);

marker_x = IMUL_T;
marker_m = IMUL_M;
marker_b = IMUL_B;

[vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b);
plot_markers_and_transforms(vx,vy,vz,marker_x,marker_m,marker_b,'L')
T_L = [vx', vy', vz', marker_m'; 0 0 0 1];

T_RL = inv(T_R)*T_L


%%
function pt = get_pt(data, row, cols)

x = str2double(cell2mat(data{row,cols(1)}));
y = str2double(cell2mat(data{row,cols(2)}));
z = str2double(cell2mat(data{row,cols(3)}));

pt = [x,y,z]/1000;

end

function [vx,vy,vz] = compute_vectors(marker_x, marker_m, marker_b)
    v1 = marker_x - marker_m;

    vx = v1/norm(v1);

    v2 = marker_b - marker_m;

    v3 = cross(v2, vx);
    vz = v3/norm(v3);
    vy = -cross(vx,vz);
end

function plot_markers_and_transforms(vx, vy, vz, marker_x, marker_m, marker_b, label)
    scale = 0.1;
    hold on
    
    % Plot the transformed vectors
    plot3([marker_m(1), scale*vx(1) + marker_m(1)], [marker_m(2), scale*vx(2) + marker_m(2)], [marker_m(3), scale*vx(3) + marker_m(3)], 'r')
    plot3([marker_m(1), scale*vy(1) + marker_m(1)], [marker_m(2), scale*vy(2) + marker_m(2)], [marker_m(3), scale*vy(3) + marker_m(3)], 'g')
    plot3([marker_m(1), scale*vz(1) + marker_m(1)], [marker_m(2), scale*vz(2) + marker_m(2)], [marker_m(3), scale*vz(3) + marker_m(3)], 'b')
    
    % Plot the markers
    plot3(marker_m(1), marker_m(2), marker_m(3), 'ko')
    plot3(marker_x(1), marker_x(2), marker_x(3), 'ro')
    plot3(marker_b(1), marker_b(2), marker_b(3), 'ko')
    
    % Add the label underneath marker_m
    % Offset the label position slightly along the Z-axis to place it below the marker
    text(marker_m(1), marker_m(2), marker_m(3) - 0.1, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10)
    
    % Plot settings
    axis equal
    grid on
end
