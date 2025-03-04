%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nathaniel Luttmer 
% Combine one subject/object data
% Library Paper
% Version: January 3, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: NATE L Wheelbarrow LIFT in paper
% This shows the sample code plot and the 
% 10 signals put together
% Also exports the final dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NATE NATE NATE - CHECK WHAT KNEE ANGLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YOU ARE USING RIGHT NOW I AM SWITCHING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINE 341 !!!!!!!!!!!!!!!!!!!!!!!!!!

clear all 
close all 
clc
% This should be fully automated for syncing data, only requiring files as
% an input from the user

[Devices_VICON,  Trajectories_VICON, Trajectories_Cut, AddBio_Data, DSpace_Data, IMU_Data] = importExperimentData(5,0,31,1);

%%%%%%%%%% VICON DATA SHIFT
Devices_time = Devices_VICON(1:end,1)/100;
Devices_Impulse = Devices_VICON(1:end,3);

% Shift to the impulses in VICON
for i = 1:length(Devices_time)
    if Devices_Impulse(i) > 1
        Vicon_shift = Devices_time(i); % VICON SHIFT VARIABLE
        break
    end
end

Devices_time = Devices_time - Vicon_shift;

Traj_Vicon_Time = Trajectories_VICON(1:end,1)/100;
Traj_Vicon_Time = Traj_Vicon_Time - Vicon_shift;
Traj_Vicon_RTOE = Trajectories_VICON(1:end,200);

Trajectories_Cut_time = Trajectories_Cut(1:end,1)/100 - Vicon_shift;
Trajectories_Cut_data = Trajectories_Cut(1:end,200);

% Shift to where the VICON data was cut
Cut_shift_Start = Trajectories_Cut_time(1,1);
Cut_shift_End = Trajectories_Cut_time(end,1);
Traj_Cut_Shift_Time_Start = Trajectories_Cut_time - Cut_shift_Start;

%%%%%%%%% ADDBIOMECHANICS DATA
AddBio_Time = AddBio_Data.data(:,1);
Delta_T = AddBio_Time(3,1)-AddBio_Time(2,1);

Right_Knee_pos = AddBio_Data.data(:,11)*180/pi;
Left_Knee_pos = AddBio_Data.data(:,18)*180/pi;

Right_Knee_vel = zeros(length(Right_Knee_pos)-1,1);
Left_Knee_vel = zeros(length(Left_Knee_pos)-1,1);

for i = 1:length(Right_Knee_pos)-1
    Right_Knee_vel(i,1) = (Right_Knee_pos(i+1) - Right_Knee_pos(i))/Delta_T;
    Left_Knee_vel(i,1) = (Left_Knee_pos(i+1) - Left_Knee_pos(i))/Delta_T;
end

Right_Knee_acc = zeros(length(Right_Knee_vel)-1,1);
Left_Knee_acc = zeros(length(Left_Knee_vel)-1,1);

for i = 1:length(Right_Knee_vel)-1
    Right_Knee_acc(i,1) = (Right_Knee_vel(i+1) - Right_Knee_vel(i))/Delta_T;
    Left_Knee_acc(i,1) = (Left_Knee_vel(i+1) - Left_Knee_vel(i))/Delta_T;
end

%%%%%%% DSPACE DATA SHIFT
fields = fieldnames(DSpace_Data);
DSpace_time = DSpace_Data.(fields{1}).X(1).Data;
DSpace_Impulse = DSpace_Data.(fields{1}).Y(20).Data;

for i = 1:length(DSpace_time)
    if DSpace_Impulse(i) > 1
        dSpace_shift = DSpace_time(i);
        break
    end
end
DSpace_time = DSpace_time - dSpace_shift;

Fy1 = DSpace_Data.(fields{1}).Y(8).Data;

%%%%%%%%% IMU Shift
IMU_Acc_Data = IMU_Data(:,1:3);
IMU_Gyro_Data = IMU_Data(:,4:6);

dSpace_Acc_time1 = DSpace_Data.(fields{1}).X(1).Data;

IMU_Acc_y = nonzeros(IMU_Acc_Data(:,2));
IMU_Acc_time = ((1:1:length(IMU_Acc_y))')/1000;
IMU_Gyro_x = nonzeros(IMU_Gyro_Data(:,1));
IMU_Gyro_time = ((1:1:length(IMU_Gyro_x))')/1000;

Fy1_new = Fy1 - Fy1(1);

IMU_Acc_y_new = IMU_Acc_y - IMU_Acc_y(1);
IMU_Acc_y_new = IMU_Acc_y_new*10;

count_F = 0;

for i = 1:1:length(Fy1_new)
    if Fy1_new(i) > 3
        Fx1_time = dSpace_Acc_time1(i);
        count_F = count_F + 1;
        break
    end
end

count_G = 0;

for i = 1:1:length(IMU_Acc_y_new)
    if IMU_Acc_y_new(i) < -2
        Acc_time = IMU_Acc_time(i);
        count_G = count_G + 1;
        break
    end
end

dSpace_Acc_time1 = dSpace_Acc_time1 - dSpace_shift;% +0.4020;
IMU_Shift = diff([Fx1_time, dSpace_shift]); %% Magic happening

IMU_Gyro_time = IMU_Gyro_time - Acc_time - IMU_Shift;
IMU_Acc_time = IMU_Acc_time - Acc_time - IMU_Shift;


%%%% PLOTTING
figure(1)
subplot(4,1,1)
hold on
plot(Devices_time,Devices_Impulse,'k-')
plot(DSpace_time,DSpace_Impulse,'r-')
title('VICON Dspace Impulse')
hold off

subplot(4,1,2)
hold on
plot(Traj_Vicon_Time,Traj_Vicon_RTOE,'k-')
plot(Trajectories_Cut_time,Trajectories_Cut_data,'r-')
title('VICON/VICON Cut')
hold off

subplot(4,1,3)
plot(Traj_Cut_Shift_Time_Start,Trajectories_Cut_data)
title('VICON Cut Shift')

subplot(4,1,4)
hold on 
plot(IMU_Acc_time,IMU_Acc_y_new,'r-')
plot(dSpace_Acc_time1,Fy1_new,'k-')
legend('Accel','Load Cell 1')
xlabel('Time')
ylabel('Signal')
title('dSpace/IMU Sync')
grid on 
grid minor
hold off


% figure(2)
% subplot(3,1,1)
% plot(Time,Right_Knee_pos,'k-',Time,Left_Knee_pos,'r-')
% xlabel('Time [s]')
% ylabel('Postion [deg]')
% title('RKnee Position')
% legend('Right Knee','Left Knee')
% 
% subplot(3,1,2)
% plot(Time(1:end-1),Right_Knee_vel,'k-',Time(1:end-1),Left_Knee_vel,'r-')
% xlabel('Time [s]')
% ylabel('Velocity [deg/s]')
% title('RKnee Velocity')
% legend('Right Knee','Left Knee')
% 
% subplot(3,1,3)
% plot(Time(1:end-2),Right_Knee_acc,'k-',Time(1:end-2),Left_Knee_acc,'r-')
% xlabel('Time [s]')
% ylabel('Acceleration [deg/s2]')
% title('RKnee Acceleration')
% legend('Right Knee','Left Knee')

figure
hold on
plot(IMU_Acc_time-Cut_shift_Start,IMU_Acc_y_new,'r-')
plot(dSpace_Acc_time1-Cut_shift_Start,Fy1_new,'k-')
plot(Traj_Vicon_Time-Cut_shift_Start,Traj_Vicon_RTOE,'b-')
plot(Trajectories_Cut_time-Cut_shift_Start,Trajectories_Cut_data,'g-')
legend('IMU','Load Cell','VICON','VICON Cut')
hold off

%%%%%% Data to extract/check plots

New_IMU_Time = round((IMU_Acc_time-Cut_shift_Start)*1000)/1000;
New_LC_Time = round((dSpace_Acc_time1-Cut_shift_Start)*1000)'/1000;
New_Cut_Time = Trajectories_Cut_time-Cut_shift_Start;
New_AddBio_Time = AddBio_Time;
AddBio_Start = New_AddBio_Time(1,1);
AddBio_End = New_AddBio_Time(end,1);

[IMU_Row_Start,IMU_Col_Start] = find(New_IMU_Time == AddBio_Start,1);
[IMU_Row_End,IMU_Col_End] = find(New_IMU_Time == AddBio_End,1);
[LC_Row_Start,LC_Col_Start] = find(New_LC_Time == AddBio_Start,1);
[LC_Row_End,LC_Col_End] = find(New_LC_Time == AddBio_End,1);

IMU_Cut_Time = New_IMU_Time(IMU_Row_Start:IMU_Row_End)';
IMU_Acc_Data_x = nonzeros(IMU_Acc_Data(:,1));
IMU_Acc_Data_y = nonzeros(IMU_Acc_Data(:,2));
IMU_Acc_Data_z = nonzeros(IMU_Acc_Data(:,3));
IMU_AccCut_Data_x = IMU_Acc_Data_x(IMU_Row_Start:IMU_Row_End);
IMU_AccCut_Data_y = IMU_Acc_Data_y(IMU_Row_Start:IMU_Row_End);
IMU_AccCut_Data_z = IMU_Acc_Data_z(IMU_Row_Start:IMU_Row_End);


IMU_Gyro_Data_x = nonzeros(IMU_Gyro_Data(:,1));
IMU_Gyro_Data_y = nonzeros(IMU_Gyro_Data(:,2));
IMU_Gyro_Data_z = nonzeros(IMU_Gyro_Data(:,3));
IMU_GyroCut_Data_x = IMU_Gyro_Data_x(IMU_Row_Start:IMU_Row_End);
IMU_GyroCut_Data_y = IMU_Gyro_Data_y(IMU_Row_Start:IMU_Row_End);
IMU_GyroCut_Data_z = IMU_Gyro_Data_z(IMU_Row_Start:IMU_Row_End);

LC_Cut_Time = New_LC_Time(LC_Row_Start:LC_Row_End);
Fx1 = DSpace_Data.(fields{1}).Y(5).Data;
Fy1 = DSpace_Data.(fields{1}).Y(8).Data;
Fz1 = DSpace_Data.(fields{1}).Y(11).Data;
Tx1 = DSpace_Data.(fields{1}).Y(22).Data;
Ty1 = DSpace_Data.(fields{1}).Y(25).Data;
Tz1 = DSpace_Data.(fields{1}).Y(28).Data;

Fx_Cut_Data = Fx1(LC_Row_Start:LC_Row_End); % CHANGED THIS TO EASILY GET THE RESULTANT - CHANGE BACK
Fy1_Cut_Data = Fy1(LC_Row_Start:LC_Row_End);
Fz1_Cut_Data = Fz1(LC_Row_Start:LC_Row_End);

Fx1_Cut_Data = sqrt(Fx_Cut_Data.^2 + Fy1_Cut_Data.^2 + Fz1_Cut_Data.^2); % CHANGED THIS TO EASILY GET THE RESULTANT - CHANGE BACK

Tx_Cut_Data = Tx1(LC_Row_Start:LC_Row_End); % CHANGED THIS TO EASILY GET THE RESULTANT - CHANGE BACK
Ty1_Cut_Data = Ty1(LC_Row_Start:LC_Row_End);
Tz1_Cut_Data = Tz1(LC_Row_Start:LC_Row_End);

Tx1_Cut_Data = sqrt(Tx_Cut_Data.^2 + Ty1_Cut_Data.^2 + Tz1_Cut_Data.^2); % CHANGED THIS TO EASILY GET THE RESULTANT - CHANGE BACK


figure % SUBPLOT ATTEMPT
ax1 = subplot(4,1,1);
plot(AddBio_Time,Right_Knee_pos,'k-',AddBio_Time,Left_Knee_pos,'r-')
xlabel('Time [s]')
ylabel('Postion [deg]')
title('RKnee Position')
legend('Right Knee','Left Knee')%,'Location','eastoutside')

ax2 = subplot(4,1,2);
plot(IMU_Cut_Time,IMU_AccCut_Data_x,'k-',IMU_Cut_Time,IMU_AccCut_Data_y,'r-', ...
    IMU_Cut_Time,IMU_AccCut_Data_z,'b-')
xlabel('Time [s]')
ylabel('Acc [g]')
title('IMU Acceleration')
legend('AccX','AccY','AccZ')%,'Location','eastoutside')

ax3 = subplot(4,1,3);
hold on
plot(LC_Cut_Time,Fx1_Cut_Data,'k-')
plot(LC_Cut_Time,Fy1_Cut_Data,'r-')
plot(LC_Cut_Time,Fz1_Cut_Data,'b-')
xlabel('Time [s]')
ylabel('Force [N]')
title('LC Force')
legend('Fx','Fy','Fz')%,'Location','eastoutside')
hold off

ax4 = subplot(4,1,4);
hold on
plot(LC_Cut_Time,Tx1_Cut_Data,'k-')
plot(LC_Cut_Time,Ty1_Cut_Data,'r-')
plot(LC_Cut_Time,Tz1_Cut_Data,'b-')
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('LC Torque')
legend('Tx','Ty','Tz')%,'Location','eastoutside')
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');


figure % NEXTTILE ATTEMPT
% Create a tiled layout for uniform subplots
t = tiledlayout(4, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Plot in each tile
ax21 = nexttile; 
plot(AddBio_Time,Right_Knee_pos,'k-',AddBio_Time,Left_Knee_pos,'r-')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Postion [deg]')
title('Knee Position')
legend('Right Knee','Left Knee','Location','eastoutside')

ax22 = nexttile; 
plot(IMU_Cut_Time,IMU_AccCut_Data_x,'k-',IMU_Cut_Time,IMU_AccCut_Data_y,'r-', ...
    IMU_Cut_Time,IMU_AccCut_Data_z,'b-')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Acc [g]')
title('IMU Acceleration')
legend('IMU AccX','IMU AccY','IMU AccZ','Location','eastoutside')

ax23 = nexttile; 
hold on
plot(LC_Cut_Time,Fx1_Cut_Data,'k-',LC_Cut_Time,Fy1_Cut_Data,'r-',LC_Cut_Time,Fz1_Cut_Data,'b-')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Force [N]')
title('LC Force')
legend('Load Cell Fx','Load Cell Fy','Load Cell Fz','Location','eastoutside')
hold off

ax24 = nexttile; 
hold on
plot(LC_Cut_Time,Tx1_Cut_Data,'k-',LC_Cut_Time,Ty1_Cut_Data,'r-',LC_Cut_Time,Tz1_Cut_Data,'b-')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Torque [Nm]')
title('LC Torque')
legend('Load Cell Tx','Load Cell Ty','Load Cell Tz','Location','eastoutside')
hold off

linkaxes([ax21,ax22,ax23,ax24],'x');

%% Cutting the full data set into each lift and combining the plots

Cut_time = Trajectories_Cut_time-Cut_shift_Start; % Starts at 0
Trajectories_Cut_data = Trajectories_Cut(1:end,200);

LKnee_Marker_x = Trajectories_Cut(1:end,156);
LKnee_Marker_y = Trajectories_Cut(1:end,157);
LKnee_Marker_z = Trajectories_Cut(1:end,158);

StartLook = find(AddBio_Time == AddBio_Time(1,1));
EndLook = find(AddBio_Time == AddBio_Time(end,1));

Left_Knee_pos = Right_Knee_pos; % SWITCHES TO RIGHT KNEE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% getting a x axis for slope
x = (1:length(Left_Knee_pos)).';

% Making a variable to change whenever we get a negative slope
hit = 0;

% Start time vector ie whenever slope is negative for 25 points
Start_vect = [];
Start_Index = 1;

for j = StartLook:1:EndLook
    
    % Breaks once we find all the negative slopes
    if length(Start_vect) >= 20

        break;

    end

    if j + 24 ~= EndLook
        y1 = Left_Knee_pos(j,1);
        y2 = Left_Knee_pos(j+1,1);
        y3 = Left_Knee_pos(j+2,1);
        y4 = Left_Knee_pos(j+3,1);
        y5 = Left_Knee_pos(j+4,1);
        y6 = Left_Knee_pos(j+5,1);
        y7 = Left_Knee_pos(j+6,1);
        y8 = Left_Knee_pos(j+7,1);
        y9 = Left_Knee_pos(j+8,1);
        y10 = Left_Knee_pos(j+9,1);
        y11 = Left_Knee_pos(j+10,1);
        y12 = Left_Knee_pos(j+11,1);
        y13 = Left_Knee_pos(j+12,1);
        y14 = Left_Knee_pos(j+13,1);
        y15 = Left_Knee_pos(j+14,1);
        y16 = Left_Knee_pos(j+15,1);
        y17 = Left_Knee_pos(j+16,1);
        y18 = Left_Knee_pos(j+17,1);
        y19 = Left_Knee_pos(j+18,1);
        y20 = Left_Knee_pos(j+19,1);
        y21 = Left_Knee_pos(j+20,1);
        y22 = Left_Knee_pos(j+21,1);
        y23 = Left_Knee_pos(j+22,1);
        y24 = Left_Knee_pos(j+23,1);
        y25 = Left_Knee_pos(j+24,1);

        x1 = x(j,1);
        x2 = x(j+1,1);
        x3 = x(j+2,1);
        x4 = x(j+3,1);
        x5 = x(j+4,1);
        x6 = x(j+5,1);
        x7 = x(j+6,1);
        x8 = x(j+7,1);
        x9 = x(j+8,1);
        x10 = x(j+9,1);
        x11 = x(j+10,1);
        x12 = x(j+11,1);
        x13 = x(j+12,1);
        x14 = x(j+13,1);
        x15 = x(j+14,1);
        x16 = x(j+15,1);
        x17 = x(j+16,1);
        x18 = x(j+17,1);
        x19 = x(j+18,1);
        x20 = x(j+19,1);
        x21 = x(j+20,1);
        x22 = x(j+21,1);
        x23 = x(j+22,1);
        x24 = x(j+23,1);
        x25 = x(j+24,1);

        slope1 = (y25 - y24) / (x25 - x24);
        slope2 = (y24 - y23) / (x24 - x23);
        slope3 = (y23 - y22) / (x23 - x22);
        slope4 = (y22 - y21) / (x22 - x21);
        slope5 = (y21 - y20) / (x21 - x20);
        slope6 = (y20 - y19) / (x20 - x19);
        slope7 = (y19 - y18) / (x19 - x18);
        slope8 = (y18 - y17) / (x18 - x17);
        slope9 = (y17 - y16) / (x17 - x16);
        slope10 = (y16 - y15) / (x16 - x15);
        slope11 = (y15 - y14) / (x15 - x14);
        slope12 = (y14 - y13) / (x14 - x13);
        slope13 = (y13 - y12) / (x13 - x12);
        slope14 = (y12 - y11) / (x12 - x11);
        slope15 = (y11 - y10) / (x11 - x10);
        slope16 = (y10 - y9) / (x10 - x9);
        slope17 = (y9 - y8) / (x9 - x8);
        slope18 = (y8 - y7) / (x8 - x7);
        slope19 = (y7 - y6) / (x7 - x6);
        slope20 = (y6 - y5) / (x6 - x5);
        slope21 = (y5 - y4) / (x5 - x4);
        slope22 = (y4 - y3) / (x4 - x3);
        slope23 = (y3 - y2) / (x3 - x2);
        slope24 = (y2 - y1) / (x2 - x1);
    
        
        if (hit == 0) && (slope1 < 0) && (slope2 < 0) && (slope3 < 0) && (slope4 < 0)...
                && (slope5 < 0) && (slope6 < 0) && (slope7 < 0) && (slope8 < 0)...
                && (slope9 < 0) && (slope10 < 0) && (slope11 < 0) && (slope12 < 0) && (slope13 < 0)...
                && (slope14 < 0) && (slope15 < 0) && (slope16 < 0) && (slope17 < 0)...
                && (slope18 < 0) && (slope19 < 0) && (slope20 < 0) && (slope21 < 0) && (slope22 < 0)...
                && (slope23 < 0) && (slope24 < 0)
    
            Start_vect(Start_Index) = x1;
            Start_Index = Start_Index + 1;
            hit = 1;
        end
        
    % found a positive slope and will be looking for the next negative slope 
        if (hit == 1) && (slope1 > 0) && (slope2 > 0) && (slope3 > 0) && (slope4 > 0)...
                && (slope5 > 0) && (slope6 > 0) && (slope7 > 0) && (slope8 > 0)...
                && (slope9 > 0) && (slope10 > 0) && (slope11 > 0) && (slope12 > 0) && (slope13 > 0)...
                && (slope14 > 0) && (slope15 > 0) && (slope16 > 0) && (slope17 > 0)...
                && (slope18 > 0) && (slope19 > 0) && (slope20 > 0) && (slope21 > 0) && (slope22 > 0)...
                && (slope23 > 0) && (slope24 > 0)
            
            hit = 0;
    
        end

    else 
        break;
    end

end

clear x x1 x10 x2 x3 x4 x5 x6 x7 x8 x9 y1 y10 y2 y3 y4 y5 y6 y7 y8 y9...
    slope9 slope1 slope2 slope3 slope4 slope5 slope6 slope7 slope8...
    j hit slope10 slope11 slope12 slope13 slope14 slope15 slope16...
    slope17 slope18 slope19 slope20 slope21 slope22 slope23 slope24...
    Start_Index x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24...
    x25 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25


%%%%% Repeating processs to cut out the holding portion on data
% getting a x axis for slope
x = (1:length(Left_Knee_pos)).';

% Making a variable to change whenever we get a negative slope
hit = 0;

% Start time vector ie whenever slope is negative for 25 points
End_vect = [];
End_Index = 1;
reverse = flip(StartLook:1:EndLook);  

for j = reverse
    
    % Breaks once we find all the negative slopes
    if length(End_vect) >= 20

        break;

    end


    if j + 24 ~= StartLook
        y1 = Left_Knee_pos(j,1);
        y2 = Left_Knee_pos(j-1,1);
        y3 = Left_Knee_pos(j-2,1);
        y4 = Left_Knee_pos(j-3,1);
        y5 = Left_Knee_pos(j-4,1);
        y6 = Left_Knee_pos(j-5,1);
        y7 = Left_Knee_pos(j-6,1);
        y8 = Left_Knee_pos(j-7,1);
        y9 = Left_Knee_pos(j-8,1);
        y10 = Left_Knee_pos(j-9,1);
        y11 = Left_Knee_pos(j-10,1);
        y12 = Left_Knee_pos(j-11,1);
        y13 = Left_Knee_pos(j-12,1);
        y14 = Left_Knee_pos(j-13,1);
        y15 = Left_Knee_pos(j-14,1);
        y16 = Left_Knee_pos(j-15,1);
        y17 = Left_Knee_pos(j-16,1);
        y18 = Left_Knee_pos(j-17,1);
        y19 = Left_Knee_pos(j-18,1);
        y20 = Left_Knee_pos(j-19,1);
        y21 = Left_Knee_pos(j-20,1);
        y22 = Left_Knee_pos(j-21,1);
        y23 = Left_Knee_pos(j-22,1);
        y24 = Left_Knee_pos(j-23,1);
        y25 = Left_Knee_pos(j-24,1);

        x1 = x(j,1);
        x2 = x(j-1,1);
        x3 = x(j-2,1);
        x4 = x(j-3,1);
        x5 = x(j-4,1);
        x6 = x(j-5,1);
        x7 = x(j-6,1);
        x8 = x(j-7,1);
        x9 = x(j-8,1);
        x10 = x(j-9,1);
        x11 = x(j-10,1);
        x12 = x(j-11,1);
        x13 = x(j-12,1);
        x14 = x(j-13,1);
        x15 = x(j-14,1);
        x16 = x(j-15,1);
        x17 = x(j-16,1);
        x18 = x(j-17,1);
        x19 = x(j-18,1);
        x20 = x(j-19,1);
        x21 = x(j-20,1);
        x22 = x(j-21,1);
        x23 = x(j-22,1);
        x24 = x(j-23,1);
        x25 = x(j-24,1);

        slope1 = (y25 - y24) / (x25 - x24);
        slope2 = (y24 - y23) / (x24 - x23);
        slope3 = (y23 - y22) / (x23 - x22);
        slope4 = (y22 - y21) / (x22 - x21);
        slope5 = (y21 - y20) / (x21 - x20);
        slope6 = (y20 - y19) / (x20 - x19);
        slope7 = (y19 - y18) / (x19 - x18);
        slope8 = (y18 - y17) / (x18 - x17);
        slope9 = (y17 - y16) / (x17 - x16);
        slope10 = (y16 - y15) / (x16 - x15);
        slope11 = (y15 - y14) / (x15 - x14);
        slope12 = (y14 - y13) / (x14 - x13);
        slope13 = (y13 - y12) / (x13 - x12);
        slope14 = (y12 - y11) / (x12 - x11);
        slope15 = (y11 - y10) / (x11 - x10);
        slope16 = (y10 - y9) / (x10 - x9);
        slope17 = (y9 - y8) / (x9 - x8);
        slope18 = (y8 - y7) / (x8 - x7);
        slope19 = (y7 - y6) / (x7 - x6);
        slope20 = (y6 - y5) / (x6 - x5);
        slope21 = (y5 - y4) / (x5 - x4);
        slope22 = (y4 - y3) / (x4 - x3);
        slope23 = (y3 - y2) / (x3 - x2);
        slope24 = (y2 - y1) / (x2 - x1);

        if (hit == 0) && (slope1 > 0) && (slope2 > 0) && (slope3 > 0) && (slope4 > 0)...
            && (slope5 > 0) && (slope6 > 0) && (slope7 > 0) && (slope8 > 0)...
            && (slope9 > 0) && (slope10 > 0) && (slope11 > 0) && (slope12 > 0) && (slope13 > 0)...
            && (slope14 > 0) && (slope15 > 0) && (slope16 > 0) && (slope17 > 0)...
            && (slope18 > 0) && (slope19 > 0) && (slope20 > 0) && (slope21 > 0) && (slope22 > 0)...
            && (slope23 > 0) && (slope24 > 0)

            End_vect(End_Index) = x25;
            End_Index = End_Index + 1;
            hit = 1;
        end
    
    % found a positive slope and will be looking for the next negative slope 
        if (hit == 1) && (slope1< 0) && (slope2< 0) && (slope3< 0) && (slope4< 0)...
                && (slope5< 0) && (slope6< 0) && (slope7< 0) && (slope8< 0)...
                && (slope9< 0) && (slope10< 0) && (slope11< 0) && (slope12< 0) && (slope13< 0)...
                && (slope14< 0) && (slope15< 0) && (slope16< 0) && (slope17< 0)...
                && (slope18< 0) && (slope19< 0) && (slope20< 0) && (slope21< 0) && (slope22< 0)...
                && (slope23< 0) && (slope24 < 0)
            
            hit = 0;
    
        end


    else 
        break;
    end
end

Start_vect_time = [];
End_vect_time = [];
for k = 1:length(Start_vect)
    Start_vect_time(k) = AddBio_Time(Start_vect(k));
    End_vect_time(k) = AddBio_Time(End_vect(k));
end

% Fliping end vect so its the same as start
End_vect = flip(End_vect);
End_vect_time = flip(End_vect_time);

Start_cut = [];
Start_cut_time = [];
End_cut = [];
End_cut_time = [];
for i = 2:2:length(End_vect)
    Start_cut(i) = Start_vect(i-1);
    Start_cut_time(i) = End_vect_time(i-1);
    End_cut(i) = Start_vect(i);
    End_cut_time(i) = End_vect_time(i);
end
    
Start_cut = nonzeros(Start_cut)';
Start_cut_time = nonzeros(Start_cut_time)';
End_cut = nonzeros(End_cut)';
End_cut_time = nonzeros(End_cut_time)';
%%
% Verify Plot
figure % NEXTTILE ATTEMPT
% Create a tiled layout for uniform subplots
t = tiledlayout(4, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Plot in each tile
ax21 = nexttile; 
hold on
plot(AddBio_Time,Right_Knee_pos,'k-',AddBio_Time,Left_Knee_pos,'r-')
xline(Start_cut_time,'--r')
xline(End_cut_time,'--b')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Postion [deg]')
title('Knee Position')
legend('Right Knee','Left Knee','Location','eastoutside')
hold off

ax22 = nexttile; 
hold on
plot(IMU_Cut_Time,IMU_AccCut_Data_x,'k-',IMU_Cut_Time,IMU_AccCut_Data_y,'r-', ...
    IMU_Cut_Time,IMU_AccCut_Data_z,'b-')
xline(Start_cut_time,'--r')
xline(End_cut_time,'--b')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Acc [g]')
title('IMU Acceleration')
legend('IMU AccX','IMU AccY','IMU AccZ','Location','eastoutside')
hold off

ax23 = nexttile; 
hold on
plot(LC_Cut_Time,Fx1_Cut_Data,'k-',LC_Cut_Time,Fy1_Cut_Data,'r-',LC_Cut_Time,Fz1_Cut_Data,'b-')
xline(Start_cut_time,'--r')
xline(End_cut_time,'--b')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Force [N]')
title('LC Force')
legend('Load Cell Fx','Load Cell Fy','Load Cell Fz','Location','eastoutside')
hold off

ax24 = nexttile; 
hold on
plot(LC_Cut_Time,Tx1_Cut_Data,'k-',LC_Cut_Time,Ty1_Cut_Data,'r-',LC_Cut_Time,Tz1_Cut_Data,'b-')
xline(Start_cut_time,'--r')
xline(End_cut_time,'--b')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Torque [Nm]')
title('LC Torque')
legend('Load Cell Tx','Load Cell Ty','Load Cell Tz','Location','eastoutside')
hold off

linkaxes([ax21,ax22,ax23,ax24],'x');
%%
% Creating Cuts
Time1 = AddBio_Time(Start_cut(1):End_cut(1));
Wavelength1 = Left_Knee_pos(Start_cut(1):End_cut(1));

Right_Hip = AddBio_Data.data(:,8)*180/pi;
Right_Ankle = AddBio_Data.data(:,12)*180/pi;

% Cut AddBio
AddBio_Cut_Time1 = AddBio_Time(Start_cut(1):End_cut(1));
AddBio_Cut_Time2 = AddBio_Time(Start_cut(2):End_cut(2));
AddBio_Cut_Time3 = AddBio_Time(Start_cut(3):End_cut(3));
AddBio_Cut_Time4 = AddBio_Time(Start_cut(4):End_cut(4));
AddBio_Cut_Time5 = AddBio_Time(Start_cut(5):End_cut(5));
AddBio_Cut_Time6 = AddBio_Time(Start_cut(6):End_cut(6));
AddBio_Cut_Time7 = AddBio_Time(Start_cut(7):End_cut(7));
AddBio_Cut_Time8 = AddBio_Time(Start_cut(8):End_cut(8));
AddBio_Cut_Time9 = AddBio_Time(Start_cut(9):End_cut(9));
AddBio_Cut_Time10 = AddBio_Time(Start_cut(10):End_cut(10));

AddBio_LKnee1 = Left_Knee_pos(Start_cut(1):End_cut(1));
AddBio_LKnee2 = Left_Knee_pos(Start_cut(2):End_cut(2));
AddBio_LKnee3 = Left_Knee_pos(Start_cut(3):End_cut(3));
AddBio_LKnee4 = Left_Knee_pos(Start_cut(4):End_cut(4));
AddBio_LKnee5 = Left_Knee_pos(Start_cut(5):End_cut(5));
AddBio_LKnee6 = Left_Knee_pos(Start_cut(6):End_cut(6));
AddBio_LKnee7 = Left_Knee_pos(Start_cut(7):End_cut(7));
AddBio_LKnee8 = Left_Knee_pos(Start_cut(8):End_cut(8));
AddBio_LKnee9 = Left_Knee_pos(Start_cut(9):End_cut(9));
AddBio_LKnee10 = Left_Knee_pos(Start_cut(10):End_cut(10));

AddBio_RHip1 = Right_Hip(Start_cut(1):End_cut(1));
AddBio_RHip2 = Right_Hip(Start_cut(2):End_cut(2));
AddBio_RHip3 = Right_Hip(Start_cut(3):End_cut(3));
AddBio_RHip4 = Right_Hip(Start_cut(4):End_cut(4));
AddBio_RHip5 = Right_Hip(Start_cut(5):End_cut(5));
AddBio_RHip6 = Right_Hip(Start_cut(6):End_cut(6));
AddBio_RHip7 = Right_Hip(Start_cut(7):End_cut(7));
AddBio_RHip8 = Right_Hip(Start_cut(8):End_cut(8));
AddBio_RHip9 = Right_Hip(Start_cut(9):End_cut(9));
AddBio_RHip10 = Right_Hip(Start_cut(10):End_cut(10));

AddBio_RAnkle1 = Right_Ankle(Start_cut(1):End_cut(1));
AddBio_RAnkle2 = Right_Ankle(Start_cut(2):End_cut(2));
AddBio_RAnkle3 = Right_Ankle(Start_cut(3):End_cut(3));
AddBio_RAnkle4 = Right_Ankle(Start_cut(4):End_cut(4));
AddBio_RAnkle5 = Right_Ankle(Start_cut(5):End_cut(5));
AddBio_RAnkle6 = Right_Ankle(Start_cut(6):End_cut(6));
AddBio_RAnkle7 = Right_Ankle(Start_cut(7):End_cut(7));
AddBio_RAnkle8 = Right_Ankle(Start_cut(8):End_cut(8));
AddBio_RAnkle9 = Right_Ankle(Start_cut(9):End_cut(9));
AddBio_RAnkle10 = Right_Ankle(Start_cut(10):End_cut(10));

AddBio_LKnee_size1 = size(AddBio_LKnee1);
AddBio_LKnee_size2 = size(AddBio_LKnee2);
AddBio_LKnee_size3 = size(AddBio_LKnee3);
AddBio_LKnee_size4 = size(AddBio_LKnee4);
AddBio_LKnee_size5 = size(AddBio_LKnee5);
AddBio_LKnee_size6 = size(AddBio_LKnee6);
AddBio_LKnee_size7 = size(AddBio_LKnee7);
AddBio_LKnee_size8 = size(AddBio_LKnee8);
AddBio_LKnee_size9 = size(AddBio_LKnee9);
AddBio_LKnee_size10 = size(AddBio_LKnee10);
AddBio_MaxSize = max([AddBio_LKnee_size1,AddBio_LKnee_size2,AddBio_LKnee_size3,AddBio_LKnee_size4,AddBio_LKnee_size5,...
    AddBio_LKnee_size6,AddBio_LKnee_size7,AddBio_LKnee_size8,AddBio_LKnee_size9,AddBio_LKnee_size10]);
AddBio_MinSize = min([AddBio_LKnee_size1(1),AddBio_LKnee_size2(1),AddBio_LKnee_size3(1),AddBio_LKnee_size4(1),AddBio_LKnee_size5(1),...
    AddBio_LKnee_size6(1),AddBio_LKnee_size7(1),AddBio_LKnee_size8(1),AddBio_LKnee_size9(1),AddBio_LKnee_size10(1)]);
AddBio_Matrix_Knee(AddBio_MaxSize,10) = 0;
AddBio_Matrix_Hip(AddBio_MaxSize,10) = 0;
AddBio_Matrix_Ankle(AddBio_MaxSize,10) = 0;

for i = 1:1:AddBio_LKnee_size1(1)
    AddBio_Matrix_Knee(i,1) = AddBio_LKnee1(i);
    AddBio_Matrix_Hip(i,1) = AddBio_RHip1(i);
    AddBio_Matrix_Ankle(i,1) = AddBio_RAnkle1(i);
end
for i = 1:1:AddBio_LKnee_size2(1)
    AddBio_Matrix_Knee(i,2) = AddBio_LKnee2(i);
    AddBio_Matrix_Hip(i,2) = AddBio_RHip2(i);
    AddBio_Matrix_Ankle(i,2) = AddBio_RAnkle2(i);
end
for i = 1:1:AddBio_LKnee_size3(1)
    AddBio_Matrix_Knee(i,3) = AddBio_LKnee3(i);
    AddBio_Matrix_Hip(i,3) = AddBio_RHip3(i);
    AddBio_Matrix_Ankle(i,3) = AddBio_RAnkle3(i);
end
for i = 1:1:AddBio_LKnee_size4(1)
    AddBio_Matrix_Knee(i,4) = AddBio_LKnee4(i);
    AddBio_Matrix_Hip(i,4) = AddBio_RHip4(i);
    AddBio_Matrix_Ankle(i,4) = AddBio_RAnkle4(i);
end
for i = 1:1:AddBio_LKnee_size5(1)
    AddBio_Matrix_Knee(i,5) = AddBio_LKnee5(i);
    AddBio_Matrix_Hip(i,5) = AddBio_RHip5(i);
    AddBio_Matrix_Ankle(i,5) = AddBio_RAnkle5(i);
end
for i = 1:1:AddBio_LKnee_size6(1)
    AddBio_Matrix_Knee(i,6) = AddBio_LKnee6(i);
    AddBio_Matrix_Hip(i,6) = AddBio_RHip6(i);
    AddBio_Matrix_Ankle(i,6) = AddBio_RAnkle6(i);
end
for i = 1:1:AddBio_LKnee_size7(1)
    AddBio_Matrix_Knee(i,7) = AddBio_LKnee7(i);
    AddBio_Matrix_Hip(i,7) = AddBio_RHip7(i);
    AddBio_Matrix_Ankle(i,7) = AddBio_RAnkle7(i);
end
for i = 1:1:AddBio_LKnee_size8(1)
    AddBio_Matrix_Knee(i,8) = AddBio_LKnee8(i);
    AddBio_Matrix_Hip(i,8) = AddBio_RHip8(i);
    AddBio_Matrix_Ankle(i,8) = AddBio_RAnkle8(i);
end
for i = 1:1:AddBio_LKnee_size9(1)
    AddBio_Matrix_Knee(i,9) = AddBio_LKnee9(i);
    AddBio_Matrix_Hip(i,9) = AddBio_RHip9(i);
    AddBio_Matrix_Ankle(i,9) = AddBio_RAnkle9(i);
end
for i = 1:1:AddBio_LKnee_size10(1)
    AddBio_Matrix_Knee(i,10) = AddBio_LKnee10(i);
    AddBio_Matrix_Hip(i,10) = AddBio_RHip10(i);
    AddBio_Matrix_Ankle(i,10) = AddBio_RAnkle10(i);
end
AddBio_Matrix_Knee = AddBio_Matrix_Knee(1:1:AddBio_MinSize,:);

zeroMask = (AddBio_Matrix_Knee == 0);
numZeros = sum(zeroMask(:));



Matrix_time_AddBio = AddBio_Cut_Time1-AddBio_Cut_Time1(1);
AddBio_Mean = zeros(1,size(AddBio_Matrix_Knee,1)).';
AddBio_Std = zeros(1,size(AddBio_Matrix_Knee,1)).';
for ii = 1:size(AddBio_Matrix_Knee,1)        
    AddBio_mean_calc = size(nonzeros(AddBio_Matrix_Knee(ii,:)));
    AddBio_Std_calc = size(nonzeros(AddBio_Matrix_Knee(ii,:)));

    AddBio_Mean(ii,1) = mean(nonzeros(AddBio_Matrix_Knee(ii,:)));
    AddBio_Std(ii,1) = std(nonzeros(AddBio_Matrix_Knee(ii,:)));
end
AddBio_Upper = AddBio_Mean + AddBio_Std*2;
AddBio_Lower = AddBio_Mean - AddBio_Std*2;
Matrix_time_plot_AddBio = Matrix_time_AddBio(1:1:AddBio_MinSize,:);
AddBio_Upper_plot = AddBio_Upper(1:1:AddBio_MinSize,:);
AddBio_Lower_plot = AddBio_Lower(1:1:AddBio_MinSize,:);
AddBio_LowerNeg_plot = zeros(length(AddBio_Lower_plot),1);
for i = 1:length(AddBio_Lower_plot)
    if AddBio_Lower_plot(i) < 0 
        AddBio_LowerNeg_plot(i) = AddBio_Lower_plot(i);
    else
        AddBio_LowerNeg_plot(i) = 0;
    end
end
[row,col] = find(AddBio_LowerNeg_plot<0);
AddBio_LowerNeg_plot = AddBio_LowerNeg_plot(1:1:length(row));
Matrix_timeNeg_plot_AddBio = Matrix_time_plot_AddBio(1:1:length(row));

% Cut IMU
IMU_AccRes = sqrt(IMU_AccCut_Data_x.^2 + IMU_AccCut_Data_y.^2 + IMU_AccCut_Data_z.^2);

IMU_Start_Cut_times = [];
IMU_End_Cut_times = [];
for i = 1:1:length(Start_cut)
    IMU_Start_Cut_times(i) = find(IMU_Cut_Time == AddBio_Time(Start_cut(i)));
    IMU_End_Cut_times(i) = find(IMU_Cut_Time == AddBio_Time(End_cut(i)));
end

IMU_Cut_Time1 = IMU_Cut_Time(IMU_Start_Cut_times(1):IMU_End_Cut_times(1));
IMU_Wave1 = IMU_AccRes(IMU_Start_Cut_times(1):IMU_End_Cut_times(1));
IMU_Cut_Time2 = IMU_Cut_Time(IMU_Start_Cut_times(2):IMU_End_Cut_times(2));
IMU_Wave2 = IMU_AccRes(IMU_Start_Cut_times(2):IMU_End_Cut_times(2));
IMU_Cut_Time3 = IMU_Cut_Time(IMU_Start_Cut_times(3):IMU_End_Cut_times(3));
IMU_Wave3 = IMU_AccRes(IMU_Start_Cut_times(3):IMU_End_Cut_times(3));
IMU_Cut_Time4 = IMU_Cut_Time(IMU_Start_Cut_times(4):IMU_End_Cut_times(4));
IMU_Wave4 = IMU_AccRes(IMU_Start_Cut_times(4):IMU_End_Cut_times(4));
IMU_Cut_Time5 = IMU_Cut_Time(IMU_Start_Cut_times(5):IMU_End_Cut_times(5));
IMU_Wave5 = IMU_AccRes(IMU_Start_Cut_times(5):IMU_End_Cut_times(5));
IMU_Cut_Time6 = IMU_Cut_Time(IMU_Start_Cut_times(6):IMU_End_Cut_times(6));
IMU_Wave6 = IMU_AccRes(IMU_Start_Cut_times(6):IMU_End_Cut_times(6));
IMU_Cut_Time7 = IMU_Cut_Time(IMU_Start_Cut_times(7):IMU_End_Cut_times(7));
IMU_Wave7 = IMU_AccRes(IMU_Start_Cut_times(7):IMU_End_Cut_times(7));
IMU_Cut_Time8 = IMU_Cut_Time(IMU_Start_Cut_times(8):IMU_End_Cut_times(8));
IMU_Wave8 = IMU_AccRes(IMU_Start_Cut_times(8):IMU_End_Cut_times(8));
IMU_Cut_Time9 = IMU_Cut_Time(IMU_Start_Cut_times(9):IMU_End_Cut_times(9));
IMU_Wave9 = IMU_AccRes(IMU_Start_Cut_times(9):IMU_End_Cut_times(9));
IMU_Cut_Time10 = IMU_Cut_Time(IMU_Start_Cut_times(10):IMU_End_Cut_times(10));
IMU_Wave10 = IMU_AccRes(IMU_Start_Cut_times(10):IMU_End_Cut_times(10));

figure
hold on

plot(IMU_Cut_Time,IMU_AccCut_Data_x,'k-',IMU_Cut_Time,IMU_AccCut_Data_y,'r-', ...
    IMU_Cut_Time,IMU_AccCut_Data_z,'b-')

plot(IMU_Cut_Time1,IMU_Wave1,'m-',IMU_Cut_Time2,IMU_Wave2,'m-', ...
    IMU_Cut_Time3,IMU_Wave3,'m-',IMU_Cut_Time4,IMU_Wave4,'m-', ...
    IMU_Cut_Time5,IMU_Wave5,'m-',IMU_Cut_Time6,IMU_Wave6,'m-', ...
    IMU_Cut_Time7,IMU_Wave7,'m-',IMU_Cut_Time7,IMU_Wave7,'m-', ...
    IMU_Cut_Time8,IMU_Wave8,'m-',IMU_Cut_Time9,IMU_Wave9,'m-', ...
    IMU_Cut_Time10,IMU_Wave10,'m-')

xline(IMU_Cut_Time(IMU_Start_Cut_times),'--r')
xline(IMU_Cut_Time(IMU_End_Cut_times),'--b')
xlabel('Time [s]')
%xlim([0,6])
ylabel('Acc [g]')
title('IMU Acceleration')
legend('IMU AccX','IMU AccY','IMU AccZ','Location','eastoutside')
hold off

IMU_size1 = size(IMU_Wave1);
IMU_size2 = size(IMU_Wave2);
IMU_size3 = size(IMU_Wave3);
IMU_size4 = size(IMU_Wave4);
IMU_size5 = size(IMU_Wave5);
IMU_size6 = size(IMU_Wave6);
IMU_size7 = size(IMU_Wave7);
IMU_size8 = size(IMU_Wave8);
IMU_size9 = size(IMU_Wave9);
IMU_size10 = size(IMU_Wave10);
IMU_MaxSize = max([IMU_size1,IMU_size2,IMU_size3,IMU_size4,IMU_size5,...
    IMU_size6,IMU_size7,IMU_size8,IMU_size9,IMU_size10]);
IMU_MinSize = min([IMU_size1(1),IMU_size2(1),IMU_size3(1),IMU_size4(1),IMU_size5(1),...
    IMU_size6(1),IMU_size7(1),IMU_size8(1),IMU_size9(1),IMU_size10(1)]);
IMU_Matrix(IMU_MaxSize,10) = 0;

for i = 1:1:IMU_size1(1)
    IMU_Matrix(i,1) = IMU_Wave1(i);
end
for i = 1:1:IMU_size2(1)
    IMU_Matrix(i,2) = IMU_Wave2(i);
end
for i = 1:1:IMU_size3(1)
    IMU_Matrix(i,3) = IMU_Wave3(i);
end
for i = 1:1:IMU_size4(1)
    IMU_Matrix(i,4) = IMU_Wave4(i);
end
for i = 1:1:IMU_size5(1)
    IMU_Matrix(i,5) = IMU_Wave5(i);
end
for i = 1:1:IMU_size6(1)
    IMU_Matrix(i,6) = IMU_Wave6(i);
end
for i = 1:1:IMU_size7(1)
    IMU_Matrix(i,7) = IMU_Wave7(i);
end
for i = 1:1:IMU_size8(1)
    IMU_Matrix(i,8) = IMU_Wave8(i);
end
for i = 1:1:IMU_size9(1)
    IMU_Matrix(i,9) = IMU_Wave9(i);
end
for i = 1:1:IMU_size10(1)
    IMU_Matrix(i,10) = IMU_Wave10(i);
end
IMU_Matrix = IMU_Matrix(1:1:IMU_MinSize,:);

zeroMask = (IMU_Matrix == 0);
numZeros = sum(zeroMask(:));


Matrix_time_IMU = IMU_Cut_Time1-IMU_Cut_Time1(1);
IMU_Mean = zeros(1,size(IMU_Matrix,1)).';
IMU_Std = zeros(1,size(IMU_Matrix,1)).';
for ii = 1:size(IMU_Matrix,1)        
    IMU_mean_calc = size(nonzeros(IMU_Matrix(ii,:)));
    IMU_Std_calc = size(nonzeros(IMU_Matrix(ii,:)));

    IMU_Mean(ii,1) = mean(nonzeros(IMU_Matrix(ii,:)));
    IMU_Std(ii,1) = std(nonzeros(IMU_Matrix(ii,:)));
end
IMU_Upper = IMU_Mean + IMU_Std*2;
IMU_Lower = IMU_Mean - IMU_Std*2;
Matrix_time_IMU = Matrix_time_IMU';
Matrix_time_plot_IMU = Matrix_time_IMU(1:1:IMU_MinSize,:);
IMU_Upper_plot = IMU_Upper(1:1:IMU_MinSize,:);
IMU_Lower_plot = IMU_Lower(1:1:IMU_MinSize,:);
IMU_LowerNeg_plot = zeros(length(IMU_Lower_plot),1);
for i = 1:length(IMU_Lower_plot)
    if IMU_Lower_plot(i) < 0 
        IMU_LowerNeg_plot(i) = IMU_Lower_plot(i);
    else
        IMU_LowerNeg_plot(i) = 0;
    end
end
[row,col] = find(IMU_LowerNeg_plot<0);
IMU_LowerNeg_plot = IMU_LowerNeg_plot(1:1:length(row));
Matrix_timeNeg_plot_IMU = Matrix_time_plot_IMU(1:1:length(row));

% Cut dSpace
dSpace_Start_Cut_times = [];
dSpace_End_Cut_times = [];
for i = 1:1:length(Start_cut)
    dSpace_Start_Cut_times(i) = find(LC_Cut_Time == AddBio_Time(Start_cut(i)));
    dSpace_End_Cut_times(i) = find(LC_Cut_Time == AddBio_Time(End_cut(i)));
end

LC_Cut_Wavelength_time1 = LC_Cut_Time(dSpace_Start_Cut_times(1):dSpace_End_Cut_times(1));
Fx1_Cut_Wavelength1 = Fx1_Cut_Data(dSpace_Start_Cut_times(1):dSpace_End_Cut_times(1));
LC_Cut_Wavelength_time2 = LC_Cut_Time(dSpace_Start_Cut_times(2):dSpace_End_Cut_times(2));
Fx1_Cut_Wavelength2 = Fx1_Cut_Data(dSpace_Start_Cut_times(2):dSpace_End_Cut_times(2));
LC_Cut_Wavelength_time3 = LC_Cut_Time(dSpace_Start_Cut_times(3):dSpace_End_Cut_times(3));
Fx1_Cut_Wavelength3 = Fx1_Cut_Data(dSpace_Start_Cut_times(3):dSpace_End_Cut_times(3));
LC_Cut_Wavelength_time4 = LC_Cut_Time(dSpace_Start_Cut_times(4):dSpace_End_Cut_times(4));
Fx1_Cut_Wavelength4 = Fx1_Cut_Data(dSpace_Start_Cut_times(4):dSpace_End_Cut_times(4));
LC_Cut_Wavelength_time5 = LC_Cut_Time(dSpace_Start_Cut_times(5):dSpace_End_Cut_times(5));
Fx1_Cut_Wavelength5 = Fx1_Cut_Data(dSpace_Start_Cut_times(5):dSpace_End_Cut_times(5));
LC_Cut_Wavelength_time6 = LC_Cut_Time(dSpace_Start_Cut_times(6):dSpace_End_Cut_times(6));
Fx1_Cut_Wavelength6 = Fx1_Cut_Data(dSpace_Start_Cut_times(6):dSpace_End_Cut_times(6));
LC_Cut_Wavelength_time7 = LC_Cut_Time(dSpace_Start_Cut_times(7):dSpace_End_Cut_times(7));
Fx1_Cut_Wavelength7 = Fx1_Cut_Data(dSpace_Start_Cut_times(7):dSpace_End_Cut_times(7));
LC_Cut_Wavelength_time8 = LC_Cut_Time(dSpace_Start_Cut_times(8):dSpace_End_Cut_times(8));
Fx1_Cut_Wavelength8 = Fx1_Cut_Data(dSpace_Start_Cut_times(8):dSpace_End_Cut_times(8));
LC_Cut_Wavelength_time9 = LC_Cut_Time(dSpace_Start_Cut_times(9):dSpace_End_Cut_times(9));
Fx1_Cut_Wavelength9 = Fx1_Cut_Data(dSpace_Start_Cut_times(9):dSpace_End_Cut_times(9));
LC_Cut_Wavelength_time10 = LC_Cut_Time(dSpace_Start_Cut_times(10):dSpace_End_Cut_times(10));
Fx1_Cut_Wavelength10 = Fx1_Cut_Data(dSpace_Start_Cut_times(10):dSpace_End_Cut_times(10));

% figure
% hold on
% plot(LC_Cut_Wavelength_time1-LC_Cut_Wavelength_time1(1),Fx1_Cut_Wavelength1)
% plot(LC_Cut_Wavelength_time2-LC_Cut_Wavelength_time2(1),Fx1_Cut_Wavelength2)
% plot(LC_Cut_Wavelength_time3-LC_Cut_Wavelength_time3(1),Fx1_Cut_Wavelength3)
% plot(LC_Cut_Wavelength_time4-LC_Cut_Wavelength_time4(1),Fx1_Cut_Wavelength4)
% plot(LC_Cut_Wavelength_time5-LC_Cut_Wavelength_time5(1),Fx1_Cut_Wavelength5)
% plot(LC_Cut_Wavelength_time6-LC_Cut_Wavelength_time6(1),Fx1_Cut_Wavelength6)
% plot(LC_Cut_Wavelength_time7-LC_Cut_Wavelength_time7(1),Fx1_Cut_Wavelength7)
% plot(LC_Cut_Wavelength_time8-LC_Cut_Wavelength_time8(1),Fx1_Cut_Wavelength8)
% plot(LC_Cut_Wavelength_time9-LC_Cut_Wavelength_time9(1),Fx1_Cut_Wavelength9)
% plot(LC_Cut_Wavelength_time10-LC_Cut_Wavelength_time10(1),Fx1_Cut_Wavelength10)
% hold off

% Combining Matrices
Fx1_size1 = size(Fx1_Cut_Wavelength1);
Fx1_size2 = size(Fx1_Cut_Wavelength2);
Fx1_size3 = size(Fx1_Cut_Wavelength3);
Fx1_size4 = size(Fx1_Cut_Wavelength4);
Fx1_size5 = size(Fx1_Cut_Wavelength5);
Fx1_size6 = size(Fx1_Cut_Wavelength6);
Fx1_size7 = size(Fx1_Cut_Wavelength7);
Fx1_size8 = size(Fx1_Cut_Wavelength8);
Fx1_size9 = size(Fx1_Cut_Wavelength9);
Fx1_size10 = size(Fx1_Cut_Wavelength10);
Fx1_MaxSize = max([Fx1_size1,Fx1_size2,Fx1_size3,Fx1_size4,Fx1_size5,...
    Fx1_size6,Fx1_size7,Fx1_size8,Fx1_size9,Fx1_size10]);
Fx1_MinSize = min([Fx1_size1(2),Fx1_size2(2),Fx1_size3(2),Fx1_size4(2),Fx1_size5(2),...
    Fx1_size6(2),Fx1_size7(2),Fx1_size8(2),Fx1_size9(2),Fx1_size10(2)]);
Fx1_Matrix(Fx1_MaxSize,10) = 0;

for i = 1:1:Fx1_size1(2)
    Fx1_Matrix(i,1) = Fx1_Cut_Wavelength1(i);
end
for i = 1:1:Fx1_size2(2)
    Fx1_Matrix(i,2) = Fx1_Cut_Wavelength2(i);
end
for i = 1:1:Fx1_size3(2)
    Fx1_Matrix(i,3) = Fx1_Cut_Wavelength3(i);
end
for i = 1:1:Fx1_size4(2)
    Fx1_Matrix(i,4) = Fx1_Cut_Wavelength4(i);
end
for i = 1:1:Fx1_size5(2)
    Fx1_Matrix(i,5) = Fx1_Cut_Wavelength5(i);
end
for i = 1:1:Fx1_size6(2)
    Fx1_Matrix(i,6) = Fx1_Cut_Wavelength6(i);
end
for i = 1:1:Fx1_size7(2)
    Fx1_Matrix(i,7) = Fx1_Cut_Wavelength7(i);
end
for i = 1:1:Fx1_size8(2)
end
    Fx1_Matrix(i,8) = Fx1_Cut_Wavelength8(i);
for i = 1:1:Fx1_size9(2)
    Fx1_Matrix(i,9) = Fx1_Cut_Wavelength9(i);
end
for i = 1:1:Fx1_size10(2)
    Fx1_Matrix(i,10) = Fx1_Cut_Wavelength10(i);
end
Fx1_Matrix = Fx1_Matrix(1:1:Fx1_MinSize,:);

Matrix_time_Fx1 = LC_Cut_Wavelength_time1-LC_Cut_Wavelength_time1(1);
Fx1_Mean = zeros(1,size(Fx1_Matrix,1)).';
Fx1_Std = zeros(1,size(Fx1_Matrix,1)).';
for i = 1:size(Fx1_Matrix,1)        
        Fx1_Mean(i,1) = mean(nonzeros(Fx1_Matrix(i,:)));
        Fx1_Std(i,1) = std(nonzeros(Fx1_Matrix(i,:)));
end
Fx1_Upper = Fx1_Mean + Fx1_Std*2;
Fx1_Lower = Fx1_Mean - Fx1_Std*2;
Matrix_time_plot_Fx1 = Matrix_time_Fx1(1:1:Fx1_MinSize,:);
Fx1_Upper_plot = Fx1_Upper(1:1:Fx1_MinSize,:);
Fx1_Lower_plot = Fx1_Lower(1:1:Fx1_MinSize,:);
Fx1_LowerNeg_plot = zeros(length(Fx1_Lower_plot),1);
for i = 1:length(Fx1_Lower_plot)
    if Fx1_Lower_plot(i) < 0 
        Fx1_LowerNeg_plot(i) = Fx1_Lower_plot(i);
    else
        Fx1_LowerNeg_plot(i) = 0;
    end
end
[row,col] = find(Fx1_LowerNeg_plot<0);
Fx1_LowerNeg_plot = Fx1_LowerNeg_plot(1:1:length(row));
Matrix_timeNeg_plot_Fx1 = Matrix_time_plot_Fx1(1:1:length(row));

% dSpace cut Torque
dSpace_Start_Cut_times = [];
dSpace_End_Cut_times = [];
for i = 1:1:length(Start_cut)
    dSpace_Start_Cut_times(i) = find(LC_Cut_Time == AddBio_Time(Start_cut(i)));
    dSpace_End_Cut_times(i) = find(LC_Cut_Time == AddBio_Time(End_cut(i)));
end

LC_Cut_Wavelength_time1 = LC_Cut_Time(dSpace_Start_Cut_times(1):dSpace_End_Cut_times(1));
Tx1_Cut_Wavelength1 = Tx1_Cut_Data(dSpace_Start_Cut_times(1):dSpace_End_Cut_times(1));
LC_Cut_Wavelength_time2 = LC_Cut_Time(dSpace_Start_Cut_times(2):dSpace_End_Cut_times(2));
Tx1_Cut_Wavelength2 = Tx1_Cut_Data(dSpace_Start_Cut_times(2):dSpace_End_Cut_times(2));
LC_Cut_Wavelength_time3 = LC_Cut_Time(dSpace_Start_Cut_times(3):dSpace_End_Cut_times(3));
Tx1_Cut_Wavelength3 = Tx1_Cut_Data(dSpace_Start_Cut_times(3):dSpace_End_Cut_times(3));
LC_Cut_Wavelength_time4 = LC_Cut_Time(dSpace_Start_Cut_times(4):dSpace_End_Cut_times(4));
Tx1_Cut_Wavelength4 = Tx1_Cut_Data(dSpace_Start_Cut_times(4):dSpace_End_Cut_times(4));
LC_Cut_Wavelength_time5 = LC_Cut_Time(dSpace_Start_Cut_times(5):dSpace_End_Cut_times(5));
Tx1_Cut_Wavelength5 = Tx1_Cut_Data(dSpace_Start_Cut_times(5):dSpace_End_Cut_times(5));
LC_Cut_Wavelength_time6 = LC_Cut_Time(dSpace_Start_Cut_times(6):dSpace_End_Cut_times(6));
Tx1_Cut_Wavelength6 = Tx1_Cut_Data(dSpace_Start_Cut_times(6):dSpace_End_Cut_times(6));
LC_Cut_Wavelength_time7 = LC_Cut_Time(dSpace_Start_Cut_times(7):dSpace_End_Cut_times(7));
Tx1_Cut_Wavelength7 = Tx1_Cut_Data(dSpace_Start_Cut_times(7):dSpace_End_Cut_times(7));
LC_Cut_Wavelength_time8 = LC_Cut_Time(dSpace_Start_Cut_times(8):dSpace_End_Cut_times(8));
Tx1_Cut_Wavelength8 = Tx1_Cut_Data(dSpace_Start_Cut_times(8):dSpace_End_Cut_times(8));
LC_Cut_Wavelength_time9 = LC_Cut_Time(dSpace_Start_Cut_times(9):dSpace_End_Cut_times(9));
Tx1_Cut_Wavelength9 = Tx1_Cut_Data(dSpace_Start_Cut_times(9):dSpace_End_Cut_times(9));
LC_Cut_Wavelength_time10 = LC_Cut_Time(dSpace_Start_Cut_times(10):dSpace_End_Cut_times(10));
Tx1_Cut_Wavelength10 = Tx1_Cut_Data(dSpace_Start_Cut_times(10):dSpace_End_Cut_times(10));

% figure
% hold on
% plot(LC_Cut_Wavelength_time1-LC_Cut_Wavelength_time1(1),Tx1_Cut_Wavelength1)
% plot(LC_Cut_Wavelength_time2-LC_Cut_Wavelength_time2(1),Tx1_Cut_Wavelength2)
% plot(LC_Cut_Wavelength_time3-LC_Cut_Wavelength_time3(1),Tx1_Cut_Wavelength3)
% plot(LC_Cut_Wavelength_time4-LC_Cut_Wavelength_time4(1),Tx1_Cut_Wavelength4)
% plot(LC_Cut_Wavelength_time5-LC_Cut_Wavelength_time5(1),Tx1_Cut_Wavelength5)
% plot(LC_Cut_Wavelength_time6-LC_Cut_Wavelength_time6(1),Tx1_Cut_Wavelength6)
% plot(LC_Cut_Wavelength_time7-LC_Cut_Wavelength_time7(1),Tx1_Cut_Wavelength7)
% plot(LC_Cut_Wavelength_time8-LC_Cut_Wavelength_time8(1),Tx1_Cut_Wavelength8)
% plot(LC_Cut_Wavelength_time9-LC_Cut_Wavelength_time9(1),Tx1_Cut_Wavelength9)
% plot(LC_Cut_Wavelength_time10-LC_Cut_Wavelength_time10(1),Tx1_Cut_Wavelength10)
% hold off

% Combining Matrices
Tx1_size1 = size(Tx1_Cut_Wavelength1);
Tx1_size2 = size(Tx1_Cut_Wavelength2);
Tx1_size3 = size(Tx1_Cut_Wavelength3);
Tx1_size4 = size(Tx1_Cut_Wavelength4);
Tx1_size5 = size(Tx1_Cut_Wavelength5);
Tx1_size6 = size(Tx1_Cut_Wavelength6);
Tx1_size7 = size(Tx1_Cut_Wavelength7);
Tx1_size8 = size(Tx1_Cut_Wavelength8);
Tx1_size9 = size(Tx1_Cut_Wavelength9);
Tx1_size10 = size(Tx1_Cut_Wavelength10);
Tx1_MaxSize = max([Tx1_size1,Tx1_size2,Tx1_size3,Tx1_size4,Tx1_size5,...
    Tx1_size6,Tx1_size7,Tx1_size8,Tx1_size9,Tx1_size10]);
Tx1_MinSize = min([Tx1_size1(2),Tx1_size2(2),Tx1_size3(2),Tx1_size4(2),Tx1_size5(2),...
    Tx1_size6(2),Tx1_size7(2),Tx1_size8(2),Tx1_size9(2),Tx1_size10(2)]);
Tx1_Matrix(Tx1_MaxSize,10) = 0;

for i = 1:1:Tx1_size1(2)
    Tx1_Matrix(i,1) = Tx1_Cut_Wavelength1(i);
end
for i = 1:1:Tx1_size2(2)
    Tx1_Matrix(i,2) = Tx1_Cut_Wavelength2(i);
end
for i = 1:1:Tx1_size3(2)
    Tx1_Matrix(i,3) = Tx1_Cut_Wavelength3(i);
end
for i = 1:1:Tx1_size4(2)
    Tx1_Matrix(i,4) = Tx1_Cut_Wavelength4(i);
end
for i = 1:1:Tx1_size5(2)
    Tx1_Matrix(i,5) = Tx1_Cut_Wavelength5(i);
end
for i = 1:1:Tx1_size6(2)
    Tx1_Matrix(i,6) = Tx1_Cut_Wavelength6(i);
end
for i = 1:1:Tx1_size7(2)
    Tx1_Matrix(i,7) = Tx1_Cut_Wavelength7(i);
end
for i = 1:1:Tx1_size8(2)
end
    Tx1_Matrix(i,8) = Tx1_Cut_Wavelength8(i);
for i = 1:1:Tx1_size9(2)
    Tx1_Matrix(i,9) = Tx1_Cut_Wavelength9(i);
end
for i = 1:1:Tx1_size10(2)
    Tx1_Matrix(i,10) = Tx1_Cut_Wavelength10(i);
end
Tx1_Matrix = Tx1_Matrix(1:1:Tx1_MinSize,:);

Matrix_time_Tx1 = LC_Cut_Wavelength_time1-LC_Cut_Wavelength_time1(1);
Tx1_Mean = zeros(1,size(Tx1_Matrix,1)).';
Tx1_Std = zeros(1,size(Tx1_Matrix,1)).';
for i = 1:size(Tx1_Matrix,1)        
        Tx1_Mean(i,1) = mean(nonzeros(Tx1_Matrix(i,:)));
        Tx1_Std(i,1) = std(nonzeros(Tx1_Matrix(i,:)));
end
Tx1_Upper = Tx1_Mean + Tx1_Std*2;
Tx1_Lower = Tx1_Mean - Tx1_Std*2;
Matrix_time_plot_Tx1 = Matrix_time_Tx1(1:1:Tx1_MinSize,:);
Tx1_Upper_plot = Tx1_Upper(1:1:Tx1_MinSize,:);
Tx1_Lower_plot = Tx1_Lower(1:1:Tx1_MinSize,:);
Tx1_LowerNeg_plot = zeros(length(Tx1_Lower_plot),1);
for i = 1:length(Tx1_Lower_plot)
    if Tx1_Lower_plot(i) < 0 
        Tx1_LowerNeg_plot(i) = Tx1_Lower_plot(i);
    else
        Tx1_LowerNeg_plot(i) = 0;
    end
end
[row,col] = find(Tx1_LowerNeg_plot<0);
Tx1_LowerNeg_plot = Tx1_LowerNeg_plot(1:1:length(row));
Matrix_timeNeg_plot_Tx1 = Matrix_time_plot_Tx1(1:1:length(row));


AddBio_percent_vect = Matrix_time_AddBio(1:1:AddBio_MinSize,:);
AddBio_percent = AddBio_percent_vect/AddBio_percent_vect(end)*100;

End_Matrix = Matrix_time_AddBio(1:1:AddBio_MinSize,:);
End_Martix = End_Matrix(end);

%%
figure
ax1 = subplot(4,1,1);
hold on
area(Matrix_time_plot_AddBio, max(AddBio_Lower_plot,AddBio_Upper_plot),...
    'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
area(Matrix_time_plot_AddBio, min(AddBio_Lower_plot,AddBio_Upper_plot),...
    'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none');
area(Matrix_timeNeg_plot_AddBio, AddBio_LowerNeg_plot,'FaceColor', 'r',...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Matrix_time_AddBio(1:1:AddBio_MinSize,:), AddBio_Matrix_Knee,'Color',[0.5, 0.5, 0.5], 'LineWidth', 1)
plot(Matrix_time_AddBio(1:1:AddBio_MinSize,:), AddBio_Mean(1:1:AddBio_MinSize,:),'k-','LineWidth',3)
plot(Matrix_time_AddBio(1:1:AddBio_MinSize,:), AddBio_Upper(1:1:AddBio_MinSize,:),'r-','LineWidth',2)
plot(Matrix_time_AddBio(1:1:AddBio_MinSize,:), AddBio_Lower(1:1:AddBio_MinSize,:),'r-','LineWidth',2)
xline(0,'k-')
yline(0,'k-')
xlabel('Time [s]')
%xlim([0,Matrix_time_plot_AddBio(AddBio_MinSize)])
xlim([0,1.2])
ylabel('Postion [\circ]')%,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
title('AddBio Right Knee')
%legend('Fx','Load Cell Fy','Load Cell Fz','Location','eastoutside')
hold off

ax2 = subplot(4,1,2);
hold on
area(Matrix_time_plot_IMU, max(IMU_Lower_plot,IMU_Upper_plot),...
    'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
area(Matrix_time_plot_IMU, min(IMU_Lower_plot,IMU_Upper_plot),...
    'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none');
area(Matrix_timeNeg_plot_IMU, IMU_LowerNeg_plot,'FaceColor', 'r',...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Matrix_time_IMU(1:1:IMU_MinSize,:),IMU_Matrix,'Color',[0.5, 0.5, 0.5], 'LineWidth', 1)
plot(Matrix_time_IMU(1:1:IMU_MinSize,:),IMU_Mean(1:1:IMU_MinSize,:),'k-','LineWidth',3)
plot(Matrix_time_IMU(1:1:IMU_MinSize,:),IMU_Upper(1:1:IMU_MinSize,:),'r-','LineWidth',2)
plot(Matrix_time_IMU(1:1:IMU_MinSize,:),IMU_Lower(1:1:IMU_MinSize,:),'r-','LineWidth',2)
xline(0,'k-')
yline(0,'k-')
xlabel('Time [s]')
% xlim([0,Matrix_time_plot_IMU(IMU_MinSize)])
xlim([0,1.2])
ylabel('Accel [g]')%,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
title('IMU Acceleration')
%legend('Fx','Load Cell Fy','Load Cell Fz','Location','eastoutside')
hold off

ax3 = subplot(4,1,3);
hold on
area(Matrix_time_plot_Fx1, max(Fx1_Lower_plot,Fx1_Upper_plot),...
    'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
area(Matrix_time_plot_Fx1, min(Fx1_Lower_plot,Fx1_Upper_plot),...
    'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none');
area(Matrix_timeNeg_plot_Fx1, Fx1_LowerNeg_plot,'FaceColor', 'r',...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Matrix_time_Fx1(1:1:Fx1_MinSize,:),Fx1_Matrix,'Color',[0.5, 0.5, 0.5], 'LineWidth', 1)
plot(Matrix_time_Fx1(1:1:Fx1_MinSize,:),Fx1_Mean(1:1:Fx1_MinSize,:),'k-','LineWidth',3)
plot(Matrix_time_Fx1(1:1:Fx1_MinSize,:),Fx1_Upper(1:1:Fx1_MinSize,:),'r-','LineWidth',2)
plot(Matrix_time_Fx1(1:1:Fx1_MinSize,:),Fx1_Lower(1:1:Fx1_MinSize,:),'r-','LineWidth',2)
xline(0,'k-')
yline(0,'k-')
xlabel('Time [s]')
% xlim([0,Matrix_time_plot_Fx1(Fx1_MinSize)])
xlim([0,1.2])
ylabel('Force [N]')%,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
title('LC Force')
%legend('Fx','Load Cell Fy','Load Cell Fz','Location','eastoutside')
hold off

ax4 = subplot(4,1,4);
hold on
area(Matrix_time_plot_Tx1, max(Tx1_Lower_plot,Tx1_Upper_plot),...
    'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
area(Matrix_time_plot_Tx1, min(Tx1_Lower_plot,Tx1_Upper_plot),...
    'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'none');
area(Matrix_timeNeg_plot_Tx1, Tx1_LowerNeg_plot,'FaceColor', 'r',...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Matrix_time_Tx1(1:1:Tx1_MinSize,:),Tx1_Matrix,'Color',[0.5, 0.5, 0.5], 'LineWidth', 1)
plot(Matrix_time_Tx1(1:1:Tx1_MinSize,:),Tx1_Mean(1:1:Tx1_MinSize,:),'k-','LineWidth',3)
plot(Matrix_time_Tx1(1:1:Tx1_MinSize,:),Tx1_Upper(1:1:Tx1_MinSize,:),'r-','LineWidth',2)
plot(Matrix_time_Tx1(1:1:Tx1_MinSize,:),Tx1_Lower(1:1:Tx1_MinSize,:),'r-','LineWidth',2)
xline(0,'k-')
yline(0,'k-')
xlabel('Time [s]')
% xlim([0,Matrix_time_plot_Tx1(Tx1_MinSize)])
xlim([0,1.2])
ylabel('Torque [Nm]')%,'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
title('LC Torque')
%legend('Tx','Load Cell Fy','Load Cell Fz','Location','eastoutside')
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');


%% EXPORT
IMU_labels = {'Time','Acc x','Acc y','Acc z','Gyro x','Gyro y','Gyro z'};
IMU_Export = [IMU_Cut_Time' IMU_AccCut_Data_x IMU_AccCut_Data_y IMU_AccCut_Data_z...
    IMU_GyroCut_Data_x IMU_GyroCut_Data_y IMU_GyroCut_Data_z];
IMU_output = [IMU_labels; num2cell(IMU_Export)];
writecell(IMU_output,'IMU TEST.csv')

LC_labels = {'Time','Force x','Force y','Force z','Torque x','Torque y','Torque z'};
LC_Export = [LC_Cut_Time Fx1_Cut_Data' Fy1_Cut_Data' Fz1_Cut_Data' Tx1_Cut_Data'...
    Ty1_Cut_Data' Tz1_Cut_Data'];
LC_output = [LC_labels; num2cell(LC_Export)];
writecell(LC_output,'LC TEST.csv')

AddBio_labels = AddBio_Data.textdata(end,:);
AddBio_Export = AddBio_Data.data;
AddBio_output = [AddBio_labels; num2cell(AddBio_Export)];
writecell(AddBio_output,'AddBio TEST.csv')



%% ROM calcs

Knee_Data_Load = load('RightKneeROM.mat');
Hip_Data_Load = load('RightHipROM.mat');
Ankle_Data_Load = load('RightAnkleROM.mat');

Knee_Data = Knee_Data_Load.AddBio_Matrix;
Hip_Data = Hip_Data_Load.AddBio_Matrix_Hip(1:length(Knee_Data),:);
Ankle_Data = Ankle_Data_Load.AddBio_Matrix_Ankle(1:length(Knee_Data),:);

Knee_max_columns = zeros(1,width(Knee_Data));
Knee_min_columns = zeros(1,width(Knee_Data));
ROM_Knee = zeros(1,width(Knee_Data));
Hip_max_columns = zeros(1,width(Hip_Data));
Hip_min_columns = zeros(1,width(Hip_Data));
ROM_Hip = zeros(1,width(Hip_Data));
Ankle_max_columns = zeros(1,width(Ankle_Data));
Ankle_min_columns = zeros(1,width(Ankle_Data));
ROM_Ankle = zeros(1,width(Ankle_Data));
for i = 1:1:width(Knee_Data)
    Knee_max_columns(i) = max(Knee_Data(:,i));
    Knee_min_columns(i) = min(Knee_Data(:,i));
    Hip_max_columns(i) = max(Hip_Data(:,i));
    Hip_min_columns(i) = min(Hip_Data(:,i));
    Ankle_max_columns(i) = max(Ankle_Data(:,i));
    Ankle_min_columns(i) = min(Ankle_Data(:,i));
    ROM_Knee(i) = Knee_max_columns(i) - Knee_min_columns(i);
    ROM_Hip(i) = Hip_max_columns(i) - Hip_min_columns(i);
    ROM_Ankle(i) = Ankle_max_columns(i) - Ankle_min_columns(i);
end

Knee_max_Mean = mean(Knee_max_columns)
Knee_max_STD = std(Knee_max_columns)
Knee_min_Mean = mean(Knee_min_columns)
Knee_min_STD = std(Knee_min_columns)
Knee_Rom_Mean = mean(ROM_Knee)
Knee_Rom_STD = std(ROM_Knee)

Hip_max_Mean = mean(Hip_max_columns)
Hip_max_STD = std(Hip_max_columns)
Hip_min_Mean = mean(Hip_min_columns)
Hip_min_STD = std(Hip_min_columns)
Hip_Rom_Mean = mean(ROM_Hip)
Hip_Rom_STD = std(ROM_Hip)

Ankle_max_Mean = mean(Ankle_max_columns)
Ankle_max_STD = std(Ankle_max_columns)
Ankle_min_Mean = mean(Ankle_min_columns)
Ankle_min_STD = std(Ankle_min_columns)
Ankle_Rom_Mean = mean(ROM_Ankle)
Ankle_Rom_STD = std(ROM_Ankle)






%% Helper Methods
%%%%%%%%%%%%%%%%%%%%
% Jaden Hunter
% November 22, 2024
%%%%%%%%%%%%%%%%%%%%



% input: user selected *.mot file(s)
% return: a 1x1 structure containing all files in the order as selected

function arr = motFilesToArray()
       arr = [];
    [filenames, folderPath] = uigetfile("MultiSelect","on", "*.mot");

   if(iscell(filenames))
        arr = importdata("" + folderPath + filenames(1));
        for i = 2:length(filenames)
           if(cell2mat(filenames(i)) ~= 0)
               filePath = "" + folderPath + filenames(i);
               arr = [importdata(filePath); arr];
           else
               error("error during csvread: file at input" + i + "was not input properly");
           end

        end
   else
       if(filenames ~= 0)
       arr = [arr, importdata("" + folderPath + filenames)];
       else
           error("input required");
       end
   end
end

% input: user selected *.csv file(s)
% input: C - parameter for csvread() (row offset R)
% input R - parameter for csvread() (column offset C)
% return: a matrix containing all files in order as selected

function arr = csvFilesToAray(C, R)
    arr = [];
    [filenames, folderPath] = uigetfile("MultiSelect","on", "*.csv");

   if(iscell(filenames))
        arr = csvread("" + folderPath + filenames(1), C, R);
        for i = 2:length(filenames)
           if(cell2mat(filenames(i)) ~= 0)
               filePath = "" + folderPath + filenames(i);
               arr = [csvread(filePath, C, R); arr];
           else
               error("error during csvread: file at input " + i + " was not recieved");
           end

        end
   else
       if(filenames ~= 0)
            arr = [arr, csvread("" + folderPath + filenames, C, R)];
       else
           error("input required");
       end
   end
end


function [Devices_VICON,  Trajectories_VICON, Trajectories_Cut, AddBio_Data, DSpace_Data, IMU_Data] = importExperimentData(VICON_C, VICON_R, IMU_C, IMU_R)
disp("IMPORTANT:")
disp("1. if there are multiple files, please select files in descending numerical order (e.g. 3, 2, 1, 0)");
disp("2. ensure all necessary files are in the same folder");
disp(" ")
disp("select Data for VICON Devices");
Devices_VICON = csvFilesToAray(VICON_C,VICON_R);

disp("select VICON Trajectories");
Trajectories_VICON = csvFilesToAray(VICON_C,VICON_R);

disp("select sample VICON trajectory snapshot");
Trajectories_Cut = csvFilesToAray(VICON_C,VICON_R);

disp("select AddBiomechanics motion data");
AddBio_Temp = motFilesToArray();
dataCells = {AddBio_Temp(:).data};
dataArr = cell2mat(dataCells(1));
for i = 2:length(dataCells)
    dataArr = [dataArr; cell2mat(dataCells(i))];
end
AddBio_Data.data = dataArr;
AddBio_Data.textdata = AddBio_Temp(1).textdata;
AddBio_Data.colheaders = AddBio_Temp(1).colheaders;


disp("select Dspace data")
[filename, filepath] = uigetfile(".mat");
DSpace_Data = load("" + filepath + filename);

disp("Select IMU data")
IMU_Data = csvFilesToAray(IMU_C, IMU_R); 


end

