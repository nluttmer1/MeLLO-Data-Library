%%%%%%%%%%%%%%%%%%%
% Nathaniel Luttmer
% 1/17/25 
% Cyclogram
%%%%%%%%%%%%%%%%%%%

clc
clear all 
close all 

Reg_Gait_markers = csvread('Subject 1 Gait 5.csv',6,0);
Table_Reg = readtable('Subject 1 Gait 5.csv');
Reg_Gait_AddBio = importdata('Subject 1 Gait 5_segment_0_ik.mot');

Reg_Gait_Time_index = Reg_Gait_markers(:,1);
Reg_Gait_Toe_z = Reg_Gait_markers(:,200);% 184 - z
Reg_Gait_Heel_z = Reg_Gait_markers(:,197); % 181 - z
Reg_Gait_Ankle_z = Reg_Gait_markers(:,191); % 175 - z

Reg_AddBio_Time = Reg_Gait_AddBio.data(:,1);
Reg_AddBio_knee_Right = Reg_Gait_AddBio.data(:,11)*180/pi();
Reg_AddBio_knee_Left = Reg_Gait_AddBio.data(:,18)*180/pi();
Reg_AddBio_hip_right_flexion = Reg_Gait_AddBio.data(:,8)*180/pi();
Reg_AddBio_hip_right_adduction = Reg_Gait_AddBio.data(:,9)*180/pi();
Reg_AddBio_hip_right_rotation = Reg_Gait_AddBio.data(:,10)*180/pi();

Reg_AddBio_hip_left_flexion = Reg_Gait_AddBio.data(:,15)*180/pi();
Reg_AddBio_hip_left_adduction = Reg_Gait_AddBio.data(:,16)*180/pi();
Reg_AddBio_hip_left_rotation = Reg_Gait_AddBio.data(:,17)*180/pi();

Reg_AddBio_hip_resultant = zeros(size(Reg_AddBio_hip_right_flexion));
for i = 1:1:length(Reg_AddBio_hip_right_flexion)
    Reg_AddBio_hip_resultant(i) = sqrt(Reg_AddBio_hip_right_adduction(i)^2 + Reg_AddBio_hip_right_flexion(i)^2 +...
        Reg_AddBio_hip_right_rotation(i)^2);
end
% 
% figure
% hold on
% plot(Reg_Gait_Time_index,Reg_Gait_Toe_z,'k-')
% plot(Reg_Gait_Time_index,Reg_Gait_Heel_z,'r-')
% plot(Reg_Gait_Time_index,Reg_Gait_Ankle_z,'b-')
% legend('Toes','Ankle', 'Heel')

% Plot toe off / mid stance
Reg_zToesPos_L_NL_1 = Reg_Gait_Toe_z;
Reg_zAnklePos_L_NL_1 = Reg_Gait_Ankle_z;
Reg_zHeelPos_L_NL_1 = Reg_Gait_Heel_z;

Reg_Heel_strike = [8,125];

figure
hold on
plot(Reg_zToesPos_L_NL_1)
plot(Reg_zAnklePos_L_NL_1)
plot(Reg_zHeelPos_L_NL_1)
legend('Toes','Ankle', 'Heel')
% xline(toe_offs,'--r')
xline(Reg_Heel_strike,'b--')

Reg_Marker_time = (Reg_Gait_Time_index - Reg_Gait_Time_index(1));
Reg_Cut_time_start = Reg_Marker_time(Reg_Heel_strike(1));
Reg_Cut_time_end = Reg_Marker_time(Reg_Heel_strike(2));

% figure
% hold on
% plot(Reg_AddBio_Time,Reg_AddBio_knee_Right)
% xline(Reg_AddBio_Time(Reg_Cut_time_start),'--r')
% xline(Reg_AddBio_Time(Reg_Cut_time_end),'--r')
% xlabel('Time [s]')
% ylabel('Knee Angle [Degrees]')
% title('AddBio Knee Rotation')
% hold off

Reg_Cyclogram_RKnee = Reg_AddBio_knee_Right(Reg_Cut_time_start:Reg_Cut_time_end);
Reg_Cyclogram_LKnee = Reg_AddBio_knee_Left(Reg_Cut_time_start:Reg_Cut_time_end);
Reg_Cyclogram_hip_resultant = Reg_AddBio_hip_resultant(Reg_Cut_time_start:Reg_Cut_time_end);

% Whellbarrow
WB_Gait_markers = csvread('Subject 1 Wheelbarrow Walk 5.csv',6,0);
Table_WB = readtable('Subject 1 Wheelbarrow Walk 5.csv');
WB_Gait_AddBio = importdata('Subject 1 Wheelbarrow Walk 5_segment_0_ik.mot');

WB_Gait_Time_index = WB_Gait_markers(:,1);
WB_Gait_Toe_z = WB_Gait_markers(:,233);% 184 - z
WB_Gait_Heel_z = WB_Gait_markers(:,230); % 181 - z
WB_Gait_Ankle_z = WB_Gait_markers(:,224); % 175 - z

WB_AddBio_Time = WB_Gait_AddBio.data(:,1);
WB_AddBio_knee_Right = WB_Gait_AddBio.data(:,11)*180/pi();
WB_AddBio_knee_Left = WB_Gait_AddBio.data(:,18)*180/pi();
Wb_AddBio_hip_right_flexion = WB_Gait_AddBio.data(:,8)*180/pi();
Wb_AddBio_hip_right_adduction = WB_Gait_AddBio.data(:,9)*180/pi();
Wb_AddBio_hip_right_rotation = WB_Gait_AddBio.data(:,10)*180/pi();

WB_AddBio_hip_resultant = zeros(size(Wb_AddBio_hip_right_flexion));
for i = 1:1:length(Wb_AddBio_hip_right_flexion)
    WB_AddBio_hip_resultant(i) = sqrt(Wb_AddBio_hip_right_adduction(i)^2 + Wb_AddBio_hip_right_flexion(i)^2 +...
        Wb_AddBio_hip_right_rotation(i)^2);
end

% figure 
% hold on
% plot(Wb_AddBio_hip_right_flexion)
% plot(Wb_AddBio_hip_right_adduction)
% plot(Wb_AddBio_hip_right_rotation)
% legend('Flexion','Adduction','Rotation')
% hold off

% figure
% hold on
% plot(WB_Gait_Time_index,WB_Gait_Toe_z,'k-')
% plot(WB_Gait_Time_index,WB_Gait_Heel_z,'r-')
% plot(WB_Gait_Time_index,WB_Gait_Ankle_z,'b-')
% legend('Toes','Ankle', 'Heel')

% Plot toe off / mid stance
WB_zToesPos_L_NL_1 = WB_Gait_Toe_z;
WB_zAnklePos_L_NL_1 = WB_Gait_Ankle_z;
WB_zHeelPos_L_NL_1 = WB_Gait_Heel_z;

WB_Heel_strike = [26,141];

figure
hold on
plot(WB_zToesPos_L_NL_1)
plot(WB_zAnklePos_L_NL_1)
plot(WB_zHeelPos_L_NL_1)
legend('Toes','Ankle', 'Heel')
% xline(toe_offs,'--r')
xline(WB_Heel_strike,'b--')

WB_Marker_time = (WB_Gait_Time_index - WB_Gait_Time_index(1));
WB_Cut_time_start = WB_Marker_time(WB_Heel_strike(1));
WB_Cut_time_end = WB_Marker_time(WB_Heel_strike(2));

% figure
% hold on
% plot(WB_AddBio_Time,WB_AddBio_knee_Right)
% xline(WB_AddBio_Time(WB_Cut_time_start),'--r')
% xline(WB_AddBio_Time(WB_Cut_time_end),'--r')
% xlabel('Time [s]')
% ylabel('Knee Angle [Degrees]')
% title('AddBio Knee Rotation')
% hold off

WB_Cyclogram_RKnee = WB_AddBio_knee_Right(WB_Cut_time_start:WB_Cut_time_end);
WB_Cyclogram_LKnee = WB_AddBio_knee_Left(WB_Cut_time_start:WB_Cut_time_end);
WB_Cyclogram_hip_resultant = WB_AddBio_hip_resultant(WB_Cut_time_start:WB_Cut_time_end);

figure
hold on
plot(Reg_Cyclogram_RKnee,Reg_Cyclogram_LKnee,'r-o','LineWidth',1.5)
plot(WB_Cyclogram_RKnee,WB_Cyclogram_LKnee,'b-o','LineWidth',1.5)
plot(Reg_Cyclogram_RKnee(1),Reg_Cyclogram_LKnee(1),'go',Reg_Cyclogram_RKnee(end),Reg_Cyclogram_LKnee(end),'ko','LineWidth',3)
plot(WB_Cyclogram_RKnee(1),WB_Cyclogram_LKnee(1),'go',WB_Cyclogram_RKnee(end),WB_Cyclogram_LKnee(end),'ko','LineWidth',3)
xlabel('Right Knee Angle [deg]')
ylabel('Left Knee Angle [deg]')
title('Knee Cyclogram')
% legend('Regular Gait','Wheelbarrow Gait','Start Point','End Point')
grid on
grid minor
hold off

figure
hold on
plot(Reg_Cyclogram_RKnee,Reg_Cyclogram_hip_resultant,'r-o','LineWidth',1.5)
plot(WB_Cyclogram_RKnee, WB_Cyclogram_hip_resultant,'b-o','LineWidth',1.5)
plot(Reg_Cyclogram_RKnee(1),Reg_Cyclogram_hip_resultant(1),'go',Reg_Cyclogram_RKnee(end),Reg_Cyclogram_hip_resultant(end),'ko','LineWidth',3)
plot(WB_Cyclogram_RKnee(1),WB_Cyclogram_hip_resultant(1),'go',WB_Cyclogram_RKnee(end),WB_Cyclogram_hip_resultant(end),'ko','LineWidth',3)
xlabel('Right Knee Angle [deg]')
ylabel('Right Hip Angle [deg]')
title('Knee vs. Hip Cyclogram')
legend('Regular Gait','Wheelbarrow Gait','Start Point','End Point','Location','northeastoutside')
grid on
grid minor
hold off
