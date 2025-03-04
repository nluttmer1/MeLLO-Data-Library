%% Process the full marker data
clear all
close all 
clc 

Full_marker = readtable('One Handed Test - With Case04.csv', 'PreserveVariableNames', true);
Full_marker_changed2 = zeros(size(Full_marker));
for i = 2:height(Full_marker)
    for j = 1:width(Full_marker)
        if iscellstr(Full_marker{i,j})
            Full_marker_changed2(i,j) = str2double(Full_marker{i,j});
        end
    end
end
save(' Reg Gait.mat', 'Full_marker_changed2');

%% 
clear all 
close all
clc

Full_marker_processed = load("Reg Gait.mat");
Full_marker_fix = fillmissing(Full_marker_processed.Full_marker_changed2(2:end,:),'constant',0);
Time = (1:1:length(Full_marker_fix))'/100;

Marker1 = Full_marker_fix;

Marker1_RToe_x = Marker1(:,198)/1000; % WheelBarrow 231 - Gait 198
Marker1_RToe_y = Marker1(:,199)/1000; % WheelBarrow 232 - Gait 199
Marker1_RToe_z = Marker1(:,200)/1000; % WheelBarrow 233 - Gait 200

Marker1_RHeel_x = Marker1(:,195)/1000; % WheelBarrow 228 - Gait 195
Marker1_RHeel_y = Marker1(:,196)/1000; % WheelBarrow 229 - Gait 196
Marker1_RHeel_z = Marker1(:,197)/1000; % WheelBarrow 230 - Gait 197

Marker1_RAnkle_x = Marker1(:,189)/1000; % WheelBarrow 222 - Gait 189
Marker1_RAnkle_y = Marker1(:,190)/1000; % WheelBarrow 223 - Gait 190
Marker1_RAnkle_z = Marker1(:,191)/1000; % WheelBarrow 224 - Gait 191

Marker1_LToe_x = Marker1(:,183)/1000; % WheelBarrow 183 - Gait 183 %%%%% CHECK!!!!
Marker1_LToe_y = Marker1(:,184)/1000; % WheelBarrow 184 - Gait 184
Marker1_LToe_z = Marker1(:,185)/1000; % WheelBarrow 185 - Gait 185

Marker1_LHeel_x = Marker1(:,147)/1000; % WheelBarrow 180 - Gait 147
Marker1_LHeel_y = Marker1(:,148)/1000; % WheelBarrow 181 - Gait 148
Marker1_LHeel_z = Marker1(:,149)/1000; % WheelBarrow 182 - Gait 149

Marker1_LAnkle_x = Marker1(:,141)/1000; % WheelBarrow 174 - Gait 141
Marker1_LAnkle_y = Marker1(:,142)/1000; % WheelBarrow 175 - Gait 142
Marker1_LAnkle_z = Marker1(:,143)/1000; % WheelBarrow 176 - Gait 143

RHeel_Start_Index1 = 534;
RHeel_End_Index1 = 657;
RHeel_Start_Index2 = 1420;
RHeel_End_Index2 = 1542;
RHeel_Start_Index3 = 2217;
RHeel_End_Index3 = 2242;
RHeel_Start_Index4 = 3013;
RHeel_End_Index4 = 3135;
RHeel_Start_Index5 = 3680;
RHeel_End_Index5 = 3800;
RHeel_Start_Index6 = 4542;
RHeel_End_Index6 = 4665;
RHeel_Start_Index7 = 5186;
RHeel_End_Index7 = 5307;
RHeel_Start_Index8 = 6060;
RHeel_End_Index8 = 6179;
RHeel_Start_Index9 = 6701;
RHeel_End_Index9 = 6823;
RHeel_Start_Index10 = 7437;
RHeel_End_Index10 = 7557;
RHeel_Start_Index = [RHeel_Start_Index1 RHeel_Start_Index2 RHeel_Start_Index3 RHeel_Start_Index4...
    RHeel_Start_Index5 RHeel_Start_Index6 RHeel_Start_Index7 RHeel_Start_Index8 RHeel_Start_Index9...
    RHeel_Start_Index10];
RHeel_End_Index = [RHeel_End_Index1 RHeel_End_Index2 RHeel_End_Index3 RHeel_End_Index4...
    RHeel_End_Index5 RHeel_End_Index6 RHeel_End_Index7 RHeel_End_Index8 RHeel_End_Index9...
    RHeel_End_Index10];

LHeel_Start_Index1 = 596;
LHeel_End_Index1 = 715;
LHeel_Start_Index2 = 1358;
LHeel_End_Index2 = 1483;
LHeel_Start_Index3 = 2180;
LHeel_End_Index3 = 2304;
LHeel_Start_Index4 = 2952;
LHeel_End_Index4 = 3075;
LHeel_Start_Index5 = 3741;
LHeel_End_Index5 = 3858;
LHeel_Start_Index6 = 4482;
LHeel_End_Index6 = 4603;
LHeel_Start_Index7 = 5246;
LHeel_End_Index7 = 5367;
LHeel_Start_Index8 = 6003;
LHeel_End_Index8 = 6120;
LHeel_Start_Index9 = 6761;
LHeel_End_Index9 = 6881;
LHeel_Start_Index10 = 7498;
LHeel_End_Index10 = 7618;
LHeel_Start_Index = [LHeel_Start_Index1 LHeel_Start_Index2 LHeel_Start_Index3 LHeel_Start_Index4...
    LHeel_Start_Index5 LHeel_Start_Index6 LHeel_Start_Index7 LHeel_Start_Index8 LHeel_Start_Index9...
    LHeel_Start_Index10];
LHeel_End_Index = [LHeel_End_Index1 LHeel_End_Index2 LHeel_End_Index3 LHeel_End_Index4...
    LHeel_End_Index5 LHeel_End_Index6 LHeel_End_Index7 LHeel_End_Index8 LHeel_End_Index9...
    LHeel_End_Index10];

figure
ax1 = subplot(4,1,1);
hold on
plot(Time*100,Marker1_RAnkle_z,'r',Time*100,Marker1_RHeel_z,'b',Time*100,Marker1_RToe_z,'k')
xline(RHeel_Start_Index,'k-')
xline(RHeel_End_Index,'r-')
legend('Ankle','Heel','Toe')
hold off

ax2 = subplot(4,1,2);
hold on
plot(Time*100,Marker1_RHeel_x,'r',Time*100,Marker1_RHeel_y,'b',Time*100,Marker1_RHeel_z,'k')
xline(RHeel_Start_Index,'k-')
xline(RHeel_End_Index,'r-')
legend('Ankle','Heel','Toe')
hold off

ax3 = subplot(4,1,3);
hold on
plot(Time*100,Marker1_LAnkle_z,'r',Time*100,Marker1_LHeel_z,'b',Time*100,Marker1_LToe_z,'k')
xline(LHeel_Start_Index,'k-')
xline(LHeel_End_Index,'r-')
legend('Ankle','Heel','Toe')
hold off

ax4 = subplot(4,1,4);
hold on
plot(Time*100,Marker1_LHeel_x,'r',Time*100,Marker1_LHeel_y,'b',Time*100,Marker1_LHeel_z,'k')
xline(LHeel_Start_Index,'k-')
xline(LHeel_End_Index,'r-')
legend('Ankle','Heel','Toe')
hold off

linkaxes([ax1,ax2,ax3,ax4],'x')

Step_Length = zeros([1,length(RHeel_Start_Index)]);
Stride_Length = zeros([1,length(RHeel_Start_Index)]);
Step_Time = zeros([1,length(RHeel_Start_Index)]);
Stride_Time = zeros([1,length(RHeel_Start_Index)]);
Walking_Speed = zeros([1,length(RHeel_Start_Index)]);
Cadence = zeros([1,length(RHeel_Start_Index)]);
Stride_Freq = zeros([1,length(RHeel_Start_Index)]);
for i = 1:1:length(RHeel_Start_Index)
    Step_Length(1,i) = abs(Marker1_LHeel_y(LHeel_Start_Index(1,i)) - Marker1_RHeel_y(RHeel_Start_Index(1,i)));%m * 3.28 % ft
    Stride_Length(1,i) = abs(Marker1_LHeel_y(LHeel_Start_Index(1,i)) - Marker1_LHeel_y(LHeel_End_Index(1,i)));%m * 3.28 % ft
    Step_Time(1,i) = abs(Time(RHeel_Start_Index(1,i)) - Time(LHeel_Start_Index(1,i))); % sec
    Stride_Time(1,i) = abs(Time(LHeel_Start_Index(1,i)) - Time(LHeel_End_Index(1,i)));% sec
    Walking_Speed(1,i) = Stride_Length(1,i) / Stride_Time(1,i); %m/s
    Cadence(1,i) = 1/Stride_Time(1,i)*60; %spm
    Stride_Freq(1,i) = Cadence(1,i)/2; %Hz
end

Step_Length_Mean = mean(Step_Length);
Step_Length_STD = std(Step_Length);
Step_Length_SE = Step_Length_STD/ sqrt(length(Step_Length));
Step_Length_df = length(Step_Length) - 1;
Step_Length_t_value = tinv(0.975, Step_Length_df);
Step_Length_Margin_of_Error = Step_Length_t_value * Step_Length_SE;
Step_Length_CI_Lower = Step_Length_Mean - Step_Length_Margin_of_Error;
Step_Length_CI_Upper = Step_Length_Mean + Step_Length_Margin_of_Error;

Stride_Length_Mean = mean(Stride_Length);
Stride_Length_STD = std(Stride_Length);
Stride_Length_SE = Stride_Length_STD/ sqrt(length(Stride_Length));
Stride_Length_df = length(Stride_Length) - 1;
Stride_Length_t_value = tinv(0.975, Stride_Length_df);
Stride_Length_Margin_of_Error = Stride_Length_t_value * Stride_Length_SE;
Stride_Length_CI_Lower = Stride_Length_Mean - Stride_Length_Margin_of_Error;
Stride_Length_CI_Upper = Stride_Length_Mean + Stride_Length_Margin_of_Error;

Step_Time_Mean = mean(Step_Time);
Step_Time_STD = std(Step_Time);
Step_Time_SE = Step_Time_STD/ sqrt(length(Step_Time));
Step_Time_df = length(Step_Time) - 1;
Step_Time_t_value = tinv(0.975, Step_Time_df);
Step_Time_Margin_of_Error = Step_Time_t_value * Step_Time_SE;
Step_Time_CI_Lower = Step_Time_Mean - Step_Time_Margin_of_Error;
Step_Time_CI_Upper = Step_Time_Mean + Step_Time_Margin_of_Error;

Stride_Time_Mean = mean(Stride_Time);
Stride_Time_STD = std(Stride_Time);
Stride_Time_SE = Stride_Time_STD/ sqrt(length(Stride_Time));
Stride_Time_df = length(Stride_Time) - 1;
Stride_Time_t_value = tinv(0.975, Stride_Time_df);
Stride_Time_Margin_of_Error = Stride_Time_t_value * Stride_Time_SE;
Stride_Time_CI_Lower = Stride_Time_Mean - Stride_Time_Margin_of_Error;
Stride_Time_CI_Upper = Stride_Time_Mean + Stride_Time_Margin_of_Error;

Walking_Speed_Mean = mean(Walking_Speed);
Walking_Speed_STD = std(Walking_Speed);
Walking_Speed_SE = Walking_Speed_STD/ sqrt(length(Walking_Speed));
Walking_Speed_df = length(Walking_Speed) - 1;
Walking_Speed_t_value = tinv(0.975, Walking_Speed_df);
Walking_Speed_Margin_of_Error = Walking_Speed_t_value * Walking_Speed_SE;
Walking_Speed_CI_Lower = Walking_Speed_Mean - Walking_Speed_Margin_of_Error;
Walking_Speed_CI_Upper = Walking_Speed_Mean + Walking_Speed_Margin_of_Error;

Cadence_Mean = mean(Cadence);
Cadence_STD = std(Cadence);
Cadence_SE = Cadence_STD/ sqrt(length(Cadence));
Cadence_df = length(Cadence) - 1;
Cadence_t_value = tinv(0.975, Cadence_df);
Cadence_Margin_of_Error = Cadence_t_value * Cadence_SE;
Cadence_CI_Lower = Cadence_Mean - Cadence_Margin_of_Error;
Cadence_CI_Upper = Cadence_Mean + Cadence_Margin_of_Error;

Stride_Freq_Mean = mean(Stride_Freq);
Stride_Freq_STD = std(Stride_Freq);
Stride_Freq_SE = Stride_Freq_STD/ sqrt(length(Stride_Freq));
Stride_Freq_df = length(Stride_Freq) - 1;
Stride_Freq_t_value = tinv(0.975, Stride_Freq_df);
Stride_Freq_Margin_of_Error = Stride_Freq_t_value * Stride_Freq_SE;
Stride_Freq_CI_Lower = Stride_Freq_Mean - Stride_Freq_Margin_of_Error;
Stride_Freq_CI_Upper = Stride_Freq_Mean + Stride_Freq_Margin_of_Error;


% 
% Para_labels = {'Step Length','Stride Length','Step Time','Stride Time','Walking Speed','Cadence','Stride Freq'};
% Para_Export = [Step_Length Stride_Length Step_Time Stride_Time Walking_Speed Cadence Stride_Freq];
% Para_output = [Para_labels; num2cell(Para_Export)];
% writecell(Para_output,'Parameters1.csv')
% 
% Index_labels = {'Left Heel Start','Left Heel End','Right Heel Start','Right Heel End'};
% Index_Export = [LHeel_Start_Index LHeel_End_Index RHeel_Start_Index RHeel_End_Index];
% Index_output = [Index_labels; num2cell(Index_Export)];
% writecell(Index_output,'Index2.csv')
% 
%% Joint Data

AddBio_Temp = importdata('Subject 1 Wheelbarrow Walk 1_segment_0_ik.mot');