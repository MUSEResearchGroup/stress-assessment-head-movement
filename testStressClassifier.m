%% Classification Accuracy: Test
    clear all
    clc

%This scripts contains the training and validation of the Stress Classifier. 
%Please refer to the paper: Ferrarotti et al., "Stress Assessment for Augmented Reality Applications Based on Head Movement Features",  IEEE TVCG, 2024.
%for further details on the proposed architecture. 
%% Load trained models
% trained_models_unw.mat contains all the SVMs trained on the head
% movement data of users performing the Stroop Color Word Test (SCWT).
    
load trained_models_unw.mat
%% Data Extraction: Stroop Color Word Test
%Extraction of recorded head movement data from the 20 participants who
%performed the Stroop Color Word Test.

%stroopNames.txt is a file containing all the directory paths of the txt file
%containing the head movements recordings of the 20 users performing the
%SCWT, representing the test sample
fid = fopen('.\file_paths\stroopNames.txt','r');
tline = fgetl(fid);
nomeFile=[];

while ischar(tline)
    nome=string(tline);
    nomeFile= [nomeFile; nome];
    tline = fgetl(fid);
end
fclose(fid);

coor_disp_phase1 = [];
coor_disp_phase2 = [];
coor_disp_phase3 = [];

for i=1:size(nomeFile,1)
    [disp_phase1, disp_phase2, disp_phase3] = getFeatures_fixed(nomeFile(i));

    if i == 1
       coor_disp_phase1 = disp_phase1; 
    elseif size(coor_disp_phase1,1) < size(disp_phase1,1)
        coor_disp_phase1 = [coor_disp_phase1; coor_disp_phase1(end,:,:).*ones(size(disp_phase1,1)-size(coor_disp_phase1,1), size(coor_disp_phase1,2),size(coor_disp_phase1,3))];
        coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    elseif size(coor_disp_phase1,1) > size(disp_phase1,1)
        disp_phase1 = [disp_phase1; disp_phase1(end,:).*ones(size(coor_disp_phase1,1)-size(disp_phase1,1),size(coor_disp_phase1,2))];
        coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    else
         coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    end

    if i == 1
       coor_disp_phase2 = disp_phase2;  
    elseif size(coor_disp_phase2,1) < size(disp_phase2,1)
        coor_disp_phase2 = [coor_disp_phase2; coor_disp_phase2(end,:,:).*ones(size(disp_phase2,1)-size(coor_disp_phase2,1), size(coor_disp_phase2,2),size(coor_disp_phase2,3))];
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    elseif size(coor_disp_phase2,1) > size(disp_phase2,1)
        disp_phase2 = [disp_phase2; disp_phase2(end,:).*ones(size(coor_disp_phase2,1)-size(disp_phase2,1),size(coor_disp_phase2,2))];
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    else
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    end


    if i == 1
       coor_disp_phase3 = disp_phase3; 
    elseif size(coor_disp_phase3,1) < size(disp_phase3,1)
        coor_disp_phase3 = [coor_disp_phase3; coor_disp_phase3(end,:,:).*ones(size(disp_phase3,1)-size(coor_disp_phase3,1), size(coor_disp_phase3,2),size(coor_disp_phase3,3))];
        coor_disp_phase3 =  cat(3, coor_disp_phase3, disp_phase3);
    elseif size(coor_disp_phase3,1) > size(disp_phase3,1)
        disp_phase3 = [disp_phase3; disp_phase3(end,:).*ones(size(coor_disp_phase3,1)-size(disp_phase3,1),size(coor_disp_phase3,2))];
        coor_disp_phase3 =  cat(3, coor_disp_phase3, disp_phase3);
    else
        coor_disp_phase3 =  cat(3, coor_disp_phase3, disp_phase3);
    end
end

speed_phase1 = 10*diff(coor_disp_phase1);
speed_phase2 = 10*diff(coor_disp_phase2);
speed_phase3 = 10*diff(coor_disp_phase3);

tot_disp_phase1=[];
tot_disp_phase2=[];
tot_disp_phase3=[];

for j=1:size(coor_disp_phase1,3)
    for i=1:length(coor_disp_phase1)
        tot_disp_phase1(i,j) = norm(coor_disp_phase1(i,:,j));
    end
end

for j=1:size(coor_disp_phase2,3)
    for i=1:length(coor_disp_phase2)
        tot_disp_phase2(i,j) = norm(coor_disp_phase2(i,:,j));
    end
end

for j=1:size(coor_disp_phase3,3)
    for i=1:length(coor_disp_phase3)
        tot_disp_phase3(i,j) = norm(coor_disp_phase3(i,:,j));
    end
end

speed_phase1 = 10*diff(coor_disp_phase1);
magn_speed_phase1 = zeros(length(speed_phase1), size(speed_phase1,3));

for j = 1:size(speed_phase1,3)
    for i=1:length(speed_phase1)
        magn_speed_phase1(i,j) = norm(speed_phase1(i,:,j));
    end
end

speed_phase2 = 10*diff(coor_disp_phase2);
magn_speed_phase2 = zeros(length(speed_phase2), 1);

for j = 1:size(speed_phase2,3)
    for i=1:size(speed_phase2,1)
        magn_speed_phase2(i,j) = norm(speed_phase2(i,:,j));
    end
end

speed_phase3 = 10*diff(coor_disp_phase3);
magn_speed_phase3 = zeros(length(speed_phase3), 1);

for j = 1:size(speed_phase3,3)
    for i=1:size(speed_phase3,1)
        magn_speed_phase3(i,j) = norm(speed_phase3(i,:,j));
    end
end

angle_yz_1 = zeros(size(coor_disp_phase1,1),size(coor_disp_phase1,3));
angle_yz_2 = zeros(size(coor_disp_phase2,1),size(coor_disp_phase2,3));
angle_yz_3 = zeros(size(coor_disp_phase3,1),size(coor_disp_phase3,3));

angle_xz_1 = zeros(size(coor_disp_phase1,1),size(coor_disp_phase1,3));
angle_xz_2 = zeros(size(coor_disp_phase2,1),size(coor_disp_phase2,3));
angle_xz_3 = zeros(size(coor_disp_phase3,1),size(coor_disp_phase3,3));

angle_xy_1 = zeros(size(coor_disp_phase1,1),size(coor_disp_phase1,3));
angle_xy_2 = zeros(size(coor_disp_phase2,1),size(coor_disp_phase2,3));
angle_xy_3 = zeros(size(coor_disp_phase3,1),size(coor_disp_phase3,3));

for i=1:size(coor_disp_phase1,3)
    angle_yz_1(:,i) = atan(coor_disp_phase1(:,2,i)./coor_disp_phase1(:,3,i));
    angle_yz_2(:,i) = atan(coor_disp_phase2(:,2,i)./coor_disp_phase2(:,3,i));
    angle_yz_3(:,i) = atan(coor_disp_phase3(:,2,i)./coor_disp_phase3(:,3,i));
    
    angle_xz_1(:,i) = atan(coor_disp_phase1(:,1,i)./coor_disp_phase1(:,3,i));
    angle_xz_2(:,i) = atan(coor_disp_phase2(:,1,i)./coor_disp_phase2(:,3,i));
    angle_xz_3(:,i) = atan(coor_disp_phase3(:,1,i)./coor_disp_phase3(:,3,i));
    
    angle_xy_1(:,i) = atan(coor_disp_phase1(:,2,i)./coor_disp_phase1(:,1,i));
    angle_xy_2(:,i) = atan(coor_disp_phase2(:,2,i)./coor_disp_phase2(:,1,i));
    angle_xy_3(:,i) = atan(coor_disp_phase3(:,2,i)./coor_disp_phase3(:,1,i));
    
end
for i = 1:size(angle_xy_3,2)
    angle_yz_1(isnan(angle_yz_1(:,i)),i) = 0;
    angle_xz_1(isnan(angle_xz_1(:,i)),i) = 0;
    angle_xy_1(isnan(angle_xy_1(:,i)),i) = 0;

    angle_yz_2(isnan(angle_yz_2(:,i)),i) = 0;
    angle_xz_2(isnan(angle_xz_2(:,i)),i) = 0;
    angle_xy_2(isnan(angle_xy_2(:,i)),i) = 0;

    angle_yz_3(isnan(angle_yz_3(:,i)),i) = 0;
    angle_xz_3(isnan(angle_xz_3(:,i)),i) = 0;
    angle_xy_3(isnan(angle_xy_3(:,i)),i) = 0;
end
%% STFT ANALYSIS: Stroop Color Word Test
freq_1 = squeeze(sum(unwrap(angle(stft(tot_disp_phase1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2 = squeeze(sum(unwrap(angle(stft(tot_disp_phase2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_3 = squeeze(sum(unwrap(angle(stft(tot_disp_phase3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

freq_1_xy = squeeze(sum(unwrap(angle(stft(angle_xy_1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xy = squeeze(sum(unwrap(angle(stft(angle_xy_2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_3_xy = squeeze(sum(unwrap(angle(stft(angle_xy_3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

freq_1_xz = squeeze(sum(unwrap(angle(stft(angle_xz_1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xz = squeeze(sum(unwrap(angle(stft(angle_xz_2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_3_xz = squeeze(sum(unwrap(angle(stft(angle_xz_3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

freq_1_yz = squeeze(sum(unwrap(angle(stft(angle_yz_1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_yz = squeeze(sum(unwrap(angle(stft(angle_yz_2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_3_yz = squeeze(sum(unwrap(angle(stft(angle_yz_3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

freq_1_speed = squeeze(sum(unwrap(angle(stft(magn_speed_phase1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_speed = squeeze(sum(unwrap(angle(stft(magn_speed_phase2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_3_speed = squeeze(sum(unwrap(angle(stft(magn_speed_phase3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

%% Data Extraction: Mental Arithmetic Test
%Extraction of recorded head movement data from the 20 participants who
%performed the Mental Arithmetic (MA) Test.


%maNames.txt is a file containing all the directory paths of the txt file
%containing the head movements recordings of the 20 users performing the
%MA, representing the test sample

fid = fopen('.\file_paths\maNames.txt','r');
tline = fgetl(fid);
nomeFile_MA=[];

while ischar(tline)
    nome=string(tline);
    nomeFile_MA= [nomeFile_MA; nome];
    tline = fgetl(fid);
end
fclose(fid);

coor_disp_phase1 = [];
coor_disp_phase2 = [];
coor_disp_phase3 = [];

for i=1:size(nomeFile_MA,1)
    [disp_phase1, disp_phase2] = getFeatures_MA(nomeFile_MA(i));

    if i == 1
       coor_disp_phase1 = disp_phase1; 
    elseif size(coor_disp_phase1,1) < size(disp_phase1,1)
        coor_disp_phase1 = [coor_disp_phase1; coor_disp_phase1(end,:,:).*ones(size(disp_phase1,1)-size(coor_disp_phase1,1), size(coor_disp_phase1,2),size(coor_disp_phase1,3))];
        coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    elseif size(coor_disp_phase1,1) > size(disp_phase1,1)
        disp_phase1 = [disp_phase1; disp_phase1(end,:).*ones(size(coor_disp_phase1,1)-size(disp_phase1,1),size(coor_disp_phase1,2))];
        coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    else
         coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    end

    if i == 1
       coor_disp_phase2 = disp_phase2;  
    elseif size(coor_disp_phase2,1) < size(disp_phase2,1)
        coor_disp_phase2 = [coor_disp_phase2; coor_disp_phase2(end,:,:).*ones(size(disp_phase2,1)-size(coor_disp_phase2,1), size(coor_disp_phase2,2),size(coor_disp_phase2,3))];
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    elseif size(coor_disp_phase2,1) > size(disp_phase2,1)
        disp_phase2 = [disp_phase2; disp_phase2(end,:).*ones(size(coor_disp_phase2,1)-size(disp_phase2,1),size(coor_disp_phase2,2))];
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    else
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    end
end

speed_phase1_MA = 10*diff(coor_disp_phase1);
speed_phase2_MA = 10*diff(coor_disp_phase2);

tot_disp_phase1_MA=[];
tot_disp_phase2_MA=[];

for j=1:size(coor_disp_phase1,3)
    for i=1:length(coor_disp_phase1)
        tot_disp_phase1_MA(i,j) = norm(coor_disp_phase1(i,:,j));
    end
end

for j=1:size(coor_disp_phase2,3)
    for i=1:length(coor_disp_phase2)
        tot_disp_phase2_MA(i,j) = norm(coor_disp_phase2(i,:,j));
    end
end

magn_speed_phase1_MA = zeros(length(speed_phase1_MA), size(speed_phase1_MA,3));

for j = 1:size(speed_phase1_MA,3)
    for i=1:length(speed_phase1_MA)
        magn_speed_phase1_MA(i,j) = norm(speed_phase1_MA(i,:,j));
    end
end

magn_speed_phase2_MA = zeros(length(speed_phase2_MA), 1);

for j = 1:size(speed_phase2_MA,3)
    for i=1:size(speed_phase2_MA,1)
        magn_speed_phase2_MA(i,j) = norm(speed_phase2_MA(i,:,j));
    end
end

angle_yz_1_MA = zeros(size(coor_disp_phase1,1),size(coor_disp_phase1,3));
angle_yz_2_MA = zeros(size(coor_disp_phase2,1),size(coor_disp_phase2,3));

angle_xz_1_MA = zeros(size(coor_disp_phase1,1),size(coor_disp_phase1,3));
angle_xz_2_MA = zeros(size(coor_disp_phase2,1),size(coor_disp_phase2,3));

angle_xy_1_MA = zeros(size(coor_disp_phase1,1),size(coor_disp_phase1,3));
angle_xy_2_MA = zeros(size(coor_disp_phase2,1),size(coor_disp_phase2,3));

for i=1:size(coor_disp_phase1,3)
    angle_yz_1_MA(:,i) = atan(coor_disp_phase1(:,2,i)./coor_disp_phase1(:,3,i));
    angle_yz_2_MA(:,i) = atan(coor_disp_phase2(:,2,i)./coor_disp_phase2(:,3,i));
    
    angle_xz_1_MA(:,i) = atan(coor_disp_phase1(:,1,i)./coor_disp_phase1(:,3,i));
    angle_xz_2_MA(:,i) = atan(coor_disp_phase2(:,1,i)./coor_disp_phase2(:,3,i));
    
    angle_xy_1_MA(:,i) = atan(coor_disp_phase1(:,2,i)./coor_disp_phase1(:,1,i));
    angle_xy_2_MA(:,i) = atan(coor_disp_phase2(:,2,i)./coor_disp_phase2(:,1,i));
    
end
for i = 1:size(angle_xy_1_MA,2)
    angle_yz_1_MA(isnan(angle_yz_1_MA(:,i)),i) = 0;
    angle_xz_1_MA(isnan(angle_xz_1_MA(:,i)),i) = 0;
    angle_xy_1_MA(isnan(angle_xy_1_MA(:,i)),i) = 0;

    angle_yz_2_MA(isnan(angle_yz_2(:,i)),i) = 0;
    angle_xz_2_MA(isnan(angle_xz_2(:,i)),i) = 0;
    angle_xy_2_MA(isnan(angle_xy_2(:,i)),i) = 0;
end

%% STFT ANALYSIS: Mental Arithmetic Test

%The frequency analysis is performed for each subsequence of data. Then,
%the mean of the results is evaluated.
freq_1_MA_1 = squeeze(sum(unwrap(angle(stft(tot_disp_phase1_MA(1:450,:), 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_1_MA_2 = squeeze(sum(unwrap(angle(stft(tot_disp_phase1_MA(451:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_1_MA = (freq_1_MA_1+freq_1_MA_2)/2;

freq_2_MA_1 = squeeze(sum(unwrap(angle(stft(tot_disp_phase2_MA(1:650,:), 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2_MA_2 = squeeze(sum(unwrap(angle(stft(tot_disp_phase2_MA(651:1301,:), 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2_MA_3 = squeeze(sum(unwrap(angle(stft(tot_disp_phase2_MA(1302:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2_MA = (freq_2_MA_1+freq_2_MA_2+freq_2_MA_3)/3;

freq_1_xy_MA_1 = squeeze(sum(unwrap(angle(stft(angle_xy_1_MA(1:450,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_xy_MA_2 = squeeze(sum(unwrap(angle(stft(angle_xy_1_MA(451:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_xy_MA = (freq_1_xy_MA_1+freq_1_xy_MA_2)/2;

freq_2_xy_MA_1 = squeeze(sum(unwrap(angle(stft(angle_xy_2_MA(1:650,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xy_MA_2 = squeeze(sum(unwrap(angle(stft(angle_xy_2_MA(651:1301,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xy_MA_3 = squeeze(sum(unwrap(angle(stft(angle_xy_2_MA(1302:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xy_MA = (freq_2_xy_MA_1+freq_2_xy_MA_2+freq_2_xy_MA_3)/3;

freq_1_xz_MA_1 = squeeze(sum(unwrap(angle(stft(angle_xz_1_MA(1:450,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_xz_MA_2 = squeeze(sum(unwrap(angle(stft(angle_xz_1_MA(451:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_xz_MA = (freq_1_xz_MA_1+freq_1_xz_MA_2)/2;

freq_2_xz_MA_1 = squeeze(sum(unwrap(angle(stft(angle_xz_2_MA(1:650,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xz_MA_2 = squeeze(sum(unwrap(angle(stft(angle_xz_2_MA(651:1301,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xz_MA_3 = squeeze(sum(unwrap(angle(stft(angle_xz_2_MA(1302:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xz_MA = (freq_2_xz_MA_1+freq_2_xz_MA_2+freq_2_xz_MA_3)/3;

freq_1_yz_MA_1 = squeeze(sum(unwrap(angle(stft(angle_yz_1_MA(1:450,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_yz_MA_2 = squeeze(sum(unwrap(angle(stft(angle_yz_1_MA(451:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_yz_MA = (freq_1_xy_MA_1+freq_1_xy_MA_2)/2;

freq_2_yz_MA_1 = squeeze(sum(unwrap(angle(stft(angle_yz_2_MA(1:650,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_yz_MA_2 = squeeze(sum(unwrap(angle(stft(angle_yz_2_MA(651:1301,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_yz_MA_3 = squeeze(sum(unwrap(angle(stft(angle_yz_2_MA(1302:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_yz_MA = (freq_2_yz_MA_1+freq_2_yz_MA_2+freq_2_yz_MA_3)/3;

freq_1_speed_MA_1 = squeeze(sum(unwrap(angle(stft(magn_speed_phase1_MA(1:450,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_speed_MA_2 = squeeze(sum(unwrap(angle(stft(magn_speed_phase1_MA(451:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_1_speed_MA = (freq_1_speed_MA_1+freq_1_speed_MA_2)/2;

freq_2_speed_MA_1 = squeeze(sum(unwrap(angle(stft(magn_speed_phase2_MA(1:650,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_speed_MA_2 = squeeze(sum(unwrap(angle(stft(magn_speed_phase2_MA(651:1301,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_speed_MA_3 = squeeze(sum(unwrap(angle(stft(magn_speed_phase2_MA(1302:end,:), 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_speed_MA = (freq_2_speed_MA_1+freq_2_speed_MA_2+freq_2_speed_MA_3)/3;
%% ACCURACY: Stroop Color Word Test
%The accuracy of the trained classifiers is evaluated for each of the
%computed features.
labels_test = [zeros(20*2,1); ones(20,1)];

test = [freq_1; freq_2; freq_3];
pred_freq_disp = predict(Mdl_freq_disp,test);
accuracy_disp_S = sum(labels_test==pred_freq_disp)/size(labels_test,1)

test = [freq_1_xy; freq_2_xy; freq_3_xy];
pred_freq_xy = predict(Mdl_freq_xy,test);
accuracy_xy_S = sum(labels_test==pred_freq_xy)/size(labels_test,1)

test = [freq_1_xz; freq_2_xz; freq_3_xz];
pred_freq_xz = predict(Mdl_freq_xz,test);
accuracy_xz_S = sum(labels_test==pred_freq_xz)/size(labels_test,1)

test = [freq_1_yz; freq_2_yz; freq_3_yz];
pred_freq_yz = predict(Mdl_freq_yz,test);
accuracy_yz_S = sum(labels_test==pred_freq_yz)/size(labels_test,1)

pred_freq_angles_MA = nan(size(pred_freq_yz,1), size(pred_freq_yz,2));

for i=1:size(pred_freq_yz,1)
    predictions = [pred_freq_xy(i) pred_freq_xz(i) pred_freq_yz(i)];
    if length(find(predictions==1)) >= 2
       pred_freq_angles_MA(i) = 1;
    else
       pred_freq_angles_MA(i) = 0;
    end
end

accuracy_angle_S = sum(labels_test==pred_freq_angles_MA)/size(labels_test,1)

test = [freq_1_speed; freq_2_speed; freq_3_speed];
pred_freq_speed = predict(Mdl_freq_speed,test);
accuracy_speed_S = sum(labels_test==pred_freq_speed)/size(labels_test,1)
%% ACCURACY: Mental Arithmetic Test
%The accuracy of the trained classifiers is evaluated for each of the
%computed features.
labels_test_MA = [zeros(20,1); ones(20,1)];

test = [freq_1_MA; freq_2_MA];
pred_freq_disp_MA = predict(Mdl_freq_disp,test);
accuracy_disp_MA = sum(labels_test_MA==pred_freq_disp_MA)/size(labels_test_MA,1)

test = [freq_1_xy_MA; freq_2_xy_MA];
pred_freq_xy_MA = predict(Mdl_freq_xy,test);
accuracy_xy_MA = sum(labels_test_MA==pred_freq_xy_MA)/size(labels_test_MA,1)

test = [freq_1_xz_MA; freq_2_xz_MA];
pred_freq_xz_MA = predict(Mdl_freq_xz,test);
accuracy_xz_MA = sum(labels_test_MA==pred_freq_xz_MA)/size(labels_test_MA,1)

test = [freq_1_yz_MA; freq_2_yz_MA];
pred_freq_yz_MA = predict(Mdl_freq_yz,test);
accuracy_yz_MA = sum(labels_test_MA==pred_freq_yz_MA)/size(labels_test_MA,1)

pred_freq_angles_MA = nan(size(pred_freq_yz_MA,1), size(pred_freq_yz_MA,2));

for i=1:size(pred_freq_yz_MA,1)
    predictions = [pred_freq_xy_MA(i) pred_freq_xz_MA(i) pred_freq_yz_MA(i)];
    if length(find(predictions==1)) >= 2
       pred_freq_angles_MA(i) = 1;
    else
       pred_freq_angles_MA(i) = 0;
    end
end

accuracy_angle__MA = sum(labels_test_MA==pred_freq_angles_MA)/size(labels_test_MA,1)

test = [freq_1_speed_MA; freq_2_speed_MA];
pred_freq_speed_MA = predict(Mdl_freq_speed,test);
accuracy_speed_MA = sum(labels_test_MA==pred_freq_speed_MA)/size(labels_test_MA,1)
%% Majority voting angles
pred_freq_angles_MA = nan(size(pred_freq_yz_MA,1), size(pred_freq_yz_MA,2));

for i=1:size(pred_freq_yz_MA,1)
    predictions = [pred_freq_xy_MA(i) pred_freq_xz_MA(i) pred_freq_yz_MA(i)];
    if length(find(predictions==1)) >= 2
       pred_freq_angles_MA(i) = 1;
    else
       pred_freq_angles_MA(i) = 0;
    end
end

accuracy_angle_mv_MA = sum(labels_test_MA==pred_freq_angles_MA)/size(labels_test_MA,1)

pred_freq_angles_S = nan(size(pred_freq_yz,1), size(pred_freq_yz,2));

for i=1:size(pred_freq_yz,1)
    predictions = [pred_freq_xy(i) pred_freq_xz(i) pred_freq_yz(i)];
    if length(find(predictions==1)) >= 2
       pred_freq_angles_S(i) = 1;
    else
       pred_freq_angles_S(i) = 0;
    end
end

accuracy_angle_mv_S = sum(labels_test==pred_freq_angles_S)/size(labels_test,1)
%% Combined classifier accuracy
%The accuracy of the proposed combined classifier is evaluated.
val_disp = 1;
val_angles = 0.98;
val_speed = 0.98;

w1 = (val_disp/val_angles)/((val_disp/val_angles)+(val_angles/val_angles)+(val_speed/val_angles));
w2 = (val_angles/val_angles)/((val_disp/val_angles)+(val_angles/val_angles)+(val_speed/val_angles));
w3 = (val_speed/val_angles)/((val_disp/val_angles)+(val_angles/val_angles)+(val_speed/val_angles));

c1_MA = ones(40,1);
c1_MA(pred_freq_disp_MA==0) = -1;
c2_MA = ones(40,1);
c2_MA(pred_freq_angles_MA==0) = -1;
c3_MA = ones(40,1);
c3_MA(pred_freq_speed_MA==0) = -1;

d_MA = w1*c1_MA + w2*c2_MA + w3*c3_MA;

predictions_MA = nan(40,1);
predictions_MA(d_MA>0) = 1;
predictions_MA(d_MA<0) = 0;

c1_S = ones(60,1);
c1_S(pred_freq_disp==0) = -1;
c2_S = ones(60,1);
c2_S(pred_freq_angles_S==0) = -1;
c3_S = ones(60,1);
c3_S(pred_freq_speed==0) = -1;

d_S = w1*c1_S + w2*c2_S + w3*c3_S;
predictions_S = nan(60,1);
predictions_S(d_S>0) = 1;
predictions_S(d_S<0) = 0;

accuracy_combined_MA = sum(labels_test_MA==predictions_MA)/size(labels_test_MA,1)
accuracy_combined_S = sum(labels_test==predictions_S)/size(labels_test,1)

