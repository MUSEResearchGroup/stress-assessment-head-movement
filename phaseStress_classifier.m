%% Stress classifier
%This script contains the training and validation of the Stress Classifier. 
%Please refer to the paper: Ferrarotti et al., "Stress Assessment for Augmented Reality Applications Based on Head Movement Features",  IEEE TVCG, 2024.
%for further details on the proposed architecture. 
%% Data Extraction (Training sample)
%This section extracts the head movement recorded data from the training
%sample.

%Training is performed on the first group of 60 users performing the Stroop
%Color Word Test (SCWT).
clear all
clc

%fileNames.txt is a file containing all the directory paths of the txt file
%containing the head movements recordings

fid = fopen('.\file_paths\fileNames.txt','r');
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

%% Feature Extraction (Training sample)

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

%% Validation data extraction

%validationFiles.txt is a file containing all the directory paths of the txt file
%containing the head movements recordings of the 20 users performing the
%SCWT, representing the validation sample
fid = fopen('.\file_paths\validationFiles.txt','r');
tline = fgetl(fid);
nomeFile=[];

while ischar(tline)
    nome=string(tline);
    nomeFile= [nomeFile; nome];
    tline = fgetl(fid);
end
fclose(fid);

val_coor_disp_phase1 = [];
val_coor_disp_phase2 = [];
val_coor_disp_phase3 = [];

for i=1:size(nomeFile,1)
    [disp_phase1, disp_phase2, disp_phase3] = getFeatures_fixed(nomeFile(i));

    if i == 1
       val_coor_disp_phase1 = disp_phase1; 
    elseif size(val_coor_disp_phase1,1) < size(disp_phase1,1)
        val_coor_disp_phase1 = [val_coor_disp_phase1; val_coor_disp_phase1(end,:,:).*ones(size(disp_phase1,1)-size(val_coor_disp_phase1,1), size(val_coor_disp_phase1,2),size(val_coor_disp_phase1,3))];
        val_coor_disp_phase1 =  cat(3, val_coor_disp_phase1, disp_phase1);
    elseif size(val_coor_disp_phase1,1) > size(disp_phase1,1)
        disp_phase1 = [disp_phase1; disp_phase1(end,:).*ones(size(val_coor_disp_phase1,1)-size(disp_phase1,1),size(val_coor_disp_phase1,2))];
        val_coor_disp_phase1 =  cat(3, val_coor_disp_phase1, disp_phase1);
    else
         val_coor_disp_phase1 =  cat(3, val_coor_disp_phase1, disp_phase1);
    end

    if i == 1
       val_coor_disp_phase2 = disp_phase2;  
    elseif size(val_coor_disp_phase2,1) < size(disp_phase2,1)
        val_coor_disp_phase2 = [val_coor_disp_phase2; val_coor_disp_phase2(end,:,:).*ones(size(disp_phase2,1)-size(val_coor_disp_phase2,1), size(val_coor_disp_phase2,2),size(val_coor_disp_phase2,3))];
        val_coor_disp_phase2 =  cat(3, val_coor_disp_phase2, disp_phase2);
    elseif size(val_coor_disp_phase2,1) > size(disp_phase2,1)
        disp_phase2 = [disp_phase2; disp_phase2(end,:).*ones(size(val_coor_disp_phase2,1)-size(disp_phase2,1),size(val_coor_disp_phase2,2))];
        val_coor_disp_phase2 =  cat(3, val_coor_disp_phase2, disp_phase2);
    else
        val_coor_disp_phase2 =  cat(3, val_coor_disp_phase2, disp_phase2);
    end


    if i == 1
       val_coor_disp_phase3 = disp_phase3; 
    elseif size(val_coor_disp_phase3,1) < size(disp_phase3,1)
        val_coor_disp_phase3 = [val_coor_disp_phase3; val_coor_disp_phase3(end,:,:).*ones(size(disp_phase3,1)-size(val_coor_disp_phase3,1), size(val_coor_disp_phase3,2),size(val_coor_disp_phase3,3))];
        val_coor_disp_phase3 =  cat(3, val_coor_disp_phase3, disp_phase3);
    elseif size(val_coor_disp_phase3,1) > size(disp_phase3,1)
        disp_phase3 = [disp_phase3; disp_phase3(end,:).*ones(size(val_coor_disp_phase3,1)-size(disp_phase3,1),size(val_coor_disp_phase3,2))];
        val_coor_disp_phase3 =  cat(3, val_coor_disp_phase3, disp_phase3);
    else
        val_coor_disp_phase3 =  cat(3, val_coor_disp_phase3, disp_phase3);
    end
end

%% Features Extraction (Validation sample)

val_speed_phase1 = 10*diff(val_coor_disp_phase1);
val_speed_phase2 = 10*diff(val_coor_disp_phase2);
val_speed_phase3 = 10*diff(val_coor_disp_phase3);

val_tot_disp_phase1=[];
val_tot_disp_phase2=[];
val_tot_disp_phase3=[];

for j=1:size(val_coor_disp_phase1,3)
    for i=1:length(val_coor_disp_phase1)
        val_tot_disp_phase1(i,j) = norm(val_coor_disp_phase1(i,:,j));
    end
end

for j=1:size(val_coor_disp_phase2,3)
    for i=1:length(val_coor_disp_phase2)
        val_tot_disp_phase2(i,j) = norm(val_coor_disp_phase2(i,:,j));
    end
end

for j=1:size(val_coor_disp_phase3,3)
    for i=1:length(val_coor_disp_phase3)
        val_tot_disp_phase3(i,j) = norm(val_coor_disp_phase3(i,:,j));
    end
end

magn_speed_phase1 = zeros(length(speed_phase1), size(speed_phase1,3));

for j = 1:size(val_speed_phase1,3)
    for i=1:length(val_speed_phase1)
        val_magn_speed_phase1(i,j) = norm(val_speed_phase1(i,:,j));
    end
end

val_magn_speed_phase2 = zeros(length(val_speed_phase2), 1);

for j = 1:size(val_speed_phase2,3)
    for i=1:size(val_speed_phase2,1)
        val_magn_speed_phase2(i,j) = norm(val_speed_phase2(i,:,j));
    end
end

val_magn_speed_phase3 = zeros(length(val_speed_phase3), 1);

for j = 1:size(val_speed_phase3,3)
    for i=1:size(val_speed_phase3,1)
        val_magn_speed_phase3(i,j) = norm(val_speed_phase3(i,:,j));
    end
end

val_angle_yz_1 = zeros(size(val_coor_disp_phase1,1),size(val_coor_disp_phase1,3));
val_angle_yz_2 = zeros(size(val_coor_disp_phase2,1),size(val_coor_disp_phase2,3));
val_angle_yz_3 = zeros(size(val_coor_disp_phase3,1),size(val_coor_disp_phase3,3));

val_angle_xz_1 = zeros(size(val_coor_disp_phase1,1),size(val_coor_disp_phase1,3));
val_angle_xz_2 = zeros(size(val_coor_disp_phase2,1),size(val_coor_disp_phase2,3));
val_angle_xz_3 = zeros(size(val_coor_disp_phase3,1),size(val_coor_disp_phase3,3));

val_angle_xy_1 = zeros(size(val_coor_disp_phase1,1),size(val_coor_disp_phase1,3));
val_angle_xy_2 = zeros(size(val_coor_disp_phase2,1),size(val_coor_disp_phase2,3));
val_angle_xy_3 = zeros(size(val_coor_disp_phase3,1),size(val_coor_disp_phase3,3));

for i=1:size(val_coor_disp_phase1,3)
    val_angle_yz_1(:,i) = atan(val_coor_disp_phase1(:,2,i)./val_coor_disp_phase1(:,3,i));
    val_angle_yz_2(:,i) = atan(val_coor_disp_phase2(:,2,i)./val_coor_disp_phase2(:,3,i));
    val_angle_yz_3(:,i) = atan(val_coor_disp_phase3(:,2,i)./val_coor_disp_phase3(:,3,i));
    
    val_angle_xz_1(:,i) = atan(val_coor_disp_phase1(:,1,i)./val_coor_disp_phase1(:,3,i));
    val_angle_xz_2(:,i) = atan(val_coor_disp_phase2(:,1,i)./val_coor_disp_phase2(:,3,i));
    val_angle_xz_3(:,i) = atan(val_coor_disp_phase3(:,1,i)./val_coor_disp_phase3(:,3,i));
    
    val_angle_xy_1(:,i) = atan(val_coor_disp_phase1(:,2,i)./val_coor_disp_phase1(:,1,i));
    val_angle_xy_2(:,i) = atan(val_coor_disp_phase2(:,2,i)./val_coor_disp_phase2(:,1,i));
    val_angle_xy_3(:,i) = atan(val_coor_disp_phase3(:,2,i)./val_coor_disp_phase3(:,1,i));
    
end
for i = 1:size(val_angle_xy_3,2)
    val_angle_yz_1(isnan(val_angle_yz_1(:,i)),i) = 0;
    val_angle_xz_1(isnan(val_angle_xz_1(:,i)),i) = 0;
    val_angle_xy_1(isnan(val_angle_xy_1(:,i)),i) = 0;

    val_angle_yz_2(isnan(val_angle_yz_2(:,i)),i) = 0;
    val_angle_xz_2(isnan(val_angle_xz_2(:,i)),i) = 0;
    val_angle_xy_2(isnan(val_angle_xy_2(:,i)),i) = 0;

    val_angle_yz_3(isnan(val_angle_yz_3(:,i)),i) = 0;
    val_angle_xz_3(isnan(val_angle_xz_3(:,i)),i) = 0;
    val_angle_xy_3(isnan(val_angle_xy_3(:,i)),i) = 0;
end

%% SVM training for total displacement

rng(7)
idx_train = randperm(180);
idx_val = randperm(60);

%STFT analysis of total displacement

freq_1 = squeeze(sum(unwrap(angle(stft(tot_disp_phase1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2 = squeeze(sum(unwrap(angle(stft(tot_disp_phase2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_3 = squeeze(sum(unwrap(angle(stft(tot_disp_phase3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

val_freq_1 = squeeze(sum(unwrap(angle(stft(val_tot_disp_phase1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
val_freq_2 = squeeze(sum(unwrap(angle(stft(val_tot_disp_phase2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
val_freq_3 = squeeze(sum(unwrap(angle(stft(val_tot_disp_phase3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

training = [freq_1; freq_2; freq_3];
validation = [val_freq_1; val_freq_2; val_freq_3];

labels_training = [zeros(60*2,1); ones(60,1)];
labels_validation = [zeros(20*2,1); ones(20,1)];

training = training(idx_train,:);
validation = validation(idx_val,:);
labels_training = labels_training(idx_train);
labels_validation = labels_validation(idx_val);

%Training of classifier and accuracy on validation set
Mdl_freq_disp = fitcsvm(training, labels_training, 'OptimizeHyperparameters','all', 'Standardize', true);
pred_freq_disp = predict(Mdl_freq_disp,validation);
accuracy_freq_disp = sum(labels_validation==pred_freq_disp)/size(labels_validation,1)

%% SVM training for displacement of angle xy

%STFT analysis of displacement of angle xy

freq_1_xy = squeeze(sum(unwrap(angle(stft(angle_xy_1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_xy = squeeze(sum(unwrap(angle(stft(angle_xy_2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_3_xy = squeeze(sum(unwrap(angle(stft(angle_xy_3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

val_freq_1_xy = squeeze(sum(unwrap(angle(stft(val_angle_xy_1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
val_freq_2_xy = squeeze(sum(unwrap(angle(stft(val_angle_xy_2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
val_freq_3_xy = squeeze(sum(unwrap(angle(stft(val_angle_xy_3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

training_xy = [freq_1_xy; freq_2_xy; freq_3_xy];
validation_xy = [val_freq_1_xy; val_freq_2_xy; val_freq_3_xy];

labels_training = [zeros(60*2,1); ones(60,1)];
labels_validation = [zeros(20*2,1); ones(20,1)]; 

%Training of classifier and accuracy on validation set
Mdl_freq_xy = fitcsvm(training_xy, labels_training, 'OptimizeHyperparameters','all');
pred_freq_xy = predict(Mdl_freq_xy,validation_xy);
accuracy_freq_xy = sum(labels_validation==pred_freq_xy)/size(labels_validation,1)

%% SVM training for displacement of angle xz

%STFT analysis of displacement of angle xz

freq_1_xz = squeeze(sum(unwrap(angle(stft(angle_xz_1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2_xz = squeeze(sum(unwrap(angle(stft(angle_xz_2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_3_xz = squeeze(sum(unwrap(angle(stft(angle_xz_3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

val_freq_1_xz = squeeze(sum(unwrap(angle(stft(val_angle_xz_1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
val_freq_2_xz = squeeze(sum(unwrap(angle(stft(val_angle_xz_2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
val_freq_3_xz = squeeze(sum(unwrap(angle(stft(val_angle_xz_3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

training_xz = [freq_1_xz; freq_2_xz; freq_3_xz];
validation_xz = [val_freq_1_xz; val_freq_2_xz; val_freq_3_xz];

%Training of classifier and accuracy on validation set

Mdl_freq_xz = fitcsvm(training_xz, labels_training, 'OptimizeHyperparameters','all');
pred_freq_xz = predict(Mdl_freq_xz,validation_xz);
accuracy_freq_xz = sum(labels_validation==pred_freq_xz)/size(labels_validation,1)
%% SVM training for displacement of angle yz

%STFT analysis of displacement of angle yz

freq_1_yz = squeeze(sum(unwrap(angle(stft(angle_yz_1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_2_yz = squeeze(sum(unwrap(angle(stft(angle_yz_2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
freq_3_yz = squeeze(sum(unwrap(angle(stft(angle_yz_3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

val_freq_1_yz = squeeze(sum(unwrap(angle(stft(val_angle_yz_1, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
val_freq_2_yz = squeeze(sum(unwrap(angle(stft(val_angle_yz_2, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';
val_freq_3_yz = squeeze(sum(unwrap(angle(stft(val_angle_yz_3, 'Window', flattopwin(16), 'OverlapLength', 8))),2))';

%Training of classifier and accuracy on validation set

training_yz = [freq_1_yz; freq_2_yz; freq_3_yz];
validation_yz = [val_freq_1_yz; val_freq_2_yz; val_freq_3_yz];

Mdl_freq_yz = fitcsvm(training_yz, labels_training, 'OptimizeHyperparameters','all');
pred_freq_yz = predict(Mdl_freq_yz,validation_yz);
accuracy_freq_yz = sum(labels_validation==pred_freq_yz)/size(labels_validation,1)

%% SVM training for displacement of angle yz 

decision = nan(size(pred_freq_yz,1), size(pred_freq_yz,2));

for i=1:size(pred_freq_yz,1)
    predictions = [pred_freq_xy(i) pred_freq_xz(i) pred_freq_yz(i)];
    if length(find(predictions==1)) >= 2
       decision(i) = 1;
    else
       decision(i) = 0;
    end
end

accuracy_angle_mv = sum(labels_validation==decision)/size(labels_validation,1)

%% SVM training for displacement of speed

%STFT analysis of total speed

freq_1_speed = squeeze(sum(unwrap(angle(stft(magn_speed_phase1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_2_speed = squeeze(sum(unwrap(angle(stft(magn_speed_phase2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
freq_3_speed = squeeze(sum(unwrap(angle(stft(magn_speed_phase3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

val_freq_1_speed = squeeze(sum(unwrap(angle(stft(val_magn_speed_phase1, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
val_freq_2_speed = squeeze(sum(unwrap(angle(stft(val_magn_speed_phase2, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';
val_freq_3_speed = squeeze(sum(unwrap(angle(stft(val_magn_speed_phase3, 'Window', flattopwin(16), 'OverlapLength', 8))), 2))';

%Training of classifier and accuracy on validation set

training_speed = [freq_1_speed; freq_2_speed; freq_3_speed];
validation_speed = [val_freq_1_speed; val_freq_2_speed; val_freq_3_speed];

labels_training = [zeros(60*2,1); ones(60,1)];
labels_validation = [zeros(20*2,1); ones(20,1)];

Mdl_freq_speed = fitcsvm(training_speed, labels_training, 'OptimizeHyperparameters','all', 'Standardize',true);
pred_freq_speed = predict(Mdl_freq_speed,validation_speed);
accuracy_freq_speed = sum(labels_validation==pred_freq_speed)/size(labels_validation,1)
