%Stress Assessment through Head Movements
%This script allows to replicate the ANOVA tests performed on the data of
%the Stroop Color Word Test (SCWT).


%% Data Extraction
%This section extracts the head movement recorded data.
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

%This section read the recordings and saved in the coor_disp_phase1,
%coor_disp_phase2 and coor_disp_phase3 variables the x,y,z variation of the
%head movement
%These variables are array with dimension:
%number of samples of phase x 3 (x,y,z) x 60 (number of users)

coor_disp_phase1 = [];
coor_disp_phase2 = [];
coor_disp_phase3 = [];

for i=1:size(nomeFile,1)
    [disp_phase1, disp_phase2, disp_phase3] = getFeatures_fixed(nomeFile(i));

    if i == 1
       coor_disp_phase1 = disp_phase1; 
    elseif size(coor_disp_phase1,1) < size(disp_phase1,1)
        coor_disp_phase1 = [coor_disp_phase1; zeros(size(disp_phase1,1)-size(coor_disp_phase1,1), size(coor_disp_phase1,2),size(coor_disp_phase1,3))];
        coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    elseif size(coor_disp_phase1,1) > size(disp_phase1,1)
        disp_phase1 = [disp_phase1; zeros(size(coor_disp_phase1,1)-size(disp_phase1,1),size(coor_disp_phase1,2))];
        coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    else
         coor_disp_phase1 =  cat(3, coor_disp_phase1, disp_phase1);
    end

    if i == 1
       coor_disp_phase2 = disp_phase2;  
    elseif size(coor_disp_phase2,1) < size(disp_phase2,1)
        coor_disp_phase2 = [coor_disp_phase2; zeros(size(disp_phase2,1)-size(coor_disp_phase2,1), size(coor_disp_phase2,2),size(coor_disp_phase2,3))];
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    elseif size(coor_disp_phase2,1) > size(disp_phase2,1)
        disp_phase2 = [disp_phase2; zeros(size(coor_disp_phase2,1)-size(disp_phase2,1),size(coor_disp_phase2,2))];
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    else
        coor_disp_phase2 =  cat(3, coor_disp_phase2, disp_phase2);
    end


    if i == 1
       coor_disp_phase3 = disp_phase3; 
    elseif size(coor_disp_phase3,1) < size(disp_phase3,1)
        coor_disp_phase3 = [coor_disp_phase3; zeros(size(disp_phase3,1)-size(coor_disp_phase3,1), size(coor_disp_phase3,2),size(coor_disp_phase3,3))];
        coor_disp_phase3 =  cat(3, coor_disp_phase3, disp_phase3);
    elseif size(coor_disp_phase3,1) > size(disp_phase3,1)
        disp_phase3 = [disp_phase3; zeros(size(coor_disp_phase3,1)-size(disp_phase3,1),size(coor_disp_phase3,2))];
        coor_disp_phase3 =  cat(3, coor_disp_phase3, disp_phase3);
    else
        coor_disp_phase3 =  cat(3, coor_disp_phase3, disp_phase3);
    end
end

%% Feature Extraction

%Total Displacement
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

%Entity of displacement

abs_diff_phase1 = abs(diff(coor_disp_phase1));
abs_diff_phase2 = abs(diff(coor_disp_phase2));
abs_diff_phase3 = abs(diff(coor_disp_phase3));

%Sign of displacement

sign_disp_phase1 = sign(diff(coor_disp_phase1));
sign_disp_phase2 = sign(diff(coor_disp_phase2));
sign_disp_phase3 = sign(diff(coor_disp_phase3));

%Speed

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

%Angles

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

%Angular entity of displacement

diff_xy_1 = abs(diff(angle_xy_1(2:end,:)));
diff_xy_2 = abs(diff(angle_xy_2(2:end,:)));
diff_xy_3 = abs(diff(angle_xy_3(2:end,:)));

diff_xz_1 = abs(diff(angle_xz_1(2:end,:)));
diff_xz_2 = abs(diff(angle_xz_2(2:end,:)));
diff_xz_3 = abs(diff(angle_xz_3(2:end,:)));

diff_yz_1 = abs(diff(angle_yz_1(2:end,:)));
diff_yz_2 = abs(diff(angle_yz_2(2:end,:)));
diff_yz_3 = abs(diff(angle_yz_3(2:end,:)));


%Angular speed

speed_angle_xz_1 = 10*diff(angle_xz_1(2:end,:));
speed_angle_xz_2 = 10*diff(angle_xz_2(2:end,:));
speed_angle_xz_3 = 10*diff(angle_xz_3(2:end,:));

speed_angle_xy_1 = 10*diff(angle_xy_1(2:end,:));
speed_angle_xy_2 = 10*diff(angle_xy_2(2:end,:));
speed_angle_xy_3 = 10*diff(angle_xy_3(2:end,:));

speed_angle_yz_1 = 10*diff(angle_yz_1(2:end,:));
speed_angle_yz_2 = 10*diff(angle_yz_2(2:end,:));
speed_angle_yz_3 = 10*diff(angle_yz_3(2:end,:));

%Sign of angles

sign_angle_xz_1 = sign(10*diff(angle_xz_1(2:end,:)));
sign_angle_xz_2 = sign(10*diff(angle_xz_2(2:end,:)));
sign_angle_xz_3 = sign(10*diff(angle_xz_3(2:end,:)));

sign_angle_xy_1 = sign(10*diff(angle_xy_1(2:end,:)));
sign_angle_xy_2 = sign(10*diff(angle_xy_2(2:end,:)));
sign_angle_xy_3 = sign(10*diff(angle_xy_3(2:end,:)));

sign_angle_yz_1 = sign(10*diff(angle_yz_1(2:end,:)));
sign_angle_yz_2 = sign(10*diff(angle_yz_2(2:end,:)));
sign_angle_yz_3 = sign(10*diff(angle_yz_3(2:end,:)));

%% ANOVA TEST: HYPOTHESIS CHECK
%This section evaluates if the selected features respects the hypothesis
%for the ANOVA
Features = ["X coordinate", "Y coordinate", "Z coordinate", "Total displacement", "X Entity of displacement", "Y Entity of displacement", "Z Entity of displacement", "X sign", "Y sign", "Z sign", "Speed along x axis", "Speed along y axis", "Speed along z axis", "Total speed", "Angle xy", "Angle yz", "Angle xz", "Entity of displacement xy", "Entity of displacement yz", "Entity of displacement xz", "Angular speed xy", "Angular speed yz", "Angular speed xz", "Sign xy", "Sign yz", "Sign xz"]';
Homoschedasticity = nan(size(Features,1),1);
Normality= nan(size(Features,1),1);

[Normality(1), Homoschedasticity(1)] = hypothesis_check(coor_disp_phase1(:,1,:), coor_disp_phase2(:,1,:), coor_disp_phase3(:,1,:));
[Normality(2), Homoschedasticity(2)] = hypothesis_check(coor_disp_phase1(:,2,:), coor_disp_phase2(:,2,:), coor_disp_phase3(:,2,:));
[Normality(3), Homoschedasticity(3)] = hypothesis_check(coor_disp_phase1(:,3,:), coor_disp_phase2(:,3,:), coor_disp_phase3(:,3,:));

[Normality(4), Homoschedasticity(4)] = hypothesis_check(tot_disp_phase1, tot_disp_phase2, tot_disp_phase3);

[Normality(5), Homoschedasticity(5)] = hypothesis_check(abs_diff_phase1(:,1,:), abs_diff_phase2(:,1,:), abs_diff_phase3(:,1,:));
[Normality(6), Homoschedasticity(6)] = hypothesis_check(abs_diff_phase1(:,2,:), abs_diff_phase2(:,2,:), abs_diff_phase3(:,2,:));
[Normality(7), Homoschedasticity(7)] = hypothesis_check(abs_diff_phase1(:,3,:), abs_diff_phase2(:,3,:), abs_diff_phase3(:,3,:));

[Normality(8), Homoschedasticity(8)] = hypothesis_check(sign_disp_phase1(:,1,:), sign_disp_phase2(:,1,:), sign_disp_phase3(:,1,:));
[Normality(9), Homoschedasticity(9)] = hypothesis_check(sign_disp_phase1(:,2,:), sign_disp_phase2(:,2,:), sign_disp_phase3(:,2,:));
[Normality(10), Homoschedasticity(10)] = hypothesis_check(sign_disp_phase1(:,3,:), sign_disp_phase2(:,3,:), sign_disp_phase3(:,3,:));

[Normality(11), Homoschedasticity(11)] = hypothesis_check(speed_phase1(:,1,:), speed_phase2(:,1,:), speed_phase3(:,1,:));
[Normality(12), Homoschedasticity(12)] = hypothesis_check(speed_phase1(:,2,:), speed_phase2(:,2,:), speed_phase3(:,2,:));
[Normality(13), Homoschedasticity(13)] = hypothesis_check(speed_phase1(:,3,:), speed_phase2(:,3,:), speed_phase3(:,3,:));

[Normality(14), Homoschedasticity(14)] = hypothesis_check(magn_speed_phase1, magn_speed_phase2, magn_speed_phase3);

[Normality(15), Homoschedasticity(15)] = hypothesis_check(angle_xy_1, angle_xy_2, angle_xy_3);
[Normality(16), Homoschedasticity(16)] = hypothesis_check(angle_yz_1, angle_yz_2, angle_yz_3);
[Normality(17), Homoschedasticity(17)] = hypothesis_check(angle_xz_1, angle_xz_1, angle_xz_3);

[Normality(18), Homoschedasticity(18)] = hypothesis_check(diff_xy_1, diff_xy_2, diff_xy_3);
[Normality(19), Homoschedasticity(19)] = hypothesis_check(diff_yz_1, diff_yz_2, diff_yz_3);
[Normality(20), Homoschedasticity(20)] = hypothesis_check(diff_xz_1, diff_xz_1, diff_xz_3);

[Normality(21), Homoschedasticity(21)] = hypothesis_check(speed_angle_xy_1, speed_angle_xy_2, speed_angle_xy_3);
[Normality(22), Homoschedasticity(22)] = hypothesis_check(speed_angle_yz_1, speed_angle_yz_2, speed_angle_yz_3);
[Normality(23), Homoschedasticity(23)] = hypothesis_check(speed_angle_xz_1, speed_angle_xz_1, speed_angle_xz_3);

[Normality(24), Homoschedasticity(24)] = hypothesis_check(sign_angle_xy_1, sign_angle_xy_2, sign_angle_xy_3);
[Normality(25), Homoschedasticity(25)] = hypothesis_check(sign_angle_yz_1, sign_angle_yz_2, sign_angle_yz_3);
[Normality(26), Homoschedasticity(26)] = hypothesis_check(sign_angle_xz_1, sign_angle_xz_1, sign_angle_xz_3);

Hypothesis = [Features Normality Homoschedasticity];
%% ANOVA TEST: Displacement along X, Y, Z coordinates

Observations_x = [];
Observations_y = [];
Observations_z = [];

for i=1:size(coor_disp_phase1,3)
    Observations_x = [Observations_x, coor_disp_phase1(:,1,i)', coor_disp_phase2(:,1,i)', coor_disp_phase3(:,1,i)'];
    Observations_y = [Observations_y, coor_disp_phase1(:,2,i)', coor_disp_phase2(:,2,i)', coor_disp_phase3(:,2,i)'];
    Observations_z = [Observations_z, coor_disp_phase1(:,3,i)', coor_disp_phase2(:,3,i)', coor_disp_phase3(:,3,i)'];
end

grp1 = repmat("Fase1", [1 size(coor_disp_phase1,1)]);
grp2 = repmat("Fase2", [1 size(coor_disp_phase2,1)]);
grp3 = repmat("Fase3", [1 size(coor_disp_phase3,1)]);
grp = repmat([grp1 grp2 grp3], [1 size(coor_disp_phase3,3)]);

Partecipants = string(zeros(length(coor_disp_phase1)+length(coor_disp_phase2)+length(coor_disp_phase3), size(coor_disp_phase1,3)));

for i=1:size(coor_disp_phase1,3)
    Partecipants(:,i) = string(i);
end


Partecipants = Partecipants(:);

[p_x,tbl_x,stats_x] = anovan(Observations_x,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Test Phase','Subject'},'alpha',0.05);
fx=figure;
[c_x,m_x,h_x,gnames_x] = multcompare(stats_x);

[p_y,tbl_y,stats_y] = anovan(Observations_y,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Test Phase','Subject'},'alpha',0.05);
figure
[c_y,m_y,h_y,gnames_y] = multcompare(stats_y);

[p_z,tbl_z,stats_z] = anovan(Observations_z,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Test Phase','Subject'},'alpha',0.05);
fz=figure;
[c_z,m_z,h_z,gnames_z] = multcompare(stats_z);

%% Total displacement
Observations_disp=[];

for i=1:size(coor_disp_phase1,3)
    Observations_disp = [Observations_disp,tot_disp_phase1(:,i)',tot_disp_phase2(:,i)',tot_disp_phase3(:,i)'];
end

grp1 = repmat("Ph1", [1 size(coor_disp_phase1,1)]);
grp2 = repmat("Ph2", [1 size(coor_disp_phase2,1)]);
grp3 = repmat("Ph3", [1 size(coor_disp_phase3,1)]);
grp = repmat([grp1 grp2 grp3], [1 size(coor_disp_phase3,3)]);

Partecipants = string(zeros(length(tot_disp_phase1)+length(tot_disp_phase2)+length(tot_disp_phase3), size(tot_disp_phase1,3)));

for i=1:size(tot_disp_phase1,2)
    Partecipants(:,i) = string(i);
end

Partecipants = reshape(Partecipants, [1 length(Observations_disp)]);
[p_disp,tbl_disp,stats_disp] = anovan(Observations_disp,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Test Phase','Subject'},'alpha',0.05);

f_disp=figure;
[c_disp,m_disp,h_disp,gnames_disp] = multcompare(stats_disp);

%% Entity of displacement

abs_diff_x = [];
abs_diff_y = [];
abs_diff_z = [];

Partecipants_x = [];
Partecipants_y = [];
Partecipants_z = [];

grp1 = repmat("Ph1", 1, size(abs_diff_phase1,1));
grp2 = repmat("Ph2", 1, size(abs_diff_phase2,1));
grp3 = repmat("Ph3", 1, size(abs_diff_phase3,1));
grp= repmat([grp1 grp2 grp3],1,size(abs_diff_phase3,3));

for i=1:size(abs_diff_phase1,3)
    abs_diff_x = [abs_diff_x, abs_diff_phase1(:,1,i)', abs_diff_phase2(:,1,i)', abs_diff_phase3(:,1,i)'];
    Partecipants_x = [Partecipants_x, repmat(string(i),1, size(abs_diff_phase1,1)+size(abs_diff_phase2,1)+size(abs_diff_phase3,1))];

    abs_diff_y = [abs_diff_y, abs_diff_phase1(:,2,i)', abs_diff_phase2(:,2,i)', abs_diff_phase3(:,2,i)'];
    Partecipants_y = [Partecipants_y, repmat(string(i),1, size(abs_diff_phase1,1)+size(abs_diff_phase2,1)+size(abs_diff_phase3,1))];

    abs_diff_z = [abs_diff_z, abs_diff_phase1(:,3,i)', abs_diff_phase2(:,3,i)', abs_diff_phase3(:,3,i)'];
    Partecipants_z = [Partecipants_z, repmat(string(i),1, size(abs_diff_phase1,1)+size(abs_diff_phase2,1)+size(abs_diff_phase3,1))];
end

[p_abs_x,tbl_abs_x,stats_abs_x] = anovan(abs_diff_x,{grp,Partecipants_x},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_abs_y,tbl_abs_y,stats_abs_y] = anovan(abs_diff_y,{grp,Partecipants_y},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_abs_z,tbl_abs_z,stats_abs_z] = anovan(abs_diff_z,{grp,Partecipants_z},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);

fx=figure;
[c_x,m1_x,h_x,gnames_x] = multcompare(stats_abs_x);

fy=figure;
[c_y,m1_y,h_y,gnames_y] = multcompare(stats_abs_y);

fz=figure;
[c_z,m1_z,h_z,gnames_z] = multcompare(stats_abs_z);

%% ANOVA TEST: Sign of displacement
sign_diff_x = [];
sign_diff_y = [];
sign_diff_z = [];

Partecipants_x = [];
Partecipants_y = [];
Partecipants_z = [];

for i=1:size(sign_disp_phase1,3)
    sign_diff_x = [sign_diff_x, sign_disp_phase1(:,1,i)', sign_disp_phase2(:,1,i)', sign_disp_phase3(:,1,i)'];
    Partecipants_x = [Partecipants_x, repmat(string(i),1, size(sign_disp_phase1,1)+size(sign_disp_phase2,1)+size(sign_disp_phase3,1))];

    sign_diff_y = [sign_diff_y, sign_disp_phase1(:,2,i)', sign_disp_phase2(:,2,i)', sign_disp_phase3(:,2,i)'];
    Partecipants_y = [Partecipants_y, repmat(string(i),1, size(sign_disp_phase1,1)+size(sign_disp_phase2,1)+size(sign_disp_phase3,1))];

    sign_diff_z = [sign_diff_z, sign_disp_phase1(:,3,i)', sign_disp_phase2(:,3,i)', sign_disp_phase3(:,3,i)'];
    Partecipants_z = [Partecipants_z, repmat(string(i),1, size(sign_disp_phase1,1)+size(sign_disp_phase2,1)+size(sign_disp_phase3,1))];
end

grp1 = repmat("Ph1", 1, size(sign_disp_phase1,1));
grp2 = repmat("Ph2", 1, size(sign_disp_phase2,1));
grp3 = repmat("Ph3", 1, size(sign_disp_phase3,1));
grp= repmat([grp1 grp2 grp3],1,size(sign_disp_phase3,3));

[p_sign_x,tbl_sign_x,stats_sign_x] = anovan(sign_diff_x,{grp,Partecipants_x},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_sign_y,tbl_sign_y,stats_sign_y] = anovan(sign_diff_y,{grp,Partecipants_y},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_sign_z,tbl_sign_z,stats_sign_z] = anovan(sign_diff_z,{grp,Partecipants_z},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);

fx = figure;
[c_xs,m1_xs,h_xs,gnames_xs] = multcompare(stats_sign_x);

fy=figure;
[c_ys,m1_ys,h_ys,gnames_ys] = multcompare(stats_sign_y);

fz=figure; 
[c_zs,m1_zs,h_zs,gnames_zs] = multcompare(stats_sign_z);

%% ANOVA TEST: Speed along x, y, z axes

Partecipants_x = [];
Partecipants_y = [];
Partecipants_z = [];

speed_x = [];
speed_y = [];
speed_z = [];

for i=1:size(speed_phase1,3)
    speed_x = [speed_x, speed_phase1(:,1,i)', speed_phase2(:,1,i)', speed_phase3(:,1,i)'];
    Partecipants_x = [Partecipants_x, repmat(string(i),1, size(speed_phase1,1)+size(speed_phase2,1)+size(speed_phase3,1))];

    speed_y = [speed_y, speed_phase1(:,2,i)', speed_phase2(:,2,i)', speed_phase3(:,2,i)'];
    Partecipants_y = [Partecipants_y, repmat(string(i),1, size(speed_phase1,1)+size(speed_phase2,1)+size(speed_phase3,1))];

    speed_z = [speed_z, speed_phase1(:,3,i)', speed_phase2(:,3,i)', speed_phase3(:,3,i)'];
    Partecipants_z = [Partecipants_z, repmat(string(i),1, size(speed_phase1,1)+size(speed_phase2,1)+size(speed_phase3,1))];
end

grp1 = repmat("Fase1", 1, size(speed_phase1,1));
grp2 = repmat("Fase2", 1, size(speed_phase2,1));
grp3 = repmat("Fase3", 1, size(speed_phase3,1));
grp= repmat([grp1 grp2 grp3],1,size(speed_phase1,3));

[p_speed_x,tbl_speed_x,stats_speed_x] = anovan(speed_x,{grp,Partecipants_x},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_speed_y,tbl_speed_y,stats_speed_y] = anovan(speed_y,{grp,Partecipants_y},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_speed_z,tbl_speed_z,stats_speed_z] = anovan(speed_z,{grp,Partecipants_z},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);

fx=figure
[c_x_speed,m1_x_speed,h_x_speed,gnames_x_speed] = multcompare(stats_speed_x);

fy=figure;
[c_y_speed,m1_y_speed,h_y_speed,gnames_y_speed] = multcompare(stats_speed_y);

fz=figure;
[c_z_speed,m1_z_speed,h_z_speed,gnames_z_speed] = multcompare(stats_speed_z);

%% ANOVA TEST: Total speed

Partecipants = [];

speed_1 = magn_speed_phase1(:)';
speed_2 = magn_speed_phase2(:)';
speed_3 = magn_speed_phase3(:)';

Observations_speed = [speed_1 speed_2 speed_3];

grp_1 = repmat("Fase1", 1, length(speed_1));
grp_2 = repmat("Fase2", 1, length(speed_2));
grp_3 = repmat("Fase3", 1, length(speed_3));
grp = [grp_1 grp_2 grp_3];

Partecipants1 = string(zeros(size(magn_speed_phase1,1),size(magn_speed_phase1,2)));
Partecipants2 = string(zeros(size(magn_speed_phase2,1),size(magn_speed_phase2,2)));
Partecipants3 = string(zeros(size(magn_speed_phase3,1),size(magn_speed_phase3,2)));


for i=1:size(magn_speed_phase1,2)
    Partecipants1(:,i) = string(i);
    Partecipants2(:,i) = string(i);
    Partecipants3(:,i) = string(i);
end

Partecipants1 = Partecipants1(:)';
Partecipants2 = Partecipants2(:)';
Partecipants3 = Partecipants3(:)';

Partecipants = [Partecipants1 Partecipants2 Partecipants3];

[p_speed,tbl_speed,stats_speed] = anovan(Observations_speed,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);

f = figure;
[c_speed,m1_speed,h_speed,gnames_speed] = multcompare(stats_speed);

%% ANOVA TEST: Angular displacement

reshaped_xy_1 = angle_xy_1(:)';
reshaped_xy_2 = angle_xy_2(:)';
reshaped_xy_3 = angle_xy_3(:)';

reshaped_xz_1 = angle_xz_1(:)';
reshaped_xz_2 = angle_xz_2(:)';
reshaped_xz_3 = angle_xz_3(:)';

reshaped_yz_1 = angle_yz_1(:)';
reshaped_yz_2 = angle_yz_2(:)';
reshaped_yz_3 = angle_yz_3(:)';

% reshaped_xy_1 = reshape(angle_xy_1(2:end,:), [1 (size(angle_xy_1,1)-1)*size(angle_xy_1,2)]);
% reshaped_xy_2 = reshape(angle_xy_2(2:end,:), [1 (size(angle_xy_2,1)-1)*size(angle_xy_2,2)]);
% reshaped_xy_3 = reshape(angle_xy_3(2:end,:), [1 (size(angle_xy_3,1)-1)*size(angle_xy_3,2)]);
% 
% reshaped_xz_1 = reshape(angle_xz_1(2:end,:), [1 (size(angle_xz_1,1)-1)*size(angle_xz_1,2)]);
% reshaped_xz_2 = reshape(angle_xz_2(2:end,:), [1 (size(angle_xz_2,1)-1)*size(angle_xz_2,2)]);
% reshaped_xz_3 = reshape(angle_xz_3(2:end,:), [1 (size(angle_xz_3,1)-1)*size(angle_xz_3,2)]);
% 
% reshaped_yz_1 = reshape(angle_yz_1(2:end,:), [1 (size(angle_yz_1,1)-1)*size(angle_yz_1,2)]);
% reshaped_yz_2 = reshape(angle_yz_2(2:end,:), [1 (size(angle_yz_2,1)-1)*size(angle_yz_2,2)]);
% reshaped_yz_3 = reshape(angle_yz_3(2:end,:), [1 (size(angle_yz_3,1)-1)*size(angle_yz_3,2)]);

Observations_xy =  [reshaped_xy_1 reshaped_xy_2 reshaped_xy_3];
Observations_xz =  [reshaped_xz_1 reshaped_xz_2 reshaped_xz_3];
Observations_yz =  [reshaped_yz_1 reshaped_yz_2 reshaped_yz_3];
 
grp1 = repmat("Ph1", 1, size(angle_xy_1,1)*size(angle_xy_1,2));
grp2 = repmat("Ph2", 1, size(angle_xy_2,1)*size(angle_xy_2,2));
grp3 = repmat("Ph3", 1, size(angle_xy_3,1)*size(angle_xy_1,2));
grp = [grp1 grp2 grp3];


Partecipants_1 = [];
Partecipants_2 = [];
Partecipants_3 = [];

for i=1:size(angle_xy_1,2)
    Partecipants_1 = [Partecipants_1 repmat(string(i), [1 length(angle_xy_1)])];
    Partecipants_2 = [Partecipants_2 repmat(string(i), [1 length(angle_xy_2)])];
    Partecipants_3 = [Partecipants_3 repmat(string(i), [1 length(angle_xy_3)])];
end

Partecipants = [Partecipants_1 Partecipants_2 Partecipants_3];

[p_angle_xz,tbl_angle_xz,stats_angle_xz] = anovan(Observations_xz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_xy,tbl_angle_xy,stats_angle_xy] = anovan(Observations_xy,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_yz,tbl_angle_yz,stats_angle_yz] = anovan(Observations_yz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);

fxz = figure;
[c_xz_angle,m1_xz_angle,h_xz_angle,gnames_xz_angle] = multcompare(stats_angle_xz);

fyz=figure;
[c_yz_angle,m1_yz_angle,h_yz_angle,gnames_yz_angle] = multcompare(stats_angle_yz);

fxy=figure;
[c_xy_angle,m1_xy_angle,h_xy_angle,gnames_xy_angle] = multcompare(stats_angle_xy);

%% ANOVA TEST: Angular entity of displacement

diff_xy_1 = diff_xy_1(:)';
diff_xy_2 = diff_xy_2(:)';
diff_xy_3 = diff_xy_3(:)';
Observations_xy = [diff_xy_1 diff_xy_2 diff_xy_3];

diff_xz_1 = diff_xz_1(:)';
diff_xz_2 = diff_xz_2(:)';
diff_xz_3 = diff_xz_3(:)';
Observations_xz = [diff_xz_1 diff_xz_2 diff_xz_3];

diff_yz_1 = diff_yz_1(:)';
diff_yz_2 = diff_yz_2(:)';
diff_yz_3 = diff_yz_3(:)';
Observations_yz = [diff_yz_1 diff_yz_2 diff_yz_3];

grp1 = repmat("Ph1", 1, size(diff_xy_1,2));
grp2 = repmat("Ph2", 1, size(diff_xy_2,2));
grp3 = repmat("Ph3", 1, size(diff_xy_3,2));
grp = [grp1 grp2 grp3];

Partecipants_1 = [];
Partecipants_2 = [];
Partecipants_3 = [];

for i=1:size(angle_xy_1,2)
    Partecipants_1 = [Partecipants_1 repmat(string(i), [1 length(angle_xy_1)-2])];
    Partecipants_2 = [Partecipants_2 repmat(string(i), [1 length(angle_xy_2)-2])];
    Partecipants_3 = [Partecipants_3 repmat(string(i), [1 length(angle_xy_3)-2])];
end

Partecipants = [Partecipants_1 Partecipants_2 Partecipants_3];

[p_angle_xz,tbl_angle_xz,stats_angle_xz] = anovan(Observations_xz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_xy,tbl_angle_xy,stats_angle_xy] = anovan(Observations_xy,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_yz,tbl_angle_yz,stats_angle_yz] = anovan(Observations_yz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);

fxz = figure;
[c_xz_angle,m1_xz_angle,h_xz_angle,gnames_xz_angle] = multcompare(stats_angle_xz);

fyz=figure;
[c_yz_angle,m1_yz_angle,h_yz_angle,gnames_yz_angle] = multcompare(stats_angle_yz);

fxy=figure;
[c_xy_angle,m1_xy_angle,h_xy_angle,gnames_xy_angle] = multcompare(stats_angle_xy);

%% ANOVA TEST: Angular speed

reshaped_xy_1 = speed_angle_xy_1(:)';
reshaped_xy_2 = speed_angle_xy_2(:)';
reshaped_xy_3 = speed_angle_xy_3(:)';

reshaped_xz_1 = speed_angle_xz_1(:)';
reshaped_xz_2 = speed_angle_xz_2(:)';
reshaped_xz_3 = speed_angle_xz_3(:)';

reshaped_yz_1 = speed_angle_yz_1(:)';
reshaped_yz_2 = speed_angle_yz_2(:)';
reshaped_yz_3 = speed_angle_yz_3(:)';

Observations_xy =  [reshaped_xy_1 reshaped_xy_2 reshaped_xy_3];
Observations_xz =  [reshaped_xz_1 reshaped_xz_2 reshaped_xz_3];
Observations_yz =  [reshaped_yz_1 reshaped_yz_2 reshaped_yz_3];

grp1 = repmat("Ph1", 1, size(speed_angle_xy_1,1)*size(speed_angle_xy_1,2));
grp2 = repmat("Ph2", 1, size(speed_angle_xy_2,1)*size(speed_angle_xy_2,2));
grp3 = repmat("Ph3", 1, size(speed_angle_xy_3,1)*size(speed_angle_xy_3,2));
grp= [grp1 grp2 grp3];

Partecipants_1 = [];
Partecipants_2 = [];
Partecipants_3 = [];

for i=1:size(angle_xy_1,2)
    Partecipants_1 = [Partecipants_1 repmat(string(i), [1 length(speed_angle_xy_1)])];
    Partecipants_2 = [Partecipants_2 repmat(string(i), [1 length(speed_angle_xy_2)])];
    Partecipants_3 = [Partecipants_3 repmat(string(i), [1 length(speed_angle_xy_3)])];
end

Partecipants = [Partecipants_1 Partecipants_2 Partecipants_3];

[p_angle_xz,tbl_angle_xz,stats_angle_xz] = anovan(Observations_xz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_xy,tbl_angle_xy,stats_angle_xy] = anovan(Observations_xy,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_yz,tbl_angle_yz,stats_angle_yz] = anovan(Observations_yz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);


[c_xz_angle,m1_xz_angle,h_xz_angle,gnames_xz_angle] = multcompare(stats_angle_xz);
figure
[c_yz_angle,m1_yz_angle,h_yz_angle,gnames_yz_angle] = multcompare(stats_angle_yz);
figure
[c_xy_angle,m1_xy_angle,h_xy_angle,gnames_xy_angle] = multcompare(stats_angle_xy);

%% ANOVA TEST: Sign of angles

reshaped_xy_1 = sign_angle_xy_1(:)';
reshaped_xy_2 = sign_angle_xy_2(:)';
reshaped_xy_3 = sign_angle_xy_3(:)';

reshaped_xz_1 = sign_angle_xz_1(:)';
reshaped_xz_2 = sign_angle_xz_2(:)';
reshaped_xz_3 = sign_angle_xz_3(:)';

reshaped_yz_1 = sign_angle_yz_1(:)';
reshaped_yz_2 = sign_angle_yz_2(:)';
reshaped_yz_3 = sign_angle_yz_3(:)';

Observations_xy =  [reshaped_xy_1 reshaped_xy_2 reshaped_xy_3];
Observations_xz =  [reshaped_xz_1 reshaped_xz_2 reshaped_xz_3];
Observations_yz =  [reshaped_yz_1 reshaped_yz_2 reshaped_yz_3];

grp1 = repmat("Ph1", 1, size(sign_angle_xy_1,1)*size(sign_angle_xy_1,2));
grp2 = repmat("Ph2", 1, size(sign_angle_xy_2,1)*size(sign_angle_xy_1,2));
grp3 = repmat("Ph3", 1, size(sign_angle_xy_3,1)*size(sign_angle_xy_1,2));
grp= [grp1 grp2 grp3];

Partecipants_1 = [];
Partecipants_2 = [];
Partecipants_3 = [];

for i=1:size(angle_xy_1,2)
    Partecipants_1 = [Partecipants_1 repmat(string(i), [1 length(sign_angle_xy_1)])];
    Partecipants_2 = [Partecipants_2 repmat(string(i), [1 length(sign_angle_xy_2)])];
    Partecipants_3 = [Partecipants_3 repmat(string(i), [1 length(sign_angle_xy_3)])];
end

Partecipants = [Partecipants_1 Partecipants_2 Partecipants_3];

[p_angle_xz,tbl_angle_xz,stats_angle_xz] = anovan(Observations_xz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_xy,tbl_angle_xy,stats_angle_xy] = anovan(Observations_xy,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);
[p_angle_yz,tbl_angle_yz,stats_angle_yz] = anovan(Observations_yz,{grp,Partecipants},'model',[1 0;0 1],'varnames',{'Phase','Subject'},'alpha',0.05);


[c_xz_angle,m1_xz_angle,h_xz_angle,gnames_xz_angle] = multcompare(stats_angle_xz);
% figure
fyz = figure;
[c_yz_angle,m1_yz_angle,h_yz_angle,gnames_yz_angle] = multcompare(stats_angle_yz);

figure
[c_xy_angle,m1_xy_angle,h_xy_angle,gnames_xy_angle] = multcompare(stats_angle_xy);

%% Plot statistical analysis of total displacement
f__disp = figure;
axes1 = axes;
hold(axes1,'on');
limits = [0.0559115 0.0564547; 0.0553244 0.0557361; 0.0614097 0.061779];

for i = 1:3
    y = i;

    if i == 1 
        color = [0 0 1 0.2];
    elseif i == 2 
        color = [1 0 0 0.2];
    else
        color = [0 1 0 0.2];
    end
    rectangle('Position', [0,limits(i,1), 4, limits(i,2)-limits(i,1)], 'FaceColor', color, 'EdgeColor', [1 1 1])
    

end

for i = 1:3
    y = i;
    
    if i == 1 
        color = 'blue';
    elseif i == 2
        color = 'red';
    else
        color = '#154734';
    end
    
    hold on
    line([y,y],limits(i,:), 'Color',color, 'LineWidth', 6)
    hold on
    h=stem(y, m_disp(i,1), 'filled', 'LineStyle','None', 'MarkerFaceColor',color, 'MarkerEdgeColor', color, 'MarkerSize', 20);
    hold on
    hbase = h.BaseLine; 
    hbase.Visible = 'off';
    
end
% Create ylabel
% ylabel('Phase');
% Create xlabel
ylabel('Total Displacement (m)');
%xlim(axes1,[2 6.3]);
xlim(axes1,[0.5 3.5]);
ylim(axes1,[0.055 0.062])
hold(axes1,'off');
set(axes1,'FontName','Times','FontSize',80);
% f_disp.yticklabels = ["Phase 1" "Phase 2" "Phase 3"];
xticklabels(["Phase 1" "Phase 2" "Phase 3"])

