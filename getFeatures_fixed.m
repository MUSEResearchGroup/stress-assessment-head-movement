function [disp_phase1, disp_phase2, disp_phase3] = getFeatures_fixed(fileName)
%This function extracts the temporal variation of head movement along the
%x,y,z axes. fileName is a file .txt conteining all the recorded data of
%one user. 
%The outputs are the coordinates temporal variation divided in the three
%phases of the Stroop Color Word Test (SCWT).

fid = fopen(fileName,'r');
tline = fgetl(fid);
coord_array=[];
timestamp = [];

while ischar(tline)
    time = tline(1:10);
    timestamp = [timestamp;time];
    coord=tline(13:end-1);
    coord_new=replace(coord,' ','');
    coord_num=str2double(split(coord_new,','))';
    coord_array=[coord_array;coord_num];
    tline = fgetl(fid);
end
fclose(fid);

timestamp_array = [];

    for i=1:size(timestamp,1)
        timestamp_array = [timestamp_array; strcat(timestamp(i,1:2), ':', timestamp(i,3:4), ':', timestamp(i,5:6), '.', timestamp(i,7:10))];
    end
%Phase 1
    seconds = str2double(timestamp_array(1,7:8))+23;

    if seconds > 59
        seconds_new = seconds-60;

        if seconds_new < 10
           seconds_new = strcat('0', string(seconds_new));
        else
           seconds_new = string(seconds_new); 
        end

        minutes = str2double(timestamp_array(1,4:5))+1;
        if minutes > 59
            minutes_new = minutes-60;
            if minutes_new < 10
                minutes_new = strcat('0', string(minutes_new));
            else
                minutes_new = string(minutes_new);
            end
            hours = str2double(timestamp_array(1,1:2))+1;
            end_phase1 = strcat(string(hours), ':', minutes_new, ':', seconds_new, '.', timestamp_array(1,10:13));
        else
            if minutes < 10
                minutes = strcat('0', string(minutes));
            end

            end_phase1 = strcat(timestamp_array(1,1:2), ':', string(minutes), ':', seconds_new, '.', timestamp_array(1,10:13));
        end
    else
       end_phase1 = strcat(timestamp_array(1,1:2), ':', timestamp_array(1,4:5), ':', string(seconds), '.', timestamp_array(1,10:13));
    end
    
    idx = find(end_phase1<string(timestamp_array));
    
    coord_phase1 = coord_array(1:idx(1),:);
    

%Phase 2
    seconds = str2double(timestamp_array(idx(1)+1,7:8))+37;

    if seconds > 59
        seconds_new = seconds-60;

        if seconds_new < 10
           seconds_new = strcat('0', string(seconds_new));
        else
           seconds_new = string(seconds_new); 
        end

        minutes = str2double(timestamp_array(idx(1)+1,4:5))+1;
        if minutes > 59
            minutes_new = minutes-60;
            if minutes_new < 10
                minutes_new = strcat('0', string(minutes_new));
            else
                minutes_new = string(minutes_new);
            end
            hours = str2double(timestamp_array(1,1:2))+1;
            end_phase2 = strcat(string(hours), ':', minutes_new, ':', seconds_new, '.', timestamp_array(1,10:13));
        else
            if minutes < 10
                minutes = strcat('0', string(minutes));
            end

            end_phase2 = strcat(timestamp_array(idx(1)+1,1:2), ':', string(minutes), ':', seconds_new, '.', timestamp_array(1,10:13));
        end
    else
        end_phase2 = strcat(timestamp_array(idx(1)+1,1:2), ':', timestamp_array(idx(1)+1,4:5), ':', string(seconds), '.', timestamp_array(idx(1)+1,10:13));
    end
    
    idx_2 = find(end_phase2<string(timestamp_array));
    
    coord_phase2 = coord_array(idx(1)+1:idx_2(1),:);


%Phase 3
coord_phase3 = coord_array(idx_2(1)+1:end,:);

disp_phase1 = [];
disp_phase2 = [];
disp_phase3 = [];

for i=1:5:size(coord_phase1,1)
    if i == 1
       disp_phase1 = [disp_phase1; mean(coord_phase1(i:i+4,:))];
    elseif i+4 < size(coord_phase1,1)
       disp_phase1 = [disp_phase1; mean(coord_phase1(i:i+4,:))];
   elseif i == size(coord_phase1,1)
       disp_phase1  = [disp_phase1; coord_phase1(i,:)];
    else
       disp_phase1  = [disp_phase1; mean(coord_phase1(i:end,:))];
    end

end

for i=1:5:size(coord_phase2,1)
    if i == 1
       disp_phase2 = [disp_phase2; mean(coord_phase2(i:i+4,:))];
    elseif i+4 < size(coord_phase2,1)
       disp_phase2 = [disp_phase2; mean(coord_phase2(i:i+4,:))];     
    elseif i == size(coord_phase2,1)
       disp_phase2  = [disp_phase2; coord_phase2(i,:)];
    else
       disp_phase2  = [disp_phase2; mean(coord_phase2(i:end,:))];
    end
end

for i=1:5:size(coord_phase3,1)
    if i == 1
       disp_phase3 = [disp_phase3; mean(coord_phase3(i:i+4,:))];
    elseif i+4 < size(coord_phase3,1)
       disp_phase3 = [disp_phase3; mean(coord_phase3(i:i+4,:))];
   elseif i == size(coord_phase3,1)
       disp_phase3  = [disp_phase3; coord_phase3(i,:)];
    else
       disp_phase3  = [disp_phase3; mean(coord_phase3(i:end,:))];
    end

end

end