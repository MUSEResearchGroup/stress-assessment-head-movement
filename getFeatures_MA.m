function [disp_phase1,disp_phase2] = getFeatures_MA(fileName)
%This function extracts the temporal variation of head movement along the
%x,y,z axes. fileName is a file .txt conteining all the recorded data of
%one user. 
%The outputs are the coordinates temporal variation divided in the two
%phases of the Mental Arithmetic Test (MA).

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

    final_time = string(duration(timestamp_array(1,1:8))+duration("00:01:15"));
    end_phase1 = strcat(final_time, timestamp_array(1,9:13));

    idx = find(string(timestamp_array) >= end_phase1);
    coord_phase1 = coord_array(1:idx(1),:);
    coord_phase2 = coord_array(idx(1)+1:end,:);

%Mean

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


end