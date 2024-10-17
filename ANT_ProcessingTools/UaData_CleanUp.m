function UaData_CleanUp

exptype = "Inverse";
domain = "ANT_nsmbl";
idrange = [10000,18000];
delims = ["_It","_Yrs"];
suffix_to_delete = "-RestartFile"+delims;% "NewMeshFile.mat"

% find folders
casefolder = "../ANT_"+exptype+"/cases/";
folders = dir(casefolder+domain+"*");
for ii=1:numel(folders)
    folder_name_full = folders(ii).name;
    folder_name_split = strsplit(folder_name_full,"_");
    expid = double(string(folder_name_split{end}));
    for nn=1:size(idrange,1)
        if expid >= idrange(nn,1) && expid <= idrange(nn,2)
            % find files to delete
            for ss=1:numel(suffix_to_delete)
                files_to_delete = dir(folders(ii).folder+"/"+folders(ii).name+...
                    "/"+domain+"_"+exptype+"_"+string(expid)+suffix_to_delete(ss)+"*.mat");
                % sort according to iteration and output year  
                filenums=[];
                for tt=1:numel(files_to_delete)
                    tmp = files_to_delete(tt).name;                    
                    tmp = erase(tmp,[domain+"_"+exptype+"_"+string(expid)+suffix_to_delete(ss),".mat"]);
                    % replace k with .
                    tmp = strrep(tmp,"k",".");
                    filenums(tt) = str2double(tmp);
                end
                [~,I] = sort(filenums,"ascend");
                if numel(I)>1
                    for jj=I(1):I(end-1) % do not delete the last restart file
                        tmp = files_to_delete(I(jj)).folder+"/"+files_to_delete(I(jj)).name;
                        fprintf("Deleting %s\n",tmp);
                        delete(tmp);
                        pause(3);
                    end   
                end
            end
        end
    end

end