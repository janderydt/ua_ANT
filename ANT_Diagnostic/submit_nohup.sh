#!/bin/sh

pid=$$;
pgid=$(ps -p $pid -o pgid --no-headers);
pgid=$(echo $pgid | sed 's/ //g');
echo $pgid;

$(

start="addpath /mnt/md0/Ua/cases/ANT; cd /mnt/md0/Ua/cases/ANT/ANT_Diagnostic; ANT_MatlabWrapper($pgid,\"Diagnostic\")";

/usr/local/MATLAB/R2023b/bin/matlab -nodisplay -r "$start; quit;";

)&
