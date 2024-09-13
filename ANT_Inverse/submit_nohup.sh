#!/bin/sh

pid=$$;
pgid=$(ps -p $pid -o pgid --no-headers);
pgid=$(echo $pgid | sed 's/ //g');
echo $pgid;

host=$(hostname)

if [ "$host" = "C23000100" ]; then
    start="addpath(\"/mnt/md0/Ua/cases/ANT\"); cd /mnt/md0/Ua/cases/ANT/ANT_Inverse; ANT_UaWrapper("","$pgid",\"Inverse\")";    
fi

if [ "$host" = "sauron" ]; then
    start="addpath(\"/home/wchm8/Documents/Ua/cases/ANT\"); cd /home/wchm8/Documents/Ua/cases/ANT/ANT_Inverse; ANT_UaWrapper("","$pgid",\"Inverse\")";
fi

$(

/usr/local/MATLAB/R2023b/bin/matlab -nodisplay -r "$start; quit;";

)&
