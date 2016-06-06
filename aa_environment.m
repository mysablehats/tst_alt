function env = aa_environment()
global logpath
if ismac
    env.wheretosavestuff = '/Volumes/Seagate';
    env.homepath = '~/matlabprogs/';
    %disp('reached ismac')
elseif isunix
    env.wheretosavestuff = '/media/fbklein/Elements/fall_detection_datasets/var';
    env.homepath = '/home/fbklein/Documents/classifier/';
    %disp('reached isunix')
else
    %addpath(genpath(homepath))
    
    disp('oh-oh')
end
%open dialog box?
%have to see how to do it
logpath = strcat(env.homepath,'tst/var/log.txt');
[env.SLASH, env.pathtodata] = OS_VARS();
end