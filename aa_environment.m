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
elseif ispc
    env.wheretosavestuff = 'd:\';
    env.homepath = 'C:\\Users\\Frederico\\Documents\\GitHub\\classifier';
else
    disp('oh-oh')
end
%addpath(genpath(homepath))
%open dialog box?
%have to see how to do it
[env.SLASH, env.pathtodata] = OS_VARS();
if ispc
    logpath = strcat(env.homepath, env.SLASH , env.SLASH ,'var',env.SLASH ,env.SLASH,'log.txt');
else
    logpath = strcat(env.homepath, env.SLASH ,'var',env.SLASH,'log.txt');
end
if ~exist(logpath,'file')
    %if this fails there may be a broken path,,, usually happens when I
    %change the computer I am using
    fid = fopen( logpath, 'wt' );
    fprintf( fid, 'New logfile created.');
    fclose(fid);
end
end