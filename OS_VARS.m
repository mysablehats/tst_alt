function [SLASH, pathtodata] = OS_VARS()
if ispc
    SLASH = '\\'; % windows 
    pathtodata = 'D:\fall_detection_datasets\TST Fall detection database ver. 2\';
elseif ismac
    pathtodata = '/Volumes/Elements/fall_detection_datasets/TST Fall detection database ver. 2/';
    SLASH = '/'; %
elseif isunix
    %pathtodata = '/media/Elements/fall_detection_datasets/TST Fall detection database ver. 2/';
    pathtodata = '/media/fbklein/Elements/fall_detection_datasets/TST Fall detection database ver. 2/';
    SLASH = '/'; %
else
    error('Cant determine OS')
end
if ~exist(pathtodata,'dir')    %%%this should be generalized so that it makes sense...
    pathtodata = 'E:\fall_detection_datasets\TST Fall detection database ver. 2\';
    if ~exist(pathtodata,'dir')
        pathtodata = 'G:\fall_detection_datasets\TST Fall detection database ver. 2\';
        if ~exist(pathtodata,'dir')
            error('data not there!')
        end
    end
end
end