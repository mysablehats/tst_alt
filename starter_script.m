function simvar = starter_script(varargin)
%%%% STARTING MESSAGES PART FOR THIS RUN
global VERBOSE LOGIT TEST
VERBOSE = true;
LOGIT = false;
TEST = false; % set to false to actually run it
dbgmsg('=======================================================================================================================================================================================================================================')
dbgmsg('Running starter script')
dbgmsg('=======================================================================================================================================================================================================================================')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Each trial is trained on freshly partitioned/ generated data, so that we
% have an unbiased understanding of how the chained-gas is classifying.
%
% They are generated in a way that you can use nnstart to classify them and
% evaluated how much better (or worse) a neural network or some other
% algorithm can separate these datasets. Also, the data for each action
% example has different length, so the partition of datapoints is not
% equitative (there will be some fluctuation in the performance of putting
% every case in one single bin) and it will not be the same in validation
% and training sets. So in case this is annoying to you and you want to run
% always with a similar dataset, set
% simvar.generatenewdataset = false


env = aa_environment; % load environment variables

%creates a structure with the results of different trials
env.cstfilename=strcat(env.wheretosavestuff,env.SLASH,'cst.mat');
if exist(env.cstfilename,'file')
    loadst = load(env.cstfilename,'simvar');
    try
        simvar = loadst.simvar;
    catch
        simvar = struct();
    end
end
if ~exist('simvar','var')
    simvar = struct();
else
    simvar(end+1).nodes = [];%cst(1);
end

%% Choose dataset
if isempty(varargin)
    simvar(end).featuresall = 1;
    simvar(end).realtimeclassifier = false;
    simvar(end).generatenewdataset = false;
    simvar(end).datasettype = 'CAD60'; % datasettypes are 'CAD60', 'tstv2' and 'stickman'
    simvar(end).sampling_type = 'type1';
    simvar(end).activity_type = 'act_type'; %'act_type' or 'act'
    simvar(end).prefilter = 'none'; % 'filter', 'none', 'median?'
    simvar(end).labels_names = []; % necessary so that same actions keep their order number
    simvar(end).TrainSubjectIndexes = 'loo';%[9,10,11,4,8,5,3,6]; %% comment these out to have random new samples
    simvar(end).ValSubjectIndexes = [];%[1,2,7];%% comment these out to have random new samples
    simvar(end).randSubjEachIteration = true;
    simvar(end).extract = {'rand', 'wantvelocity'};
    simvar(end).preconditions = {'nohips', 'normal', 'norotatehips','mirrorx'};% {'nohips', 'norotatehips','mirrorx'}; %
    simvar(end).trialdataname = strcat('skel',simvar(end).datasettype,'_',simvar(end).sampling_type,simvar(end).activity_type,'_',simvar(end).prefilter, [simvar(end).extract{:}],[simvar(end).preconditions{:}]);
    simvar(end).trialdatafile = strcat(env.wheretosavestuff,env.SLASH,simvar(end).trialdataname,'.mat');
else
    simvar(end).featuresall = 3;%size(varargin{1},2);
    simvar(end).generatenewdataset = false;
    simvar(end).datasettype = 'Ext!';
    simvar(end).sampling_type = '';
    simvar(end).activity_type = ''; %'act_type' or 'act'
    simvar(end).prefilter = 'none'; % 'filter', 'none', 'median?'
    simvar(end).labels_names = []; % necessary so that same actions keep their order number
    simvar(end).TrainSubjectIndexes = [];%[9,10,11,4,8,5,3,6]; %% comment these out to have random new samples
    simvar(end).ValSubjectIndexes = [];%[1,2,7];%% comment these out to have random new samples
    simvar(end).randSubjEachIteration = true;
    simvar(end).extract = {''};
    simvar(end).preconditions = {''};
    simvar(end).trialdataname = strcat('other',simvar(end).datasettype,'_',simvar(end).sampling_type,simvar(end).activity_type,'_',simvar(end).prefilter, [simvar(end).extract{:}],[simvar(end).preconditions{:}]);
    simvar(end).trialdatafile = strcat(env.wheretosavestuff,env.SLASH,simvar(end).trialdataname,'.mat');
end

%% Setting up runtime variables

% set other additional simulation variables
simvar(end).TEST = TEST; %change this in the beginning of the program
simvar(end).PARA = 0;
simvar(end).P = 1;
simvar(end).NODES_VECT = 1000;
simvar(end).MAX_EPOCHS_VECT = [20];
simvar(end).ARCH_VECT = [1];
simvar(end).MAX_NUM_TRIALS = 1;
simvar(end).MAX_RUNNING_TIME = 1;%3600*10; %%% in seconds, will stop after this

% set parameters for gas:

params.layertype = '';
params.MAX_EPOCHS = [];
params.removepoints = true;
params.PLOTIT = true;
params.RANDOMSTART = true; % if true it overrides the .startingpoint variable
params.RANDOMSET = true; % if true, each sample (either alone or sliding window concatenated sample) will be presented to the gas at random
params.savegas.resume = false; % do not set to true. not working
params.savegas.save = false;
params.savegas.path = env.wheretosavestuff;
params.savegas.parallelgases = true;
params.savegas.parallelgasescount = 0;
params.savegas.accurate_track_epochs = true;
params.savegas.P = simvar(end).P;
params.startingpoint = [1 2];
params.amax = 50; %greatest allowed age
params.nodes = []; %maximum number of nodes/neurons in the gas
params.en = 0.006; %epsilon subscript n
params.eb = 0.2; %epsilon subscript b
params.gamma = 1; % for the denoising function
params.plottingstep = 0; % zero will make it plot only the end-gas

%Exclusive for gwr
params.STATIC = true;
params.at = 0.95; %activity threshold
params.h0 = 1;
params.ab = 0.95;
params.an = 0.95;
params.tb = 3.33;
params.tn = 3.33;

%Exclusive for gng
params.age_inc                  = 1;
params.lambda                   = 3;
params.alpha                    = .5;     % q and f units error reduction constant.
params.d                           = .99;   % Error reduction factor.


%% Begin loop
for architectures = simvar(end).ARCH_VECT
    for NODES = simvar(end).NODES_VECT
        for MAX_EPOCHS = simvar(end).MAX_EPOCHS_VECT
            for featuress = 1:simvar(end).featuresall
                if NODES ==100000 && (simvar(end).MAX_EPOCHS==1||simvar(end).MAX_EPOCHS==1)
                    dbgmsg('Did this already',1)
                    break
                end
                simvar(end).arch = architectures;
                simvar(end).NODES =  NODES;
                simvar(end).MAX_EPOCHS = MAX_EPOCHS;
                
                params.MAX_EPOCHS = simvar(end).MAX_EPOCHS;
                params.nodes = simvar(end).NODES; %maximum number of nodes/neurons in the gas
                
                 %% Classifier structure definitions
                
                simvar(end).allconn = allconnset(simvar(end).arch, params);
                
                
                %% Loading data
                data =  makess(length(baq(simvar(end).allconn))); % this breaks exectution core changes... 
                if isempty(varargin)
                    datasetmissing = false;
                    if ~exist(simvar(end).trialdatafile, 'file')&&~simvar(end).generatenewdataset
                        dbgmsg('There is no data on the specified location. Will generate new dataset.',1)
                        datasetmissing = true;
                    end
                    if simvar(end).generatenewdataset||datasetmissing
                        [allskel1, allskel2, simvar(end).TrainSubjectIndexes, simvar(end).ValSubjectIndexes] = generate_skel_data(simvar(end).datasettype, simvar(end).sampling_type, simvar(end).TrainSubjectIndexes, simvar(end).ValSubjectIndexes, simvar(end).randSubjEachIteration);
                        [allskel1, allskel2] = conformactions(allskel1,allskel2, simvar(end).prefilter); %%% should be split in 2
                        [data.train, simvar(end).labels_names] = extractdata(allskel1, simvar(end).activity_type, simvar(end).labels_names,simvar(end).extract{:});
                        [data.train, params.skelldef] = conformskel(data.train, simvar(end).preconditions{:});
                        %does same for validation data
                        [data.val, simvar(end).labels_names] = extractdata(allskel2, simvar(end).activity_type, simvar(end).labels_names,simvar(end).extract{:});
                        [data.val, ~                ] = conformskel(data.val,   simvar(end).preconditions{:}); 

                        simvar(end).trialdatafile = savefilesave(simvar(end).trialdataname, {data, simvar,params},env);
                        %save(simvar(end).trialdataname,'data', 'simvar','params');
                        dbgmsg('Training and Validation data saved.')
                        clear datasetmissing
                    else
                        loadedtrial = loadfileload(simvar(end).trialdataname,env);
                        data = loadedtrial.data;
                        params.skelldef = loadedtrial.params.skelldef;
                        simvar(end).generatenewdataset = false;
                    end
                    simvar(end).datainputvectorsize = size(data.train.data,1);
                else
                    data = varargin{1};
                    data = data(featuress);                    
                    simvar(end).datainputvectorsize = size(data.inputs,2);
                    params.skelldef = struct('length', simvar(end).datainputvectorsize, 'notskeleton', true, 'awk', struct('pos', [],'vel',[]), 'pos', simvar(end).datainputvectorsize, 'vel', []);
                    data.train.data = data.inputs'; % not empty so that the algorithm doesnt complain
                    data.train.y = data.labelsM;
                    data.train.ends = ones(1,size(data.inputs,1));
                    data.val.data = data.inputs';
                    data.val.y = data.labelsM;
                    data.val.ends = ones(1,size(data.inputs,1));                    
                end
                %% Classifier structure definitions
                
                simvar(end).allconn = allconnset(simvar(end).arch, params);
            
  %%%%%does this look like good programming?           
                
                
                %% Setting up different parameters for each of parallel tria
                % Maybe you want to do that; in case you don't, then we just
                % use a for to put the same parameters for each.
                a = struct([]);
                for i = 1:simvar(end).P
                    simvar(end).paramsZ(i) = params;
                    a(i).a = struct([]);
                end
                
                clear a
                b = [];
                
                if ~TEST % there are so many different ways I want to test it, that this definition is pretty much pointless.
                    starttime = tic;
                    while toc(starttime)< simvar(end).MAX_RUNNING_TIME
                        if length(b)> simvar(end).MAX_NUM_TRIALS
                            break
                        end
                        if simvar(end).PARA
                            spmd(simvar(end).P)
                                a(labindex).a = executioncore_in_starterscript(simvar(end),labindex, data);
                            end
                            %b = cat(2,b,a.a);
                            for i=1:length(a)
                                c = a{i};
                                a{i} = [];
                                b = [c.a b];
                            end
                            clear a c
                            a(1:simvar(end).P) = struct();
                        else
                            for i = 1:simvar(end).P
                                a(i).a = executioncore_in_starterscript(simvar(end),i, data);
                            end
                            b = cat(2,b,a.a);
                            clear a
                            a(1:simvar(end).P) = struct();
                        end
                    end
                else
                    b = executioncore_in_starterscript(simvar(end),1, data);
                end
                
                simvar(end).metrics = gen_cst(b); %%% it takes the important stuff from b;;; hopefully
                if isempty(varargin)
                    save(strcat(env.wheretosavestuff,env.SLASH,'cst.mat'),'simvar')
                    
                    savevar = strcat('b',num2str(simvar(end).NODES),'_', num2str(params.MAX_EPOCHS),'epochs',num2str(size(b,2)), simvar(end).sampling_type, simvar(end).datasettype, simvar(end).activity_type);
                    eval(strcat(savevar,'=simvar(end);'))
                    simvar(end).savesave = savefilesave(savevar, simvar(end),env);
                    %             simvar(end).savesave = strcat(env.wheretosavestuff,env.SLASH,savevar,'.mat');
                    %             ver = 1;
                    %
                    %             while exist(simvar(end).savesave,'file')
                    %                 simvar(end).savesave = strcat(env.wheretosavestuff,env.SLASH,savevar,'[ver(',num2str(ver),')].mat');
                    %                 ver = ver+1;
                    %             end
                    %             save(simvar(end).savesave,savevar)
                    dbgmsg('Trial saved in: ',simvar(end).savesave,1)
                end
                simvar(end+1) = simvar(end);
            end
            clear b
            clock
        end
    end
end
simvar(end) = [];
disp('hello')
end
function allconn = allconnset(n, params)
allconn_set = {...
    {... %%%% ARCHITECTURE 1
    {'gwr1layer',   'gwr',{'pos'},                    'pos',[1 0],params}...
    {'gwr2layer',   'gwr',{'vel'},                    'vel',[1 0],params}...
    {'gwr3layer',   'gwr',{'gwr1layer'},              'pos',[3 2],params}...
    {'gwr4layer',   'gwr',{'gwr2layer'},              'vel',[3 2],params}...
    {'gwrSTSlayer', 'gwr',{'gwr3layer','gwr4layer'},  'all',[3 2],params}...
    }...
    {...%%%% ARCHITECTURE 2
    {'gng1layer',   'gng',{'pos'},                    'pos',[1 0],params}...
    {'gng2layer',   'gng',{'vel'},                    'vel',[1 0],params}...
    {'gng3layer',   'gng',{'gng1layer'},              'pos',[3 2],params}...
    {'gng4layer',   'gng',{'gng2layer'},              'vel',[3 2],params}...
    {'gngSTSlayer', 'gng',{'gng4layer','gng3layer'},  'all',[3 2],params}...
    }...
    {...%%%% ARCHITECTURE 3
    {'gng1layer',   'gng',{'pos'},                    'pos',[1 0],params}...
    {'gng2layer',   'gng',{'vel'},                    'vel',[1 0],params}...
    {'gng3layer',   'gng',{'gng1layer'},              'pos',[3 0],params}...
    {'gng4layer',   'gng',{'gng2layer'},              'vel',[3 0],params}...
    {'gngSTSlayer', 'gng',{'gng4layer','gng3layer'},  'all',[3 0],params}...
    }...
    {...%%%% ARCHITECTURE 4
    {'gwr1layer',   'gwr',{'pos'},                    'pos',[1 0],params}...
    {'gwr2layer',   'gwr',{'vel'},                    'vel',[1 0],params}...
    {'gwr3layer',   'gwr',{'gwr1layer'},              'pos',[3 0],params}...
    {'gwr4layer',   'gwr',{'gwr2layer'},              'vel',[3 0],params}...
    {'gwrSTSlayer', 'gwr',{'gwr3layer','gwr4layer'},  'all',[3 0],params}...
    }...
    {...%%%% ARCHITECTURE 5
    {'gwr1layer',   'gwr',{'pos'},                    'pos',[1 2 3],params}...
    {'gwr2layer',   'gwr',{'vel'},                    'vel',[1 2 3],params}...
    {'gwr3layer',   'gwr',{'gwr1layer'},              'pos',[3 2],params}...
    {'gwr4layer',   'gwr',{'gwr2layer'},              'vel',[3 2],params}...
    {'gwrSTSlayer', 'gwr',{'gwr3layer','gwr4layer'},  'all',[3 2],params}...
    }...
    {...%%%% ARCHITECTURE 6
    {'gwr1layer',   'gwr',{'pos'},                    'pos',[3 4 2],params}...
    {'gwr2layer',   'gwr',{'vel'},                    'vel',[3 4 2],params}...
    {'gwrSTSlayer', 'gwr',{'gwr1layer','gwr2layer'},  'all',[3 2],params}...
    }...
    {...%%%% ARCHITECTURE 7
    {'gwr1layer',   'gwr',{'all'},                    'all',[3 2], params}...
    {'gwr2layer',   'gwr',{'gwr1layer'},              'all',[3 2], params}...
    }...
    {...%%%% ARCHITECTURE 8
    {'gwr1layer',   'gwr',{'pos'},                    'pos',[1 0], params}... %% now there is a vector where q used to be, because we have the p overlap variable...
    }...
    {...%%%% ARCHITECTURE 9
    {'gwr1layer',   'gwr',{'pos'},                    'pos',3,params}...
    {'gwr2layer',   'gwr',{'vel'},                    'vel',3,params}...
    {'gwr3layer',   'gwr',{'gwr1layer'},              'pos',3,params}...
    {'gwr4layer',   'gwr',{'gwr2layer'},              'vel',3,params}...
    {'gwr5layer',   'gwr',{'gwr3layer'},              'pos',3,params}...
    {'gwr6layer',   'gwr',{'gwr4layer'},              'vel',3,params}...
    {'gwrSTSlayer', 'gwr',{'gwr6layer','gwr5layer'},  'all',3,params}
    }...
    {... %%%% ARCHITECTURE 10
    {'gwr1layer',   'gwr',{'pos'},                    'pos',[1 0],params}...
    {'gwr2layer',   'gwr',{'vel'},                    'vel',[1 0],params}...
    {'gwrSTSlayer', 'gwr',{'gwr1layer','gwr2layer'},  'all',[3 2],params}...
    }...
    };
allconn = allconn_set{n};
end
function a = executioncore_in_starterscript(simvar, i, data)
paramsZ = simvar.paramsZ(i);
allconn = simvar.allconn;

global TEST
n = randperm(size(data.train.data,2)-3,2); % -(q-1) necessary because concatenation reduces the data size!
paramsZ.startingpoint = [n(1) n(2)];
pallconn = allconn;
pallconn{1,1}{1,6} = paramsZ; % I only change the initial points of the position gas
%pallconn{1,3}{1,6} = paramsZ; %but I want the concatenation to reflect the same position that I randomized. actually this is not going to happen because of the sliding window scheme
%pallconn{1,4}{1,6} = pallconn{1,2}{1,6};

%[a.sv, a.mt] = starter_sc(data, pallconn, 1);
if TEST
    dbgmsg('TEST RUN. Generating sham output data. Data will not be saved.',1)
    confconf = struct('val','val', 'train', '');
    ouout = struct('accumulatedepochs',0);
    for i =1:4
        for j =1:5
            a.mt(i,j) = struct('conffig', 'hello','confusions', confconf,'conffvig', 'hello','outparams',ouout);
        end
    end
    
else
    [outstruct, a.mt, ~] = starter_sc(data, pallconn);
    %%disp('hello')
    %%chunk = makechunk(data);
    %%load('chunk.mat');
    save('realclassifier.mat', 'outstruct', 'pallconn', 'simvar')
 %   [~, something_to_classify] = realvideo(outstruct, baq(pallconn), simvar,0);   
    % realvideo(outstruct, baq(pallconn), simvar);

    %%online_classifier(chunkchunk, outstruct, baq(pallconn), simvar);
end
end

