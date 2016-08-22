function simvar = starter_script(varargin)
%%%% STARTING MESSAGES PART FOR THIS RUN
dbgmsg('=======================================================================================================================================================================================================================================')
dbgmsg('Running starter script')
dbgmsg('=======================================================================================================================================================================================================================================')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global VERBOSE LOGIT TEST
VERBOSE = true;
LOGIT = true;
TEST = false; % set to false to actually run it

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
    simvar = loadst.simvar;
end
if ~exist('simvar','var')
    simvar = struct();
else
    simvar(end+1).nodes = [];%cst(1);
end

%% Choose dataset
if isempty(varargin)
    simvar(end).featuresall = 1;
    simvar(end).generatenewdataset = false;
    simvar(end).datasettype = 'tstv2'; % datasettypes are 'CAD60', 'tstv2' and 'stickman'
    simvar(end).sampling_type = 'type2';
    simvar(end).activity_type = 'act'; %'act_type' or 'act'
    simvar(end).prefilter = 'none'; % 'filter', 'none', 'median?'
    simvar(end).labels_names = []; % necessary so that same actions keep their order number
    simvar(end).TrainSubjectIndexes = [];%[9,10,11,4,8,5,3,6]; %% comment these out to have random new samples
    simvar(end).ValSubjectIndexes = [];%[1,2,7];%% comment these out to have random new samples
    simvar(end).randSubjEachIteration = true;
    simvar(end).extract = {'rand', 'wantvelocity'};
    simvar(end).preconditions = {'nohips', 'normal', 'norotatehips','mirrorx'};
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
simvar(end).TEST = TEST;
simvar(end).PARA = 1;
simvar(end).P = 4;
simvar(end).NODES_VECT = 1000;
simvar(end).MAX_EPOCHS_VECT = [1];
simvar(end).ARCH_VECT = [8];
simvar(end).MAX_NUM_TRIALS = 1;
simvar(end).MAX_RUNNING_TIME = 1;%3600*10; %%% in seconds, will stop after this

% set parameters for gas:

params.MAX_EPOCHS = [];
params.removepoints = true;
params.PLOTIT = false;
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
                
                %% Loading data
                if isempty(varargin)
                    datasetmissing = false;
                    if ~exist(simvar(end).trialdatafile, 'file')&&~simvar(end).generatenewdataset
                        dbgmsg('There is no data on the specified location. Will generate new dataset.',1)
                        datasetmissing = true;
                    end
                    if simvar(end).generatenewdataset||datasetmissing
                        [allskel1, allskel2, simvar(end).TrainSubjectIndexes, simvar(end).ValSubjectIndexes] = generate_skel_data(simvar(end).datasettype, simvar(end).sampling_type, simvar(end).TrainSubjectIndexes, simvar(end).ValSubjectIndexes, simvar(end).randSubjEachIteration);
                        [allskel1, allskel2] = conformactions(allskel1,allskel2, simvar(end).prefilter);
                        [data.train, simvar(end).labels_names] = extractdata(allskel1, simvar(end).activity_type, simvar(end).labels_names,simvar(end).extract{:});
                        [data.val, simvar(end).labels_names] = extractdata(allskel2, simvar(end).activity_type, simvar(end).labels_names,simvar(end).extract{:});
                        [data, params.skelldef] = conformskel(data, simvar(end).preconditions{:});
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
                                a(labindex).a = executioncore_in_starterscript(simvar(end).paramsZ(labindex),simvar(end).allconn, data);
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
                                a(i).a = executioncore_in_starterscript(simvar(end).paramsZ(i),simvar(end).allconn, data);
                            end
                            b = cat(2,b,a.a);
                            clear a
                            a(1:simvar(end).P) = struct();
                        end
                    end
                else
                    b = executioncore_in_starterscript(simvar(end).paramsZ(1),simvar(end).allconn, data);
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
end
function metrics = gen_cst(b)
metrics = struct;
metrics(length(b),length(b(1).mt)) = struct; %,2,size(b(1).mt(1).confusions.val,1),size(b(1).mt(1).confusions.val,2));
for ii = 1:length(b)
    for jj= 1:length(b(ii).mt)
        metrics(ii,jj).conffig = b(ii).mt(jj).conffig;
        metrics(ii,jj).val = b(ii).mt(jj).confusions.val;
        metrics(ii,jj).train = b(ii).mt(jj).confusions.train;
        if ~isfield(b(ii).mt(jj).outparams, 'accumulatedepochs')
            metrics(ii,jj).accumulatedepochs = paramsZ(ii).MAX_EPOCHS;
        else
            metrics(ii,jj).accumulatedepochs = b(ii).mt(jj).outparams.accumulatedepochs;
        end
    end
end
% function cst = gen_cst(b)
% cst = struct([]);
% cst(end).metrics(length(b),length(b(1).mt)) = struct; %,2,size(b(1).mt(1).confusions.val,1),size(b(1).mt(1).confusions.val,2));
% for ii = 1:length(b)
%     for jj= 1:length(b(ii).mt)
%         cst(end).metrics(ii,jj).conffig = b(ii).mt(jj).conffig;
%         cst(end).metrics(ii,jj).val = b(ii).mt(jj).confusions.val;
%         cst(end).metrics(ii,jj).train = b(ii).mt(jj).confusions.train;
%         if ~isfield(b(ii).mt(jj).outparams, 'accumulatedepochs')
%             cst(end).metrics(ii,jj).accumulatedepochs = paramsZ(ii).MAX_EPOCHS;
%         else
%             cst(end).metrics(ii,jj).accumulatedepochs = b(ii).mt(jj).outparams.accumulatedepochs;
%         end
%     end
% end
end
function savesave = savefilesave(filename, savevar,env)
global TEST
ver = 1;
savesave = strcat(env.wheretosavestuff,env.SLASH,filename,'.mat');
while exist(savesave,'file')
    savesave = strcat(env.wheretosavestuff,env.SLASH,filename,'[ver(',num2str(ver),')].mat');
    ver = ver+1;
end
if ~TEST
    if iscell(savevar)&&(length(savevar)==3) % hack! I know the only thing I will save that has 3 cells is the dataset
        data = savevar{1}; %#ok<*NASGU>
        simvar = savevar{2};
        params = savevar{3};
        save(savesave, 'data', 'simvar','params')
    else
        save(savesave,'savevar')
    end
    dbgmsg('Saved file:',savesave,1)
end
end
function loadload = loadfileload(filename,env)
ver = 0;
loadfile = strcat(env.wheretosavestuff,env.SLASH,filename,'.mat');
while exist(loadfile,'file')
    ver = ver+1;
    loadfile = strcat(env.wheretosavestuff,env.SLASH,filename,'[ver(',num2str(ver),')].mat');
end
if ver == 1
    loadfile = strcat(env.wheretosavestuff,env.SLASH,filename,'.mat');
else
    ver = ver - 1;
    loadfile = strcat(env.wheretosavestuff,env.SLASH,filename,'[ver(',num2str(ver),')].mat');
end
loadload = load(loadfile);
dbgmsg('Loaded file:',loadfile,1)
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
function a = executioncore_in_starterscript(paramsZ,allconn, data)
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
    [~, a.mt] = starter_sc(data, pallconn);
end
end
function [savestructure, metrics] = starter_sc(savestructure, allconn)
%% starter_sc
% This is the main function to run the chained classifier, label and
% generate confusion matrices and recall, precision and F1 values for the
% skeleton classifier of activities using an architecture of chained neural
% gases on skeleton activities data (the STS V2 Dataset). This work is an
% attempt to implement Parisi, 2015's paper.

%%
global VERBOSE LOGIT
VERBOSE = true;
LOGIT = true;
dbgmsg('ENTERING MAIN LOOP')
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% making metrics structure

metrics = struct('confusions',[],'conffig',[],'outparams',[]);
%%%% building arq_connect
arq_connect = struct;
%arq_connect(1:length(allconn)) = struct('name','','method','','sourcelayer','', 'layertype','','q',[1 0],'params',struct([]));
for i = 1:length(allconn)
    arq_connect(i).name = allconn{i}{1};
    arq_connect(i).method = allconn{i}{2};
    arq_connect(i).sourcelayer = allconn{i}{3};
    arq_connect(i).layertype = allconn{i}{4};
    arq_connect(i).q = allconn{i}{5};
    arq_connect(i).params = allconn{i}{6};
    %%% hack I need this info in params as well
    arq_connect(i).params.q = arq_connect(i).q;
end
inputs = struct('input_clip',[],'input',[],'input_ends',[],'oldwhotokill',struct([]), 'index', []);
gas_data = struct('name','','class',[],'y',[],'inputs',inputs,'bestmatch',[],'bestmatchbyindex',[],'whotokill',struct([]));
gas_methodss = struct;
gas_methodss = struct('name','','edges',[],'nodes',[],'fig',[],'nodesl',[]); %bestmatch will have the training matrix for subsequent layers
for i =1:length(arq_connect)
    gas_methods(i) = gas_methodss;
end
savestructure.gas = gas_methods;
savestructure.train.indexes = [];
savestructure.train.gas = gas_data;
savestructure.val.indexes = [];
savestructure.val.gas = gas_data;

for i = 1:length(savestructure) % oh, I don't know how to do it elegantly
    savestructure.figset = {};
end

%%% end of gas structures region


%% Gas-chain Classifier
% This part executes the chain of interlinked gases. Each iteration is one
% gas, and currently it works as follows:
% 1. Function setinput() chooses based on the input defined in allconn
% 2. Run either a Growing When Required (gwr) or Growing Neural Gas (GNG)
% on this data
% 3. Generate matrix of best matching units to be used by the next gas
% architecture

dbgmsg('Starting chain structure for GWR and GNG for nodes:',num2str(labindex),1)
dbgmsg('###Using multilayer GWR and GNG ###',1)

for j = 1:length(arq_connect)
    [savestructure, savestructure.train] = gas_method(savestructure, savestructure.train,'train', arq_connect(j),j, size(savestructure.train.data,1)); % I had to separate it to debug it.
    metrics(j).outparams = savestructure.gas(j).outparams;
    [savestructure, savestructure.val ]= gas_method(savestructure, savestructure.val,'val', arq_connect(j),j, size(savestructure.train.data,1));
end


%% Gas Outcomes
% This should be made into a function... It is a nice graph to perhaps
% debug the gas...
figure
for j = 1:length(arq_connect)
    if arq_connect(j).params.PLOTIT
        subplot (1,length(arq_connect),j)
        hist(savestructure.gas(j).outparams.graph.errorvect)
        title((savestructure.gas(j).name))
    end
end
%% Labelling
% The current labelling procedure for both the validation and training datasets. As of this moment I label all the gases
% to see how adding each part increases the overall performance of the
% structure, but since this is slow, the variable whatIlabel can be changed
% to contain only the last gas.
%
% The labelling procedure is simple. It basically picks the label of the
% closest point and assigns to that. In a sense the gas can be seen as a
% dimensional (as opposed to temporal) filter, encoding prototypical
% action-lets.
dbgmsg('Labelling',num2str(labindex),1)

whatIlabel = 1:length(savestructure.gas); %change this series for only the last value to label only the last gas

%%
% Specific part on what I want to label
for j = whatIlabel
    dbgmsg('Applying labels for gas: ''',savestructure.gas(j).name,''' (', num2str(j),') for process:',num2str(i),1)
    [savestructure.train.gas(j).class, savestructure.val.gas(j).class,savestructure.gas(j).nodesl ] = labeller(savestructure.gas(j).nodes, savestructure.train.gas(j).bestmatchbyindex,  savestructure.val.gas(j).bestmatchbyindex, savestructure.train.gas(j).inputs.input, savestructure.train.gas(j).y);
    %%%% I dont understand what I did, so I will code this again.
    %%% why did I write .bestmatch when it should be nodes??? what was I thinnking [savestructure.train.gas(j).class, savestructure.val.gas(j).class] = untitled6(savestructure.gas(j).bestmatch, savestructure.train.gas(j).inputs.input,savestructure.val.gas(j).inputs.input, y_train);
    
end

%% Displaying multiple confusion matrices for GWR and GNG for nodes
% This part creates the matrices that can later be shown with the
% plotconfusion() function.

savestructure.figset = {}; %% you should clear the set first if you want to rebuild them
dbgmsg('Displaying multiple confusion matrices for GWR and GNG for nodes:',num2str(labindex),1)

for j = whatIlabel
    [~,metrics(j).confusions.val,~,~] = confusion(savestructure.val.gas(j).y,savestructure.val.gas(j).class);
    [~,metrics(j).confusions.train,~,~] = confusion(savestructure.train.gas(j).y,savestructure.train.gas(j).class);
    
    dbgmsg(savestructure.gas(j).name,' Confusion matrix on this validation set:',writedownmatrix(metrics(j).confusions.val),1)
    savestructure.gas(j).fig.val =   {savestructure.val.gas(j).y,     savestructure.val.gas(j).class,  strcat(savestructure.gas(j).name,savestructure.gas(j).method,'V')};
    savestructure.gas(j).fig.train = {savestructure.train.gas(j).y,   savestructure.train.gas(j).class,strcat(savestructure.gas(j).name,savestructure.gas(j).method,'T')};
    %savestructure.figset = [savestructure.figset, savestructure.gas(j).fig.val, savestructure.gas(j).fig.train];
    %%%
    metrics(j).conffig = savestructure.gas(j).fig;
end

%% Actual display of the confusion matrices:
metitems = [];
for j = whatIlabel
    if arq_connect(j).params.PLOTIT
        metitems = [metitems j*arq_connect(j).params.PLOTIT];
    end
end
if ~isempty(metitems)
    figure
    plotconf(metrics(metitems))
end
%plotconf(savestructure.figset{:})
figure
plotconfusion(savestructure.gas(end).fig.val{:})

end
function [sst, sstv] = gas_method(sst, sstv, vot, arq_connect,j, dimdim)
%% Gas Method
% This is a function to go over a gas of the classifier, populate it with the apropriate input and generate the best matching units for the next layer.
%% Setting up some labels
sst.gas(j).name = arq_connect.name;
sst.gas(j).method = arq_connect.method;
sst.gas(j).layertype = arq_connect.layertype;
arq_connect.params.layertype = arq_connect.layertype;

%% Choosing the right input for this layer
% This calls the function set input that chooses what will be written on the .inputs variable. It also handles the sliding window concatenations and saves the .input_ends properties, so that this can be done recursevely.
% After some consideration, I have decided that all of the long inputing
% will be done inside setinput, because it it would be easier.
dbgmsg('Working on gas: ''',sst.gas(j).name,''' (', num2str(j),') with method: ',sst.gas(j).method ,' for process:',num2str(labindex),1)

[sstv.gas(j).inputs.input_clip, sstv.gas(j).inputs.input, sstv.gas(j).inputs.input_ends, sstv.gas(j).y, sstv.gas(j).inputs.oldwhotokill, sstv.gas(j).inputs.index, sstv.gas(j).inputs.awk ]  = setinput(arq_connect, sst, dimdim, sstv); %%%%%%

%%
% After setting the input, we can actually run the gas, either a GNG or the
% GWR function we wrote.
if strcmp(vot, 'train')
    %DO GNG OR GWR
    [sst.gas(j).nodes, sst.gas(j).edges, sst.gas(j).outparams] = gas_wrapper(sstv.gas(j).inputs.input_clip,arq_connect);
end
%%%% POS-MESSAGE
dbgmsg('Finished working on gas: ''',sst.gas(j).name,''' (', num2str(j),') with method: ',sst.gas(j).method ,'.Num of nodes reached:',num2str(sst.gas(j).outparams.graph.nodesvect(end)),' for process:',num2str(labindex),1)

%% Best-matching units
% The last part is actually finding the best matching units for the gas.
% This is a simple procedure where we just find from the gas units (nodes
% or vectors, as you wish to call them), which one is more like our input.
% It is a filter of sorts, and the bestmatch matrix is highly repetitive.

% I questioned if I actually need to compute this matrix here or maybe
% inside the setinput function. But I think this doesnt really matter.
% Well, for the last gas it does make a difference, since these units will
% not be used... Still I will  not fix it unless I have to.
%PRE MESSAGE
dbgmsg('Finding best matching units for gas: ''',sst.gas(j).name,''' (', num2str(j),') for process:',num2str(labindex),1)
[~, sstv.gas(j).bestmatchbyindex] = genbestmmatrix(sst.gas(j).nodes, sstv.gas(j).inputs.input, arq_connect.layertype, arq_connect.q); %assuming the best matching node always comes from initial dataset!

%% Post-conditioning function
%This will be the noise removing function. I want this to be optional or allow other things to be done to the data and I
%am still thinking about how to do it. Right now I will just create the
%whattokill property and let setinput deal with it.
if arq_connect.params.removepoints
    dbgmsg('Flagging noisy input for removal from gas: ''',sst.gas(j).name,''' (', num2str(j),') with points with more than',num2str(arq_connect.params.gamma),' standard deviations, for process:',num2str(labindex),1)
    sstv.gas(j).whotokill = removenoise(sst.gas(j).nodes, sstv.gas(j).inputs.input, sstv.gas(j).inputs.oldwhotokill, arq_connect.params.gamma, sstv.gas(j).inputs.index);
else
    dbgmsg('Skipping removal of noisy input for gas:',sst.gas(j).name)
end
end
function allskel = LoadDataBase(idx_folderi)

%%%%%%Messages part. Provides feedback for the user about what is being
%%%%%%done
dbgmsg('Loading Database and building allskel structure. Please wait, this takes a while.')
%%%%%%


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Choose:
%startFoldIdx = 1;   % Index of the first test folder to load (1->'Data1')
%stopFoldIdx = 5;    % Index of the last test folder to load (11->'Data11')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%The script loads:
% - PckVal_*:       packet from wearable device
% - TStamp_QPC_*:   timestamp assigned by the PC to the packet, when it arrives to the PC [us]
% - TStamp_Int_*:   internal timestamp of the wearable device [ms]
% - *_acc_raw_*:    acceleration data
% - M:              depth frame
% - KinectTime:     timestamps assigned by the PC to the frame, when it is available to the PC [100ns]/[us]
% - jMatDep:        skeleton joints in depth space
% - jMatSkl:        skeleton joints in skeleton space
% - KinectTimeBody: timestamps assigned by the PC to the skeleton, when it is available to the PC [100ns]/[us]

[SLASH, pathtodata] = OS_VARS(); % adjust environment variables to system in question (PC, MAC, UNIX)

ADLFolderName = {'sit','grasp','walk','lay'};
FallFolderName = {'front','back','side','EndUpSit'};
%*****Wearable device*****
device1 = '35EE'; %Number of device1
device2 = '36F9'; %Number of device2

%*******Kinect - V2*******
rowPixel = 424; %[pixel] number of row
columnPixel = 512;  %[pixel] number of column
frameIdx = 1; %depth frame number

% Load selected folders
for idx_folder = idx_folderi
    for groupName = {'ADL' 'Fall'}
        if strcmp(groupName,'ADL')
            subfolder = ADLFolderName;
        else
            subfolder = FallFolderName;
        end
        for name_Subfolder = subfolder
            for idx_test = 1:3
                Folder = strcat(pathtodata, 'Data',num2str(idx_folder),SLASH,cell2mat(groupName),SLASH,cell2mat(name_Subfolder),SLASH,num2str(idx_test)); %Folder where are stored the data
                %*************************
                %Load wearable device data
                %*************************
                %**35EE**
                PckVal_35EE = loadPackets(device1,Folder);
                PckValSize_35EE = size(PckVal_35EE);
                TStamp_QPC_35EE = csvread(strcat(Folder,SLASH,'Time',SLASH,'TimeStamps',device1,'.csv'));
                TStamp_Int_35EE = (PckVal_35EE(:,3)+256.*PckVal_35EE(:,4)+256.*256.*PckVal_35EE(:,5));
                %accelerometer values
                X_acc_raw_35EE = PckVal_35EE(:,8)+256*PckVal_35EE(:,9);
                Y_acc_raw_35EE = PckVal_35EE(:,10)+256*PckVal_35EE(:,11);
                Z_acc_raw_35EE = PckVal_35EE(:,12)+256*PckVal_35EE(:,13);
                %36F9
                PckVal_36F9 = loadPackets(device2,Folder);
                PckValSize_36F9 = size(PckVal_36F9);
                TStamp_QPC_36F9 = csvread(strcat(Folder,SLASH,'Time',SLASH,'TimeStamps',device2,'.csv'));
                TStamp_Int_36F9 = (PckVal_36F9(:,3)+256.*PckVal_36F9(:,4)+256.*256.*PckVal_36F9(:,5));
                %accelerometer values
                X_acc_raw_36F9 = PckVal_36F9(:,8)+256*PckVal_36F9(:,9);
                Y_acc_raw_36F9 = PckVal_36F9(:,10)+256*PckVal_36F9(:,11);
                Z_acc_raw_36F9 = PckVal_36F9(:,12)+256*PckVal_36F9(:,13);
                
                %*************************
                %****Load depth frame*****
                %*************************
                fid = fopen(strcat(Folder,SLASH,'Depth',SLASH,'Filedepth_',num2str(frameIdx-1),'.bin'),'r');
                arrayFrame = fread(fid,'uint16');
                fclose(fid);
                M = zeros(rowPixel, columnPixel);
                for r=1:rowPixel
                    M(r,:) = arrayFrame((r-1)*columnPixel+1:r*columnPixel);
                end
                %Load time information
                KinectTime = csvread(strcat(Folder,SLASH,'Time',SLASH,'DepthTime.csv'));
                
                %*************************
                %******Load skeleton******
                %*************************
                fileNameSk1DS = strcat(Folder,SLASH,'Body',SLASH,'Fileskeleton.csv'); %joint in the depth frame
                fileNameSk1SS = strcat(Folder,SLASH,'Body',SLASH,'FileskeletonSkSpace.csv'); %joint in 3D space
                Sk1SkDepth = csvread(fileNameSk1DS);
                Sk1SkSpace = csvread(fileNameSk1SS);
                %Find player number
                for idx_player = 1:6
                    if Sk1SkDepth(25*(idx_player-1)+1,1) ~= 0
                        break;
                    end
                end
                %Find row index of the specific player in Sk1SkDepth
                row_idx = find(Sk1SkDepth(:,5) == idx_player-1);
                Sk1SkDepth = Sk1SkDepth(row_idx,:);
                NumFrameSkelDepth = fix(length(Sk1SkDepth(:,1))/25);
                Sk1SkSpace = Sk1SkSpace(row_idx,:);
                NumFrameSkelSpace = fix(length(Sk1SkSpace(:,1))/25);
                %restore 25 joints in groups
                jMatDep = zeros(25, 3, NumFrameSkelDepth);
                jMatSkl = zeros(25, 3, NumFrameSkelSpace);
                for n = 1:NumFrameSkelDepth
                    jMatDep(:,:,n) = Sk1SkDepth(((n-1)*25+1):n*25,1:3);
                end
                for n = 1:NumFrameSkelSpace
                    jMatSkl(:,:,n) = Sk1SkSpace(((n-1)*25+1):n*25,1:3);
                end
                %Load time information
                KinectTimeBody = csvread(strcat(Folder,SLASH,'Time',SLASH,'BodyTime.csv'));
                
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % Put here your code!
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                deltatimes = diff(KinectTimeBody(:,1));
                meany = mean(deltatimes);
                steady = std(deltatimes);
                if any(((deltatimes)-meany)>.8*meany)||steady>6e3||any((deltatimes-meany).^2/steady^2>3)||meany<2e5||meany>6e5
                    figure
                    hist(deltatimes)
                    %error('You will have to handle sample rate!!!')
                end
                %%% Calculate velocities!!
                vel = zeros(size(jMatSkl)); % initial velocity is zero
                ncnc = meany*30; %normalizing constant since the sensor is hopefully always at 30 fps
                for i = 2:size(jMatSkl,3)
                    vel(:,:,i) = ncnc*(jMatSkl(:,:,i) - jMatSkl(:,:,i-1))/(KinectTimeBody(i,1)-KinectTimeBody(i-1,1));
                end
                jskelstruc = struct('skel',jMatSkl, 'act',groupName,'act_type', name_Subfolder, 'index', idx_test, 'subject', idx_folder,'time',KinectTimeBody,'vel',vel);
                %%%%%% size(jskelstruc.skel) % this was here for debugging
                if exist('allskel','var') % initialize my big matrix of skeletons
                    allskel = cat(2,allskel,jskelstruc);
                else
                    allskel = jskelstruc;
                end
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % Plot raw Acceleration data of each test
                %figure;
                %subplot 121; title('35EE'); hold on;
                %plot(X_acc_raw_35EE,'b'); plot(Y_acc_raw_35EE,'r'); plot(Z_acc_raw_35EE,'k');
                %subplot 122; title('36F9'); hold on;
                %plot(X_acc_raw_36F9,'b'); plot(Y_acc_raw_36F9,'r'); plot(Z_acc_raw_36F9,'k');
                
                % Plot Depth and Skeleton data of each test
                %figure;
                %subplot 131;
                %imagesc(M); title('depth frame');
                %subplot 132;  hold on;
                %plot3(jMatSkl(:,1,1),jMatSkl(:,3,1),jMatSkl(:,2,1),'.r','markersize',20); view(0,0); axis equal;
                %title('skeleton in skeleton space');
                %subplot 133;
                %plot3(jMatDep(:,1,1),jMatDep(:,3,1),jMatDep(:,2,1),'.b','markersize',20); view(0,0); axis equal; set(gca,'ZDir','Reverse');
                %title('skeleton in depth space');
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            end
        end
    end
end

end
function PckVal = loadPackets(device,path)
[SLASH, ~] = OS_VARS();
n = 22; % No. of columns of T
fid = fopen(strcat(path,SLASH,'Shimmer',SLASH,'Packets',device,'.bin'),'r');
B = fread(fid,'uint8');
fclose(fid);
BB = reshape(B, n,[]);
PckVal = permute(BB,[2,1]);
end
function [allskel1, allskel2] = conformactions(allskel1,allskel2, whichfilter)
%%%disp('hello')
%whichfilter = 'median';

switch whichfilter
    case 'filter'
        %%% this filter is likely bad because it introduces phase shift!!
        windowSize = 5;
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
        filterfun = @(x)filter(b,a,x);
    case 'median'
        medianmedian = 5;
        filterfun = @(x)medfilt1(x,medianmedian);
    case 'none'
        return
    otherwise
        error('Unknown filter.')
end

for i = 1:size(allskel1,2)
    for j = 1:size(allskel1(i).vel,1)
        for k = 1:size(allskel1(i).vel,2)
            allskel1(i).vel(j,k,:) = filterfun(allskel1(i).vel(j,k,:));
        end
    end
end

for i = 1:size(allskel2,2)
    for j = 1:size(allskel2(i).vel,1)
        for k = 1:size(allskel2(i).vel,2)
            allskel2(i).vel(j,k,:) = filterfun(allskel2(i).vel(j,k,:));
        end
    end
end


% create the matrix from the structure to access the database and run the
% classification....
end
function [data, lab] = extractdata(structure, typetype, inputlabels,varargin)
for i = 1:length(varargin)
    switch varargin{i}
        case 'wantvelocity'
            WANTVELOCITY = true;
        case 'rand'
            RANDSQE = true;
        case 'novelocity'
            WANTVELOCITY = true;
        case 'seq'
            RANDSQE = false;
        otherwise
            error('unexpected argument')
    end
end

%%%%%%%%Messages part. Feedback for the user about the algorithm
dbgmsg('Extracting data from skeleton structure',1)
if WANTVELOCITY
    dbgmsg('Constructing long vectors with velocity data as well',1)
end
%%%%%%%%
%typetype= 'act_type';
%typetype= 'act';

Data = [];
ends = [];
% approach
if strcmp(typetype,'act')
    [labelZ,~] = alllabels(structure,inputlabels);
elseif strcmp(typetype,'act_type')
    [~, labelZ] = alllabels(structure,inputlabels);
else
    error('weird typetype!')
end
lab = sort(labelZ);

Y = [];

if RANDSQE
    randseq = randperm(length(structure));
else
    randseq = 1:length(structure);
end

for i = randseq % I think each iteration is one action
    Data = cat(3, Data, structure(i).skel);
    Y = cat(2, Y, repmat(whichlab(structure(i),lab,typetype),1,size(structure(i).skel,3)));
    ends = cat(2, ends, size(structure(i).skel,3));
end
if WANTVELOCITY
    Data_vel = [];
    for i = randseq % I think each iteration is one action
        Data_vel = cat(3, Data_vel, structure(i).vel);
    end
    Data = cat(1,Data, Data_vel);
end
% It will also construct data for a clustering analysis, whatever the hell
% that might mean in this sense
vectordata = [];
for i = 1:length(Data)
    vectordata = cat(2,vectordata, [Data(:,1,i); Data(:,2,i); Data(:,3,i)]);
end
%Y = Y';
%Data, vectordata, Y, ends, lab]

data.data = vectordata;
%data.tensor = Data;
data.y = Y;
data.ends = ends;

end
function [lab, biglab] = alllabels(st,lab)
%lab = cell(0);
biglab = lab;

if isfield(st,'act')&&isfield(st,'act_type')
    for i = 1:length(st) % I think each iteration is one action
        cu = strfind(lab, st(i).act);
        if isempty(lab)||isempty(cell2mat(cu))
            lab = [{st(i).act}, lab];
            biglab = [{[st(i).act st(i).act_type]}, biglab];
        end
        bgilab = [st(i).act st(i).act_type];
        cu = strfind(biglab, bgilab);
        if isempty(biglab)||isempty(cell2mat(cu))
            biglab = [{bgilab}, biglab];
        end
    end
end
if isfield(st,'act_type')
    for i = 1:length(st) % I think each iteration is one action
        cu = strfind(biglab, st(i).act_type);
        if isempty(biglab)||isempty(cell2mat(cu))
            biglab = [{st(i).act_type}, biglab];
        end
        
    end
elseif isfield(st,'act')
    for i = 1:length(st) % I think each iteration is one action
        cu = strfind(lab, st(i).act);
        if isempty(lab)||isempty(cell2mat(cu))
            lab = [{st(i).act}, lab];
        end
    end
    
else
    error('No action fields in data structure.')
end

end
function outlab = whichlab(st,lb,tt)
numoflabels = size(lb,2);
switch tt
    case 'act_type'
        if isfield(st, 'act')
            comp_act = [st.act st.act_type];
        else
            comp_act = st.act_type;
        end
    case 'act'
        comp_act = st.act;
    otherwise
        error('Unknown classification type!')
end

for i = 1:numoflabels
    if strcmp(lb{i},comp_act)
        lab = i;%i-1;
    end
end


%I thought lab was a good choice, but matlab
outlab = zeros(numoflabels,1);
outlab(lab) = 1;
end
function allskel = readcad60(subject)

if length(subject)>1
    allskel = [];
    for i = subject
        allskel = [allskel readcad60(i)];
    end
    return
end

[slash, pathtodata] = OS_VARS;
%%% read
realpath = strcat(pathtodata, '..',slash,'CAD60',slash,'data', num2str(subject), slash);
ftoopen = strcat(realpath, 'activityLabel.txt');
fidfid = fopen(ftoopen, 'r');
Axr = textscan(fidfid, '%f,%s','Delimiter', '\n', 'CollectOutput', true);
fclose(fidfid);

% allskel = Axr

for i=1:length(Axr{1})
    activitytoopen = strcat(realpath, num2str(Axr{1}(i)), '.txt');
    fidfid = fopen(activitytoopen, 'r');
    if fidfid == -1
        activitytoopen = strcat(realpath, '0', num2str(Axr{1}(i)), '.txt');
        fidfid = fopen(activitytoopen, 'r');
        if fidfid == -1
            error(strcat('could not open file',activitytoopen))
        end
    end
    ori_def =  strcat(repmat('%f,',1,9),'%u8,');
    j_def = strcat(repmat('%f,',1,3),'%u8,');
    massofformlessdata = textscan(fidfid, strcat('%u,',repmat(strcat(ori_def, j_def),1,11),repmat(j_def,1,4)),'Delimiter', '\n', 'CollectOutput', true);
    fclose(fidfid);
    actlen = size(massofformlessdata{1},1);
    skel = nan(15,3,actlen); %initalize it with nans, because why not?
    vel = zeros(15,3,actlen);
    skelori = nan(11,9,actlen);
    jskelstruc = struct('skel',skel,'skelori', skelori,'act_type', Axr{2}{i}(1:end-1), 'index', Axr{1}(i), 'subject', subject,'time',[],'vel',vel);
    for j = 1:actlen
        for k = 1:11
            skel(k,:,j) = massofformlessdata{4*k}(j,:);
            skelori(k,:,j) = massofformlessdata{4*k-2}(j,:);
        end
        for k = 12:15
            skel(k,:,j) = massofformlessdata{22+2*k}(j,:);
        end
        if j>1
            for k = 1:15
                vel(k,:,j) = skel(k,:,j)-skel(k,:,j-1);
            end
        end
        
    end
    jskelstruc.skel = skel;
    jskelstruc.skelori = skelori;
    jskelstruc.vel = vel;
    % size(jskelstruc.skel) % this was here for debugging
    if exist('allskel','var') % initialize my big matrix of skeletons
        allskel = cat(2,allskel,jskelstruc);
    else
        allskel = jskelstruc;
    end
    
    
end
end
function [allskel1, allskel2, allskeli1, allskeli2] = generate_skel_data(varargin)
% generates the .mat files that have the generated datasets for training and
% validation
% based on 2 types of random sampling
%
%%%% type 1 data sampling: known subjects, unknown individual activities from them:
%
%%%% type 2 data sampling: unknown subjects, unknown individual activities from them:
%
%
% data_train, y_train
% data_val, y_val
%

%aa_environment
if nargin<2
    error('I need at least the name of the data set (either ''CAD60'' or ''tstv2'') and the data sampling type (either ''type1'' or ''type2'')')
end

dataset = varargin{1};
sampling_type = varargin{2};
if nargin>2
    allskeli1 = varargin{3};
end
if nargin>3
    allskeli2 = varargin{4};
end

if isempty(allskeli1)||isempty(allskeli2)||varargin{5}
    %%% maybe it would make sense here to just generate the complementary
    %%% group with setdiff, but I decided for: if any subject definition is
    %%% empty, or simvar(end).randSubjEachIteration is true, then generate
    %%% new sample.
    clear allskeli1
    clear allskeli2
end

switch dataset
    case 'CAD60'
        loadfun = @readcad60;
        datasize = 4;
    case 'tstv2'
        loadfun = @LoadDataBase;
        datasize = 11;
    case 'stickman'
        loadfun = @generate_falling_stick;
        datasize = 20;
        if strcmp(sampling_type,'type1')
            dbgmsg('Sampling type1 not implemented for falling_stick!!! Using type2',1)
            sampling_type = 'type2';
        end
    otherwise
        error('Unknown database.')
end

%%% checks to see if indices will be within array range
if exist('allskeli1','var')
    if any(allskeli1>datasize)
        error('Index 1 is out of range for selected dataset.')
    end
end
if exist('allskeli2','var')
    if any(allskeli2>datasize)
        error('Index 2 is out of range for selected dataset.')
    end
end


%%%%%%%%Messages part: provides feedback for the user
dbgmsg('Generating random datasets for training and validation')
if exist('allskeli1','var')
    dbgmsg('Variable allskeli1 is defined. Will skip randomization.')
end
%%%%%%%%%%%%

%%%% type 1 data sampling: known subjects, unknown individual activities from them:
if strcmp(sampling_type,'type1')
    allskel = loadfun(1:datasize); %main data
    if ~exist('allskeli1','var')
        allskeli1 = randperm(length(allskel),fix(length(allskel)*.8)); % generates the indexes for sampling the dataset
    end
    allskel1 = allskel(allskeli1);
    if ~exist('allskeli2','var')
        allskeli2 = setdiff(1:length(allskel),allskeli1); % use the remaining data as validation set
    end
    allskel2 = allskel(allskeli2);
end

%%%% type 2 data sampling: unknown subjects, unknown individual activities from them:
if strcmp(sampling_type,'type2')
    if ~exist('allskeli1','var')
        allskeli1 = randperm(datasize,fix(datasize*.8)); % generates the indexes for sampling the dataset
    end
    allskel1 = loadfun(allskeli1(1)); %initializes the training dataset
    for i=2:length(allskeli1)
        allskel1 = cat(2,loadfun(allskeli1(i)),allskel1 );
    end
    
    allskeli2 = setdiff(1:datasize,allskeli1); % use the remaining data as validation set
    
    allskel2 = loadfun(allskeli2(1)); %initializes the training dataset
    for i=2:length(allskeli2)
        allskel2 = cat(2,loadfun(allskeli2(i)),allskel2 );
    end
end
%%%%%%
% saves data
%%%%%%

end
function allskel = generate_falling_stick(numsticks)

% if length(numsticks)>1
%     allskel = [];
%     for i = subject
%         allskel = [allskel generate_falling_stick(i)];
%     end
%     return
% end
stickmodel = 'skel25'; % 'rod', 'skel25', 'skel15'
thetanoisepower = 0.01;
translationnoisepower = 1;
floorsize = 1500;

for kk = 1:numsticks
    num_of_points_in_stick = 25;
    
    %%% random size
    %%% random start location
    %%% random initial velocity for falling stick
    
    l = (1.6+rand()*.3)*1000;
    
    for akk = [0,1]
        
        startlocation = floorsize*[rand(), 0.01*rand(), rand()];
        phi = 2*pi*rand();
        phinoisepower = 0.1*rand;
        initial_velocity = -1*rand();
        
        if akk
            act = 'Fall';
            % from http://www.chem.mtu.edu/~tbco/cm416/MatlabTutorialPart2.pdf and
            % from http://ocw.mit.edu/courses/mechanical-engineering/2-003j-dynamics-and-control-i-spring-2007/lecture-notes/lec10.pdf
            % the equation from the falling stick is
            %
            % -m*g*l/2*cos(t) = (Ic+m*l^2/4*cos(t)^2)*tdd - m*l^2/4*cos(t)*sin(t*td^2)
            %
            % tdd = diff(td)
            % td = diff(t)
            
            testOptions = odeset('RelTol',1e-3,'AbsTol', [1e-4; 1e-2]);
            notsatisfiedwithmodel = true;
            while(notsatisfiedwithmodel)
                [t,x] =     ode45(@stickfall, [0 2+1.5*rand()], [pi/2;initial_velocity], testOptions);
                
                
                
                %%% resample to kinect average sample rate, i.e. 15 or 30 hz
                x = resample(x(:,1),t,30);
                
                %%% upon visual inpection we see that the beginnings and ends
                %%% of the results from ode45 sequence are very noisy, possibly
                %%% due to the unnatural constraints on the differential
                %%% equation mode. Some doctoring is required:
                
                x(1,:) = x(2,:)+0.01*rand(); %%% fixing initial point
                %             %now to fix the ends it is even dirtier:
                %             plot(x(:,1))
                %             hold on
                diffx = abs(diff(x));
                for i = length(diffx):-1:3
                    if all([diffx(i-2), diffx(i-1), diffx(i)]<0.07) %this was done on trial and error basis
                        x = x(1:i);
                        if abs(x(end))<0.1||abs(x(end)-pi)<0.1 %%it has to end in a natural resting angle
                            notsatisfiedwithmodel = false;
                        end
                        break
                    end
                end
            end
            padding = x(end).*ones(round(130*rand()),1);
            padding = padding+thetanoisepower*rand(size(padding));
            x = [x;padding];
            %             plot(x(:,1))
            %             hold off
            phinoise = zeros(size(x)); %%% a freefalling stick would never change angle, ok maybe if it were spinning, or hit by something...
            translationnoise = translationnoisepower*zeros(size(x,1),3);
            [skel, vel] = construct_skel(x,l,num_of_points_in_stick ,startlocation, translationnoise, phi, phinoise, stickmodel);
            
        else
            act = 'Walk';
            %%% non falling activity
            x = pi/2*ones(90+round(25*rand()),1)+0*thetanoisepower; % does it even make sense if the velocity changes and the position does not?
            phinoise = 0*phinoisepower*rand(size(x));
            %%% the walk now
            t = linspace(0,1,length(x));
            realtransx = 2*sin(phi)*t*floorsize;
            realtransy = 0*t*floorsize;
            realtransz = 2*cos(phi)*t*floorsize;
            translationnoise = 0*translationnoisepower*rand(size(x,1),3)+[realtransx; realtransy; realtransz]';
            [skel, vel] = construct_skel(x, l, num_of_points_in_stick, startlocation, translationnoise , phi, phinoise, stickmodel);
            
        end
        %for i = 1:s
        construct_sk_struct = struct('x',x,'l',l,'num_points', num_of_points_in_stick,'startlocation',startlocation,'phi',phi);
        jskelstruc = struct('skel',skel,'act_type', act, 'index', kk, 'subject', kk,'time',[],'vel', vel, 'construct_sk_struct', construct_sk_struct);
        
        %plot(stickstick(:,1), stickstick(:,2), '*')
        if exist('allskel','var') % initialize my big matrix of skeletons
            allskel = cat(2,allskel,jskelstruc);
        else
            allskel = jskelstruc;
        end
    end
end
end
function dx = stickfall(t,x)

m = 1;
g = 9.8;
l = 1.6;
t = 0;

Ic = 1/3*m*l^2;%% for a slender rod rotating on one end

x1 = x(1);
x2 = x(2);

% state equations
% I can put any nonlinearilty I want, so I will make so that when  x1 == 0
% or x1 = pi, then the velocity changes sign, i.e. it bounces

if ((x1 < 0)&&x2<0)||((x1> pi)&&x2>0)
    dx1 = -0.5*x2; %so that it dampens as well
else
    dx1 = x2;
end

dx2 = (m*l^2/4*cos(x1*x2^2) -m*g*l/2*cos(x1))/(Ic+m*l^2/4*cos(x1)^2);

dx = [dx1;dx2];

end
function [stickstick,stickvel] = construct_skel(thetha, l, num_of_points_in_stick, displacement, tn, phi, phinoise, stickmodel)

bn = [rand(), rand(), rand()]/10; %% I need a bit of noise or the classifier gets insane

simdim = size(thetha,1); % 4;%

stickstick = zeros(num_of_points_in_stick,3,simdim);
stickvel = stickstick; %%%
switch stickmodel
    case 'rod'
        for i=1:simdim
            dpdp =  displacement+tn(i,:);
            tt=thetha(i);
            PHI = phi+phinoise(i);
            
            stickstick(1,:,i) = dpdp-l/2*cos(tt)*[cos(PHI) 0 sin(PHI)];
            for j = 2:num_of_points_in_stick
                stickstick(j,:,i) = ([cos(tt)*cos(PHI), sin(tt), cos(tt)*sin(PHI)]+bn)*l*j/num_of_points_in_stick+dpdp-l/2*cos(tt)*[cos(PHI) 0 sin(PHI)];
            end
        end
    case 'skel25'
        num_of_points_in_stick = 25;
        %%%% skeleton model was obtained after using generate_skel_data
        %protoskel = allskel1(3).skel(:,:,1) -repmat(mean(allskel1(3).skel(:,:,1)),25,1);
        prot = load('protoskel.mat');
        height = 1552;
        if ~isfield(prot,'protoskel')
            error('Problems to load protoskel.mat: skeleton model not found.')
        end
        for i=1:simdim
            dpdp =  displacement+tn(i,:);
            tt=thetha(i);
            PHI = phi+phinoise(i);
            %%% will need to make 2 rotation matrices
            % PHI is first, it is a rotation around the y axis
            %             PHIMAT = [[cos(PHI) 0 sin(PHI)];[0 1 0];[-sin(PHI) 0 cos(PHI)]];
            %             % now, is a fall a rotation on x or z axis?
            %             % I will do both and comment out the one I dislike; or maybe I
            %             % can make it a different kind of fall.
            %             %around x
            %             if 1
            %                 ttmatx = [[1 0 0 ];[0 cos(tt) -sin(tt)];[0 sin(tt) cos(tt)]];
            %                 ttmatz = eye(3);
            %             else
            %                 %around z
            %                 ttmatx = eye(3);
            %                 ttmatz = [[cos(tt) -sin(tt) 0];[sin(tt) cos(tt) 0 ];[0 0 1]];
            %             end
            %stickstick(1,:,i) = dpdp;
            askel = rotskel(prot.protoskel2+repmat([0 height/2 0],num_of_points_in_stick,1),0,PHI,pi/2-tt);
            for j = 1:num_of_points_in_stick
                %                 stickstick(j,:,i) = (ttmatz*PHIMAT*ttmatx*prot.protoskel2(j,:).').'/height*l+dpdp-l/2*cos(tt)*[cos(PHI) 0 sin(PHI)]+[0 height/2 0];
                stickstick(j,:,i) = askel(j,:)/height*l+dpdp-l/2*cos(tt)*[cos(PHI) 0 sin(PHI)];
                
            end
            %there is something wrong with my matrix multiplication I think,   but I have no time to be elegant about it
            %or probably means I am wrong about the kinect data's axis; in any case, this would only affect conformskel mirror functions,
            %I should plot database(posei), hold and plot
            %mirrordatabase(posei) to see if the skeleton is upsidedown...
            %2/3 chances it isnt...
            %quick and dirty
            %stickstick(:,1,i) = - stickstick(:,1,i);
            
        end
end
%%% the initial velocities are not zero, but I had to artificially change
%%% the x point because of limitations of the simulated stick. They are
%%% already zeroed, so I just need to start from the second point

%%% for the next ones I will calculate the points
for i = 2:simdim
    stickvel(:,:,i-1) = stickstick(:,:,i)-stickstick(:,:,i-1);
end

%stickstick; % ok, I forgot that the reshape happens latter %reshape(stickstick,[],simdim);


end
function [extinput_clipped, extinput, inputends,y, removeremove, indexes, awko] = setinput(arq_connect, savestruc,data_size, svst_t_v) %needs to receive the correct data size so that generateidx may work well
%%%%%% this is the place to get long inputs actually.
%arqconnect has only the current layer, so it is flat
%inputends need to be the same for everything to work out fine
%theoretically it makes sense that they are not the same size and then the
%smallest should be used and the end bits of each action set, discarded and
%then rematched to fit the real correspondent. I don't really know how to
%do that, maybe I need to match each action with some extra indexing, or
%make a clever indexing function that will discard the right amount of end
%bits at the right part
% Or I perhaps should carry each action separatelly, because it would make
% things easier, but this would require a major rewrite of everything I did
% so far. So in short: same "ends" for every component.
extinput = [];
midremove = [];
awko = [];
inputinput = cell(length(arq_connect.sourcelayer),1);
removeremove = struct([]);
awk = inputinput; %it is also the same size...
%if inputends dont coincide this function will give a strange error
inputends = [];
[posidx, velidx] = generateidx(data_size, arq_connect.params.skelldef);
for j = 1:length(arq_connect.sourcelayer)
    foundmysource = false;
    for i = 1:length(savestruc.gas)
        if strcmp(arq_connect.sourcelayer{j}, savestruc.gas(i).name)
            if isempty( svst_t_v.gas(i).bestmatchbyindex)
                error('wrong computation order. bestmatch field not yet defined.')
            end
            oldinputends = inputends;
            [inputinput{j},inputends,y, indexes] = longinput( savestruc.gas(i).nodes(:,svst_t_v.gas(i).bestmatchbyindex), arq_connect.q, svst_t_v.gas(i).inputs.input_ends, svst_t_v.gas(i).y,svst_t_v.gas(i).inputs.index);
            
            %%%check for misalignments of inputends
            if ~isempty(oldinputends)
                if ~all(oldinputends==inputends)
                    error('Misaligned layers! Alignment not yet implemented.')
                end
            end
            %%% old  longinput call. I will no longer create .bestmatch, so
            %%% I need to create it on the fly from gasnodes
            %            [inputinput{j},inputends,y] = longinput( svst_t_v.gas(i).bestmatch, arq_connect.q, svst_t_v.gas(i).inputs.input_ends, svst_t_v.gas(i).y);
            
            %inputinput{j} = longinput(savestruc.gas(i).bestmatch; %
            removeremove = structcat(svst_t_v.gas(i).whotokill, removeremove);
            %%% this part will construct my newly created awk vector out of
            %%% initial awk vectors
            awk{j} = makeawk(arq_connect.q, svst_t_v.gas(i).inputs.awk);
            foundmysource = true;
        end
    end
    if ~foundmysource
        if strcmp(arq_connect.layertype, 'pos')
            [inputinput{j},inputends,y, indexes] = longinput(svst_t_v.data(posidx,:), arq_connect.q, svst_t_v.ends, svst_t_v.y, (1:size(svst_t_v.data,2)));
            %inputinput{j} = svst_t_v.data(posidx,:); %
            %ends is savestructure.train.ends
            awk{j} = makeawk(arq_connect.q, arq_connect.params.skelldef.awk.pos);
        elseif strcmp(arq_connect.layertype, 'vel')
            [inputinput{j},inputends,y, indexes] = longinput(svst_t_v.data(velidx,:), arq_connect.q, svst_t_v.ends, svst_t_v.y, (1:size(svst_t_v.data,2)));
            %inputinput{j} = svst_t_v.data(velidx,:); %
            %ends is savestructure.train.ends
            awk{j} = makeawk(arq_connect.q, arq_connect.params.skelldef.awk.vel);
        elseif strcmp(arq_connect.layertype, 'all')
            [inputinput{j},inputends,y, indexes] = longinput(svst_t_v.data, arq_connect.q, svst_t_v.ends, svst_t_v.y, (1:size(svst_t_v.data,2)));
            %inputinput{j} = svst_t_v.data; %
            %ends is savestructure.train.ends
            awk{j} = makeawk(arq_connect.q, [arq_connect.params.skelldef.awk.pos;  arq_connect.params.skelldef.awk.vel] );
        end
    end
    if isempty(inputinput)
        error(strcat('Unknown layer type:', arq_connect.layertype,'or sourcelayer:',arq_connect.sourcelayer))
    end
end
if length(inputinput)>1
    for i = 1:length(inputinput)
        extinput = cat(1,extinput,inputinput{i}); % this part should check for the right ends, ends should also be a cell array, and they should be concatenated properly
        %%% dealing with possible empty sets::
        if ~isempty(removeremove)
            midremove = structcat(midremove,removeremove); %{i}
            %midremove = cat(1,midremove,removeremove{i}); %{i}
        end
        awko = cat(1,awko,awk{i});
    end
else
    extinput = inputinput{:};
    midremove = removeremove; %no turtles!
    awko = awk{:};
    %oh, it may be wrong,,, have
    %to check wrappung and unwrapping %%it was wrong, I will check
    %unwrapping inside removebaddata
    %[extinput_clipped, inputends_clipped, y_clipped]= removebaddata(extinput, inputends, y, removeremove); % this wrong, since I will only avoid putting bad sets into the gas. so it is fortunately simpler than I thought!
    
end
extinput_clipped= removebaddata(extinput, indexes, midremove, arq_connect.q);
end
function awk = makeawk(q,inawk)
awk = repmat(inawk,q(1),1);
end
function [linput,newends, newy, indexes] = longinput(shortinput, qp, ends, y, iindex)
% this function was getting messy, so I decided to recreate the structure
% that generated her, so to make easier debugging
% It is very disellegant of me. I apologise.
q = qp(1);
switch length(qp)
    case 1
        p = 0;
        r = 1;
    case 2
        p = qp(2);
        r = 1;
    case 3
        p = qp(2);
        r = qp(3);
end

realends = cumsum(ends,2);
actionstructure = struct;
actionstructure(1:size(ends,2)) = struct();
actionstructure(1).pose = shortinput(:,1:realends(1));
actionstructure(1).end = ends(1);
actionstructure(1).y = y(:,realends(1));
actionstructure(1).index = iindex(1:ends(1));

for i = 2:size(ends,2)
    actionstructure(i).pose = shortinput(:,realends(i-1)+1:realends(i));
    actionstructure(i).end = ends(i);
    actionstructure(i).y = y(:,realends(i));
    actionstructure(i).index = iindex(realends(i-1):realends(i));
end
shortdim = size(shortinput,1);
for i = 1:length(actionstructure)
    m = 1;
    for j = 1:1+p:actionstructure(i).end
        a = zeros(shortdim*q,1);
        indexx = zeros(1,q);
        if j+q*r-1>actionstructure(i).end
            %cant complete the whole vector!
            break
        else
            k = 1;
            for lop = 1:q
                a(1+(k-1)*shortdim:k*shortdim) = actionstructure(i).pose(:,j+lop*r-1);
                indexx(lop) = actionstructure(i).index(j+lop*r-1);
                k = k+1;
            end
        end
        %have to save a somewhere
        actionstructure(i).long(m).vec = a;
        actionstructure(i).long(m).index = indexx;
        m = m+1;
    end
    %should concatenate long now
    actionstructure(i).newend = length(actionstructure(i).long);
    actionstructure(i).longinput = zeros(q*shortdim,actionstructure(i).newend);
    actionstructure(i).longy = repmat(actionstructure(i).y,1,actionstructure(i).newend);
    actionstructure(i).longindex = zeros(q,size(actionstructure(i).longy,2));
    for j = 1:actionstructure(i).newend
        actionstructure(i).longinput(:,j) = actionstructure(i).long(j).vec;
        actionstructure(i).longindex(:,j) = actionstructure(i).long(j).index;
    end
end
linput = actionstructure(1).longinput;
newy = actionstructure(1).longy;
newends = zeros(1,length(actionstructure));
newends(1) = actionstructure(1).newend;
indexes = actionstructure(1).longindex;
for i = 2:length(actionstructure)
    linput = cat(2,linput,actionstructure(i).longinput);
    newy = cat(2,newy, actionstructure(i).longy);
    newends(i) = actionstructure(i).newend;
    indexes = cat(2,indexes, actionstructure(i).longindex);
end
end
function stringmatrix = writedownmatrix(matrix)
stringmatrix = '';
for i =1:size(matrix,1)
    for j = 1:size(matrix,2)
        stringmatrix = strcat(stringmatrix,'(', num2str(i),',',num2str(j),'):', num2str(matrix(i,j)),'\t');
    end
end
end
function [ matmat, matmat_byindex] = genbestmmatrix(nodes, data, ~,~)
%matmat = zeros(size(nodes,1),size(data,2));
%matmat_byindex = zeros(1,size(data,2));
[~,matmat_byindex] = pdist2(nodes',data','euclidean','Smallest',1);
matmat = data(matmat_byindex);
%
% for i = 1:size(data,2)
%        [ matmat(:,i), matmat_byindex(i)] = bestmatchingunit(data(:,i),gwr_nodes,whichisit,q);
% end

% % old genbestmatch. I figured out it is the other way around, which makes
% more sense, so I should just filter spatially my data, so the
% bestmatching data should have the same dimension as the initial dataset,
% duh...
% matmat = zeros(size(data,1),size(gwr_nodes,2));
% for i = 1:size(gwr_nodes,2) %%%%%%%%%% it is actually the other way around!!!!!
%         matmat(:,i) = bestmatchingunit(gwr_nodes(:,i),data,whichisit);
% end
% end
end
function labels = labeling(nodes, data, y)
% the labeling function for the GWR and GNG
% cf parisi, adopt the label of the nearest node

%so first we need to know the label of each node
% we do exactly as parisi, with the labeling after the learning phase
% this he calls the training phase (but everything is the training phase)
% go through the nodeset and see in the data_set what is more similar
%%%%%MESSAGES
dbgmsg('Applying labels to the prototypical nodes.',1)
%%%%%
try
    [~,ni] = pdist2(data',nodes', 'euclidean', 'Smallest',1);
    labels = y(:,ni);
catch
    [labels, ~ ]= labelling(nodes, data, y);
end


end
function [labels, ni1 ]= labelling(nodes, data, y)
maxmax = size(nodes,2);
labels = zeros(1,maxmax);

for i = 1:maxmax
    [~, ~, ni1 , ~ , ~] = findnearest(nodes(:,i), data); % s1 is the index of the nearest point in data
    labels(i) = y(ni1);
end
end
function [train_class, val_class, nodesl] = labeller(nodes, train_bestmatchbyindex, val_bestmatchbyindex, train_input, y_train)

nodesl = labeling(nodes,train_input,y_train);

tcsize = size(train_bestmatchbyindex,2);
vcsize = size(val_bestmatchbyindex,2);
ycsize = size(y_train,1);

train_class = zeros(ycsize,tcsize);
val_class = zeros(ycsize,vcsize);

for i =1:tcsize
    train_class(:,i) = nodesl(:,train_bestmatchbyindex(i));
end
for i =1:vcsize
    val_class(:,i) = nodesl(:,val_bestmatchbyindex(i));
end
end
