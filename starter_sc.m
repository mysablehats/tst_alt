function [savestructure, metrics, gases] = starter_sc(savestructure, allconn)
gases = [];
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
    if size(savestructure.val.data,1) > 0
        [savestructure, savestructure.val ]= gas_method(savestructure, savestructure.val,'val', arq_connect(j),j, size(savestructure.train.data,1));
    end
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
    [savestructure.train.gas(j).class, savestructure.val.gas(j).class, savestructure.gas(j).nodesl ] = labeller(savestructure.gas(j).nodes, savestructure.train.gas(j).bestmatchbyindex,  savestructure.val.gas(j).bestmatchbyindex, savestructure.train.gas(j).inputs.input, savestructure.train.gas(j).y);
    %%%% I dont understand what I did, so I will code this again.
    %%% why did I write .bestmatch when it should be nodes??? what was I thinnking [savestructure.train.gas(j).class, savestructure.val.gas(j).class] = untitled6(savestructure.gas(j).bestmatch, savestructure.train.gas(j).inputs.input,savestructure.val.gas(j).inputs.input, y_train);
    
end

%% Displaying multiple confusion matrices for GWR and GNG for nodes
% This part creates the matrices that can later be shown with the
% plotconfusion() function.

savestructure.figset = {}; %% you should clear the set first if you want to rebuild them
dbgmsg('Displaying multiple confusion matrices for GWR and GNG for nodes:',num2str(labindex),1)

for j = whatIlabel
    [~,metrics(j).confusions.train,~,~] = confusion(savestructure.train.gas(j).y,savestructure.train.gas(j).class);
    savestructure.gas(j).fig.train = {savestructure.train.gas(j).y,   savestructure.train.gas(j).class,strcat(savestructure.gas(j).name,savestructure.gas(j).method,'T')};
    
    if isfield(savestructure.val, 'data')&&~isempty(savestructure.val.data)%&&~isempty(savestructure.val.data)
        [~,metrics(j).confusions.val,~,~] = confusion(savestructure.val.gas(j).y,savestructure.val.gas(j).class);      
        
        dbgmsg(savestructure.gas(j).name,' Confusion matrix on this validation set:',writedownmatrix(metrics(j).confusions.val),1)
        savestructure.gas(j).fig.val =   {savestructure.val.gas(j).y,     savestructure.val.gas(j).class,  strcat(savestructure.gas(j).name,savestructure.gas(j).method,'V')};
        %savestructure.figset = [savestructure.figset, savestructure.gas(j).fig.val, savestructure.gas(j).fig.train];
        %%%
    else
        metrics(j).confusions.val = [];
        savestructure.gas(j).fig.val = [];
    end
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
if isfield(savestructure.gas(end).fig, 'val')&&~isempty(savestructure.gas(end).fig.val)
    figure
    plotconfusion(savestructure.gas(end).fig.val{:})
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