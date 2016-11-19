function data = online_classifier(data, allskel, arq_connect, simvar)
global VERBOSE LOGIT
VERBOSE = 1; % true;
LOGIT = false;
%a = struct();

%runs conformactions

[allskel,~] = conformactions(allskel, [], simvar.prefilter); %%% should be split in 2

%extractdata
[data.test, ~] = extractdata(allskel, simvar.activity_type, simvar.labels_names,simvar.extract{:}); %%% I dont know if I need this, or if it makes sense to use this function
%conformskel
[data.test, ~ ] = conformskel(data.test, simvar.preconditions{:}, 'test', simvar.paramsZ.skelldef);
%loop from l.65 of starter_sc -> make it into a function???
for j = 1:length(arq_connect)
    [~, data.test] = gas_method(data, data.test, 'vot', arq_connect(j),j,size(data.train.data,1));
end


for i = 1:length(arq_connect)
    %labeller
    %for i gases... like above
    [data.test.gas(i).class, ~] = find(data.gas(i).nodesl(:, data.test.gas(i).bestmatchbyindex));
end
if isfield(simvar,'labels_names')&&~isempty(simvar.labels_names)
    %size(simvar.labels_names)
    for i = 1:length(arq_connect)
        aoutout = [];
        [a, b ]= hist(data.test.gas(i).class, unique(data.test.gas(i).class));
        [~, ii] = sort(a,'descend');
        for j = ii
            aoutout = ([ aoutout num2str(a(j)) ' ' simvar.labels_names{b(j)} '-   ']);
        end, 
        disp(aoutout), 
        clear aoutout 
    end
else
    disp('something went wrong')
    disp('no label names defined!')
end
%disp('')
if simvar.paramsZ.PLOTIT
    figure
    plot_act(data, data.test.data, simvar.paramsZ.skelldef) %% data.test!!!!
end

end