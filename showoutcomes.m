function showoutcomes(simvar,data, numskells)

if isfield(simvar,'labels_names')&&~isempty(simvar.labels_names)
    %size(simvar.labels_names)
    for i = 1:length(data.gas)
        aoutout = [];
        if numskells>size(data.gas(i).class,2)
            claas = sum(data.gas(i).class(:,1:end).*(1:size(data.gas(i).class,1)).');
        else
            claas = sum(data.gas(i).class(:,1:numskells).*(1:size(data.gas(i).class,1)).');
        end
        [a, b ]= hist(claas, unique(claas));
        [~, ii] = sort(a,'descend');
        for j = ii
            if a(j) >0
                aoutout = ([ aoutout num2str(a(j)) ' ' simvar.labels_names{b(j)} '-   ']);
            end
        end,
        disp(aoutout),
        clear aoutout
    end
else
    disp('something went wrong')
    disp('no label names defined!')
end
end