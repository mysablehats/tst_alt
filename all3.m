function [data,labels_names, skelldef] = all3(allskel, simvar)
[allskel] = conformactions(allskel, simvar.prefilter);
[data, labels_names] = extractdata(allskel, simvar.activity_type, simvar.labels_names,simvar.extract{:});
[data, skelldef] = conformskel(data, simvar.preconditions{:});
end