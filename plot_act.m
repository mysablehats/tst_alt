function plot_act(datagas,data,i, skelldef, numskels)
if numskels > size(data.gas(i).bestmatchbyindex,2)
    numskels = size(data.gas(i).bestmatchbyindex,2);
end

%%% this function is strange. I am making too

%%% so many functions!!!! who can keep track of all of them without
%%% documentation?
%%%
%%% Certainly not me
%%%
%%% This function will plot an action sequence recorded from the sensor in
%%% matlab environment and then plot the bestmatching action sequence it
%%% found from the gas.
%%%
%%% This function should be called from the scope of online_classifier and
%%% it needs this kind of data.
%numskels = 3;
%maybe plot with subplot and show also bestmatching
subplot(3,1,1)
dataplotskel = skelreshape(data.gas(i).inputs.input(:,1:numskels),skelldef);

dataplotdata = skeldraw(dataplotskel,'ts',skelldef);%,'rt', skelldef); %%to plot this i need the skelldef :-(
% I will just plot the last one...

%%%listen this function should be more generic and accept parameters, the
%%%whole 405 thing is disgusting. it is a 45*9 samples long
subplot(3,1,2)
%i = 5;
%%%% here there is a funny fact. we only have 9 skeletons long so we need
%%%% to check if we have an overflow
gasplotskel = skelreshape(datagas(i).nodes(:,data.gas(i).bestmatchbyindex(1:numskels)),skelldef);
gasplotdata = skeldraw(gasplotskel,'ts',skelldef);

subplot(3,1,3)
plot3(dataplotdata(1,:),dataplotdata(2,:),dataplotdata(3,:))
hold
plot3(gasplotdata(1,:),gasplotdata(2,:),gasplotdata(3,:))
hold off