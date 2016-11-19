function plot_act(data, other, skelldef)

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

%maybe plot with subplot and show also bestmatching
subplot(3,1,1)
skeldraw(other(:,1:9),'rt', skelldef); %%to plot this i need the skelldef :-(
% I will just plot the last one...

%%%listen this function should be more generic and accept parameters, the
%%%whole 405 thing is disgusting. it is a 45*9 samples long
subplot(3,1,2)
i = 5;
skeldraw(data.gas(i).nodes(1:405,data.test.gas(i).bestmatchbyindex),'mt',9);

subplot(3,1,3)
skeldraw(other(:,1:9),'rt', skelldef);
hold
skeldraw(data.gas(i).nodes(1:405,data.test.gas(i).bestmatchbyindex),'mt',9);