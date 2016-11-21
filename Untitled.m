clear all;
close all;
load('onlineclass.mat'); 
labellabel = online_classifier(outstruct,allskel3, arc_conn, simvar);