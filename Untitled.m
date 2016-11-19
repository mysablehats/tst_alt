clear all;
close all;
load('onlineclass.mat'); 
allskel3 = generate_skel_online(chunk.chunk);
labellabel = online_classifier(outstruct,allskel3, arc_conn, simvar);