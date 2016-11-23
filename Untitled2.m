clear all;
close all;
env = aa_environment;
classfdata = loadfileload('realclassifier',env); 
realvideo(classfdata.outstruct, baq(classfdata.pallconn), classfdata.simvar,0); 