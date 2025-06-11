clear all
close all
clc

rgb_fn = 'C:\damon\data\VanOrley\A14340_F-OL.tif';
sc = 0.9592;
j01 = 3;
j02 = 4;
repcnt_max = 1;
max_shift = 12;
minpts = 64;
Npoly = 0;  
Ntran = 0;
memory_limited = 0;
morepts = 0;
useresults = 0; 
filteron = 0;
addsharp = 0;
debugmode = 0;
usegui = 0;

ir_fn1 = 'C:\damon\data\VanOrley\VanOrley_ChristwithDoctors_H.img';
ir_extract(ir_fn1,0,usegui)
rough_mosaic(ir_fn1,usegui)

ir_fn2 = 'C:\damon\data\VanOrley\VanOrley_ChristwithDoctors_K.img';
ir_extract(ir_fn2,0,usegui)
rough_mosaic(ir_fn2,usegui)

art_register(rgb_fn,sc,ir_fn1,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui)

rgb_fn = 'VanOrley_ChristwithDoctors_H_IR.tif';
sc = 1;
art_register(rgb_fn,sc,ir_fn2,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui)
