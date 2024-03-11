# exoplasim to genie regrid
processing code for genie, exoplasim

currently: 
(1) reads in 'MOST_0009X.nc' files for the last 10 years of an exoplasim model run
(2) extracts core data
    (a) wind products
    (b) planetary albedo
(3) returns genie input files
    (a) wind field input (tau, spd) files 
    (b) 1-d and 2-d planetary albedo files

in genie userconfig:
set bg_ctrl_force_windspeed=.true
for exoplasim v. 3.0.6 
  - set tau scaling (ea_11, go_13) to 2.0
  - set bg_par_gastransfer_a=0.715
for exoplasim v. 3.3.0
  - set tau scaling (ea_11, go_13) to 2.6
  - set bg_par_gastransfer_a=1.1
