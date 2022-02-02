# this bash script is for spawning new mom6 runs

source /glade/u/home/jsimkins/.bashrc

rm /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/INPUT/MOM.res*
rm /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/coupler.res
rm /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/ice_model.res.nc
rm /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/RESTART/*

cp new.input.nml input.nml

qsub /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/mom.sub


