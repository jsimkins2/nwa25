# this bash script is for spawning new mom6 runs

source /glade/u/home/jsimkins/.bashrc

rm /glade/work/jsimkins/runs-MOM6/nwa25/INPUT/MOM.res*
rm /glade/work/jsimkins/runs-MOM6/nwa25/INPUT/coupler.res
rm /glade/work/jsimkins/runs-MOM6/nwa25/INPUT/ice_model.res.nc
rm /glade/work/jsimkins/runs-MOM6/nwa25/RESTART/*

cp new.input.nml input.nml

qsub /glade/work/jsimkins/runs-MOM6/nwa25/mom.sub


