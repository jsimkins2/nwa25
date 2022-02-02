source /glade/u/home/jsimkins/.bashrc

cp /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/RESTART/* /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/INPUT/

cp restart.input.nml input.nml

qsub /glade/work/jsimkins/runs-MOM6/old_runs/nwa25_epbl_update_12_2021/mom.sub


