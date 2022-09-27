source /glade/u/home/jsimkins/.bashrc

cp /glade/work/jsimkins/runs-MOM6/nwa25/RESTART/* /glade/work/jsimkins/runs-MOM6/nwa25/INPUT/

cp restart.input.nml input.nml

qsub /glade/work/jsimkins/runs-MOM6/nwa25/mom.sub


