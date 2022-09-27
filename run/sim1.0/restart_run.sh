source /ncrc/home1/James.Simkins/.bashrc

cp $sim/RESTART/* $sim/INPUT/

cp restart.input.nml input.nml

sbatch $sim/mom.sub


