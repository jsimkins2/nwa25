# this bash script is for spawning new mom6 runs

source /ncrc/home1/James.Simkins/.bashrc

rm $sim/INPUT/MOM.res*
rm $sim/INPUT/coupler.res
rm $sim/INPUT/ice_model.res.nc
rm $sim/RESTART/*

cd $sim
cp new.input.nml input.nml

sbatch $sim/mom.sub


