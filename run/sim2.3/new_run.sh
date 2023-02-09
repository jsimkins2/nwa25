# this bash script is for spawning new mom6 runs

source /ncrc/home1/James.Simkins/.bashrc

rm INPUT/MOM.res*
rm INPUT/coupler.res
rm INPUT/ice_model.res.nc
rm RESTART/*

cp new.input.nml input.nml

sbatch mom.sub


