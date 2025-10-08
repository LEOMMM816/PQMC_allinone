PQMC_allinone_git

0. Originally intended for one specific QMC project, now it is extended to be able to manage all kinds of lattice and interactions as the repo name implies. 
One only need to edit the model_file as input. There are two example files in the model folder.
Not for quasi-crystal. Measurables still needs some case-by-case coding. 

1. To change model, edit the .nml files in ./model. 
Then change the TEMPLATE accordingly in .sh file.

2. For now, most input parameters are still only accessible through ./code/input.f90.
So is the adjustments to the measurement module, ./code/mod_meas.f90

3. The raw data is stored in ./data/data/.
More readable data are in ./data/out_files/.

4. To debug, run testlocal.sh. To put on cluster, sbatch mpi_32c.sh
