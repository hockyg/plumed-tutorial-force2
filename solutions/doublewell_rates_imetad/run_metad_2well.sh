#!/bin/bash
source /scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/GMXRC.bash.modules
export PATH=$PATH:/scratch/work/hockygroup/software/plumed2-gcc-eds-masterclass/bin/

n_runs=40

for force in 0.0 -0.2 -0.4 -0.6 -0.8 -1.0;do
    for seed in `seq 1 $n_runs`;do
        outdir=run_F$force/$seed
        outprefix=run_metad_F$force
	if [ ! -d $outdir ];then
            mkdir -p $outdir
            sed -e "s/_FORCE_/$force/" -e "s/_OUTPREFIX_/$outprefix/g" doublewell_prod.plumed.dat > $outdir/doublewell_prod.plumed.dat
    	sed -e "s/__SEED__/$seed/" ./doublewell_prod.pesmd.input > $outdir/doublewell_prod.pesmd.input
            pushd $outdir
            plumed --no-mpi pesmd < ./doublewell_prod.pesmd.input
            popd
	fi
    done
done

