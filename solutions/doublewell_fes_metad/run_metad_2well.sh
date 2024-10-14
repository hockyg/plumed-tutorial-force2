#!/bin/bash
source /scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/GMXRC.bash.modules
export PATH=$PATH:/scratch/work/hockygroup/software/plumed2-gcc-eds-masterclass/bin/

for force in 0.0 -1.0 -0.5 0.5 1.0 -2.5 2.5;do
    outdir=run_F$force
    outprefix=run_metad_F$force
    mkdir -p $outdir
    sed -e "s/_FORCE_/$force/" -e "s/_OUTPREFIX_/$outprefix/g" doublewell_prod.plumed.dat > $outdir/doublewell_prod.plumed.dat
    cat $outdir/doublewell_prod.plumed.dat
    pushd $outdir
    plumed --no-mpi pesmd < ../doublewell_prod.pesmd.input
    popd
done

