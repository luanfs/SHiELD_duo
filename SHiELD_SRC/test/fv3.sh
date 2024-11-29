#!/bin/sh
#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the SHiELD Build System.
#*
#* The SHiELD Build System free software: you can redistribute it
#* and/or modify it under the terms of the
#* GNU Lesser General Public License as published by the
#* Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* The SHiELD Build System distributed in the hope that it will be
#* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
#* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#* See the GNU General Public License for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with theSHiELD Build System
#* If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
#
#  DISCLAIMER: This script is provided as-is for the purpose of continuous
#              integration software regression testing and as such is
#              unsupported.
#
#SBATCH --ntasks=24
BUILDDIR="/SHiELD_build/"
if [ -z ${BUILDDIR} ] ; then
  echo -e "\nERROR:\tset BUILDDIR to the base path /<path>/SHiELD_build/ \n"
  exit 99
fi

# when we are running this test for CI, we run it in a container
# an expected value for container would be "--mpi=pmi2 singularity exec -B <file directory> <container>"
# this is needed within the run command
#if [ -z "$1" ] ; then
#  CONTAINER=""
#else
#  CONTAINER=$1
#  echo -e "\nThis test will be run inside of a container with the command $1"
#fi

COMPILER="gnu"
# configure the site
COMPILER=${COMPILER:-intel}
. ${BUILDDIR}/site/environment.${COMPILER}.sh

set -x

# necessary for OpenMP
export OMP_STACKSIZE=10g


# Check if three command-line arguments are provided
if [ "$#" -ne 17 ]; then
  echo "Usage: ./sw.sh tc N dt gtype duogrid tcname alpha mp adv days hours minutes seconds comp nsplit hord cw dd"
  exit 1
fi

# Assign the command-line arguments to variables
tc="$1"
N="$2"
dt_atmos="$3"
gtype="$4"
dg="$5"
tcname="$6"
alpha="$7"
days="${8}"
hours="${9}"
minutes="${10}"
seconds="${11}"
COMP="${12}"
n_split="${13}"
hord="${14}"
dd="${15}"
adv="${16}"
mf="${17}"
plotfreq=$((days * 1))
#plotfreq=1

#echo $tc $N $dt_atmos $gtype $tcname $alpha $days $hours $minutes $seconds $COMP


# case specific details
res=$N
MEMO="sw."$tcname # trying repro executable
TYPE="sw"         # choices:  nh, hydro
MODE="64bit"      # choices:  32bit, 64bit
GRID="C$res"
HYPT="on"         # choices:  on, off  (controls hyperthreading)
#COMP="debug"      # choices:  debug, repro, prod

# variables for hyperthreading
# duogrid scheme
if [ $dg = "1" ] ; then
  dgname="dg1"
  dgflag=".true."
elif [ $dg = "2" ] ; then
  dgname="dg2"
  dgflag=".true."
else
  dgname="kinked"
  dgflag=".false."
fi
echo $dgname

# midpoint formulation
#if [ $mp = ".true." ] ; then
#  mpname="c"
#elif [ $mp = ".false." ] ; then
#  mpname="s"
#fi
#echo $mpname $mp

# variables for hyperthreading
if [ "$tc" = "-2" ]; then
  WORKDIR="/SHiELD_SRC/test/CI/BATCH-CI/${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$mpname.$dgname/"
# variables for hyperthreading
elif [ "$tc" = "-10" ]; then
  WORKDIR="/SHiELD_SRC/test/CI/BATCH-CI/${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$mpname.$dgname.hord$hord/"
elif [ "$tc" -le 1 ]; then
  #WORKDIR="/SHiELD_SRC/test/CI/BATCH-CI/${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$mpname.$dgname.adv$adv.hord$hord.tf$days/"
  WORKDIR="/SHiELD_SRC/test/CI/BATCH-CI/${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$dgname.adv$adv.hord$hord.mf$mf.tf$days/"
else
  #WORKDIR="/SHiELD_SRC/test/CI/BATCH-CI/${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$mpname.$dgname.adv$adv.hord$hord.dd$dd.tf$days/"
  WORKDIR="/SHiELD_SRC/test/CI/BATCH-CI/${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$dgname.adv$adv.hord$hord.dd$dd.mf$mf.tf$days/"
fi


# directory structure
#WORKDIR=${SCRATCHDIR:-${BUILDDIR}}/CI/BATCH-CI/${GRID}.${MEMO}/

executable=${BUILDDIR}/Build/bin/SOLO_${TYPE}.${COMP}.${MODE}.${COMPILER}.x

# convert alpha from degree to rad
alpha=$(awk "BEGIN { printf \"%.10f\", $alpha * 0.01745329251 }")

# changeable parameters
# dycore definitions
npx=$(($N+1))
npy=$(($N+1))
npz="1" #Shallow water
layout_x="1"
layout_y="1"
io_layout="1,1" #Want to increase this in a production run??
nthreads="2"

# set variables in input.nml for initial run
na_init=0 # TRY 1
curr_date="0,0,0,0"

make_nh=".F."
hydrostatic=".T."
phys_hydrostatic=".F."     # will be ignored in hydro mode
use_hydro_pressure=".T."   # have to be .T. in hydro mode
consv_te="0."

# variables for hyperthreading
if [ ${HYPT} = "on" ] ; then
  hyperthread=".true."
  div=2
else
  hyperthread=".false."
  div=1
fi
skip=`expr ${nthreads} \/ ${div}`

# when running with threads, need to use the following command
npes=`expr ${layout_x} \* ${layout_y} \* 6 `
LAUNCHER=${LAUNCHER:-srun}
if [ ${LAUNCHER} = 'srun' ] ; then
    export SLURM_CPU_BIND=verbose
    run_cmd="${LAUNCHER} --label  --allow-run-as-root --ntasks=$npes --cpus-per-task=$skip $CONTAINER ./${executable##*/}"
else
    export MPICH_ENV_DISPLAY=YES
    run_cmd="${LAUNCHER}  --allow-run-as-root -np $npes $CONTAINER ./${executable##*/}"
fi

# set up the run area
\rm -rf $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR
mkdir -p RESTART
mkdir -p INPUT

# copy over the executable
cp $executable .

# copy over the tables
#
# create an empty data_table
touch data_table


#
# build the field_table
cat > field_table <<EOF

 "TRACER", "atmos_mod", "sphum"
           "longname",     "specific humidity"
           "units",        "kg/kg"
           "profile_type", "fixed", "surface_value=3.e-6" /
EOF


#
# build the diag_table with the experiment name and date stamp
cat > diag_table << EOF
${GRID}.${MODE}
0 0 0 0 0 0
"grid_spec",    -1,  "hours",  1, "days", "time",
"atmos_daily",  24,  "hours",  1, "days", "time",

"dynamics", "grid_lon", "grid_lon", "grid_spec", "all", .false.,  "none", 2,
"dynamics", "grid_lat", "grid_lat", "grid_spec", "all", .false.,  "none", 2,
"dynamics", "area", "area", "grid_spec", "all", .false.,  "none", 2,

"dynamics", "ps_ic", "ps_ic",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "ps",    "ps",      "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "ua_ic", "ua_ic",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "va_ic", "va_ic",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "ucomp", "ucomp",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "vcomp", "vcomp",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "vort", "vort",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "pv", "pv",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "delp", "delp",   "atmos_daily", "all", .false.,  "none", 2,

EOF

#
# build the input.nml
cat > input.nml <<EOF
 &fms_affinity_nml
    affinity = .false.
/
&fms_io_nml
  checksum_required   = .false.
  max_files_r = 100,
  max_files_w = 100,
/

&fms_affinity_nml
  affinity=.F.
/

&fms_nml
  clock_grain = 'ROUTINE',
  domains_stack_size = 2000000000,
  print_memory_usage = .false.
/

&fv_core_nml
  do_vort_damp=.false.
  vtdm4=0.0
  layout   = $layout_x,$layout_y
  io_layout = $io_layout
  npx      = $npx
  npy      = $npy
  ntiles   = 6
  npz    = $npz
  nwat = 0
  grid_type = $gtype
  na_init = 0
  dnats = 0
  nord = 2
  d4_bg = $dd
  mountain = .F.
  hord_mt = $hord
  hord_vt = $hord
  hord_tm = $hord
  hord_dp = $hord
  hord_tr = $hord
  n_split = $n_split
  print_freq = $plotfreq
  warm_start = .F.
  duogrid    = $dgflag
  duogrid_scheme = $dg
  adv_scheme = $adv
  mass_fixer = $mf
  do_schmidt = .false.
  adiabatic = .true.
/

&test_case_nml
    test_case = $tc
    alpha = $alpha
/

 &main_nml
       days  = $days
       hours = $hours
       minutes = $minutes
       seconds = $seconds
       dt_atmos = $dt_atmos
       current_time =  $curr_date
       atmos_nthreads = 4
       use_hyper_thread = .true.
/

EOF



# run the executable
${run_cmd} | tee fms.out || exit
#
# verification
# < add logic for verification >
