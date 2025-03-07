#!/bin/tcsh 
#SBATCH --output=/ncrc/home2/Luan.Santos/SHiELD_duo/scripts/run_gaea/stdout/%x.o%j
#SBATCH --job-name=solo.sw
#SBATCH --partition=batch
#SBATCH --qos=urgent
#SBATCH --account=gfdl_w
#SBATCH --time=1:30:00
#SBATCH --cluster=c5
#SBATCH --nodes=10


set echo
# Check if 18 command-line arguments are provided
if ($#argv != 18) then
  echo "Usage: ./sw.sh tc N dt gtype duogrid tcname alpha days hours minutes seconds comp nsplit hord adv dd vd mf"  
endif
# Assign the command-line arguments to variables
set tc="$1"
set N="$2"
set dt_atmos="$3"
set gtype="$4"
set dg="$5"
set tcname="$6"
set alpha="$7"
set days="${8}"
set hours="${9}"
set minutes="${10}"
set seconds="${11}"
set COMP="${12}"
set n_split="${13}"
set hord="${14}"
set adv="${15}"
set dd="${16}"
set vd="${17}"
set mf="${18}"
@ plotfreq = $days * 1
#set plotfreq=$((days * 1))

set BUILD_AREA = "/ncrc/home2/Luan.Santos/SHiELD_duo/SHiELD_build" 
set SCRATCHROOT = "/gpfs/f5/gfdl_w/scratch/Luan.Santos"
set SCRIPT_AREA = /ncrc/home2/Luan.Santos/SHiELD_duo/SHiELD_run

# case specific details
@ Np1 = $N + 1
set res=$N
set MEMO="sw."$tcname # trying repro executable
set TYPE="sw"         # choices:  nh, hydro
set MODE="64bit"      # choices:  32bit, 64bit
set GRID="C$res"
set HYPT="on"         # choices:  on, off  (controls hyperthreading)
set COMP = "repro"       # choices:  debug, repro, prod
set NO_SEND = "no_send"    # choices:  send, no_send
set NUM_TOT = 1         # run cycle, 1: no restart # z2: increased
set RELEASE = "solo_sw"         # run cycle, 1: no restart # z2: increased
set EXE  = "intel.x"

# variables for hyperthreading
# duogrid scheme
if ( $dg == "1" ) then
  set dgname="dg1"
  set dgflag=".true."
else if ( $dg == "2" ) then
  set dgname="dg2"
  set dgflag=".true."
else
  set dgname="kinked"
  set dgflag=".false."
endif
echo $dgname

# variables for hyperthreading
if ( "$tc" == "-2" ) then
  set OUTDIR="${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$mpname.$dgname"
else if ( "$tc" == "-10" ) then
  set OUTDIR="${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$mpname.$dgname.hord$hord"
else if ( "$tc" <= 1 ) then
  set OUTDIR="${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$dgname.adv$adv.hord$hord.mf$mf.tf$days" 
else
  set OUTDIR="${GRID}.${MEMO}.tc$tc.alpha$alpha.g$gtype.$dgname.adv$adv.hord$hord.dd$dd.vd$vd.mf$mf.tf$days"
endif
echo $OUTDIR

if ($vd == "0") then
  set do_vort_damp=.false.
else
  set do_vort_damp=.true.
endif
echo $do_vort_damp

# convert alpha to radians
set alpha = `awk 'BEGIN { printf "%.10f", '"$alpha"' * 0.01745329251 }'`

# directory structure
set WORKDIR    = ${SCRATCHROOT}/${RELEASE}/${OUTDIR}
#set executable = ${BUILD_AREA}/Build/bin/SOLO_${TYPE}.${COMP}.${MODE}.x
set executable = ${BUILD_AREA}/Build/bin/SOLO_${TYPE}.${COMP}.${MODE}.${EXE}

# input filesets
echo $WORKDIR


#changeable parameters
# dycore definitions
set npx=$Np1
set npy=$Np1
set npz = "1" #Shallow water
set layout_x = "6" 
set layout_y = "6" 
set io_layout = "1,1" #Want to increase this in a production run??
set nthreads = "2"

# run length
set days = $days
set hours = "0"
set minutes = "0"
set seconds = "0"
set dt_atmos = $dt_atmos

set make_nh = ".F."
set hydrostatic = ".T."
set phys_hydrostatic = ".F."     # will be ignored in hydro mode
set use_hydro_pressure = ".T."   # have to be .T. in hydro mode
set consv_te = "0."

# variables for hyperthreading
if (${HYPT} == "on") then
  set hyperthread = ".true."
  set div = 2
else
  set hyperthread = ".false."
  set div = 1
endif
@ skip = ${nthreads} / ${div}

# when running with threads, need to use the following command
@ npes = ${layout_x} * ${layout_y} * 6

set run_cmd = "srun --ntasks=$npes --cpus-per-task=$skip ./$executable:t"

setenv MPICH_ENV_DISPLAY
setenv MPICH_MPIIO_CB_ALIGN 2
setenv MALLOC_MMAP_MAX_ 0
setenv MALLOC_TRIM_THRESHOLD_ 536870912
setenv NC_BLKSZ 1M

# necessary for OpenMP when using Intel
setenv KMP_STACKSIZE 256m
setenv SLURM_CPU_BIND verbose

set RUN_DIR = $PWD
#if ( "$SLURM_JOB_NAME" == "sh" ) then
#  set SCRIPT = "${RUN_DIR}/$0"
#else
#  set SCRIPT = "${RUN_DIR}/$SLURM_JOB_NAME"
#endif
mkdir -p $WORKDIR/restart
set RST_COUNT = $WORKDIR/restart/rst.count
#######DEBUG
\rm -f $RST_COUNT
########END DEBUG
if ( -f ${RST_COUNT} ) then
   source ${RST_COUNT}
   if ( x"$num" == "x" || ${num} < 1 ) then
     set RESTART_RUN = "F"
   else
     set RESTART_RUN = "T"
   endif
else
   set num = 0
   set RESTART_RUN = "F"
endif

#NEED TO BE CAREFUL OF SETUP CODE WRT RESTARTS!!
if (${RESTART_RUN} == "F") then

  \rm -rf $WORKDIR/rundir

  mkdir -p $WORKDIR/rundir
  cd $WORKDIR/rundir

  mkdir -p RESTART
  mkdir -p INPUT

  # set variables in input.nml for initial run
  set na_init = 0 # TRY 1
else

  cd $WORKDIR/rundir
  \rm -rf INPUT/*

  # move the restart data into INPUT/
  #mv RESTART/* INPUT/.
  ln -s ${restart_dir}/coupler.res ${restart_dir}/[^0-9]*.nc ${restart_dir}/[^0-9]*.nc.???? INPUT/.

  # reset values in input.nml for restart run
  set na_init = 0
endif
ls INPUT/
ls RESTART/
echo ${RESTART_RUN}

# copy over the executable
cp $executable .
pwd
ls
echo $executable


#copy over the other tables and executable
cp ${SCRIPT_AREA}/common/data_table data_table 
cat >! field_table <<EOF

 "TRACER", "atmos_mod", "sphum" 
           "longname",     "specific humidity"
           "units",        "kg/kg" 
           "profile_type", "fixed", "surface_value=3.e-6" /

EOF

cp $executable .

# # build the date for curr_date from DATE
 set curr_date = "0,0,0,0"

# build the diag_table with the experiment name and date stamp
cat >! diag_table << EOF
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

cat >! input.nml <<EOF

 &fms_io_nml
       checksum_required   = .false.
       max_files_r = 100,
       max_files_w = 100,
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
  vtdm4 = $vd
  do_vort_damp = $do_vort_damp
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
       atmos_nthreads = $nthreads
       use_hyper_thread = $hyperthread
/

EOF

# run the executable
${run_cmd} | tee fms.out || exit
@ num ++
echo "set num = ${num}" >! ${RST_COUNT}

