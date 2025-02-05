#!/bin/tcsh
# Gaea 
#SBATCH --output=/ncrc/home2/Luan.Santos/SHiELD_duo/scripts/run/stdout/%x.o%j
#SBATCH --job-name=solo.nh
#SBATCH --partition=batch
#SBATCH --qos=urgent
#SBATCH --account=gfdl_w
#SBATCH --time=3:00:00
#SBATCH --cluster=c5
#SBATCH --ntasks=1000

set BUILD_AREA = "/ncrc/home2/Luan.Santos/SHiELD_duo/SHiELD_build" 
set SCRATCHROOT = "/gpfs/f5/gfdl_w/scratch/Luan.Santos"
set SCRIPT_AREA = /ncrc/home2/Luan.Santos/SHiELD_duo/SHiELD_build

# Stellar
##SBATCH --output=/home/ls9640/SHiELD_duo/scripts/run/stdout/%x.o%j
##SBATCH --job-name=solo.nh
##SBATCH --partition=batch
##SBATCH --qos=urgent
##SBATCH --account=cimes2
##SBATCH --time=3:00:00
##SBATCH --cluster=stellar
##SBATCH --ntasks=1000

#set BUILD_AREA = "/home/ls9640/SHiELD_duo/SHiELD_build" 
#set SCRATCHROOT = "/scratch/cimes/ls9640"
#set SCRIPT_AREA = /home/ls9640/SHiELD_duo/SHiELD_build
#set stellar environement
#source /home/ls9640/SHiELD_duo/SHiELD_build/site/environment.stellar.sh_ok

##################################################################################
# Simulation parameters
set adv=2             # 1-Putman and Lin 2007 scheme; 2-LT2
set dg=1              # duogrid (always 1)
set gtype=0           # grid type(0-equiedge; 2-equiangular)
set hord=5            # PPM scheme
set N=128             # N
set npz="5"

set TYPE="nh"           # choices:  nh, hydro
set TYPE="hydro"           # choices:  nh, hydro

set Tf="100"
set dt_atmos="1200"    # atmos time step
set n_split="17"
set dgflag=".true."
set test_case="45"
set testname="modons3d"
set layout=5
##################################################################################

# set divergence/vorticity damping coefficient
if ($hord == "5") then
   if ($adv == "1") then
      set div_damp=0.12
      set vort_damp=0.12
   else
      set div_damp=0.15
      set vort_damp=0.12
   endif
else if ($hord == "6") then
   if ($adv == "1") then
      set div_damp=0.12
      set vort_damp=0.06
   else
      set div_damp=0.12
      set vort_damp=0.12
   endif
else if ($hord == "8") then
   if ($adv == "1") then
      set div_damp=0.12
      set vort_damp=0.12
   else
      set div_damp=0.18
      set vort_damp=0.06
   endif
else
   set vort_damp=0
endif


##################################################################################
# rotation angles
set alpha_deg=0
set alpha = `awk 'BEGIN { printf "%.10f", '"$alpha_deg"' * 0.01745329251 }'`
##################################################################################

# case specific details
set res=$N
set MEMO=$TYPE"."$testname # trying repro executable
set MODE="64bit"        # choices:  32bit, 64bit
set GRID="C$res"
set HYPT="on"         # choices:  on, off  (controls hyperthreading)
set COMP="repro"       # choices:  debug, repro, prod
set RELEASE = "solo_"$TYPE         # run cycle, 1: no restart # z2: increased
set EXE  = "intel.x"


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

if ($vort_damp == "0") then
  set do_vort_damp=.false.
else
  set do_vort_damp=.true.
endif

set OUTDIR="${GRID}.${MEMO}.g$gtype.$dgname.adv$adv.hord$hord.dd$div_damp.vd$vort_damp"
#set OUTDIR="${GRID}.${MEMO}.g$gtype.$dgname.adv$adv.hord$hord.vd$vort_damp"

# directory structure
set WORKDIR =  ${SCRATCHROOT}/${RELEASE}/${OUTDIR}
#set executable = ${BUILD_AREA}/Build/bin/SOLO_${TYPE}.${COMP}.${MODE}.x
set executable = ${BUILD_AREA}/Build/bin/SOLO_${TYPE}.${COMP}.${MODE}.${EXE}



# changeable parameters
# dycore definitions
@ Np1 = $res + 1
set npx=$Np1
set npy=$Np1
set layout_x=$layout
set layout_y=$layout
set io_layout="1,1"
set nthreads="2"

# run length
set days=$Tf
set hours="0"
set minutes="0"
set seconds="0"

# set variables in input.nml for initial run
set na_init=0
set curr_date="0,0,0,0"


if (${TYPE} == "nh")then
   # non-hydrostatic options
   set make_nh=".T."
   set hydrostatic=".F."
   set phys_hydrostatic=".F."     # can be tested
   set use_hydro_pressure=".F."   # can be tested
   set consv_te="1."
else
   # hydrostatic options
   set make_nh=".F."
   set hydrostatic=".T."
   set phys_hydrostatic=".F."     # will be ignored in hydro mode
   set use_hydro_pressure=".T."   # have to be .T. in hydro mode
   set consv_te="0."
endif


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
setenv KMP_STACKSIZE 4g
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


#copy over the other tables and executable
cp ${SCRIPT_AREA}/tables/data_table data_table 
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
"dynamics", "u500", "u500",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "v500", "v500",   "atmos_daily", "all", .false.,  "none", 2,
"dynamics", "vort500", "vort500",   "atmos_daily", "all", .false.,  "none", 2,


EOF

cat >! input.nml <<EOF

 &diag_manager_nml 
    !flush_nc_files = .false.
    prepend_date = .F.
! this diag table creates a lot of files
! next three lines are necessary
    max_num_axis_sets = 100 
    max_files = 100
    max_axes = 240
    !do_diag_field_log = .T.
/

 &fms_io_nml
       checksum_required   = .false.
       max_files_r = 100,
       max_files_w = 100,
/

 &fms_nml
       clock_grain = 'ROUTINE',
       domains_stack_size = 16000000,
       print_memory_usage = .false.
/

 &fv_core_nml
        layout   = $layout_x,$layout_y
        mountain = .false.
        npx      = 129
        npy      = 129
        ntiles   = 6
        npz      = 5
  grid_type = 0
 reset_eta = .F.

 k_split = 1
 n_split = 16
     nwat = 0
     fill = .false.
   n_sponge = -1
 d2_bg_k1 = 0.0
 d2_bg_k2 = 0.00
 d4_bg = 0.08
           d2_bg = 0.0
	   d_ext = 0.0
 nord = 2
	   d_con = 0.
 do_vort_damp = .F.
    vtdm4 = 0.00
 consv_te = 0.0
        hord_mt = 8
        hord_vt = 8
        hord_tm = 8
        hord_dp = 8
        hord_tr = 8
      kord_tm = -9
      kord_mt =  9
      kord_wz =  9
      kord_tr =  9
        adiabatic = .true.  

        print_freq = 24
        warm_start = .false.

	!fv_debug = .true.
 check_negative = .F.

hydrostatic = .F.
   a_imp = 1.0
    beta = 0.
   p_fac = 0.25
 do_schmidt = .T.
    target_lat = 0.
    target_lon = 180
    !duogrid    = $dgflag
    !duogrid_scheme = $dg
    !adv_scheme = $adv


/

&test_case_nml ! cold start
    test_case = 45
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
