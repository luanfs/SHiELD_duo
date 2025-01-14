#!/bin/tcsh
# Gaea 
#SBATCH --output=/ncrc/home2/Luan.Santos/SHiELD_duo/scripts/run/stdout/%x.o%j
#SBATCH --job-name=solo.nh
#SBATCH --partition=batch
#SBATCH --qos=urgent
#SBATCH --account=gfdl_w
#SBATCH --time=1:00:00
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
##SBATCH --time=1:00:00
##SBATCH --cluster=stellar
##SBATCH --ntasks=1000

#set BUILD_AREA = "/home/ls9640/SHiELD_duo/SHiELD_build" 
#set SCRATCHROOT = "/scratch/cimes/ls9640"
#set SCRIPT_AREA = /home/ls9640/SHiELD_duo/SHiELD_build
#set stellar environement
#source /home/ls9640/workspace_stellar/site/environment.stellar.sh_ok

##################################################################################
# Simulation parameters
set adv=2             # 1-Putman and Lin 2007 scheme; 2-LT2
set dg=1              # duogrid (always 1)
set hord=5            # PPM scheme
set N=128             # N
set npz="63"
set gtype=0

#set TYPE="nh"           # choices:  nh, hydro
set TYPE="hydro"           # choices:  nh, hydro

set Tf="5"
set dt_atmos="600"    # atmos time step
set n_split="8"       # 
set div_damp=0.25     # divergence damping coefficient
set dgflag=".true."
set test_case="55"
set testname="tc"
set layout=5
##################################################################################

# set vorticity damping coefficient
set vort_damp=0.12


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

set OUTDIR="${GRID}.${MEMO}.g$gtype.$dgname.adv$adv.hord$hord"
#set OUTDIR="${GRID}.${MEMO}.g$gtype.$dgname.adv$adv.hord$hord.vd$vort_damp"
#set adv=1             # 1-Putman and Lin 2007 scheme; 2-LT2

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

# set variables in input.nml for initial run
set ecmwf_ic=".F."
set mountain=".F."
set external_ic=".F."
set warm_start=".F."
set na_init=1
set curr_date="0,0,0,0"

# fms yaml
set use_yaml=".F." #if True, requires data_table.yaml and field_table.yaml

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


#copy over the other tables and executable
cp ${SCRIPT_AREA}/tables/data_table data_table 

# build the date for curr_date from DATE
set curr_date = "0,0,0,0"

# build the diag_table with the experiment name and date stamp
cat > diag_table << EOF
${GRID}.${MODE}
0 0 0 0 0 0
"grid_spec",              -1,  "months",   1, "days",  "time"
"atmos_static",           -1,  "hours",    1, "hours", "time"
"atmos_4xdaily",           6, "hours",  1, "days",  "time"
"atmos_4xdaily_ave",       6, "hours",  1, "days",  "time"
#output variables
#=======================
# ATMOSPHERE DIAGNOSTICS
#=======================
###
# grid_spec
###
 "dynamics", "grid_lon", "grid_lon", "grid_spec", "all", .false.,  "none", 2,
 "dynamics", "grid_lat", "grid_lat", "grid_spec", "all", .false.,  "none", 2,
 "dynamics", "grid_lont", "grid_lont", "grid_spec", "all", .false.,  "none", 2,
 "dynamics", "grid_latt", "grid_latt", "grid_spec", "all", .false.,  "none", 2,
 "dynamics", "area",     "area",     "grid_spec", "all", .false.,  "none", 2,
### 850mb
 "dynamics",  "vort850",        "vort850",       "atmos_4xdaily", "all", .false.,  "none", 2
 "dynamics",  "u850",        "u850",       "atmos_4xdaily", "all", .false.,  "none", 2
 "dynamics",  "v850",        "v850",       "atmos_4xdaily", "all", .false.,  "none", 2

#### Lowest-layer fields
 "dynamics",  "us",          "us",        "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "vs",          "vs",        "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "tb",          "tb",        "atmos_4xdaily", "all", .false., "none", 2


#### Integrated fields
 "dynamics",  "tq",          "PWAT",        "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "lw",          "VIL",         "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "iw",          "iw",          "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "ps",          "PRESsfc",     "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "ctp",         "PRESctp",      "atmos_4xdaily", "all", .false., "none", 2
 "dynamics",  "tm",          "TMP500_300", "atmos_4xdaily", "all", .false.,  "none", 2
 "dynamics",  "cond",        "condensation", "atmos_4xdaily", "all", .false.,  "none", 2
 "dynamics",  "reevap",      "evaporation", "atmos_4xdaily", "all", .false.,  "none", 2

 "dynamics",  "prec",        "prec", "atmos_4xdaily_ave", "all", .true.,  "none", 2
 "dynamics",  "wmaxup",      "wmaxup", "atmos_4xdaily_ave", "all", "max",  "none", 2

###
# gfs static data
###
 "dynamics",      "pk",          "pk",           "atmos_static",      "all", .false.,  "none", 2
 "dynamics",      "bk",          "bk",           "atmos_static",      "all", .false.,  "none", 2
 "dynamics",      "hyam",        "hyam",         "atmos_static",      "all", .false.,  "none", 2
 "dynamics",      "hybm",        "hybm",         "atmos_static",      "all", .false.,  "none", 2
 "dynamics",      "zsurf",       "HGTsfc",          "atmos_static",      "all", .false.,  "none", 2
EOF


#
# build the field_table
cat > field_table <<EOF
# added by FRE: sphum must be present in atmos
# specific humidity for moist runs
 "TRACER", "atmos_mod", "sphum"
           "longname",     "specific humidity"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=1.e30" /
# prognostic cloud water mixing ratio
 "TRACER", "atmos_mod", "liq_wat"
           "longname",     "cloud water mixing ratio"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=1.e30" /
 "TRACER", "atmos_mod", "rainwat"
           "longname",     "rain mixing ratio"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=1.e30" /
 "TRACER", "atmos_mod", "ice_wat"
           "longname",     "cloud ice mixing ratio"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=1.e30" /
 "TRACER", "atmos_mod", "snowwat"
           "longname",     "snow mixing ratio"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=1.e30" /
 "TRACER", "atmos_mod", "graupel"
           "longname",     "graupel mixing ratio"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=1.e30" /

# non-prognostic cloud amount
 "TRACER", "atmos_mod", "cld_amt"
           "longname",     "cloud amount"
           "units",        "1"
       "profile_type", "fixed", "surface_value=1.e30" /
EOF


#
# build the input.nml
cat > input.nml <<EOF
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

 &fms_affinity_nml
    affinity = .false.
/

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
       layout   = $layout_x,$layout_y
       io_layout = $io_layout
       npx      = $npx
       npy      = $npy
       ntiles   = 6
       npz    = $npz
       npz_type = 'gfs'
       grid_type = 0
       make_nh = $make_nh
       range_warn = .T.
       reset_eta = .F.
       sg_cutoff = 150.e2 !replaces old "n_sponge"
       fv_sg_adj = 3600
       nudge_qv = .T.
       RF_fast = .F.
       tau_h2o = 0.
       tau = 10.
       rf_cutoff = 30.e2
       d2_bg_k1 = 0.20
       d2_bg_k2 = 0.10 ! z2: increased
       kord_tm = -8
       kord_mt = 8
       kord_wz = 8
       kord_tr = 8
       hydrostatic = $hydrostatic
       phys_hydrostatic = $phys_hydrostatic
       use_hydro_pressure = $use_hydro_pressure
       beta = 0.
       a_imp = 1.
       p_fac = 0.05
       k_split = 2
       n_split = 8
       nwat = 6
       na_init = $na_init
       d_ext = 0.0
       dnats = 1
       d2_bg = 0.
       nord = 3 ! z8: 2 --> 3
       dddmp = 0.5
       d4_bg = 0.15 ! z8: increased with nord change
       vtdm4 = 0.06 ! z10: increased
       delt_max = 0.002
       convert_ke = .true.
       ke_bg = 0.
       do_vort_damp = .T.
       external_ic = .F. !COLD START
       !is_ideal_case = .T.
       mountain = .F.
       d_con = 1.
       hord_mt = 5
       hord_vt = 5
       hord_tm = 5
       hord_dp = -5
       hord_tr = -5 ! z2: changed
       adjust_dry_mass = .F.
       consv_te = 0.0
       fill = .F.
       dwind_2d = .F.
       print_freq = 6
       warm_start = $warm_start
       z_tracer = .T.
       fill_dp = .T.
       adiabatic = .F.

       !do_cube_transform = .true. !Replaces do_schmidt
       !target_lat = 17.5
       !target_lon = 172.5
       !stretch_fac = 3.0

       duogrid    = $dgflag
       duogrid_scheme = $dg
       adv_scheme = $adv


/

 &integ_phys_nml
       do_inline_mp = .T.
       do_sat_adj = .F.
/

&fv_diag_plevs_nml
    nplev=7
    levs = 50,200,300,500,700,850,1000
    levs_ave = 10,100,300,500,700,850,1000
/

 &sim_phys_nml
  do_reed_sim_phys = .T.
!  do_gfdl_sim_phys = .T.
/

&reed_sim_phys_nml
   do_reed_cond = .F.
   reed_test = 0 !uniform SST
/


&gfdl_sim_phys_nml
  sst0 = 302.15
  !mixed_layer = .T.
  uniform_sst = .T.
  !do_abl = .T.
  do_mon_obkv = .T.
/


 &ocean_rough_nml
        do_highwind = .T.
        do_cap40 = .T.
        v10m = 32.5
        v10n = 17.5
        do_init = .T.
        zcoh1 = 0.
        zcoq1 = 0. !  1.3e-4
        rough_scheme = 'beljaars'
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

 &external_ic_nml
       filtered_terrain = .T.
       levp = 64
       gfs_dwinds = .T.
       checker_tr = .F.
       nt_checker = 0
/

 &gfdl_mp_nml
       do_sedi_heat = .F.
       do_sedi_w = .F.
       rad_snow = .true.
       rad_graupel = .true.
       rad_rain = .true.
       const_vi = .F.
       const_vs = .F.
       const_vg = .F.
       const_vr = .F.
       vi_max = 1.
       vs_max = 2.
       vg_max = 16.
       vr_max = 16.
       qi_lim = 1.
       prog_ccn = .false.
       do_qa = .true.
       tau_l2v = 300.
       tau_v2l = 90. ! z7: enabled
       do_cond_timescale = .true. ! z7: enabled
       rthresh = 10.e-6  !   10.e-6  ! This is a key parameter for cloud water
      dw_land  = 0.15
      dw_ocean = 0.10
       ql_gen = 1.0e-3
    ql_mlt = 2.0e-3
    qs_mlt = 1.e-6
       qi0_crt = 8.E-5
       qs0_crt = 3.0e-3
       tau_i2s = 1000.
       c_psaci = 0.05
       c_pgacs = 0.01
       rh_inc = 0.0
       rh_inr = 0.0
       rh_ins = 0.0
       ccn_l = 300.
       ccn_o = 100.
       c_paut =  0.55
       z_slope_liq  = .T.
       z_slope_ice  = .T.
       fix_negative = .true.
       irain_f = 0
       icloud_f = 0
/

!# From LJZ mar 2019
 &cld_eff_rad_nml
       reiflag        = 4
        rewmin = 5.0
        rewmax = 10.0
        reimin = 10.0
        reimax = 150.0
/

 &test_case_nml
    test_case = $test_case
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
