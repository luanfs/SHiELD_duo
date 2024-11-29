#!/bin/bash
# Luan Santos - 2023
# Script to run sw test cases for different parameters


#----------------------------------------------------------------------------------------------
# Let us compile the code first
COMP="repro"      # choices:  debug, repro, prod
#COMP="debug"
clear
cd ../compilation
./compile_sw.sh $COMP noclean
cd -

#----------------------------------------------------------------------------------------------
# Modify this part for your purposes
# The output may plotted using the python scripts
# plot_scalar_field.py, plot_converge.py and plot_error_graph.py

#test case
tc="2"

# 1d advection scheme
#hords=(5 8)
hords=(8)

# 2d advection scheme
#advs=(1 1 1 1 1 1)
advs=(1 1 2 2)
#advs=(1 2)
#advs=(1 1 1 1)
#advs=(2 2 2 2)
#advs=(2 2)
#advs=(1)

# Simulation arrays
#grid_type (0-equiedge; 2-equiangular)
#gtypes=(2)
#gtypes=(2 2 2 0 0 0)
gtypes=(0 2 0 2)
#gtypes=(0 0 2 2)
#gtypes=(2 2 2 2 2 2)
#gtypes=(2 2)
#gtypes=(0 2)
#gtypes=(2)
#gtypes=(0 0 0 0)
#gtypes=(2 2 2 2)

# duogrid schemes (0-none; 1-interpolation based on geodesic distances; 2-interpolation base on cube distances)
dgs=(2 2 2 2 2 2)
#dgs=(1 1 1 1 1 1)
#dgs=(0 0 0 0 0 0)

#mass fixer
mfs=(1 1 1 1)
#mfs=(0 1 0 1)

#divergence damping (only for sw)
# variables for hyperthreading
if [ $tc = "2" ]; then
   dds=(0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12)
else
   dds=(0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12)
   #dds=(0 0 0 0 0 0 0 0 0 0)
   #dds=(0 0.12 0 0.12)
fi
#dds=(0 0.12 0 0.12)

# rotation angles
alphas=(45 45 45 45 45 45 45 45)
#alphas=(0 0 0 0 0 0)

# Initial value of N
N=48

#number of grids to be tested (we double the values N for each new grid and divide dt by 2)
Ng=4

#----------------------------------------------------------------------------------------------


case $tc in
    #-------------SW tests-------------
    2)
    tcname="geobalance"
    days="3"
    dt=3600
    n_split=7
    ;;

    5)
    tcname="mountain"
    days="14"
    dt=1800
    n_split=7
    ;;

    6)
    tcname="rhwave"
    days="105"
    dt=1200
    n_split=0
    ;;

    7)
    tcname="galewski"
    days="8"
    dt=1800
    n_split=8
    ;;

    110)
    tcname="divwind"
    days="7"
    dt=3200
    n_split=7
    ;;
esac

#----------------------------------------------------------------------------------------------
# run length
hours="0"
minutes="0"
seconds="0"
#----------------------------------------------------------------------------------------------

# Loop over all schemes
size=${#gtypes[@]}
h=${#hords[@]}

# Loop on adv scheme
# Perform the loop from 1 to N
for ((j=1; j<=$Ng; j++)); do
    for ((k=0; k<=$h-1; k++)); do
        hord=${hords[k]}
        for ((i=0; i<=size-1; i++)); do
            gtype=${gtypes[i]}
            dg=${dgs[i]}
            alpha=${alphas[i]}
            dd=${dds[i]}
            adv=${advs[i]}
            mf=${mfs[i]}
 
            if [[ $tc -ne 1 && $tc -ne 2 ]]; then
                alpha=0
            fi
            # Run the code shallow water script
            sbatch ./fv3_sw.csh $tc $N $dt $gtype $dg $tcname $alpha $days $hours $minutes $seconds $COMP $n_split $hord $dd $adv $mf
            #echo N$N g$gtype dg$dg a$alpha dt$dt tc$tc
            #sleep 10
        done
    done
    # Update the values of N and dt for the next iteration
    N=$((N * 2))
    dt=$((dt / 2))
done
