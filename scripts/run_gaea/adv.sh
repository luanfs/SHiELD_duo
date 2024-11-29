#!/bin/bash
# Luan Santos - 2023
# Script to run advection/sw test cases for different parameters


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
tc="-6"

# 1d advection scheme
hords=(0 8)
#hords=(8)

# 2d advection scheme
#advs=(1 2 1 2 1 2 1 2)
advs=(1 2 1 2)
#advs=(2)

# Simulation arrays
#grid_type (0-equiedge; 2-equiangular)
#gtypes=(2 2 0 0 2 2 0 0)
gtypes=(2 2 0 0)
#gtypes=(0 2 0 2)
#gtypes=(2 0)
#gtypes=(0 0 0)
#gtypes=(2)

#mass fixer
#mfs=(1 1 1 1 0 0 0 0)
mfs=(1 1 1 1) #always true

#midpoint cube
#mps=(".false." ".false." ".true.")
mps=(".false." ".false." ".false." ".false." ".false." ".false." ".false." ".false." ".false." ".false." ".false." ".false.")

# duogrid schemes (0-none; 1-interpolation based on geodesic distances; 2-interpolation base on cube distances)
dgs=(2 2 2 2 2 2 2 2)
#dgs=(1 1 2 2)
#dgs=(0 2 2)
#dgs=(2 2)

#divergence damping (only for sw)
# variables for hyperthreading
if [ $tc = "2"  ]; then
   dds=(0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12)
else
   dds=(0 0 0 0 0 0 0 0 0 0 0 0)
fi

# rotation angles
alphas=(45 45 45 45 45 45 45 45)

# Initial value of N
N=48

#number of grids to be tested (we double the values N for each new grid and divide dt by 2)
Ng=4

#----------------------------------------------------------------------------------------------


case $tc in
    #------------- advection tests-------------
    1)
    tcname="cosine-zonal"
    days="12"
    dt=3600
    n_split=1
    ;;

    -3)
    tcname="gaussian-zonal"
    days="12"
    dt=3600
    n_split=1
    ;;

    -4)
    tcname="geobalance"
    days="12"
    dt=3600
    n_split=1
    ;;

    -5)
    tcname="twogaussians-ndiv"
    days="12"
    dt=1600
    n_split=1
    ;;

    -6)
    tcname="twogaussians-div"
    days="12"
    dt=6400
    #dt=2592
    n_split=1
    ;;

    -7)
    tcname="cylinder-zonal"
    days="12"
    dt=3600
    n_split=1
    ;;

    -8)
    tcname="cylinder-ndiv"
    days="12"
    dt=1600
    n_split=1
    ;;

    -9)
    tcname="cylinder-div"
    days="12"
    dt=6400
    n_split=1
    ;;
 
    -2)
    tcname="duogrid-test"
    days="12"
    dt=3600
    n_split=1
    ;;

    -10)
    tcname="recon-test"
    days="12"
    dt=3600
    n_split=1
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
   
            if [[ $tc == -6 || $tc == -5  || $tc == -8 || $tc == -9 ]]; then
                alpha=0
            fi
            # Run the code shallow water script
            echo N$N g$gtype dg$dg a$alpha  tc$tc dt$dt
            sbatch ./fv3_sw.csh $tc $N $dt $gtype $dg $tcname $alpha $days $hours $minutes $seconds $COMP $n_split $hord $dd $adv $mf
	    #exit

            #sleep 10
        done
    done
    # Update the values of N and dt for the next iteration
    N=$((N * 2))
    dt=$((dt / 2))
done
