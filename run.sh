#!/bin/sh

# default values
vBgT=-1
vBgEv="sin"
compile=false
run_main=true
run_plots=false
make_mp4=true
parfile="parameters/parameters_example.dat"
vortfile="parameters/inVort_example.dat"
sim_folder="sim"

# help message
show_help() {
    echo "Usage: ./run.sh [--vBgT] [--vBgEv] [-c] [--no_main] [--parfile] [--vortfile] [-p] [-t] [-o] [-h]"
    echo
    echo "Options:"
    echo "  -c                          Compile the code before running"
    echo "  --vBgT <vBgT>               First argument passed to gintonic (optional)"
    echo "  --vBgEv <vBgEv>             Second argument passed to gintonic (optional between 'sin', 'cost', 'limit')"
    echo "  --no_main                   Do not run the main gintonic exe (only run post-process)"
    echo "  --parfile <path/to/file>    Path to the parameter file; default=parameters/parameters_example.dat"
    echo "  --vortfile <path/to/file>   Path to the vortex parameter file; default=parameters/inVort_example.dat"
    echo "  -p                          Run also plots.py"
    echo "  --keep_jpg                  Do not combine jpgs into mp4"
    echo "  -o <dir>                    Path to the simulation folder (creates the folder, cp exe and utils and run there); default=sim"
    echo "  -h                          Show this help message"
    exit 0
}

# parse arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        -c)
            compile=true
            shift
            ;;
        --vBgT)
            vBgT="$2"
            shift 2
            ;;
        --vBgEv)
            vBgEv="$2"
            shift 2
            ;;
        --parfile)
            parfile="$2"
            shift 2
            ;;
        --vortfile)
            vortfile="$2"
            shift 2
            ;;
        --no_main)
            run_main=false
            shift
            ;;
        -p)
            run_plots=true
            shift
            ;;
        --keep_jpg)
            make_mp4=false
            shift
            ;;
        -o)
            if [ -n "$2" ]; then
                sim_folder="$2"
                shift 2
            else
                echo "Error: -o requires a directory path argument."
                exit 1
            fi
            ;;
        -h)
            show_help
            ;;
        *)
            shift
            ;;
    esac
done

# save original working dir
og_pwd=$(pwd)

# create sim folder
if ! [ -d "$sim_folder" ]; then
    mkdir "$sim_folder"
fi
cp run.sh $sim_folder/run.sh

# compile gintonic and copy the exe to the sim folder
if [ "$run_main" = true ]; then
    if [ "$compile" = true ]; then
        if ! make; then
            echo
            echo "Compilation failed! Exiting."
            exit 1
        fi
    fi
    cp gin $sim_folder/
fi

# go to sim folder
cp $parfile $sim_folder/parameters.dat
cp $vortfile $sim_folder/inVort.dat
cd $sim_folder
echo "pwd: " $(pwd)

# run gin (main code)
if [ "$run_main" = true ]; then
    # create out folder
    out_folder="out"
    if [ -d "$out_folder" ]; then
        rm -r "$out_folder"
    fi
    mkdir "$out_folder"

    # copy parfiles
    cp $og_pwd/$parfile $out_folder/parameters.dat
    if [ -f "$og_pwd/$vortfile" ]; then
        cp $og_pwd/$vortfile $out_folder/inVort.dat
    fi

    # std output
    echo
    echo
    echo "Running gintonic in a Cartesian cell with parameter.dat as input file"
    echo

    # run with timer
    start=`date +%s`
    ./gin $vBgT $vBgEv
    end=`date +%s`
    runtime="$((end-start))"
    echo
    echo "Running time = $runtime s"
fi

# run plotting script
if [ "$run_plots" = true ]; then
    # create jpg folder
    jpg_folder="jpg"
    if [ -d "$jpg_folder" ]; then
        rm -r "$jpg_folder"
    fi
    mkdir "$jpg_folder"

    # run plotting script
    echo
    cp $og_pwd/plots.py .
    echo "Running plots.py"
    python3 plots.py -i $og_pwd/$sim_folder/"out" -o $og_pwd/$sim_folder/"jpg"
    echo

    # combining jpgs into mp4 and deleting jpgs
    if [ "$make_mp4" = true ]; then
        cd jpg
        ffmpeg -framerate 24 -pattern_type glob -i "phase_*.jpg" -c:v libx264 -pix_fmt yuv420p phase.mp4
        mv phase_ini.jpg phaseIni.jpg
        rm phase_*.jpg
        cd ..
    fi
fi

cd $og_pwd
echo "pwd: " $(pwd)
