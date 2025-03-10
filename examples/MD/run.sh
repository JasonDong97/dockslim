#!/bin/bash
#source /usr/local/gromacs/bin/GMXRC
MAIN_DIR='/home/gromacs-2024.5/md/MD'


start_time=$(date +%s)

for SUBDIR in "$MAIN_DIR"/*/ ; do
    if [[ -f "$SUBDIR/md.tpr" ]]; then
        echo "Processing directory: $SUBDIR"
	cd "$SUBDIR" || continue
	/usr/local/gromacs/bin/gmx mdrun -v -deffnm md -maxh 0.01667
	cd "$MAIN_DIR"
    else
        echo "Skipping $SUBDIR as it does not contain necessary files."
    fi
done


end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Total elapsed time: $elapsed_time seconds"

hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total elapsed time: ${hours}h ${minutes}m ${seconds}s"



