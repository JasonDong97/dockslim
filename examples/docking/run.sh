#!/bin/bash
WORK_DIR="$(dirname "$(readlink -f "$0")")"
cd $WORK_DIR

# 设置脚本目录
SCRIPTS_HOME='../../dockslim'

python $SCRIPTS_HOME/docking.py cli -c 2F0Z1.pdb -l ligands/ -o outdir
