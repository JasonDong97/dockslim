source /usr/local/gromacs/bin/GMXRC
service ssh start
service shellinabox start

python dockslim/docking.py mq &
python mdslim/md.py &
python ecal/ecal_demo.py &

sleep infinity