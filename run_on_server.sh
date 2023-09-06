## RSYNC downward

USER=quimortiz
SERVER=hal-9000.lis.tu-berlin.de


echo user:$USER
echo server:$SERVER


# time rsync -avrL  \
#   --exclude '*.pdf' \
#   --exclude '*.dat' \
#   --exclude '*.dot' \
#   --exclude '*.pddl' \
#   --exclude '*.o' \
#   --exclude '*.a' \
#   --exclude '*.png' \
#   --exclude '*.so' \
#   --exclude 'x.exe' \
#   --exclude '.git/' \
#   --exclude '.png' \
#   --exclude '01-joint-vs-conditional/tmp' \
#   /home/quim/stg/pddlss_gnlp/ $USER@$SERVER:~/stg/pddlss_gnlp
#
#
#
# cmd="python3 run_experiments_thesis.py"
#
# ssh  $USER@$SERVER  << EOF
#   cd ~/stg/pddlss_gnlp/01-joint-vs-conditional
#   source ~/kino/bin/activate.fish
#   make  -j20
#   $cmd | tee from_laptop.log
# EOF


# $USER@$SERVER:~/stg/pddlss_gnlp/01-joint-vs-conditional/results \
time rsync -avr  \
  --exclude '*.dat' \
  --exclude '*.dot' \
  $USER@$SERVER:~/stg/pddlss_gnlp/01-joint-vs-conditional/results_plots \
  $USER@$SERVER:~/stg/pddlss_gnlp/01-joint-vs-conditional/results \
  /home/quim/stg/pddlss_gnlp/01-joint-vs-conditional





