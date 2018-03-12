#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# This script runs the standard swizzle fcl file.

input=$UBOONE_EXAMPLE_DATA_DIR/ubdaq/PhysicsRun-2018_2_26_13_49_57-0015224-00007.ubdaq
fcl=swizzle.fcl
out=`basename $fcl .fcl`.out
err=`basename $fcl .fcl`.err
cmd="lar --rethrow-all -c $fcl -s $input -n 5"
echo $cmd
$cmd > $out 2> $err
stat=$?
echo "Command finished with status $stat"
if [ $stat -ne 0 ]; then
  exit $stat
fi
