#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# Exit if stash cache isn't mounted.

UBOONE_EXAMPLE_DATA_DIR=/cvmfs/uboone.osgstorage.org/stash/uboone_example_overlay
if [ ! -d $UBOONE_EXAMPLE_DATA_DIR ]; then
  echo "Quittig because stash cache isn't available."
  exit
fi

# This script runs the full mc+reco chain using standard released fcl files.

input=$UBOONE_EXAMPLE_DATA_DIR/pnfs/uboone/data/uboone/raw/swizzle_trigger_streams/mergeext_unbiased/prod_v04_26_04_05/00/00/83/07/PhysicsRun-2016_10_15_0_6_5-0008307-00029_20161015T141511_ext_unbiased_20161020T160531_merged.root
for fcl in standard_overlay_gen_driver.fcl standard_g4_uboone.fcl standard_detsim_nonoise_overlay_uboone.fcl standard_overlay_uboone.fcl reco_uboone_mcc8_driver_overlay_stage1a.fcl reco_uboone_mcc8_driver_overlay_stage1b.fcl reco_uboone_mcc8_driver_overlay_stage1c.fcl reco_uboone_mcc8_driver_overlay_stage2.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 5"
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi
  input=$output
done
