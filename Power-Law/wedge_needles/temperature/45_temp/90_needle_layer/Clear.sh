#!/bin/sh
#. $WM_PROJECT_DIR/bin/tools/CleanFunctions
#cleanTimeDirectories
foamListTimes -rm
rm -r processor* 2>/dev/null
rm error.log 2>/dev/null
rm plot_residuals.png 2>/dev/null
rm sorted_residuals.csv 2>/dev/null