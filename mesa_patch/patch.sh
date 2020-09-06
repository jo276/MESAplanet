#!/bin/bash/

# script to patch v6208 of MESA to run planetary evolution
# it does two things 1) change the calculation of KH to a proper integral
# 2) puts a timestep criterion on the envelope mass fraction not total mass

\cp *.f $MESA_DIR/star/private/.
\cp history_columns.list $MESA_DIR/star/defaults/.