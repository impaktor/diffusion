#! /bin/bash

BASE_PATH=${HOME}"/crowding/write_up/correlation_data_fit"

SCRIPT_PATH=${HOME}"/crowding/src"
DATA_TARGET_PATH=${BASE_PATH}"/data_corr/seed1"
INPUT=${BASE_PATH}"/data_corr/input.dat"

# DATA_TARGET_PATH=${BASE_PATH}"/data_corr/small"
# INPUT=${BASE_PATH}"/data_corr/input_small.dat"

FILE_NAME_BASE=${DATA_TARGET_PATH}"/renormalization.dat"
PLOT_PATH=${BASE_PATH}"/fig"

# nice -19 ./prog -r $INPUT -w $FILE_NAME_BASE

echo "Summing diagonals"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix1"   > ${DATA_TARGET_PATH}/"diag_sum1"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix2"   > ${DATA_TARGET_PATH}/"diag_sum2"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix2b"  > ${DATA_TARGET_PATH}/"diag_sum2b"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix3"   > ${DATA_TARGET_PATH}/"diag_sum3"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix3b"  > ${DATA_TARGET_PATH}/"diag_sum3b"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix3c"   > ${DATA_TARGET_PATH}/"diag_sum3c"

${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix1"   10 > ${DATA_TARGET_PATH}/"diag_sum1_small"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix2"   10 > ${DATA_TARGET_PATH}/"diag_sum2_small"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix2b"  10 > ${DATA_TARGET_PATH}/"diag_sum2b_small"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix3"   10 > ${DATA_TARGET_PATH}/"diag_sum3_small"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix3b"  10 > ${DATA_TARGET_PATH}/"diag_sum3b_small"
${SCRIPT_PATH}/diag_sum.py ${FILE_NAME_BASE}"_matrix3c"   10 > ${DATA_TARGET_PATH}/"diag_sum3c_small"

echo "Processing matrices1"
nice -19 ${SCRIPT_PATH}/lomholt.py $FILE_NAME_BASE ${FILE_NAME_BASE}"_matrix1"   > $DATA_TARGET_PATH/converge_method1.dat
echo "Processing matrices2"
nice -19 ${SCRIPT_PATH}/lomholt.py $FILE_NAME_BASE ${FILE_NAME_BASE}"_matrix2"   > $DATA_TARGET_PATH/converge_method2.dat
echo "Processing matrices2b"
nice -19 ${SCRIPT_PATH}/lomholt.py $FILE_NAME_BASE ${FILE_NAME_BASE}"_matrix2b" > $DATA_TARGET_PATH/converge_method2b.dat
echo "Processing matrices3"
nice -19 ${SCRIPT_PATH}/lomholt.py $FILE_NAME_BASE ${FILE_NAME_BASE}"_matrix3"   > $DATA_TARGET_PATH/converge_method3.dat
echo "Processing matrices3b"
nice -19 ${SCRIPT_PATH}/lomholt.py $FILE_NAME_BASE ${FILE_NAME_BASE}"_matrix3b" > $DATA_TARGET_PATH/converge_method3b.dat
echo "Processing matrices3c"
nice -19 ${SCRIPT_PATH}/lomholt.py $FILE_NAME_BASE ${FILE_NAME_BASE}"_matrix3c" > $DATA_TARGET_PATH/converge_method3c.dat

# echo "Generating plots"
# cd $PLOT_PATH; for f in *gnu; do gnuplot $f; done
