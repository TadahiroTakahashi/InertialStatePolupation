
The purpose of this repository is to calculate the population ratio within the atomic internal states. I plan to create the following two codes:

(1) Load atomic transition rates from the NIST database and generate a matrix in Excel.
(2) Calculate the population ratio within the atomic internal states using the matrix described in Excel.
At this time, only (2) will be uploaded.

When running the function in pumpingRatio_calc.py from (2), please ensure that "grotorian.xlsx" and the /dataSet/ folder are in the same directory. Structure the files in the /dataSet/ folder as illustrated in "0_test", which means:

(1) Create a directory (you can name it as you wish).
(2) Inside that directory, prepare paramSet.xlsx (please use the provided format by copying it).
With paramSet.xlsx, you can set various parameters. It also allows you to add more lasers.

When executing pumpingRatio_calc.py, input the name of the directory you created as an argument.as bellow:

> python pumpingRatio_calc.py 0_test

Thank you!
