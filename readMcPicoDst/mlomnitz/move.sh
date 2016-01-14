#!/bin/csh
set tops=$2
set i=$1
while ( $i <= $tops )
    mv KFHE_Out_wPXLIST/testSIMU_$i.root KF_21_11_2014/.
    @ i = $i + 1
end
