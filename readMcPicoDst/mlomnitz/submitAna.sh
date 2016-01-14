#!/bin/csh
set tops=$2
set i=$1
while ( $i <= $tops )
    writeJob.sh $i
    @ i = $i + 1
end
