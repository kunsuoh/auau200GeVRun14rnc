#!/bin/bash
rm -rf LocalLibraries.* 
date
nrun=$1
path=/star/data01/pwg/kunsu/Pi0Dalitz/
star-submit-template -template submit.xml -entities nRuns=${nrun},base=${path}
#Fix privileges for report, log, etc.
find ./ -user $USER -exec chgrp rhstar {} \;
find ./ -user $USER -exec chmod g+rw {} \;
