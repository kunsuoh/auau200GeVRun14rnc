#!/bin/bash
date
particle=pi0Dalitz
nrun=$1
path=/star/data01/pwg/kunsu/$particle/
rm -rf LocalLibraries_$particle*
star-submit-template -template submit.xml -entities nRuns=${nrun},base=${path},particle=${particle}
#Fix privileges for report, log, etc.
find ./ -user $USER -exec chgrp rhstar {} \;
find ./ -user $USER -exec chmod g+rw {} \;
