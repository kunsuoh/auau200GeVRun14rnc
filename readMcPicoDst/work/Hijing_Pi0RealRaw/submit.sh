#!/bin/bash
date
particle=pi0real
nrun=$1
path=/star/data01/pwg/kunsu/$particle/
tag=_wrongRow_0.4_0.2
star-submit-template -template submit.xml -entities nRuns=${nrun},base=${path},particle=${particle},tag=${tag}
#Fix privileges for report, log, etc.
find ./ -user $USER -exec chgrp rhstar {} \;
find ./ -user $USER -exec chmod g+rw {} \;
