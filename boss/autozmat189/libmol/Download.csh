#!/bin/csh

if (${#argv} > 0) then
	set LIB_MOL_DIR = $1
endif

foreach file (local_include/*.h src/*.c)
	set src = $LIB_MOL_DIR/$file
	diff $src $file >& /dev/null
	if ($status == 0) continue;

	echo /bin/cp $src $file
	/bin/cp $src $file
end

