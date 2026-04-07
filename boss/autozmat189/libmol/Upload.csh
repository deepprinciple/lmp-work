#!/bin/csh

if (${#argv} > 0) then
	set LIB_MOL_DIR = $1
endif

foreach file (local_include/*.h src/*.c)
	set dst = $LIB_MOL_DIR/$file
	diff $dst $file >& /dev/null
	if ($status == 0) continue;

	echo /bin/cp $file $dst
	/bin/cp $file $dst
end

