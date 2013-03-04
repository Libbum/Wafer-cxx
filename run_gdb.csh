#!/bin/csh -f

if ("$LAMRANK" == "0") then
  gdb $*
else
  $*
endif
exit 0
