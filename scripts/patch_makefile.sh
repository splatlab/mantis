#!/bin/bash
MAKEFILE=$1
LIBPATHS="$(./libpaths.sh)"
LIBPATHS="$(echo "$LIBPATHS" | sed '$ ! s/$/\\/')"    # escape new lines
sed -i "/TARGETS=/aMYLIBS=$LIBPATHS" $MAKEFILE        # add MYLIBS var with all libraries path
sed -i "s/LDFLAGS += /LDFLAGS += \$\(MYLIBS\) /g" $MAKEFILE  # add MYLIBS to LDFLAGS
