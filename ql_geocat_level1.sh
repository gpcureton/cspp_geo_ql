#!/bin/bash
# Wrapper script for GVAR components from SSEC
#
# Copyright 2011, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_GEO_GEOCAT_HOME" ]; then
    echo "CSPP_GEO_GEOCAT_HOME must be set to the path where the CSPP software was installed."
    echo "export CSPP_GEO_GEOCAT_HOME=/home/me/GEOCAT/CSPP_GEOCAT/cspp-geo-geocat-x.xx"
    exit 1
fi

. ${CSPP_GEO_GEOCAT_HOME}/geocat_runtime.sh

usage() {
    $PY $CSPP_GEO_GEOCAT_HOME/l2/ql_geocat_level1.py --help
}

if [ -z "$1" ]; then
    usage
    exit 3
fi

$PY ${CSPP_GEO_GEOCAT_HOME}/l2/ql_geocat_level1.py "$@"

