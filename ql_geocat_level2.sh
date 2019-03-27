#!/usr/bin/env bash
#
# ql_geocat_level2.sh
#
# * DESCRIPTION: Wrapper script for ql_geocat_level2.py
#
# Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-04-01.
# Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
# Licensed under GNU GPLv3.

if [ -z "$CSPP_GEO_GEOCAT_HOME" ]; then
    echo "CSPP_GEO_GEOCAT_HOME must be set to the path where the CSPP software was installed."
    echo "export CSPP_GEO_GEOCAT_HOME=/home/me/GEOCAT/CSPP_GEOCAT/cspp-geo-geocat-x.xx"
    exit 1
fi

. ${CSPP_GEO_GEOCAT_HOME}/cspp_geo_geocat_runtime.sh

usage() {
    $PY_QL $CSPP_GEO_GEOCAT_HOME/scripts/ql_geocat_level2.py --help
}

if [ -z "$1" ]; then
    usage
    exit 3
fi

$PY_QL ${CSPP_GEO_GEOCAT_HOME}/scripts/ql_geocat_level2.py "$@"
