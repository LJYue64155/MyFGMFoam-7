#!/bin/bash
. /OpenFOAM/of7/OpenFOAM-7/etc/bashrc WM_NCOMPROCS=2;
/usr/bin/gdb "$@"