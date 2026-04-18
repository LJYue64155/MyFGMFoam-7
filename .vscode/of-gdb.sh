#!/bin/bash
# 获取传入的完整命令（包括求解器路径和参数）
# 废了，不能调试！！
CMD="$@"
/usr/bin/gdb "$@"

. /home/liujiayue/OpenFOAM/of7/OpenFOAM-7/etc/bashrc WM_NCOMPROCS=2;
. /home/liujiayue/OpenFOAM/of7/liujiayue-7/applications/MyFGMFoam/src/bashrc;

exec $CMD