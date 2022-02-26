#!/usr/bin/env sh

if [ $# -eq 1 ] && ([ "$1"="Debug" ] || [ "$1"="Release" ]); then
    build_type=$1
else
    build_type="Release"
fi

cmake -S . -B build -DCMAKE_BUILD_TYPE=$build_type\
    -DCMAKE_VERBOSE_MAKEFILE=ON && cmake --build build\
    && mv build/compile_commands.json .

