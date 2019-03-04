#!/bin/bash
find src/mesodyn -depth -name "*.cpp" -exec sh -c 'mv "$1" "${1%.cpp}.cu"' _ {} \;
find src -name "tools.cpp" -exec sh -c 'mv "$1" "${1%.cpp}.cu"' _ {} \;
find src/ -name "mesodyn.cpp" -exec sh -c 'mv "$1" "${1%.cpp}.cu"' _ {} \;