#!/bin/bash
find src/ -depth -name "*.cu" -exec sh -c 'mv "$1" "${1%.cu}.cpp"' _ {} \;