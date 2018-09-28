#!/bin/bash
make
cp ./bin/namics ./inputs/
cd inputs/
./namics tag5
#valgrind --track-origins=yes ./namics tag_x52.in
