#!/bin/bash

#root -l "RunSEAview.c(1.0, true, false, false, true, false, 1.0, 0, 100, 0.9, 100, true, 0.5, 2.0)"

rm *o
g++ -g RunSEAview.c plothelper.c  RecoOpAng1.c SEAview_lite_epem.c  SetRealAspectRatio2.c -c $(root-config --libs --cflags) -v
g++ -g *.o $(root-config --libs --cflags) -o RunSEA -v

#valgrind --leak-check=full --show-leak-kinds=all --log-file=mylog ./RunSEA


