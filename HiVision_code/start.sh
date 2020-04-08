#!/bin/bash
indexpath="./indexes/"
shppath="./datasets/alaska_OSM/"
patternpath="./patterns/"
process=4
redishost="127.0.0.1"
redisport=6379
serviceport=10080
nohup mpirun -np $process ./build/hivision_engine $indexpath $redishost $redisport >hivision_engine.log 2>&1 &
nohup ./build/hivision_server $shppath $indexpath $patternpath $redishost $redisport $serviceport > ./hivision_server.log 2>&1 &
