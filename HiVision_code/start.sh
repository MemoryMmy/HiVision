indexpath="./indexes/"
shppath="./datasets/alaska_OSM/"
patternpath="./patterns/"
process=4
redishost="127.0.0.1"
redisport=6379
serviceport=10080
nohup mpirun -np $process ./Vision $indexpath $redishost $redisport >mpi.log 2>&1 &
nohup ./visioncrowserver/Crowserver $shppath $indexpath $patternpath $redishost $redisport $serviceport > ./crowserver.log 2>&1 &
