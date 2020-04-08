#redis-cli flushall
ps -ef | grep ./Vision | grep -v grep| awk '{print "kill -9 " $2}'| sh
ps -ef | grep ./visioncrowserver/Crowserver | grep -v grep| awk '{print "kill -9 " $2}'| sh
