#!/bin/bash
#redis-cli flushall
ps -ef | grep hivision_engine | grep -v grep| awk '{print "kill -9 " $2}'| sh
ps -ef | grep hivision_server | grep -v grep| awk '{print "kill -9 " $2}'| sh

