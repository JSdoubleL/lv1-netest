#!/bin/bash

dir=$1
for i in $1/polytomy*; do 
	./paup4a168_centos64 $i
done

