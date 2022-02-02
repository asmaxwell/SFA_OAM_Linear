#!/bin/sh
#g++ -c -Wall -Werror -fpic HpgForms.cpp
#g++ -shared -o HpgForms.so HpgForms.o
python3 setup.py build_ext --inplace
