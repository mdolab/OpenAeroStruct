#!/bin/bash
set -e

testflo -n 2 -v . --coverage --coverpkg openaerostruct
