#!/bin/bash
set -e
sed -i '/numpy/d; /openmdao/d' "$HOME/.config/pip/constraints.txt"  # Remove the pip constraint on the numpy and openmdao version
pip install .[testing,mphys]
