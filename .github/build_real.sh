#!/bin/bash
set -e
sed -i '/openmdao/d' "$HOME/.config/pip/constraints.txt" # Remove the pip constraint on the openmdao version
pip install .[testing,mphys]
