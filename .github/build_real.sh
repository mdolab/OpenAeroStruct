#!/bin/bash
set -e
sed -i '/mphys/d' "$HOME/.config/pip/constraints.txt" # Remove the pip constraint on the mphys version
pip install .[testing,mphys]
