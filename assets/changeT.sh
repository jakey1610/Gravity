#!/bin/sh
export OMP_NUM_THREADS=1;
python3 gridGen6.py
export OMP_NUM_THREADS=2;
python3 gridGen6.py
export OMP_NUM_THREADS=3;
python3 gridGen6.py
export OMP_NUM_THREADS=4;
python3 gridGen6.py
export OMP_NUM_THREADS=5;
python3 gridGen6.py
export OMP_NUM_THREADS=6;
python3 gridGen6.py
export OMP_NUM_THREADS=7;
python3 gridGen6.py
export OMP_NUM_THREADS=8;
python3 gridGen6.py
