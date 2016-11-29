For now, model is independent from CATS.

To compile:
make -j 10

Test run with 400 time steps and 5 fishes:
rm aa; ./model -s 42  -n 100 -F 5 -i ../arena/arena2roomsLargeSmall.png > aa ; ./visualizer.py aa 5 ../arena/arena2roomsLargeSmall.png

If the display is slow on your computer, try with a smaller arena file (e.g with resolution 100x100).


