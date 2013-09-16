## What is this?
A brute-force implementation of Expectation-Maximation (EM) training.
Based on Kevin Knight's classnotes for CSCI 662 (Advanced NLP), p73:
http://d.pr/i/ADyl

### Purpose?
Research, understanding the algorithms, etc. So it's a bit hack-y.

### Code Organization
Main file: TwoTagExample.cc.

## Running
Compile using cmake. Example:

    mkdir build && cd build
    cmake ..
    make

Run the generated executable. Alternatively run any test executables, but make
sure you do `cmake -Dtest=ON ..` first.

## TODO
Implement the Forward-Backward algorithm - the more efficient version.
Then do the Viterbi algorithm, or maybe Dijkstra's on the -log of edges.

