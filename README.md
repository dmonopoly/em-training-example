## What is this?
Implementations of brute force EM, forward-backward EM, and Viterbi to find
the best matching tag sequence for a simple observation of sequence.
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

### Examples
Forward-backward and Viterbi:

    ./two_tag_example

Brute force (pass any param):

    ./two_tag_example -b

