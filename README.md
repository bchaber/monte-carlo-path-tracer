# Monte Carlo Path Tracer
A simple Monte Carlo path tracer written in Julia.

# Motivation
My way of learning new things is to find a minimalistic examples of a method, rewrite it to my taste in order to get some insight how and why something works. To understand how Monte Carlo Path Tracer works I've found an amazing code by Keving Beason (https://www.kevinbeason.com/smallpt/). It is written in plain-C++.

The first attempt to rewrite it to Julia will be recorded in `unoptimized.jl`. It uses vectors. However, the performance was much worse than the corresponding C++ code. The second version uses tuples instead of vectors, as they are stack allocated and compiler can reason about their size and unroll some loops. The idea about using tuples comes from `parenthephobia` at https://www.reddit.com/r/Julia/comments/adkr9i/ray_tracer_needs_optimising/.
