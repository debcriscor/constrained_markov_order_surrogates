# Contrained Markov Order Surrogates

This is a surrogate algorithm based on swapping that provides an exact test for Markov order. 

The surrogates are constrained realisations which preserve Markov properties of a given order. 

Example

Create surrogate that exactly preserves second-order Markov properties of
the original time series:

>> ts = [1 2 1 1 2 1 2 1 1 2];
>> surrogate = generate_constrained_surrogate(ts,2);

Relevant paper:

Corrêa, D. C., Moore, J. M., Jüngling, T., & Small, M. (2020). Constrained Markov order surrogates. Physica D: Nonlinear Phenomena, Vol 406, 132437.

https://doi.org/10.1016/j.physd.2020.132437