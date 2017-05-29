import random


def directed_gen(alpha, beta, fg, n):
    """
    a function that generates same-sum sample pair of bi-degree using the algorithm
    the function returns a tuple  (in-degree, out-degree)

    Input Parameters
    -----------------
    alpha: the power for W+ 
    beta: the power for W-
    fg: the distribution of bi-degree 
    n: simulation size
    """

    # derive kappa
    kappa = min(1 - 1/alpha, 1 - 1/beta, 1/2)
    delta_0 = (1 - 1 / n) * kappa  # take a fixed delta0
    tol = n ** (1 - kappa + delta_0) # take the tolerance limit for delta_n

    # derive the sample bi-degree sequence
    bi_seq = fg.iid(n)
    in_seq = bi_seq[0] # the sample in-degree sequence
    out_seq = bi_seq[1] # the sample out-degree sequence

    # derive the Delta_n
    in_sum = sum(in_seq)
    out_sum = sum(out_seq)
    delta_n = in_sum - out_sum

    # repeat sample generation until Delta_n is small enough to be "negligible"
    while abs(delta_n) > tol:
        # repeat the sample working
        bi_seq = fg.iid(n)
        in_seq = bi_seq[0]
        out_seq = bi_seq[1]

        in_sum = sum(in_seq)
        out_sum = sum(out_seq)
        delta_n = in_sum - out_sum

    if abs(delta_n) > 0:
        # derive random sample node's index  delta_i
        delta_i = random.sample(range(0, n), int(abs(delta_n)))
        for i in delta_i:
            if delta_n > 0: # in-degree sequence sum is larger , increase the random out-degree
                out_seq[i] += 1
            else: # out-degree is larger, increase the random in-degree
                in_seq[i] += 1

    return in_seq, out_seq
