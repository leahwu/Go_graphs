import random


def directed_gen(alpha, beta, fg, n):
    """
    a function that generates same-sum sample pair of bi-degree using the algorithm
    the function returns a tuple  (in-degree, out-degree)

    Input Parameters
    -----------------
    alpha : 
    beta: 
    n: simulation graph size
    """

    # derive kappa
    kappa = min(1 - 1/alpha, 1 - 1/beta, 1/2)
    delta_0 = (1 - 1 / n) * kappa  # kappa fixed
    tol = n ** (1 - kappa + delta_0)

    bi_seq = fg.iid(n)
    in_seq = bi_seq[0]
    out_seq = bi_seq[1]
    bi_seq = 0
    in_sum = sum(in_seq)
    out_sum = sum(out_seq)
    delta_n = in_sum - out_sum

    # repeat until Delta_n is small enough to be "negligible"
    while abs(delta_n) > tol:
        bi_seq = fg.iid(n)
        in_seq = bi_seq[0]
        out_seq = bi_seq[1]
        bi_seq = 0
        in_sum = sum(in_seq)
        out_sum = sum(out_seq)
        delta_n = in_sum - out_sum

    if abs(delta_n) > 0:
        # derive random sample nodes deltai
        delta_i = random.sample(range(0, n), int(abs(delta_n)))
        for i in delta_i:
            if delta_n > 0:
                out_seq[i] += 1
            else:
                in_seq[i] += 1

    return in_seq, out_seq
