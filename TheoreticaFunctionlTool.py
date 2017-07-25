import math
import matplotlib.pyplot as plt
import PowerLawDistribution as pld


def corr(alpha, beta, E, d):
    s = beta/alpha
    c = E * (beta - 1) / beta

    a = c **((alpha - beta) / alpha) * (beta * (alpha - 1)) / (alpha * (beta - 1))

    b = E * (alpha - 1) / alpha


    E_W_1s = beta / (beta - 1 - s) * c ** (1 + s)
    Var_W_s = beta / (beta - 2 * s) * c ** (2 * s) - (beta * c ** s / (beta - s)) ** 2


    Cov_Dplus_Dminus = d * (a * E_W_1s - E ** 2)

    Var_Dplus = a ** 2 * (d ** 2 + (1 - d) ** 2) * Var_W_s + E
    Var_Dminus = beta / (beta - 2 ) * c ** 2 - E ** 2 + E


    corr = Cov_Dplus_Dminus / ((Var_Dplus) ** (1/2) * (Var_Dminus) ** (1/2))


    return corr


alpha = 3
beta = 3
E = 3


def test_corr(alpha, beta, E):
    d_lst = range(0,101,1)
    corr_lst = []

    for d in d_lst:
        d_i = d/100
        cor_d = corr(alpha, beta, E, d_i)
        corr_lst.append(cor_d)

    return corr_lst

