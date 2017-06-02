import ValidDegree as vd
import DCMRevised as dcm_r
import PowerLawDistribution as pld


def generate_simple(a, alpha, beta, n, b = 2):
    #    1. Generate bi-degree-sequence according to Section 2.1, with F and G having finite variance.
    #    2. (Optional) Verify graphicality using Theorem 2.2.
    #    3. Randomly pair the in-degrees and out-degrees.
    #    4. If the resulting graph is not simple, repeat from step 3 (or from step
    #    1 if skipping step 2).

    flag = True
    while flag:
        fg = pld.PowerLaw(a, alpha, beta, b)
        bi_seq = vd.directed_gen(alpha, beta, fg, n)
        (model, flag) = dcm_r.directed_configuration_model_revised(bi_seq[0].tolist(), bi_seq[1].tolist())
    return model


a = 3
alpha = 4
beta = 3
n = 2000
b = 2

fg = pld.PowerLaw(a, alpha, beta)
(din, dout, original) = vd.directed_gen(alpha, beta, fg, n)
din = din.tolist()
dout = dout.tolist()
din_original = original[0].tolist()
dout_original = original[1].tolist()


#dcm = generate_simple(a, alpha, beta, n)

