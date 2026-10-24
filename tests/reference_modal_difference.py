"""Independent 100-digit recurrence reference; uses only Python's standard library.

Print CSV fixtures for check_higher_orders.R. Evaluate the ordinary exact and
backward-Euler formulas at high precision, then subtract them. This deliberately
does not reproduce the cancellation-free identity used by the production code.
Inputs are converted from their binary64 values exactly.
"""
import csv
from decimal import Decimal as D, localcontext
import math
import sys


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def mul(a, b):
    return a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]


def div(a, b):
    denominator = b[0]*b[0] + b[1]*b[1]
    return ((a[0]*b[0] + a[1]*b[1])/denominator,
            (a[1]*b[0] - a[0]*b[1])/denominator)


def exp_i(x):
    cosine, sine = D(1), x
    ct, st = cosine, sine
    for n in range(1, 200):
        ct *= -x*x / D((2*n-1)*(2*n))
        st *= -x*x / D((2*n)*(2*n+1))
        cosine += ct
        sine += st
        if max(abs(ct), abs(st)) < D('1e-110'):
            return cosine, sine
    raise ArithmeticError('Taylor series did not converge')


def reference(values):
    k, fraction, dt, rate, numerical_rate, initial, forcing, amplitude, omega = (
        D.from_float(float(v)) for v in values)
    time = (k-1+fraction)*dt
    theta = omega*dt
    q = 1/(1+dt*numerical_rate)
    qk = (k*q.ln()).exp()
    decay = (-rate*time).exp()
    exact = initial*decay
    numerical = initial*qk
    if omega and forcing and amplitude:
        exact += forcing*amplitude*div(sub(exp_i(omega*time), (decay, D(0))),
                                       (rate, omega))[1]
        sinc = exp_i(theta/2)[1]/(theta/2)
        response = mul(exp_i(theta/2),
                       div(sub(exp_i(theta*k), (qk, D(0))),
                           sub(exp_i(theta), (q, D(0)))))[1]
        numerical += forcing*amplitude*q*dt*sinc*response
    return exact-numerical


writer = csv.writer(sys.stdout)
writer.writerow(['k', 'fraction', 'time_step', 'exact_rate', 'numerical_rate',
                 'initial', 'forcing', 'A', 'omega', 'difference'])
with localcontext() as context:
    context.prec = 100
    for dt, k in [(0.01, 31), (2.25e-16, 2e15)]:
        for fraction in [0.0, 0.37, 1.0]:
            for initial, forcing in [(10, 0), (0, 10), (7, -3)]:
                for relative_change in [0, 1e-10]:
                    row = [k, fraction, dt, 0.00025,
                           0.00025*(1+relative_change), initial, forcing, 0.7, math.pi]
                    writer.writerow(row + [str(reference(row))])
