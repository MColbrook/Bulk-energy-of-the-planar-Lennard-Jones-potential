#!/usr/bin/env python3
"""Exact rational enclosures for the triangular Lennard-Jones lattice constants.

All mathematical arithmetic uses fractions.Fraction and integers. Decimal is
used only to display bounds with outward rounding. No special functions or
numerical quadrature are used.

For integer s>1, rational a>0, z=N+a, and m>=1, Euler--Maclaurin gives

  zeta(s,a) = sum_{n=0}^{N-1}(n+a)^(-s) + z^(1-s)/(s-1)
             + z^(-s)/2
             + sum_{k=1}^m B_(2k)/(2k)! (s)_(2k-1) z^(1-s-2k) + R.

Here

  |R| <= C_m (s)_(2m-1)/(2m)! z^(1-s-2m),
  C_m = sum_{j=0}^{2m} binomial(2m,j)*abs(B_j).

Proof of the bound: write B_n(t)=sum_j binomial(n,j) B_j t^(n-j).
For 0<=t<=1 this gives |B_(2m)(t)|<=C_m by the triangle inequality.
Integration by parts on successive unit intervals, using
  B_n'(t)=n B_(n-1)(t), B_n(1)=B_n(0) for n>=2,
gives the exact Euler--Maclaurin remainder
  R = -(2m)!^(-1) integral_N^infinity B_(2m)({x}) f^(2m)(x) dx,
where f(x)=(x+a)^(-s). Since f^(2m)=(s)_(2m)(x+a)^(-s-2m),
the claimed bound follows by integrating this positive power function.
Every quantity in the finite expression and remainder bound is rational.

The standard triangular Epstein zeta identity is
  Z_(2s) = 6*zeta(s)*L(s,chi_-3),
  L(s,chi_-3) = 3^(-s)*(zeta(s,1/3)-zeta(s,2/3)).
Then
  a_star = (Z12/Z6)^(1/6), e_triangle = -Z6^2/(2*Z12).
The sixth-root enclosure is proved by integer-power comparisons.
"""

if not __debug__:
    raise RuntimeError("Certificate checks require Python assertions; do not use -O or -OO.")

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import argparse
import json


def bernoulli_numbers(n):
    b = [F(1)]
    for k in range(1, n + 1):
        b.append(-sum(F(comb(k + 1, j)) * b[j] for j in range(k)) / (k + 1))
    return b


def rising(s, k):
    result = 1
    for j in range(k):
        result *= s + j
    return result


def hurwitz_interval(s, a, N=32, m=8):
    assert isinstance(s, int) and s > 1 and a > 0 and N >= 1 and m >= 1
    b = bernoulli_numbers(2 * m)
    z = F(N) + a
    center = sum((F(n) + a) ** (-s) for n in range(N))
    center += z ** (1 - s) / (s - 1) + z ** (-s) / 2
    for k in range(1, m + 1):
        center += b[2 * k] * F(rising(s, 2 * k - 1), factorial(2 * k)) * z ** (1 - s - 2 * k)
    C_m = sum(abs(F(comb(2 * m, j)) * b[j]) for j in range(2 * m + 1))
    error = C_m * F(rising(s, 2 * m - 1), factorial(2 * m)) * z ** (1 - s - 2 * m)
    return (center - error, center + error), error, C_m


def mul_positive(x, y):
    assert 0 < x[0] <= x[1] and 0 < y[0] <= y[1]
    return x[0] * y[0], x[1] * y[1]


def divide_positive(x, y):
    assert 0 < x[0] <= x[1] and 0 < y[0] <= y[1]
    return x[0] / y[1], x[1] / y[0]


def integer_root_floor(value, n):
    assert isinstance(value, int) and value >= 0 and n >= 1
    if value == 0:
        return 0
    lo = 0
    hi = 1 << ((value.bit_length() + n - 1) // n)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** n <= value:
            lo = mid
        else:
            hi = mid
    assert lo ** n <= value < (lo + 1) ** n
    return lo


def root_interval(x, n=6, digits=24):
    assert 0 < x[0] <= x[1]
    scale = 10 ** digits
    lower = integer_root_floor((x[0].numerator * scale ** n) // x[0].denominator, n)
    upper = integer_root_floor((x[1].numerator * scale ** n) // x[1].denominator, n) + 1
    lo, hi = F(lower, scale), F(upper, scale)
    assert lo ** n <= x[0] and hi ** n >= x[1]
    return lo, hi


def decimal_outward(value, lower, digits=18):
    # Both the division and final quantization round in the required direction.
    mode = ROUND_FLOOR if lower else ROUND_CEILING
    with localcontext() as context:
        context.prec = max(80, digits + 20)
        context.rounding = mode
        number = Decimal(value.numerator) / Decimal(value.denominator)
        return format(number.quantize(Decimal(1).scaleb(-digits), rounding=mode), 'f')


def interval_record(x, digits=18):
    return {
        'lower_exact': str(x[0]),
        'upper_exact': str(x[1]),
        'lower_decimal_outward': decimal_outward(x[0], True, digits),
        'upper_decimal_outward': decimal_outward(x[1], False, digits),
        'width_exact': str(x[1] - x[0]),
    }


def compute(N=32, m=8, certificate=None):
    hurwitz = {}
    z = {}
    C_m = None
    for s in (3, 6):
        for a in (F(1), F(1, 3), F(2, 3)):
            interval, error, C_m = hurwitz_interval(s, a, N, m)
            z[s, a] = interval
            hurwitz[f'{s},{a}'] = {
                **interval_record(interval),
                'remainder_bound_exact': str(error),
            }
    Z = {}
    L = {}
    for s in (3, 6):
        x, y = z[s, F(1, 3)], z[s, F(2, 3)]
        L[s] = ((x[0] - y[1]) / 3 ** s, (x[1] - y[0]) / 3 ** s)
        product = mul_positive(z[s, F(1)], L[s])
        Z[2 * s] = 6 * product[0], 6 * product[1]
    ratio = divide_positive(Z[12], Z[6])
    a_star = root_interval(ratio)
    energy_magnitude = divide_positive(mul_positive(Z[6], Z[6]),
                                       (2 * Z[12][0], 2 * Z[12][1]))
    energy = -energy_magnitude[1], -energy_magnitude[0]
    result = {
        'method': 'Exact rational Euler--Maclaurin with explicit Bernoulli polynomial remainder bound',
        'N': N, 'm': m,
        'bernoulli_polynomial_sup_bound_exact': str(C_m),
        'hurwitz_zeta_intervals': hurwitz,
        'L_chi_minus3_at_3': interval_record(L[3]),
        'L_chi_minus3_at_6': interval_record(L[6]),
        'Z6': interval_record(Z[6]),
        'Z12': interval_record(Z[12]),
        'a_star': interval_record(a_star),
        'e_triangle': interval_record(energy),
        'minus_e_triangle': interval_record(energy_magnitude),
    }
    if certificate is not None:
        data = json.loads(Path(certificate).read_text(encoding='utf-8'))
        bound = F(data['bound'])
        gap = bound - energy_magnitude[1], bound - energy_magnitude[0]
        result['certificate_bound_exact'] = str(bound)
        result['certificate_gap_to_triangle'] = interval_record(gap)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='triangular_energy_intervals.json')
    parser.add_argument('--certificate')
    args = parser.parse_args()
    result = compute(certificate=args.certificate)
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    for key in ('Z6', 'Z12', 'a_star', 'e_triangle', 'certificate_gap_to_triangle'):
        if key in result:
            interval = result[key]
            print(f"{key}: [{interval['lower_decimal_outward']}, {interval['upper_decimal_outward']}]")


if __name__ == '__main__':
    main()
