__author__ = 'sebasanper'

def simpson_integrate(f, a, b, n):
    h = (b - a) / float(n)
    res = 0.0
    for i in range(n):
        res += abs(f(a + i * h) + 4.0 * f(a + i * h + h / 2.0) + f(a + (i + 1.0) * h))
    res *= h / 6.0
    return res

if __name__ == '__main__':
    from math import sin, pi
    def func(x):
        return 4.0

    print simpson_integrate(func, 0., 2.0, 500)
