from scipy.optimize import golden
from tabulate import tabulate
from math import sqrt
from sympy import *

x, y = symbols('x,y')
f = 4 * pow(x, 2) + 3 * pow(y, 2) - 4 * x * y + x

def getAgrad(x1, x2):
    return [agrad(diff(f, x), x1, x2), agrad(diff(f, y), x1, x2)]

def agrad(x, x1, x2):
    return -1 * eval(str(x).replace(" ", "").replace("x", str(x1)).replace("y", str(x2)))

def processing(v, eps, x1, x2):  # v - выбор метода : 1 метод наискорейшего спуска, 2 метод сопреженных гардиентов
    k, b1, b2 = 0, 0, 0
    agradx = getAgrad(x1, x2)
    table = [['N', 'x(n)', 'S(n)'], [k, f'({x1};{x2})', f'({agradx[0]};{agradx[1]})']]
    while sqrt(pow(agradx[0], 2) + pow(agradx[1], 2)) >= eps:
        min, k = golden(lambda z: lambdify((x, y), f)(x1 + agradx[0] * z, x2 + agradx[1] * z)), k + 1
        x1, x2 = x1 + agradx[0] * min, x2 + agradx[1] * min
        if v == 1:
            agradx = getAgrad(x1, x2)
        else:
            b = pow(pow(getAgrad(x1,x2)[0],2) + pow(getAgrad(x1,x2)[1],2), 2) / pow(pow(agradx[0],2) + pow(agradx[1],2), 2)
            agradx[0], agradx[1] = getAgrad(x1, x2)[0] + b * agradx[0], getAgrad(x1, x2)[1] + b * agradx[1]
        table.append([k, f'({x1};{x2})', f'({agradx[0]};{agradx[1]})'])
    print(tabulate(table, headers='keys') + '\n')

processing(1, 0.0000001, 0, 0), processing(2, 0.0001, 0, 0)
