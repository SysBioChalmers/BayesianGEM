#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/21/20

"""ode_model.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from sklearn.metrics import mean_squared_error as MSE


# from tqdm import tqdm


def dydt_func(t, y, theta):
    '''
    y:
        y[0]: m1,
        y[1]: p1,
        y[2]: m2,
        y[3]: p2,
        y[4]: m3,
        y[5]: p3
    theta (parameters):
        theta[0]: alpha0: a0
        theta[1]: n: n
        theta[2]: beta: b,
        theta[3]: alpha: a,
    in math:
        dm1/dt = -m1 + a/(1+p3**n) + a0
        dp1/dt = -b*(p1-m1)
        dm2/dt = -m2 + a/(1+p1**n) + a0
        dp2/dt = -b*(p2-m2)
        dm3/dt = -m3 + a/(1+p2**n) + a0
        dp2/dt = -b*(p2-m2)

    '''
    m1, p1, m2, p2, m3, p3 = y

    a0, n, b, a = theta

    dydt = np.zeros(6)
    dydt[0] = -m1 + a / (1 + p3 ** n) + a0
    dydt[1] = -b * (p1 - m1)
    dydt[2] = -m2 + a / (1 + p1 ** n) + a0
    dydt[3] = -b * (p2 - m2)
    dydt[4] = -m3 + a / (1 + p2 ** n) + a0
    dydt[5] = -b * (p3 - m3)
    #
    # if dydt_func.pbar is not None:
    #     dydt_func.pbar.update(1)
    #     dydt_func.pbar.set_description('t = {:.3f}'.format(t))

    return dydt


def do_solve_ivp(t_span, y0, theta):
    # dydt_func.pbar = None

    # with tqdm() as pbar:
    #     dydt_func.pbar = pbar
    sol = solve_ivp(
        fun=lambda t, y: dydt_func(t, y, theta),
        t_span=(t_span.min(), t_span.max()),
        y0=y0,
        t_eval=t_span,
    )

    return sol


def distance(sol_exp_y, sol_y):
    #print(sol_exp_y['data'].shape,sol_y['data'].shape)
    try: return MSE(sol_exp_y['data'].T[:, [0, 2, 4]], sol_y['data'].T[:, [0, 2, 4]])
    except: return 100000
 

def distance2(sol_exp_y, sol_y):
    #print(sol_exp_y['data'].shape,sol_y['data'].shape)
    try: return MSE(sol_exp_y['data'].T, sol_y['data'].T)
    except: return 100000
    


def plot_sol(t,sol,outname=None):
    fig, axs = plt.subplots(3, 2)
    ylables = ['$m_1$', '$p_1$', '$m_2$', '$p_2$', '$m_3$', '$p_3$']
    for i in range(0, 6):
        ax = axs[i // 2, i % 2]
        ax.plot(t, sol.T[:, i])
        ax.set_xlabel("Time ")
        ax.set_ylabel(ylables[i])
    plt.tight_layout()
    if outname is not None: plt.savefig(outname)
    plt.show()


def dydt_func_p53(t, y, theta):
    '''

    '''
    A, B, C, D = y

    k1, k2, k3, k4, k5, k6, k7, k8, k9 = theta

    dydt = np.zeros(4)
    dydt[0] = k1 - k2 * A - k3 * A * B + k4 * C + k5 * C
    dydt[1] = k6 * D - k3 * A * B + k4 * C - k5 * B
    dydt[2] = k3 * A * B - k4 * C - k7 * C - k5 * C
    dydt[3] = k8 * (A ** 2) - k9 * D

    return dydt


def do_solve_ivp_p53(t_span, y0, theta):
    sol = solve_ivp(
        fun=lambda t, y: dydt_func_p53(t, y, theta),
        t_span=(t_span.min(), t_span.max()),
        y0=y0,
        #dense_output=True,
        t_eval=t_span,
        method='BDF',
    )

    return sol


def plot_sol_p53(t, sol, outname=None):
    fig, axs = plt.subplots(2, 2)
    ylables = ['A (nuclear-p53)', 'B (Mdm2)', 'C (p53-Mdm2)', 'D (Mdm2 mRNA)']
    for i in range(0, 4):
        ax = axs[i // 2, i % 2]
        ax.plot(t, sol.T[:, i])
        ax.set_xlabel("Time ")
        ax.set_ylabel(ylables[i])
    plt.tight_layout()
    if outname is not None: plt.savefig(outname)
    plt.show()


def info():
    print('\n ---------- ode toy model information ----------')
    print('''
        ODE Parameters Description: theta
        theta[0]: alpha0: a0
        theta[1]: n: n
        theta[2]: beta: b,
        theta[3]: alpha: a,

        Mathematical Description:
        dm1/dt = -m1 + a/(1+p3**n) + a0
        dp1/dt = -b*(p1-m1)
        dm2/dt = -m2 + a/(1+p1**n) + a0
        dp2/dt = -b*(p2-m2)
        dm3/dt = -m3 + a/(1+p2**n) + a0
        dp2/dt = -b*(p2-m2)''')

    print(''' default:
    t_span = np.linspace(0, 10, 101)
    y0 = [0, 2, 0, 1, 0, 3]     # ref
    default_parameters = [1, 2, 5, 1000]  # ref  
    ''')


def get_experiment_data(t_span, y0=[0, 2, 0, 1, 0, 3], default_parameters=[1, 2, 5, 1000]):
    sol_exp = do_solve_ivp(t_span, y0, theta=default_parameters)

    # print('add noise as reference (N(0,5**2)) Note: please check')
    noise = np.random.normal(0, 5, sol_exp.y.shape)
    sol_exp.y = sol_exp.y + noise
    return sol_exp


if __name__ == '__main__':
    t_span_p53 = np.linspace(0, 3, 101)
    y0_p53 = [0, 0, 0, 0, ]  # ref
    default_parameters_p53 = np.array([10, 2, 30, 40, 5, 60, 7, 80, 9])*0.1  # ref

    sol_exp_p53 = do_solve_ivp_p53(t_span_p53, y0_p53, theta=default_parameters_p53)
    plot_sol_p53(t_span_p53, sol_exp_p53.y, outname=None)

    # info()
    # t_span = np.linspace(0, 10, 101)
    # y0 = [0, 2, 0, 1, 0, 3]  # ref
    # default_parameters = [1, 2, 5, 1000]  # ref
    #
    # print('under default_parameters,the sol could be sol_exp: ')
    # sol_exp = do_solve_ivp(t_span, y0, theta=default_parameters)
    # plot_sol(sol_exp)
    #
    # print('add noise as reference (N(0,5**2)) Note: please check')
    # noise = np.random.normal(0, 5, sol_exp.y.shape)
    # sol_exp.y = sol_exp.y + noise
    # plot_sol(sol_exp)
    #
    # print('''
    # set new_parameters by: new_parameters = [1, 2, 5, 1000]
    # new_sol = do_solve_ivp(t_span, y0, new_parameters)
    # ''')
    # new_parameters = [1, 2, 5, 1000]
    # new_sol = do_solve_ivp(t_span, y0, new_parameters)
    # plot_sol(new_sol)
    # print('distance: only consider m1,m2,m3')
    #
    # print('distance(sol_exp.y, new_sol.y): ', distance(sol_exp.y, new_sol.y))
    #
    # print('r2_score(sol_exp.y.T[:,[0,2,4]], new_sol.y.T[:,[0,2,4]]): ',
    #       r2_score(sol_exp.y.T[:, [0, 2, 4]], new_sol.y.T[:, [0, 2, 4]]))
