# SZO-QQ
Safe-Zeroth-Order-using-sequential-Quadratic-Constrained-Quadratic-Programming

Code for the paper "Safe Zeroth-Order Optimization Using Quadratic Local Approximations", co-authored by Baiwei Guo, Yuning Jiang, Giancarlo Ferrari-Trecate and Maryam Kamgarpour

####################

In the folder **6_1_Toy_Example**, the plots in the paper can be generated by plot/plot_1.m, plot/plot_2.m and plot/plot_3_4.m. The source for SZO-QQ is qp_optimizer.m. Sources for ES and LB-SGD can be found in the folders **ES_TOOL** and **log_barrier_code**.
To run the code in this folder, the user needs Mosek and Yalmip.

Mosek: https://www.mosek.com
Yalmip: https://yalmip.github.io/

####################
In the folder **6_2_Optimal_Control**, the source is Example_6_2.m
To run ths code in this folder, the user needs Casadi and IPOPT

Casadi: https://web.casadi.org
IPOPT: https://coin-or.github.io/Ipopt/

###################

In the folder **6_3_OPF**, the plots in the paper can be generated by plot_5.m. The sources for three methods can be found in the folder **code**
To run the code in this folder, the user needs Matpower (for generating the reference).