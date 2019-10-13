| [Linux][lin-link] |  [Codecov][cov-link]  |
| :---------------: | :-------------------: |
| ![lin-badge]      | ![cov-badge]          |

[lin-badge]: https://github.com/phillyfan1138/second_order_ode_rust/workflows/Rust/badge.svg
[lin-link]:  https://github.com/phillyfan1138/second_order_ode_rust/actions
[cov-badge]: https://codecov.io/gh/phillyfan1138/second_order_ode_rust/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/phillyfan1138/second_order_ode


## ODE Solver for second order differential equations

This repository can solve second order ODEs: h(x)f''(x)+g(x)f'(x)+c(x)f(x)=0.  The user must specify 
* the functions h(x), g(x), and c(x)
* two boundary conditions (at f(xmin) and f(xmax)).
* the domain of x (xmin and xmax)
* the number of discrete points to evaluate the function at.  