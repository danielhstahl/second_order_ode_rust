#[macro_use]
#[cfg(test)]
extern crate approx;


fn compute_dx(x_discrete:usize, x_min:f64, x_max:f64)->f64{
    (x_max-x_min)/((x_discrete as f64)-1.0) 
}

fn compute_x(index:usize, dx:f64, x_min:f64)-> f64 {
    x_min+(index as f64)*dx
}

fn thomas_algorithm<'a, Id>(
    diag:Id
)->Vec<f64> 
where 
    Id:Iterator<Item=(Option<f64>, f64, Option<f64>, f64)>, //lower, main, upper, solution
{
    let mut upper_v:Vec<f64>=vec![];
    let mut solve_v:Vec<f64>=vec![];
    for (index, (lower, main, upper, sol)) in diag.enumerate(){
        if lower.is_some()  {
            let upper_v_prev=upper_v[index-1];
            let solve_v_prev=solve_v[index-1];
            if upper.is_some() {
                upper_v.push(upper.unwrap()/(main-upper_v_prev*lower.unwrap()));
            }
            solve_v.push((sol-lower.unwrap()*solve_v_prev)/(main-upper_v_prev*lower.unwrap()))
        }
        else if lower.is_none(){
            upper_v.push(upper.unwrap()/main);
            solve_v.push(sol/main);
        }
    }

    for (index, cprime) in upper_v.iter().enumerate().rev(){
        solve_v[index]=solve_v[index]-cprime*solve_v[index+1];
    }
    solve_v  
}

/// Solves ODEs of the form fn2(x)*f''(x)+fn1(x)*f'(x)+fn*f(x)=0
/// # Examples
/// ```
/// let fn2=|_|1.5;
/// let fn1=|_|5.0;
/// let fnc=|_|1.5;
/// let init_cond_lower=0.0;
/// let init_cond_upper=1.0;
/// let x_min=0.0;
/// let x_max=1.0;
/// let n=100;
/// let result=second_order_ode::solve_ode(
///     &fn2, &fn1, 
///     &fnc, init_cond_lower, init_cond_upper, 
///     x_min, x_max, n
/// );
/// ```
pub fn solve_ode(
    second_deriv_coef:&Fn(f64)->f64,
    first_deriv_coef:&Fn(f64)->f64,
    fn_coef:&Fn(f64)->f64,
    initial_condition_lower:f64,
    initial_condition_upper:f64,
    x_min:f64,
    x_max:f64,
    num_steps:usize
)->Vec<f64>{
    let dx=compute_dx(num_steps+2, x_min, x_max);
    let dx_sq=dx.powi(2);
    let dx2=dx*2.0;
    let get_upper_coef=|index:usize|{
        let x=compute_x(index, dx, x_min);
        second_deriv_coef(x)/dx_sq+first_deriv_coef(x)/dx2
    };
    let get_lower_coef=|index:usize|{
        let x=compute_x(index, dx, x_min);
        second_deriv_coef(x)/dx_sq-first_deriv_coef(x)/dx2
    };
    let get_main_coef=|index:usize|{
        let x=compute_x(index, dx, x_min);
        fn_coef(x)-second_deriv_coef(x)*2.0/dx_sq
    };
    let diag=(1..num_steps+1).map(|index|{
        if index==1 {
            (
                None, 
                get_main_coef(index), 
                Some(get_upper_coef(index+1)), 
                -initial_condition_lower*get_lower_coef(index)
            )
        }
        else if index==num_steps {
            (
                Some(get_lower_coef(index-1)), 
                get_main_coef(index), 
                None, 
                -initial_condition_upper*get_upper_coef(index)
            )
        }
        else {
            (
                Some(get_lower_coef(index-1)), 
                get_main_coef(index), 
                Some(get_upper_coef(index+1)),
                0.0
            )
        }
    });
    thomas_algorithm(diag)

}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_thomas_algorithm() {
        let mut diag:Vec<(Option<f64>, f64, Option<f64>, f64)>=vec![];
        diag.push((None, 0.3, Some(0.9), 0.3));
        diag.push((Some(0.4), 0.4, Some(0.2), -0.5));
        diag.push((Some(0.6), 0.2, None, 0.3));
        let expected:Vec<f64>=vec![
            -1.5714286,
            0.8571429,
            -1.0714286
        ];
        let result=thomas_algorithm(diag.iter().map(|v|*v));
        for (res, ex) in result.iter().zip(expected.iter()){
            assert_abs_diff_eq!(res, ex, epsilon=0.000001);
        }
    }

    #[test]
    fn test_solve_ode(){
        let fn2=|_|1.5;
        let fn1=|_|5.0;
        let fnc=|_|1.5;
        let init_cond_lower=0.0;
        let init_cond_upper=1.0;
        let x_min=0.0;
        let x_max=1.0;
        let n=100;
        let expected_fnc=|x:f64|{
            let coef1:f64=-1.0/3.0;
            let coef2:f64=-3.0;

            let c2=1.0/(coef1.exp()-coef2.exp());
            let c1=-c2;
            c1*(coef2*x).exp()+c2*(coef1*x).exp()
        };
        let result=solve_ode(
            &fn2, &fn1, 
            &fnc, init_cond_lower, init_cond_upper, 
            x_min, x_max, n
        );
        let dx=compute_dx(n+2, x_min, x_max);
        for (index, res) in result.iter().enumerate(){
            assert_abs_diff_eq!(*res, expected_fnc(compute_x(index+1, dx, x_min)), epsilon=0.001);
        }
    }
}
