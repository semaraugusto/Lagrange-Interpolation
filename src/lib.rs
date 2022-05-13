#![allow(dead_code)]
// TODO: Change prime field from bls12 to F_5
// use ark_bls12_381::Fr;
use ark_ff::PrimeField;
use itertools::Itertools;

struct Lagrange<P: PrimeField, F: Fn(&Vec<bool>) -> P> {
    func: F,
    domain: usize,
}

impl<P: PrimeField, F: Fn(&Vec<bool>) -> P> Lagrange<P, F> {
    fn new(func: F, domain: usize) -> Self {
        Self { func, domain }
    }

    fn slow_interpolate(&self) -> impl Fn(Vec<P>) -> P {
        let (w_vals, f_evals) = self.get_inputs_and_evaluations();
        let one = P::from(1u8);

        move |input: Vec<P>| {
            print!(" ");
            let chi = w_vals.iter().map(|w| {
                w.iter()
                    .zip(input.iter())
                    .map(|(w_i, x_i)| if *w_i { *x_i } else { one - *x_i })
                    .product::<P>()
            });

            f_evals.iter().zip(chi).map(|(a, b)| *a * b).sum::<P>()
        }
    }
    // pretty sure this is not optimized yet and is just a rewritten version of the above function
    fn interpolate(&self) -> impl Fn(Vec<P>) -> P {
        let (w_vals, f_evals) = self.get_inputs_and_evaluations();
        let one = P::from(1u8);
        let mut num_mult = 0;
        move |input: Vec<P>| {
            let chi = w_vals.iter().fold(Vec::new(), |mut chi, w| {
                let mut chi_r = one;
                for (w_i, x_i) in w.iter().zip(input.iter()) {
                    let val = if *w_i { *x_i } else { one - *x_i };
                    chi_r *= val;
                }
                chi.push(chi_r);
                chi
            });

            f_evals.iter().zip(chi).map(|(a, b)| *a * b).sum::<P>()
        }
    }

    fn get_inputs_and_evaluations(&self) -> (Vec<Vec<bool>>, Vec<P>) {
        let mut w_vals = Vec::new();
        let mut f_evals = Vec::new();
        for w in (0..self.domain)
            .map(|_| [false, true])
            .multi_cartesian_product()
        {
            f_evals.push((self.func)(&w));
            w_vals.push(w);
        }
        (w_vals, f_evals)
    }
}

// Function f mapping {0,1}^2 to prime field F_5 (see Figure 3.1)
fn f<P: PrimeField>(x: &Vec<bool>) -> P {
    assert_eq!(x.len(), 2, "invalid input domain size");
    P::from(match (x[0] as u8, x[1] as u8) {
        (0, 0) => 1,
        (0, 1) => 2,
        (1, 0) => 1,
        (1, 1) => 4,
        _ => panic!(),
    } as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{fields::Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "5"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    fn test_func(a: u32, b: u32) -> Fq {
        Fq::from(match (a, b) {
            (0, 0) => 1,
            (0, 1) => 2,
            (0, 2) => 3,
            (0, 3) => 4,
            (0, 4) => 5,

            (1, 0) => 1,
            (1, 1) => 4,
            (1, 2) => 2,
            (1, 3) => 0,
            (1, 4) => 3,

            (2, 0) => 1,
            (2, 1) => 1,
            (2, 2) => 1,
            (2, 3) => 1,
            (2, 4) => 1,

            (3, 0) => 1,
            (3, 1) => 3,
            (3, 2) => 0,
            (3, 3) => 2,
            (3, 4) => 4,

            (4, 0) => 1,
            (4, 1) => 0,
            (4, 2) => 4,
            (4, 3) => 3,
            (4, 4) => 2,
            _ => panic!(),
        } as u32)
    }
    #[test]
    fn test_slow() {
        for x in (0..2).map(|_| (0..5)).multi_cartesian_product() {
            let a = test_func(x[0], x[1]);
            let lagrange = Lagrange::new(f, 2);
            let f_tilde = lagrange.slow_interpolate();
            let b = f_tilde(vec![Fq::from(x[0]), Fq::from(x[1])]);
            println!("{x:?}{a:?} {b:?}");
            assert_eq!(
                a, b,
                "incorrect evaluation of multilinear extension f_tilde"
            );
        }
    }
    #[test]
    fn test_fast() {
        for x in (0..2).map(|_| (0..5)).multi_cartesian_product() {
            let a = test_func(x[0], x[1]);
            let lagrange = Lagrange::new(f, 2);
            let f_tilde = lagrange.interpolate();
            let b = f_tilde(vec![Fq::from(x[0]), Fq::from(x[1])]);
            println!("{x:?}{a:?} {b:?}");
            assert_eq!(
                a, b,
                "incorrect evaluation of multilinear extension f_tilde"
            );
        }
    }
}
