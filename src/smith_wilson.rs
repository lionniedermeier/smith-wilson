//! This module implements the Smith-Wilson extrapolation method with intensities,
//! which is used for yield curve construction and extrapolation
//! beyond the last liquid point as specified in the EIOPA guidelines.

use nalgebra::{DMatrix, DVector};

/// Implements the Smith-Wilson method with intensities
///
/// # Parameters
/// - `maturities`: liquid maturities
/// - `rates`: liquid rates
/// - `ufr`: ultimate forward rate.
/// - `instrument`: the type of instrument zero, swap.
/// - `alpha_min`: lower bound for alpha.
/// - `num_of_coupon`: coupon payments per year.
/// - `tau`: Convergence tolerance acceptable deviation from the convergence point.
/// - `convergence_period`: year to convergence point from last-liquid-point.
/// - `llp`: Last-liquid-point.
pub struct SmithWilson {
    pub maturities: Vec<u8>,
    pub rates: Vec<f64>,
    pub ufr: f64,
    pub instrument: Instrument,
    pub alpha_min: f64,
    pub num_of_coupon: u8,
    pub tau: f64,
    pub convergence_period: u8,
    pub llp: u8,
}

impl SmithWilson {
    pub fn new(
        maturities: Vec<u8>,
        rates: Vec<f64>,
        instrument: Instrument,
        ufr: f64,
        alpha_min: f64,
        num_of_coupon: u8,
        tau: f64,
        convergence_period: u8,
        llp: u8,
    ) -> Self {
        Self {
            maturities,
            rates,
            instrument,
            ufr,
            alpha_min,
            num_of_coupon,
            tau,
            convergence_period,
            llp,
        }
    }

    pub fn fit(&self) {
        let ln_ufr = (1.0 + self.ufr).ln();

        let q_mat = q_matrix(
            &self.maturities,
            &self.rates,
            ln_ufr,
            &self.instrument,
            self.num_of_coupon,
        );

        let precision = 6;
        let t2 = (self.llp + self.convergence_period).max(60) as f64;
        let num_of_rates = self.rates.len() as u16;
        let umax = *self.maturities.last().unwrap();
        let tau = self.tau / 10000.0;

        // let (error, _gamma) = g_alpha(
        //     self.alpha_min,
        //     &q_mat,
        //     num_of_rates,
        //     umax,
        //     self.num_of_coupon,
        //     t2,
        //     tau,
        // );

        let (alpha, gamma) = optimize_alpha(
            self.alpha_min,
            &q_mat,
            num_of_rates,
            umax,
            self.num_of_coupon,
            t2,
            tau,
            precision,
        );

        let tenors: u8 = 150;

        let mut h: DMatrix<f64> = DMatrix::zeros(150, (umax * self.num_of_coupon).into());
        let mut g: DMatrix<f64> = DMatrix::zeros(150, (umax * self.num_of_coupon).into());

        for i in 0..tenors {
            for j in 0..u16::from(umax * self.num_of_coupon) {
                let i_f: f64 = f64::from(i);
                let j_f: f64 = f64::from(j);
                let num_of_coupon_f = f64::from(self.num_of_coupon);

                h[(i as usize, j as usize)] =
                    h_mat(alpha * i_f, alpha * (j_f + 1.0) / num_of_coupon_f);
                if ((j + 1) / u16::from(self.num_of_coupon)) > i.into() {
                    g[(i as usize, j as usize)] = alpha
                        * (1.0 - (-alpha * (j_f + 1.0) / num_of_coupon_f).exp())
                        * (alpha * i_f).cosh();
                } else {
                    g[(i as usize, j as usize)] = alpha
                        * (-alpha * i_f).exp()
                        * (alpha * (j_f + 1.0) / num_of_coupon_f).sinh();
                }
            }
        }

        let tt_discount: DVector<f64>;
        let tt_intensity: DVector<f64>;
        tt_discount = h * gamma.clone();
        tt_intensity = g * gamma.clone();
        // Debug
        println!(""); // DEBUG
        for i in 0..19 {
            println!("Gamma {}", gamma[i as usize]);
        }
        println!(""); // DEBUG

        let mut temp: f64 = 0.0;
        for i in 0..u16::from(umax * self.num_of_coupon) {
            temp = temp
                + (1.0 - ((-alpha * f64::from(i + 1) / f64::from(self.num_of_coupon)).exp()))
                    * gamma[(i as usize, 0)];
        }

        let mut discount: DVector<f64> = DVector::zeros(150);
        let mut fw_intensity: DVector<f64> = DVector::zeros(150);
        let mut yld_intensity: DVector<f64> = DVector::zeros(150);
        let mut forward_ac: DVector<f64> = DVector::zeros(150);
        let mut zero_ac: DVector<f64> = DVector::zeros(150);
        // let mut forward_cc: DVector<f64> = DVector::zeros(150);
        // let mut zero_cc: DVector<f64> = DVector::zeros(150);

        yld_intensity[0] = ln_ufr - alpha * temp;
        fw_intensity[0] = yld_intensity[0];
        discount[0] = 1.0;

        for i in 1..tenors {
            yld_intensity[i as usize] = ln_ufr - (1.0 + tt_discount[i as usize] / f64::from(i));
            fw_intensity[i as usize] =
                ln_ufr - tt_intensity[i as usize] / (1.0 + tt_discount[i as usize]);
            discount[i as usize] = (-ln_ufr * f64::from(i)).exp() * (1.0 + tt_discount[i as usize]);
            zero_ac[i as usize] = (1.0 / discount[i as usize]).powf(1.0 / f64::from(i)) - 1.0;
            forward_ac[i as usize] = discount[(i - 1) as usize] / discount[i as usize] - 1.0;
        }

        // DEBUG
        // for i in 0..tenors {
        //     println!("{}", zero_ac[i as usize]);
        // }
    }
}

pub enum Instrument {
    ZERO,
    SWAP,
}

struct SmithWilsonResult {
    pub zero_coupon_rate: Vec<f64>,
    pub yield_intensity: Vec<f64>,
    pub forward_intensity: Vec<f64>,

}



/// Calculates the Q matrix
///
/// # Parameters
/// - `u`: Vector of liquid maturities
/// - `r`: Vector of liquid rates
/// - `ufr`: ultimate forward rate
/// - `instrument`: type of instrument either zero or swap
/// - `num_of_coupon`: coupon payments per year
pub fn q_matrix(
    u: &Vec<u8>,
    r: &Vec<f64>,
    ufr: f64,
    instrument: &Instrument,
    num_of_coupon: u8,
) -> DMatrix<f64> {
    let num_of_rates: u8 = r.len().try_into().unwrap();
    let mut q_mat: DMatrix<f64>;
    let umax = u[u.len() - 1];

    match instrument {
        Instrument::ZERO => {
            q_mat = DMatrix::zeros(num_of_rates.into(), (umax * num_of_coupon).into());

            for i in 0..num_of_rates {
                let i = i as usize;
                let ui = (u[i] - 1) as usize;
                let uf = f64::from(u[i]);
                q_mat[(i, ui)] = (-ufr * uf).exp() * ((1.0 + r[i]).powf(uf));
            }
        }
        _ => {
            q_mat = DMatrix::zeros(num_of_rates.into(), (umax * num_of_coupon).into());

            for i in 0..num_of_rates {
                let i = i as usize;
                let ui = u[i];
                let num_of_coupon_f = f64::from(num_of_coupon);
                let mut j: u16 = 0;

                for _ in 0..(num_of_coupon * ui - 1) {
                    let jf = f64::from(j + 1);
                    let ju = j as usize;
                    q_mat[(i, ju)] = (-ufr * jf / num_of_coupon_f).exp() * (r[i] / num_of_coupon_f);
                    j += 1;
                }
                q_mat[(i, usize::from(j))] = (-ufr * f64::from(j + 1) / num_of_coupon_f).exp()
                    * (1.0 + r[i] / num_of_coupon_f);
            }
        }
    }
    q_mat
}

/// Helper function
fn _h(z: f64) -> f64 {
    (z + (-z).exp()) / 2.0
}

/// Calculates an element from the H matrix
pub fn h_mat(u: f64, v: f64) -> f64 {
    _h(u + v) - _h((u - v).abs())
}

/// Computes the absolute difference to the foward-intensity function at the convergence point
///
/// # Parameters
///
/// - `alpha`: A scalar parameter influencing the convergence behavior.
/// - `q_mat`: A matrix (DMatrix) containing transition or rate values.
/// - `num_of_rates`: The number of rate categories or transition states.
/// - `umax`: The upper limit for iterations or model parameters (e.g., maximum usage or level).
/// - `num_of_coupon`: coupon payments per year
/// - `t2`: Convergence point — the target or expected value the model should approach.
/// - `tau`: Convergence tolerance — acceptable deviation from the convergence point.
///
/// # Returns
///
/// A tuple containing:
/// - the difference g(a) - tau
/// - the vector Qb
fn g_alpha(
    alpha: f64,
    q_mat: &DMatrix<f64>,
    num_of_rates: u16,
    umax: u8,
    num_of_coupon: u8,
    t2: f64,
    tau: f64,
) -> (f64, DVector<f64>) {
    let mut h: DMatrix<f64> =
        DMatrix::zeros((umax * num_of_coupon).into(), (umax * num_of_coupon).into());

    let num_of_coupon_f = f64::from(num_of_coupon);

    for i in 0..umax * num_of_coupon {
        for j in 0..umax * num_of_coupon {
            h[(i as usize, j as usize)] = h_mat(
                alpha * f64::from(i + 1) / num_of_coupon_f,
                alpha * f64::from(j + 1) / num_of_coupon_f,
            )
        }
    }
    
    let temp1 = DVector::from_element(num_of_rates.into(), 1.0) - q_mat.column_sum();
    let temp_mat = ((q_mat * h) * q_mat.transpose()).try_inverse();
    if temp_mat.is_none() {
        // return Err("Matrix is not invertible");
        panic!("Matrix is not invertible");

    }
    let b = temp_mat.unwrap() * temp1;
    let q_b: DVector<f64> = q_mat.transpose() * b;

    let mut temp2: f64 = 0.0;
    let mut temp3: f64 = 0.0;
    for i in 0..umax * num_of_coupon {
        temp2 = temp2 + q_b[i as usize] * f64::from(i + 1) / num_of_coupon_f;
        temp3 = temp3 + q_b[i as usize] * (alpha * f64::from(i + 1) / num_of_coupon_f).sinh();
    }

    let kappa = (1.0 + alpha * temp2) / temp3;
    let g_a = alpha / (1.0 - kappa * (t2 * alpha).exp()).abs() - tau;
    (g_a, q_b)
}

fn scan_for_alpha(
    prev_alpha: f64,
    stepsize: f64,
    q_mat: &DMatrix<f64>,
    num_of_rates: u16,
    umax: u8,
    num_of_coupon: u8,
    t2: f64,
    tau: f64,
    digit: u8
) -> (f64, DVector<f64>) {
    let mut alpha = prev_alpha - stepsize + stepsize / 10.0;
    let mut g_alpha_out: (f64, DVector<f64>) = (0.0, DVector::zeros(1));
    while alpha <= prev_alpha + stepsize / 10.0 {
        // alpha = round_to_decimal(alpha, stepsize * 1000.0);
        // let factor = 100.0 * 10f64.powi(digit as i32); Rounding 
        // println!("Alpha D {:.12} Times {:.4}", trunc_to_decimal(alpha, i32::from(digit + 2)), factor);
        g_alpha_out = g_alpha(alpha, &q_mat, num_of_rates, umax, num_of_coupon, t2, tau);
        println!("Alpha Scan: {}, error: {}", alpha, g_alpha_out.0);
        if g_alpha_out.0 <= 0.0 {
            break;
        }
        alpha += stepsize / 10.0;
    }
    (alpha, g_alpha_out.1)
}

fn optimize_alpha(
    alpha_min: f64,
    q_mat: &DMatrix<f64>,
    num_of_rates: u16,
    umax: u8,
    num_of_coupon: u8,
    t2: f64,
    tau: f64,
    precision: u8,
) -> (f64, DVector<f64>) {
    let g_alpha_out = g_alpha(
        alpha_min,
        &q_mat,
        num_of_rates,
        umax,
        num_of_coupon,
        t2,
        tau,
    );
    let mut alpha: f64;
    let mut gamma = DVector::zeros((umax * num_of_coupon) as usize);

    if g_alpha_out.0 <= 0.0 {
        alpha = alpha_min;
        gamma = g_alpha_out.1;
    } else {
        let mut stepsize = 0.1;
        alpha = alpha_min;

        while alpha < 20.0 {
            println!("Alpha 1: {}", alpha); // DEBUG
            if g_alpha(alpha, &q_mat, num_of_rates, umax, num_of_coupon, t2, tau).0 <= 0.0 {
                break;
            }
            alpha += stepsize;
        }

        for digit in 0..precision {
            let (alpha_out, gamma_out) = scan_for_alpha(
                alpha,
                stepsize,
                &q_mat,
                num_of_rates,
                umax,
                num_of_coupon,
                t2,
                tau,
                digit
            );
            alpha = alpha_out;
            gamma = gamma_out;
            stepsize /= 10.0;
        }
    }
    println!("Alpha Ende: {}", alpha); // DEBUG

    (alpha, gamma)
}

fn trunc_to_decimal(mut v: f64, digits: i32) -> f64 {
    v = (v * 10_f64.powi(digits)).trunc();
    v = v / 10_f64.powi(digits);
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_u() -> Vec<u8> {
        let u: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 19, 20];
        u
    }
    fn get_r() -> Vec<f64> {
        let r: Vec<f64> = vec![
            0.9250, 1.0590, 1.1032, 1.1044, 1.3788, 1.6944, 1.8429, 1.8642, 2.0102, 2.1872, 2.2981,
            2.4615, 2.5868, 2.6940, 2.7786, 2.7962, 2.8432, 2.8665, 2.8957, 2.9492,
        ]
        .iter()
        .map(|x| x / 100.0)
        .collect();
        r
    }

    // #[test]
    // fn test_get_liquid_rates() {
    //     let u = get_u();
    //     let r = get_r();
    //     let liq_r: Vec<f64> = get_liquid_rates(&u, &r);
    //     println!("{:?}", liq_r);
    // }

    // #[test]
    // fn test_qmat() {
    //     let u = get_u();
    //     let r = get_r();
    //     let rl = get_liquid_rates(&u, &r);
    //     let ufr = (1.0 + f64::from(3.3 / 100.0)).ln();
    //     let instrument = Instrument::SWAP;
    //     let q_mat = q_matrix(&u, &rl, ufr, &instrument, 1);
    //     // Print matrix
    //     for row in q_mat.row_iter() {
    //         let row_out = row
    //             .iter()
    //             .map(|v| format!("{:.5}", v))
    //             .collect::<Vec<_>>()
    //             .join(" ");
    //         println!("{}", row_out);
    //     }
    // }

    // #[test]
    // fn test_g_alpha() {
    //     let u = get_u();
    //     let r = get_r();
    //     let rl = get_liquid_rates(&u, &r);
    //     let ufr = (1.0 + f64::from(3.3)).ln();
    //     let instrument = Instrument::SWAP;
    //     let q_mat = q_matrix(&u, &rl, ufr, &instrument, 1);
    //     let umax = u[(u.len() - 1) as usize];
    //     println!("{}", rl.len());
    //     let num_of_rates = rl.len().try_into().unwrap();
    //     let aout = g_alpha(0.05, &q_mat, num_of_rates, umax, 1, 60.0, 1.0 / 10000.0);
    //     println!("{:?}", aout.1)
    // }
}
