/// Returns Vector of liquid rates
pub fn get_liquid_rates(u: &Vec<u8>, r: &Vec<f64>) -> Vec<f64> {
    let mut liquid_rates: Vec<f64> = vec![0.0; u.len()];

    for i in 0..u.len() {
        let ui = (u[i] - 1) as usize;
        liquid_rates[i] = r[ui];
    }
    liquid_rates
}