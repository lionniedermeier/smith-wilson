mod utils;
mod smith_wilson;

#[cfg(test)]
mod tests {
    use crate::smith_wilson::{Instrument, SmithWilson};
    use crate::utils::get_liquid_rates;

    fn get_u() -> Vec<u8> {
        let u: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 19, 20];
        u
    }
    fn get_r() -> Vec<f64> {
        let r: Vec<f64> = vec![
            0.9250, 1.0590, 1.1032, 1.1044, 1.3788, 1.6944, 1.8429, 1.8642, 2.0102, 2.1872, 2.2981,
            2.4615, 2.5868, 2.6940, 2.7786, 2.7962, 2.8432, 2.8665, 2.8957, 2.9492,
        ].iter().map(|x| x / 100.0).collect();
        r
    }

    #[test]
    fn test_smith_wilson() {
        let ufr = 3.3 / 100.0;
        let u = get_u();
        let r = get_r();
        let mut r_liq = get_liquid_rates(&u, &r);
        let instrument = Instrument::SWAP;
        let cra = 10;

        for i in 0..r_liq.len() {
            r_liq[i] =  r_liq[i] - f64::from(cra / 10000);
        }

        let smith_wilson = SmithWilson::new(
            u,
            r_liq,
            instrument,
            ufr,
            0.05,
            1,
            1.0,
            40,
            20
        );
        smith_wilson.fit();
    }
}
