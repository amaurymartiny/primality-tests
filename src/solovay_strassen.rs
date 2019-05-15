extern crate rug;

use rug::rand::RandState;
use rug::Integer;

pub fn solovay_strassen(p: u128, iterations: u32) -> bool {
	if p <= 2 {
		return false;
	}

	let mut rand = RandState::new();
	for _ in 0..iterations {
		// Take a random integer in [2, p-1]
		let a: Integer = (Integer::from(p - 2)).random_below(&mut rand) + 2;
		// Get the jacobian (a/p) mod p
		let (_, jacobian) =
			Integer::from(a.jacobi(&Integer::from(p))).div_rem_floor(Integer::from(p));
		// Calculate a^((p-1)/2) mod p
		let modulo = a
			.pow_mod(
				&((Integer::from(p - 1)) / Integer::from(2)),
				&Integer::from(p),
			)
			.unwrap();

		// If p is prime, then jacobian == modulo
		if jacobian.eq(&0) || !jacobian.eq(&modulo) {
			return false;
		}
	}

	return true;
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_composite_2143() {
		assert_eq!(solovay_strassen(1024, 20), false);
	}

	#[test]
	fn test_prime_2143() {
		assert_eq!(solovay_strassen(2143, 20), true);
	}
}
