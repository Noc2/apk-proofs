use ark_bw6_761::{Fr as F, BW6_761};
use ark_ec::ProjectiveCurve;
use ark_ec::short_weierstrass_jacobian::GroupAffine;
use ark_ff::{FftField, Field, One, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial, UVPolynomial, Radix2EvaluationDomain};
use ark_poly::univariate::DensePolynomial;
use merlin::Transcript;

use crate::{KZG_BW6, Proof, point_in_g1_complement, Bitmask};
use crate::transcript::ApkTranscript;
use crate::signer_set::SignerSetCommitment;
use crate::kzg::ProverKey;
use crate::bls::PublicKey;
use crate::utils;

fn mul<F: Field>(s: F, p: &DensePolynomial<F>) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs.iter().map(|c| s * c).collect()
    )
}

fn mul_by_x<F: Field>(p: &DensePolynomial<F>) -> DensePolynomial<F> {
    let mut px = vec![F::zero()];
    px.extend_from_slice(&p.coeffs);
    DensePolynomial::from_coefficients_vec(px)
}

fn add_constant<F: FftField, D: EvaluationDomain<F>>(p: &Evaluations<F, D>, c: F, d: D) ->  Evaluations<F, D> {
    Evaluations::from_vec_and_domain(p.evals.iter().map(|x| c + x).collect(), d)
}

pub struct Prover<'a> {
    domain_size: usize,
    kzg_pk: ProverKey<BW6_761>,
    pks: &'a[PublicKey],
    h: ark_bls12_377::G1Affine,
    preprocessed_transcript: Transcript,
}

impl<'a> Prover<'a> {
    pub fn new(
        domain_size: usize,
        kzg_pk: ProverKey<BW6_761>,
        signer_set_comm: &SignerSetCommitment,
        pks: &'a[PublicKey],
        mut empty_transcript: Transcript,
    ) -> Self {
        assert!(domain_size.is_power_of_two(), "domain size should be a power of 2");
        assert!(domain_size <= kzg_pk.max_coeffs(), "domain size shouldn't exceed srs length");
        // empty_transcript.set_protocol_params(); //TODO
        empty_transcript.set_signer_set(&signer_set_comm);
        Self { domain_size, kzg_pk, pks, h: point_in_g1_complement(), preprocessed_transcript: empty_transcript }
    }

    #[allow(non_snake_case)]
    pub fn prove(&self, b: &Bitmask) -> Proof {
        let m = self.pks.len();

        assert_eq!(b.size(), m);
        assert!(b.count_ones() > 0);

        let apk = b.to_bits().iter()
            .zip(self.pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p.0)
            .sum::<ark_bls12_377::G1Projective>();

        let mut transcript = self.preprocessed_transcript.clone();
        transcript.append_public_input(&apk.into(), b);

        let (pks_x, pks_y): (Vec<F>, Vec<F>) = self.pks.iter()
            .map(|p| p.0.into_affine())
            .map(|p| (p.x, p.y))
            .unzip();

        let mut acc = vec![self.h; m+1];
        for (i, (b, p)) in b.to_bits().iter().zip(self.pks.iter()).enumerate() {
            acc[i+1] = if *b {
                acc[i] + p.0.into_affine()
            } else {
                acc[i]
            }
        }

        let (mut acc_x, mut acc_y): (Vec<F>, Vec<F>) = acc.iter()
            .map(|p| (p.x, p.y))
            .unzip();


        assert_eq!(b.size(), m);
        assert_eq!(pks_x.len(), m);
        assert_eq!(pks_y.len(), m);
        assert_eq!(acc_x.len(), m+1);
        assert_eq!(acc_y.len(), m+1);
        assert_eq!(GroupAffine::new(acc_x[0], acc_y[0], false), self.h);
        assert_eq!(GroupAffine::new(acc_x[m], acc_y[m], false), apk.into_affine() + self.h);

        let mut b = b.to_bits().iter()
            .map(|b| if *b { F::one() } else { F::zero() })
            .collect::<Vec<_>>();

        let n = self.domain_size;
        let subdomain = Radix2EvaluationDomain::<F>::new(n).unwrap();

        // Extend the computation to the whole domain
        b.resize_with(n, || F::zero());
        // So we don't care about pks, but
        let apk_plus_h_x = acc_x[m];
        let apk_plus_h_y = acc_y[m];
        acc_x.resize_with(n, || apk_plus_h_x);
        acc_y.resize_with(n, || apk_plus_h_y);

        let mut acc_x_shifted = acc_x.clone();
        let mut acc_y_shifted = acc_y.clone();
        acc_x_shifted.rotate_left(1);
        acc_y_shifted.rotate_left(1);

        let mut l1 = vec![F::zero(); n];
        let mut ln = vec![F::zero(); n];
        l1[0] = F::one();
        ln[n-1] = F::one();

        let b_poly = Evaluations::from_vec_and_domain(b, subdomain).interpolate();
        let pks_x_poly = Evaluations::from_vec_and_domain(pks_x, subdomain).interpolate();
        let pks_y_poly = Evaluations::from_vec_and_domain(pks_y, subdomain).interpolate();
        let acc_x_poly = Evaluations::from_vec_and_domain(acc_x, subdomain).interpolate();
        let acc_y_poly = Evaluations::from_vec_and_domain(acc_y, subdomain).interpolate();

        let b_comm = KZG_BW6::commit(&self.kzg_pk, &b_poly);
        let acc_x_comm = KZG_BW6::commit(&self.kzg_pk, &acc_x_poly);
        let acc_y_comm = KZG_BW6::commit(&self.kzg_pk, &acc_y_poly);

        transcript.append_proof_point(b"b_comm", &b_comm);
        transcript.append_proof_point(b"acc_x_comm", &acc_x_comm);
        transcript.append_proof_point(b"acc_y_comm", &acc_y_comm);
        let r = transcript.get_128_bit_challenge(b"r");

        let powers_of_2 = utils::powers(F::from(2u8), 255);
        assert_eq!(self.domain_size % 256, 0); //TODO: 256 is the highest power of 2 that fits field bit capacity
        let powers_of_r = utils::powers(r, self.domain_size / 256 - 1);
        // tensor product (powers_of_r X powers_of_2)
        let c = powers_of_r.iter().flat_map(|rj|
            powers_of_2.iter().map(move |_2k| *rj * _2k)
        ).collect::<Vec<F>>();

        let c_poly = Evaluations::from_vec_and_domain(c, subdomain).interpolate();
        let c_comm = KZG_BW6::commit(&self.kzg_pk, &c_poly);
        transcript.append_proof_point(b"c_comm", &c_comm);
        let phi = transcript.get_128_bit_challenge(b"phi");

        let acc_x_shifted_poly = Evaluations::from_vec_and_domain(acc_x_shifted, subdomain).interpolate();
        let acc_y_shifted_poly = Evaluations::from_vec_and_domain(acc_y_shifted, subdomain).interpolate();
        let l1_poly = Evaluations::from_vec_and_domain(l1, subdomain).interpolate();
        let ln_poly = Evaluations::from_vec_and_domain(ln, subdomain).interpolate();

        assert_eq!(b_poly.coeffs.len(), n);
        assert_eq!(b_poly.degree(), n-1);

        let domain = GeneralEvaluationDomain::<F>::new(4*n).unwrap();
        assert_eq!(domain.size(), 4*n);

        let B = b_poly.evaluate_over_domain_by_ref(domain);
        let x1 = acc_x_poly.evaluate_over_domain_by_ref(domain);
        let y1 = acc_y_poly.evaluate_over_domain_by_ref(domain);
        let x2 = pks_x_poly.evaluate_over_domain_by_ref(domain);
        let y2 = pks_y_poly.evaluate_over_domain_by_ref(domain);
        let x3 = acc_x_shifted_poly.evaluate_over_domain(domain);
        let y3 = acc_y_shifted_poly.evaluate_over_domain(domain);
        let L1 = l1_poly.evaluate_over_domain(domain);
        let Ln = ln_poly.evaluate_over_domain(domain);

        let nB = Evaluations::from_vec_and_domain(
            B.evals.iter().map(|x| F::one() - x).collect(),
            domain
        );

        let a1 =
            &(
                &B *
                    &(
                        &(
                            &(
                                &(&x1 - &x2) * &(&x1 - &x2)
                            ) *
                                &(
                                    &(&x1 + &x2) + &x3
                                )
                        ) -
                            &(
                                &(&y2 - &y1) * &(&y2 - &y1)
                            )
                    )
            ) +
                &(
                    &nB * &(&y3 - &y1)
                );

        let a2 =
            &(
                &B *
                    &(
                        &(
                            &(&x1 - &x2) * &(&y3 + &y1)
                        ) -
                            &(
                                &(&y2 - &y1) * &(&x3 - &x1)
                            )
                    )
            ) +
                &(
                    &nB * &(&x3 - &x1)
                );

        let a3 = &B * &nB;




        let acc_minus_h_x = add_constant(&x1, -self.h.x, domain);
        let acc_minus_h_y = add_constant(&y1, -self.h.y, domain);

        let acc_minus_h_plus_apk_x = add_constant(&x1, -apk_plus_h_x, domain);
        let acc_minus_h_plus_apk_y = add_constant(&y1, -apk_plus_h_y, domain);

        let a4 = &(&acc_minus_h_x * &L1) + &(&acc_minus_h_plus_apk_x * &Ln);
        let a5 = &(&acc_minus_h_y * &L1) + &(&acc_minus_h_plus_apk_y * &Ln);

        let a1_poly = a1.interpolate();
        let a2_poly = a2.interpolate();
        let a3_poly = a3.interpolate();
        let a4_poly = a4.interpolate();
        let a5_poly = a5.interpolate();

        assert_eq!(a1_poly.degree(), 4*(n-1));
        assert_eq!(a2_poly.degree(), 3*(n-1));
        assert_eq!(a3_poly.degree(), 2*(n-1));
        assert_eq!(a4_poly.degree(), 2*(n-1));
        assert_eq!(a5_poly.degree(), 2*(n-1));

        let a1_poly_ = &mul_by_x(&a1_poly) - &mul(subdomain.group_gen_inv, &a1_poly);
        let a2_poly_ = &mul_by_x(&a2_poly) - &mul(subdomain.group_gen_inv, &a2_poly);
        assert_eq!(a1_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
        assert_eq!(a2_poly_.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
        assert_eq!(a3_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
        assert_eq!(a4_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());
        assert_eq!(a5_poly.divide_by_vanishing_poly(subdomain).unwrap().1, DensePolynomial::zero());


        let mut curr = phi;
        let mut powers_of_phi = vec![curr];
        for _ in 0..3 {
            curr *= &phi;
            powers_of_phi.push(curr);
        }

        let mut w = &a1_poly + &mul(powers_of_phi[0], &a2_poly); // a1 + phi a2
        w = &mul_by_x(&w) - &mul(subdomain.group_gen_inv, &w); // X w - omega_inv w = w (X - omega_inv)
        w = &w + &mul(powers_of_phi[1], &a3_poly);
        w = &w + &mul(powers_of_phi[2], &a4_poly);
        w = &w + &mul(powers_of_phi[3], &a5_poly);

        let (q_poly, r) = w.divide_by_vanishing_poly(subdomain).unwrap();
        assert_eq!(r, DensePolynomial::zero());
        assert_eq!(q_poly.degree(), 3*n-3);

        assert_eq!(self.kzg_pk.max_degree(), q_poly.degree()); //TODO: check at the prover creation
        let q_comm = KZG_BW6::commit(&self.kzg_pk, &q_poly);

        transcript.append_proof_point(b"q_comm", &q_comm);
        let zeta = transcript.get_128_bit_challenge(b"zeta");

        let b_zeta = b_poly.evaluate(&zeta);
        let pks_x_zeta = pks_x_poly.evaluate(&zeta);
        let pks_y_zeta = pks_y_poly.evaluate(&zeta);
        let acc_x_zeta = acc_x_poly.evaluate(&zeta);
        let acc_y_zeta = acc_y_poly.evaluate(&zeta);
        let q_zeta = q_poly.evaluate(&zeta);

        let zeta_omega = zeta * subdomain.group_gen;
        let acc_x_zeta_omega = acc_x_poly.evaluate(&zeta_omega);
        let acc_y_zeta_omega = acc_y_poly.evaluate(&zeta_omega);

        transcript.append_proof_scalar(b"b_zeta", &b_zeta);
        transcript.append_proof_scalar(b"pks_x_zeta", &pks_x_zeta);
        transcript.append_proof_scalar(b"pks_y_zeta", &pks_y_zeta);
        transcript.append_proof_scalar(b"acc_x_zeta", &acc_x_zeta);
        transcript.append_proof_scalar(b"acc_y_zeta", &acc_y_zeta);
        transcript.append_proof_scalar(b"q_zeta", &q_zeta);
        let nu: F = transcript.get_128_bit_challenge(b"nu");

        let mut curr = nu;
        let mut powers_of_nu = vec![curr];
        for _ in 0..5 {
            curr *= &nu;
            powers_of_nu.push(curr);
        }

        let w2 = &acc_x_poly + &mul(powers_of_nu[0], &acc_y_poly);
        let w2_proof = KZG_BW6::open(&self.kzg_pk, &w2, zeta_omega);

        let mut w1 = &pks_x_poly + &mul(powers_of_nu[0], &pks_y_poly);
        w1 = &w1 + &mul(powers_of_nu[1], &b_poly);
        w1 = &w1 + &mul(powers_of_nu[2], &q_poly);
        w1 = &w1 + &mul(powers_of_nu[3], &w2);
        let w1_proof = KZG_BW6::open(&self.kzg_pk, &w1, zeta);

        Proof {
            b_comm,
            acc_x_comm,
            acc_y_comm,
            q_comm,

            c_comm,

            w1_proof,
            w2_proof,

            b_zeta,
            pks_x_zeta,
            pks_y_zeta,
            acc_x_zeta,
            acc_y_zeta,
            acc_x_zeta_omega,
            acc_y_zeta_omega,
            q_zeta,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_larger_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<F>::new(4*n).unwrap();

        let p_evals1 = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let p_poly1 = Evaluations::from_vec_and_domain(p_evals1, d1).interpolate();

        let p_evals4 = p_poly1.evaluate_over_domain_by_ref(d4);
        let p_poly4 = p_evals4.interpolate();

        assert_eq!(p_poly1, p_poly4);
    }

    #[test]
    fn test_mul_domain() {
        let rng = &mut test_rng();
        let n = 2;
        let d1 = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let d4 = GeneralEvaluationDomain::<F>::new(4*n).unwrap();

        let a_evals1 = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let b_evals1 = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let a_poly1 = Evaluations::from_vec_and_domain(a_evals1, d1).interpolate();
        let b_poly1 = Evaluations::from_vec_and_domain(b_evals1, d1).interpolate();
        assert_eq!(a_poly1.degree(), n-1);
        assert_eq!(b_poly1.degree(), n-1);

        let a_evals4 = a_poly1.evaluate_over_domain_by_ref(d4);
        let b_evals4 = b_poly1.evaluate_over_domain_by_ref(d4);

        let c_evals4 = &a_evals4 * &b_evals4;
        let c_poly4 = c_evals4.interpolate();

        assert_eq!(c_poly4.degree(), 2*(n-1));
    }

    #[test]
    fn test_shift() {
        let rng = &mut test_rng();
        let n = 2;
        let d = GeneralEvaluationDomain::<F>::new(n).unwrap();

        let p_evals = (0..n).map(|_| F::rand(rng)).collect::<Vec<_>>();
        let mut p_evals_shifted = p_evals.clone();
        p_evals_shifted.rotate_left(1);

        let p = Evaluations::from_vec_and_domain(p_evals, d).interpolate();
        let p_shifted = Evaluations::from_vec_and_domain(p_evals_shifted, d).interpolate();

        if let ark_poly::GeneralEvaluationDomain::Radix2(d) = d {
            let omega = d.group_gen;
            assert_eq!(p.evaluate(&omega), p_shifted.evaluate(&F::one()));
            let x  =  F::rand(rng);
            assert_eq!(p.evaluate(&(x * omega)), p_shifted.evaluate(&x));
        } else {
            assert_eq!(0, 1);
        }
    }
}