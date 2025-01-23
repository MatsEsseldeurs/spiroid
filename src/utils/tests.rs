use super::*;
use pretty_assertions::assert_eq;

#[test]
fn _factorial() {
    let results = [1.0, 1.0, 2.0, 6.0, 24.0];
    for (i, result) in results.iter().enumerate() {
        assert_eq!(factorial(i), *result);
    }
}

#[test]
fn _kronecker_delta() {
    let expected = [1., 0., 0., 1.];
    let result = [
        kronecker_delta(0, 0),
        kronecker_delta(1, 0),
        kronecker_delta(0, 1),
        kronecker_delta(1, 1),
    ];
    assert_eq!(expected, result);
}

#[test]
fn _map_3d_to_1d() {
    let mut result = vec![];
    for q in 0..=14 {
        for p in 0..3 {
            for m in 0..3 {
                result.push(map_3d_to_1d(m, 3, p, 3, q));
            }
        }
    }
    let expected = (0..=134).collect::<Vec<usize>>();
    assert_eq!(expected, result);
}
