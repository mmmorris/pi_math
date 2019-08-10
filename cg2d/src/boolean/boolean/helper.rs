use std::cmp::Ordering;

#[inline]
pub fn less_if(condition: bool) -> Ordering {
    if condition {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

#[cfg(test)]
pub mod test {
    use geo2d::Point2;

    pub fn xy<X: Into<f64>, Y: Into<f64>>(x: X, y: Y) -> Point2<f64> {
        Point2 {
            x: x.into(),
            y: y.into(),
        }
    }
}
