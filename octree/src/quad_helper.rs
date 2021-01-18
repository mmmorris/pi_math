use std::mem;

use cgmath::{BaseFloat, Point2, Vector2};
use collision::{Aabb, Aabb2};

pub struct QuadHelper {}

impl QuadHelper {
    pub fn get_deap<S: BaseFloat>(
        d: &mut Vector2<S>,
        loose_layer: usize,
        max_loose: &Vector2<S>,
        deep: usize,
        min_loose: &Vector2<S>,
    ) -> usize {
        let dd = Vector2::new(d.x.to_f64().unwrap(), d.y.to_f64().unwrap());
        let two = S::one() + S::one();
        let x = ((max_loose.x.to_f64().unwrap() / dd.x + 1.0) / 2.0).powf(loose_layer as f64);
        let y = ((max_loose.y.to_f64().unwrap() / dd.y + 1.0) / 2.0).powf(loose_layer as f64);
        d.x = S::from(x * dd.x).unwrap();
        d.y = S::from(y * dd.y).unwrap();
        let deep = if loose_layer < deep {
            // 高于该层的节点，松散值都是用最小值， 也可计算其下每层的八叉节点的大小
            // 八叉节点的大小如果小于最小松散值的2倍， 应该停止向下划分， 因为最小松散值占据了八叉节点的大部分
            // 最大层由设置值和该停止划分的层的最小值
            let mut calc_deep = loose_layer;
            let min = min_loose * two;
            while calc_deep < deep && d.x >= min.x && d.y >= min.y {
                *d = (*d + min_loose) / two;
                calc_deep += 1;
            }
            calc_deep
        } else {
            deep
        };
        deep
    }

    #[inline]
    pub fn smaller_than_min_loose<S: BaseFloat>(d: &Vector2<S>, min_loose: &Vector2<S>) -> bool {
        if d.x <= min_loose.x && d.y <= min_loose.y {
            return true;
        };
        return false;
    }

    #[inline]
    pub fn calc_layer<S: BaseFloat>(loose: &Vector2<S>, el: &Vector2<S>) -> usize {
        let x = if el.x == S::zero() {
            usize::max_value()
        } else {
            (loose.x / el.x).to_usize().unwrap()
        };
        let y = if el.y == S::zero() {
            usize::max_value()
        } else {
            (loose.y / el.y).to_usize().unwrap()
        };

        let min = x.min(y);
        if min == 0 {
            return 0;
        }
        (mem::size_of::<usize>() << 3) - (min.leading_zeros() as usize) - 1
    }

    #[inline]
    pub fn get_child<S: BaseFloat>(point: Point2<S>, aabb: &Aabb2<S>) -> usize {
        let mut i: usize = 0;
        if aabb.max.x > point.x {
            i += 1;
        }
        if aabb.max.y > point.y {
            i += 2;
        }
        i
    }

    #[inline]
    pub fn get_contain_child<S: BaseFloat>(
        parent_aabb: &Aabb2<S>,
        parent_loose: &Vector2<S>,
        aabb: &Aabb2<S>,
    ) -> usize {
        Self::get_child(
            Self::get_aabb_center_add_half_loose(parent_aabb, parent_loose),
            aabb,
        )
    }

    #[inline]
    pub fn get_aabb_center_add_half_loose<S: BaseFloat>(
        aabb: &Aabb2<S>,
        loose: &Vector2<S>,
    ) -> Point2<S> {
        let two = S::one() + S::one();
        let x = (aabb.min.x + aabb.max.x + loose.x) / two;
        let y = (aabb.min.y + aabb.max.y + loose.y) / two;
        Point2::new(x, y)
    }

    #[inline]
    pub fn get_aabb_center_min_half_loose<S: BaseFloat>(
        aabb: &Aabb2<S>,
        loose: &Vector2<S>,
    ) -> Point2<S> {
        let two = S::one() + S::one();
        let x = (aabb.min.x + aabb.max.x - loose.x) / two;
        let y = (aabb.min.y + aabb.max.y - loose.y) / two;
        Point2::new(x, y)
    }

    pub fn get_aabbs_childs<S: BaseFloat>(aabb: &Aabb2<S>, loose: &Vector2<S>) -> Vec<Aabb2<S>> {
        let mut result = Vec::new();
        for index in 0..4 {
            result.push(Self::create_child(aabb, loose, index))
        }

        result
    }

    pub fn create_child<S: BaseFloat>(
        aabb: &Aabb2<S>,
        loose: &Vector2<S>,
        index: usize,
    ) -> Aabb2<S> {
        let p1 = Self::get_aabb_center_min_half_loose(&aabb, &loose);
        let p2 = Self::get_aabb_center_add_half_loose(&aabb, &loose);
        match index {
            0 => Aabb2::new(aabb.min(), p2),
            1 => Aabb2::new(
                Point2::new(p1.x, aabb.min.y),
                Point2::new(aabb.max.x, p2.y),
            ),
            2 => Aabb2::new(
                Point2::new(aabb.min.x, p1.y),
                Point2::new(p2.x, aabb.max.y),
            ),
            _ => Aabb2::new(p1, aabb.max()),
        }
    }
}
