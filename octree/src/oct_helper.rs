use std::mem;

use cgmath::{BaseFloat, Point3, Vector3};
use collision::{Aabb, Aabb3};

pub struct OctHelper {}

impl OctHelper {
    pub fn get_deap<S: BaseFloat>(
        d: &mut Vector3<S>,
        loose_layer: usize,
        max_loose: &Vector3<S>,
        deep: usize,
        min_loose: &Vector3<S>,
    ) -> usize {
        let dd = Vector3::new(
            d.x.to_f64().unwrap(),
            d.y.to_f64().unwrap(),
            d.z.to_f64().unwrap(),
        );
        let two = S::one() + S::one();
        let x = ((max_loose.x.to_f64().unwrap() / dd.x + 1.0) / 2.0).powf(loose_layer as f64);
        let y = ((max_loose.y.to_f64().unwrap() / dd.y + 1.0) / 2.0).powf(loose_layer as f64);
        let z = ((max_loose.z.to_f64().unwrap() / dd.z + 1.0) / 2.0).powf(loose_layer as f64);
        d.x = S::from(x * dd.x).unwrap();
        d.y = S::from(y * dd.y).unwrap();
        d.z = S::from(z * dd.z).unwrap();
        let deep = if loose_layer < deep {
            // 高于该层的节点，松散值都是用最小值， 也可计算其下每层的八叉节点的大小
            // 八叉节点的大小如果小于最小松散值的2倍， 应该停止向下划分， 因为最小松散值占据了八叉节点的大部分
            // 最大层由设置值和该停止划分的层的最小值
            let mut calc_deep = loose_layer;
            let min = min_loose * two;
            while calc_deep < deep && d.x >= min.x && d.y >= min.y && d.z >= min.z {
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
    pub fn smaller_than_min_loose<S: BaseFloat>(d: &Vector3<S>, min_loose: &Vector3<S>) -> bool {
        if d.x <= min_loose.x && d.y <= min_loose.y && d.z <= min_loose.z {
            return true;
        };
        return false;
    }

    #[inline]
    pub fn calc_layer<S: BaseFloat>(loose: &Vector3<S>, el: &Vector3<S>) -> usize {
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
        let z = if el.z == S::zero() {
            usize::max_value()
        } else {
            (loose.z / el.z).to_usize().unwrap()
        };
        let min = x.min(y).min(z);
        if min == 0 {
            return 0;
        }
        (mem::size_of::<usize>() << 3) - (min.leading_zeros() as usize) - 1
    }

    #[inline]
    pub fn get_child<S: BaseFloat>(point: Point3<S>, aabb: &Aabb3<S>) -> usize {
        let mut i: usize = 0;
        if aabb.max.x > point.x {
            i += 1;
        }
        if aabb.max.y > point.y {
            i += 2;
        }
        if aabb.max.z > point.z {
            i += 4;
        }
        i
    }

    #[inline]
    pub fn get_contain_child<S: BaseFloat>(
        parent_aabb: &Aabb3<S>,
        parent_loose: &Vector3<S>,
        aabb: &Aabb3<S>,
    ) -> usize {
        Self::get_child(
            Self::get_aabb_center_add_half_loose(parent_aabb, parent_loose),
            aabb,
        )
    }

    #[inline]
    pub fn get_aabb_center_add_half_loose<S: BaseFloat>(
        aabb: &Aabb3<S>,
        loose: &Vector3<S>,
    ) -> Point3<S> {
        let two = S::one() + S::one();
        let x = (aabb.min.x + aabb.max.x + loose.x) / two;
        let y = (aabb.min.y + aabb.max.y + loose.y) / two;
        let z = (aabb.min.z + aabb.max.z + loose.z) / two;
        Point3::new(x, y, z)
    }

    #[inline]
    pub fn get_aabb_center_min_half_loose<S: BaseFloat>(
        aabb: &Aabb3<S>,
        loose: &Vector3<S>,
    ) -> Point3<S> {
        let two = S::one() + S::one();
        let x = (aabb.min.x + aabb.max.x - loose.x) / two;
        let y = (aabb.min.y + aabb.max.y - loose.y) / two;
        let z = (aabb.min.z + aabb.max.z - loose.z) / two;
        Point3::new(x, y, z)
    }

    pub fn get_aabbs_childs<S: BaseFloat>(aabb: &Aabb3<S>, loose: &Vector3<S>) -> Vec<Aabb3<S>> {
        let mut result = Vec::new();
        for index in 0..8 {
            result.push(Self::create_child(aabb, loose, index))
        }

        result
    }

    pub fn create_child<S: BaseFloat>(
        aabb: &Aabb3<S>,
        loose: &Vector3<S>,
        index: usize,
    ) -> Aabb3<S> {
        let p1 = Self::get_aabb_center_min_half_loose(&aabb, &loose);
        let p2 = Self::get_aabb_center_add_half_loose(&aabb, &loose);
        match index {
            0 => Aabb3::new(aabb.min(), p2),
            1 => Aabb3::new(
                Point3::new(p1.x, aabb.min.y, aabb.min.z),
                Point3::new(aabb.max.x, p2.y, p2.z),
            ),
            2 => Aabb3::new(
                Point3::new(aabb.min.x, p1.y, aabb.min.z),
                Point3::new(p2.x, aabb.max.y, p2.z),
            ),
            3 => Aabb3::new(
                Point3::new(p1.x, p1.y, aabb.min.z),
                Point3::new(aabb.max.x, aabb.max.y, p2.z),
            ),
            4 => Aabb3::new(
                Point3::new(aabb.min.x, aabb.min.y, p1.z),
                Point3::new(p2.x, p2.y, aabb.max.z),
            ),
            5 => Aabb3::new(
                Point3::new(p1.x, aabb.min.y, p1.z),
                Point3::new(aabb.max.x, p2.y, aabb.max.z),
            ),
            6 => Aabb3::new(
                Point3::new(aabb.min.x, p1.y, p1.z),
                Point3::new(p2.x, aabb.max.y, aabb.max.z),
            ),
            _ => Aabb3::new(p1, aabb.max()),
        }
    }
}
