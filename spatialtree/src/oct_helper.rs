//! 四叉相关接口
//! 为抽象叉树（四叉或八叉通用一套代码）而抽离的接口struct
use std::mem;

use nalgebra::{Point3, RealField, Scalar, Vector3};
use ncollide3d::bounding_volume::AABB;
use num_traits::{Float, FromPrimitive};
pub struct OctHelper {}

impl OctHelper {
    /// 计算八叉树的深度
    pub fn get_deap<S: Scalar + RealField + FromPrimitive>(
        d: &mut Vector3<S>,
        loose_layer: usize,
        max_loose: &Vector3<S>,
        deep: usize,
        min_loose: &Vector3<S>,
    ) -> usize {
        let two = S::one() + S::one();
        let x = ((max_loose.x / d.x + S::one()) / two)
            .powf(FromPrimitive::from_usize(loose_layer).unwrap());
        let y = ((max_loose.y / d.y + S::one()) / two)
            .powf(FromPrimitive::from_usize(loose_layer).unwrap());
        let z = ((max_loose.z / d.z + S::one()) / two)
            .powf(FromPrimitive::from_usize(loose_layer).unwrap());
        d.x *= x;
        d.y *= y;
        d.z *= z;
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
    /// 判定指定向量是否小于最小“松散”尺寸
    pub fn smaller_than_min_loose<S: Scalar + RealField>(
        d: &Vector3<S>,
        min_loose: &Vector3<S>,
    ) -> bool {
        if d.x <= min_loose.x && d.y <= min_loose.y && d.z <= min_loose.z {
            return true;
        };
        return false;
    }

    #[inline]
    /// 指定向量以及最大松散尺寸计算对应的层
    pub fn calc_layer<S: Float + Scalar + RealField>(loose: &Vector3<S>, el: &Vector3<S>) -> usize {
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
    /// 判断所在的子节点
    pub fn get_child<S: Scalar + RealField>(point: Point3<S>, aabb: &AABB<S>) -> usize {
        let mut i: usize = 0;
        if aabb.maxs.x > point.x {
            i += 1;
        }
        if aabb.maxs.y > point.y {
            i += 2;
        }
        if aabb.maxs.z > point.z {
            i += 4;
        }
        i
    }

    #[inline]
    ///  判断所在的子节点
    pub fn get_contain_child<S: Scalar + RealField>(
        parent_aabb: &AABB<S>,
        parent_loose: &Vector3<S>,
        aabb: &AABB<S>,
    ) -> usize {
        Self::get_child(
            Self::get_aabb_center_add_half_loose(parent_aabb, parent_loose),
            aabb,
        )
    }

    #[inline]
    pub fn get_aabb_center_add_half_loose<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector3<S>,
    ) -> Point3<S> {
        let two = S::one() + S::one();
        let x = (aabb.mins.x + aabb.maxs.x + loose.x) / two;
        let y = (aabb.mins.y + aabb.maxs.y + loose.y) / two;
        let z = (aabb.mins.z + aabb.maxs.z + loose.z) / two;
        Point3::new(x, y, z)
    }

    #[inline]
    pub fn get_aabb_center_min_half_loose<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector3<S>,
    ) -> Point3<S> {
        let two = S::one() + S::one();
        let x = (aabb.mins.x + aabb.maxs.x - loose.x) / two;
        let y = (aabb.mins.y + aabb.maxs.y - loose.y) / two;
        let z = (aabb.mins.z + aabb.maxs.z - loose.z) / two;
        Point3::new(x, y, z)
    }

    /// 创建ab的子节点集合
    pub fn get_aabbs_childs<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector3<S>,
    ) -> Vec<AABB<S>> {
        let mut result = Vec::new();
        for index in 0..8 {
            result.push(Self::create_child(aabb, loose, index))
        }

        result
    }

    /// 指定创建ab的子节点
    pub fn create_child<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector3<S>,
        index: usize,
    ) -> AABB<S> {
        let p1 = Self::get_aabb_center_min_half_loose(&aabb, &loose);
        let p2 = Self::get_aabb_center_add_half_loose(&aabb, &loose);
        match index {
            0 => AABB::new(aabb.mins, p2),
            1 => AABB::new(
                Point3::new(p1.x, aabb.mins.y, aabb.mins.z),
                Point3::new(aabb.maxs.x, p2.y, p2.z),
            ),
            2 => AABB::new(
                Point3::new(aabb.mins.x, p1.y, aabb.mins.z),
                Point3::new(p2.x, aabb.maxs.y, p2.z),
            ),
            3 => AABB::new(
                Point3::new(p1.x, p1.y, aabb.mins.z),
                Point3::new(aabb.maxs.x, aabb.maxs.y, p2.z),
            ),
            4 => AABB::new(
                Point3::new(aabb.mins.x, aabb.mins.y, p1.z),
                Point3::new(p2.x, p2.y, aabb.maxs.z),
            ),
            5 => AABB::new(
                Point3::new(p1.x, aabb.mins.y, p1.z),
                Point3::new(aabb.maxs.x, p2.y, aabb.maxs.z),
            ),
            6 => AABB::new(
                Point3::new(aabb.mins.x, p1.y, p1.z),
                Point3::new(p2.x, aabb.maxs.y, aabb.maxs.z),
            ),
            _ => AABB::new(p1, aabb.maxs),
        }
    }
}

#[cfg(test)]
mod octtests {
    use crate::*;
    use map::vecmap::VecMap;
    use map::Map;
    use nalgebra::{Point3, RealField, Scalar, Vector3};
    use ncollide3d::bounding_volume::{BoundingVolume, AABB};
    use num_traits::Float;
    use rand::Rng;
    use slab::Slab;
    use std::mem;
    use time::Time;
    crate::custom_dimension!(Point3 { x, y, z }, Vector3 { x, y, z }, AABB, OctHelper, 8);

    struct AbQueryArgs<S: Scalar + RealField + Float, T> {
        aabb: AABB<S>,
        result: Vec<(usize, T)>,
    }

    #[inline]
    fn intersects<S: Scalar + RealField + Float>(a: &AABB<S>, b: &AABB<S>) -> bool {
        a.intersects(b)
    }
    impl<S: Scalar + RealField + Float, T: Clone> AbQueryArgs<S, T> {
        fn new(aabb: AABB<S>) -> AbQueryArgs<S, T> {
            AbQueryArgs {
                aabb: aabb,
                result: Vec::with_capacity(100),
            }
        }

        fn result(&mut self) -> Vec<(usize, T)> {
            mem::replace(&mut self.result, Vec::new())
        }
    }

    fn ab_query_func<S: Scalar + RealField + Float, T: Clone>(
        arg: &mut AbQueryArgs<S, T>,
        id: usize,
        aabb: &AABB<S>,
        bind: &T,
    ) {
        if intersects(&arg.aabb, aabb) {
            arg.result.push((id, bind.clone()));
        }
    }

    #[test]
    fn test_oct() {
        let max = Vector3::new(100f32, 100f32, 100f32);
        let min = max / 100f32;
        let mut tree = Tree::new(
            AABB::new(
                Point3::new(0f32, 0f32, 0f32),
                Point3::new(1000f32, 1000f32, 1000f32),
            ),
            
            max,
            min,
            6,
            8,
            50,
        );
        let ab_num = 100000;
        let query_times = 100;
        let small_ab_num = (ab_num as f32 * 0.8) as usize;

        let mut aabbs = Vec::with_capacity(ab_num);
        for id in 0..small_ab_num as usize {
            let point1 = Point3::new(
                rand::thread_rng().gen_range(0, 999) as f32,
                rand::thread_rng().gen_range(0, 999) as f32,
                rand::thread_rng().gen_range(0, 999) as f32,
            );
            let point2 = Point3::new(point1.x + 1.0, point1.y + 1.0, point1.z + 1.0);
            aabbs.push(AABB::new(point1, point2));
            // tree.add(id + 1, AABB::new(point1, point2), 1);
        }

        for id in small_ab_num..ab_num {
            let point1 = Point3::new(
                rand::thread_rng().gen_range(0, 900) as f32,
                rand::thread_rng().gen_range(0, 900) as f32,
                rand::thread_rng().gen_range(0, 900) as f32,
            );
            let point2 = Point3::new(
                point1.x + rand::thread_rng().gen_range(1, 100) as f32,
                point1.y + rand::thread_rng().gen_range(1, 100) as f32,
                point1.z + rand::thread_rng().gen_range(1, 100) as f32,
            );
            aabbs.push(AABB::new(point1, point2));
            // tree.add(id + 1, AABB::new(point1, point2), 1);
        }
        let start = Time::now(); //获取开始时间
        for idx in 0..ab_num {
            tree.add(idx + 1, aabbs[idx], 1);
        }
        let end = Time::now();
        println!("insert time:{:?},  ", end - start,);

        let start = Time::now(); //获取开始时间
        tree.collect();
        let end = Time::now();
        println!("collect time:{:?}  ", end - start);

        let point1 = Point3::new(500.0, 500.0, 500.0);
        let point2 = Point3::new(point1.x + 1.0, point1.y + 1.0, point1.z + 1.0);
        let aabb = AABB::new(point1, point2);
        let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());
        let start = Time::now(); //获取开始时间
        for _ in 0..query_times {
            args.result.clear();
            tree.query(&aabb, intersects, &mut args, ab_query_func);
        }
        let end = Time::now();
        println!("narrow query time:{:?}  ,find num: {}", end - start, args.result.len());

        let point1 = Point3::new(500.0, 500.0, 500.0);
        let point2 = Point3::new(point1.x + 100.0, point1.y + 100.0, point1.z + 100.0);
        let aabb = AABB::new(point1, point2);
        let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());
        let start = Time::now(); //获取开始时间
        for _ in 0..query_times {
            args.result.clear();
            tree.query(&aabb, intersects, &mut args, ab_query_func);
        }
        let end = Time::now();
        println!("wide query time:{:?}  ,find num: {}", end - start, args.result.len());

        aabbs.clear();
        for id in 0..small_ab_num as usize {
            let point1 = Point3::new(
                rand::thread_rng().gen_range(0, 999) as f32,
                rand::thread_rng().gen_range(0, 999) as f32,
                rand::thread_rng().gen_range(0, 999) as f32,
            );
            let point2 = Point3::new(point1.x + 1.0, point1.y + 1.0, point1.z + 1.0);
            aabbs.push(AABB::new(point1, point2));
            // tree.add(id + 1, AABB::new(point1, point2), 1);
        }
        
        let start = Time::now(); //获取开始时间
        for idx in 0..small_ab_num {
            tree.update(idx + 1, aabbs[idx]);
        }
        let end = Time::now();
        println!(
            "update time:{:?},  ",
            end - start,
        );

        let start = Time::now(); //获取开始时间
        for id in 0..ab_num {
            tree.remove(id + 1);
        }
        
        let end = Time::now();
        println!("remove time:{:?}  ", end - start);

        let start = Time::now(); //获取开始时间
        tree.collect();
        let end = Time::now();
        println!("collect time:{:?}  ", end - start);
    }
}
