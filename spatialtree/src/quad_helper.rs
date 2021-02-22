//! 四叉相关接口
//! 为抽象叉树（四叉或八叉通用一套代码）而抽离的接口struct

use std::mem;

use nalgebra::{Point2, RealField, Scalar, Vector2};
use ncollide2d::bounding_volume::AABB;
use num_traits::{Float, FromPrimitive};
pub struct QuadHelper {}

impl QuadHelper {
    /// 计算四叉树的深度
    pub fn get_deap<S: Scalar + RealField + FromPrimitive>(
        d: &mut Vector2<S>,
        loose_layer: usize,
        max_loose: &Vector2<S>,
        deep: usize,
        min_loose: &Vector2<S>,
    ) -> usize {
        let two = S::one() + S::one();
        let x = ((max_loose.x / d.x + S::one()) / two)
            .powf(FromPrimitive::from_usize(loose_layer).unwrap());
        let y = ((max_loose.y / d.y + S::one()) / two)
            .powf(FromPrimitive::from_usize(loose_layer).unwrap());
        d.x = x;
        d.y = y;
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
    /// 判定指定向量是否小于最小“松散”尺寸
    pub fn smaller_than_min_loose<S: Scalar + RealField>(
        d: &Vector2<S>,
        min_loose: &Vector2<S>,
    ) -> bool {
        if d.x <= min_loose.x && d.y <= min_loose.y {
            return true;
        };
        return false;
    }

    #[inline]
    /// 指定向量以及最大松散尺寸计算对应的层
    pub fn calc_layer<S: Float + Scalar + RealField>(loose: &Vector2<S>, el: &Vector2<S>) -> usize {
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
    /// 判断所在的子节点
    pub fn get_child<S: Scalar + RealField>(point: Point2<S>, aabb: &AABB<S>) -> usize {
        let mut i: usize = 0;
        if aabb.maxs.x > point.x {
            i += 1;
        }
        if aabb.maxs.y > point.y {
            i += 2;
        }
        i
    }

    #[inline]
    ///  判断所在的子节点
    pub fn get_contain_child<S: Scalar + RealField>(
        parent_aabb: &AABB<S>,
        parent_loose: &Vector2<S>,
        aabb: &AABB<S>,
    ) -> usize {
        Self::get_child(
            Self::get_aabb_center_add_half_loose(parent_aabb, parent_loose),
            aabb,
        )
    }

    #[inline]
    ///
    pub fn get_aabb_center_add_half_loose<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector2<S>,
    ) -> Point2<S> {
        let two = S::one() + S::one();
        let x = (aabb.mins.x + aabb.maxs.x + loose.x) / two;
        let y = (aabb.mins.y + aabb.maxs.y + loose.y) / two;
        Point2::new(x, y)
    }

    #[inline]
    ///
    pub fn get_aabb_center_min_half_loose<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector2<S>,
    ) -> Point2<S> {
        let two = S::one() + S::one();
        let x = (aabb.mins.x + aabb.maxs.x - loose.x) / two;
        let y = (aabb.mins.y + aabb.maxs.y - loose.y) / two;
        Point2::new(x, y)
    }

    /// 创建ab的子节点集合
    pub fn get_aabbs_childs<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector2<S>,
    ) -> Vec<AABB<S>> {
        let mut result = Vec::new();
        for index in 0..4 {
            result.push(Self::create_child(aabb, loose, index))
        }

        result
    }

    /// 指定创建ab的子节点
    pub fn create_child<S: Scalar + RealField>(
        aabb: &AABB<S>,
        loose: &Vector2<S>,
        index: usize,
    ) -> AABB<S> {
        let p1 = Self::get_aabb_center_min_half_loose(&aabb, &loose);
        let p2 = Self::get_aabb_center_add_half_loose(&aabb, &loose);
        match index {
            0 => AABB::new(aabb.mins, p2),
            1 => AABB::new(
                Point2::new(p1.x, aabb.mins.y),
                Point2::new(aabb.maxs.x, p2.y),
            ),
            2 => AABB::new(
                Point2::new(aabb.mins.x, p1.y),
                Point2::new(p2.x, aabb.maxs.y),
            ),
            _ => AABB::new(p1, aabb.maxs),
        }
    }
}

#[cfg(test)]
mod quadtests {
    use crate::*;
    use map::vecmap::VecMap;
    use map::Map;
    use nalgebra::{Point2, RealField, Scalar, Vector2};
    use ncollide2d::bounding_volume::{BoundingVolume, AABB};
    use num_traits::Float;
    use rand::Rng;
    use slab::Slab;
    use std::mem;
    use time::Time;
    crate::custom_dimension!(Point2 { x, y }, Vector2 { x, y }, AABB, QuadHelper, 4);

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
                result: Vec::new(),
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
    fn test() {
        let max = Vector2::new(100f32, 100f32);
        let min = max / 100f32;
        let mut tree = Tree::new(
            AABB::new(Point2::new(0f32, 0f32), Point2::new(1000f32, 1000f32)),
            max,
            min,
            6,
            8,
            9999,
        );
        let ab_num = 5000;
        let small_ab_num = (ab_num as f32 * 0.8) as usize;

        for id in 0..small_ab_num as usize {
            let point1 = Point2::new(
                rand::thread_rng().gen_range(0, 999) as f32,
                rand::thread_rng().gen_range(0, 999) as f32,
            );
            let point2 = Point2::new(point1.x + 1.0, point1.y + 1.0);

            tree.add(id + 1, AABB::new(point1, point2), 1);
        }

        for id in small_ab_num + 1..ab_num {
            let point1 = Point2::new(
                rand::thread_rng().gen_range(1, 900) as f32,
                rand::thread_rng().gen_range(1, 900) as f32,
            );
            let point2 = Point2::new(
                point1.x + rand::thread_rng().gen_range(1, 100) as f32,
                point1.y + rand::thread_rng().gen_range(1, 100) as f32,
            );

            tree.add(id + 1, AABB::new(point1, point2), 1);
        }
        tree.collect();

        println!(
            "AABBs Num: {}, ab_1: {:?}",
            ab_num,
            tree.ab_map.get(1).unwrap()
        );

        let point1 = Point2::new(
            rand::thread_rng().gen_range(0, 999) as f32,
            rand::thread_rng().gen_range(0, 999) as f32,
        );
        let point2 = Point2::new(point1.x + 1.0, point1.y + 1.0);
        let start = Time::now(); //获取开始时间
        tree.update(1, AABB::new(point1, point2));
        let end = Time::now();
        println!(
            "update time:{:?}, ab_1: {:?},  ",
            end - start,
            tree.ab_map.get(1).unwrap()
        );

        let point1 = Point2::new(
            rand::thread_rng().gen_range(0, 999) as f32,
            rand::thread_rng().gen_range(0, 999) as f32,
        );
        let point2 = Point2::new(point1.x + 1.0, point1.y + 1.0);
        let start = Time::now(); //获取开始时间
        tree.add(ab_num + 2, AABB::new(point1, point2), 1);
        let end = Time::now();
        println!(
            "add time:{:?}, ab_new: {:?},  ",
            end - start,
            tree.ab_map.get(ab_num + 2).unwrap()
        );

        let start = Time::now(); //获取开始时间
        tree.remove(ab_num + 2);
        let end = Time::now();
        println!("remove time:{:?}  ", end - start);

        let start = Time::now(); //获取开始时间
        tree.collect();
        let end = Time::now();
        println!("collect time:{:?}  ", end - start);

        let point1 = Point2::new(
            rand::thread_rng().gen_range(1, 999) as f32,
            rand::thread_rng().gen_range(1, 999) as f32,
        );
        let point2 = Point2::new(point1.x + 1.0, point1.y + 1.0);
        let aabb = AABB::new(point1, point2);
        let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());
        let start = Time::now(); //获取开始时间
        tree.query(&aabb, intersects, &mut args, ab_query_func);
        let end = Time::now();
        println!("narrow query time:{:?}  ", end - start);

        let point1 = Point2::new(
            rand::thread_rng().gen_range(1, 900) as f32,
            rand::thread_rng().gen_range(1, 900) as f32,
        );
        let point2 = Point2::new(point1.x + 100.0, point1.y + 100.0);
        let aabb = AABB::new(point1, point2);
        let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());
        let start = Time::now(); //获取开始时间
        tree.query(&aabb, intersects, &mut args, ab_query_func);
        let end = Time::now();
        println!("wide query time:{:?}  ", end - start);

    }
}
