use num_traits::{Float};

pub use cgmath::{Vector2, Point2};
pub use collision::{Ray2, Line2, Aabb2 as Rectangle};

// 三角形，逆时针的点
#[derive(Clone)]
pub struct Triangle<S: Float> {
    pub vertices: [Point2<S>; 3],
}

// 一般多边形
#[derive(Clone)]
pub struct CommonPolygon<S: Float> {
    pub vertices: Vec<Point2<S>>,
    pub hole_indices: Vec<usize>, // 每个洞的起点在vertices数组的索引
}

// 凸多边形
#[derive(Clone)]
pub struct ConvexPolygon<S: Float> {
    pub vertices: Vec<Point2<S>>,
}

// 凹多边形
// TODO: 需要找凸分割算法。
// TODO：需要指明哪条边是内边
#[derive(Clone)]
pub struct ConcavePolygon<S: Float> {
    pub vertices: Vec<Point2<S>>,
    pub converx_indices: Vec<Vec<usize>>,
}

// 简单多边形？
#[derive(Clone)]
pub struct SimplePolygon<S: Float> {
    pub vertices: Vec<Point2<S>>,
}

// 多边形
#[derive(Clone)]
pub enum Polygon<S: Float> {
    Convex(ConvexPolygon<S>),
    Concave(ConcavePolygon<S>),
    Simple(SimplePolygon<S>),
    Common(CommonPolygon<S>),
}