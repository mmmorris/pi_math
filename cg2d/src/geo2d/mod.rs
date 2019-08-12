use num_traits::{Float};

pub use cgmath::{Vector2, Point2};
pub use collision::{Ray2, Line2, Aabb2 as Rectangle};

// 三角形，逆时针的点
#[derive(Clone)]
pub struct Triangle<F: Float> {
    pub vertices: [Point2<F>; 3],
}

/** 
 * 多边形：一堆线段组成的区域
 *    + 简单多边形：全连通；每点2条边；顶点数=边数；任意两条边除了顶点外没有额外的交点；
 *       * 凸多边形。
 *       * 凹多边形：可以看成一堆凸多边形组成，需要找到Convex Decomposition算法
 *    + 弱简单多边形：可以带孔(不全连通)；其他条件同简单多边形
 *    + 非简单多边形：边除了顶点外，还可以自相交。
 */

// 一般多边形
#[derive(Clone)]
pub struct Polygon<F: Float> {
    pub vertices: Vec<Point2<F>>,
    pub hole_indices: Vec<usize>, // 每个洞的起点在vertices数组的索引
}

// 闭合的一堆点集
pub struct LineString<F: Float> {
    pub vertices: Vec<Point2<F>>, 
}

impl<F: Float> Polygon<F> {

    // 设置外围多边形
    pub fn set_exterior(&mut self, string: &LineString<F>) {
        debug_assert!(self.vertices.len() == 0, "invalid time to set");
        self.vertices.extend_from_slice(string.vertices.as_slice());
    }

    // 插入一个孔
    pub fn push_hole(&mut self, string: &LineString<F>) {
        self.hole_indices.push(self.vertices.len());
        self.vertices.extend_from_slice(string.vertices.as_slice());
    }

    // 外围多边形：点
    pub fn exterior(&self) -> PolygonIter<F> {
        let hole_num = self.hole_indices.len();
        if hole_num > 0 {
            PolygonIter::<F>::new(&self.vertices[0..self.hole_indices[0]])
        } else {
            PolygonIter::<F>::new(&self.vertices[..])
        }
    }

    // 孔多边形的数量
    pub fn hole_num(&self) -> usize {
        self.hole_indices.len()
    }

    // 第i个孔多边形的数量
    pub fn hole(&self, i: usize) -> PolygonIter<F> {
        let len = self.hole_indices.len();
        if i >= len {
            PolygonIter::<F>::new(&self.vertices[0..0])
        } else {
            let start = self.hole_indices[i];
            let end = if i + 1 == len {
                self.vertices.len()
            } else {
                self.hole_indices[i + 1]
            };

            PolygonIter::<F>::new(&self.vertices[start..end])
        }
    }
}

impl<F: Float> Default for Polygon<F> {
    fn default() -> Self {
        Self {
            vertices: vec![],
            hole_indices: vec![],
        }
    }
}

impl<F: Float> Default for LineString<F> {
    fn default() -> Self {
        Self {
            vertices: vec![],
        }
    }
}

impl<F: Float> LineString<F> {

    pub fn new(vertices: Vec<Point2<F>>) -> Self {
        Self {
            vertices,
        }
    }

    pub fn push_point(&mut self, vertex: Point2<F>) {
        self.vertices.push(vertex)
    }

    pub fn points(&self) -> PointIterator<F> {
        PointIterator::<F> {
            cursor: 0,
            vertices: &self.vertices,
        }
    }

    pub fn lines(&self) -> LineIterator<F> {
        LineIterator::<F> {
            cursor: 0,
            vertices: &self.vertices
        }
    }

    pub fn triangles(&self) -> TriIterator<F> {
        TriIterator::<F> {
            cursor: 0,
            vertices: &self.vertices
        }
    }
}

pub struct PolygonIter<'a, F: Float> {
    vertices: &'a [Point2<F>],
}

pub struct PointIterator<'a, F: Float> {
    cursor: usize,
    vertices: &'a [Point2<F>],
}

pub struct LineIterator<'a, F: Float> {
    cursor: usize,
    vertices: &'a [Point2<F>],
}

pub struct TriIterator<'a, F: Float> {
    cursor: usize,
    vertices: &'a [Point2<F>],
}

impl<'a, F: Float> PolygonIter<'a, F> {

    pub fn new(vertices: &'a [Point2<F>]) -> Self {
        Self {
            vertices,
        }
    }

    pub fn points(&self) -> PointIterator<F> {
        PointIterator::<F> {
            cursor: 0,
            vertices: &self.vertices,
        }
    }

    pub fn lines(&self) -> LineIterator<F> {
        LineIterator::<F> {
            cursor: 0,
            vertices: &self.vertices,
        }
    }

    pub fn triangles(&self) -> TriIterator<F> {
        TriIterator::<F> {
            cursor: 0,
            vertices: &self.vertices,
        }
    }
}

impl<'a, F: Float> Iterator for PointIterator<'a, F> {
    type Item = &'a Point2<F>;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor < self.vertices.len() {
            self.cursor += 1;
            Some(&self.vertices[self.cursor - 1])   
        } else {
            None
        }
    }
}

impl<'a, F: Float> Iterator for LineIterator<'a, F> {
    type Item = (&'a Point2<F>, &'a Point2<F>);
    
    fn next(&mut self) -> Option<Self::Item> {
        let cursor = self.cursor;
        let len = self.vertices.len();
        if len < 2 || cursor >= len {
            return None;
        }
        
        self.cursor += 1;
        let pts = self.vertices;
        if 1 + cursor < len {
            Some((&pts[cursor], &pts[cursor + 1]))
        } else if len > 2 {
            Some((&pts[cursor], &pts[0]))
        } else {
            None
        }
    }
}

impl<'a, F: Float> Iterator for TriIterator<'a, F> {
    type Item = (&'a Point2<F>, &'a Point2<F>, &'a Point2<F>);
    
    fn next(&mut self) -> Option<Self::Item> {
        let cursor = self.cursor;
        let len = self.vertices.len();
        if len < 3 || cursor >= len {
            return None;
        }

        self.cursor += 1;
        let pts = self.vertices;
        if 2 + cursor < len {
            Some((&pts[cursor], &pts[cursor + 1], &pts[cursor + 2]))
        } else if 2 + cursor == len && len > 3 {
            Some((&pts[cursor], &pts[cursor + 1], &pts[0]))
        } else if 1 + cursor == len && len > 3 {
            Some((&pts[cursor], &pts[0], &pts[1]))
        } else {
            None
        }
    }
}