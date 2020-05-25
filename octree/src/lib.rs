//! 高性能的松散八叉树
extern crate core;

extern crate cgmath;
extern crate collision;
extern crate map;
extern crate slab;
use std::mem;

use cgmath::{BaseNum, Point3, Vector3};
use collision::{Aabb, Aabb3, Contains};

use map::vecmap::VecMap;
use map::Map;
use slab::Slab;

// aabb是否相交
#[inline]
pub fn intersects<S: BaseNum>(a: &Aabb3<S>, b: &Aabb3<S>) -> bool {
    a.min.x < b.max.x
        && a.max.x > b.min.x
        && a.min.y < b.max.y
        && a.max.y > b.min.y
        && a.min.z < b.max.z
        && a.max.z > b.min.z
}

/// aabb的查询函数的参数
pub struct AbQueryArgs<S: BaseNum, T> {
    aabb: Aabb3<S>,
    result: Vec<(usize, T)>,
}
impl<S: BaseNum, T: Clone> AbQueryArgs<S, T> {
    pub fn new(aabb: Aabb3<S>) -> AbQueryArgs<S, T> {
        AbQueryArgs {
            aabb: aabb,
            result: Vec::new(),
        }
    }
    pub fn result(&mut self) -> Vec<(usize, T)> {
        mem::replace(&mut self.result, Vec::new())
    }
}

/// aabb的ab查询函数, aabb的oct查询函数应该使用intersects
pub fn ab_query_func<S: BaseNum, T: Clone>(
    arg: &mut AbQueryArgs<S, T>,
    id: usize,
    aabb: &Aabb3<S>,
    bind: &T,
) {
    if intersects(&arg.aabb, aabb) {
        arg.result.push((id, bind.clone()));
    }
}

/// OctTree
pub struct Tree<S: BaseNum, T> {
    oct_slab: Slab<OctNode<S>>,
    ab_map: VecMap<AbNode<S, T>>,
    loose_ratio: usize,                     //松散系数，0-10000之间， 默认3000
    adjust: (usize, usize),                 //小于min，节点收缩; 大于max，节点分化。默认(4, 5)
    deep: usize,                            // 最大深度
    loose: Vector3<S>,                      //第一层的松散大小
    outer: NodeList, // 和根节点不相交的ab节点列表，及节点数量。 相交的放在root的nodes上了。 该AbNode的parent为0
    dirty: (Vec<Vec<usize>>, usize, usize), // 脏的OctNode节点, 及脏节点数量，及脏节点的起始层
}

impl<S: BaseNum, T> Tree<S, T> {
    pub fn new(
        root: Aabb3<S>,
        loose_ratio: usize,
        adjust_min: usize,
        adjust_max: usize,
        deep: usize,
    ) -> Tree<S, T> {
        let d = root.dim();
        let loose_ratio = if loose_ratio == 0 { LOOSE } else { loose_ratio };
        let loose_ratio = if loose_ratio > LOOSE_MAX {
            LOOSE_MAX
        } else {
            loose_ratio
        };
        let lr = S::from(loose_ratio).unwrap();
        let w = S::from(10000).unwrap();
        let loose = Vector3::new(d.x * lr / w, d.y * lr / w, d.z * lr / w);
        let mut s = Slab::new();
        s.insert(OctNode::new(root, loose.clone(), 0, 0, 0));
        let deep = if deep == 0 || deep > DEEP_MAX {
            DEEP_MAX
        } else {
            deep
        };
        let two = S::one() + S::one();
        let deep = if S::from(1usize).unwrap() / two == S::zero() {
            // 如果是整数空间，则必须计算最大深度，否则会出现物件不在子节点上
            let mut i = 0;
            let mut l = loose / two;
            while i < deep && l.x > S::zero() && l.y > S::zero() && l.z > S::zero() {
                l = l / two;
                i += 1;
            }
            i
        } else {
            deep
        };
        let adjust_min = if adjust_min == 0 {
            ADJUST_MIN
        } else {
            adjust_min
        };
        let adjust_max = if adjust_max == 0 {
            ADJUST_MAX
        } else {
            adjust_max
        };
        let adjust_max = if adjust_min > adjust_max {
            adjust_min
        } else {
            adjust_max
        };
        Tree {
            oct_slab: s,
            ab_map: VecMap::default(),
            loose_ratio: loose_ratio,
            adjust: (adjust_min, adjust_max),
            deep: deep,
            loose: loose,
            outer: NodeList::new(),
            dirty: (Vec::new(), 0, usize::max_value()),
        }
    }
    pub fn mem_size(&self) -> usize {
        let mut r = self.oct_slab.mem_size()
            + self.ab_map.mem_size()
            + self.outer.len() * std::mem::size_of::<usize>();
        for v in self.dirty.0.iter() {
            r += v.capacity() * std::mem::size_of::<usize>();
        }
        r
    }
    // 获得松散系数
    pub fn get_loose_ratio(&self) -> usize {
        self.loose_ratio
    }
    // 获得节点收缩和分化的阈值
    pub fn get_adjust(&self) -> (usize, usize) {
        (self.adjust.0, self.adjust.1)
    }
    // 获得该aabb对应的层
    pub fn get_layer(&self, aabb: &Aabb3<S>) -> usize {
        calc_layer(&self.loose, &aabb.dim())
    }
    // 添加一个aabb及其绑定
    pub fn add(&mut self, id: usize, aabb: Aabb3<S>, bind: T) {
        let layer = calc_layer(&self.loose, &aabb.dim());
        match self.ab_map.insert(id, AbNode::new(aabb, bind, layer)) {
            Some(_) => panic!("duplicate id: {}", id),
            _ => (),
        }
        let next = {
            let node = unsafe { self.ab_map.get_unchecked_mut(id) };
            let root = unsafe { self.oct_slab.get_unchecked_mut(1) };
            if root.aabb.contains(&node.aabb) {
                set_tree_dirty(
                    &mut self.dirty,
                    down(&mut self.oct_slab, self.adjust.1, self.deep, 1, node, id),
                );
            } else if intersects(&root.aabb, &node.aabb) {
                // 相交的放在root的nodes上
                node.parent = 1;
                node.next = root.nodes.head;
                root.nodes.push(id);
            } else {
                // 和根节点不相交的ab节点, 该AbNode的parent为0
                node.next = self.outer.head;
                self.outer.push(id);
            }
            node.next
        };
        if next > 0 {
            let n = unsafe { self.ab_map.get_unchecked_mut(next) };
            n.prev = id;
        }
    }
    // 获取指定id的aabb及其绑定
    pub fn get(&self, id: usize) -> Option<(&Aabb3<S>, &T)> {
        match self.ab_map.get(id) {
            Some(node) => Some((&node.aabb, &node.bind)),
            _ => None,
        }
    }
    // 获取指定id的aabb及其绑定
    pub unsafe fn get_unchecked(&self, id: usize) -> (&Aabb3<S>, &T) {
        let node = self.ab_map.get_unchecked(id);
        (&node.aabb, &node.bind)
    }
    // 更新指定id的aabb
    pub fn update(&mut self, id: usize, aabb: Aabb3<S>) -> bool {
        let r = match self.ab_map.get_mut(id) {
            Some(node) => {
                node.layer = calc_layer(&self.loose, &aabb.dim());
                node.aabb = aabb;
                update(
                    &mut self.oct_slab,
                    &self.adjust,
                    self.deep,
                    &mut self.outer,
                    &mut self.dirty,
                    id,
                    node,
                )
            }
            _ => return false,
        };
        remove_add(self, id, r);
        true
    }
    // 移动指定id的aabb，性能比update要略好
    pub fn shift(&mut self, id: usize, distance: Vector3<S>) -> bool {
        let r = match self.ab_map.get_mut(id) {
            Some(node) => {
                node.aabb = Aabb3::new(node.aabb.min + distance, node.aabb.max + distance);
                update(
                    &mut self.oct_slab,
                    &self.adjust,
                    self.deep,
                    &mut self.outer,
                    &mut self.dirty,
                    id,
                    node,
                )
            }
            _ => return false,
        };
        remove_add(self, id, r);
        true
    }
    // 获取指定id的可写绑定
    pub unsafe fn get_mut(&mut self, id: usize) -> Option<&mut T> {
        match self.ab_map.get_mut(id) {
            Some(n) => Some(&mut n.bind),
            _ => None,
        }
    }
    // 获取指定id的可写绑定
    pub unsafe fn get_unchecked_mut(&mut self, id: usize) -> &mut T {
        let node = self.ab_map.get_unchecked_mut(id);
        &mut node.bind
    }
    // 更新指定id的绑定
    pub fn update_bind(&mut self, id: usize, bind: T) -> bool {
        match self.ab_map.get_mut(id) {
            Some(node) => {
                node.bind = bind;
                true
            }
            _ => false,
        }
    }
    // 移除指定id的aabb及其绑定
    pub fn remove(&mut self, id: usize) -> Option<(Aabb3<S>, T)> {
        let node = match self.ab_map.remove(id) {
            Some(n) => n,
            _ => return None,
        };
        if node.parent > 0 {
            let (p, c) = {
                let parent = unsafe { self.oct_slab.get_unchecked_mut(node.parent) };
                if node.parent_child < 8 {
                    // 在节点的childs上
                    match parent.childs[node.parent_child] {
                        ChildNode::Ab(ref mut ab) => {
                            ab.remove(&mut self.ab_map, node.prev, node.next)
                        }
                        _ => panic!("invalid state"),
                    }
                } else {
                    // 在节点的nodes上
                    parent.nodes.remove(&mut self.ab_map, node.prev, node.next);
                }
                (parent.parent, parent.parent_child)
            };
            remove_up(&mut self.oct_slab, self.adjust.0, &mut self.dirty, p, c);
        } else {
            // 表示在outer上
            self.outer.remove(&mut self.ab_map, node.prev, node.next);
        }
        Some((node.aabb, node.bind))
    }
    // 整理方法，只有整理方法才会创建或销毁OctNode
    pub fn collect(&mut self) {
        let mut count = self.dirty.1;
        if count == 0 {
            return;
        }
        for i in self.dirty.2..self.dirty.0.len() {
            let vec = unsafe { self.dirty.0.get_unchecked_mut(i) };
            let c = vec.len();
            if c == 0 {
                continue;
            }
            for j in 0..c {
                let oct_id = unsafe { vec.get_unchecked(j) };
                collect(
                    &mut self.oct_slab,
                    &mut self.ab_map,
                    &self.adjust,
                    self.deep,
                    *oct_id,
                );
            }
            vec.clear();
            if count <= c {
                break;
            }
            count -= c;
        }
        self.dirty.1 = 0;
        self.dirty.2 = usize::max_value();
    }

    // 查询空间内及相交的ab节点
    pub fn query<A, B>(
        &self,
        oct_arg: &A,
        oct_func: fn(arg: &A, aabb: &Aabb3<S>) -> bool,
        ab_arg: &mut B,
        ab_func: fn(arg: &mut B, id: usize, aabb: &Aabb3<S>, bind: &T),
    ) {
        query(
            &self.oct_slab,
            &self.ab_map,
            1,
            oct_arg,
            oct_func,
            ab_arg,
            ab_func,
        )
    }
    // 查询空间外的ab节点
    pub fn query_outer<B>(
        &self,
        arg: &mut B,
        func: fn(arg: &mut B, id: usize, aabb: &Aabb3<S>, bind: &T),
    ) {
        let mut id = self.outer.head;
        while id > 0 {
            let ab = unsafe { self.ab_map.get_unchecked(id) };
            func(arg, id, &ab.aabb, &ab.bind);
            id = ab.next;
        }
    }

    // 检查碰撞对，不会检查outer的aabb。一般arg包含1个hashset，用(big, little)做键，判断是否已经计算过。
    pub fn collision<A>(
        &self,
        id: usize,
        _limit_layer: usize,
        arg: &mut A,
        func: fn(
            arg: &mut A,
            a_id: usize,
            a_aabb: &Aabb3<S>,
            a_bind: &T,
            b_id: usize,
            b_aabb: &Aabb3<S>,
            b_bind: &T,
        ) -> bool,
    ) {
        let a = match self.ab_map.get(id) {
            Some(ab) => ab,
            _ => return,
        };
        // 先判断root.nodes是否有节点，如果有则遍历root的nodes
        let node = unsafe { self.oct_slab.get_unchecked(1) };
        collision_list(
            &self.ab_map,
            id,
            &a.aabb,
            &a.bind,
            arg,
            func,
            node.nodes.head,
        );
        // 和同列表节点碰撞
        collision_list(&self.ab_map, id, &a.aabb, &a.bind, arg, func, a.next);
        let mut prev = a.prev;
        while prev > 0 {
            let b = unsafe { self.ab_map.get_unchecked(prev) };
            func(arg, id, &a.aabb, &a.bind, prev, &b.aabb, &b.bind);
            prev = b.prev;
        }
        // 需要计算是否在重叠区，如果在，则需要上溯检查重叠的兄弟节点。不在，其实也需要上溯检查父的匹配节点，但可以提前计算ab节点的最小层
        //}
    }
}

#[derive(Debug, Clone)]
struct NodeList {
    head: usize,
    len: usize,
}
impl NodeList {
    #[inline]
    pub fn new() -> NodeList {
        NodeList { head: 0, len: 0 }
    }
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }
    #[inline]
    pub fn push(&mut self, id: usize) {
        self.head = id;
        self.len += 1;
    }
    #[inline]
    pub fn remove<S: BaseNum, T>(
        &mut self,
        map: &mut VecMap<AbNode<S, T>>,
        prev: usize,
        next: usize,
    ) {
        if prev > 0 {
            let node = unsafe { map.get_unchecked_mut(prev) };
            node.next = next;
        } else {
            self.head = next;
        }
        if next > 0 {
            let node = unsafe { map.get_unchecked_mut(next) };
            node.prev = prev;
        }
        self.len -= 1;
    }
}

const LOOSE: usize = 3000;
const LOOSE_MAX: usize = 5000;
const DEEP_MAX: usize = 24;
const ADJUST_MIN: usize = 4;
const ADJUST_MAX: usize = 5;

#[derive(Debug, Clone)]
struct OctNode<S: BaseNum> {
    aabb: Aabb3<S>,         // 包围盒
    loose: Vector3<S>,      // 本层的松散值
    parent: usize,          // 父八叉节点
    parent_child: usize,    // 对应父八叉节点childs的位置
    childs: [ChildNode; 8], // 子八叉节点
    layer: usize,           // 表示第几层， 根据aabb大小，决定最低为第几层
    nodes: NodeList,        // 匹配本层大小的ab节点列表，及节点数量
    dirty: usize, // 脏标记, 1-128对应节点被修改。添加了节点，并且某个子八叉节点(AbNode)的数量超过阈值，可能分化。删除了节点，并且自己及其下ab节点的数量超过阈值，可能收缩
}
impl<S: BaseNum> OctNode<S> {
    #[inline]
    pub fn new(
        aabb: Aabb3<S>,
        loose: Vector3<S>,
        parent: usize,
        child: usize,
        layer: usize,
    ) -> OctNode<S> {
        OctNode {
            aabb: aabb,
            loose: loose,
            parent: parent,
            parent_child: child,
            childs: [
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
                ChildNode::Ab(NodeList::new()),
            ],
            layer: layer,
            nodes: NodeList::new(),
            dirty: 0,
        }
    }
}
#[derive(Debug, Clone)]
enum ChildNode {
    Oct(usize, usize), // 对应的OctNode, 及其下ab节点的数量
    Ab(NodeList),      // ab节点列表，及节点数量
}

#[derive(Debug, Clone)]
struct AbNode<S: BaseNum, T> {
    aabb: Aabb3<S>,      // 包围盒
    bind: T,             // 绑定
    layer: usize,        // 表示第几层， 根据aabb大小，决定最低为第几层
    parent: usize,       // 父八叉节点
    parent_child: usize, // 父八叉节点所在的子八叉节点， 8表示不在子八叉节点上
    prev: usize,         // 前ab节点
    next: usize,         // 后ab节点
}
impl<S: BaseNum, T> AbNode<S, T> {
    pub fn new(aabb: Aabb3<S>, bind: T, layer: usize) -> AbNode<S, T> {
        AbNode {
            aabb: aabb,
            bind: bind,
            layer: layer,
            parent: 0,
            parent_child: 8,
            prev: 0,
            next: 0,
        }
    }
}

// 计算该aabb对应的层
#[inline]
fn calc_layer<S: BaseNum>(loose: &Vector3<S>, el: &Vector3<S>) -> usize {
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
    (mem::size_of::<usize>() << 3) - (min.leading_zeros() as usize)
}
// 检查自己、父亲及兄弟节点
fn check_contain<S: BaseNum>(
    parent: &Aabb3<S>,
    loose: &Vector3<S>,
    node: &Aabb3<S>,
    child: usize,
) -> usize {
    let two = S::one() + S::one();
    let x1 = (parent.min.x + parent.max.x - loose.x) / two;
    let y1 = (parent.min.y + parent.max.y - loose.y) / two;
    let z1 = (parent.min.z + parent.max.z - loose.z) / two;
    let x2 = (parent.min.x + parent.max.x + loose.x) / two;
    let y2 = (parent.min.y + parent.max.y + loose.y) / two;
    let z2 = (parent.min.z + parent.max.z + loose.z) / two;
    let a = match child {
        0 => Aabb3::new(Point3::new(x1, y1, z1), parent.max()),
        1 => Aabb3::new(
            Point3::new(parent.min.x, y1, z1),
            Point3::new(x2, parent.max.y, parent.max.z),
        ),
        2 => Aabb3::new(
            Point3::new(x1, parent.min.y, z1),
            Point3::new(parent.max.x, y2, parent.max.z),
        ),
        3 => Aabb3::new(
            Point3::new(x1, y1, parent.min.z),
            Point3::new(parent.max.x, parent.max.y, z2),
        ),
        4 => Aabb3::new(
            Point3::new(parent.min.x, parent.min.y, z1),
            Point3::new(x2, y2, parent.max.z),
        ),
        5 => Aabb3::new(
            Point3::new(parent.min.x, y1, parent.min.z),
            Point3::new(x2, parent.max.y, z2),
        ),
        6 => Aabb3::new(
            Point3::new(x1, parent.min.y, parent.min.z),
            Point3::new(parent.max.x, y2, z2),
        ),
        _ => Aabb3::new(parent.min(), Point3::new(x2, y2, z2)),
    };
    if a.contains(node) {
        return child;
    }
    if !parent.contains(node) {
        return 9;
    }
    let a = Aabb3::new(Point3::new(x1, y1, z1), parent.max());
    if a.contains(node) {
        return 0;
    }
    let a = Aabb3::new(
        Point3::new(parent.min.x, y1, z1),
        Point3::new(x2, parent.max.y, parent.max.z),
    );
    if a.contains(node) {
        return 1;
    }
    let a = Aabb3::new(
        Point3::new(x1, parent.min.y, z1),
        Point3::new(parent.max.x, y2, parent.max.z),
    );
    if a.contains(node) {
        return 2;
    }
    let a = Aabb3::new(
        Point3::new(parent.min.x, parent.min.y, z1),
        Point3::new(x2, y2, parent.max.z),
    );
    if a.contains(node) {
        return 4;
    }
    let a = Aabb3::new(
        Point3::new(parent.min.x, y1, parent.min.z),
        Point3::new(x2, parent.max.y, z2),
    );
    if a.contains(node) {
        return 5;
    }
    let a = Aabb3::new(
        Point3::new(x1, parent.min.y, parent.min.z),
        Point3::new(parent.max.x, y2, z2),
    );
    if a.contains(node) {
        return 6;
    }
    let a = Aabb3::new(parent.min(), Point3::new(x2, y2, z2));
    if a.contains(node) {
        return 7;
    }
    return 3;
}
// ab节点下降
fn down<S: BaseNum, T>(
    slab: &mut Slab<OctNode<S>>,
    adjust: usize,
    deep: usize,
    oct_id: usize,
    node: &mut AbNode<S, T>,
    id: usize,
) -> (usize, usize) {
    let parent = unsafe { slab.get_unchecked_mut(oct_id) };
    if parent.layer >= node.layer {
        node.parent = oct_id;
        node.next = parent.nodes.head;
        parent.nodes.push(id);
        return (0, 0);
    }
    let two = S::one() + S::one();
    #[macro_use()]
    macro_rules! child_macro {
        ($a:ident, $i:tt) => {
            if $a.contains(&node.aabb) {
                match parent.childs[$i] {
                    ChildNode::Oct(oct, ref mut num) => {
                        *num += 1;
                        return down(slab, adjust, deep, oct, node, id);
                    }
                    ChildNode::Ab(ref mut list) => {
                        node.parent = oct_id;
                        node.parent_child = $i;
                        node.next = list.head;
                        list.push(id);
                        if list.len > adjust && parent.layer < deep {
                            return set_dirty(&mut parent.dirty, $i, parent.layer, oct_id);
                        }
                        return (0, 0);
                    }
                }
            }
        };
    }
    let x1 = (parent.aabb.min.x + parent.aabb.max.x - parent.loose.x) / two;
    let y1 = (parent.aabb.min.y + parent.aabb.max.y - parent.loose.y) / two;
    let z1 = (parent.aabb.min.z + parent.aabb.max.z - parent.loose.z) / two;
    let x2 = (parent.aabb.min.x + parent.aabb.max.x + parent.loose.x) / two;
    let y2 = (parent.aabb.min.y + parent.aabb.max.y + parent.loose.y) / two;
    let z2 = (parent.aabb.min.z + parent.aabb.max.z + parent.loose.z) / two;
    let a = Aabb3::new(Point3::new(x1, y1, z1), parent.aabb.max());
    child_macro!(a, 0);
    let a = Aabb3::new(
        Point3::new(parent.aabb.min.x, y1, z1),
        Point3::new(x2, parent.aabb.max.y, parent.aabb.max.z),
    );
    child_macro!(a, 1);
    let a = Aabb3::new(
        Point3::new(x1, parent.aabb.min.y, z1),
        Point3::new(parent.aabb.max.x, y2, parent.aabb.max.z),
    );
    child_macro!(a, 2);
    let a = Aabb3::new(
        Point3::new(x1, y1, parent.aabb.min.z),
        Point3::new(parent.aabb.max.x, parent.aabb.max.y, z2),
    );
    child_macro!(a, 3);
    let a = Aabb3::new(
        Point3::new(parent.aabb.min.x, parent.aabb.min.y, z1),
        Point3::new(x2, y2, parent.aabb.max.z),
    );
    child_macro!(a, 4);
    let a = Aabb3::new(
        Point3::new(parent.aabb.min.x, y1, parent.aabb.min.z),
        Point3::new(x2, parent.aabb.max.y, z2),
    );
    child_macro!(a, 5);
    let a = Aabb3::new(
        Point3::new(x1, parent.aabb.min.y, parent.aabb.min.z),
        Point3::new(parent.aabb.max.x, y2, z2),
    );
    child_macro!(a, 6);
    let a = Aabb3::new(parent.aabb.min(), Point3::new(x2, y2, z2));
    child_macro!(a, 7);
    (0, 0)
}
// 更新aabb
fn update<S: BaseNum, T>(
    slab: &mut Slab<OctNode<S>>,
    adjust: &(usize, usize),
    deep: usize,
    outer: &mut NodeList,
    dirty: &mut (Vec<Vec<usize>>, usize, usize),
    id: usize,
    node: &mut AbNode<S, T>,
) -> Option<(usize, usize, usize, usize, usize)> {
    let old_p = node.parent;
    if old_p > 0 {
        let old_c = node.parent_child;
        let mut parent = unsafe { slab.get_unchecked_mut(old_p) };
        if node.layer > parent.layer {
            // ab节点能在当前Oct节点的容纳范围
            // 获得新位置
            let child = check_contain(&parent.aabb, &parent.loose, &node.aabb, old_c);
            if old_c == child {
                return None;
            }
            if child < 8 {
                let prev = node.prev;
                let next = node.next;
                node.prev = 0;
                // 移动到兄弟节点
                match parent.childs[child] {
                    ChildNode::Oct(oct, ref mut num) => {
                        *num += 1;
                        node.parent_child = 8;
                        set_tree_dirty(dirty, down(slab, adjust.1, deep, oct, node, id));
                        return Some((old_p, old_c, prev, next, node.next));
                    }
                    ChildNode::Ab(ref mut list) => {
                        node.parent_child = child;
                        node.next = list.head;
                        list.push(id);
                        if list.len > adjust.1 && node.layer < deep {
                            set_dirty(&mut parent.dirty, child, parent.layer, id);
                        }
                        return Some((old_p, old_c, prev, next, node.next));
                    }
                }
            }
        // 需要向上
        } else if node.layer == parent.layer {
            if parent.aabb.contains(&node.aabb) {
                if old_c == 8 {
                    return None;
                }
                let prev = node.prev;
                let next = node.next;
                node.prev = 0;
                // 从child 移到 nodes
                node.parent_child = 8;
                node.next = parent.nodes.head;
                parent.nodes.push(id);
                return Some((old_p, old_c, prev, next, node.next));
            }
        // 在当前节点外
        } else {
            // 比当前节点大
        };
        let prev = node.prev;
        let next = node.next;
        if old_p > 1 {
            // 向上移动
            let mut p = parent.parent;
            let mut c = parent.parent_child;
            loop {
                parent = unsafe { slab.get_unchecked_mut(p) };
                match parent.childs[c] {
                    ChildNode::Oct(_, ref mut num) => {
                        *num -= 1;
                        if *num < adjust.0 {
                            let d = set_dirty(&mut parent.dirty, c, parent.layer, p);
                            if d.1 > 0 {
                                set_tree_dirty(dirty, d);
                            }
                        }
                    }
                    _ => panic!("invalid state"),
                }
                if parent.layer <= node.layer && parent.aabb.contains(&node.aabb) {
                    node.prev = 0;
                    node.parent_child = 8;
                    set_tree_dirty(dirty, down(slab, adjust.1, deep, p, node, id));
                    return Some((old_p, old_c, prev, next, node.next));
                }
                p = parent.parent;
                c = parent.parent_child;
                if p == 0 {
                    break;
                }
            }
        }
        // 判断根节点是否相交
        if intersects(&parent.aabb, &node.aabb) {
            if old_p == 1 && old_c == 8 {
                return None;
            }
            // 相交的放在root的nodes上
            node.parent = 1;
            node.next = parent.nodes.head;
            parent.nodes.push(id);
        } else {
            node.parent = 0;
            node.next = outer.head;
            outer.push(id);
        }
        node.prev = 0;
        node.parent_child = 8;
        return Some((old_p, old_c, prev, next, node.next));
    } else {
        // 边界外物体更新
        let root = unsafe { slab.get_unchecked_mut(1) };
        if intersects(&root.aabb, &node.aabb) {
            // 判断是否相交或包含
            let prev = node.prev;
            let next = node.next;
            node.prev = 0;
            node.parent_child = 8;
            if root.aabb.contains(&node.aabb) {
                set_tree_dirty(dirty, down(slab, adjust.1, deep, 1, node, id));
            } else {
                // 相交的放在root的nodes上
                node.parent = 1;
                node.next = root.nodes.head;
                root.nodes.push(id);
            }
            Some((0, 0, prev, next, node.next))
        } else {
            // 表示还在outer上
            None
        }
    }
}
// 从NodeList中移除，并可能添加
pub fn remove_add<S: BaseNum, T>(
    tree: &mut Tree<S, T>,
    id: usize,
    r: Option<(usize, usize, usize, usize, usize)>,
) {
    // 从NodeList中移除
    if let Some((rid, child, prev, next, cur_next)) = r {
        if rid > 0 {
            let oct = unsafe { tree.oct_slab.get_unchecked_mut(rid) };
            if child < 8 {
                match oct.childs[child] {
                    ChildNode::Ab(ref mut ab) => ab.remove(&mut tree.ab_map, prev, next),
                    _ => panic!("invalid state"),
                }
            } else {
                oct.nodes.remove(&mut tree.ab_map, prev, next);
            }
        } else {
            tree.outer.remove(&mut tree.ab_map, prev, next);
        }
        if cur_next > 0 {
            let n = unsafe { tree.ab_map.get_unchecked_mut(cur_next) };
            n.prev = id;
        }
    }
}

// 移除时，向上修改数量，并可能设脏
#[inline]
fn remove_up<S: BaseNum>(
    slab: &mut Slab<OctNode<S>>,
    adjust: usize,
    dirty: &mut (Vec<Vec<usize>>, usize, usize),
    parent: usize,
    child: usize,
) {
    if parent == 0 {
        return;
    }
    let (p, c) = {
        let node = unsafe { slab.get_unchecked_mut(parent) };
        match node.childs[child] {
            ChildNode::Oct(_, ref mut num) => {
                *num -= 1;
                if *num < adjust {
                    let d = set_dirty(&mut node.dirty, child, node.layer, parent);
                    if d.1 > 0 {
                        set_tree_dirty(dirty, d);
                    }
                }
            }
            _ => panic!("invalid state"),
        }
        (node.parent, node.parent_child)
    };
    remove_up(slab, adjust, dirty, p, c);
}

#[inline]
fn set_dirty(dirty: &mut usize, index: usize, layer: usize, rid: usize) -> (usize, usize) {
    if *dirty == 0 {
        *dirty |= 1 << index;
        return (layer, rid);
    }
    *dirty |= 1 << index;
    return (0, 0);
}
// 设置脏标记
#[inline]
fn set_tree_dirty(dirty: &mut (Vec<Vec<usize>>, usize, usize), (layer, rid): (usize, usize)) {
    if rid == 0 {
        return;
    }
    dirty.1 += 1;
    if dirty.2 > layer {
        dirty.2 = layer;
    }
    if dirty.0.len() <= layer {
        for _ in dirty.0.len()..layer + 1 {
            dirty.0.push(Vec::new())
        }
    }
    let vec = unsafe { dirty.0.get_unchecked_mut(layer) };
    vec.push(rid);
}

// 创建指定的子节点
fn create_child<S: BaseNum>(
    aabb: &Aabb3<S>,
    loose: &Vector3<S>,
    layer: usize,
    parent_id: usize,
    child: usize,
) -> OctNode<S> {
    let two = S::one() + S::one();
    let x1 = (aabb.min.x + aabb.max.x - loose.x) / two;
    let y1 = (aabb.min.y + aabb.max.y - loose.y) / two;
    let z1 = (aabb.min.z + aabb.max.z - loose.z) / two;
    let x2 = (aabb.min.x + aabb.max.x + loose.x) / two;
    let y2 = (aabb.min.y + aabb.max.y + loose.y) / two;
    let z2 = (aabb.min.z + aabb.max.z + loose.z) / two;
    let a = match child {
        0 => Aabb3::new(Point3::new(x1, y1, z1), aabb.max()),
        1 => Aabb3::new(
            Point3::new(aabb.min.x, y1, z1),
            Point3::new(x2, aabb.max.y, aabb.max.z),
        ),
        2 => Aabb3::new(
            Point3::new(x1, aabb.min.y, z1),
            Point3::new(aabb.max.x, y2, aabb.max.z),
        ),
        3 => Aabb3::new(
            Point3::new(x1, y1, aabb.min.z),
            Point3::new(aabb.max.x, aabb.max.y, z2),
        ),
        4 => Aabb3::new(
            Point3::new(aabb.min.x, aabb.min.y, z1),
            Point3::new(x2, y2, aabb.max.z),
        ),
        5 => Aabb3::new(
            Point3::new(aabb.min.x, y1, aabb.min.z),
            Point3::new(x2, aabb.max.y, z2),
        ),
        6 => Aabb3::new(
            Point3::new(x1, aabb.min.y, aabb.min.z),
            Point3::new(aabb.max.x, y2, z2),
        ),
        _ => Aabb3::new(aabb.min(), Point3::new(x2, y2, z2)),
    };
    return OctNode::new(a, loose / two, parent_id, child, layer + 1);
}

// 整理方法，只有整理方法才会创建或销毁OctNode
fn collect<S: BaseNum, T>(
    oct_slab: &mut Slab<OctNode<S>>,
    ab_map: &mut VecMap<AbNode<S, T>>,
    adjust: &(usize, usize),
    deep: usize,
    parent_id: usize,
) {
    let (dirty, childs, ab, loose, layer) = {
        let parent = match oct_slab.get_mut(parent_id) {
            Some(oct) => oct,
            _ => return,
        };
        let dirty = parent.dirty;
        if parent.dirty == 0 {
            return;
        }
        parent.dirty = 0;
        (
            dirty,
            parent.childs.clone(),
            parent.aabb.clone(),
            parent.loose.clone(),
            parent.layer,
        )
    };
    for i in 0..8 {
        if dirty & (1 << i) != 0 {
            match childs[i] {
                ChildNode::Oct(oct, num) if num < adjust.0 => {
                    let mut list = NodeList::new();
                    if num > 0 {
                        shrink(oct_slab, ab_map, parent_id, i, oct, &mut list);
                    }
                    let parent = unsafe { oct_slab.get_unchecked_mut(parent_id) };
                    parent.childs[i] = ChildNode::Ab(list);
                }
                ChildNode::Ab(ref list) if list.len > adjust.1 => {
                    let child_id = split(
                        oct_slab, ab_map, adjust, deep, list, &ab, &loose, layer, parent_id, i,
                    );
                    let parent = unsafe { oct_slab.get_unchecked_mut(parent_id) };
                    parent.childs[i] = ChildNode::Oct(child_id, list.len);
                }
                _ => (),
            }
        }
    }
}
// 收缩OctNode
fn shrink<S: BaseNum, T>(
    oct_slab: &mut Slab<OctNode<S>>,
    ab_map: &mut VecMap<AbNode<S, T>>,
    parent: usize,
    parent_child: usize,
    oct_id: usize,
    result: &mut NodeList,
) {
    let node = oct_slab.remove(oct_id);
    if node.nodes.len > 0 {
        shrink_merge(ab_map, parent, parent_child, &node.nodes, result);
    }
    #[macro_use()]
    macro_rules! child_macro {
        ($i:tt) => {
            match node.childs[$i] {
                ChildNode::Ab(ref list) if list.len > 0 => {
                    shrink_merge(ab_map, parent, parent_child, &list, result);
                }
                ChildNode::Oct(oct, len) if len > 0 => {
                    shrink(oct_slab, ab_map, parent, parent_child, oct, result);
                }
                _ => (),
            }
        };
    }
    child_macro!(0);
    child_macro!(1);
    child_macro!(2);
    child_macro!(3);
    child_macro!(4);
    child_macro!(5);
    child_macro!(6);
    child_macro!(7);
}
// 合并ab列表到结果列表中
#[inline]
fn shrink_merge<S: BaseNum, T>(
    ab_map: &mut VecMap<AbNode<S, T>>,
    parent: usize,
    parent_child: usize,
    list: &NodeList,
    result: &mut NodeList,
) {
    let old = result.head;
    result.head = list.head;
    result.len += list.len;
    let mut id = list.head;
    loop {
        let ab = unsafe { ab_map.get_unchecked_mut(id) };
        ab.parent = parent;
        ab.parent_child = parent_child;
        if ab.next == 0 {
            ab.next = old;
            break;
        }
        id = ab.next;
    }
    if old > 0 {
        let ab = unsafe { ab_map.get_unchecked_mut(old) };
        ab.prev = id;
    }
}

// 分裂出OctNode
#[inline]
fn split<S: BaseNum, T>(
    oct_slab: &mut Slab<OctNode<S>>,
    ab_map: &mut VecMap<AbNode<S, T>>,
    adjust: &(usize, usize),
    deep: usize,
    list: &NodeList,
    parent_ab: &Aabb3<S>,
    parent_loose: &Vector3<S>,
    parent_layer: usize,
    parent_id: usize,
    child: usize,
) -> usize {
    let oct = create_child(parent_ab, parent_loose, parent_layer, parent_id, child);
    let oct_id = oct_slab.insert(oct);
    let oct = unsafe { oct_slab.get_unchecked_mut(oct_id) };
    if split_down(ab_map, adjust.1, deep, oct, oct_id, list) > 0 {
        collect(oct_slab, ab_map, adjust, deep, oct_id);
    }
    oct_id
}
// 将ab节点列表放到分裂出来的八叉节点上
fn split_down<S: BaseNum, T>(
    map: &mut VecMap<AbNode<S, T>>,
    adjust: usize,
    deep: usize,
    parent: &mut OctNode<S>,
    parent_id: usize,
    list: &NodeList,
) -> usize {
    let two = S::one() + S::one();
    let x1 = (parent.aabb.min.x + parent.aabb.max.x - parent.loose.x) / two;
    let y1 = (parent.aabb.min.y + parent.aabb.max.y - parent.loose.y) / two;
    let z1 = (parent.aabb.min.z + parent.aabb.max.z - parent.loose.z) / two;
    let x2 = (parent.aabb.min.x + parent.aabb.max.x + parent.loose.x) / two;
    let y2 = (parent.aabb.min.y + parent.aabb.max.y + parent.loose.y) / two;
    let z2 = (parent.aabb.min.z + parent.aabb.max.z + parent.loose.z) / two;
    #[macro_use()]
    macro_rules! child_macro {
        ($a:ident, $node:ident, $id:tt, $i:tt) => {
            if $a.contains(&$node.aabb) {
                match parent.childs[$i] {
                    ChildNode::Ab(ref mut list) => {
                        $node.parent = parent_id;
                        $node.parent_child = $i;
                        $node.next = list.head;
                        list.push($id);
                        if list.len > adjust && parent.layer < deep {
                            set_dirty(&mut parent.dirty, $i, parent.layer, parent_id);
                        }
                        continue;
                    }
                    _ => panic!("invalid state"),
                }
            }
        };
    }
    let mut id = list.head;
    while id > 0 {
        let node = unsafe { map.get_unchecked_mut(id) };
        let nid = id;
        id = node.next;
        node.prev = 0;
        if parent.layer >= node.layer {
            node.parent = parent_id;
            node.parent_child = 8;
            node.next = parent.nodes.head;
            parent.nodes.push(nid);
            continue;
        }
        id = node.next;
        let a = Aabb3::new(Point3::new(x1, y1, z1), parent.aabb.max());
        child_macro!(a, node, nid, 0);
        let a = Aabb3::new(
            Point3::new(parent.aabb.min.x, y1, z1),
            Point3::new(x2, parent.aabb.max.y, parent.aabb.max.z),
        );
        child_macro!(a, node, nid, 1);
        let a = Aabb3::new(
            Point3::new(x1, parent.aabb.min.y, z1),
            Point3::new(parent.aabb.max.x, y2, parent.aabb.max.z),
        );
        child_macro!(a, node, nid, 2);
        let a = Aabb3::new(
            Point3::new(x1, y1, parent.aabb.min.z),
            Point3::new(parent.aabb.max.x, parent.aabb.max.y, z2),
        );
        child_macro!(a, node, nid, 3);
        let a = Aabb3::new(
            Point3::new(parent.aabb.min.x, parent.aabb.min.y, z1),
            Point3::new(x2, y2, parent.aabb.max.z),
        );
        child_macro!(a, node, nid, 4);
        let a = Aabb3::new(
            Point3::new(parent.aabb.min.x, y1, parent.aabb.min.z),
            Point3::new(x2, parent.aabb.max.y, z2),
        );
        child_macro!(a, node, nid, 5);
        let a = Aabb3::new(
            Point3::new(x1, parent.aabb.min.y, parent.aabb.min.z),
            Point3::new(parent.aabb.max.x, y2, z2),
        );
        child_macro!(a, node, nid, 6);
        let a = Aabb3::new(parent.aabb.min(), Point3::new(x2, y2, z2));
        child_macro!(a, node, nid, 7);
    }
    fix_prev(map, parent.nodes.head);
    for i in 0..8 {
        match parent.childs[i] {
            ChildNode::Ab(ref list) => fix_prev(map, list.head),
            _ => (), // panic
        }
    }
    parent.dirty
}
// 修复prev
#[inline]
fn fix_prev<S: BaseNum, T>(map: &mut VecMap<AbNode<S, T>>, mut head: usize) {
    if head == 0 {
        return;
    }
    let node = unsafe { map.get_unchecked(head) };
    let mut next = node.next;
    while next > 0 {
        let node = unsafe { map.get_unchecked_mut(next) };
        node.prev = head;
        head = next;
        next = node.next;
    }
}

// 查询空间内及相交的ab节点
fn query<S: BaseNum, T, A, B>(
    oct_slab: &Slab<OctNode<S>>,
    ab_map: &VecMap<AbNode<S, T>>,
    oct_id: usize,
    oct_arg: &A,
    oct_func: fn(arg: &A, aabb: &Aabb3<S>) -> bool,
    ab_arg: &mut B,
    ab_func: fn(arg: &mut B, id: usize, aabb: &Aabb3<S>, bind: &T),
) {
    let node = unsafe { oct_slab.get_unchecked(oct_id) };
    let mut id = node.nodes.head;
    while id > 0 {
        let ab = unsafe { ab_map.get_unchecked(id) };
        ab_func(ab_arg, id, &ab.aabb, &ab.bind);
        id = ab.next;
    }
    #[macro_use()]
    macro_rules! child_macro {
        ($a:ident, $i:tt) => {
            match node.childs[$i] {
                ChildNode::Oct(oct, ref num) if *num > 0 => {
                    if oct_func(oct_arg, &$a) {
                        query(oct_slab, ab_map, oct, oct_arg, oct_func, ab_arg, ab_func);
                    }
                }
                ChildNode::Ab(ref list) if list.head > 0 => {
                    if oct_func(oct_arg, &$a) {
                        let mut id = list.head;
                        loop {
                            let ab = unsafe { ab_map.get_unchecked(id) };
                            ab_func(ab_arg, id, &ab.aabb, &ab.bind);
                            id = ab.next;
                            if id == 0 {
                                break;
                            }
                        }
                    }
                }
                _ => (),
            }
        };
    }
    let two = S::one() + S::one();
    let x1 = (node.aabb.min.x + node.aabb.max.x - node.loose.x) / two;
    let y1 = (node.aabb.min.y + node.aabb.max.y - node.loose.y) / two;
    let z1 = (node.aabb.min.z + node.aabb.max.z - node.loose.z) / two;
    let x2 = (node.aabb.min.x + node.aabb.max.x + node.loose.x) / two;
    let y2 = (node.aabb.min.y + node.aabb.max.y + node.loose.y) / two;
    let z2 = (node.aabb.min.z + node.aabb.max.z + node.loose.z) / two;
    let a = Aabb3::new(Point3::new(x1, y1, z1), node.aabb.max());
    child_macro!(a, 0);
    let a = Aabb3::new(
        Point3::new(node.aabb.min.x, y1, z1),
        Point3::new(x2, node.aabb.max.y, node.aabb.max.z),
    );
    child_macro!(a, 1);
    let a = Aabb3::new(
        Point3::new(x1, node.aabb.min.y, z1),
        Point3::new(node.aabb.max.x, y2, node.aabb.max.z),
    );
    child_macro!(a, 2);
    let a = Aabb3::new(
        Point3::new(x1, y1, node.aabb.min.z),
        Point3::new(node.aabb.max.x, node.aabb.max.y, z2),
    );
    child_macro!(a, 3);
    let a = Aabb3::new(
        Point3::new(node.aabb.min.x, node.aabb.min.y, z1),
        Point3::new(x2, y2, node.aabb.max.z),
    );
    child_macro!(a, 4);
    let a = Aabb3::new(
        Point3::new(node.aabb.min.x, y1, node.aabb.min.z),
        Point3::new(x2, node.aabb.max.y, z2),
    );
    child_macro!(a, 5);
    let a = Aabb3::new(
        Point3::new(x1, node.aabb.min.y, node.aabb.min.z),
        Point3::new(node.aabb.max.x, y2, z2),
    );
    child_macro!(a, 6);
    let a = Aabb3::new(node.aabb.min(), Point3::new(x2, y2, z2));
    child_macro!(a, 7);
}

// 和指定的列表进行碰撞
fn collision_list<S: BaseNum, T, A>(
    map: &VecMap<AbNode<S, T>>,
    id: usize,
    aabb: &Aabb3<S>,
    bind: &T,
    arg: &mut A,
    func: fn(
        arg: &mut A,
        a_id: usize,
        a_aabb: &Aabb3<S>,
        a_bind: &T,
        b_id: usize,
        b_aabb: &Aabb3<S>,
        b_bind: &T,
    ) -> bool,
    mut head: usize,
) {
    while head > 0 {
        let b = unsafe { map.get_unchecked(head) };
        func(arg, id, aabb, bind, head, &b.aabb, &b.bind);
        head = b.next;
    }
}

// 和指定的节点进行碰撞
// fn collision_node<S: BaseNum, T, A>(
//   oct_slab: &Slab<OctNode<S>>,
//   ab_map: &Slab<AbNode<S, T>>,
//   id: usize,
//   aabb: &Aabb3<S>,
//   bind: &T,
//   arg: &mut A,
//   func: fn(arg: &mut A, a_id: usize, a_aabb: &Aabb3<S>, a_bind: &T, b_id: usize, b_aabb: &Aabb3<S>, b_bind: &T) -> bool,
//   parent: usize,
//   parent_child: usize,
// ) {

// }

#[test]
fn test1() {
    println!("test1-----------------------------------------");

    let mut tree = Tree::new(
        Aabb3::new(
            Point3::new(-1024f32, -1024f32, -4194304f32),
            Point3::new(3072f32, 3072f32, 4194304f32),
        ),
        0,
        0,
        0,
        0,
    );
    for i in 0..1 {
        tree.add(
            i + 1,
            Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
            i + 1,
        );
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("00000, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    for i in 1..tree.ab_map.len() + 1 {
        println!("00000, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.collect();
    for i in 1..tree.ab_map.len() + 1 {
        println!("00000, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    for i in 1..5 {
        tree.add(
            i + 1,
            Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
            i + 3,
        );
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("00001, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );

    tree.update(
        4,
        Aabb3::new(
            Point3::new(0.0, 700.0, 0.0),
            Point3::new(1000.0, 1400.0, 1.0),
        ),
    );

    tree.update(
        5,
        Aabb3::new(
            Point3::new(0.0, 1400.0, 0.0),
            Point3::new(1000.0, 1470.0, 1.0),
        ),
    );
    tree.update(
        6,
        Aabb3::new(
            Point3::new(0.0, 1470.0, 0.0),
            Point3::new(1000.0, 1540.0, 1.0),
        ),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.collect();
    for i in 1..tree.ab_map.len() + 1 {
        println!("00002, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    //   tree.update(1, Aabb3::new(Point3::new(0.0,0.0,0.0), Point3::new(1000.0, 800.0, 1.0)));
    //   tree.update(2, Aabb3::new(Point3::new(0.0,0.0,0.0), Point3::new(1000.0, 800.0, 1.0)));
    //   tree.update(3, Aabb3::new(Point3::new(0.0,0.0,0.0), Point3::new(1000.0, 800.0, 1.0)));
    //   tree.update(4, Aabb3::new(Point3::new(0.0,0.0,0.0), Point3::new(1000.0, 800.0, 1.0)));

    //   tree.update(5, Aabb3::new(Point3::new(0.0,800.0,0.0), Point3::new(1000.0, 1600.0, 1.0)));

    //    tree.update(6, Aabb3::new(Point3::new(0.0,1600.0,0.0), Point3::new(1000.0, 2400.0, 1.0)));
    //   tree.update(7, Aabb3::new(Point3::new(0.0,2400.0,0.0), Point3::new(1000.0, 3200.0, 1.0)));
    //   for i in 1..tree.ab_map.len() + 1 {
    //   println!("22222, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    //  }
    // tree.collect();
    for i in 1..tree.oct_slab.len() + 1 {
        println!(
            "000000 000000, id:{}, oct: {:?}",
            i,
            tree.oct_slab.get(i).unwrap()
        );
    }
    println!("outer:{:?}", tree.outer);
    let aabb = Aabb3::new(
        Point3::new(500f32, 500f32, -4194304f32),
        Point3::new(500f32, 500f32, 4194304f32),
    );
    let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());
    tree.query(&aabb, intersects, &mut args, ab_query_func);
    //assert_eq!(args.result(), [1, 3, 4]);
}

#[test]
fn test2() {
    println!("test2-----------------------------------------");
    let mut tree = Tree::new(
        Aabb3::new(
            Point3::new(0f32, 0f32, 0f32),
            Point3::new(1000f32, 1000f32, 1000f32),
        ),
        0,
        0,
        0,
        0,
    );
    for i in 0..9 {
        tree.add(
            i + 1,
            Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
            i + 1,
        );
    }
    println!("loose:{:?} deep:{}", tree.loose, tree.deep);
    for i in 1..tree.oct_slab.len() + 1 {
        println!("000000, id:{}, oct: {:?}", i, tree.oct_slab.get(i).unwrap());
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("00000, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0)),
    );
    tree.collect();
    for i in 1..tree.oct_slab.len() + 1 {
        println!(
            "000000 000000, id:{}, oct: {:?}",
            i,
            tree.oct_slab.get(i).unwrap()
        );
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!(
            "000000 000000, id:{}, ab: {:?}",
            i,
            tree.ab_map.get(i).unwrap()
        );
    }
    println!("tree -new ------------------------------------------");
    let mut tree = Tree::new(
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0)),
        0,
        0,
        0,
        0,
    );
    for i in 0..6 {
        tree.add(
            i + 1,
            Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.1, 0.1, 0.1)),
            i + 1,
        );
    }
    for i in 1..tree.oct_slab.len() + 1 {
        println!("test1, id:{}, oct: {:?}", i, tree.oct_slab.get(i).unwrap());
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("test1, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.collect();
    for i in 1..tree.oct_slab.len() + 1 {
        println!("test2, id:{}, oct: {:?}", i, tree.oct_slab.get(i).unwrap());
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("test2, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.shift(4, Vector3::new(2.0, 2.0, 1.0));
    tree.shift(5, Vector3::new(4.0, 4.0, 1.0));
    tree.shift(6, Vector3::new(10.0, 10.0, 1.0));
    for i in 1..tree.oct_slab.len() + 1 {
        println!("test3, id:{}, oct: {:?}", i, tree.oct_slab.get(i).unwrap());
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("test3, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    tree.collect();
    for i in 1..tree.oct_slab.len() + 1 {
        println!("test4, id:{}, oct: {:?}", i, tree.oct_slab.get(i).unwrap());
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("test4, id:{}, ab: {:?}", i, tree.ab_map.get(i).unwrap());
    }
    println!("outer:{:?}", tree.outer);
    let aabb = Aabb3::new(
        Point3::new(0.05f32, 0.05f32, 0f32),
        Point3::new(0.05f32, 0.05f32, 1000f32),
    );
    let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());
    tree.query(&aabb, intersects, &mut args, ab_query_func);
    //assert_eq!(args.result(), [1, 2, 3]);
}

#[test]
fn test3() {
    let z_max: f32 = 4194304.0;
    let aabb = Aabb3::new(
        Point3::new(700.0, 100.0, -z_max),
        Point3::new(700.0, 100.0, z_max),
    );
    let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());

    // let mut tree = Tree::new(Aabb3::new(Point3::new(0f32,0f32,0f32), Point3::new(1000f32,1000f32,1000f32)),
    // 	0,
    // 	0,
    // 	0,
    // 	0,
    // );
    let mut tree = Tree::new(
        Aabb3::new(
            Point3::new(-1024f32, -1024f32, -z_max),
            Point3::new(3072f32, 3072f32, z_max),
        ),
        0,
        0,
        0,
        0,
    );
    tree.add(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();

    tree.update(
        0,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );

    tree.add(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        8,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        9,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 350.0, 1.0)),
    );
    tree.update(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(800.0, 175.0, 1.0)),
    );
    tree.update(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(640.0, 140.0, 1.0)),
    );
    tree.update(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(512.0, 112.0, 1.0)),
    );
    tree.update(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(410.0, 90.0, 1.0)),
    );
    tree.update(
        8,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1600.0, 175.0, 1.0),
        ),
    );
    tree.update(
        9,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1440.0, 140.0, 1.0),
        ),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.query(&aabb, intersects, &mut args, ab_query_func);
    tree.remove(7);
    tree.remove(6);
    tree.remove(5);

    tree.add(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();
    tree.update(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(640.0, 140.0, 1.0)),
    );
    tree.update(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(512.0, 112.0, 1.0)),
    );
    tree.update(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(410.0, 90.0, 1.0)),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 350.0, 1.0)),
    );
    tree.update(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(800.0, 175.0, 1.0)),
    );
    tree.update(
        8,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1600.0, 175.0, 1.0),
        ),
    );
    tree.update(
        9,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1440.0, 140.0, 1.0),
        ),
    );
    for i in 1..tree.oct_slab.len() + 1 {
        println!("test||||||, id:{}, oct: {:?}", i, tree.oct_slab.get(i));
    }
    for i in 1..tree.ab_map.len() + 1 {
        println!("test----------, id:{}, ab: {:?}", i, tree.ab_map.get(i));
    }
    println!(
        "-------------------------------------------------------dirtys:{:?}",
        tree.dirty
    );
    tree.query(&aabb, intersects, &mut args, ab_query_func);

    tree.remove(7);
    tree.remove(6);
    tree.remove(5);

    tree.add(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();

    tree.update(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(640.0, 140.0, 1.0)),
    );
    tree.update(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(512.0, 112.0, 1.0)),
    );
    tree.update(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(410.0, 90.0, 1.0)),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 350.0, 1.0)),
    );
    tree.update(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(800.0, 175.0, 1.0)),
    );
    tree.update(
        8,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1600.0, 175.0, 1.0),
        ),
    );
    tree.update(
        9,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1440.0, 140.0, 1.0),
        ),
    );
    tree.query(&aabb, intersects, &mut args, ab_query_func);
    for i in 1..tree.oct_slab.len() + 1 {
        println!("test||||||, id:{}, oct: {:?}", i, tree.oct_slab.get(i));
    }
    for i in 1..tree.ab_map.len() + 10 {
        println!("test----------, id:{}, ab: {:?}", i, tree.ab_map.get(i));
    }
    println!(
        "-------------------------------------------------------outer:{:?}",
        tree.outer
    );
}

#[cfg(test)]
extern crate rand;
#[test]
fn test4() {
    use rand;
    use rand::Rng;
    let z_max: f32 = 4194304.0;
    let aabb = Aabb3::new(
        Point3::new(700.0, 100.0, -z_max),
        Point3::new(700.0, 100.0, z_max),
    );
    let mut args: AbQueryArgs<f32, usize> = AbQueryArgs::new(aabb.clone());

    let mut tree = Tree::new(
        Aabb3::new(
            Point3::new(-1024f32, -1024f32, -z_max),
            Point3::new(3072f32, 3072f32, z_max),
        ),
        0,
        0,
        0,
        0,
    );
    tree.add(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();

    tree.update(
        0,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );

    tree.add(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        8,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        9,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 350.0, 1.0)),
    );
    tree.update(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(800.0, 175.0, 1.0)),
    );
    tree.update(
        5,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(640.0, 140.0, 1.0)),
    );
    tree.update(
        6,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(512.0, 112.0, 1.0)),
    );
    tree.update(
        7,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(410.0, 90.0, 1.0)),
    );
    tree.update(
        8,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1600.0, 175.0, 1.0),
        ),
    );
    tree.update(
        9,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1440.0, 140.0, 1.0),
        ),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );

    tree.add(
        10,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        11,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        12,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();
    tree.update(
        10,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(640.0, 140.0, 1.0)),
    );
    tree.update(
        11,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(512.0, 112.0, 1.0)),
    );
    tree.update(
        12,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(410.0, 90.0, 1.0)),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 350.0, 1.0)),
    );
    tree.update(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(800.0, 175.0, 1.0)),
    );
    tree.update(
        8,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1600.0, 175.0, 1.0),
        ),
    );
    tree.update(
        9,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1440.0, 140.0, 1.0),
        ),
    );

    tree.add(
        13,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        14,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.add(
        15,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)),
        1,
    );
    tree.collect();

    log(&tree.oct_slab, &tree.ab_map, 10000);
    tree.update(
        13,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(640.0, 140.0, 1.0)),
    );
    log(&tree.oct_slab, &tree.ab_map, 13);
    println!("oct========================");
    tree.update(
        14,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(512.0, 112.0, 1.0)),
    );
    tree.update(
        15,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(410.0, 90.0, 1.0)),
    );
    tree.update(
        1,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        2,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 700.0, 1.0)),
    );
    tree.update(
        3,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1000.0, 350.0, 1.0)),
    );
    tree.update(
        4,
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(800.0, 175.0, 1.0)),
    );
    tree.update(
        8,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1600.0, 175.0, 1.0),
        ),
    );
    tree.update(
        9,
        Aabb3::new(
            Point3::new(800.0, 0.0, 0.0),
            Point3::new(1440.0, 140.0, 1.0),
        ),
    );
    let mut rng = rand::thread_rng();
    for _ in 0..1000000 {
        let r = rand::thread_rng().gen_range(1, 16);
        let x = rand::thread_rng().gen_range(0, 128) - 64;
        let y = rand::thread_rng().gen_range(0, 128) - 64;
        let old = clone_tree(&tree.oct_slab, &tree.ab_map);
        tree.shift(r, Vector3::new(x as f32, y as f32, 0.0));
        assert_eq!(check_tree(&tree.oct_slab, &tree.ab_map, old, r), false);
        if x + y > 120 {
            let old = clone_tree(&tree.oct_slab, &tree.ab_map);
            tree.collect();
            assert_eq!(check_tree(&tree.oct_slab, &tree.ab_map, old, r), false);
        }
    }
}
#[cfg(test)]
fn clone_tree(
    oct_slab: &Slab<OctNode<f32>>,
    ab_map: &VecMap<AbNode<f32, usize>>,
) -> (Slab<OctNode<f32>>, VecMap<AbNode<f32, usize>>) {
    (oct_slab.clone(), ab_map.clone())
}

#[cfg(test)]
fn check_tree(
    oct_slab: &Slab<OctNode<f32>>,
    ab_map: &VecMap<AbNode<f32, usize>>,
    old: (Slab<OctNode<f32>>, VecMap<AbNode<f32, usize>>),
    r: usize,
) -> bool {
    for (id, _n) in oct_slab.iter() {
        if check_ab(oct_slab, ab_map, id) {
            log(&old.0, &old.1, r);
            log(oct_slab, ab_map, id);
            return true;
        }
    }
    false
}
#[cfg(test)]
fn check_ab(
    oct_slab: &Slab<OctNode<f32>>,
    ab_map: &VecMap<AbNode<f32, usize>>,
    oct_id: usize,
) -> bool {
    let node = unsafe { oct_slab.get_unchecked(oct_id) };
    let mut old = VecMap::default();
    for c in 0..8 {
        match &node.childs[c] {
            ChildNode::Ab(list) => {
                if check_list(ab_map, oct_id, c, &list, &mut old) {
                    return true;
                };
            }
            _ => (),
        }
    }
    check_list(ab_map, oct_id, 8, &node.nodes, &mut old)
}
#[cfg(test)]
fn check_list(
    ab_map: &VecMap<AbNode<f32, usize>>,
    parent: usize,
    parent_child: usize,
    list: &NodeList,
    old: &mut VecMap<usize>,
) -> bool {
    let mut id = list.head;
    let mut prev = 0;
    let mut i = 0;
    while id > 0 {
        old.insert(id, id);
        let ab = unsafe { ab_map.get_unchecked(id) };
        if ab.prev != prev {
            println!("------------0-oct_id: {}, ab_id: {}", parent, id);
            return true;
        }
        if ab.parent != parent {
            println!("------------1-oct_id: {}, ab_id: {}", parent, id);
            return true;
        }
        if ab.parent_child != parent_child {
            println!("------------2-oct_id: {}, ab_id: {}", parent, id);
            return true;
        }
        if old.contains(ab.next) {
            println!("------------3-oct_id: {}, ab_id: {}", parent, id);
            return true;
        }
        prev = id;
        id = ab.next;
        i += 1;
    }
    if i != list.len {
        println!("------------4-oct_id: {}, ab_id: {}", parent, id);
        return true;
    }
    return false;
}

#[cfg(test)]
fn log(oct_slab: &Slab<OctNode<f32>>, ab_map: &VecMap<AbNode<f32, usize>>, oct_id: usize) {
    println!("oct_id----------, id:{}", oct_id);
    let mut i = 0;
    for or in ab_map.iter() {
        i += 1;
        if let Some(r) = or {
            //let r = ab_map.get(r).unwrap();
            println!("ab----------, id:{}, ab: {:?}", i, r);
        }
    }
    for (id, n) in oct_slab.iter() {
        //let r = ab_map.get(r).unwrap();
        println!("oct=========, id:{}, oct: {:?}", id, n);
    }
}
