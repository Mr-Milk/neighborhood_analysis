mod utils;
use utils::*;

use itertools::Itertools;
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::HashMap;

use kdbush::KDBush;
use counter::Counter;
use rayon::prelude::*;
use rstar::{RTree, RTreeObject, AABB};
use spade::BoundingRect;

// pyo3 dependencies
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pymodule]
fn neighborhood_analysis(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(get_bbox))?;
    m.add_wrapped(wrap_pyfunction!(get_point_neighbors))?;
    m.add_wrapped(wrap_pyfunction!(get_bbox_neighbors))?;
    m.add_class::<CellCombs>()?;
    m.add_wrapped(wrap_pyfunction!(comb_bootstrap))?;
    m.add_wrapped(wrap_pyfunction!(neighbor_components))?;
    Ok(())
}



#[pyfunction]
pub fn neighbor_components(neighbors: HashMap<usize, Vec<usize>>, types: HashMap<usize, &str>)
    -> (Vec<usize>, Vec<&str>, Vec<Vec<usize>>) {

    let mut uni_types: HashMap<&str, i64> = HashMap::new();
    for (_, t) in &types {
        uni_types.entry(*t).or_insert(0);
    }
    let uni_types: Vec<&str> = uni_types.keys().map(|k| *k).collect_vec();
    let mut cent_order: Vec<usize> = vec![];
    let result: Vec<Vec<usize>> = neighbors.iter().map(|(cent, neigh)| {
        let count: HashMap<&&str, usize> = neigh.iter().map(|i| &types[i]).collect::<Counter<_>>().into_map();
        let mut result_v: Vec<usize> = vec![];
        for t in &uni_types {
            let v = match count.get(t) {
                Some(v) => *v,
                None => {0}
            };
            result_v.push(v);
        }
        cent_order.push(*cent);
        result_v
        }).collect();

    (cent_order, uni_types, result)

}



/// get_bbox(points_collections)
/// --
///
/// A utility function to return minimum bounding box list of polygons
///
/// Args:
///     points_collections: List[List[(float, float)]]; List of 2d points collections
///
/// Return:
///     A list of bounding box

#[pyfunction]
pub fn get_bbox(points_collections: Vec<Vec<(f64, f64)>>) -> Vec<(f64, f64, f64, f64)> {
    let bbox: Vec<(f64, f64, f64, f64)> = points_collections
        .par_iter()
        .map(|p| {
            let points: Vec<[f64; 2]> = p.iter().map(|ps| [ps.0, ps.1]).collect();
            let rect = BoundingRect::from_points(points);
            let lower: [f64; 2] = rect.lower();
            let upper: [f64; 2] = rect.upper();
            (lower[0], lower[1], upper[0], upper[1])
        })
        .collect();

    bbox
}

/// get_point_neighbors(points, r)
/// --
///
/// A utility function to search for point neighbors using kd-tree
///
/// Args:
///     points: List[tuple(float, float)]; Two dimension points
///     r: float; The search radius
///
/// Return:
///     A list of neighbors' index, return as the order of the input
///
#[pyfunction]
pub fn get_point_neighbors(points: Vec<(f64, f64)>, r: f64, labels: Option<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut has_labels = false;
    let labels: Vec<usize> = match labels {
        Some(data) => {
            has_labels = true;
            data
        },
        None => vec![0],
    };

    let tree = KDBush::create(points.to_owned(), kdbush::DEFAULT_NODE_SIZE); // make an index
    let result: HashMap<usize, Vec<usize>> = points
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let mut neighbors: Vec<usize> = vec![];
            tree.within(p.0, p.1, r, |id| neighbors.push(id));
            (i, neighbors)
        })
        .collect();

    let count = points.len();
    let mut neighbors = vec![];
    if has_labels {
        for i in 0..count {
            neighbors.push(result.get(&i).unwrap().clone().iter().map(
                |t| {labels[*t]}
            ).collect())
        }
    } else {
        for i in 0..count {
            neighbors.push(result.get(&i).unwrap().clone())
        }
    }

    neighbors
}

// customize object to insert in to R-tree
struct Rect {
    minx: f64,
    miny: f64,
    maxx: f64,
    maxy: f64,
    index: usize,
}

impl Rect {
    fn new(bbox: (f64, f64, f64, f64), index: usize) -> Rect {
        Rect {
            minx: bbox.0,
            miny: bbox.1,
            maxx: bbox.2,
            maxy: bbox.3,
            index,
        }
    }
}

impl RTreeObject for Rect {
    type Envelope = AABB<[f64; 2]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_corners([self.minx, self.miny], [self.maxx, self.maxy])
    }
}

/// get_bbox_neighbors(bbox_list, expand=1.0, scale=1.0)
/// --
///
/// A utility function to search for bbox neighbors using r-tree
///
/// Args:
///     bbox_list: List[tuple(float, float, float, float)]; The minimum bounding box of any polygon
///               (minx, miny, maxx, maxy)
///     expand: float; The expand unit
///     scale: float; The scale fold number
///
/// Return:
///     A list of neighbors' index, return as the order of the input
///
#[pyfunction]
pub fn get_bbox_neighbors(
    bbox_list: Vec<(f64, f64, f64, f64)>,
    expand: Option<f64>,
    scale: Option<f64>,
    labels: Option<Vec<usize>>
) -> Vec<Vec<usize>> {
    let mut expand_na: bool = true;
    let expand: f64 = match expand {
        Some(data) => {
            expand_na = false;
            data
        }
        None => 0.0,
    };

    let scale: f64 = match scale {
        Some(data) => data,
        None => 1.0,
    };

    let mut has_labels = false;
    let labels: Vec<usize> = match labels {
        Some(data) => {
            has_labels = true;
            data
        },
        None => vec![0],
    };

    let aabb: Vec<Rect> = bbox_list
        .par_iter()
        .enumerate()
        .map(|(i, b)| Rect::new(b.to_owned(), i))
        .collect();
    let tree: RTree<Rect> = RTree::<Rect>::bulk_load(aabb);
    let search_aabb: Vec<Rect> = {
        if expand_na == false {
            let expand_aabb: Vec<Rect> = bbox_list
                .par_iter()
                .enumerate()
                .map(|(i, b)| {
                    Rect::new((b.0 - expand, b.1 - expand, b.2 + expand, b.3 + expand), i)
                })
                .collect();
            expand_aabb
        } else {
            let scale_aabb: Vec<Rect> = bbox_list
                .par_iter()
                .enumerate()
                .map(|(i, b)| {
                    let xexpand: f64 = (b.2 - b.0) * (scale - 1.0);
                    let yexpand: f64 = (b.3 - b.1) * (scale - 1.0);
                    Rect::new(
                        (b.0 - xexpand, b.1 - yexpand, b.2 + xexpand, b.3 + yexpand),
                        i,
                    )
                })
                .collect();
            scale_aabb
        }
    };
    let result: HashMap<usize, Vec<usize>> = search_aabb
        .par_iter()
        .map(|rect| {
            let envelop = rect.envelope();
            let search_result: Vec<&Rect> =
                tree.locate_in_envelope_intersecting(&envelop).collect();
            let neighbors: Vec<usize> = search_result.iter().map(|r| r.index).collect();
            (rect.index, neighbors)
        })
        .collect();

    let count = bbox_list.len();
    let mut neighbors = vec![];
    if has_labels {
        for i in 0..count {
            neighbors.push(result.get(&i).unwrap().clone().iter().map(
                |t| {labels[*t]}
            ).collect())
        }
    } else {
        for i in 0..count {
            neighbors.push(result.get(&i).unwrap().clone())
        }
    }

    neighbors
}

/// comb_bootstrap(x_status, y_status, neighbors, times=500, ignore_self=False)
/// --
///
/// Bootstrap between two types
///
/// If you want to test co-localization between protein X and Y, first determine if the cell is X-positive
/// and/or Y-positive. True is considered as positive and will be counted.
///
/// Args:
///     x_status: List[bool]; If cell is type x
///     y_status: List[bool]; If cell is type y
///     neighbors: Dict[int, List[int]]; eg. {1:[4,5], 2:[6,7]}, cell at index 1 has neighbor cells from index 4 and 5
///     times: int (500); How many times to perform bootstrap
///     ignore_self: bool (False); Whether to consider self as a neighbor
///
/// Return:
///     The z-score for the spatial relationship between X and Y
///
#[pyfunction]
fn comb_bootstrap(
    py: Python,
    x_status: PyObject,
    y_status: PyObject,
    neighbors: PyObject,
    times: Option<usize>,
    ignore_self: Option<bool>,
) -> PyResult<f64> {
    let x: Vec<bool> = match x_status.extract(py) {
        Ok(data) => data,
        Err(_) => {
            return Err(PyTypeError::new_err(
                "Can't resolve `x_status`, should be list of bool.",
            ))
        }
    };

    let y: Vec<bool> = match y_status.extract(py) {
        Ok(data) => data,
        Err(_) => {
            return Err(PyTypeError::new_err(
                "Can't resolve `y_status`, should be list of bool.",
            ))
        }
    };

    let neighbors_data: Vec<Vec<usize>> = match neighbors.extract(py) {
        Ok(data) => data,
        Err(_) => {
            return Err(PyTypeError::new_err(
                "Can't resolve `neighbors`, should be a dict.",
            ))
        }
    };

    let times = match times {
        Some(data) => data,
        None => 500,
    };

    let ignore_self = match ignore_self {
        Some(data) => data,
        None => false,
    };
    let neighbors = utils::remove_rep_neighbors(neighbors_data, ignore_self);
    let real: f64 = comb_count_neighbors(&x, &y, &neighbors) as f64;

    let perm_counts: Vec<usize> = (0..times)
        .into_par_iter()
        .map(|_| {
            let mut rng = thread_rng();
            let mut shuffle_y = y.to_owned();
            shuffle_y.shuffle(&mut rng);
            let perm_result = comb_count_neighbors(&x, &shuffle_y, &neighbors);
            perm_result
        })
        .collect();

    let m = mean(&perm_counts);
    let sd = std(&perm_counts);

    Ok((real - m) / sd)
}

/// Constructor function
///
/// Args:
///     types: List[str]; All the type of cells in your research
///     order: bool (False); If False, A->B and A<-B is the same
///
#[pyclass]
struct CellCombs {
    #[pyo3(get)]
    cell_types: PyObject,
    #[pyo3(get)]
    cell_combs: PyObject,
    #[pyo3(get)]
    order: bool,
}

unsafe impl Send for CellCombs {}

#[pymethods]
impl CellCombs {
    #[new]
    fn new(py: Python, types: PyObject, order: Option<bool>) -> PyResult<Self> {
        let types_data: Vec<&str> = match types.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyTypeError::new_err(
                    "Can't resolve `types`, should be list of string.",
                ))
            }
        };

        let order_data: bool = match order {
            Some(data) => data,
            None => false,
        };

        let uni: Vec<&str> = types_data.into_iter().unique().collect();
        let mut combs = vec![];

        if order_data {
            for i1 in uni.to_owned() {
                for i2 in uni.to_owned() {
                    combs.push((i1, i2));
                }
            }
        } else {
            for (i1, e1) in uni.to_owned().iter().enumerate() {
                for (i2, e2) in uni.to_owned().iter().enumerate() {
                    if i2 >= i1 {
                        combs.push((e1, e2));
                    }
                }
            }
        }

        Ok(CellCombs {
            cell_types: uni.to_object(py),
            cell_combs: combs.to_object(py),
            order: order_data,
        })
    }

    /// Bootstrap functions
    ///
    /// If method is 'pval', 1.0 means association, -1.0 means avoidance, 0.0 means insignificance.
    /// If method is 'zscore', results is the exact z-score value.
    ///
    /// Args:
    ///     types: List[str]; The type of all the cells
    ///     neighbors: List[List[int]]; eg. {1:[4,5], 2:[6,7]}, cell at index 1 has neighbor cells from index 4 and 5
    ///     times: int (500); How many times to perform bootstrap
    ///     pval: float (0.05); The threshold of p-value
    ///     method: str ('pval'); 'pval' or 'zscore'
    ///     ignore_self: bool (False); Whether to consider self as a neighbor
    ///
    /// Return:
    ///     List of tuples, eg.(('a', 'b'), 1.0), the type a and type b has a relationship as association
    ///
    fn bootstrap(
        &self,
        py: Python,
        types: PyObject,
        neighbors: PyObject,
        times: Option<usize>,
        pval: Option<f64>,
        method: Option<&str>,
        ignore_self: Option<bool>,
    ) -> PyResult<PyObject> {
        let types_data: Vec<&str> = match types.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyTypeError::new_err(
                    "Can't resolve `types`, should be list of string.",
                ))
            }
        };
        let neighbors_data: Vec<Vec<usize>> = match neighbors.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyTypeError::new_err(
                    "Can't resolve `neighbors`, should be a list.",
                ))
            }
        };

        let times = match times {
            Some(data) => data,
            None => 500,
        };

        let pval = match pval {
            Some(data) => data,
            None => 0.05,
        };

        let method = match method {
            Some(data) => data,
            None => "pval",
        };

        let ignore_self = match ignore_self {
            Some(data) => data,
            None => false,
        };

        let cellcombs: Vec<(&str, &str)> = match self.cell_combs.extract(py) {
            Ok(data) => data,
            Err(_) => return Err(PyTypeError::new_err("Resolve cell_combs failed.")),
        };

        let neighbors = utils::remove_rep_neighbors(neighbors_data, ignore_self);

        let real_data = count_neighbors(&types_data, &neighbors, &cellcombs, self.order);

        let mut simulate_data = cellcombs
            .iter()
            .map(|comb| (comb.to_owned(), vec![]))
            .collect::<HashMap<(&str, &str), Vec<f64>>>();

        let all_data: Vec<HashMap<(&str, &str), f64>> = (0..times)
            .into_par_iter()
            .map(|_| {
                let mut rng = thread_rng();
                let mut shuffle_types = types_data.to_owned();
                shuffle_types.shuffle(&mut rng);
                let perm_result =
                    count_neighbors(&shuffle_types, &neighbors, &cellcombs, self.order);
                perm_result
            })
            .collect();

        for perm_result in all_data {
            for (k, v) in perm_result.iter() {
                simulate_data.get_mut(k).unwrap().push(*v);
            }
        }

        let mut results: Vec<((&str, &str), f64)> = vec![];

        for (k, v) in simulate_data.iter() {
            let real = real_data[k];

            if method == "pval" {
                let mut gt: f64 = 0.0;
                let mut lt: f64 = 0.0;
                for i in v.iter() {
                    if i >= &real {
                        gt += 1.0
                    }
                    if i <= &real {
                        lt += 1.0
                    }
                }
                let gt: f64 = gt as f64 / (times.to_owned() as f64 + 1.0);
                let lt: f64 = lt as f64 / (times.to_owned() as f64 + 1.0);
                let dir: f64 = (gt < lt) as i32 as f64;
                let udir: f64 = !(gt < lt) as i32 as f64;
                let p: f64 = gt * dir + lt * udir;
                let sig: f64 = (p < pval) as i32 as f64;
                let sigv: f64 = sig * (dir - 0.5).signum();
                results.push((k.to_owned(), sigv));
            } else {
                let m = mean_f(v);
                let sd = std_f(v);
                if sd != 0.0 {
                    results.push((k.to_owned(), (real - m) / sd));
                } else {
                    results.push((k.to_owned(), 0.0));
                }
            }
        }

        let results_py = results.to_object(py);

        Ok(results_py)
    }
}
