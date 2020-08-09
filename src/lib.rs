mod utils;
use utils::*;

use rand::thread_rng;
use rand::seq::SliceRandom;
use itertools::Itertools;
use kdbush::KDBush;
use std::collections::HashMap;
use counter::Counter;
use rayon::prelude::*;

// pyo3 dependencies
use pyo3::prelude::*;
use pyo3::exceptions::ValueError as PyValueError;
use pyo3::wrap_pyfunction;

#[pymodule]
fn neighborhood_analysis(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<CellCombs>()?;
    m.add_wrapped(wrap_pyfunction!(get_neighbors))?;
    Ok(())
}

/// A utility function to search for neighbors
///
/// Args:
///     points: List[tuple(float, float)]; Two dimension points
///     r: float; The search radius
///
/// Return:
///     A dictionary of the index of every points, with the index of its neighbors
#[pyfunction]
fn get_neighbors(points: Vec<(f64, f64)>, r: f64)
    -> HashMap<usize, Vec<usize>>{
    let tree = KDBush::create(points.to_owned(), kdbush::DEFAULT_NODE_SIZE); // make an index
    let mut result: HashMap<usize, Vec<usize>> = (0..points.len()).map(|i| (i, vec![])).collect();
    for (i, p) in points.iter().enumerate() {
        tree.within(p.0, p.1, r, |id| result.get_mut(&i).unwrap().push(id));
    }

    result
}

/// Constructor function
///
/// Args:
///     types: List[str]; All the type of cells in your research
///     order: bool (False); If False, the ('A', 'B') and ('A', 'B') is the same.
///
/// Return:
///     self
///
#[pyclass]
struct CellCombs {
    #[pyo3(get)]
    cell_types: PyObject,
    #[pyo3(get)]
    cell_combs: PyObject,
    #[pyo3(get)]
    cell_relationships: PyObject,
}

unsafe impl Send for CellCombs {}

#[pymethods]
impl CellCombs {
    #[new]
    fn new(py:Python, types: PyObject, order: Option<bool>) -> PyResult<Self> {

        let types_data: Vec<&str> = match types.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyValueError::py_err("Can't resolve `types`, should be list of string."))
            }
        };

        let order_data: bool = match order {
            Some(data) => data,
            None => false,
        };

        let uni: Vec<&str> = types_data.into_iter().unique().collect();
        let mut combs = vec![];
        let mut relationships = HashMap::new();

        for i1 in uni.to_owned() {
            relationships.insert(i1, vec![]);
            for i2 in uni.to_owned() {
                if !&order_data {
                    combs.push(vec![i1, i2]);
                }
                relationships.get_mut(i1).unwrap().push(vec![i1, i2]);
            }
        }

        if order_data {
            let end: usize = uni.len();
            for (i, e) in uni.to_owned().iter().enumerate() {
                for t in i..end {
                    combs.push(vec![e, uni[t]]);
                }
            }
        }

        let uni_py = uni.to_object(py);
        let combs_py = combs.to_object(py);
        let relationships_py = relationships.to_object(py);

        Ok(CellCombs {
            cell_types: uni_py,
            cell_combs: combs_py,
            cell_relationships: relationships_py,
        })
    }

    /// Bootstrap functions
    ///
    /// If method is 'pval', 1.0 means association, -1.0 means avoidance.
    /// If method is 'zscore', results is the exact z-score value.
    ///
    /// Args:
    ///     types: List[str]; The type of all the cells
    ///     neighbors: Dict[int, List[int]]; eg. {1:[4,5], 2:[6,7]}, cell at index 1 has neighbor cells from index 4 and 5
    ///     times: int (500); How many times to perform bootstrap
    ///     pval: float (0.05); The threshold of p-value
    ///     method: str ('pval'); 'pval' or 'zscore'
    ///
    /// Return:
    ///     List of tuples, eg.(['a', 'b'], 1.0), the type a and type b has a relationship as association
    ///
    fn bootstrap(&self,
                 py: Python,
                 types: PyObject,
                 neighbors: PyObject,
                 times: Option<usize>,
                 pval: Option<f64>,
                 method: Option<&str>,
    ) -> PyResult<PyObject> {

        let types_data: Vec<&str> = match types.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyValueError::py_err("Can't resolve `types`, should be list of string."))
            }
        };
        let neighbors_data: HashMap<usize, Vec<usize>> = match neighbors.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyValueError::py_err("Can't resolve `neighbors`, should be a dict."))
            }
        };

        let times = match times {
            Some(data) => data,
            None => 500
        };

        let pval = match pval {
            Some(data) => data,
            None => 0.05
        };

        let method = match method {
            Some(data) => data,
            None => "pval"
        };

        let cellcombs: Vec<Vec<&str>> = match self.cell_combs.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyValueError::py_err("Resolve cell_combs failed."))
            }
        };
        let cellrelatetionship: HashMap<&str, Vec<Vec<&str>>> = match self.cell_relationships.extract(py) {
            Ok(data) => data,
            Err(_) => {
                return Err(PyValueError::py_err("Resolve cell_relationship failed."))
            }
        };

        let real_data = count_neighbors(
            &types_data,
            &neighbors_data,
            &cellcombs,
            &cellrelatetionship
        );

        let mut simulate_data = cellcombs.iter()
            .map(|comb| (comb.to_owned(), vec![]))
            .collect::<HashMap<Vec<&str>, Vec<f64>>>();

        let all_data: Vec<HashMap<Vec<&str>, f64>> = (0..times).into_par_iter().map(|_| {
            let mut rng = thread_rng();
            let mut shuffle_types = types_data.to_owned();
            shuffle_types.shuffle(&mut rng);
            let perm_result = count_neighbors(
                &shuffle_types,
                &neighbors_data,
            &cellcombs,
            &cellrelatetionship);
            perm_result
            })
            .collect();

        for perm_result in all_data {
            for (k, v) in perm_result.iter() {
                simulate_data.get_mut(k).unwrap().push(*v);
        }
        };
/*
        let mut results = cellcombs.iter()
            .map(|comb| (comb.to_owned(), 0.0))
            .collect::<HashMap<Vec<&str>, f64>>();*/
        let mut results: Vec<(Vec<&str>, f64)> = vec![];

        for (k, v) in simulate_data.iter() {
            let real = real_data[k];

            if method == "pval" {
                let mut gt: i32 = 0;
                let mut lt: i32 = 0;
                for i in v.iter() {
                    if i >= &real {
                        gt += 1
                    } else if i <= &real {
                        lt += 1
                    }
                }
                let gt: f64 = gt as f64 / (times.to_owned() as i32 + 1) as f64;
                let lt: f64 = lt as f64 / (times.to_owned() as i32 + 1) as f64;
                let dir: f64 = (gt < lt) as i32 as f64;
                let udir: f64 = -dir;
                let p: f64 = gt * dir + lt * udir;
                let sig: f64 = (p < pval) as i32 as f64;
                let sigv: f64 = sig * (dir - 0.5).signum();
                //println!("{:?} {:?} {:?}", dir, p, sig);
                //println!("{:?} {:?}", k, sigv);
                // *results.get_mut(k).unwrap() += sigv;
                results.push((k.to_owned(), sigv));
            } else {
                let m = mean_f(v);
                let sd = std(v);
                if sd != 0.0 {
                    results.push((k.to_owned(), (real - m) / sd));
                    //*results.get_mut(k).unwrap() += (real - m) / sd;
                } else {
                    results.push((k.to_owned(), 0.0));
                    //*results.get_mut(k).unwrap() += 0.0;
                }

            }

        }

        let results_py = results.to_object(py);

        Ok(results_py)
    }
}

fn count_neighbors<'a>(
                   types: &Vec<&str>,
                   neighbors: &HashMap<usize, Vec<usize>>,
                   cell_combs: &Vec<Vec<&'a str>>,
                   cell_relationships: &HashMap<&'a str, Vec<Vec<&'a str>>>,)
               -> HashMap<Vec<&'a str>, f64> {
    let mut storage = cell_combs.iter()
        .map(|comb| (comb.to_owned(), vec![]))
        .collect::<HashMap<Vec<&str>, Vec<usize>>>();

    for (k, v) in neighbors.iter() {
        let cent_type = types[*k];
        let neigh_type: Counter<_> = v.iter().map(|i| types[*i]).collect::<Counter<_>>();
        for comb in cell_relationships[cent_type].iter() {
            let counts = neigh_type.get(comb[1]).unwrap_or(&0);
            match storage.get_mut(comb) {
                None => {
                    storage.get_mut(&vec![comb[1], comb[0]]).unwrap().push(*counts)
                },
                Some(s) => s.push(*counts),
            };
        }
    }

    let mut results: HashMap<Vec<&'a str>, f64>= HashMap::new();
    for (k, v) in storage.iter() {
        results.insert(k.to_owned(), mean(&v));
    }

    results
}
