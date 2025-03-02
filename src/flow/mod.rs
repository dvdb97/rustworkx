use std::u128;

use crate::digraph;
use crate::edge_weights_from_callable;
use crate::weight_callable;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::EdgeIndex;
use rustworkx_core::flow;

use hashbrown::HashMap;

use pyo3::prelude::*;
use pyo3::Python;
use pyo3::exceptions::PyIndexError;



#[pyfunction]
#[pyo3(
    signature=(graph, source, sink, cap_fn),
    text_signature = "(graph, source, sink, cap_fn)"
)]
pub fn ford_fulkerson(
    py: Python,
    graph: &digraph::PyDiGraph,
    source: usize,
    sink: usize,
    cap_fn: Option<PyObject>
) -> PyResult<HashMap<(usize, usize), u64>> {
    if !graph.graph.contains_node(NodeIndex::new(source)) {
        return Err(PyIndexError::new_err(format!(
            "Node source index \"{source}\" out of graph bound"
        )));
    }

    if !graph.graph.contains_node(NodeIndex::new(sink)) {
        return Err(PyIndexError::new_err(format!(
            "Node sink index \"{sink}\" out of graph bound"
        )));
    }

    let source = NodeIndex::new(source);
    let sink = NodeIndex::new(sink);

    let edge_caps: Vec<Option<u64>> = edge_weights_from_callable(py, &graph.graph, &cap_fn, u64::MAX)?;
    let int_cap_fn = |e: EdgeIndex| -> Result<u64, PyErr> {
        match edge_caps[e.index()] {
            Some(weight) => Ok(weight),
            None => Err(PyIndexError::new_err("No edge found for index")),
        }
    };

    let edges = graph.edges();
    
    flow::ford_fulkerson(&graph.graph, source, sink, |e| weight_callable(py, &cap_fn, e.weight(), 0u64))
}