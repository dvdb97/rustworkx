use std::u128;

use crate::digraph;
use crate::edge_weights_from_callable;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::EdgeIndex;
use petgraph::visit::EdgeRef;
use rustworkx_core::flow;

use hashbrown::HashMap;

use pyo3::prelude::*;
use pyo3::Python;
use pyo3::exceptions::PyIndexError;

use crate::weight_callable;


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
) -> PyResult<HashMap<(usize, usize), u128>> {
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

    let edge_caps: Vec<Option<u128>> = edge_weights_from_callable(py, &graph.graph, &cap_fn, u128::MAX)?;
    let int_cap_fn = |e: EdgeIndex| -> PyResult<u128> {
        match edge_caps[e.index()] {
            Some(weight) => Ok(weight),
            None => Err(PyIndexError::new_err("No edge found for index")),
        }
    };
    
    flow::ford_fulkerson(&graph.graph, source, sink, |e| int_cap_fn(e.id()))
}