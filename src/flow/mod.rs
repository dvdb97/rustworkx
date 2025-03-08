use crate::digraph;
use crate::weight_callable;
use petgraph::graph::NodeIndex;
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
    flow::ford_fulkerson(&graph.graph, source, sink, |e| weight_callable(py, &cap_fn, e.weight(), 0u64))
}


#[pyfunction]
#[pyo3(
    signature=(graph, cap_fn, cost_fn),
    text_signature = "(graph, cap_fn, cost_fn)"
)]
pub fn cycle_canceling(
    py: Python,
    graph: &digraph::PyDiGraph,
    cap_fn: Option<PyObject>,
    cost_fn: Option<PyObject>
) -> PyResult<HashMap<(usize, usize), u64>> {
    let flow = HashMap::new();   
    flow::cycle_canceling(&graph.graph, flow, |e| weight_callable(py, &cap_fn, e.weight(), u64::MAX), |e| weight_callable(py, &cost_fn, e.weight(), 0i64))
}