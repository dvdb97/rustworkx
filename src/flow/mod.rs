use crate::digraph;
use crate::weight_callable;
use crate::CostFn;
use petgraph::graph::NodeIndex;
use petgraph::graph::NodeWeightsMut;
use petgraph::visit::NodeRef;
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
    signature=(graph, flow, cap_fn, cost_fn),
    text_signature = "(graph, flow, cap_fn, cost_fn)"
)]
pub fn cycle_canceling(
    py: Python,
    graph: &digraph::PyDiGraph,
    flow: HashMap<(usize, usize), u64>,
    cap_fn: Option<PyObject>,
    cost_fn: Option<PyObject>
) -> PyResult<HashMap<(usize, usize), u64>> { 
    flow::cycle_canceling(&graph.graph, flow, |e| weight_callable(py, &cap_fn, e.weight(), u64::MAX), |e| weight_callable(py, &cost_fn, e.weight(), 0i64))
}


/// Find a minimum cost flow satisfying all the demands in a directed graph
///
/// This is a primal network simplex algorithm that uses the leaving arc rule
/// to prevent cycling
///
/// ``graph`` is treated as a digraph (the trait bound don't force it but
/// for undirected inputs it will be treated as one) where edges have cost and
/// capacity and noes have demand, in other words nodes want to send or receive
/// some flow. A negative demands means that the node wants to send flow,
/// a positive demand means the node wants to receive flow. A flow on the
/// digraph ``graph`` satisfies all demand if the net flow into each node is
/// equal
///
/// Note this function uses signed integers for edge weights, capacity, and
/// node demands. This is to avoid floating point errors preventing causing
/// the algorithm to fail as it is particularly sensitive to this. If you need
/// to use a floating point number it is recommend to multiply the float by
/// a constant factor and cast it to an int (such as
/// ``demand=lambda x: int(x * 100)``) to approximate the floating point value.
///
/// This implementation is based on the NetworkX ``network_simplex()``
/// function [1], which is based on [2] and [3].
///
/// [1] https://github.com/networkx/networkx/blob/6971d845e25aa4b76879831bb85ac2fa23ae0c9e/networkx/algorithms/flow/networksimplex.py#L330
///
/// [2] Z. Kiraly, P. Kovacs.
///     Efficient implementation of minimum-cost flow algorithms.
///     Acta Universitatis Sapientiae, Informatica 4(1):67--118. 2012.
/// [3] R. Barr, F. Glover, D. Klingman.
///     Enhancement of spanning tree labeling procedures for network
///     optimization.
///     INFOR 17(1):16--34. 1979.
///
/// :param PyDiGraph graph: The input graph object to run the algorithm on
/// :param demand: A callable which will receive the node data payload
///     from ``graph`` and is expected to return the demand to indicate how
///     much flow that node wants to send (a negative value) or receive
///     (a positive value). Note that the sum of all demands in the graph
///     should be equal to 0 or the problem is not feasible.
/// :param capacity: A callable which will receive the edge data payload
///     from ``graph`` and is expected to the capacity, or how much flow
///     the edge can support.
/// * `weight` - A function which will receive the edge data payload from
///     ``graph`` and is expected to return the weight of the edge (i.e. the
///     cost incurred by sending flow on the edge).
#[pyfunction]
#[pyo3(
    signature=(graph, demand, capacity, weight),
    text_signature = "(graph, demand, capacity, weight/)"
)]
pub fn network_simplex(
    py: Python,
    graph: &digraph::PyDiGraph,
    demand: PyObject,
    capacity: Option<PyObject>,
    weight: Option<PyObject>,
) -> PyResult<(i64, HashMap<usize, HashMap<usize, i64>>)> {
    // let flow = flow::network_simplex(
    //     &graph.graph, 
    //     |n| demand.call1(py, (graph.graph.node_weight(n),))?.extract::<i64>(py), 
    //     |e| weight_callable(py, &capacity, e.weight(), i64::MAX), 
    //     |e| weight_callable(py, &weight, e.weight(), 0))?;
    let flow = flow::network_simplex(
        &graph.graph, 
        |n| demand.call1(py, (graph.graph.node_weight(n),))?.extract::<i64>(py), 
        |e| weight_callable(py, &capacity, e.weight(), i64::MAX), 
        |e| weight_callable(py, &weight, e.weight(), 0))?;

    match flow {
        Some((cost, flow)) => Ok((cost, flow)),
        None => Ok((0, HashMap::new()))
    }
}