use std::collections::VecDeque;
use std::hash::Hash;
use std::ops::Sub;
use hashbrown::HashMap;

use petgraph::adj::EdgeReference;
use petgraph::csr::IndexType;
use petgraph::graph::NodeIndex;
use petgraph::visit::{
    EdgeCount, EdgeIndexable, EdgeRef, IntoEdgeReferences, IntoEdges, IntoEdgesDirected, IntoNodeIdentifiers, NodeCount, NodeIndexable, NodeRef, VisitMap, Visitable
};
use petgraph::Direction::{Incoming, Outgoing};


enum ResArcType {
    Forward(u128),
    Backward(u128),
}


fn residual_capacity<G, F, E>(
    graph: G,
    v: G::NodeId, 
    w: G::NodeId,
    cap_fn: F,
    flow: &HashMap<(G::NodeId, G::NodeId), u128>,
    edges: Vec<G::EdgeRef>
) -> ResArcType
where
    G: IntoEdges + Visitable + NodeIndexable + NodeCount + IntoNodeIdentifiers + IntoEdgesDirected + EdgeIndexable,
    G::NodeId: Eq + IndexType + Hash,
    F: FnMut(G::EdgeRef) -> Result<u128, E>
{
    match (flow.contains_key(&(v, w)), flow.contains_key(&(w, v))) {
        (true, true)    => {
            let fwd = cap_fn((v, w)) - flow[&(v, w)];
            let bwd = flow[&(v, w)];

            if fwd >= bwd {
                ResArcType::Forward(fwd)
            } else {
                ResArcType::Backward(bwd)
            }
        }
        (true, false)   => ResArcType::Backward(),
        (false, true)   => 0,
        (false, false)  => 0
    }
}

fn find_augmenting_path<G, F, E>(
    graph: G,
    source: G::NodeId,
    sink: G::NodeId,
    mut cap_fn: F,
    flow: &HashMap<(G::NodeId, G::NodeId), u128>,
    pred: &mut Vec<Option<G::EdgeRef>>
) -> bool
where 
    G: IntoEdges + Visitable + NodeIndexable + NodeCount + IntoNodeIdentifiers + IntoEdgesDirected + EdgeIndexable,
    G::NodeId: Eq + IndexType + Hash,
    F: FnMut(G::EdgeRef) -> Result<u128, E>
{
    let mut queue = VecDeque::from([source]);
    let mut visited = graph.visit_map();

    while let Some(vertex) = queue.pop_front() {
        if vertex == sink {
            return true;
        }

        visited.visit(vertex);

        for edge in graph.edges_directed(vertex, Outgoing) {
            let succ = edge.target();

            if cap_fn(edge).is_ok_and(|u| u - flow[&(edge.source(), edge.target())] <= 0) {
                continue;
            }

            if !visited.is_visited(&succ) {
                queue.push_back(succ);
                pred[NodeIndexable::to_index(&graph, succ)] = Some(edge);
            }
        }

        for edge in graph.edges_directed(vertex, Incoming) {
            
        }
    }

    false
}

pub fn ford_fulkerson<G, F, E>(
    graph: G,
    source: G::NodeId,
    sink: G::NodeId,
    mut cap_fn: F,
) -> Result<HashMap<(G::NodeId, G::NodeId), u128>, E>
where
    G: IntoEdges + Visitable + NodeIndexable + NodeCount + EdgeCount + IntoNodeIdentifiers + IntoEdgesDirected + EdgeIndexable,
    G::NodeId: Eq + IndexType + Hash,
    F: FnMut(G::EdgeRef) -> Result<u128, E>
{
    let mut flow = HashMap::new();
    let mut pred = vec![None; graph.node_count()];
    let mut edges = graph.edge_references().collect::<Vec<G::EdgeRef>>();

    while find_augmenting_path(graph, source, sink, &mut cap_fn, &flow, &pred) {
        let mut vertex = source;
        let mut vertex_index = NodeIndexable::to_index(&graph, vertex);

        let mut bottleneck = 0;

        while let Some(edge) = pred[vertex_index] {

        }
    }

    Ok(flow)
}