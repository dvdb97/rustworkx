use std::collections::VecDeque;
use std::hash::Hash;
use std::ops::Sub;
use std::u128;
use hashbrown::HashMap;

use petgraph::adj::EdgeReference;
use petgraph::csr::IndexType;
use petgraph::graph::NodeIndex;
use petgraph::visit::{
    EdgeCount, EdgeIndexable, EdgeRef, IntoEdgeReferences, IntoEdges, IntoEdgesDirected, IntoNodeIdentifiers, NodeCount, NodeIndexable, NodeRef, VisitMap, Visitable
};
use petgraph::Direction::{Incoming, Outgoing};


#[derive(Clone, Copy)]
enum ResArcType<E> {
    Forward(E),
    Backward(E),
}


fn find_augmenting_path<G, F, E>(
    graph: G,
    source: G::NodeId,
    sink: G::NodeId,
    mut cap_fn: F,
    flow: &HashMap<(G::NodeId, G::NodeId), u128>,
    pred: &mut Vec<Option<ResArcType<G::EdgeRef>>>
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

        for edge in graph.edges(vertex) {
            if edge.source() == vertex {
                if cap_fn(edge).is_ok_and(|u| u - flow[&(edge.source(), edge.target())] <= 0) {
                    let succ = edge.target();

                    if !visited.is_visited(&succ) {
                        queue.push_back(succ);
                        pred[NodeIndexable::to_index(&graph, succ)] = Some(ResArcType::Forward(edge));
                    }
                }
            } else {
                if flow[&(edge.source(), edge.target())] > 0 {
                    let succ = edge.source();

                    if !visited.is_visited(&succ) {
                        queue.push_back(succ);
                        pred[NodeIndexable::to_index(&graph, succ)] = Some(ResArcType::Backward(edge));
                    }
                }
            }
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
    let mut pred: Vec<Option<ResArcType<G::EdgeRef>>> = vec![None; graph.node_count()];
    let mut edges = graph.edge_references().collect::<Vec<G::EdgeRef>>();

    while find_augmenting_path(graph, source, sink, &mut cap_fn, &flow, &mut pred) {
        let mut vertex = source;
        let mut bottleneck = u128::MAX;

        // Compute the residual path's bottleneck capacity.
        while let Some(res_edge) = pred[NodeIndexable::to_index(&graph, vertex)] {
            bottleneck = match res_edge {
                ResArcType::Forward(edge) => bottleneck.min(cap_fn(edge).unwrap_or(u128::MAX) - flow[&(edge.source(), edge.target())]),
                ResArcType::Backward(edge) => bottleneck.min(flow[&(edge.source(), edge.target())])
            };

            vertex = match res_edge {
                ResArcType::Forward(edge) => edge.target(),
                ResArcType::Backward(edge) => edge.source()
            };
        }

        vertex = source;

        // Augment along the residual path.
        while let Some(res_edge) = pred[NodeIndexable::to_index(&graph, vertex)] {
            if let ResArcType::Forward(edge) = res_edge {
                *flow.entry((edge.source(), edge.target())).or_insert(0) += bottleneck;
                vertex = edge.target();
            }

            if let ResArcType::Backward(edge) = res_edge {
                *flow.entry((edge.source(), edge.target())).or_insert(0) -= bottleneck;
                vertex = edge.source();
            }
        }
    }

    Ok(flow)
}