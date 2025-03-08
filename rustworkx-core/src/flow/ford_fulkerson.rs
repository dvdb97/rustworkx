use std::collections::VecDeque;
use std::hash::Hash;
use hashbrown::HashMap;

use petgraph::csr::IndexType;
use petgraph::visit::{
    Control, ControlFlow, EdgeCount, EdgeRef, IntoEdges, IntoEdgesDirected, NodeCount, NodeIndexable, VisitMap, Visitable
};
use petgraph::Direction;

use crate::Result;
use crate::traversal::{breadth_first_search, BfsEvent};
use crate::flow::residual::ResArcType;


fn find_augmenting_path<G>(
    graph: G,
    source: G::NodeId,
    sink: G::NodeId,
    caps: &HashMap<(usize, usize), u64>,
    flow: &HashMap<(usize, usize), u64>,
    preds: &mut Vec<Option<ResArcType<(usize, usize)>>>
) -> bool
where 
    G: IntoEdgesDirected + Visitable + NodeIndexable,
    G::NodeId: Eq + IndexType + Hash
{
    let mut queue = VecDeque::from([source]);
    let mut visited = graph.visit_map();

    while let Some(vertex) = queue.pop_front() {
        if vertex == sink {
            return true;
        }

        visited.visit(vertex);

        for edge_ref in graph.edges_directed(vertex, Direction::Outgoing) {
            let edge = (edge_ref.source().index(), edge_ref.target().index());

            if caps[&edge] - flow[&edge] > 0 {
                let succ = edge_ref.target();

                if !visited.is_visited(&succ) {
                    queue.push_back(succ);
                    preds[edge.1] = Some(ResArcType::Forward(edge));
                }
            }
        }

        for edge_ref in graph.edges_directed(vertex, Direction::Incoming) {
            let edge = (edge_ref.source().index(), edge_ref.target().index());

            if flow[&edge] > 0 {
                let succ = edge_ref.source();

                if !visited.is_visited(&succ) {
                    queue.push_back(succ);
                    preds[edge.0] = Some(ResArcType::Backward(edge));
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
) -> Result<HashMap<(usize, usize), u64>, E>
where
    G: IntoEdgesDirected + Visitable + NodeIndexable + EdgeCount + NodeCount,
    G::NodeId: Eq + IndexType + Hash,
    F: FnMut(G::EdgeRef) -> Result<u64, E>
{
    let mut flow: HashMap<(usize, usize), u64> = HashMap::with_capacity(graph.edge_count());
    let mut caps: HashMap<(usize, usize), u64> = HashMap::with_capacity(graph.edge_count());

    for edge_ref in graph.edge_references() {
        let edge = (edge_ref.source().index(), edge_ref.target().index());
        *flow.entry(edge).or_default() = 0u64;
        *caps.entry(edge).or_default() = cap_fn(edge_ref)?;
    }

    // Predecessor mapping for shortest augmenting paths.
    let mut pred: Vec<Option<ResArcType<(usize, usize)>>> = vec![None; graph.node_count()];

    while find_augmenting_path(graph, source, sink, &caps, &flow, &mut pred) {
        let mut vertex = NodeIndexable::to_index(&graph, sink);
        let mut bottleneck = u64::MAX;

        // Compute the residual path's bottleneck capacity.
        while let Some(res_edge) = pred[vertex] {
            bottleneck = match res_edge {
                ResArcType::Forward(edge) => bottleneck.min(caps[&edge] - flow[&edge]),
                ResArcType::Backward(edge) => bottleneck.min(flow[&edge])
            };

            vertex = match res_edge {
                ResArcType::Forward(edge) => edge.0,
                ResArcType::Backward(edge) => edge.1
            };
        }

        let mut vertex = NodeIndexable::to_index(&graph, sink);

        // Augment along the residual path.
        while let Some(res_edge) = pred[vertex] {
            if let ResArcType::Forward(edge) = res_edge {
                *flow.entry(edge).or_insert(0) += bottleneck;
                vertex = edge.0;
            }

            if let ResArcType::Backward(edge) = res_edge {
                *flow.entry(edge).or_insert(0) -= bottleneck;
                vertex = edge.1;
            }
        }
    }

    Ok(flow)
}


#[cfg(test)]
mod tests {
    use crate::flow::ford_fulkerson;

    use hashbrown::HashMap;
    use petgraph::prelude::*;
    use rand::Error;

    #[test]
    fn test_simple_network() {
        let mut g = Graph::new();
        let a = g.add_node("A");
        let b = g.add_node("B");
        let c = g.add_node("C");
        let d = g.add_node("D");
        let e = g.add_node("E");

        g.add_edge(a, b, ());
        g.add_edge(a, d, ());
        g.add_edge(b, c, ());
        g.add_edge(c, d, ());
        g.add_edge(d, e, ());
        g.add_edge(b, e, ());

        let flow: Result<HashMap<(usize, usize), u64>, Error> = ford_fulkerson(&g, a, e, |_| Ok(2u64));

        assert!(flow.is_ok());

        let flow = flow.unwrap();

        assert!(flow[&(0, 1)] == 2);
        assert!(flow[&(1, 4)] == 2);
        assert!(flow[&(0, 3)] == 2);
        assert!(flow[&(3, 4)] == 2);
        assert!(flow[&(1, 2)] == 0);
        assert!(flow[&(2, 3)] == 0);
    }

    #[test]
    fn test_residual_augmentations() {
        let mut graph: Graph<&str, u64> = Graph::new();
        let a = graph.add_node("A");
        let b = graph.add_node("B");
        let c = graph.add_node("C");
        let d = graph.add_node("D");
        let e = graph.add_node("E");
        let f = graph.add_node("F");
        let g = graph.add_node("G");
        let h = graph.add_node("H");
        let i = graph.add_node("I");

        graph.add_edge(a, b, 5);
        graph.add_edge(a, c, 10);
        graph.add_edge(b, e, 5);
        graph.add_edge(c, d, 5);
        graph.add_edge(c, f, 10);
        graph.add_edge(e, f, 5);
        graph.add_edge(f, i, 10);
        graph.add_edge(d, g, 5);
        graph.add_edge(g, h, 5);
        graph.add_edge(h, i, 5);

        let flow: Result<HashMap<(usize, usize), u64>, Error> = ford_fulkerson(&graph, a, i, |e| Ok(*e.weight()));
        assert!(flow.is_ok());

        let flow = flow.unwrap();
        println!("{:#?}", flow);

        assert_eq!(flow[&(2, 5)], 5);
    }

    #[test]
    fn test_non_connected_network() {
        let mut graph: Graph<&str, u64> = Graph::new();
        let a = graph.add_node("A");
        let b = graph.add_node("B");
        let c = graph.add_node("C");
        let d = graph.add_node("D");

        graph.add_edge(a, b, 10);
        graph.add_edge(c, d, 8);

        let flow: Result<HashMap<(usize, usize), u64>, Error> = ford_fulkerson(&graph, a, d, |e| Ok(*e.weight()));
        assert!(flow.is_ok());

        let flow = flow.unwrap();
        assert_eq!(flow[&(0, 1)], 0);
        assert_eq!(flow[&(2, 3)], 0);
    }
}