use std::{cmp::Ordering, hash::Hash, io::Sink};
use hashbrown::HashMap;

use petgraph::{adj::NodeIndex, csr::IndexType, data::Build, prelude::StableDiGraph, visit::{EdgeCount, EdgeRef, IntoEdges, IntoNodeIdentifiers, NodeCount, Visitable}};

use crate::graph_ext::EdgeRemovable;
use crate::flow::ford_fulkerson;


fn compute_balanced_flow<G, B, F, E>(
    graph: G,
    mut bal_fn: B,
    mut cap_fn: F,
) -> Result<HashMap<(usize, usize), u64>, E>
where
    G: IntoEdges + Visitable + EdgeRemovable + NodeCount + EdgeCount + IntoNodeIdentifiers,
    G::NodeId: Eq + IndexType + Hash,
    B: FnMut(G::NodeId)  -> Result<i64, E>,
    F: FnMut(G::EdgeRef) -> Result<u64, E>,
{
    let mut balances = vec![0; graph.node_count()];
    let mut all_zeros = true;

    for node in graph.node_identifiers() {
        balances[node.index()] = bal_fn(node)?; 
        all_zeros = all_zeros && (balances[node.index()] == 0);
    }

    if all_zeros {
        Ok(HashMap::new())
    } else {
        let mut ext_graph: StableDiGraph<(), ()> = StableDiGraph::with_capacity(
            graph.node_count()+2, 
            graph.edge_count() + graph.node_count()
        );

        let mut node_mapping = vec![NodeIndex::default(); graph.node_count()];

        for node in graph.node_identifiers() {
            node_mapping[node.index()] = ext_graph.add_node(());
        }

        for edge in graph.edge_references() {
            ext_graph.add_edge(node_mapping[edge.source().index()], node_mapping[edge.target().index()], ());
        }

        let source = ext_graph.add_node(());
        let sink = ext_graph.add_node(());

        for idx in 0..balances.len() {
            if balances[idx] > 0 {
                ext_graph.add_edge(source, node_mapping[idx], ());
            }

            if balances[idx] < 0 {
                ext_graph.add_edge(node_mapping[idx], sink, ());
            }
        }

        ford_fulkerson(&ext_graph, source, sink, |edge| {
            Ok(0u64)
        })
    }
}


pub fn cycle_canceling<G, B, F, C, E>(
    graph: G,
    mut bal_fn: B,
    mut cap_fn: F,
    mut cost_fn: C
) -> Result<HashMap<(usize, usize), u64>, E>
where
    G: IntoEdges + Visitable + EdgeRemovable + NodeCount + EdgeCount + IntoNodeIdentifiers,
    G::NodeId: Eq + IndexType + Hash,
    B: FnMut(G::NodeId)  -> Result<i64, E>,
    F: FnMut(G::EdgeRef) -> Result<u64, E>,
    C: FnMut(G::EdgeRef) -> Result<i64, E>
{
    let mut flow = compute_balanced_flow(graph, bal_fn, cap_fn)?;

    Ok(flow)
}