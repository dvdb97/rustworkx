use std::hash::Hash;
use hashbrown::HashMap;

use petgraph::csr::IndexType;
use petgraph::visit::{
    Control, ControlFlow, EdgeCount, EdgeRef, IntoEdges, NodeCount, NodeIndexable, Visitable
};

use crate::Result;
use crate::traversal::{breadth_first_search, BfsEvent};
use crate::flow::residual::ResArcType;


pub fn find_augmenting_path<G>(
    graph: G,
    source: G::NodeId,
    sink: G::NodeId,
    caps: &HashMap<(usize, usize), u64>,
    flow: &HashMap<(usize, usize), u64>,
    preds: &mut Vec<Option<ResArcType<(usize, usize)>>>
) -> bool
where 
    G: IntoEdges + Visitable + NodeIndexable,
    G::NodeId: Eq + IndexType + Hash
{
    let rslt = breadth_first_search(&graph, Some(source), |event| {
        if let BfsEvent::TreeEdge(u, v, _) = event {
            let edge = (u.index(), v.index());
            let bwd_edge = (v.index(), u.index());

            if flow.contains_key(&edge) && caps[&edge] - flow[&edge] > 0 {
                preds[v.index()] = Some(ResArcType::Forward(edge));
            } else if flow.contains_key(&bwd_edge) && flow[&bwd_edge] > 0 {
                preds[v.index()] = Some(ResArcType::Backward(bwd_edge));
            }

            if v == sink {
                return Control::Break(v);
            }
        }
        Control::Continue
    });

    rslt.should_break()
}


pub fn ford_fulkerson<G, F, E>(
    graph: G,
    source: G::NodeId,
    sink: G::NodeId,
    mut cap_fn: F,
) -> Result<HashMap<(usize, usize), u64>, E>
where
    G: IntoEdges + Visitable + NodeIndexable + EdgeCount + NodeCount,
    G::NodeId: Eq + IndexType + Hash,
    F: FnMut(G::EdgeRef) -> Result<u64, E>
{
    let mut flow: HashMap<(usize, usize), u64> = HashMap::with_capacity(graph.edge_count());
    let mut caps: HashMap<(usize, usize), u64> = HashMap::with_capacity(graph.edge_count());

    for edge_ref in graph.edge_references() {
        let edge = (edge_ref.source().index(), edge_ref.target().index());
        flow.entry(edge).insert(0u64);
        caps.entry(edge).insert(cap_fn(edge_ref)?);
    }

    // Predecessor mapping for shortest augmenting paths.
    let mut pred: Vec<Option<ResArcType<(usize, usize)>>> = vec![None; graph.node_count()];

    while find_augmenting_path(&graph, source, sink, &caps, &flow, &mut pred) {
        let mut vertex = NodeIndexable::to_index(&graph, sink);
        let mut bottleneck = u64::MAX;

        // Compute the residual path's bottleneck capacity.
        while let Some(res_edge) = pred[vertex] {
            bottleneck = match res_edge {
                ResArcType::Forward(edge) => bottleneck.min(caps[&edge] - flow[&edge]),
                ResArcType::Backward(edge) => bottleneck.min(flow[&(edge.1, edge.0)])
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
    use crate::flow::{find_augmenting_path, ford_fulkerson, residual::ResArcType};
    use crate::petgraph::graph::DiGraph;
    use crate::petgraph::visit::{IntoEdges};

    use hashbrown::HashMap;
    use petgraph::prelude::*;

    #[test]
    fn test_find_augmenting_path() {
        let mut g: DiGraph<&str, ()> = DiGraph::new();
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

        let flow: HashMap<(usize, usize), u64> = HashMap::from([
            ((a.index(), b.index()), 2), 
            ((b.index(), c.index()), 2), 
            ((c.index(), d.index()), 2),
            ((d.index(), e.index()), 2)
        ]);
        
        let caps: HashMap<(usize, usize), u64> = HashMap::from_iter(g.edge_references().map(|e| ((e.source().index(), e.target().index()), 2u64)));
        let mut preds: Vec<Option<ResArcType<(usize, usize)>>> = vec![None; g.node_count()];

        let rslt = find_augmenting_path(&g, a, e, &caps, &flow, &mut preds); 

        assert!(rslt);       
    }
}