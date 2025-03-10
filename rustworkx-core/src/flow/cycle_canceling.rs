use std::{cmp::{max_by, Ordering}, hash::Hash, io::Sink, u64};
use hashbrown::HashMap;

use petgraph::{adj::NodeIndex, csr::IndexType, data::Build, prelude::{DiGraphMap, GraphMap, StableDiGraph}, visit::{EdgeCount, EdgeRef, IntoEdgeReferences, IntoEdges, IntoNodeIdentifiers, NodeCount, NodeIndexable, Visitable}, Graph};

use crate::shortest_path::negative_cycle_finder;

use super::residual::ResArcType;


fn residual_edge_costs<E>(
    arc: &(usize, usize),
    flow: &HashMap<(usize, usize), u64>,
    caps: &HashMap<(usize, usize), u64>,
    costs: &HashMap<(usize, usize), i64>
) -> Result<i64, E>
{
    let fwd_arc = *arc;
    let bwd_arc = (arc.1, arc.0);

    let fwd_costs = if caps.contains_key(&fwd_arc) && caps[&fwd_arc] - *flow.get(&fwd_arc).unwrap_or(&0) > 0 {
        costs[&fwd_arc]
    } else {
        i64::MAX
    };

    let bwd_costs = if caps.contains_key(&bwd_arc) && *flow.get(&bwd_arc).unwrap_or(&0) > 0 {
        -costs[&bwd_arc]
    } else {
        i64::MAX
    };

    Ok(fwd_costs.min(bwd_costs))
}


fn find_augmenting_cycle<G, E>(
    resnet: G,
    flow: &HashMap<(usize, usize), u64>,
    caps: &HashMap<(usize, usize), u64>,
    costs: &HashMap<(usize, usize), i64>
) -> Result<Option<(Vec<ResArcType<(usize, usize)>>, u64)>, E> 
where
    G: IntoEdges + Visitable + NodeCount + IntoNodeIdentifiers + IntoEdgeReferences + NodeIndexable,
    G::NodeId: Eq + IndexType + Hash
{
    let cycle = negative_cycle_finder(
        &resnet, 
        |arc| residual_edge_costs::<E>(&(arc.source().index(), arc.target().index()), &flow, &caps, &costs)
    )?;

    match cycle {
        Some(cycle) => {
            let mut res_cycle = Vec::new();
            let mut bottleneck = u64::MAX;

            for i in 0..cycle.len()-1 {
                let arc = (cycle[i].index(), cycle[i+1].index());
                let rev_arc = (arc.1, arc.0);

                let has_fwd_arc = caps.contains_key(&arc) && caps[&arc] - *flow.get(&arc).unwrap_or(&0) > 0;
                let has_bwd_arc = caps.contains_key(&rev_arc) && *flow.get(&rev_arc).unwrap_or(&0) > 0;

                let res_arc = match (has_fwd_arc, has_bwd_arc) {
                    (true, true)    => {
                        if costs[&arc] <= -costs[&rev_arc] {
                            ResArcType::Forward(arc)
                        } else {
                            ResArcType::Backward(arc)
                        }
                    },
                    (true, false)   => ResArcType::Forward(arc),
                    (false, true)   => ResArcType::Backward(arc),
                    (false, false)  => ResArcType::Forward(arc)
                };

                res_cycle.push(res_arc);

                bottleneck = match res_arc {
                    ResArcType::Forward(arc)    => bottleneck.min(caps[&arc]),
                    ResArcType::Backward(arc)   => bottleneck.min(flow[&(arc.1, arc.0)])
                };
            }

            Ok(Some((res_cycle, bottleneck)))
        },
        None => Ok(None)
    }
}


pub fn cycle_canceling<G, F, C, E>(
    graph: G,
    flow: HashMap<(usize, usize), u64>,
    mut cap_fn: F,
    mut cost_fn: C
) -> Result<HashMap<(usize, usize), u64>, E>
where
    G: IntoEdges + Visitable + NodeCount + EdgeCount + IntoNodeIdentifiers + IntoEdgeReferences,
    G::NodeId: Eq + IndexType + Hash,
    F: FnMut(G::EdgeRef) -> Result<u64, E>,
    C: FnMut(G::EdgeRef) -> Result<i64, E>
{
    let mut flow = flow.clone();
    let mut caps: HashMap<(usize, usize), u64>  = HashMap::with_capacity(graph.edge_count());
    let mut costs: HashMap<(usize, usize), i64> = HashMap::with_capacity(graph.edge_count());

    for edge_ref in graph.edge_references() {
        let edge = (edge_ref.source().index(), edge_ref.target().index());
        *caps.entry(edge).or_default() = cap_fn(edge_ref)?;
        *costs.entry(edge).or_default() = cost_fn(edge_ref)?;
    }

    let mut resnet = Graph::new();
    let mut node_mapping = vec![NodeIndex::default(); graph.node_count()];

    for node in graph.node_identifiers() {
        node_mapping[node.index()] = resnet.add_node(());
    }

    for edge_ref in graph.edge_references() {
        let edge = (node_mapping[edge_ref.source().index()], node_mapping[edge_ref.target().index()]);
        resnet.update_edge(edge.0, edge.1, ());
        resnet.update_edge(edge.1, edge.0, ());
    }

    while let Some((cycle, bottleneck)) = find_augmenting_cycle(&resnet, &flow, &caps, &costs)?
    {
        for res_arc in cycle {
            if let ResArcType::Forward(arc) = res_arc {
                *flow.entry(arc).or_default() += bottleneck;
            }

            if let ResArcType::Backward(arc) = res_arc {
                *flow.entry((arc.1, arc.0)).or_default() -= bottleneck;
            }
        }
    }

    Ok(flow)
}


mod tests {
    use crate::flow::cycle_canceling;
    use hashbrown::HashMap;
    use petgraph::prelude::*;
    use rand::Error;

    #[test]
    fn test_single_cycle() {
        let mut graph = Graph::new();
        let a = graph.add_node("A");
        let b = graph.add_node("B");
        graph.add_edge(a, b, (4, -1));
        graph.add_edge(b, a, (3, 0));

        let flow: HashMap<(usize, usize), u64> = HashMap::from_iter([((0, 1), 0), ((1, 0), 0)]);
        let flow: Result<HashMap<(usize, usize), u64>, Error> = cycle_canceling(&graph, flow, |e| Ok(e.weight().0), |e| Ok(e.weight().1));

        assert!(flow.is_ok());
        let flow = flow.unwrap();

        assert_eq!(flow[&(0, 1)], 3);
        assert_eq!(flow[&(1, 0)], 3);
    }

    #[test]
    fn test_residual_cycle() {
        let mut graph = Graph::new();
        let a = graph.add_node("A");
        let b = graph.add_node("B");
        graph.add_edge(a, b, (10, 2));
        graph.add_edge(b, a, (4, -1));

        let flow: HashMap<(usize, usize), u64> = HashMap::from_iter([((0, 1), 3), ((1, 0), 3)]);
        let flow: Result<HashMap<(usize, usize), u64>, Error> = cycle_canceling(&graph, flow, |e| Ok(e.weight().0), |e| Ok(e.weight().1));

        assert!(flow.is_ok());
        let flow = flow.unwrap();

        assert_eq!(flow[&(0, 1)], 0);
        assert_eq!(flow[&(1, 0)], 0);     
    }

    #[test]
    fn test_empty_flow() {
        let mut graph = Graph::new();
        let a = graph.add_node("A");
        let b = graph.add_node("B");
        graph.add_edge(a, b, (5, 2));
        graph.add_edge(b, a, (10, 2));

        let flow = HashMap::new();
        let flow: Result<HashMap<(usize, usize), u64>, Error> = cycle_canceling(&graph, flow, |e| Ok(e.weight().0), |e| Ok(e.weight().1));
    
        assert!(flow.is_ok());

        let flow = flow.unwrap();
        assert_eq!(*flow.get(&(0, 1)).unwrap_or(&0), 0);
        assert_eq!(*flow.get(&(1, 0)).unwrap_or(&0), 0);
    }

    #[test]
    fn test_balance_flow() {
        let mut graph = Graph::new();
        let a = graph.add_node("A");
        let b = graph.add_node("B");
        let c = graph.add_node("C");
        let d = graph.add_node("D");
        let e = graph.add_node("E");
        let f = graph.add_node("F");

        graph.add_edge(a, b, (4, 2));
        graph.add_edge(c, b, (10, 2));
        graph.add_edge(c, f, (4, 0));
        graph.add_edge(b, e, (10, 3));
        graph.add_edge(a, d, (4, 0));
        graph.add_edge(e, d, (6, 2));
        graph.add_edge(e, f, (4, 2));

        let flow = HashMap::from_iter([((0, 1), 4), ((2, 1), 6), ((2, 5), 0), ((4, 5), 4), ((1, 4), 10), ((4, 3), 6), ((0, 3), 0)]);
        let flow: Result<HashMap<(usize, usize), u64>, Error> = cycle_canceling(&graph, flow, |e| Ok(e.weight().0), |e| Ok(e.weight().1));

        assert!(flow.is_ok());
        let flow = flow.unwrap();

        println!("{:#?}", flow);

        assert_eq!(flow[&(0, 3)], 4);
        assert_eq!(flow[&(2, 5)], 4);
        assert_eq!(flow[&(1, 4)], 2);
    }
}