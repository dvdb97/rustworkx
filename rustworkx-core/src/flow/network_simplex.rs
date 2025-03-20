use std::hash::Hash;
use hashbrown::HashMap;

use petgraph::{csr::IndexType, visit::{IntoEdges, Visitable}};


pub fn network_simplex<G, B, F, C, E>(
    graph: G,
    mut bal_fn: B,
    mut cap_fn: F,
    mut cost_fn: C
) -> Result<HashMap<(usize, usize), u64>, E>
where
G: IntoEdges + Visitable,
    G::NodeId: Eq + IndexType + Hash,
    B: FnMut(G::NodeId)  -> Result<i64, E>,
    F: FnMut(G::EdgeRef) -> Result<u64, E>,
    C: FnMut(G::EdgeRef) -> Result<i64, E>
{
    Ok(HashMap::new())
}