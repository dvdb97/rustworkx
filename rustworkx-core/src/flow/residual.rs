use std::{hash::Hash};

use ahash::HashMap;
use petgraph::{adj::NodeIndex, csr::IndexType, data::Build, prelude::{DiGraphMap, GraphMap, StableDiGraph}, visit::{EdgeCount, EdgeRef, IntoEdgeReferences, IntoEdges, IntoEdgesDirected, IntoNodeIdentifiers, NodeCount, NodeIndexable, Visitable}, Graph};


struct ResNet<G>
where
    G: IntoEdgesDirected + Visitable + NodeIndexable + NodeCount,
    G::NodeId: Eq + IndexType + Hash
{
    graph: G,
    flow: HashMap<(usize, usize), i64>,

}