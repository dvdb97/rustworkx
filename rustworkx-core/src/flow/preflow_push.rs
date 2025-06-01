use std::{collections::{HashSet, VecDeque}, hash::Hash};

use hashbrown::HashMap;
use petgraph::visit::{
    EdgeCount, EdgeRef, GraphBase, GraphProp, IntoEdgeReferences, IntoEdgesDirected, IntoNodeReferences, NodeCount, NodeIndexable, NodeRef, VisitMap, Visitable
};


pub struct MaxFlowReturn<G: GraphProp> {
    pub value: u64,
    pub flow_edges: HashMap<G::NodeId, HashMap<G::NodeId, usize>>,
}


#[derive(Clone)]
pub struct Level {
    active: HashSet<usize>,
    inactive: HashSet<usize>
}

impl Level {
    pub fn new() -> Level {
        Level { active: HashSet::new(), inactive: HashSet::new() }
    }

    pub fn is_active(self, v: usize) -> bool {
        return self.active.contains(&v);
    }

    pub fn activate(&mut self, v: usize) {
        if !self.active.contains(&v) {
            self.active.insert(v);
            self.inactive.remove(&v);
        }
    }

    pub fn deactivate(&mut self, v: usize) {
        if !self.inactive.contains(&v) {
            self.inactive.insert(v);
            self.active.remove(&v);
        }
    }
}


struct PreflowPushState<G> 
where
    G: GraphProp,
    G::NodeId: Eq + Hash,
{
    s: usize,
    t: usize,
    num_nodes: usize,
    num_edges: usize,
    node_map: HashMap<G::NodeId, usize>,
    edge_map: HashMap<G::EdgeId, usize>,
    edge_capacities: Vec<u64>,
    edge_sources: Vec<Option<G::NodeId>>,
    edge_targets: Vec<Option<G::NodeId>>,
    edge_flow: Vec<u64>,
    excess: Vec<u64>,
    heights: Vec<u64>,
    levels: Vec<Level>
}

impl<G> PreflowPushState<G>
where
    G: GraphProp,
    G::NodeId: Eq + Hash,
{
    pub fn push(&mut self, edge: usize, flow: u64) {
        self.edge_flow[edge] += flow;

        let source = self.edge_sources[edge].unwrap();
        self.excess[self.node_map[&source]] -= flow;
    }

    pub fn activate(&mut self, v: usize) {
        if v != self.s && v != self.t {
            let height = self.heights[v];

            self.levels[height as usize].activate(v);
        }
    }

    pub fn deactivate(&mut self, v: usize) {
        if v != self.s && v != self.t {
            let height = self.heights[v];

            self.levels[height as usize].deactivate(v);
        }
    }

    pub fn relabel(&self) {

    }

    pub fn discharge(&self) {

    }
}

pub fn preflow_push<G, F, C, W, E>(
    graph: G,
    s: G::NodeId,
    t: G::NodeId,
    mut capacity: C,
) -> Result<Option<MaxFlowReturn<G>>, E>
where
    G: NodeIndexable
        + IntoNodeReferences
        + IntoEdgeReferences
        + IntoEdgesDirected
        + NodeCount
        + EdgeCount
        + GraphProp
        + GraphBase
        + std::fmt::Debug
        + Visitable,
    <G as GraphBase>::NodeId: Eq + Hash + std::fmt::Debug,
    <G as GraphBase>::EdgeId: Eq + Hash + std::fmt::Debug,
    C: FnMut(&G::EdgeWeight) -> Result<u64, E>
{
    let num_edges = graph.edge_count();
    let num_nodes = graph.node_count();

    if num_nodes == 0 || num_edges == 0 {
        return Ok(None)
    }

    // Look-up tables for nodes and their supply / demand.
    let node_indices: Vec<G::NodeId> = graph.node_identifiers().collect();
    let node_map: HashMap<G::NodeId, usize> = node_indices
        .iter()
        .enumerate()
        .map(|(index, val)| (*val, index))
        .collect();

    // Look-up tables for the arcs and their attributes.
    let mut edge_indices: Vec<G::EdgeId> = Vec::with_capacity(num_edges);
    let mut edge_capacities: Vec<u64> = Vec::with_capacity(num_edges);
    let mut edge_sources: Vec<Option<G::NodeId>> = Vec::with_capacity(num_edges);
    let mut edge_targets: Vec<Option<G::NodeId>> = Vec::with_capacity(num_edges);
    
    for edge in graph.edge_references() {
        let capacity = capacity(edge.weight())?;
        let source = edge.source();
        let target = edge.target();

        // TODO: Shouldn't a negative capacity just result in an exception?
        if capacity < 0 {
            return Ok(None);
        }

        if source != target && capacity != 0 {
            edge_indices.push(edge.id());
            edge_capacities.push(capacity);
            edge_sources.push(Some(source));
            edge_targets.push(Some(target));
        }
    }

    let edge_map: HashMap<G::EdgeId, usize> = edge_indices
        .iter()
        .enumerate()
        .map(|(index, val)| (*val, index))
        .collect();

    let mut state: PreflowPushState<G> = PreflowPushState{
        s: node_map[&s],
        t: node_map[&t],
        num_nodes: num_nodes,
        num_edges: num_edges,
        node_map: node_map,
        edge_map: edge_map,
        edge_capacities: edge_capacities,
        edge_sources: edge_sources,
        edge_targets: edge_targets,
        edge_flow: vec![0; num_edges],
        excess: vec![0; num_nodes],
        heights: vec![0u64; num_nodes],
        levels: vec![Level::new(); 2 * num_nodes]
    };

    // Initialize the preflow for all edges.
    for edge in graph.edges_directed(s, petgraph::Direction::Outgoing) {
        let idx = state.edge_map[&edge.id()];
        state.push(idx, state.edge_capacities[idx]);
    }

    // Perform a reversed breadth-first search to determine the node heights.
    let mut queue: VecDeque<(<G as GraphBase>::NodeId, u64)> = VecDeque::from([(s, 0u64)]);
    let mut visited = graph.visit_map();

    while !queue.is_empty() {
        let (w, height) = queue.pop_front().unwrap();
        let height = height + 1;
        visited.visit(w);

        for edge in graph.edges_directed(w, petgraph::Direction::Incoming) {
            let v = edge.source();

            if !visited.is_visited(&v) {
                state.heights[state.node_map[&v]] = height;
                queue.push_back((v, height));
            }
        }
    }

    // If t is not reachable from s in the residual network, the maximum flow has to be a zero flow.
    if !visited.is_visited(&t) {
        return Ok(None);
    }

    let max_height = state.heights
        .iter()
        .enumerate().filter(|(n, h)| *n != state.node_map[&s])
        .map(|(_, h)| h)
        .max();
    state.heights[state.node_map[&s]] = state.num_nodes as u64;

    // Populate the level data structure for efficiently finding active and inactive nodes.
    for v in 0..num_nodes {
        if v != state.node_map[&s] && v != state.node_map[&t] {
            if state.excess[v] > 0 {
                state.activate(v);
            } else {
                state.deactivate(v);
            }
        }
    }

    Ok(None)
}