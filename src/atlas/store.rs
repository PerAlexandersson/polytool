//! The JSONL storage of the atlas ([`AtlasStore`]).
//!
//! A data directory holds `nodes.jsonl`, `edges.jsonl` and `links.jsonl`.
//! Every line is one JSON record; files are only appended to.  A node id may
//! occur on several lines (revisions); the last line of an id is its current
//! record.  Edges and links are written once (a repeated id is not written
//! again).  [`AtlasStore::open`] reads everything; [`AtlasStore::commit`]
//! appends the records added since, in a deterministic order: new and
//! revised nodes in the order they were first touched, then edges and links
//! in the order they were added.

use super::error::AtlasError;
use super::types::{AtlasEdge, AtlasLink, AtlasNode};
use std::collections::{BTreeMap, HashMap};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// File name of the node records.
pub const NODES_FILE: &str = "nodes.jsonl";
/// File name of the edge records.
pub const EDGES_FILE: &str = "edges.jsonl";
/// File name of the link records.
pub const LINKS_FILE: &str = "links.jsonl";

/// The records of an atlas directory, with the pending additions.
#[derive(Debug, Clone, Default)]
pub struct AtlasStore {
    directory: PathBuf,
    nodes: Vec<AtlasNode>,
    node_position: HashMap<String, usize>,
    edges: Vec<AtlasEdge>,
    edge_position: HashMap<String, usize>,
    links: Vec<AtlasLink>,
    link_position: HashMap<String, usize>,
    node_lines: usize,
    pending_nodes: Vec<String>,
    pending_edges: usize,
    pending_links: usize,
}

fn read_lines<T: serde::de::DeserializeOwned>(path: &Path) -> Result<Vec<T>, AtlasError> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let file = File::open(path)?;
    let mut out = Vec::new();
    for (number, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let record = serde_json::from_str(&line).map_err(|error| {
            AtlasError::InvalidRecord(format!("{}:{}: {error}", path.display(), number + 1))
        })?;
        out.push(record);
    }
    Ok(out)
}

fn append_lines<T: serde::Serialize>(path: &Path, records: &[&T]) -> Result<(), AtlasError> {
    if records.is_empty() {
        return Ok(());
    }
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    let mut text = String::new();
    for record in records {
        text.push_str(&serde_json::to_string(record)?);
        text.push('\n');
    }
    file.write_all(text.as_bytes())?;
    file.flush()?;
    Ok(())
}

impl AtlasStore {
    /// Read the records of `directory` (missing files are empty).
    ///
    /// # Errors
    ///
    /// [`AtlasError::Io`] or [`AtlasError::InvalidRecord`] for unreadable
    /// files or malformed lines.
    pub fn open(directory: impl AsRef<Path>) -> Result<Self, AtlasError> {
        let directory = directory.as_ref().to_path_buf();
        let node_records: Vec<AtlasNode> = read_lines(&directory.join(NODES_FILE))?;
        let edges: Vec<AtlasEdge> = read_lines(&directory.join(EDGES_FILE))?;
        let links: Vec<AtlasLink> = read_lines(&directory.join(LINKS_FILE))?;
        let mut store = Self {
            directory,
            node_lines: node_records.len(),
            ..Self::default()
        };
        for node in node_records {
            match store.node_position.get(&node.id) {
                Some(&i) => {
                    if node.revision <= store.nodes[i].revision {
                        return Err(AtlasError::InvalidRecord(format!(
                            "node {} revision {} does not increase",
                            node.id, node.revision
                        )));
                    }
                    store.nodes[i] = node;
                }
                None => {
                    store
                        .node_position
                        .insert(node.id.clone(), store.nodes.len());
                    store.nodes.push(node);
                }
            }
        }
        for edge in edges {
            if store.edge_position.contains_key(&edge.id) {
                return Err(AtlasError::InvalidRecord(format!(
                    "edge {} occurs twice",
                    edge.id
                )));
            }
            store
                .edge_position
                .insert(edge.id.clone(), store.edges.len());
            store.edges.push(edge);
        }
        for link in links {
            if store.link_position.contains_key(&link.id) {
                return Err(AtlasError::InvalidRecord(format!(
                    "link {} occurs twice",
                    link.id
                )));
            }
            store
                .link_position
                .insert(link.id.clone(), store.links.len());
            store.links.push(link);
        }
        Ok(store)
    }

    /// The data directory.
    pub fn directory(&self) -> &Path {
        &self.directory
    }

    /// The current record of every node, in the order of first appearance.
    pub fn nodes(&self) -> &[AtlasNode] {
        &self.nodes
    }

    /// Every edge, in file order (then pending ones).
    pub fn edges(&self) -> &[AtlasEdge] {
        &self.edges
    }

    /// Every link, in file order (then pending ones).
    pub fn links(&self) -> &[AtlasLink] {
        &self.links
    }

    /// Number of node lines in the file when it was opened (revisions
    /// included).
    pub fn node_lines(&self) -> usize {
        self.node_lines
    }

    /// The current record of a node.
    pub fn node(&self, id: &str) -> Option<&AtlasNode> {
        self.node_position.get(id).map(|&i| &self.nodes[i])
    }

    /// An edge by id.
    pub fn edge(&self, id: &str) -> Option<&AtlasEdge> {
        self.edge_position.get(id).map(|&i| &self.edges[i])
    }

    /// Whether a link id is present.
    pub fn has_link(&self, id: &str) -> bool {
        self.link_position.contains_key(id)
    }

    /// Insert a new node or a new revision of an existing one (its revision
    /// is set to one more than the current record's).
    pub fn put_node(&mut self, mut node: AtlasNode) {
        match self.node_position.get(&node.id) {
            Some(&i) => {
                if self.nodes[i] == node {
                    return;
                }
                node.revision = self.nodes[i].revision + 1;
                if !self.pending_nodes.contains(&node.id) {
                    self.pending_nodes.push(node.id.clone());
                } else {
                    // Several changes before one commit make one revision.
                    node.revision = self.nodes[i].revision;
                }
                self.nodes[i] = node;
            }
            None => {
                node.revision = 0;
                self.pending_nodes.push(node.id.clone());
                self.node_position.insert(node.id.clone(), self.nodes.len());
                self.nodes.push(node);
            }
        }
    }

    /// Add an edge; returns `false` when its id is already present.
    pub fn put_edge(&mut self, edge: AtlasEdge) -> bool {
        if self.edge_position.contains_key(&edge.id) {
            return false;
        }
        self.edge_position.insert(edge.id.clone(), self.edges.len());
        self.edges.push(edge);
        self.pending_edges += 1;
        true
    }

    /// Add a link; returns `false` when its id is already present.
    pub fn put_link(&mut self, link: AtlasLink) -> bool {
        if self.link_position.contains_key(&link.id) {
            return false;
        }
        self.link_position.insert(link.id.clone(), self.links.len());
        self.links.push(link);
        self.pending_links += 1;
        true
    }

    /// Number of records waiting to be committed (nodes, edges, links).
    pub fn pending(&self) -> (usize, usize, usize) {
        (
            self.pending_nodes.len(),
            self.pending_edges,
            self.pending_links,
        )
    }

    /// Append the pending records to the files (creating the directory and
    /// files when needed).
    ///
    /// # Errors
    ///
    /// [`AtlasError::Io`] when a file cannot be written.
    pub fn commit(&mut self) -> Result<(), AtlasError> {
        std::fs::create_dir_all(&self.directory)?;
        let nodes: Vec<&AtlasNode> = self
            .pending_nodes
            .iter()
            .filter_map(|id| self.node(id))
            .collect();
        append_lines(&self.directory.join(NODES_FILE), &nodes)?;
        let edges: Vec<&AtlasEdge> = self.edges[self.edges.len() - self.pending_edges..]
            .iter()
            .collect();
        append_lines(&self.directory.join(EDGES_FILE), &edges)?;
        let links: Vec<&AtlasLink> = self.links[self.links.len() - self.pending_links..]
            .iter()
            .collect();
        append_lines(&self.directory.join(LINKS_FILE), &links)?;
        self.node_lines += self.pending_nodes.len();
        self.pending_nodes.clear();
        self.pending_edges = 0;
        self.pending_links = 0;
        Ok(())
    }

    /// Nodes indexed by fingerprint.
    pub fn fingerprint_index(&self) -> BTreeMap<String, Vec<String>> {
        let mut index: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for node in &self.nodes {
            index
                .entry(node.fingerprint.clone())
                .or_default()
                .push(node.id.clone());
        }
        index
    }
}
