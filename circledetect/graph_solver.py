import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

class EccDNAGraphSolver:
    def __init__(self):
        self.adj_matrix = None
        self.id_map = {}
        
    def build_graph(self, validated_junctions: pd.DataFrame, discordant: pd.DataFrame):
        """
        Constructs a directed graph.
        Nodes: Validated Breakpoints.
        Edges: 
          1. Self-loops weighted by Posterior Prob (Single Segment evidence).
          2. Directed edges between nodes supported by Discordant Reads (Multi Segment evidence).
        """
        unique_sites = validated_junctions[['chr', 'pos']].drop_duplicates()
        unique_sites['id'] = range(len(unique_sites))
        
        self.id_map = {row['id']: (row['chr'], row['pos']) for _, row in unique_sites.iterrows()}
        reverse_map = {(row['chr'], row['pos']): row['id'] for _, row in unique_sites.iterrows()}
        
        n_nodes = len(unique_sites)
        if n_nodes == 0:
            return
            
        row_inds = []
        col_inds = []
        data = []
        
        # 1. Add Self-Loops (Evidence for Single-Segment Circles)
        for _, row in validated_junctions.iterrows():
            if row['is_valid']:
                key = (row['chr'], row['pos'])
                if key in reverse_map:
                    nid = reverse_map[key]
                    row_inds.append(nid)
                    col_inds.append(nid)
                    data.append(row['posterior_prob'])
        
        # 2. Add Edges from Discordant Reads (Evidence for Connections)
        for _, row in discordant.iterrows():
            k1 = (row['chr1'], row['pos1'])
            k2 = (row['chr2'], row['pos2'])
            
            if k1 in reverse_map and k2 in reverse_map:
                u = reverse_map[k1]
                v = reverse_map[k2]
                row_inds.append(u)
                col_inds.append(v)
                data.append(1.0) # Binary support from discordant pair
                
                # Also add reverse edge if strand orientation suggests it (simplified here to bidirectional support)
                # In strict mode, check strand1/strand2 to determine direction
        
        if not row_inds:
            self.adj_matrix = csr_matrix((0,0))
            return

        self.adj_matrix = csr_matrix((data, (row_inds, col_inds)), shape=(n_nodes, n_nodes))
        
    def detect_cycles(self) -> pd.DataFrame:
        """
        Identifies ecDNA circles via Strongly Connected Components.
        """
        if self.adj_matrix is None or self.adj_matrix.shape[0] == 0:
            return pd.DataFrame()
            
        n_comp, labels = connected_components(
            csgraph=self.adj_matrix, 
            directed=True, 
            connection='strong'
        )
        
        circles = []
        for i in range(n_comp):
            members = np.where(labels == i)[0]
            if len(members) == 0:
                continue
                
            # Calculate Confidence
            subgraph = self.adj_matrix[np.ix_(members, members)]
            if subgraph.nnz > 0:
                conf = subgraph.data.mean()
            else:
                conf = 0.0
            
            nodes = [self.id_map[m] for m in members]
            
            circle_type = 'Multi-Segment' if len(nodes) > 1 else 'Single-Segment'
            
            # Filter low confidence single segments
            if circle_type == 'Single-Segment' and conf < 0.9:
                continue
                
            circles.append({
                'type': circle_type,
                'segment_count': len(nodes),
                'nodes': str(nodes), # Store as string for CSV
                'confidence': conf
            })
                
        return pd.DataFrame(circles)