import pandas as pd
import re
import sys
import numpy as np
from scipy import stats
import argparse

def parse_gaf_reads(gaf_file, node_lengths, min_read_coverage=0.9, min_node_coverage=0.9, max_read_node_ratio=10.0):
    """
    Parse GAF file and count reads mapped to each node.
    For reads mapping to multiple nodes, add one count to each node.
    Tracks mapping positions for ALL nodes (single and multi-node mappings).
    
    Filters out alignments where:
    1. Neither the read nor node is well-covered
    2. A single-node mapping has read much longer than node (likely false positive)
    
    GAF columns:
    0: query name
    1: query length
    2: query start
    3: query end
    4: strand
    5: target path
    6: target length
    7: target start
    8: target end
    
    Returns: 
        node_read_counts: dict {node_id: read_count}
        node_positions: dict {node_id: list of (start, end) positions}
    """
    node_read_counts = {}
    node_positions = {}
    total_lines = 0
    filtered_coverage = 0
    filtered_length_ratio = 0
    
    with open(gaf_file, 'r') as f:
        for line in f:
            if line.strip():
                total_lines += 1
                fields = line.strip().split('\t')
                
                try:
                    # Parse alignment information
                    read_length = int(fields[1])
                    read_start = int(fields[2])
                    read_end = int(fields[3])
                    path = fields[5]
                    target_length = int(fields[6])
                    target_start = int(fields[7])
                    target_end = int(fields[8])
                    
                    # Calculate coverage fractions
                    read_aligned_length = read_end - read_start
                    read_coverage = read_aligned_length / read_length if read_length > 0 else 0
                    
                    target_aligned_length = target_end - target_start
                    target_coverage = target_aligned_length / target_length if target_length > 0 else 0
                    
                    # Filter 1: require EITHER read OR node to be well-covered
                    if read_coverage < min_read_coverage and target_coverage < min_node_coverage:
                        filtered_coverage += 1
                        continue
                    
                    # Parse nodes in path
                    clean_path = path.replace('<', ' ').replace('>', ' ')
                    nodes = clean_path.split()
                    
                    # Filter 2: For single-node mappings, check if read is much longer than node
                    if len(nodes) == 1 and target_length > 0:
                        read_to_node_ratio = read_length / target_length
                        if read_to_node_ratio > max_read_node_ratio:
                            filtered_length_ratio += 1
                            continue
                    
                except (IndexError, ValueError):
                    # Skip malformed lines
                    continue
                
                # Track positions for ALL nodes (single and multi-node mappings)
                # Calculate where the alignment falls within each node
                cumulative_end = 0
                for i, node_id in enumerate(nodes):
                    node_len = node_lengths.get(node_id, 0)
                    if node_len == 0:
                        # Skip nodes with unknown length
                        continue
                    
                    cumulative_start = cumulative_end
                    cumulative_end = cumulative_start + node_len
                    
                    # Find overlap between alignment [target_start, target_end] and this node
                    align_start_in_node = max(0, target_start - cumulative_start)
                    align_end_in_node = min(node_len, target_end - cumulative_start)
                    
                    # Only track if there's actual overlap
                    if align_end_in_node > align_start_in_node:
                        # Count the read for this node
                        node_read_counts[node_id] = node_read_counts.get(node_id, 0) + 1
                        
                        # Track position
                        if node_id not in node_positions:
                            node_positions[node_id] = []
                        node_positions[node_id].append((align_start_in_node, align_end_in_node))
    
    total_filtered = filtered_coverage + filtered_length_ratio
    print(f"  Total alignment lines: {total_lines}")
    print(f"  Filtered due to poor coverage: {filtered_coverage} ({100*filtered_coverage/total_lines:.1f}%)")
    print(f"  Filtered due to read >> node (single-node): {filtered_length_ratio} ({100*filtered_length_ratio/total_lines:.1f}%)")
    print(f"  Total filtered: {total_filtered} ({100*total_filtered/total_lines:.1f}%)")
    print(f"  Alignments retained: {total_lines - total_filtered}")
    
    return node_read_counts, node_positions

def parse_gfa_lengths(gfa_file):
    """
    Parse GFA file and extract node lengths.
    Returns dict: {node_id: length}
    """
    node_lengths = {}
    
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):
                fields = line.strip().split('\t')
                node_id = fields[1]
                
                # Find the LN:i: field
                for field in fields:
                    if field.startswith('LN:i:'):
                        length = int(field.split(':')[-1])
                        node_lengths[node_id] = length
                        break
    
    return node_lengths

def parse_gfa_edges(gfa_file):
    """
    Parse L (link) lines from GFA to build adjacency graph.
    Returns: dict {node_id: set of nodes that can follow it}
    """
    edges = {}
    
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('L'):
                fields = line.strip().split('\t')
                from_node = fields[1]
                from_orient = fields[2]
                to_node = fields[3]
                to_orient = fields[4]
                
                # For simplicity, consider + orientation (forward strand)
                # Store directed edges
                if from_node not in edges:
                    edges[from_node] = set()
                edges[from_node].add(to_node)
                
                # Also store reverse edge for bidirectional traversal
                if to_node not in edges:
                    edges[to_node] = set()
                # Note: not adding reverse here as GFA links are directional
    
    return edges

def parse_bubble_core_nodes(bubble_file):
    """
    Parse bubble file and extract core nodes (first and last node of each bubble).
    Returns set of core node IDs.
    """
    core_nodes = set()
    
    with open(bubble_file, 'r') as f:
        for line in f:
            if line.strip():
                fields = line.strip().split('\t')
                # Column 12 (index 11) contains the node path
                if len(fields) > 11:
                    node_path = fields[11]
                    nodes = node_path.split(',')
                    
                    if len(nodes) >= 2:
                        # First and last nodes are core nodes
                        core_nodes.add(nodes[0])
                        core_nodes.add(nodes[-1])
    
    return core_nodes

def calculate_coverage_uniformity(positions, node_length, bin_size=100):
    """
    Calculate coverage uniformity across a node using hybrid approach.
    
    For short nodes (<5000bp): Uses position-by-position depth calculation
    For long nodes (≥5000bp): Uses binning approach
    
    Returns:
        uniformity_score: 0-1, where 1 is perfectly uniform
        coverage_fraction: fraction of node/bins with coverage
        avg_depth: average depth across the node (reads per base)
    """
    if not positions or node_length == 0:
        return 0.0, 0.0, 0.0
    
    # For short nodes, calculate depth position-by-position
    if node_length < 5000:
        coverage_array = [0] * node_length
        for start, end in positions:
            for pos in range(max(0, start), min(node_length, end)):
                coverage_array[pos] += 1
        
        if not coverage_array:
            return 0.0, 0.0, 0.0
        
        avg_depth = np.mean(coverage_array)
        
        # Calculate uniformity
        if avg_depth > 0:
            std_coverage = np.std(coverage_array)
            cv = std_coverage / avg_depth
            uniformity_score = np.exp(-cv)
        else:
            uniformity_score = 0.0
        
        # Coverage fraction
        coverage_fraction = sum(1 for c in coverage_array if c > 0) / len(coverage_array)
        
        return uniformity_score, coverage_fraction, avg_depth
    
    # For longer nodes, use binning approach
    num_bins = node_length // bin_size
    bin_coverage = [0] * num_bins
    
    # For each position span, add coverage to overlapping bins
    for start, end in positions:
        start_bin = max(0, min(int(start / bin_size), num_bins - 1))
        end_bin = max(0, min(int((end - 1) / bin_size), num_bins - 1))
        
        for b in range(start_bin, end_bin + 1):
            if b < num_bins:
                bin_coverage[b] += 1
    
    if sum(bin_coverage) == 0:
        return 0.0, 0.0, 0.0
    
    # Average depth is mean of bin coverages (since each bin represents same number of bases)
    avg_depth = np.mean(bin_coverage)
    
    # Calculate uniformity using coefficient of variation
    std_coverage = np.std(bin_coverage)
    if avg_depth > 0:
        cv = std_coverage / avg_depth
        uniformity_score = np.exp(-cv)
    else:
        uniformity_score = 0.0
    
    # Calculate fraction of bins with coverage
    bins_with_coverage = sum(1 for c in bin_coverage if c > 0)
    coverage_fraction = bins_with_coverage / num_bins
    
    return uniformity_score, coverage_fraction, avg_depth

def filter_nodes_by_uniformity(node_read_counts, node_positions, node_lengths, 
                                min_uniformity=0.3, min_coverage_fraction=0.5):
    """
    Filter nodes based on coverage uniformity.
    Nodes with uneven coverage are set to 0 reads.
    Also calculates true average depth from coverage bins.
    
    Returns:
        filtered_counts: dict with filtered read counts
        uniformity_stats: dict with uniformity metrics per node
        node_depths: dict with average depth per node
    """
    filtered_counts = {}
    uniformity_stats = {}
    node_depths = {}
    
    for node_id, count in node_read_counts.items():
        length = node_lengths.get(node_id, 0)
        positions = node_positions.get(node_id, [])
        
        if length > 0 and positions:
            uniformity_score, coverage_fraction, avg_depth = calculate_coverage_uniformity(
                positions, length
            )
            
            node_depths[node_id] = avg_depth
            
            uniformity_stats[node_id] = {
                'uniformity_score': uniformity_score,
                'coverage_fraction': coverage_fraction,
                'avg_depth': avg_depth,
                'original_count': count,
                'node_length': length,
                'num_mapped_reads': len(positions)
            }
            
            # Check if node passes uniformity filters
            if uniformity_score >= min_uniformity and coverage_fraction >= min_coverage_fraction:
                filtered_counts[node_id] = count
                uniformity_stats[node_id]['passes_filter'] = True
            else:
                filtered_counts[node_id] = 0
                uniformity_stats[node_id]['passes_filter'] = False
        else:
            # No position data - use read count as rough estimate
            # This happens for multi-node mappings where we don't track positions
            # Assume reads span the entire node on average
            if length > 0:
                fallback_depth = count / length
            else:
                fallback_depth = 0
            
            node_depths[node_id] = fallback_depth
            filtered_counts[node_id] = count
            
            if node_id not in uniformity_stats:
                uniformity_stats[node_id] = {
                    'uniformity_score': 1.0,
                    'coverage_fraction': 1.0,
                    'avg_depth': fallback_depth,
                    'original_count': count,
                    'node_length': length,
                    'num_mapped_reads': 0,
                    'passes_filter': True
                }
    
    return filtered_counts, uniformity_stats, node_depths
    """
    Calculate RPKM (Reads Per Kilobase per Million mapped reads)
    RPKM = (reads * 10^9) / (length * total_reads)
    """
    if node_length == 0 or total_reads == 0:
        return 0.0
    return (read_count * 1e9) / (node_length * total_reads)

def calculate_rpkm(read_count, node_length, total_reads):
    """
    Calculate RPKM (Reads Per Kilobase per Million mapped reads)
    RPKM = (reads * 10^9) / (length * total_reads)
    """
    if node_length == 0 or total_reads == 0:
        return 0.0
    return (read_count * 1e9) / (node_length * total_reads)

def calculate_confidence_interval(values, confidence=0.95):
    """
    Calculate confidence interval for a set of values.
    Returns (mean, lower_ci, upper_ci, std_error)
    """
    if len(values) == 0:
        return 0, 0, 0, 0
    
    n = len(values)
    mean = np.mean(values)
    std_err = stats.sem(values)
    
    # Calculate confidence interval
    ci = std_err * stats.t.ppf((1 + confidence) / 2, n - 1)
    
    return mean, mean - ci, mean + ci, std_err

def classify_node_zygosity_by_depth(node_depth, core_avg_depth, read_count, min_reads):
    """
    Classify a node based on depth comparison to flanking core nodes.
    Uses pre-calculated average depths from coverage bins.
    
    Returns: 'het', 'hom', or 'absent'
    """
    if read_count < min_reads:
        return 'absent'
    
    if core_avg_depth == 0 or node_depth == 0:
        return 'absent'
    
    depth_ratio = node_depth / core_avg_depth
    
    # More lenient thresholds
    if depth_ratio > 0.6:  # Changed from 0.75 to 0.7
        return 'hom'
    elif depth_ratio > 0.3:  # Changed from 0.25 to 0.2
        return 'het'
    else:
        return 'absent'

def natural_sort_key(s):
    """
    Key function for natural sorting (s1, s2, s10 instead of s1, s10, s2)
    """
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]

def topological_sort_with_fallback(nodes, edges):
    """
    Sort nodes respecting topology where edges exist, 
    falling back to numerical order when topology is ambiguous.
    Uses Kahn's algorithm with proper edge tracking.
    """
    if not nodes:
        return []
    
    node_set = set(nodes)
    
    # Build subgraph with only these nodes
    subgraph = {n: [] for n in nodes}
    in_degree = {n: 0 for n in nodes}
    
    for node in nodes:
        if node in edges:
            for target in edges[node]:
                if target in node_set:
                    subgraph[node].append(target)
                    in_degree[target] += 1
    
    # Start with nodes that have no incoming edges, sorted numerically
    queue = sorted([n for n in nodes if in_degree[n] == 0], key=natural_sort_key)
    result = []
    
    while queue:
        # Pick numerically smallest from queue
        current = queue.pop(0)
        result.append(current)
        
        # Process neighbors in sorted order
        neighbors = sorted(subgraph[current], key=natural_sort_key)
        for neighbor in neighbors:
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                # Insert maintaining sorted order
                insert_pos = 0
                for i, node in enumerate(queue):
                    if natural_sort_key(neighbor) < natural_sort_key(node):
                        insert_pos = i
                        break
                    insert_pos = i + 1
                queue.insert(insert_pos, neighbor)
    
    # If there are remaining nodes (shouldn't happen with proper DAG), add them sorted
    remaining = sorted([n for n in nodes if n not in result], key=natural_sort_key)
    result.extend(remaining)
    
    return result

def partition_nodes_into_paths(alternate_nodes, edges, node_classifications, core_start, core_end):
    """
    Partition nodes into separate paths considering:
    1. Each hom node gets its own individual path (to allow proper insertion of het nodes)
    2. Het nodes connected to core nodes are separate  
    3. Other het nodes are grouped by connectivity
    
    Returns: list of path objects with 'nodes' and 'classification' keys
    """
    classified_nodes = {}
    for node in alternate_nodes:
        classification = node_classifications.get(node, 'absent')
        if classification != 'absent':
            classified_nodes[node] = classification
    
    if not classified_nodes:
        return []
    
    # Separate hom and het nodes
    hom_nodes = [n for n in classified_nodes if node_classifications[n] == 'hom']
    het_nodes = [n for n in classified_nodes if node_classifications[n] == 'het']
    
    all_paths = []
    
    # Process hom nodes - each gets its own INDIVIDUAL path
    # This allows het nodes to be inserted between them during ordering
    for node in hom_nodes:
        all_paths.append({
            'nodes': [node],
            'classification': 'hom'
        })
    core_connected = set()
    # Process het nodes
    if het_nodes:
        # Find which het nodes connect directly to core nodes
        core_connected = set()
        for node in het_nodes:
            if core_start in edges and node in edges[core_start]:
                core_connected.add(node)
            if node in edges and core_end in edges[node]:
                core_connected.add(node)
        truly_core_only = set()
        for node in core_connected:
            has_other_connections = False
            
            # Check outgoing edges
            #if node in edges:
            #    for target in edges[node]:
            #        if target in het_nodes and target != node:
            #            has_other_connections = True
            #            break
            
            # Check incoming edges (reverse direction)
            if not has_other_connections:
                for other_node in het_nodes:
                    if other_node != node and other_node in edges:
                        if node in edges[other_node]:
                            has_other_connections = True
                            break
            
            # Only mark as core-connected if it has NO other connections
            if not has_other_connections:
                truly_core_only.add(node)
        
        # Core-only connected nodes get individual paths with special classification
        for node in sorted(truly_core_only, key=natural_sort_key):
            all_paths.append({
                'nodes': [node],
                'classification': 'het_core_connected'
            })
        
        # Other het nodes (including core-connected with other connections)
        other_het = [n for n in het_nodes if n not in truly_core_only]
        # Core-connected get individual paths
        #for node in sorted(core_connected, key=natural_sort_key):
        #    all_paths.append({
        #        'nodes': [node],
        #        'classification': 'het_core_connected'
        #    })
        
        # Other het nodes - build connected components
        #other_het = [n for n in het_nodes if n not in core_connected]
        if other_het:
            het_subgraph = {}
            for node in other_het:
                het_subgraph[node] = set()
                if node in edges:
                    for target in edges[node]:
                        if target in other_het:
                            het_subgraph[node].add(target)
            
            # Build connected components
            visited = set()
            def dfs_het(node, component):
                visited.add(node)
                component.add(node)
                # Check both forward and reverse edges
                for neighbor in het_subgraph.get(node, set()):
                    if neighbor not in visited:
                        dfs_het(neighbor, component)
                # Check reverse
                for other in other_het:
                    if node in het_subgraph.get(other, set()) and other not in visited:
                        dfs_het(other, component)
            
            for node in other_het:
                if node not in visited:
                    component = set()
                    dfs_het(node, component)
                    sorted_comp = topological_sort_with_fallback(list(component), het_subgraph)
                    all_paths.append({
                        'nodes': sorted_comp,
                        'classification': 'het'
                    })
#    if core_connected:
#        print(core_connected,all_paths)
    return all_paths

def order_paths_globally(paths, edges, core_start, core_end):
    """
    Order paths by checking connections between last node of one path 
    and first node of next path. Falls back to numerical order.
    
    Returns: ordered list of all nodes
    """
    if not paths:
        return []
    
    if len(paths) == 1:
        return paths[0]
    
    # Build adjacency: which path connects to which path
    # path i connects to path j if: last node of path_i connects to first node of path_j
    path_connections = {}
    for i in range(len(paths)):
        path_connections[i] = []
        last_node_i = paths[i][-1]
        
        if last_node_i in edges:
            for j in range(len(paths)):
                if i != j:
                    first_node_j = paths[j][0]
                    if first_node_j in edges[last_node_i]:
                        path_connections[i].append(j)
    
    # Also check connections from core_start
    paths_from_core = []
    for i in range(len(paths)):
        first_node = paths[i][0]
        if core_start in edges and first_node in edges[core_start]:
            paths_from_core.append(i)
    
    # Topological sort with path connections
    in_degree = {i: 0 for i in range(len(paths))}
    for i in range(len(paths)):
        for j in path_connections[i]:
            in_degree[j] += 1
    
    # Start with paths that have no incoming connections
    # Prioritize those connected to core_start
    queue = []
    for i in range(len(paths)):
        if in_degree[i] == 0:
            if i in paths_from_core:
                queue.insert(0, i)  # Core-connected go first
            else:
                queue.append(i)
    
    # Sort queue by numerical order within priority groups
    if len(queue) > 1:
        core_count = sum(1 for i in queue if i in paths_from_core)
        if core_count > 0:
            core_queue = sorted([i for i in queue if i in paths_from_core], 
                               key=lambda i: natural_sort_key(paths[i][0]))
            other_queue = sorted([i for i in queue if i not in paths_from_core],
                                key=lambda i: natural_sort_key(paths[i][0]))
            queue = core_queue + other_queue
        else:
            queue = sorted(queue, key=lambda i: natural_sort_key(paths[i][0]))
    
    result = []
    ordered = []
    
    while queue:
        current_idx = queue.pop(0)
        ordered.append(current_idx)
        result.extend(paths[current_idx])
        
        # Add paths that this one connects to
        for next_idx in path_connections[current_idx]:
            in_degree[next_idx] -= 1
            if in_degree[next_idx] == 0:
                # Insert in queue maintaining numerical order
                insert_pos = len(queue)
                for i, q_idx in enumerate(queue):
                    if natural_sort_key(paths[next_idx][0]) < natural_sort_key(paths[q_idx][0]):
                        insert_pos = i
                        break
                queue.insert(insert_pos, next_idx)
    
    # Add any remaining paths (no connections) in numerical order
    remaining = [i for i in range(len(paths)) if i not in ordered]
    remaining.sort(key=lambda i: natural_sort_key(paths[i][0]))
    for i in remaining:
        result.extend(paths[i])
    
    return result

def assign_paths_to_haplotypes(path_objs, edges, core_start, core_end):
    """
    Assign paths to haplotypes based on their classifications.
    path_objs: list of dicts with 'nodes' and 'classification' keys
    
    Returns: (hap1_nodes, hap2_nodes) - complete ordered lists
    """
    if not path_objs:
        return [], []
    
    hap1_path_objs = []
    hap2_path_objs = []
    hap1_has_core_connected = False
    hap2_has_core_connected = False
    
    def count_connections_to_haplotype(het_path, hap_paths, edges):
        """
        Count how many connections a het path has to nodes in a haplotype.
        Checks first and last nodes of the het path against first and last nodes of all paths in haplotype.
        Returns (connection_count, occupied_connections)
        occupied_connections is a set of (from_node, to_node) tuples already used in haplotype
        """
        if not hap_paths:
            return 0, set()
        
        connection_count = 0
        het_first = het_path[0]
        het_last = het_path[-1]
        
        # Track which connections are already occupied within the haplotype
        occupied = set()
        
        # First, find all edges already used WITHIN the haplotype
        for i, hap_path_i in enumerate(hap_paths):
            for j, hap_path_j in enumerate(hap_paths):
                if i != j:
                    # Check if last of path_i connects to first of path_j
                    last_i = hap_path_i[-1]
                    first_j = hap_path_j[0]
                    if last_i in edges and first_j in edges[last_i]:
                        occupied.add((last_i, first_j))
        # Now count available connections to the het path
        for hap_path in hap_paths:
            hap_first = hap_path[0]
            hap_last = hap_path[-1]
            
            # Check if het_last connects to hap_first (het path comes before)
            if het_last in edges and hap_first in edges[het_last]:
                #if (het_last, hap_first) not in occupied:
                if not any(b == hap_first for a, b in occupied):
                    connection_count += 1
            
            # Check if hap_last connects to het_first (het path comes after)
            if hap_last in edges and het_first in edges[hap_last]:
                #if (hap_last, het_first) not in occupied:
                if not any(a == hap_last for a, b in occupied):
                    connection_count += 1
            
            # Check if het_first connects to hap_first
            if het_first in edges and hap_first in edges[het_first]:
                #if (het_first, hap_first) not in occupied:
                if not any(b == hap_first for a, b in occupied):
                    connection_count += 1
            
            # Check if hap_first connects to het_first
            if hap_first in edges and het_first in edges[hap_first]:
                #if (hap_first, het_first) not in occupied:
                if not any(a == hap_first for a, b in occupied):
                    connection_count += 1
            
            # Check if het_last connects to hap_last
            if het_last in edges and hap_last in edges[het_last]:
                #if (het_last, hap_last) not in occupied:
                if not any(b == hap_last for a, b in occupied):
                    connection_count += 1
            
            # Check if hap_last connects to het_last
            if hap_last in edges and het_last in edges[hap_last]:
                #if (hap_last, het_last) not in occupied:
                if not any(a == hap_last for a, b in occupied):
                    connection_count += 1
        
        return connection_count, occupied

    for path_obj in path_objs:
        nodes = path_obj['nodes']
        classification = path_obj['classification']
        
        if classification == 'hom':
            # Homozygous - add to BOTH haplotypes
            hap1_path_objs.append(nodes)
            hap2_path_objs.append(nodes)
        
        elif classification == 'het_core_connected':
            # Core-connected het nodes - try to distribute
            # Add to smaller haplotype
            #if len(hap1_path_objs) <= len(hap2_path_objs):
            #    hap1_path_objs.append(nodes)
            #else:
            #    hap2_path_objs.append(nodes)
            if not hap1_has_core_connected:
                hap1_path_objs.append(nodes)
                hap1_has_core_connected = True
            elif not hap2_has_core_connected:
                hap2_path_objs.append(nodes)
                hap2_has_core_connected = True
            #else:
                # Both already have one, add to smaller haplotype
                #if len(hap1_path_objs) <= len(hap2_path_objs):
                #    hap1_path_objs.append(nodes)
                #else:
                #    hap2_path_objs.append(nodes)
            else:
                # Both already have one, check connectivity
                connections_to_hap1, _ = count_connections_to_haplotype(nodes, hap1_path_objs, edges)
                connections_to_hap2, _ = count_connections_to_haplotype(nodes, hap2_path_objs, edges)

                if connections_to_hap1 > connections_to_hap2:
                    hap1_path_objs.append(nodes)
                elif connections_to_hap2 > connections_to_hap1:
                    hap2_path_objs.append(nodes)
                else:
                    # Equal connections, add to smaller haplotype
                    if len(hap1_path_objs) <= len(hap2_path_objs):
                        hap1_path_objs.append(nodes)
                    else:
                        hap2_path_objs.append(nodes)
        elif classification == 'het':
            # Regular het path - add to smaller haplotype
            #if len(hap1_path_objs) <= len(hap2_path_objs):
            #    hap1_path_objs.append(nodes)
            #else:
            #    hap2_path_objs.append(nodes)
            connections_to_hap1, _ = count_connections_to_haplotype(nodes, hap1_path_objs, edges)
            connections_to_hap2, _ = count_connections_to_haplotype(nodes, hap2_path_objs, edges)
            
            if connections_to_hap1 > connections_to_hap2:
                hap1_path_objs.append(nodes)
            elif connections_to_hap2 > connections_to_hap1:
                hap2_path_objs.append(nodes)
            else:
                # Equal connections (including both zero), add to smaller haplotype
                if len(hap1_path_objs) >= len(hap2_path_objs):
                    hap1_path_objs.append(nodes)
                else:
                    hap2_path_objs.append(nodes)

    # Order paths globally for each haplotype
    hap1_nodes = order_paths_globally(hap1_path_objs, edges, core_start, core_end)
    hap2_nodes = order_paths_globally(hap2_path_objs, edges, core_start, core_end)
    
    return hap1_nodes, hap2_nodes

def build_haplotype_paths(bubble_nodes, node_depths_dict, node_length_dict, node_read_dict, 
                          edges, min_reads):
    """
    Build paths for each haplotype based on depth comparison of non-core nodes.
    Uses pre-calculated average depths from coverage bins.
    Respects graph topology from GFA edges.
    
    Returns (hap1_nodes, hap2_nodes, total_length_hap1, total_length_hap2)
    """
    if len(bubble_nodes) < 2:
        return [], [], 0, 0
    
    # First and last are core nodes, middle are alternate alleles
    core_start = bubble_nodes[0]
    core_end = bubble_nodes[-1]
    alternate_nodes = bubble_nodes[1:-1]
    
    if not alternate_nodes:
        return [], [], 0, 0
    
    # Get core node depths
    core_start_depth = node_depths_dict.get(core_start, 0)
    core_end_depth = node_depths_dict.get(core_end, 0)
    core_avg_depth = (core_start_depth + core_end_depth) / 2
    
    # Classify each alternate node
    node_classifications = {}
    for node in alternate_nodes:
        read_count = node_read_dict.get(node, 0)
        node_depth = node_depths_dict.get(node, 0)
        
        classification = classify_node_zygosity_by_depth(
            node_depth, core_avg_depth, read_count, min_reads
        )
        node_classifications[node] = classification
    
    # Partition nodes into paths with classifications
    path_objs = partition_nodes_into_paths(alternate_nodes, edges, node_classifications, 
                                           core_start, core_end)
    
    # Assign paths to haplotypes
    hap1_nodes, hap2_nodes = assign_paths_to_haplotypes(path_objs, edges, core_start, core_end)
    
    # Calculate total lengths
    hap1_length = sum(node_length_dict.get(n, 0) for n in hap1_nodes)
    hap2_length = sum(node_length_dict.get(n, 0) for n in hap2_nodes)
    
    return hap1_nodes, hap2_nodes, hap1_length, hap2_length

def generate_bed_output(bubble_file, node_rpkm_dict, node_length_dict, node_read_dict, node_depths_dict,
                        edges, sample_name, ploidy, output_prefix, min_reads):
    """
    Generate BED file(s) based on bubble file and depth classifications.
    Uses bubble-specific depth thresholds for each bubble.
    Respects graph topology when building haplotype paths.
    """
    if ploidy == 'haploid':
        output_files = [f"{output_prefix}.bed"]
        file_handles = [open(output_files[0], 'w')]
    else:  # diploid
        output_files = [f"{output_prefix}_hap1.bed", f"{output_prefix}_hap2.bed"]
        file_handles = [open(output_files[0], 'w'), open(output_files[1], 'w')]
    
    # Track bubble-specific thresholds for reporting
    bubble_thresholds = []
    
    try:
        with open(bubble_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    
                    # Extract fields from bubble file
                    chrom = fields[0]
                    start = fields[1]
                    end = fields[2]
                    node_path = fields[11]
                    nodes = node_path.split(',')
                    
                    if len(nodes) < 2:
                        continue
                    
                    core_start = nodes[0]
                    core_end = nodes[-1]
                    
                    # Track bubble stats
                    bubble_thresholds.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'core_start': core_start,
                        'core_end': core_end,
                        'core_start_depth': node_depths_dict.get(core_start, 0),
                        'core_end_depth': node_depths_dict.get(core_end, 0),
                        'avg_core_depth': (node_depths_dict.get(core_start, 0) + node_depths_dict.get(core_end, 0)) / 2
                    })
                    
                    # Build haplotype paths with topology awareness
                    hap1_nodes, hap2_nodes, hap1_length, hap2_length = build_haplotype_paths(
                        nodes, node_depths_dict, node_length_dict, node_read_dict,
                        edges, min_reads
                    )
                    
                    if ploidy == 'haploid':
                        # For haploid, use hap1 path
                        if hap1_nodes:
                            alt_path = '>' + '>'.join(hap1_nodes)
                        else:
                            alt_path = '*'
                        
                        # Format: chr start end >core_start >core_end alt_path:length:+:contig:bubble_start:bubble_end
                        bed_line = f"{chrom}\t{start}\t{end}\t>{core_start}\t>{core_end}\t{alt_path}:{hap1_length}:+:{sample_name}:{start}:{end}\n"
                        file_handles[0].write(bed_line)
                    
                    else:  # diploid
                        # Haplotype 1
                        if hap1_nodes:
                            alt_path_1 = '>' + '>'.join(hap1_nodes)
                        else:
                            alt_path_1 = '*'
                        bed_line_1 = f"{chrom}\t{start}\t{end}\t>{core_start}\t>{core_end}\t{alt_path_1}:{hap1_length}:+:{sample_name}#1:{start}:{end}\n"
                        file_handles[0].write(bed_line_1)
                        
                        # Haplotype 2
                        if hap2_nodes:
                            alt_path_2 = '>' + '>'.join(hap2_nodes)
                        else:
                            alt_path_2 = '*'
                        bed_line_2 = f"{chrom}\t{start}\t{end}\t>{core_start}\t>{core_end}\t{alt_path_2}:{hap2_length}:+:{sample_name}#2:{start}:{end}\n"
                        file_handles[1].write(bed_line_2)
    
    finally:
        for fh in file_handles:
            fh.close()
    
    return output_files, bubble_thresholds

def main(gaf_file, gfa_file, bubble_file, sample_name, ploidy, output_prefix, min_reads,
         min_uniformity, min_coverage_fraction, min_read_cov, min_node_cov, max_read_node_ratio):
    """
    Main processing pipeline
    """
    print("="*60)
    print("STEP 1: Parsing input files")
    print("="*60)
    
    print("Parsing GFA file first (needed for node lengths)...")
    node_lengths = parse_gfa_lengths(gfa_file)
    
    print("Parsing GFA edges (for topology awareness)...")
    edges = parse_gfa_edges(gfa_file)
    print(f"  Total edges parsed: {sum(len(v) for v in edges.values())}")
    
    print("Parsing GAF file...")
    print(f"  Min read coverage: {min_read_cov*100:.0f}%")
    print(f"  Min node coverage: {min_node_cov*100:.0f}%")
    print(f"  (Alignments must meet at least ONE of these thresholds)")
    print(f"  Max read/node length ratio (single-node): {max_read_node_ratio}x")
    node_read_counts, node_positions = parse_gaf_reads(gaf_file, node_lengths, min_read_cov, min_node_cov, max_read_node_ratio)
    
    print("Parsing bubble file for core nodes...")
    core_nodes = parse_bubble_core_nodes(bubble_file)
    print(f"Total core nodes identified: {len(core_nodes)}")
    
    print("\n" + "="*60)
    print("STEP 2: Filtering nodes by coverage uniformity")
    print("="*60)
    print(f"Minimum uniformity score: {min_uniformity}")
    print(f"Minimum coverage fraction: {min_coverage_fraction}")
    
    # Filter nodes based on uniformity
    filtered_read_counts, uniformity_stats, node_depths = filter_nodes_by_uniformity(
        node_read_counts, node_positions, node_lengths,
        min_uniformity, min_coverage_fraction
    )
    
    # Count how many nodes were filtered
    nodes_filtered = sum(1 for nid in filtered_read_counts 
                        if filtered_read_counts[nid] == 0 and node_read_counts[nid] > 0)
    print(f"Nodes filtered due to uneven coverage: {nodes_filtered}")
    
    # Save uniformity statistics
    uniformity_df = pd.DataFrame.from_dict(uniformity_stats, orient='index')
    uniformity_df.index.name = 'node_id'
    uniformity_df.to_csv(f'{output_prefix}_uniformity_stats.csv')
    print(f"Uniformity statistics saved to: {output_prefix}_uniformity_stats.csv")
    
    # Get all unique nodes from both files
    all_nodes = set(filtered_read_counts.keys()) | set(node_lengths.keys())
    
    # Calculate total mapped reads (using filtered counts)
    total_reads = sum(filtered_read_counts.values())
    print(f"Total read counts after uniformity filtering: {total_reads}")
    
    print("\n" + "="*60)
    print("STEP 3: Calculating RPKM values")
    print("="*60)
    print(f"Minimum read filter: {min_reads}")
    
    # Create dataframe
    data = []
    node_rpkm_dict = {}
    
    for node_id in sorted(all_nodes):
        read_count = filtered_read_counts.get(node_id, 0)
        length = node_lengths.get(node_id, 0)
        
        # Calculate RPKM
        rpkm = calculate_rpkm(read_count, length, total_reads)
        
        # Apply minimum read filter
        if read_count >= min_reads:
            node_rpkm_dict[node_id] = rpkm
        else:
            node_rpkm_dict[node_id] = 0  # Treat as absent
        
        # Mark if node is a core node
        is_core = node_id in core_nodes
        
        # Get uniformity info
        uniformity_info = uniformity_stats.get(node_id, {})
        
        data.append({
            'node_id': node_id,
            'read_count': read_count,
            'original_read_count': node_read_counts.get(node_id, 0),
            'node_length': length,
            'RPKM': rpkm,
            'avg_depth': node_depths.get(node_id, 0),
            'is_core_node': is_core,
            'uniformity_score': uniformity_info.get('uniformity_score', 1.0),
            'coverage_fraction': uniformity_info.get('coverage_fraction', 1.0),
            'passes_uniformity': uniformity_info.get('passes_filter', True)
        })
    
    df = pd.DataFrame(data)
    
    print("\n" + "="*60)
    print("STEP 4: Calculating core node statistics")
    print("="*60)
    
    # Filter core nodes > 150bp for average RPKM calculation
    core_nodes_filtered = df[(df['is_core_node'] == True) & (df['node_length'] > 150)]
    print(f"Core nodes > 150bp: {len(core_nodes_filtered)}")
    
    # Calculate average RPKM and confidence interval for core nodes (for reference)
    if len(core_nodes_filtered) > 0:
        core_rpkm_values = core_nodes_filtered['RPKM'].values
        mean_rpkm, lower_ci, upper_ci, std_err = calculate_confidence_interval(core_rpkm_values)
        
        print(f"\nGlobal Core Node RPKM Statistics (for reference):")
        print(f"  Number of core nodes: {len(core_nodes_filtered)}")
        print(f"  Average RPKM: {mean_rpkm:.4f}")
        print(f"  95% Confidence Interval: [{lower_ci:.4f}, {upper_ci:.4f}]")
        print(f"  Standard Error: {std_err:.4f}")
        print(f"\nNote: Bubble-specific depth thresholds will be calculated for each bubble")
        print(f"      based on its adjacent core nodes' average depth")
        print(f"      Het threshold: depth_ratio > 0.25, Hom threshold: depth_ratio > 0.75")
        
        # Add core node stats to a summary dict
        core_stats = {
            'n_core_nodes': len(core_nodes_filtered),
            'mean_rpkm': mean_rpkm,
            'lower_ci_95': lower_ci,
            'upper_ci_95': upper_ci,
            'std_error': std_err
        }
        
        # Save core statistics
        core_stats_df = pd.DataFrame([core_stats])
        core_stats_df.to_csv(f'{output_prefix}_core_stats.csv', index=False)
        print(f"\nGlobal core node statistics saved to: {output_prefix}_core_stats.csv")
        
        # Save filtered core nodes
        core_nodes_filtered.to_csv(f'{output_prefix}_core_nodes.csv', index=False)
        print(f"Filtered core nodes saved to: {output_prefix}_core_nodes.csv")
    else:
        print("\nWarning: No core nodes > 150bp found!")
        core_stats = None
    
    # Save all node data
    output_file = f'{output_prefix}_all_nodes.csv'
    df.to_csv(output_file, index=False)
    print(f"All node data saved to: {output_file}")
    
    print("\n" + "="*60)
    print("STEP 5: Generating BED output files")
    print("="*60)
    print(f"Sample: {sample_name}")
    print(f"Ploidy: {ploidy}")
    print(f"Using bubble-specific depth thresholds")
    print(f"Respecting graph topology from GFA edges")
    
    # Generate BED output
    output_files, bubble_thresholds = generate_bed_output(
        bubble_file, 
        node_rpkm_dict, 
        node_lengths,
        filtered_read_counts,  # Use filtered counts
        node_depths,  # Pass node depths
        edges,  # Pass GFA edges for topology
        sample_name, 
        ploidy,
        output_prefix,
        min_reads
    )
    
    # Save bubble-specific thresholds
    bubble_thresh_df = pd.DataFrame(bubble_thresholds)
    bubble_thresh_df.to_csv(f'{output_prefix}_bubble_thresholds.csv', index=False)
    print(f"\nBubble-specific thresholds saved to: {output_prefix}_bubble_thresholds.csv")
    
    print(f"\nBED file(s) generated:")
    for f in output_files:
        print(f"  - {f}")
    
    # Display statistics on bubble thresholds
    if len(bubble_thresholds) > 0:
        depths = [b['avg_core_depth'] for b in bubble_thresholds]
        print(f"\nBubble core depth statistics:")
        print(f"  Min avg core depth: {np.min(depths):.4f}")
        print(f"  Max avg core depth: {np.max(depths):.4f}")
        print(f"  Mean avg core depth: {np.mean(depths):.4f}")
        print(f"  Median avg core depth: {np.median(depths):.4f}")
    
    print("\n" + "="*60)
    print("COMPLETE!")
    print("="*60)
    
    return df, core_stats
