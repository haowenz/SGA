#ifndef SGA_SEQUENCEGRAPH_H
#define SGA_SEQUENCEGRAPH_H

#include <algorithm>
#include <fstream>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "gfa.h"
#include "sequence.h"
#include "utils.h"

namespace sga {

template <class GraphSizeType = int32_t, class QueryLengthType = int16_t,
          class ScoreType = int16_t>
struct VertexWithDistanceForDijkstra {
  GraphSizeType graph_vertex_id;
  QueryLengthType query_index;
  ScoreType distance;
  bool is_reverse_complementary;
};

template <class GraphSizeType = int32_t>
struct DijkstraAlgorithmStatistics {
  GraphSizeType forward_num_cells = 0;
  GraphSizeType rc_num_cells = 0;
};

template <class GraphSizeType = int32_t, class QueryLengthType = int16_t,
          class ScoreType = int16_t>
class SequenceGraph {
 public:
  SequenceGraph() {}
  ~SequenceGraph() {}

  GraphSizeType GetNumVerticesInCompactedGraph() {
    return compacted_graph_labels_.size();
  }

  GraphSizeType GetNumEdgesInCompactedGraph() {
    GraphSizeType num_edges = 0;
    for (std::vector<GraphSizeType> &neighbors :
         compacted_graph_adjacency_list_) {
      num_edges += neighbors.size();
    }
    return num_edges;
  }

  GraphSizeType GetNumVertices() { return labels_.size(); }

  GraphSizeType GetNumEdges() {
    GraphSizeType num_edges = 0;
    for (std::vector<GraphSizeType> &neighbors : adjacency_list_) {
      num_edges += neighbors.size();
    }
    return num_edges;
  }

  void PrintLayer(const std::vector<ScoreType> &layer,
                  const std::vector<GraphSizeType> &order) {
    for (GraphSizeType i = 0; i < GetNumVertices(); ++i) {
      std::cerr << layer[order[i]] << " ";
    }
    std::cerr << std::endl;
  }

  void GenerateCompressedRepresentation() {
    GraphSizeType num_vertices = GetNumVertices();
    GraphSizeType num_edges = GetNumEdges();
    std::cerr << "# vertices: " << num_vertices << ", # edges: " << num_edges
              << std::endl;

    look_up_table_.reserve(GetNumVertices());
    look_up_table_.push_back(0);
    for (auto &neighbor_list : adjacency_list_) {
      GraphSizeType last_sum = look_up_table_.back();
      look_up_table_.push_back(neighbor_list.size());
      look_up_table_.back() += last_sum;
      neighbor_table_.insert(neighbor_table_.end(), neighbor_list.begin(),
                             neighbor_list.end());
    }
  }

  void SetAlignmentParameters(const ScoreType substitution_penalty,
                              const ScoreType deletion_penalty,
                              const ScoreType insertion_penalty) {
    substitution_penalty_ = substitution_penalty;
    deletion_penalty_ = deletion_penalty;
    insertion_penalty_ = insertion_penalty;
  }

  void AddReverseComplementaryVertexIfNecessary(
      const gfa_t *gfa_graph, uint32_t gfa_vertex_id,
      std::vector<GraphSizeType> &reverse_complementary_compacted_vertex_id) {
    const uint32_t gfa_segment_id = gfa_vertex_id >> 1;
    const GraphSizeType vertex_id = gfa_segment_id + 1;
    const bool is_gfa_segment_forward = ((gfa_vertex_id & 1) == 0);
    const bool vertex_rc_is_not_added =
        (reverse_complementary_compacted_vertex_id[vertex_id] == vertex_id);

    // std::cerr << "gfa_vi: " << gfa_vertex_id << " si: " << gfa_segment_id
    //          << " rcvi: "
    //          << reverse_complementary_compacted_vertex_id[vertex_id]
    //          << " vi: " << vertex_id << " forward:" << is_gfa_segment_forward
    //          << std::endl;
    if (!is_gfa_segment_forward && vertex_rc_is_not_added) {
      const gfa_seg_t &gfa_segment = gfa_graph->seg[gfa_segment_id];
      std::string rc_sequence_bases;
      for (QueryLengthType qi = 0; qi < (QueryLengthType)gfa_segment.len;
           ++qi) {
        rc_sequence_bases.push_back(
            base_complement_[(int)gfa_segment.seq
                                 [(QueryLengthType)gfa_segment.len - 1 - qi]]);
      }

      compacted_graph_labels_.emplace_back(rc_sequence_bases);
      compacted_graph_adjacency_list_.emplace_back(
          std::vector<GraphSizeType>());

      const uint32_t new_compacted_graph_rc_vertex_id =
          compacted_graph_adjacency_list_.size() - 1;
      reverse_complementary_compacted_vertex_id[vertex_id] =
          new_compacted_graph_rc_vertex_id;
    }
  }

  void AddEdge(const gfa_t *gfa_graph, const gfa_arc_t &gfa_arc,
               const std::vector<GraphSizeType>
                   &reverse_complementary_compacted_vertex_id) {
    const uint32_t gfa_head_vertex_id = gfa_arc_head(gfa_arc);
    const uint32_t gfa_head_segment_id = gfa_head_vertex_id >> 1;
    const bool is_gfa_head_segment_forward = ((gfa_head_vertex_id & 1) == 0);

    GraphSizeType head_vertex_id = gfa_head_segment_id + 1;
    if (!is_gfa_head_segment_forward) {
      head_vertex_id =
          reverse_complementary_compacted_vertex_id[head_vertex_id];
    }

    const uint32_t gfa_tail_vertex_id = gfa_arc_tail(gfa_arc);
    const uint32_t gfa_tail_segment_id = gfa_tail_vertex_id >> 1;
    const bool is_gfa_tail_segment_forward = ((gfa_tail_vertex_id & 1) == 0);

    GraphSizeType tail_vertex_id = gfa_tail_segment_id + 1;
    if (!is_gfa_tail_segment_forward) {
      tail_vertex_id =
          reverse_complementary_compacted_vertex_id[tail_vertex_id];
    }

    compacted_graph_adjacency_list_[head_vertex_id].emplace_back(
        tail_vertex_id);
  }

  // An algorithm to convert a GFA to a compacted sequence graph with two passes
  // on arcs and one pass on segments. We assume there is no overlap. The
  // algorithm first passes the arcs to know which node needs to be duplicated
  // for a reverse complement. Then the segments are parsed to create the
  // vertices and labels. Finally the arcs are parsed again to add the edges.
  void LoadFromGfaFile(const std::string &graph_file_path) {
    gfa_t *gfa_graph = gfa_read(graph_file_path.data());
    // gfa_print(gfa_graph, stderr, 0);

    std::vector<GraphSizeType> reverse_complementary_compacted_vertex_id;

    const uint32_t num_segments = gfa_graph->n_seg;
    reverse_complementary_compacted_vertex_id.reserve(num_segments);
    compacted_graph_labels_.reserve(num_segments);
    compacted_graph_adjacency_list_.reserve(num_segments);

    // Add a dummy vertex.
    compacted_graph_labels_.emplace_back("N");
    reverse_complementary_compacted_vertex_id.emplace_back(0);
    compacted_graph_adjacency_list_.emplace_back(std::vector<GraphSizeType>());

    for (uint32_t si = 0; si < num_segments; ++si) {
      const gfa_seg_t &current_segment = gfa_graph->seg[si];
      compacted_graph_labels_.emplace_back(std::string(current_segment.seq));

      const uint32_t compacted_graph_vertex_id = si + 1;
      reverse_complementary_compacted_vertex_id.emplace_back(
          compacted_graph_vertex_id);

      compacted_graph_adjacency_list_.emplace_back(
          std::vector<GraphSizeType>());
    }

    const uint64_t num_arcs = gfa_graph->n_arc;
    // std::cerr << "NUM ARCS: " << num_arcs << std::endl;
    const gfa_arc_t *gfa_arcs = gfa_graph->arc;

    for (uint64_t ai = 0; ai < num_arcs; ++ai) {
      if (gfa_arcs[ai].del || gfa_arcs[ai].comp) continue;
      const uint32_t gfa_head_vertex_id = gfa_arc_head(gfa_arcs[ai]);
      AddReverseComplementaryVertexIfNecessary(
          gfa_graph, gfa_head_vertex_id,
          reverse_complementary_compacted_vertex_id);

      const uint32_t gfa_tail_vertex_id = gfa_arc_tail(gfa_arcs[ai]);
      AddReverseComplementaryVertexIfNecessary(
          gfa_graph, gfa_tail_vertex_id,
          reverse_complementary_compacted_vertex_id);
    }

    for (uint64_t ai = 0; ai < num_arcs; ++ai) {
      if (gfa_arcs[ai].del || gfa_arcs[ai].comp) continue;
      AddEdge(gfa_graph, gfa_arcs[ai],
              reverse_complementary_compacted_vertex_id);
    }
  }

  void LoadFromTxtFile(const std::string &graph_file_path) {
    std::string line;
    std::ifstream infile(graph_file_path);

    GraphSizeType num_vertices = 0;
    GraphSizeType row_index = 0;

    compacted_graph_labels_.emplace_back("N");
    compacted_graph_adjacency_list_.emplace_back(std::vector<GraphSizeType>());

    while (std::getline(infile, line)) {
      std::istringstream inputString(line);
      // get count of vertices from header row
      if (row_index == 0) {
        inputString >> num_vertices;
        compacted_graph_labels_.reserve(num_vertices);
      } else {  // get out-neighbor vertex ids and vertex label
        assert(row_index <= num_vertices);
        compacted_graph_adjacency_list_.emplace_back(
            std::vector<GraphSizeType>());

        // Parse the input line
        std::vector<std::string> tokens(
            std::istream_iterator<std::string>{inputString},
            std::istream_iterator<std::string>());
        assert(tokens.size() > 0);
        compacted_graph_labels_.emplace_back(tokens.back());
        for (auto it = tokens.begin();
             it != tokens.end() && std::next(it) != tokens.end(); it++) {
          compacted_graph_adjacency_list_.back().emplace_back(stoi(*it) + 1);
        }
      }
      row_index++;
    }

    std::cerr << "# vertices in compacted graph: "
              << GetNumVerticesInCompactedGraph()
              << ", # edges in compacted graph: "
              << GetNumEdgesInCompactedGraph() << std::endl;
  }

  void OutputCompactedGraphInGFA(std::string &output_file_path) {
    std::ofstream outstrm(output_file_path);
    outstrm << "H\tVN:Z:1.0\n";
    // Both for loops start from 1 to skip the dummy vertex.
    for (uint32_t i = 1; i < compacted_graph_labels_.size(); ++i) {
      outstrm << "S\t" << i << "\t" << compacted_graph_labels_[i] << "\n";
    }

    for (uint32_t i = 1; i < compacted_graph_adjacency_list_.size(); ++i) {
      for (auto neighbor : compacted_graph_adjacency_list_[i]) {
        outstrm << "L\t" << i << "\t+\t" << neighbor << "\t+\t0M\n";
      }
    }
  }

  void GenerateCharLabeledGraph() {
    for (const std::string &compacted_graph_label : compacted_graph_labels_) {
      labels_.emplace_back(compacted_graph_label[0]);
      adjacency_list_.emplace_back(std::vector<GraphSizeType>());
    }

    // Keep the original vertex ids unchanged.
    GraphSizeType vertex_id = compacted_graph_labels_.size();
    GraphSizeType compacted_graph_vertex_id = 0;

    for (const std::string &compacted_graph_label : compacted_graph_labels_) {
      GraphSizeType compacted_graph_label_length =
          compacted_graph_label.length();

      if (compacted_graph_label_length == 1) {
        // Add the neighbors of the chain to the neighbors of the last vertex.
        adjacency_list_[compacted_graph_vertex_id].insert(
            adjacency_list_[compacted_graph_vertex_id].end(),
            compacted_graph_adjacency_list_[compacted_graph_vertex_id].begin(),
            compacted_graph_adjacency_list_[compacted_graph_vertex_id].end());
      }

      for (GraphSizeType i = 1; i < compacted_graph_label_length; ++i) {
        labels_.emplace_back(compacted_graph_label[i]);
        adjacency_list_.emplace_back(std::vector<GraphSizeType>());

        // If this is the second vertex in the chain, add the link from the
        // first vertex to the second vertex in the chain.
        if (i == 1) {
          adjacency_list_[compacted_graph_vertex_id].push_back(vertex_id);
        }

        // If this is not the last vertex in the chain, which means it has a
        // next vertex, add link from current vertex to its next.
        if (i + 1 < compacted_graph_label_length) {
          adjacency_list_[vertex_id].push_back(vertex_id + 1);
        } else {
          // If this is the last vertex in the chain, add the neighbors of the
          // chain to the neighbors of the last vertex.
          adjacency_list_[vertex_id].insert(
              adjacency_list_[vertex_id].end(),
              compacted_graph_adjacency_list_[compacted_graph_vertex_id]
                  .begin(),
              compacted_graph_adjacency_list_[compacted_graph_vertex_id].end());
        }

        ++vertex_id;
      }

      ++compacted_graph_vertex_id;
    }
    // after the loop, compacted_graph_vertex_id should be the number of
    // vertices in the compacted graph
    // assert(compacted_graph_vertex_id == GetNumVerticesInCompactedGraph());
    // add an edge between the dummy and every other vertex
    // adjacency_list_[0].reserve(vertex_id - 1);
    // for (GraphSizeType i = 1; i < vertex_id; ++i) {
    //  adjacency_list_[0].push_back(i);
    //}
  }

  inline char GetReverseComplementaryVertexLabel(GraphSizeType vertex) const {
    return base_complement_[(int)labels_[vertex]];
  }

  inline char GetVertexLabel(GraphSizeType vertex) const {
    return labels_[vertex];
  }

  void GenerateReverseComplementaryCharLabeledGraph() {
    reverse_complementary_adjacency_list_.assign(adjacency_list_.size(),
                                                 std::vector<GraphSizeType>());
    for (GraphSizeType vertex = 0; vertex < adjacency_list_.size(); ++vertex) {
      for (const GraphSizeType neighbor : adjacency_list_[vertex]) {
        reverse_complementary_adjacency_list_[neighbor].push_back(vertex);
      }
    }
  }

  void PropagateInsertions(const std::vector<ScoreType> &initialized_layer,
                           const std::vector<GraphSizeType> &initialized_order,
                           std::vector<ScoreType> &current_layer,
                           std::vector<GraphSizeType> &current_order) {
    const GraphSizeType num_vertices = GetNumVertices();
    GraphSizeType initialized_order_index = 0;
    GraphSizeType current_order_index = 0;
    visited_.assign(num_vertices, false);

    std::deque<GraphSizeType> updated_neighbors;

    while (initialized_order_index < num_vertices ||
           !updated_neighbors.empty()) {
      GraphSizeType min_vertex = num_vertices;

      if (initialized_order_index < num_vertices &&
          (updated_neighbors.empty() ||
           current_layer[initialized_order[initialized_order_index]] <
               current_layer[updated_neighbors.front()])) {
        min_vertex = initialized_order[initialized_order_index];
        ++initialized_order_index;
      } else {
        min_vertex = updated_neighbors.front();
        updated_neighbors.pop_front();
      }

      if (!visited_[min_vertex]) {
        visited_[min_vertex] = true;
        current_order[current_order_index] = min_vertex;
        ++current_order_index;
        for (const auto &neighbor : adjacency_list_[min_vertex]) {
          if (!visited_[neighbor] &&
              current_layer[neighbor] >
                  current_layer[min_vertex] + insertion_penalty_) {
            current_layer[neighbor] =
                current_layer[min_vertex] + insertion_penalty_;
            updated_neighbors.push_back(neighbor);
          }
        }
      }
    }
  }

  // Build the order look up table.
  void BuildOrderLookUpTable(const std::vector<ScoreType> &previous_layer,
                             const std::vector<GraphSizeType> &previous_order) {
    const GraphSizeType num_vertices = GetNumVertices();
    GraphSizeType match_index = 0, substitution_index = 0, deletion_index = 0;
    GraphSizeType count = 0;

    while (count < 3 * num_vertices) {
      // Find the min.
      ScoreType min_distance =
          previous_layer[previous_order[num_vertices - 1]] +
          substitution_penalty_ + deletion_penalty_ + 1;
      int min_type = -1;  // 0 for match, 1 for substitution, 2 for deletion.
      if (match_index < num_vertices &&
          previous_layer[previous_order[match_index]] < min_distance) {
        min_distance = previous_layer[previous_order[match_index]];
        min_type = 0;
      }

      if (substitution_index < num_vertices &&
          previous_layer[previous_order[substitution_index]] +
                  substitution_penalty_ <
              min_distance) {
        min_distance = previous_layer[previous_order[substitution_index]] +
                       substitution_penalty_;
        min_type = 1;
      }

      if (deletion_index < num_vertices &&
          previous_layer[previous_order[deletion_index]] + deletion_penalty_ <
              min_distance) {
        min_distance =
            previous_layer[previous_order[deletion_index]] + deletion_penalty_;
        min_type = 2;
      }

      // Put the order of min into the look up table.
      if (min_type == 0) {
        order_look_up_table_[min_type * num_vertices +
                             previous_order[match_index]] = count;
        ++match_index;
      } else if (min_type == 1) {
        order_look_up_table_[min_type * num_vertices +
                             previous_order[substitution_index]] = count;
        ++substitution_index;
      } else {
        order_look_up_table_[min_type * num_vertices +
                             previous_order[deletion_index]] = count;
        ++deletion_index;
      }

      ++count;
    }
  }

  void InitializeDistancesWithSorting(
      const char sequence_base, const std::vector<ScoreType> &previous_layer,
      const std::vector<GraphSizeType> &previous_order,
      std::vector<ScoreType> &initialized_layer,
      std::vector<GraphSizeType> &initialized_order) {
    const GraphSizeType num_vertices = GetNumVertices();

    // Initialize the layer.
    initialized_layer[0] = previous_layer[0] + deletion_penalty_;
    for (GraphSizeType j = 1; j < num_vertices; ++j) {
      ScoreType cost = 0;
      if (sequence_base != labels_[j]) {
        cost = substitution_penalty_;
      }
      initialized_layer[j] = previous_layer[0] + cost;
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      if (initialized_layer[i] > previous_layer[i] + deletion_penalty_) {
        initialized_layer[i] = previous_layer[i] + deletion_penalty_;
      }

      for (const auto &neighbor : adjacency_list_[i]) {
        ScoreType cost = 0;

        if (sequence_base != labels_[neighbor]) {
          cost = substitution_penalty_;
        }

        if (initialized_layer[neighbor] > previous_layer[i] + cost) {
          initialized_layer[neighbor] = previous_layer[i] + cost;
        }
      }
    }

    // Use sorting to get order.
    for (GraphSizeType vertex = 0; vertex < num_vertices; ++vertex) {
      distances_with_vertices_[vertex] =
          std::make_pair((*initialized_layer)[vertex], vertex);
    }
    std::sort(distances_with_vertices_.begin(), distances_with_vertices_.end());
    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      initialized_order[i] = distances_with_vertices_[i].second;
    }
  }

  void InitializeDistances(const char sequence_base,
                           const std::vector<ScoreType> &previous_layer,
                           const std::vector<GraphSizeType> &previous_order,
                           std::vector<ScoreType> &initialized_layer,
                           std::vector<GraphSizeType> &initialized_order) {
    const GraphSizeType num_vertices = GetNumVertices();
    BuildOrderLookUpTable(previous_layer, previous_order);

    // Initialize the layer
    initialized_layer[0] = previous_layer[0] + deletion_penalty_;
    parents_[0] = 0;
    types_[0] = 2;

    for (GraphSizeType j = 1; j < num_vertices; ++j) {
      ScoreType cost = 0;
      int type = 0;
      if (sequence_base != labels_[j]) {
        cost = substitution_penalty_;
        type = 1;
      }
      initialized_layer[j] = previous_layer[0] + cost;
      parents_[j] = 0;
      types_[j] = type;
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      if (initialized_layer[i] > previous_layer[i] + deletion_penalty_) {
        initialized_layer[i] = previous_layer[i] + deletion_penalty_;
        parents_[i] = i;
        types_[i] = 2;
      }

      for (const auto &neighbor : adjacency_list_[i]) {
        ScoreType cost = 0;
        int type = 0;

        if (sequence_base != labels_[neighbor]) {
          cost = substitution_penalty_;
          type = 1;
        }

        if (initialized_layer[neighbor] > previous_layer[i] + cost) {
          initialized_layer[neighbor] = previous_layer[i] + cost;
          parents_[neighbor] = i;
          types_[neighbor] = type;
        }
      }
    }

    // Get the order. One should notice there can be multiple vertices share the
    // same parent and type.
    order_offsets_.assign(3 * num_vertices + 1, 0);
    order_counts_.assign(3 * num_vertices, 0);

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      order_offsets_
          [order_look_up_table_[types_[i] * num_vertices + parents_[i]] + 1]++;
    }

    for (GraphSizeType i = 1; i < 3 * num_vertices + 1; ++i) {
      order_offsets_[i] += order_offsets_[i - 1];
    }

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      const GraphSizeType order =
          order_look_up_table_[types_[i] * num_vertices + parents_[i]];
      initialized_order[order_offsets_[order] + order_counts_[order]] = i;
      ++order_counts_[order];
    }
  }

  ScoreType AlignUsingLinearGapPenalty(const sga::Sequence &sequence) {
    ScoreType max_cost = std::max(
        std::max(substitution_penalty_, deletion_penalty_), insertion_penalty_);
    const GraphSizeType num_vertices = GetNumVertices();
    const QueryLengthType sequence_length = sequence.GetLength();
    const std::string &sequence_bases = sequence.GetSequence();

    std::vector<ScoreType> previous_layer(num_vertices);
    std::vector<GraphSizeType> previous_order(num_vertices);
    std::vector<ScoreType> initialized_layer(num_vertices,
                                             sequence_length * max_cost + 1);
    std::vector<GraphSizeType> initialized_order(num_vertices);
    std::vector<ScoreType> current_layer(num_vertices, 0);
    std::vector<GraphSizeType> current_order;
    current_order.reserve(num_vertices);

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      current_order.push_back(i);
    }

    order_look_up_table_.assign(3 * num_vertices, 0);
    parents_.assign(num_vertices, 0);
    types_.assign(num_vertices, 0);
    // distances_with_vertices_.assign(num_vertices, std::make_pair(0,0));

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      std::swap(previous_order, current_order);
      InitializeDistances(sequence_bases[i], previous_layer, previous_order,
                          initialized_layer, initialized_order);
      // InitializeDistancesWithSorting(sequence_bases[i], previous_layer,
      // previous_order, &initialized_layer, &initialized_order);
      current_layer = initialized_layer;
      PropagateInsertions(initialized_layer, initialized_order, current_layer,
                          current_order);
    }

    const ScoreType forward_alignment_cost = current_layer[current_order[0]];

    // For reverse complement.
    initialized_layer.assign(num_vertices, sequence_length * max_cost + 1);
    current_layer.assign(num_vertices, 0);

    for (GraphSizeType i = 0; i < num_vertices; ++i) {
      current_order[i] = i;
    }

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      std::swap(previous_order, current_order);
      InitializeDistances(
          base_complement_[(int)sequence_bases[sequence_length - 1 - i]],
          previous_layer, previous_order, initialized_layer, initialized_order);
      // InitializeDistancesWithSorting(base_complement_[sequence_bases[sequence_length
      // - 1 - i]], previous_layer, previous_order, &initialized_layer,
      // &initialized_order);
      current_layer = initialized_layer;
      PropagateInsertions(initialized_layer, initialized_order, current_layer,
                          current_order);
    }

    const ScoreType reverse_complement_alignment_cost =
        current_layer[current_order[0]];

    const ScoreType min_alignment_cost =
        std::min(forward_alignment_cost, reverse_complement_alignment_cost);

    std::cerr << "Sequence length: " << sequence_length
              << ", forward alignment cost:" << forward_alignment_cost
              << ", reverse complement alignment cost:"
              << reverse_complement_alignment_cost
              << ", alignment cost:" << min_alignment_cost << std::endl;
    return min_alignment_cost;
  }

  void PropagateWithNavarroAlgorithm(
      const GraphSizeType from, const GraphSizeType to,
      GraphSizeType &num_propagations,
      std::vector<QueryLengthType> &current_layer) {
    num_propagations += 1;
    if (current_layer[to] > insertion_penalty_ + current_layer[from]) {
      current_layer[to] = insertion_penalty_ + current_layer[from];
      for (const auto &neighbor : adjacency_list_[to]) {
        PropagateWithNavarroAlgorithm(to, neighbor, num_propagations,
                                      current_layer);
      }
    }
  }

  void ComputeLayerWithNavarroAlgorithm(
      const char sequence_base,
      const std::vector<QueryLengthType> &previous_layer,
      GraphSizeType &num_propagations,
      std::vector<QueryLengthType> &current_layer) {
    const GraphSizeType num_vertices = GetNumVertices();

    // Initialize current layer
    current_layer[0] = previous_layer[0] + deletion_penalty_;
    for (GraphSizeType j = 1; j < num_vertices; ++j) {
      QueryLengthType cost = 0;

      if (sequence_base != labels_[j]) {
        cost = substitution_penalty_;
      }

      current_layer[j] = previous_layer[0] + cost;
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      if (current_layer[i] > previous_layer[i] + deletion_penalty_) {
        current_layer[i] = previous_layer[i] + deletion_penalty_;
      }

      for (const auto &neighbor : adjacency_list_[i]) {
        QueryLengthType cost = 0;

        if (sequence_base != labels_[neighbor]) {
          cost = substitution_penalty_;
        }

        if (current_layer[neighbor] > previous_layer[i] + cost) {
          current_layer[neighbor] = previous_layer[i] + cost;
        }
      }
    }

    for (GraphSizeType i = 1; i < num_vertices; ++i) {
      for (const auto &neighbor : adjacency_list_[i]) {
        PropagateWithNavarroAlgorithm(i, neighbor, num_propagations,
                                      current_layer);
      }
    }
  }

  QueryLengthType AlignUsingLinearGapPenaltyWithNavarroAlgorithm(
      const sga::Sequence &sequence) {
    QueryLengthType max_cost = std::max(
        std::max(substitution_penalty_, deletion_penalty_), insertion_penalty_);
    const GraphSizeType num_vertices = GetNumVertices();
    const QueryLengthType sequence_length = sequence.GetLength();
    const std::string &sequence_bases = sequence.GetSequence();
    std::vector<QueryLengthType> previous_layer(num_vertices,
                                                sequence_length * max_cost + 1);
    std::vector<QueryLengthType> current_layer(num_vertices, 0);

    GraphSizeType num_propagations = 0;

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      ComputeLayerWithNavarroAlgorithm(sequence_bases[i], previous_layer,
                                       num_propagations, current_layer);
    }

    const QueryLengthType forward_alignment_cost =
        *std::min_element(current_layer.begin(), current_layer.end());

    // For reverse complement.
    previous_layer.assign(num_vertices, sequence_length * max_cost + 1);
    current_layer.assign(num_vertices, 0);

    for (QueryLengthType i = 0; i < sequence_length; ++i) {
      std::swap(previous_layer, current_layer);
      ComputeLayerWithNavarroAlgorithm(
          base_complement_[(int)sequence_bases[sequence_length - 1 - i]],
          previous_layer, num_propagations, current_layer);
    }

    const QueryLengthType reverse_complement_alignment_cost =
        *std::min_element(current_layer.begin(), current_layer.end());

    const QueryLengthType min_alignment_cost =
        std::min(forward_alignment_cost, reverse_complement_alignment_cost);

    std::cerr << "Sequence length: " << sequence_length
              << ", forward alignment cost:" << forward_alignment_cost
              << ", reverse complement alignment cost:"
              << reverse_complement_alignment_cost
              << ", alignment cost:" << min_alignment_cost
              << ", num propogations: " << num_propagations << std::endl;
    return min_alignment_cost;
  }

  ScoreType AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
      const sga::Sequence &sequence,
      DijkstraAlgorithmStatistics<GraphSizeType> &stats) {
    const GraphSizeType num_vertices = GetNumVertices();
    const QueryLengthType sequence_length = sequence.GetLength();
    const std::string &sequence_bases = sequence.GetSequence();

    std::vector<std::unordered_map<QueryLengthType, ScoreType>>
        forward_vertex_distances(num_vertices);
    std::vector<std::unordered_map<QueryLengthType, ScoreType>>
        complementary_vertex_distances(num_vertices);

    // std::vector<std::unordered_map<
    //    QueryLengthType, VertexWithDistanceForDijkstra<
    //                         GraphSizeType, QueryLengthType, ScoreType>>>
    //    vertex_parent(num_vertices);

    auto compare_function =
        [](const VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
                                               ScoreType> &v1,
           const VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
                                               ScoreType> &v2) {
          if (v1.distance > v2.distance) {
            return true;
          }

          if (v1.distance == v2.distance) {
            if (v1.query_index < v2.query_index) {
              return true;
            }
            if (v1.query_index == v2.query_index) {
              if (v1.graph_vertex_id > v2.graph_vertex_id) {
                return true;
              }
              if (v1.graph_vertex_id == v2.graph_vertex_id) {
                if (!v1.is_reverse_complementary &&
                    v2.is_reverse_complementary) {
                  return true;
                }
                return false;
              }
              return false;
            }
            return false;
          }
          return false;
        };

    std::priority_queue<VertexWithDistanceForDijkstra<
                            GraphSizeType, QueryLengthType, ScoreType>,
                        std::vector<VertexWithDistanceForDijkstra<
                            GraphSizeType, QueryLengthType, ScoreType>>,
                        decltype(compare_function)>
        Q(compare_function);

    stats.forward_num_cells = 0;
    stats.rc_num_cells = 0;

    for (GraphSizeType vertex = 0; vertex < num_vertices; ++vertex) {
      // Deal with forward strand fisrt.
      ScoreType cost = 0;
      const char vertex_label = labels_[vertex];

      if (sequence_bases[0] != vertex_label) {
        cost = substitution_penalty_;
      }
      cost = std::min(cost, insertion_penalty_);

      Q.push({/*graph_vertex_id=*/vertex, /*query_index=*/0, /*distance=*/cost,
              /*is_reverse_complementary=*/false});

      ++(stats.forward_num_cells);
      forward_vertex_distances[vertex][0] = cost;

      // Now deal with reverse complementary strand.
      cost = 0;
      const char complementary_vertex_label =
          GetReverseComplementaryVertexLabel(vertex);

      if (sequence_bases[sequence_length - 1] != complementary_vertex_label) {
        cost = substitution_penalty_;
      }
      cost = std::min(cost, insertion_penalty_);

      Q.push({/*graph_vertex_id=*/vertex, /*query_index=*/0, /*distance=*/cost,
              /*is_reverse_complementary=*/true});
      complementary_vertex_distances[vertex][0] = cost;
      ++(stats.rc_num_cells);

      // vertex_parent[vertex][0] =
      //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
      //                                  ScoreType>{
      //        /*graph_vertex_id=*/vertex, /*query_index=*/0,
      //        /*distance=*/cost};
      // std::cerr << "Init PUSH: " << vertex << " " << 0 << " " << cost
      //          << std::endl;
    }

    ScoreType min_alignment_cost = 0;

    while (!Q.empty()) {
      const auto current_vertex = Q.top();
      Q.pop();

      // Check if we reach the last layer where we can stop.
      if (current_vertex.query_index + 1 == sequence_length) {
        min_alignment_cost = current_vertex.distance;
        // auto previous_it =
        //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
        //                                  ScoreType>{-1, -1, 0};
        // auto it = current_vertex;
        // while (!(it.graph_vertex_id == previous_it.graph_vertex_id &&
        //         it.query_index == previous_it.query_index)) {
        //  std::cerr << "Traceback: "
        //            << "gi: " << it.graph_vertex_id << " qi: " <<
        //            it.query_index
        //            << " d: " << it.distance
        //            << " qb: " << sequence_bases[it.query_index]
        //            << " gb: " << labels_[it.graph_vertex_id];

        //  if (it.distance + substitution_penalty_ == previous_it.distance) {
        //    std::cerr << " op: M";
        //  }

        //  if (it.graph_vertex_id == previous_it.graph_vertex_id &&
        //      it.distance + deletion_penalty_ == previous_it.distance) {
        //    std::cerr << " op: D";
        //  }

        //  if (it.query_index == previous_it.query_index) {
        //    std::cerr << " op: I";
        //  }
        //  std::cerr << std::endl;
        //  previous_it = it;
        //  it = vertex_parent[it.graph_vertex_id][it.query_index];
        //}
        // std::cerr << std::endl;
        break;
      }

      auto &vertex_distances = current_vertex.is_reverse_complementary
                                   ? complementary_vertex_distances
                                   : forward_vertex_distances;

      // Explore its neighbors.
      for (const auto &neighbor :
           adjacency_list_[current_vertex.graph_vertex_id]) {
        // Process neighbors in the same layaer.
        const ScoreType new_deletion_distance =
            current_vertex.distance + deletion_penalty_;

        if (vertex_distances[neighbor].find(current_vertex.query_index) !=
            vertex_distances[neighbor].end()) {
          if (new_deletion_distance <
              vertex_distances[neighbor][current_vertex.query_index]) {
            vertex_distances[neighbor][current_vertex.query_index] =
                new_deletion_distance;
            Q.push({/*graph_vertex_id=*/neighbor,
                    /*query_index=*/current_vertex.query_index,
                    /*distance=*/new_deletion_distance,
                    /*is_reverse_complementary=*/
                    current_vertex.is_reverse_complementary});
            // vertex_parent[neighbor][current_vertex.query_index] =
            //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
            //                                  ScoreType>{
            //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
            //        /*query_index=*/current_vertex.query_index,
            //        /*distance=*/current_vertex.distance};
            // std::cerr << "PUSH: " << neighbor << " "
            //          << current_vertex.query_index << " "
            //          << new_deletion_distance << std::endl;
            if (current_vertex.is_reverse_complementary) {
              ++(stats.rc_num_cells);
            } else {
              ++(stats.forward_num_cells);
            }
          }
        } else {
          vertex_distances[neighbor][current_vertex.query_index] =
              new_deletion_distance;
          Q.push({/*graph_vertex_id=*/neighbor,
                  /*query_index=*/current_vertex.query_index,
                  /*distance=*/new_deletion_distance,
                  /*is_reverse_complementary=*/
                  current_vertex.is_reverse_complementary});
          // vertex_parent[neighbor][current_vertex.query_index] =
          //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
          //                                  ScoreType>{
          //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
          //        /*query_index=*/current_vertex.query_index,
          //        /*distance=*/current_vertex.distance};

          // std::cerr << "PUSH: " << neighbor << " " <<
          // current_vertex.query_index
          //          << " " << new_deletion_distance << std::endl;
          if (current_vertex.is_reverse_complementary) {
            ++(stats.rc_num_cells);
          } else {
            ++(stats.forward_num_cells);
          }
        }

        // Process neighbors in the next layaer.
        const QueryLengthType query_index = current_vertex.query_index + 1;

        ScoreType cost = 0;
        const char vertex_label =
            current_vertex.is_reverse_complementary
                ? GetReverseComplementaryVertexLabel(neighbor)
                : labels_[neighbor];
        const char sequence_base =
            current_vertex.is_reverse_complementary
                ? sequence_bases[sequence_length - 1 - query_index]
                : sequence_bases[query_index];

        if (sequence_base != vertex_label) {
          cost = substitution_penalty_;
        }

        const ScoreType new_match_or_mismatch_distance =
            std::min(current_vertex.distance,
                     (ScoreType)(deletion_penalty_ * query_index)) +
            cost;

        // std::cerr << "seq base: " << sequence_bases[query_index] << " label:
        // " << labels_[neighbor] << " cost: " << cost << " d: " <<
        // new_match_or_mismatch_distance << std::endl;

        if (vertex_distances[neighbor].find(query_index) !=
            vertex_distances[neighbor].end()) {
          if (new_match_or_mismatch_distance <
              vertex_distances[neighbor][query_index]) {
            vertex_distances[neighbor][query_index] =
                new_match_or_mismatch_distance;

            Q.push({/*graph_vertex_id=*/neighbor,
                    /*query_index=*/query_index,
                    /*distance=*/new_match_or_mismatch_distance,
                    /*is_reverse_complementary=*/
                    current_vertex.is_reverse_complementary});
            // vertex_parent[neighbor][query_index] =
            //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
            //                                  ScoreType>{
            //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
            //        /*query_index=*/current_vertex.query_index,
            //        /*distance=*/current_vertex.distance};

            // std::cerr << "PUSH: " << neighbor << " " << query_index << " "
            //          << new_match_or_mismatch_distance << std::endl;

            if (current_vertex.is_reverse_complementary) {
              ++(stats.rc_num_cells);
            } else {
              ++(stats.forward_num_cells);
            }
          }
        } else {
          vertex_distances[neighbor][query_index] =
              new_match_or_mismatch_distance;

          Q.push({/*graph_vertex_id=*/neighbor,
                  /*query_index=*/query_index,
                  /*distance=*/new_match_or_mismatch_distance,
                  /*is_reverse_complementary=*/
                  current_vertex.is_reverse_complementary});
          // vertex_parent[neighbor][query_index] =
          //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
          //                                  ScoreType>{
          //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
          //        /*query_index=*/current_vertex.query_index,
          //        /*distance=*/current_vertex.distance};

          // std::cerr << "PUSH: " << neighbor << " " << query_index << " "
          //          << new_match_or_mismatch_distance << std::endl;
          if (current_vertex.is_reverse_complementary) {
            ++(stats.rc_num_cells);
          } else {
            ++(stats.forward_num_cells);
          }
        }
      }  // End for exploring neighbor.

      // Process insertions to the next layaer.
      const QueryLengthType query_index = current_vertex.query_index + 1;

      const ScoreType new_insertion_distance =
          current_vertex.distance + insertion_penalty_;

      if (vertex_distances[current_vertex.graph_vertex_id].find(query_index) !=
          vertex_distances[current_vertex.graph_vertex_id].end()) {
        if (new_insertion_distance <
            vertex_distances[current_vertex.graph_vertex_id][query_index]) {
          vertex_distances[current_vertex.graph_vertex_id][query_index] =
              new_insertion_distance;

          Q.push({/*graph_vertex_id=*/current_vertex.graph_vertex_id,
                  /*query_index=*/query_index,
                  /*distance=*/new_insertion_distance,
                  /*is_reverse_complementary=*/
                  current_vertex.is_reverse_complementary});
          // vertex_parent[current_vertex.graph_vertex_id][query_index] =
          //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
          //                                  ScoreType>{
          //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
          //        /*query_index=*/current_vertex.query_index,
          //        /*distance=*/current_vertex.distance};

          // std::cerr << "PUSH: " << current_vertex.graph_vertex_id << " "
          //          << query_index << " " << new_insertion_distance
          //          << std::endl;
          if (current_vertex.is_reverse_complementary) {
            ++(stats.rc_num_cells);
          } else {
            ++(stats.forward_num_cells);
          }
        }
      } else {
        vertex_distances[current_vertex.graph_vertex_id][query_index] =
            new_insertion_distance;

        Q.push({/*graph_vertex_id=*/current_vertex.graph_vertex_id,
                /*query_index=*/query_index,
                /*distance=*/new_insertion_distance,
                /*is_reverse_complementary=*/
                current_vertex.is_reverse_complementary});
        // vertex_parent[current_vertex.graph_vertex_id][query_index] =
        //    VertexWithDistanceForDijkstra<GraphSizeType, QueryLengthType,
        //                                  ScoreType>{
        //        /*graph_vertex_id=*/current_vertex.graph_vertex_id,
        //        /*query_index=*/current_vertex.query_index,
        //        /*distance=*/current_vertex.distance};

        // std::cerr << "PUSH: " << current_vertex.graph_vertex_id << " "
        //          << query_index << " " << new_insertion_distance <<
        //          std::endl;
        if (current_vertex.is_reverse_complementary) {
          ++(stats.rc_num_cells);
        } else {
          ++(stats.forward_num_cells);
        }
      }
    }

    // std::cerr << "Sequence length: " << sequence_length
    //          << ", forward alignment cost:" << forward_alignment_cost
    //          << ", reverse complement alignment cost:"
    //          << reverse_complement_alignment_cost
    //          << ", alignment cost:" << min_alignment_cost << std::endl;
    return min_alignment_cost;
  }

  ScoreType AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
      const sga::Sequence &sequence) {
    const QueryLengthType sequence_length = sequence.GetLength();
    const std::string &sequence_bases = sequence.GetSequence();
    DijkstraAlgorithmStatistics<GraphSizeType> stats;
    const ScoreType min_alignment_cost =
        AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(sequence, stats);

    std::cerr << "Sequence length: " << sequence_length
              << ", alignment cost:" << min_alignment_cost
              << ", forward num cells:" << stats.forward_num_cells
              << ", reverse num cells: " << stats.rc_num_cells << std::endl;
    return min_alignment_cost;
  }

 protected:
  char base_complement_[256] = {
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 'T', 4, 'G', 4, 4, 4, 'C', 4,   4, 4, 4,
      4, 4, 'N', 4, 4,   4, 4, 4, 'A', 4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 'T', 4, 'G', 4, 4, 4, 'C', 4, 4,   4, 4, 4, 4,   'N', 4, 4, 4,
      4, 4, 'A', 4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4,   4, 4,   4, 4, 4, 4,   4,   4, 4, 4,
      4, 4, 4,   4, 4,   4, 4, 4, 4};
  int8_t base_to_int_[256] = {
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

  // For graph representation
  std::vector<GraphSizeType> look_up_table_;
  std::vector<GraphSizeType> neighbor_table_;
  std::vector<std::vector<GraphSizeType>> adjacency_list_;
  std::vector<char> labels_;
  std::vector<std::vector<GraphSizeType>> compacted_graph_adjacency_list_;
  std::vector<std::string> compacted_graph_labels_;

  // For reverse complementary graph representation. We make the vertex ids
  // stable, and thus the vertex labels are just the corresponding reverse
  // complementary base in labels_. However, we construct a separate adjacency
  // list for fast neighbor queries in the reverse complementary graph.
  std::vector<std::vector<GraphSizeType>> reverse_complementary_adjacency_list_;

  // For RECOMB work
  std::vector<GraphSizeType> order_look_up_table_;
  std::vector<bool> visited_;
  std::vector<GraphSizeType> parents_;
  std::vector<int> types_;
  std::vector<GraphSizeType> order_offsets_;
  std::vector<GraphSizeType> order_counts_;

  // For alignment
  std::vector<std::pair<ScoreType, GraphSizeType>> distances_with_vertices_;
  ScoreType substitution_penalty_ = 1;
  ScoreType deletion_penalty_ = 1;
  ScoreType insertion_penalty_ = 1;
};

}  // namespace sga
#endif  // SGA_SEQUENCEGRAPH_H
