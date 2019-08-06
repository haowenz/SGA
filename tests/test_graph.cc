#include "gtest/gtest.h"
#include "sequence_batch.h"
#include "sequence_graph.h"

namespace sga_testing {

class SequenceGraphTest : public ::testing::Test {
 protected:
  SequenceGraphTest(): sequence_graph_(), sequence_batch_(5) {
  }

  ~SequenceGraphTest() override {
  }

  void SetUp() override {
    sequence_batch_.InitializeLoading(sequence_file_path_);

    sequence_graph_.LoadFromVGFile(sequence_graph_file_path_);
    sequence_graph_.GenerateCharLabeledGraph();
    sequence_graph_.GenerateCompressedRepresentation();
  }

  void TearDown() override {
    sequence_batch_.FinalizeLoading();
  }

  sga::SequenceGraph<> sequence_graph_;
  const std::string sequence_graph_file_path_ = "BRCA1_seq_graph.vg";

  sga::SequenceBatch sequence_batch_;
  const std::string sequence_file_path_ = "BRCA1_5_reads.fastq";

};

TEST_F(SequenceGraphTest, LoadFromVGFileTest) {
  uint32_t num_vertices = sequence_graph_.GetNumVerticesInCompactedGraph();
  ASSERT_EQ(num_vertices, (uint32_t)83) << "Number of vertices loaded is wrong! It should be 83 but it is " << num_vertices;
  uint32_t num_edges = sequence_graph_.GetNumEdgesInCompactedGraph();
  ASSERT_EQ(num_edges, (uint32_t)81) << "Number of edges loaded is wrong! It should be 81 but it is " << num_edges;
}

TEST_F(SequenceGraphTest, GenerateCharLabeledGraphTest) {
  uint32_t num_vertices = sequence_graph_.GetNumVertices();
  ASSERT_EQ(num_vertices, (uint32_t)81190) << "Number of vertices loaded is wrong! It should be 81190 but it is " << num_vertices;
  uint32_t num_edges = sequence_graph_.GetNumEdges();
  ASSERT_EQ(num_edges, (uint32_t)(81188)) << "Number of edges loaded is wrong! It should be (81188) but it is " << num_edges;
}

TEST_F(SequenceGraphTest, AlignUsingLinearGapPenaltyTest) {
  uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  int32_t max_alignment_scores[5] = {62, 25, 54, 9, 37};
  sequence_graph_.SetAlignmentParameters(1, 1, 1);
  for (uint32_t i = 0; i < num_loaded_sequences; ++i) {
    int32_t alignment_score = sequence_graph_.AlignUsingLinearGapPenalty(sequence_batch_.GetSequence(i));
    ASSERT_EQ(alignment_score, max_alignment_scores[i]) << "Alignment score for sequence" << i << " is wrong! It should be " << max_alignment_scores[i] << " but it is " << alignment_score;
  }
}

}  // namespace sga_testing

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
