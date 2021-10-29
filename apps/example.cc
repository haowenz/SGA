#include <string>

#include "sequence_batch.h"
#include "sequence_graph.h"

int main(int argc, char *argv[]) {
  std::string sequence_graph_file_path;
  std::string sequence_file_path;
  if (argc != 3) {
    std::cerr << "Usage:\t" << argv[0] << "\tgraph_file\tread_file\n";
    exit(-1);
  } else {
    sequence_graph_file_path = argv[1];
    sequence_file_path = argv[2];
  }
  uint32_t max_batch_size = 1000000;
  sga::SequenceGraph<> sequence_graph;
  sequence_graph.SetAlignmentParameters(2, 3, 3);
  sga::SequenceBatch sequence_batch(max_batch_size);
  sequence_graph.LoadFromTxtFile(sequence_graph_file_path);
  sequence_graph.GenerateCharLabeledGraph();
  sequence_graph.GenerateCompressedRepresentation();
  sequence_batch.InitializeLoading(sequence_file_path);
  uint32_t num_sequences = sequence_batch.LoadBatch();
  while (num_sequences > 0) {
    for (uint32_t si = 0; si < num_sequences; ++si) {
      sequence_graph.AlignUsingLinearGapPenaltyWithNavarroAlgorithm(
          sequence_batch.GetSequence(si));
    }
    num_sequences = sequence_batch.LoadBatch();
  }
  sequence_batch.FinalizeLoading();
}
