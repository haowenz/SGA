#include <benchmark/benchmark.h>

#include "sequence_batch.h"
#include "sequence_graph.h"

const static std::string txt_sequence_graph_file_path = "BRCA1_seq_graph.txt";
const static std::string sequence_file_path = "BRCA1_5_reads.fastq";
const static uint32_t max_batch_size = 1000000;
static sga::SequenceGraph<> sequence_graph;
static sga::SequenceBatch sequence_batch(max_batch_size);

static void DoSetup(const benchmark::State& state) {
  static bool is_setup_called = false;
  if (!is_setup_called) {
    std::cerr << "Setup benchmark!" << std::endl;
    sequence_graph.SetAlignmentParameters(2, 3, 3);
    sequence_graph.LoadFromTxtFile(txt_sequence_graph_file_path);
    sequence_graph.GenerateCharLabeledGraph();
    sequence_graph.GenerateCompressedRepresentation();
    sequence_batch.InitializeLoading(sequence_file_path);
    sequence_batch.LoadBatch();
    is_setup_called = true;
  }
}

static void DoTeardown(const benchmark::State& state) {
  static bool is_tear_down_called = false;
  if (!is_tear_down_called) {
    std::cerr << "Tear down benchmark!" << std::endl;
    sequence_batch.FinalizeLoading();
    is_tear_down_called = true;
  }
}

static void BM_AlignUsingLinearGapPenalty(benchmark::State& state) {
  const uint32_t num_sequences = sequence_batch.GetNumLoadedSequences();

  for (auto _ : state) {
    for (uint32_t si = 0; si < num_sequences; ++si) {
      sequence_graph.AlignUsingLinearGapPenalty(sequence_batch.GetSequence(si));
    }
  }
}

static void BM_AlignUsingLinearGapPenaltyWithNavarroAlgorithm(
    benchmark::State& state) {
  const uint32_t num_sequences = sequence_batch.GetNumLoadedSequences();
  for (auto _ : state) {
    for (uint32_t si = 0; si < num_sequences; ++si) {
      sequence_graph.AlignUsingLinearGapPenaltyWithNavarroAlgorithm(
          sequence_batch.GetSequence(si));
    }
  }
}

static void BM_AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
    benchmark::State& state) {
  const uint32_t num_sequences = sequence_batch.GetNumLoadedSequences();
  for (auto _ : state) {
    for (uint32_t si = 0; si < num_sequences; ++si) {
      sequence_graph.AlignUsingLinearGapPenaltyWithDijkstraAlgorithm(
          sequence_batch.GetSequence(si));
    }
  }
}

BENCHMARK(BM_AlignUsingLinearGapPenalty)->Setup(DoSetup)->Teardown(DoTeardown);
BENCHMARK(BM_AlignUsingLinearGapPenaltyWithNavarroAlgorithm)
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);
BENCHMARK(BM_AlignUsingLinearGapPenaltyWithDijkstraAlgorithm)
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

BENCHMARK_MAIN();
