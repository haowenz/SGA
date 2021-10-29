#include "gtest/gtest.h"
#include "sequence.h"
#include "sequence_batch.h"

namespace sga_testing {

class SequenceBatchTest : public ::testing::Test {
 protected:
  SequenceBatchTest() : sequence_batch_(2) {
    sequence_batch_.InitializeLoading(sequence_file_path_);
  }

  ~SequenceBatchTest() override { sequence_batch_.FinalizeLoading(); }

  void SetUp() override {}

  void TearDown() override {}

  sga::SequenceBatch sequence_batch_;
  const std::string sequence_file_path_ = "BRCA1_5_reads.fastq";
};

TEST_F(SequenceBatchTest, LoadBatchTest) {
  uint32_t num_loaded_sequences = sequence_batch_.LoadBatch();
  ASSERT_EQ(num_loaded_sequences, (uint32_t)2)
      << "Number of sequences loaded is wrong! It should be 2 but it is "
      << num_loaded_sequences;
  num_loaded_sequences = sequence_batch_.LoadBatch();
  ASSERT_EQ(num_loaded_sequences, (uint32_t)2)
      << "Number of sequences loaded is wrong! It should be 2 but it is "
      << num_loaded_sequences;
  num_loaded_sequences = sequence_batch_.LoadBatch();
  ASSERT_EQ(num_loaded_sequences, (uint32_t)1)
      << "Number of sequences loaded is wrong! It should be 2 but it is "
      << num_loaded_sequences;
  sga::Sequence true_sequence;
  true_sequence.Update(
      325, "S1_5",
      "TCTTGTAATTTAATTTCGATTACTAATTTCCTGAAAATTATAGAACTAGATAAAGCTATATAGTTGTGGATT"
      "ATTTATGGTATAGTTTACTTGAGAAATAATTATTAAATATTAGTGGAAAAAGCTATACTTTGGGTATGATAA"
      "GGAACTTCCTCAATTGAATTTCCTTTCCTATCTGTAAAAAGCAAGTAGGGTTAAATAGTTTTATTCCGCCAG"
      "AAGGCATCTTTTATTCTTCTCCCCTTGTCCTCACATGGGTGAATTTACCAGCACATTTAACTAAATTCAGAA"
      "ACTGGTTCCAAATGTACTGCAGATAGTAGGCAATTTC",
      ")**)**)********)*******)**********)*)****************)*******)***)******"
      "*****************)*)***)*)***)*****)****************************)*******"
      "************)****)*************)******)*))*)*******)*)***))***)****)****"
      "*)**)*******))*********)**))******)***))****)**************************)"
      "************************)*****)****)*");
  sga::Sequence loaded_sequence = sequence_batch_.GetSequence(0);
  ASSERT_EQ(loaded_sequence.GetLength(), true_sequence.GetLength())
      << "The length of the last loaded sequence is wrong. It should be 325 "
         "but it is "
      << loaded_sequence.GetLength();
  ASSERT_STREQ(loaded_sequence.GetName().c_str(),
               true_sequence.GetName().c_str());
  ASSERT_STREQ(loaded_sequence.GetSequence().c_str(),
               true_sequence.GetSequence().c_str());
}

}  // namespace sga_testing

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
