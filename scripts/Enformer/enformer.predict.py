#Imports
import tensorflow as tf
#Make sure the GPU is enabled  (use tensorflow -- GPU version)
assert tf.config.list_physical_devices('GPU'), 'Start the colab kernel with GPU: Runtime -> Change runtime type -> GPU'

import tensorflow_hub as hub

import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import sys
import argparse

#Code
# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:

  def __init__(self, tfhub_dir):
    self._model = hub.load(tfhub_dir).model
  #def __init__(self, tfhub_url):
  #  self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]

    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)
    return tf.reduce_sum(input_grad, axis=-1)

#Code
# @title `variant_centered_sequences`

class FastaStringExtractor:
    
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()
    
def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

#Variant scoring

#variant's input format
def variant_to_dict(snp_file):
    variant_dict={}
    for line in open(snp_file):
        line=line.strip()
        snpinfo=line.split('_')
        snpinfo.append(line)
        variant_dict[line] = snpinfo
        
    return variant_dict


#prediction
def enformer_predict(genome_fa,snp_file):
    #model
    model = Enformer("./moduleE")
   
    #reference genome sequences
    fasta_extractor = FastaStringExtractor(genome_fa)
    
    #SNP
    variant_dict=variant_to_dict(snp_file)
    print('SNPs were loaded',flush=True)
    
    df={}
    for pos_name in variant_dict.keys():
        print(pos_name,flush=True)
        snpinfo = variant_dict[pos_name]
        variant = kipoiseq.Variant(snpinfo[0], snpinfo[1], snpinfo[2], snpinfo[3], snpinfo[4])  # @param
        
        # Center the interval at the variant
        interval = kipoiseq.Interval(variant.chrom, variant.start, variant.start).resize(SEQUENCE_LENGTH)
        seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
        center = interval.center() - interval.start

        reference = seq_extractor.extract(interval, [], anchor=center)
        alternate = seq_extractor.extract(interval, [variant], anchor=center)

        # Make predictions for the refernece and alternate allele
        reference_prediction = model.predict_on_batch(one_hot_encode(reference)[np.newaxis])['human']#shape:1*896*5313
        alternate_prediction = model.predict_on_batch(one_hot_encode(alternate)[np.newaxis])['human']#shape:1*896*5313
        
        delta = alternate_prediction.mean(axis=1) - reference_prediction.mean(axis=1) #1*5313
        delta = delta.squeeze()
        df[pos_name]= delta
        
    diff=pd.DataFrame.from_dict(df,orient='index')
    diff=diff.reset_index().rename(columns={'index':'snp'})
    return diff


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-s", "--snp", dest="snp", required=True, help="snp information")
    parser.add_argument("-g", "--genome", dest="genome", required=True, help="genome sequence")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output file/directory (bindind difference)")
    args = parser.parse_args()
    
    #snp information
    snpfile = args.snp

    #genome sequences
    genomefa = args.genome+".fa"
    
    #output
    outfile = args.out

    #predict
    batch_result = enformer_predict(genome_fa=genomefa,snp_file=snpfile)

    #save results
    batch_result.to_csv(outfile,index=False,header=False,sep='\t')


if __name__ == "__main__":
    sys.exit(main())

