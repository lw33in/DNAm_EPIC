# DNA methylation EPIC arrays

- **Experiment**:
  * 500 ng genomic DNA was extracted from wt HEK293T cells and the three DKO clones using the Zymo Quick-DNA Miniprep kit (Zymo #D3024) and bisulfite-converted using the Zymo EZ DNA Methylation kit (Zymo #D5001) according to the manufacturer's protocol.
  * The bisulfite-converted DNA was analyzed using Illumina EPIC BeadArrays.
- **Computational analysis**:
  * The BeadArrays were scanned and the raw signal intensities were extracted from the *.IDAT files using the ‘noob’ function in the minfi R package.
  * The beta value (a measure of change in DNA methylation) was calculated as (M/(M+U)), in which M and U refer to the (pre-processed) mean methylated and unmethylated probe signal intensities, respectively.
  * Measurements in which the fluorescent intensity was not statistically significantly above background signal (detection P value > 0.05) were removed from the dataset.
  * Probes located from –1500 bp relative to the TSS and extending through the first coding exon (using the Illumina MethylationEPIC Manifest RefGene annotation) were included in the analysis as a defined set of ‘promoter’ probes for downstream analysis.
  * The cut off used for identifying hypomethylated or hypermethylated probes was 0.2 for the absolute beta value difference between the methylation level of a probe in the DKO cells versus the wt HEK293T cells.
 
Reference: W. Ni, A. A. Perez, S. Schreiner, C. M. Nicolet, and P. J. Farnham. Characterization of the ZFX family of transcription factors that bind downstream of the start site of CpG island promoters. Nucleic Acids Res, 2020. doi: 10.1093/nar/gkaa384.

  
