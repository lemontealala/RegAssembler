# RegAssembler
RegAssembler is a genome assembling tool for long reads. Sequencing reads from Pacbio and nanopore, such as raw CLS reads, hifi reads, corrected reads preprocessed by correction tools can be assembled by RegAssembler.

The workflow of RegAssembler contains four major modules. 
1. Self-mapping and preprocessing module. Long reads are self-mapped and preprocessed. Reads containing low coverage breakpoint are identified as suspicious chimeras and thus removed from self-mapping results. Overlaps with unreasonable hanging outs are also removed.
2. The assembly module. An overlap graph is constructed based on detected overlaps for long reads. Then RegAssembler orientates and positions long reads in a linear regression framework to generate contigs. A weight-decreasing breadth-first-searching (WBFS) algorithm is applied to obtain connected components whose maximum size are restricted. Furthermore, the WBFS determines the assignment of direction in order to minimize overlaps with conflicting relative orientations within each component. As for estimating the start coordinates of long reads, a two-step robust regression estimation is applied. The estimated orientations and positions generate linear arrangement of long reads, that is, contigs.
3. The splitting and iterating module. Though robust regression is applied, a few unrecognized false overlaps may bond together multiple strains from different genomic regions. RegAssembler then split regression results to obtain conflict-free linear contigs by verifing eligible overlaps between adjacent reads. By replacing the raw reads with the resulted contigs (pre-contigs), RegAssembler iterates the assembly module on contig level to obtain longer contigs. In the iteration, RegAssembler inherits the overlaps between long reads positioned at only ends of the contigs as overlaps between contigs. Thus a bigger contig graph is constructed and longer contigs are assembled by orientating and positioning pre-contigs. Consequently, long chromosomes can be reconstructed by iteration.
4. The consensus module. RegAssembler extract final consensus by a block POA algorithm. Each contig is segmented into several blocks with a fixed bandwidth. Reads in each block are multi-aligned by the Patitial order alignment(POA) algorithm. The POA algorithm is paralized implemented so that the consensus module is as effective as accuracy.

We have tested RegAssembler on both simulated and real data. A comparsion with canu and wtdgb2 on C .elegans chromosome 1 is exhibited in the following table. The quality index are from QUAST 5.2.0.

| Metric                      | RegAssembler | Canu | wtdgb2 |
|-----------------------------|------------------------|------------------------|------------------------|
| \# contigs                   | 1                      | 1                      | 1                      |
| Largest contig               | 15099728               | 15066353               | 12601506               |
| Total length                 | 15099728               | 15066353               | 15030370               |
| Reference length             | 15072434               | 15072434               | 15072434               |
| GC (%)                       | 35.72                  | 35.75                  | 35.76                  |
| Reference GC (%)             | 35.75                  | 35.75                  | 35.75                  |
| N50                          | 15099728               | 15066353               | 12601506               |
| NG50                         | 15099728               | 15066353               | 12601506               |
| N90                          | 15099728               | 15066353               | -                      |
| NG90                         | 15099728               | 15066353               | -                      |
| L50                          | 1                      | 1                      | 1                      |
| LG50                         | 1                      | 1                      | 1                      |
| L90                          | 1                      | 1                      | -                      |
| LG90                         | 1                      | 1                      | -                      |
| \# misassemblies             | 0                      | 0                      | 4                      |
| \# misassembled contigs      | 0                      | 0                      | 1                      |
| Misassembled contigs length  | 0                      | 0                      | 12601506               |
| \# local misassemblies       | 0                      | 0                      | 1                      |
| \# scaffold gap ext. mis.    | 0                      | 0                      | 0                      |
| \# scaffold gap loc. mis.    | 0                      | 0                      | 0                      |
| \# unaligned mis. contigs    | 0                      | 0                      | 0                      |
| \# unaligned contigs         | 0 + 0 part            | 0 + 0 part            | 1 + 0 part            |
| Unaligned length             | 0                      | 0                      | 4990                   |
| Genome fraction (%)          | 99.976                 | 99.961                 | 99.677                 |
| Duplication ratio            | 1.002                  | 1.000                  | 1.000                  |
| Largest alignment            | 15099728               | 15066353               | 5446807                |
| Total aligned length         | 15099728               | 15066353               | 15022949               |
| NA50                         | 15099728               | 15066353               | 4284491                |
| NGA50                        | 15099728               | 15066353               | 4284491                |
| NA90                         | 15099728               | 15066353               | -                      |
| NGA90                        | 15099728               | 15066353               | -                      |
| auNA                         | 15099728.0             | 15066353.0             | 15030370.0             |
| auNGA                        | 15127071.4             | 15060274.5             | 15002419.0             |
| LA50                         | 1                      | 1                      | 2                      |
| LGA50                        | 1                      | 1                      | 2                      |
| LA90                         | 1                      | 1                      | -                      |
| LGA90                        | 1                      | 1                      | -                      |

## Usage
The following RegAssembler.sh provides a example for runing RegAssembler, where '-r' specifies the long reads input, '-m' specifies the self-mapping results from minimap2. Usages for other arguments can be known by '-h' command.
`python RegAssembler.py -r test_data/Celegans_chr1_simul_Q10_0001.fasta -m test_data/minimap2_Q10_simu.paf -p '_' -sd 1500 -bw 3000 -HO 150 -t 40 -C -o Celegans.fa > RegAssembler.log`


