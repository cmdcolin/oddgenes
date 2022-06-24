# oddgenes

A list of weird gene annotations or things that break bioinformatics assumptions

## Gene structures

### 1bp length exon

Evidence given for a 1bp length exon in Arabadopsis and different splicing models are discussed

http://www.nature.com/articles/srep18087

Another 1bp exon is discussed here https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177959

### 0bp length exon

The phenomenon of recursive splicing can remove sequences progressively inside an intron, so there can exist "0bp exons" that are just the splice-site sequences pasted together.

"To identify potential zero nucleotide exon-type ratchet points, we parsed the RNA-Seq alignments to identify novel splice junctions where the reads mapped to an annotated 5' splice site and an unannotated 3' splice site, and the genomic sequence at the 3' splice site junction was AG/GT"

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4529404/

### Twintron

A twintron is essentially an intron-within-an-intron, which could be formed by a mobile element (TE) insertion. The original idea is that the internal intron has to be spliced first before the outer one is, but several classes have been discovered. See https://en.wikipedia.org/wiki/Twintron

### Very large introns

Satellite DNA study uncovers megabase scale introns
https://www.biorxiv.org/content/early/2018/12/11/493254

An example in this paper kl-3 spans 4.3 million bp

In human, an example is Dystrophin spanning 2.3 million bp

### Backsplicing and circRNAs

The process of "backsplicing" circularizes RNAs. There can be alternative backsplicing too

See https://academic.oup.com/nar/article/48/4/1779/5715065

### Very large number of isoforms in Dscam

"Dscam has 24 exons; exon 4 has 12 variants, exon 6 has 48 variants, exon 9 has 33 variants, and exon 17 has two variants. The combination of exons 4, 6, and 9 leads to 19,008 possible isoforms with different extracellular domains (due to differences in Ig2, Ig3 and Ig4). With two different transmembrane domains from exon 17, the total possible protein products could reach 38,016 isoforms"

Ref https://en.wikipedia.org/wiki/DSCAM https://www.wikigenes.org/e/gene/e/35652.html

### Translational frameshift

"The main distinction between frameshifts resulting from mutation and those resulting from ribosomal frameshifting is that the latter are controlled by various mechanisms found in codons...Certain codons take longer to translate, because there are not equal amounts of tRNA of that particular codon in the cytosol..." which leads to ribosomal slippage into an alternative reading frame.  

Ref https://en.wikipedia.org/wiki/Translational_frameshift

https://www.sciencedirect.com/topics/neuroscience/ribosomal-frameshifting

### Ribosome hopping
"Ribosome hopping involves ribosomes skipping over large portions of an mRNA without translating them"
Ref https://pubmed.ncbi.nlm.nih.gov/24711422/

### Internal Ribosome Entry Sites (IRES)
"Eukaryotic mRNAs are typically monocistronic and translated only a single Open Reading Frame. Some viruses can reinititate translation after translation termination using an IRES"
Ref https://en.wikipedia.org/wiki/Internal_ribosome_entry_site

### A Stop codon that is not a stop codon

In some cases a stop codon is not interpreted as such. When it is interpreted, it is sometimes called "Stop codon readthrough" and can encode for an amino acid. The amino acid Selenocysteine is coded for by a stop codon (https://en.wikipedia.org/wiki/Selenocysteine) and Pyrrolysine also is coded for by a stop codon (https://en.wikipedia.org/wiki/Pyrrolysine). Both of these lie outside the conventional 20 amino acid code

There are several other stop codon modifications described here https://www.nature.com/articles/nrg3963

Selenocysteine can be coded via a SECIS sequence https://en.wikipedia.org/wiki/SECIS_element and resulting products are called selenoproteins

Pyrolysine is coded through a pyIT tRNA gene that interprets the amber stop codon as pyrolysine

### Readthrough transcription

See also this Ensembl blog on annotating readthrough transcription which joins multiple genes http://www.ensembl.info/2019/02/11/annotating-readthrough-transcription-in-ensembl/

RNA-seq often makes extremely compelling cases for two-or-more different genes to be conjoined by splicing

Some algorithms e.g. mikado https://academic.oup.com/gigascience/article/7/8/giy093/5057872 try to avoid this calling it artifactual fusion/chimera that can be due to some tandem duplication but it does seem to be very prevalent in real data sets

### Non-canonical splice sites

The standard splice site recognition sequence is an GU in RNA (or GT in DNA) on the 5' end and AG on the 3' (remember, goes 5' to 3'). This recognition motif accounts for the large majority of splicing. If a different sequence is used it is said that a different spliceosome complex is being used "minor spliceosome"

https://en.wikipedia.org/wiki/Minor_spliceosome

### Cryptic splice sites

Some exons harbor internal splice sites (e.g. they get split) that might be unused or underused and are so called "cryptic splice sites"

Review article https://academic.oup.com/nar/article/39/14/5837/1382796

The snaptron project from Ben Langmead analyzed huge amounts of RNA-seq public data and found many types of these cryptic splicing http://snaptron.cs.jhu.edu/

### Wobble splicing

NAGNAG, GYNGYN, repeats of the splicing signal cause modified transcriptional behavior

"Another mechanism introducing small variations to protein isoforms is wobble splicing. Here, a GYN repeat at the donor splice site (5’ splice site; Y stands for C or T and N stands for A, C, G, or T) or an NAG repeat at the acceptor splice site (3’ splice site) leads to subtle length variations in the spliced transcripts and finally to alternative isoforms differing in few amino acids." ref https://onlinelibrary.wiley.com/doi/full/10.1002/bies.201900066?af=R

### Intron retention

Intron retention (IR) is a phenomenon where intron sequence is preserved, or doesn't get spliced out, in mature RNA

It can occur in both abnormal and normal biological conditions. Transcript with IR often undergo nonsense-mediated decay.

### Self-splicing RNA

Normally RNA is spliced by a specialized protein complex called a spliceosome. There is also self-splicing RNA where the splicing is done itself with RNA

The Group 1 intron type mentioned above is a "self splicing" function of RNA not requiring external spliceosome https://en.wikipedia.org/wiki/Group_I_catalytic_intron

Group 2 and group 3 with similar but different mechanisms also exist

### Introns in archaea

The only types of introns known conventionally in archaea are called "bulge-helix-bulge" but recently Group 1 introns have been discovered https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky414/4999243

### Codon tables

Many eukaryotes use the "standard genetic code" for changing codons to amino acids but frequent changes occur across the domains of life. The NCBI "genetic code" table lists several of these and contains recent additions for particular species

https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG31

One article explains how alternative genetic codes work https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6207430/

### Untranslated regions

The 5' and 3' UTR (un-translated region) is a part of the pre-mRNA at the start and end of the gene respectively that is spliced away in the mature RNA

This blog post by Ensembl shows how they annotate UTR and 19kb 3' UTR in Grin2b http://www.ensembl.info/2018/08/17/ensembl-insights-how-are-utrs-annotated/

They have many important functionality and are often targets of miRNA binding which leads to degradation.

### Poly-A tails

A poly-A tail is added to the pre-mRNA on the 3' end of the transcript to protect it from degradation. The A signal is not part of the genome https://en.wikipedia.org/wiki/Polyadenylation

A survey of poly-A using Oxford Nanopore found some transcript isoforms with 450bp ENST00000581230, with intron retention being a possible correlate of having a longer poly-A tails https://www.biorxiv.org/content/early/2018/11/09/459529.article-info

Intronic polyadenylation can also occur https://www.nature.com/articles/s41467-018-04112-z it is revealed by 3'-seq

### Circular chromosomes

Circularized chromosomes should be unsurprising to anyone working with plasmids and many prokaryotic genomes but for gene annotation formats which use linear coordinates, representing anything wrapping around the origin is challenging. 

Many genomic viewers do not do this well. For GFF format this is done by making the end go past the end of the genome. Below, the genome is 6407 bp in length, but the CDS feature extends past this and sets Is_circular=true

```
##gff-version 3.2.1
# organism Enterobacteria phage f1
# Note Bacteriophage f1, complete genome.
J02448  GenBank region  1      6407    .       +       .       ID=J02448;Name=J02448;Is_circular=true;
J02448  GenBank CDS     6006   7238    .       +       0       ID=geneII;Name=II;Note=protein II;
```
### Dynamic DNA structures in vivo

The replication of the 2 micron plasmid found in Saccharomyces cerevisiae relies on a programmed DNA rearrangement; in any population of cells two different states of the 2 micron plasmid can be expected and these will interconvert in later generations.
Reference: https://pubmed.ncbi.nlm.nih.gov/23541845/

### Overlapping genes

It is possible for gene sequences to overlap possibly in alternate coding frames

https://en.wikipedia.org/wiki/Overlapping_gene

Some articles 

- The novel EHEC gene asa overlaps the TEGT transporter gene in antisense and is regulated by NaCl and growth phase https://www.ncbi.nlm.nih.gov/m/pubmed/30552341/
- Overlapping genes in natural and engineered genomes https://www.nature.com/articles/s41576-021-00417-w
- Uncovering de novo gene birth in yeast using deep transcriptomics https://www.nature.com/articles/s41467-021-20911-3


## Flybase

### Chimeric genes

The gene Jingwei is a chimera of two genes, alcohol dehydrogenage and yellow emperor. Many chimeras are damaging but this has been selected for

http://www.pnas.org/content/101/46/16246

## Wormbase

### Adding leader sequence to mRNA

"About 70% of C. elegans mRNAs are trans-spliced to one of two 22 nucleotide spliced leaders. SL1 is used to trim off the 5' ends of pre-mRNAs and replace them with the SL1 sequence. This processing event is very closely related to cis-splicing, or intron removal."

The region that is spliced out is called an outron

http://www.wormbook.org/chapters/www_transsplicingoperons/transsplicingoperons.html

### Polycistronic transcripts/operons

Although prevalent in bacteria, operons are not common in eukaryotes. However, they are common in C. elegans specifically. "A characteristic feature of the worm genome is the existence of genes organized into operons. These polycistronic gene clusters contain two or more closely spaced genes, which are oriented in a head to tail direction. They are transcribed as a single polycistronic mRNA and separated into individual mRNAs by the process of trans-splicing"

http://www.wormbook.org/chapters/www_overviewgenestructure.2/genestructure.html

### Trans-splicing of exons on different strands

A pre-mRNA from both strands of DNA eri6 and eri7 are combined to create eri-6/7

Source http://forums.wormbase.org/index.php?topic=1225.0 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2756026/

### Exon shared across different genes

An example from drosophila, C. elegans, and rat shows a gene with a 5' exon being shared between two genes

![](twosplice.png)

Source http://forums.wormbase.org/index.php?topic=1225.0 https://www.fasebj.org/doi/full/10.1096/fj.00-0313rev

An example here shows 5'UTR exons shared across different olfactory receptor genes ("Some OR genes share 5'UTR exons")

https://www.biorxiv.org/content/biorxiv/early/2019/09/19/774612.full.pdf

## Evolution

### Possible adaptive bacteria->eukaryote HGT

A possible horizontal gene transfer from bacteria to eukaryotes is found in an insect that feeds on coffee beans. Changes that the gene had to undergo are covered (added poly-A tail, shine-dalgarno sequence deleted)

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3306691/

also https://www.cell.com/cell/fulltext/S0092-8674(19)30097-2

### Transgenerational epigenetic inheritence

This phenomena of epigenetic modifications being passed down across generations garners a lot of media attention and scientific attention. The idea of it being influenced by what "one does in life" such as experiencing famine is also very interesting.

https://en.wikipedia.org/wiki/Transgenerational_epigenetic_inheritance

There are skeptics also http://www.wiringthebrain.com/2018/07/calibrating-scientific-skepticism-wider.html but the science is hopefully what speaks for itself

## Codon usage

### Alternative start codons

"The most common start codons for known Escherichia coli genes are AUG (83% of genes), GUG (14%) and UUG (3%)"

"Here, we systematically quantified translation initiation of green fluorescent protein (GFP) from all 64 codons and nanoluciferase from 12 codons on plasmids designed to interrogate a range of translation initiation conditions."

https://www.sciencedaily.com/releases/2017/02/170221080506.htm

Testing in eukaryotes has also revealed alternative starts being viable https://en.wikipedia.org/wiki/Start_codon#Eukaryotes

## Molecular

### Complex DNA structures

The standard DNA double stranded helix is called B-DNA (https://genome.cshlp.org/content/early/2018/11/06/gr.241257.118.abstract)

Other geometries are possible https://en.wikipedia.org/wiki/Nucleic_acid_double_helix#Helix_geometries

Complex structures such as four stranded quadruplex have been found that could have biological functions

See https://news.cnrs.fr/articles/unlocking-the-secrets-of-four-strand-dna

### Polytene chromosome

Some organisms, famously insect salivary glands, create many copies of genes through multiple phases of incomplete DNA replication https://en.wikipedia.org/wiki/Polytene_chromosome

![](polytene.png)

Figure source https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5768140/

"Polytene chromosomes are produced by endoreplication, in which chromosomal DNA undergoes mitotic replication, but the strands do not separate. Ten rounds of endoreplication produces 2^10 = 1,024 DNA strands, which when arranged alongside of each other produce distinctive banding patterns. Endoreplication occurs in cells of the larval salivary glands of many species of Diptera, and increases production of mRNA for Glue Protein that the larvae use to anchor themselves to the walls of (for example) culture vials."
from https://www.mun.ca/biology/scarr/Polytene_Chromosomes.html

## RNA world

### RNA modifications

https://www.hindawi.com/journals/jna/2011/408053/tab1/

updated link on hindawi should point here http://mods.rna.albany.edu/mods/

### RNA editing

RNA editing is a post-transcriptional modification to the mRNA to change the bases. A-to-I editing is a common one in mammals which would make the RNA, when sequenced, to have a G instead of an A. So WGS would not show a SNP but the RNA-seq would appear to have A->G.

Other editing occurs also https://en.wikipedia.org/wiki/RNA_editing

### Maternal RNAs being passed down

Maternal RNAs being active against the zygote (e.g. https://en.wikipedia.org/wiki/Maternal_to_zygotic_transition) and lead to complex transgenerational effects

### Lowly expressed RNA has large effects

A lncRNA VELUCT almost flies under the radar in a lung cancer screen due to being very lowly expressed such that it is "below the detection limit in total RNA from NCI-H460 cells by RT-qPCR as well as RNA-Seq", however this study confirms it as a factor in experiments

https://www.ncbi.nlm.nih.gov/pubmed/28160600?dopt=Abstract

Note that X inactivation relies on relatively lowly expressed RNA also https://twitter.com/mitchguttman/status/1454256452990734336

### X chromosome inactivation

X chromosome inactivation is produced by a non-coding transcript called Xist is transcribed on the X that is being inactivated and actually coats the X chromosome with itself. An anti-sense transcript called Tsix regulates Xist

https://en.wikipedia.org/wiki/XIST

https://en.wikipedia.org/wiki/X-inactivation#Xist_and_Tsix_RNAs

### Types of RNA

There are many types of RNA some more weird an exotic than others, a large list https://en.wikipedia.org/wiki/List_of_RNAs

Some are named based on where they are expressed or active

Others are uniquely shaped. There are also circular RNA for example https://en.wikipedia.org/wiki/Circular_RNA

Small and long non coding RNAs often fold into important structural shapes

## Proteins

### Removal of start amino acid in proteins

This is probably obvious to many people who work on proteins but while the genome has almost all genes starting with a start codon which produces methionine, this is often post translationally removed https://en.m.wikipedia.org/wiki/Methionyl_aminopeptidase

### Inteins

An intein is like an intron but for a protein, a segment of protein that is spliced out https://en.wikipedia.org/wiki/Intein

See section here https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#pathological-cases

### Polyprotein

Viral sequences can create a polyprotein which is fully transcribed and translated before being cleaved by a protease. In some viruses (such as coronaviruses) their translation involves ribosomal frameshifting.

Dengue, HIV, flu, etc. use this

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6040172/
https://www.sciencedirect.com/science/article/abs/pii/S0959440X15000597

## Transposons

### Cross-species BovB transposon transfers

Or "How a quarter of the cow genome came from snakes" http://phenomena.nationalgeographic.com/2013/01/01/how-a-quarter-of-the-cow-genome-came-from-snakes/

Source http://www.pnas.org/content/110/3/1012.full

### LINE1 important for embryonic development

Transposon activity can mutate DNA as it will insert itself into the genome. The genome has functions for keeping transposons inactive. However, evidence shows that the LINE1 is important for embryonic development.

https://www.ucsf.edu/news/2018/06/410781/not-junk-jumping-gene-critical-early-embryo

## Immunity

### VDJ Recombination

VDJ recombination is a process of somatic recombination (using "recombination signal sequences") that is done in immune cells. Different gene segments of class "V", class "D", and class "J" exons (sometimes to exons are referred to as "genes" themselves) are somatically rearranged into coherent genes that are then transcribed and immune diversity. Splicing at the DNA level is not precise, with terminal transferase adding random nucleotides to further diversify the sequences

https://en.wikipedia.org/wiki/V(D)J_recombination


### MHC region

The MHC region is a very polymorphic region of the genome on chr6. I'm not personally aware of the intricacies of MHC beyond that it is a unique contributor of some additional hg38 alternative loci/contigs

- https://en.wikipedia.org/wiki/Major_histocompatibility_complex

- https://en.wikipedia.org/wiki/Human_leukocyte_antigen

## Structural variations

### Tandem duplication

What is a tandem duplication? Why does it occur?

Factors can include

- replication slippage
- retrotransposition
- unequal crossing over (UCO).
- imperfect repair of double-strand breaks by nonhomologous end joining (NHEJ).

Ref https://academic.oup.com/mbe/article/24/5/1190/1038942

## Pseudogenes

### A pseudogene that can protect against cancer in Elephants

The LIF gene has many copies in Elephant but many are non-functional. One copy can be "turned back on" and play a role in cancer protection. They call this a "zombie gene"

https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31145-8

https://www.sciencealert.com/lif6-pseudogene-elephant-tumour-suppression-solution-petos-paradox

## Regulation

### Intron mediated enhancement (IME)

It has been shown that some intron sequences can enhance expression similar to how promoter sequences work https://en.wikipedia.org/wiki/Intron-mediated_enhancement

The first intron of the UBQ10 gene in Arabidopsis exhibits IME, and "the sequences responsible for increasing mRNA accumulation are redundant and dispersed throughout the UBQ10 intron" http://www.plantcell.org/content/early/2017/04/03/tpc.17.00020.full.pdf+html

The classic peppered moth phenotype is a intron TE insertion https://wp.unil.ch/genomeeee/2016/12/16/peppered-moth-melanism-mutation-is-a-transposable-element/

### Bidirectional promoters

Wikipedia https://en.wikipedia.org/wiki/Promoter_(genetics)#Bidirectional_(mammalian)

"Bidirectional promoters are a common feature of mammalian genomes. About 11% of human genes are bidirectionally paired."

"The two genes are often functionally related, and modification of their shared promoter region allows them to be co-regulated and thus co-expressed"

## Chromosomal abnormalities

### Loss of Y chromosome

Older men can have a mosaic loss of the Y chromosome in blood samples

https://www.karger.com/Article/FullText/508564 (found from https://www.biostars.org/p/9482437/)

## File formats

### Non-ACGT letters in fasta files

The latest human genome, for example, downloaded from NCBI, contains a number of Non-ACGT letters in the form of IUPAC codes https://www.bioinformatics.org/sms/iupac.html These represent ambiguous bases.

Here is the incidence of non-ACGT IUPAC letters in the entire human genome GRCh38.p10 NC_000001-24

```
b: 2
k: 8
m: 8
r: 26
s: 5
w: 14
y: 35
```

Did you expect that in your bioinformatics software? Note that the mouse genome (GRCm38.p5) as far as I could tell does not contain any non-ACGT IUPAC letters

### rs SNP identifiers occuring in multiple places

Due to how dbSNP is creating, an rs SNP ID can occur in multiple places on the genome
https://www.biostars.org/p/2323/

### Weird characters in FASTA sequence names

In response to hg38 including a colon in sequence names, which conflicts with commonly used representation of a range as chr1:1-100 for example, people analyzed meta-character frequencies in sequence names https://github.com/samtools/hts-specs/issues/291

```
ENA
#   16927
*   1
,   231
-   122563947
.   521540419
/   236951
\   0
:   30181
;   72892
=   186611
@   3713
|   949

Broad(?)
     12 #
    527 *
    357 ,
1451132 -
1492749 .
  86114 /
 233731 :
   2034 =
     17 @
1735713 |

Reference sequences
 # 203
 % 203
 * 525
 + 1
 , 496
 - 154226
 . 1826561
 : 1577
 = 26
 _ 4961932
 | 1098333
```

Note that commas in FASTA names is being suggested as an illegal character because of the supplementary alignment tag in SAM/BAM using comma separated values

## Humongous chromosomes V1

Genomes such as wheat have large chromosomes averaging 806Mbp but the BAI file format is limited to 2^29-1 ~ 536Mbp in size (this is due to the binning strategy, the max bin size is listed as 2^29)

## Humongous chromosomes V2

The axolotl genome has individual chromosomes that are of size 3.14 Gbp https://genome.cshlp.org/content/29/2/317.long (2019) which is almost as big as the entire human genome

The BAM and CRAM formats can only store 2^31-1 length https://en.wikipedia.org/wiki/2,147,483,647 so bgzip/tabix SAM is used

## Largest genomes

Just some honorable mentions for largest genome

- Polychaos dubium/Amoeba dubium/Chaos chaos - ~600-1300Gbp (unsequenced, 1968 back of envelope measurement, needs confirmation) https://en.wikipedia.org/wiki/Polychaos_dubium (another ref https://bionumbers.hms.harvard.edu/bionumber.aspx?&id=117342)
- Dinoflagellates - up to 250Gbp (unsequenced, 1987 book referenced in this paper, needs confirmation, has weird chromosome "rod-like" structures) https://www.nature.com/articles/s41588-021-00841-y
- Paris japonica (canopy plant) - ~149Gbp (unsequenced)  https://en.wikipedia.org/wiki/Paris_japonica
- Tmesipteris_obliqua (fern) - ~147Gbp (unsequenced) - https://en.wikipedia.org/wiki/Tmesipteris_obliqua
- Marbled lungfish - ~133Gbp (unsequenced) https://en.wikipedia.org/wiki/Marbled_lungfish
- European mistletoe - ~90Gbp (partial sequence) https://onlinelibrary.wiley.com/doi/10.1111/tpj.15558
- Australian lungfish - ~43Gbp (sequenced)  https://www.smithsonianmag.com/smart-news/australian-lungfish-has-biggest-genome-ever-sequenced-180976837/
- Axolotl - ~32Gbp (sequenced) https://en.wikipedia.org/wiki/Axolotl
- Coastal redwood - ~26Gbp (sequenced) https://www.ucdavis.edu/climate/news/coast-redwood-and-sequoia-genome-sequences-completed
- Loblolly pine - ~22Gbp (sequenced) https://blogs.biomedcentral.com/on-biology/wp-content/uploads/sites/5/2014/03/genomelog030.jpg
- Wheat genome - ~17Gbp https://academic.oup.com/gigascience/article/6/11/gix097/4561661

Inspired by twitter thread https://twitter.com/PetrovADmitri/status/1506824610360168455

Also see http://www.genomesize.com/statistics.php?stats=entire#stats_top

## Humongous CIGAR strings

The CG tag was invented in order to store CIGAR strings longer than 64kb, since n_cigar_opt is a uint16

## Interesting gene names

- Tinman - https://en.wikipedia.org/wiki/Tinman_gene
- Sonic hedgehog (SHH) - https://en.wikipedia.org/wiki/Sonic_hedgehog
- Heart of glass (heg) - https://www.ncbi.nlm.nih.gov/pubmed/14680629
- Dracula (drc) - https://www.ncbi.nlm.nih.gov/pubmed/10985389
- Sleeping Beauty transposon - https://en.wikipedia.org/wiki/Sleeping_Beauty_transposon_system
- Skywalker protein - http://www.ebi.ac.uk/pdbe/entry/search/index?pubmed_id:27669036
- Time for coffee - http://www.plantcell.org/content/15/11/2719.abstract
- WTF - https://www.ebi.ac.uk/interpro/entry/IPR004982 https://www.sciencedaily.com/releases/2017/06/170620093209.htm
- Mothers against decapentaplegic - https://en.wikipedia.org/wiki/Mothers_against_decapentaplegic
- Saxophone (sax) - http://www.sdbonline.org/sites/fly/gene/saxophon.htm
- Beethovan (btv) - http://www.uniprot.org/uniprot/Q0E8P6
- Superman+kryptonite - https://en.wikipedia.org/wiki/Superman_(gene)
- Supervillin (SVIL) - https://www.uniprot.org/uniprot/O95425
- Wishful thinking (wit) - https://www.wikigenes.org/e/gene/e/44096.html
- Doublesex (dsx) - https://en.wikipedia.org/wiki/Doublesex
- Fruitless (fru) - https://en.wikipedia.org/wiki/Fruitless_(gene)
- Transformer (tra) - https://en.wikipedia.org/wiki/Transformer_(gene)
- Gypsy+Flamenco - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1206375/ also described in wiki https://en.wikipedia.org/wiki/Piwi-interacting_RNA#History_and_loci
- Jockey - http://flybase.org/reports/FBgn0015952.html
- Tigger - https://www.omim.org/entry/612972
- Nanog - celtic legend https://www.sciencedaily.com/releases/2003/06/030602024530.htm (source https://twitter.com/EpgntxEinstein/status/1057359656220348417)
- Jerky (jrk) - https://www.genecards.org/cgi-bin/carddisp.pl?gene=JRK
- Hippo (Hpo) - https://www.wikigenes.org/e/gene/e/37247.html
- Dishevelled (Dsh) - https://en.wikipedia.org/wiki/Dishevelled
- Glass bottom boat (gbb) - http://www.sdbonline.org/sites/fly/dbzhnsky/60a-1.htm
- Makes catepillars floppy (mcf) - https://www.pnas.org/content/99/16/10742 (source https://twitter.com/JUNIUS_64/status/1081007886560608256)
- Eyeless http://flybase.org/reports/FBgn0005558.html
- Straightjaket (stj) - http://flybase.org/reports/FBgn0261041.html
- Huluwa http://science.sciencemag.org/content/362/6417/eaat1045 ref https://twitter.com/zhouwanding/status/1065960714978897921
- frameshifts or pseudogene? - check sequence - https://www.ncbi.nlm.nih.gov/gene/?term=24562233%5Buid%5D
- Bad response to refridgeration (brr) https://twitter.com/hitenmadhani/status/1149471071675924481?s=20
- Mindbomb (mib1) - https://www.sdbonline.org/sites/fly/hjmuller/mindbomb1.htm
- β'COP http://flybase.org/reports/FBgn0025724.html (https://twitter.com/DarrenObbard/status/1260613447198412800)
- King-tubby https://www.uniprot.org/uniprot/B0XFQ9 see also tubby https://www.uniprot.org/uniprot/P50586
- fucK https://www.uniprot.org/uniprot/?query=fuck&sort=score
- Halloween genes https://en.wikipedia.org/wiki/Halloween_genes
- VANDAL21 https://www.arabidopsis.org/servlets/TairObject?type=transposon_family&id=139
- HotDog domain - superfamily of genes/proteins https://www.wikidata.org/wiki/Q24785143 https://www.ebi.ac.uk/interpro/entry/IPR029069
- Flower/fwe - https://flybase.org/reports/FBgn0261722.html
- Brahma https://www.sdbonline.org/sites/fly/polycomb/brahma.htm

### Allele names

Sometimes it is not the gene, but the allele that is named

- Bad hair day http://www.informatics.jax.org/allele/MGI:3764934
- Samba, chacha, bossa nova http://www.informatics.jax.org/allele/MGI:3764934
- Yoda http://www.informatics.jax.org/allele/MGI:3797584

Ref https://twitter.com/hmdc_mgi/status/1242893531779391496

## More reading

Great illustrations of interesting biology, including information about gene names https://twitter.com/vividbiology

Many of the stories behind fly gene nomenclature is available at https://web.archive.org/web/20110716201703/http://www.flynome.com/cgi-bin/search?source=browse including the famous ForRentApartments dot com gene (just kidding but lol https://web.archive.org/web/20110716202150/http://www.flynome.com/cgi-bin/search?storyID=180)

Musing article: "What is in a (gene) name?" https://web.archive.org/web/20180731060319/https://blogs.plos.org/toothandclaw/2012/06/17/whats-in-a-gene-name/

## Send PRs for more things!
