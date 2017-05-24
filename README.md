# oddgenes

A list of weird gene annotations or things that break bioinformatics assumptions


## Gene structures

### 1bp length exon

Evidence given for a 1bp length exon in Arabadopsis and different splicing models are discussed

http://www.nature.com/articles/srep18087

### 0bp length exon

The phenomenon of recursive splicing can remove sequences progressively inside an intron, so there can exist "0bp exons" that are just the splice-site sequences pasted together.

"To identify potential zero nucleotide exon-type ratchet points, we parsed the RNA-Seq alignments to identify novel splice junctions where the reads mapped to an annotated 5' splice site and an unannotated 3' splice site, and the genomic sequence at the 3' splice site junction was AG/GT"

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4529404/

### Very large number of isoforms in Dscam

"Dscam has 24 exons; exon 4 has 12 variants, exon 6 has 48 variants, exon 9 has 33 variants, and exon 17 has two variants. The combination of exons 4, 6, and 9 leads to 19,008 possible isoforms with different extracellular domains (due to differences in Ig2, Ig3 and Ig4). With two different transmembrane domains from exon 17, the total possible protein products could reach 38,016 isoforms"

Ref https://en.wikipedia.org/wiki/DSCAM https://www.wikigenes.org/e/gene/e/35652.html


## Wormbase

### Adding leader sequence to mRNA

"About 70% of C. elegans mRNAs are trans-spliced to one of two 22 nucleotide spliced leaders. SL1 is used to trim off the 5' ends of pre-mRNAs and replace them with the SL1 sequence. This processing event is very closely related to cis-splicing, or intron removal."

The region that is spliced out is called an outron

http://www.wormbook.org/chapters/www_transsplicingoperons/transsplicingoperons.html

### Polycistronic transcripts/operons

"A characteristic feature of the worm genome is the existence of genes organized into operons. These polycistronic gene clusters contain two or more closely spaced genes, which are oriented in a head to tail direction. They are transcribed as a single polycistronic mRNA and separated into individual mRNAs by the process of trans-splicing"

http://www.wormbook.org/chapters/www_overviewgenestructure.2/genestructure.html

## Evolution

### Possible adaptive bacteria->eukaryote HGT

A possible horizontal gene transfer from bacteria to eukaryotes is found in an insect that feeds on coffee beans. Changes that the gene had to undergo are covered (added poly-A tail, shine-dalgarno sequence deleted)

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3306691/

## Codon usage

### Alternative start codons

"The most common start codons for known Escherichia coli genes are AUG (83% of genes), GUG (14%) and UUG (3%)"

"Here, we systematically quantified translation initiation of green fluorescent protein (GFP) from all 64 codons and nanoluciferase from 12 codons on plasmids designed to interrogate a range of translation initiation conditions."

https://www.sciencedaily.com/releases/2017/02/170221080506.htm


## Transposons

### Cross-species BovB transposon transfers

Or "How a quarter of the cow genome came from snakes" http://phenomena.nationalgeographic.com/2013/01/01/how-a-quarter-of-the-cow-genome-came-from-snakes/

Source http://www.pnas.org/content/110/3/1012.full


## Regulation

### Intron mediated enhancement (IME)

It has been shown that some intron sequences can enhance expression similar to how promoter sequences work https://en.wikipedia.org/wiki/Intron-mediated_enhancement

The first intron of the UBQ10 gene in Arabidopsis exhibits IME, and "the sequences responsible for increasing mRNA accumulation are redundant and dispersed throughout the UBQ10 intron" http://www.plantcell.org/content/early/2017/04/03/tpc.17.00020.full.pdf+html


## Non-ACGT letters in fasta files

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

Full listing

```
NC_000001.11	248752513	248752514	NC_000001.11_248752513_248752514_for	1	+	M
NC_000001.11	248755121	248755122	NC_000001.11_248755121_248755122_for	1	+	R
NC_000002.12	20953874	20953875	NC_000002.12_20953874_20953875_for	1	+	Y
NC_000002.12	109493617	109493618	NC_000002.12_109493617_109493618_for	1	+	R
NC_000002.12	233142969	233142970	NC_000002.12_233142969_233142970_for	1	+	Y
NC_000002.12	233143702	233143703	NC_000002.12_233143702_233143703_for	1	+	W
NC_000002.12	233143708	233143709	NC_000002.12_233143708_233143709_for	1	+	Y
NC_000002.12	233143728	233143729	NC_000002.12_233143728_233143729_for	1	+	Y
NC_000002.12	233144846	233144847	NC_000002.12_233144846_233144847_for	1	+	K
NC_000002.12	233145085	233145086	NC_000002.12_233145085_233145086_for	1	+	W
NC_000002.12	239866290	239866291	NC_000002.12_239866290_239866291_for	1	+	M
NC_000003.12	16814154	16814155	NC_000003.12_16814154_16814155_for	1	+	R
NC_000003.12	16816463	16816464	NC_000003.12_16816463_16816464_for	1	+	Y
NC_000003.12	16894810	16894811	NC_000003.12_16894810_16894811_for	1	+	W
NC_000003.12	16902882	16902883	NC_000003.12_16902882_16902883_for	1	+	B
NC_000003.12	66186452	66186453	NC_000003.12_66186452_66186453_for	1	+	Y
NC_000003.12	66191220	66191221	NC_000003.12_66191220_66191221_for	1	+	Y
NC_000003.12	90549738	90549739	NC_000003.12_90549738_90549739_for	1	+	W
NC_000006.12	61366326	61366327	NC_000006.12_61366326_61366327_for	1	+	Y
NC_000007.14	154575458	154575459	NC_000007.14_154575458_154575459_for	1	+	R
NC_000007.14	154578613	154578614	NC_000007.14_154578613_154578614_for	1	+	Y
NC_000007.14	154578614	154578615	NC_000007.14_154578614_154578615_for	1	+	S
NC_000007.14	154578678	154578679	NC_000007.14_154578678_154578679_for	1	+	Y
NC_000009.12	89883951	89883952	NC_000009.12_89883951_89883952_for	1	+	R
NC_000009.12	89884095	89884096	NC_000009.12_89884095_89884096_for	1	+	Y
NC_000009.12	89884246	89884247	NC_000009.12_89884246_89884247_for	1	+	Y
NC_000010.11	39239442	39239443	NC_000010.11_39239442_39239443_for	1	+	R
NC_000010.11	39240084	39240085	NC_000010.11_39240084_39240085_for	1	+	R
NC_000010.11	39240481	39240482	NC_000010.11_39240481_39240482_for	1	+	R
NC_000010.11	39242507	39242508	NC_000010.11_39242507_39242508_for	1	+	Y
NC_000010.11	39244443	39244444	NC_000010.11_39244443_39244444_for	1	+	W
NC_000010.11	39245282	39245283	NC_000010.11_39245282_39245283_for	1	+	M
NC_000010.11	39248482	39248483	NC_000010.11_39248482_39248483_for	1	+	W
NC_000010.11	39250360	39250361	NC_000010.11_39250360_39250361_for	1	+	W
NC_000010.11	39284361	39284362	NC_000010.11_39284361_39284362_for	1	+	Y
NC_000010.11	39285661	39285662	NC_000010.11_39285661_39285662_for	1	+	M
NC_000010.11	39285854	39285855	NC_000010.11_39285854_39285855_for	1	+	K
NC_000010.11	39289867	39289868	NC_000010.11_39289867_39289868_for	1	+	Y
NC_000010.11	39290415	39290416	NC_000010.11_39290415_39290416_for	1	+	Y
NC_000010.11	39300295	39300296	NC_000010.11_39300295_39300296_for	1	+	W
NC_000010.11	39307179	39307180	NC_000010.11_39307179_39307180_for	1	+	Y
NC_000010.11	39307944	39307945	NC_000010.11_39307944_39307945_for	1	+	K
NC_000010.11	39308364	39308365	NC_000010.11_39308364_39308365_for	1	+	K
NC_000010.11	39311772	39311773	NC_000010.11_39311772_39311773_for	1	+	M
NC_000010.11	39311797	39311798	NC_000010.11_39311797_39311798_for	1	+	Y
NC_000010.11	39312243	39312244	NC_000010.11_39312243_39312244_for	1	+	R
NC_000010.11	39315628	39315629	NC_000010.11_39315628_39315629_for	1	+	R
NC_000010.11	39315771	39315772	NC_000010.11_39315771_39315772_for	1	+	R
NC_000010.11	39316937	39316938	NC_000010.11_39316937_39316938_for	1	+	W
NC_000010.11	39321100	39321101	NC_000010.11_39321100_39321101_for	1	+	R
NC_000010.11	124120673	124120674	NC_000010.11_124120673_124120674_for	1	+	R
NC_000010.11	124120891	124120892	NC_000010.11_124120891_124120892_for	1	+	Y
NC_000010.11	124121818	124121819	NC_000010.11_124121818_124121819_for	1	+	Y
NC_000010.11	124122079	124122080	NC_000010.11_124122079_124122080_for	1	+	B
NC_000010.11	124122245	124122246	NC_000010.11_124122245_124122246_for	1	+	R
NC_000010.11	124122267	124122268	NC_000010.11_124122267_124122268_for	1	+	R
NC_000010.11	131592410	131592411	NC_000010.11_131592410_131592411_for	1	+	R
NC_000010.11	131594305	131594306	NC_000010.11_131594305_131594306_for	1	+	R
NC_000010.11	131594830	131594831	NC_000010.11_131594830_131594831_for	1	+	K
NC_000010.11	131595664	131595665	NC_000010.11_131595664_131595665_for	1	+	W
NC_000010.11	131595827	131595828	NC_000010.11_131595827_131595828_for	1	+	R
NC_000010.11	131596142	131596143	NC_000010.11_131596142_131596143_for	1	+	S
NC_000012.12	122093259	122093260	NC_000012.12_122093259_122093260_for	1	+	Y
NC_000012.12	132222633	132222634	NC_000012.12_132222633_132222634_for	1	+	Y
NC_000012.12	132225196	132225197	NC_000012.12_132225196_132225197_for	1	+	M
NC_000013.11	100972570	100972571	NC_000013.11_100972570_100972571_for	1	+	Y
NC_000013.11	100973392	100973393	NC_000013.11_100973392_100973393_for	1	+	K
NC_000013.11	100973394	100973395	NC_000013.11_100973394_100973395_for	1	+	Y
NC_000016.10	88366824	88366825	NC_000016.10_88366824_88366825_for	1	+	R
NC_000021.9	41829695	41829696	NC_000021.9_41829695_41829696_for	1	+	M
NC_000021.9	41830238	41830239	NC_000021.9_41830238_41830239_for	1	+	M
NC_000021.9	41830507	41830508	NC_000021.9_41830507_41830508_for	1	+	R
NC_000022.11	16309655	16309656	NC_000022.11_16309655_16309656_for	1	+	R
NC_000022.11	16320468	16320469	NC_000022.11_16320468_16320469_for	1	+	W
NC_000022.11	16325089	16325090	NC_000022.11_16325089_16325090_for	1	+	Y
NC_000022.11	16326809	16326810	NC_000022.11_16326809_16326810_for	1	+	R
NC_000022.11	16327065	16327066	NC_000022.11_16327065_16327066_for	1	+	Y
NC_000023.11	222581	222582	NC_000023.11_222581_222582_for	1	+	Y
NC_000023.11	222862	222863	NC_000023.11_222862_222863_for	1	+	Y
NC_000023.11	222975	222976	NC_000023.11_222975_222976_for	1	+	R
NC_000023.11	222979	222980	NC_000023.11_222979_222980_for	1	+	S
NC_000023.11	222989	222990	NC_000023.11_222989_222990_for	1	+	W
NC_000024.10	222581	222582	NC_000024.10_222581_222582_for	1	+	Y
NC_000024.10	222862	222863	NC_000024.10_222862_222863_for	1	+	Y
NC_000024.10	222975	222976	NC_000024.10_222975_222976_for	1	+	R
NC_000024.10	222979	222980	NC_000024.10_222979_222980_for	1	+	S
NC_000024.10	222989	222990	NC_000024.10_222989_222990_for	1	+	W
```

## Interesting gene names

* Tinman - https://en.wikipedia.org/wiki/Tinman_gene
* Sonic hedgehog (SHH) - https://en.wikipedia.org/wiki/Sonic_hedgehog
* Heart of glass (heg) - https://www.ncbi.nlm.nih.gov/pubmed/14680629
* Dracula (drc) - https://www.ncbi.nlm.nih.gov/pubmed/10985389
* Sleeping Beauty transposon - https://en.wikipedia.org/wiki/Sleeping_Beauty_transposon_system
* Skywalker protin - http://www.ebi.ac.uk/pdbe/entry/search/index?pubmed_id:27669036
* Time for coffee - http://www.plantcell.org/content/15/11/2719.abstract

Great illustrations of interesting biology, including information about gene names https://twitter.com/vividbiology
