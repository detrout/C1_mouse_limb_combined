CURLOPTS:= -L --http1.1

all: 15019 15020 15084 15085 16111 16112 16134 16135 16930 16931 17298 17299

15019: 15019/15019-mm10-M4-male_genome.bam 15019/15019-mm10-M4-male_anno_rsem.genes.results
15020: 15020/15020-mm10-M4-male_genome.bam 15020/15020-mm10-M4-male_anno_rsem.genes.results
15084: 15084/15084-mm10-M4-male_genome.bam 15084/15084-mm10-M4-male_anno_rsem.genes.results
15085: 15085/15085-mm10-M4-male_genome.bam 15085/15085-mm10-M4-male_anno_rsem.genes.results
16111: 16111/16111-mm10-M4-male_genome.bam 16111/16111-mm10-M4-male_anno_rsem.genes.results
16112: 16112/16112-mm10-M4-male_genome.bam 16112/16112-mm10-M4-male_anno_rsem.genes.results
16134: 16134/16134-mm10-M4-male_genome.bam 16134/16134-mm10-M4-male_anno_rsem.genes.results
16135: 16135/16135-mm10-M4-male_genome.bam 16135/16135-mm10-M4-male_anno_rsem.genes.results
16930: 16930/16930-mm10-M4-male_genome.bam 16930/16930-mm10-M4-male_anno_rsem.genes.results
16931: 16931/16931-mm10-M4-male_genome.bam 16931/16931-mm10-M4-male_anno_rsem.genes.results
17298: 17298/17298-mm10-M4-male_genome.bam 17298/17298-mm10-M4-male_anno_rsem.genes.results
17299: 17299/17299-mm10-M4-male_genome.bam 17299/17299-mm10-M4-male_anno_rsem.genes.results

15019/15019-mm10-M4-male_genome.bam: 15019/ENCFF748RBN.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
15020/15020-mm10-M4-male_genome.bam: 15020/ENCFF271VQQ.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
15084/15084-mm10-M4-male_genome.bam: 15084/ENCFF483ZFR.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
15085/15085-mm10-M4-male_genome.bam: 15085/ENCFF793SNG.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16111/16111-mm10-M4-male_genome.bam: 16111/ENCFF943SLL.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16112/16112-mm10-M4-male_genome.bam: 16112/ENCFF808VGK.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16134/16134-mm10-M4-male_genome.bam: 16134/ENCFF445SGU.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16135/16135-mm10-M4-male_genome.bam: 16135/ENCFF223JGE.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16930/16930-mm10-M4-male_genome.bam: 16930/ENCFF766OUK.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16931/16931-mm10-M4-male_genome.bam: 16931/ENCFF188HWM.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
17298/17298-mm10-M4-male_genome.bam: 17298/ENCFF309GAT.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
17299/17299-mm10-M4-male_genome.bam: 17299/ENCFF836HAC.bam
	cd $(dir $^) ; ln -s $(notdir $^) ../$@

15019/ENCFF748RBN.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF748RBN/@@download/ENCFF748RBN.bam
15020/ENCFF271VQQ.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF271VQQ/@@download/ENCFF271VQQ.bam
15084/ENCFF483ZFR.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF483ZFR/@@download/ENCFF483ZFR.bam
15085/ENCFF793SNG.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF793SNG/@@download/ENCFF793SNG.bam
16111/ENCFF943SLL.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF943SLL/@@download/ENCFF943SLL.bam
16112/ENCFF808VGK.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF808VGK/@@download/ENCFF808VGK.bam
16134/ENCFF445SGU.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF445SGU/@@download/ENCFF445SGU.bam
16135/ENCFF223JGE.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF223JGE/@@download/ENCFF223JGE.bam
16930/ENCFF766OUK.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF766OUK/@@download/ENCFF766OUK.bam
16931/ENCFF188HWM.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF188HWM/@@download/ENCFF188HWM.bam
17298/ENCFF309GAT.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF309GAT/@@download/ENCFF309GAT.bam
17299/ENCFF836HAC.bam:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/ENCFF836HAC/@@download/ENCFF836HAC.bam


15019/ENCFF195JHC.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
15020/ENCFF457ZGF.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
15084/ENCFF677BPV.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
15085/ENCFF794QMH.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
16111/ENCFF358WYS.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
16112/ENCFF634AUL.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
16134/ENCFF532ZDE.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
16135/ENCFF003DBZ.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
16930/ENCFF879FXB.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
16931/ENCFF470WZZ.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
17298/ENCFF399HJB.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)
17299/ENCFF597VHC.tsv:
	curl $(CURLOPTS) -o $@ https://www.encodeproject.org/files/$(basename $(notdir $@))/@@download/$(notdir $@)


15019/15019-mm10-M4-male_anno_rsem.genes.results: 15019/ENCFF195JHC.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
15020/15020-mm10-M4-male_anno_rsem.genes.results: 15020/ENCFF457ZGF.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
15084/15084-mm10-M4-male_anno_rsem.genes.results: 15084/ENCFF677BPV.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
15085/15085-mm10-M4-male_anno_rsem.genes.results: 15085/ENCFF794QMH.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16111/16111-mm10-M4-male_anno_rsem.genes.results: 16111/ENCFF358WYS.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16112/16112-mm10-M4-male_anno_rsem.genes.results: 16112/ENCFF634AUL.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16134/16134-mm10-M4-male_anno_rsem.genes.results: 16134/ENCFF532ZDE.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16135/16135-mm10-M4-male_anno_rsem.genes.results: 16135/ENCFF003DBZ.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16930/16930-mm10-M4-male_anno_rsem.genes.results: 16930/ENCFF879FXB.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
16931/16931-mm10-M4-male_anno_rsem.genes.results: 16931/ENCFF470WZZ.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
17298/17298-mm10-M4-male_anno_rsem.genes.results: 17298/ENCFF399HJB.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
17299/17299-mm10-M4-male_anno_rsem.genes.results: 17299/ENCFF597VHC.tsv
	cd $(dir $^) ; ln -s $(notdir $^) ../$@
