#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* stuff about dinucleotides */

#include <float.h>
#include <math.h>
#include <jkweb/common.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/dnaLoad.h>
#include <jkweb/basicBed.h>
#include <jkweb/rangeTree.h>
#include <beato/bigs.h>

#include "prostar_consts.h"
#include "prostar.h"

void usage()
/* Explain usage and exit. */
{
errAbort(
  "dinucs - dinucleotide counting\n"
  "usage:\n"
  "   dinucs tag sequence_file output.wig\n"
  "options:\n"
  "   -prostar=word,w     get the deformation values of roll or tilt, etc with\n"
  "                       word=\"roll\", for example.  w is the window size\n"
  "   -wigtype=val        where val is fix, bg, or var\n"
    );
}

#define NANUM sqrt(-1)

static struct optionSpec options[] = {
    {"wigtype", OPTION_STRING},
    {"prostar", OPTION_STRING},
    {NULL, 0},
};

void delete_short_subsections(struct bed **pBedList, int min_len)
/* delete all the subsections shorter than 500 or -min-length bases */
{
    struct bed *section = NULL;
    struct bed *newList = NULL;
    while ((section = slPopHead(pBedList)) != NULL)
	if (section->chromStart - section->chromEnd >= min_len)
	    slAddHead(&newList, section);
	else
	    bedFree(&section);
    slReverse(&newList);
    *pBedList = newList;
}

void do_classic(struct dnaSeq *seqs, char *tag, enum wigOutType wot, FILE *out)
/* The typical algorithm */
{
    struct dnaSeq *seq;
    int tag_len = strlen(tag);
    int tag_mid = tag_len/2;
    int i, j;
    for (seq = seqs; seq != NULL; seq = seq->next)
    {
	struct bed *sec;
	struct perBaseWig *pbw = NULL;
	touppers(seq->dna);
	pbw = alloc_perBaseWig_matchingSequence(seq, TRUE);
	delete_short_subsections(&pbw->subsections, tag_len);
	for (sec = pbw->subsections; sec != NULL; sec = sec->next)
	{
	    for (i = sec->chromStart; i <= sec->chromEnd - tag_len; i++)
	    {
		boolean found = TRUE;
		for (j = 0; j < tag_len; j++)
		    if (seq->dna[i+j] != tag[j])
		    {
			found = FALSE;
			break;
		    }
		pbw->data[i+tag_mid] = (found) ? 1 : 0;
	    }
	}
	perBaseWigOutputNASkip(pbw, out, wot, 0, NULL, FALSE, TRUE);
	perBaseWigFree(&pbw);
    }
}

void do_prostar(struct dnaSeq *seqs, char *tag, enum wigOutType wot, FILE *out, int w)
/* The algorithm for finding basepair deformation values */
{
    struct dnaSeq *seq;
    struct hash *ps_hash = get_prostar_hash();
    double na = NANUM;
    int mid = w/2;
    int ps_idx = get_prostar_deform_idx(tag);
    int i, j;
    for (seq = seqs; seq != NULL; seq = seq->next)
    {
	struct bed *sec;
	struct perBaseWig *pbw = NULL;
	touppers(seq->dna);
	pbw = alloc_perBaseWig_matchingSequence(seq, TRUE);
	for (sec = pbw->subsections; sec != NULL; sec = sec->next)
	{
	    for (i = sec->chromStart; i <= sec->chromEnd - w; i++)
	    {
		double prostar_sum = 0;
		for (j = 0; j <= w-2; j++)
		{
		    char swap = seq->dna[i+j+2];
		    seq->dna[i+j+2] = '\0';
		    double ps_val = get_prostar_val(ps_hash, seq->dna + i + j, ps_idx);
		    prostar_sum += ps_val;
		    seq->dna[i+j+2] = swap;
		}
		if (w % 2 == 0)
		    pbw->data[i+mid-1] = prostar_sum/(w-1);
		else
		    pbw->data[i+mid] = prostar_sum/(w-1);
	    }
	    /* fill in zeroes with NAs */
	    int extra_for_even = (w % 2 == 0) ? 1 : 0;
	    for (i = sec->chromStart; i < mid-extra_for_even; i++)
		pbw->data[i] = na;
	    for (i = sec->chromEnd - mid; i < sec->chromEnd; i++)
		pbw->data[i] = na;
	}
	perBaseWigOutputNASkip(pbw, out, wot, 4, NULL, FALSE, TRUE);
	perBaseWigFree(&pbw);
    }
    hashFree(&ps_hash);
}

void dinucs(char *tag, char *input, char *output)
/* find nucleotide signals */
{
    char bases[5] = "ACGT";
    struct dnaSeq *seqs = dnaLoadAll(input);
    enum wigOutType wot = get_wig_out_type((char *)optionVal("wigtype", "fix"));
    boolean doo_prostar = FALSE;
    int prostar_w = 0;
    if (optionExists("prostar"))
    {
	doo_prostar = TRUE;
	prostar_w = optionInt("prostar", 2);
    }
    FILE *out = mustOpen(output, "w");
    touppers(tag);
    if (doo_prostar)
	do_prostar(seqs, tag, wot, out, prostar_w);
    else
	do_classic(seqs, tag, wot, out);
    carefulClose(&out);
    dnaSeqFreeList(&seqs);
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 4)
    usage();
dinucs(argv[1], argv[2], argv[3]);
return 0;
}
