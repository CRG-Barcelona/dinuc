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

int get_prostar_deform_idx(char *word)
/* return index needed from "roll", etc. */
/* should be from 0-5 */
{
    if (sameWord(word, "twist"))
	return 0;
    if (sameWord(word, "tilt"))
	return 1;
    if (sameWord(word, "roll"))
	return 2;
    if (sameWord(word, "shift"))
	return 3;
    if (sameWord(word, "slide"))
	return 4;
    if (sameWord(word, "rise"))
	return 5;
    errAbort("You gave a bad word");
    return -1;
}

struct hash *get_prostar_hash()
{
    struct hash *ps_hash = hashNew(4);
    hashAddInt(ps_hash, "AA", 0);
    hashAddInt(ps_hash, "AC", 6);
    hashAddInt(ps_hash, "AG", 12);
    hashAddInt(ps_hash, "AT", 18);
    hashAddInt(ps_hash, "CA", 24);
    hashAddInt(ps_hash, "CC", 30);
    hashAddInt(ps_hash, "CG", 36);
    hashAddInt(ps_hash, "GA", 42);
    hashAddInt(ps_hash, "GC", 48);
    hashAddInt(ps_hash, "TA", 54);
    /* missing ones */
    hashAddInt(ps_hash, "CT", 12);  /* CT -> AG */
    hashAddInt(ps_hash, "GG", 30);  /* GG -> CC */
    hashAddInt(ps_hash, "GT",  6);  /* GT -> AC */  
    hashAddInt(ps_hash, "TC", 42);  /* TC -> GA */
    hashAddInt(ps_hash, "TG", 24);  /* TG -> CA */
    hashAddInt(ps_hash, "TT",  0);  /* TT -> AA */
    
    return ps_hash;
}

double get_prostar_val(struct hash *ps_hash, char *dinuc, int deform_idx)
{
    int ix_off = hashIntVal(ps_hash, dinuc);
    return prostar_mat[ix_off + deform_idx];
}
