#ifndef PROSTAR_H
#define PROSTAR_H

int get_prostar_deform_idx(char *word);
/* return index needed from "roll", etc. */
/* should be from 0-5 */

struct hash *get_prostar_hash();

double get_prostar_val(struct hash *ps_hash, char *dinuc, int deform_idx);

#endif
