#include "htslib/hts.h"
#include "htslib/vcf.h"
#include <stdint.h>
#include <inttypes.h>

int32_t variant_rid(bcf1_t *b);
int64_t variant_pos(bcf1_t *b);
int64_t variant_rlen(bcf1_t *b);

char *variant_REF(bcf1_t *b);
char *variant_ALT(bcf1_t *b, int32_t i);

float variant_QUAL(bcf1_t *b);
int variant_nflt(bcf1_t *b);
char *variant_flt0(bcf1_t *b, bcf_hdr_t *h);

BGZF *fp_bgzf(htsFile *h);


int bcf_get_info_i32(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int32_t **dst, int *ndst);


