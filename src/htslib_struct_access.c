#include "htslib_struct_access.h"

int32_t variant_rid(bcf1_t *b) {
  return b->rid;
}
int64_t variant_pos(bcf1_t *b) {
  return b->pos;
}
int64_t variant_rlen(bcf1_t *b) {
  return b->rlen;
}

char *variant_REF(bcf1_t *b) {
  return b->d.allele[0];
}

char *variant_ALT(bcf1_t *b, int32_t i) {
  if(i >= b->n_allele - 1) { return NULL; }
  return b->d.allele[i+1];
}

float variant_QUAL(bcf1_t *b) {
  return b->qual;
}

int variant_nflt(bcf1_t *b) {
  return b->d.n_flt;
}
char *variant_flt0(bcf1_t *b, bcf_hdr_t *h) {
  return bcf_hdr_int2id(h, BCF_DT_ID, b->d.flt[0]);
}

int32_t variant_n_samples(bcf1_t *b) {
	return b->n_sample;
}

bool is_vcf(htsFile *hts) {
	return hts->format.format == vcf;
}

int32_t variant_errcode(bcf1_t *b) {
	return (int32_t)b->errcode;
}

int32_t header_n_samples(bcf_hdr_t *h) {
	return bcf_hdr_nsamples(h);
}

BGZF *fp_bgzf(htsFile *h) {
  return h->fp.bgzf;
}

int bcf_get_info_i32(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int32_t **dst, int *ndst) {
  return bcf_get_info_values(hdr, line, tag, (void **)dst, ndst, BCF_HT_INT);
}

