/*  convert.h -- functions for converting between VCF/BCF and related formats.

    Copyright (C) 2014-2024 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#ifndef __CONVERT_H__
#define __CONVERT_H__

#include "../samtools/htslib-1.20/htslib/vcf.h"

typedef struct _convert_t convert_t;
enum convert_option
{
    allow_undef_tags,       // see `bcftools query --allow-undef-tags`, throws an error if tag is not defined otherwise
    subset_samples,         // in bracketed expressions (e.g. [ %GT]) consider only marked samples
    header_samples,         // include sample name in bracketed tags (e.g. SAMPLE1:GT SAMPLE2:GT for [ %GT])
    force_newline,          // automatically insert a newline when not part of the formatting expression
    print_filtered,         // print the provided string instead of discarding samples not included in subset_samples
    no_hdr_indices,         // drop column indices when printing header, i.e. "#CHROM", not "#[1]CHROM"
};

convert_t *convert_init(bcf_hdr_t *hdr, int *samples, int nsamples, const char *str);
void convert_destroy(convert_t *convert);
int convert_set_option(convert_t *convert, enum convert_option opt, ...);
int convert_header(convert_t *convert, kstring_t *str);
int convert_line(convert_t *convert, bcf1_t *rec, kstring_t *str);
int convert_max_unpack(convert_t *convert);
int convert_is_tag_used(convert_t *convert, char *tag);
const char **convert_list_used_tags(convert_t *convert, int *ntags);

#endif

