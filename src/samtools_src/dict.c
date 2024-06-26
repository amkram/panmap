/*  dict.c -- create a sequence dictionary file.

    Copyright (C) 2015, 2020 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <getopt.h>
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "samtools.h"

KHASH_SET_INIT_STR(str)
KSEQ_INIT(gzFile, gzread)

typedef struct _args_t
{
    char *output_fname, *alt_fname;
    char *assembly, *species, *uri;
    int  alias, header;
    khash_t(str) *is_alt;
}
args_t;

static void write_dict(const char *fn, args_t *args)
{
    hts_md5_context *md5;
    int l, i, k;
    gzFile fp;
    kseq_t *seq;
    unsigned char digest[16];
    char hex[33];

    fp = strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        print_error_errno("dict", "Cannot open %s", fn);
        exit(1);
    }
    FILE *out = stdout;
    if (args->output_fname) {
        out = fopen(args->output_fname, "w");
        if (out == NULL) {
          print_error_errno("dict", "Cannot open %s for writing", args->output_fname);
          exit(1);
        }
    }

    if (!(md5 = hts_md5_init()))
        exit(1);

    seq = kseq_init(fp);
    if (args->header) fprintf(out, "@HD\tVN:1.0\tSO:unsorted\n");
    while ((l = kseq_read(seq)) >= 0) {
        for (i = k = 0; i < seq->seq.l; ++i) {
            if (seq->seq.s[i] >= '!' && seq->seq.s[i] <= '~')
                seq->seq.s[k++] = toupper(seq->seq.s[i]);
        }
        hts_md5_reset(md5);
        hts_md5_update(md5, (unsigned char*)seq->seq.s, k);
        hts_md5_final(digest, md5);
        hts_md5_hex(hex, digest);
        fprintf(out, "@SQ\tSN:%s\tLN:%d\tM5:%s", seq->name.s, k, hex);
        if (args->is_alt && kh_get(str, args->is_alt, seq->name.s) != kh_end(args->is_alt))
            fprintf(out, "\tAH:*");
        if (args->alias) {
            const char *name = seq->name.s;
            if (strncmp(name, "chr", 3) == 0) {
                name += 3;
                fprintf(out, "\tAN:%s", name);
            }
            else
                fprintf(out, "\tAN:chr%s", name);

            if (strcmp(name, "M") == 0)
                fprintf(out, ",chrMT,MT");
            else if (strcmp(name, "MT") == 0)
                fprintf(out, ",chrM,M");
        }
        if (args->uri)
            fprintf(out, "\tUR:%s", args->uri);
        else if (strcmp(fn, "-") != 0) {
#ifdef _WIN32
            char *real_path = _fullpath(NULL, fn, PATH_MAX);
#else
            char *real_path = realpath(fn, NULL);
#endif
            fprintf(out, "\tUR:file://%s", real_path);
            free(real_path);
        }
        if (args->assembly) fprintf(out, "\tAS:%s", args->assembly);
        if (args->species) fprintf(out, "\tSP:%s", args->species);
        fprintf(out, "\n");
    }
    kseq_destroy(seq);
    hts_md5_destroy(md5);

    if (args->output_fname) fclose(out);
    gzclose(fp);
}

static void read_alt_file(khash_t(str) *is_alt, const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) {
        print_error_errno("dict", "Cannot open %s", fname);
        exit(1);
    }

    // .alt files are in a SAM-like format, but we don't use sam_read1()
    // as these files may not have a complete set of @SQ headers.

    kstring_t str = KS_INITIALIZE;
    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l == 0 || str.s[0] == '@') continue;

        char *tab = strchr(str.s, '\t');
        if (tab) *tab = '\0';

        int ret;
        char *seqname = strdup(str.s);
        kh_put(str, is_alt, seqname, &ret);
        if (ret == 0) free(seqname); // Already present
    }

    ks_free(&str);
    hts_close(fp);
}

static int dict_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Create a sequence dictionary file from a fasta file\n");
    fprintf(stderr, "Usage:   samtools dict [options] <file.fa|file.fa.gz>\n\n");
    fprintf(stderr, "Options: -a, --assembly STR    assembly\n");
    fprintf(stderr, "         -A, --alias, --alternative-name\n");
    fprintf(stderr, "                               add AN tag by adding/removing 'chr'\n");
    fprintf(stderr, "         -H, --no-header       do not print @HD line\n");
    fprintf(stderr, "         -l, --alt FILE        add AH:* tag to alternate locus sequences\n");
    fprintf(stderr, "         -o, --output FILE     file to write out dict file [stdout]\n");
    fprintf(stderr, "         -s, --species STR     species\n");
    fprintf(stderr, "         -u, --uri STR         URI [file:///abs/path/to/file.fa]\n");
    fprintf(stderr, "\n");
    return 1;
}

int dict_main(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->header = 1;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"no-header", no_argument, NULL, 'H'},
        {"alias", no_argument, NULL, 'A'},
        {"alt", required_argument, NULL, 'l'},
        {"alternative-name", no_argument, NULL, 'A'},
        {"assembly", required_argument, NULL, 'a'},
        {"species", required_argument, NULL, 's'},
        {"uri", required_argument, NULL, 'u'},
        {"output", required_argument, NULL, 'o'},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ( (c=getopt_long(argc,argv,"?AhHa:l:s:u:o:",loptions,NULL))>0 )
    {
        switch (c)
        {
            case 'A': args->alias = 1; break;
            case 'a': args->assembly = optarg; break;
            case 'l': args->alt_fname = optarg; break;
            case 's': args->species = optarg; break;
            case 'u': args->uri = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'H': args->header = 0; break;
            case 'h': return dict_usage();
            default: return dict_usage();
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(STDIN_FILENO) ) fname = "-";  // reading from stdin
        else return dict_usage();
    }
    else fname = argv[optind];

    if (args->alt_fname) {
        args->is_alt = kh_init(str);
        read_alt_file(args->is_alt, args->alt_fname);
    }

    write_dict(fname, args);

    if (args->is_alt) {
        khint_t k;
        for (k = 0; k < kh_end(args->is_alt); ++k)
            if (kh_exist(args->is_alt, k)) free((char *) kh_key(args->is_alt, k));
        kh_destroy(str, args->is_alt);
    }

    free(args);
    return 0;
}
