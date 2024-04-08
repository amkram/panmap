#include "align.hpp"
#include "genotype.hpp"
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include "mmpriv.h"
#include "minimap.h"
#include "sample.h"

using namespace seeding;
using namespace tree;


/* from minimap2 bam_plcmd.c */
#define dummy_free(p)

#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_SMART_OVERLAPS (1<<10)

#define MPLP_PRINT_MAPQ_CHAR (1<<11)
#define MPLP_PRINT_QPOS  (1<<12)
#define MPLP_PRINT_QNAME (1<<13)
#define MPLP_PRINT_FLAG  (1<<14)
#define MPLP_PRINT_RNAME (1<<15)
#define MPLP_PRINT_POS   (1<<16)
#define MPLP_PRINT_MAPQ  (1<<17)
#define MPLP_PRINT_CIGAR (1<<18)
#define MPLP_PRINT_RNEXT (1<<19)
#define MPLP_PRINT_PNEXT (1<<20)
#define MPLP_PRINT_TLEN  (1<<21)
#define MPLP_PRINT_SEQ   (1<<22)
#define MPLP_PRINT_QUAL  (1<<23)
#define MPLP_PRINT_MODS  (1<<24)
#define MPLP_PRINT_QPOS5 (1<<25)

#define MPLP_PRINT_LAST  (1<<26) // terminator for loop

#define MPLP_MAX_DEPTH 8000
#define MPLP_MAX_INDEL_DEPTH 250

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}


typedef struct sam_global_args {
    htsFormat in;
    htsFormat out;
    char *reference;
    int nthreads;
    int write_index;
} sam_global_args;


typedef struct {
    char *ref[2];
    int ref_id[2];
    hts_pos_t ref_len[2];
} mplp_ref_t;


typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, all, rev_del;
    int rflag_require, rflag_filter;
    char *reg, *pl_list, *fai_fname, *output_fname;
    faidx_t *fai;
    void *bed, *rghash, *auxlist;
    int argc;
    char **argv;
    char sep, empty, no_ins, no_ins_mods, no_del, no_ends;
    sam_global_args ga;
} mplp_conf_t;

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    sam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;


extern "C" {
    // minimap2 symbols
    mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos);
    mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos);
    int mpileup(mplp_conf_t *conf, int n, char **fn, char **fn_idx);
    void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a);
    mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a);

    // samtools symbols
    void sam_global_args_init(sam_global_args *ga);
}

//takes as input:
//  a pointer to a header type (sam_hdr_t)
//  an array of (bam1_t *)
//  an int for length of array
//  a pointer to a reference string
//  an int for ref string length
//
//  a pointer to an empty kstring_t for output (this will be populated with the mpileup string)
//
//destroys the header and nothing else
void bam_and_ref_to_mplp(sam_hdr_t *header, bam1_t **bam_lines, int nbams, char *ref_string, int lref, kstring_t *mplp_string) {

    char *global_ref_string;
    hts_pos_t global_ref_length;
    int bam_index;
    int bams;
    bam1_t **bam_records;

    mplp_conf_t mplp;
    memset(&mplp, 0, sizeof(mplp_conf_t));
    mplp.min_baseQ  = 13;
    mplp.capQ_thres = 0;
    mplp.max_depth  = MPLP_MAX_DEPTH;
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
    mplp.argc = 4; mplp.argv = NULL;
    mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    mplp.output_fname = NULL;
    mplp.all = 0;
    mplp.rev_del = 0;
    mplp.sep = ',';
    mplp.empty = '*';
    sam_global_args_init(&mplp.ga);

    bam_index = 0;
    bams = nbams;
    bam_records = bam_lines;

    global_ref_string = ref_string;
    global_ref_length = lref;

    char *fn = {(char *) "panmap-test.bam"};  
    char *fni = {(char *) "panmap-test.bam.bai"};  

    mpileup(&mplp, nbams, &fn, &fni);
}

// mi: holds reference info and alignment flags and options
// read_length
// read_seq: sequence of read, as a C string of A C G T or U
// n_regs: output of function, number of alignment hits
// regs: output of function, array of alignment hits, must free (*regs[i])->p using free(), must free (*regs) using free()
// b: thread buffer, just initialize using mm_tbuf_init() before all alignments, then afterward destroy using mm_tbuf_destroy()
// opt: alignment options
// seed_l: length of seeds
// n_seeds: number of seed hits in this read
// r_poss: position of seeds in reference coords (the ends of the seed, the last index included in the seed)
// q_poss: position of seeds in query(read) coords (the ends of the seed, the last index included in the seed)
void align_read_given_seeds(const mm_idx_t *mi,const int read_length,const char *read_seq, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, int seed_l, int64_t n_seeds, int *r_poss, int *q_poss, uint8_t *reversed) 
{
	int i, j, rep_len = 0, qlen_sum = 0, n_regs0, n_mini_pos = 0;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;
	float chn_pen_gap, chn_pen_skip;

	n_regs[0] = 0;  //Initialize outputs
	regs[0] = 0;    //

	qlen_sum += read_length;

	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	//mini_pos = (uint64_t*)kmalloc(b->km, n_seeds * sizeof(uint64_t));
	//kfree(b->km, mv.a);
	//mv.n = n_seeds;
	//mv.m = n_seeds;
	//kv_resize(mm128_t, b->km, mv, mv.n);
	//kfree(b->km, mini_pos);

	a = (mm128_t*)kmalloc(b->km, n_seeds * sizeof(mm128_t));

	for (i = 0; i < n_seeds; ++i){
		uint32_t qpos = (uint32_t)q_poss[i];            //query position, position on read coordinates, position of the END of the seed, including it. so if this is 6 and read[6] is A then A is the last character of the seed
		uint32_t rpos = (uint32_t)r_poss[i];            //reference position of the END of the seed
		uint32_t span = (uint32_t)seed_l;               //seed length
		
		uint64_t x, y;
		
		x = ((uint64_t)reversed[i]<<63) | rpos;
		y = (uint64_t)span << 32 | qpos;

		a[i].x = x;
		a[i].y = y;
	}
	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	chn_pen_gap  = opt->chain_gap_scale * 0.01 * mi->k;
	chn_pen_skip = opt->chain_skip_scale * 0.01 * mi->k;
	if (opt->flag & MM_F_RMQ) {
		a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
						  chn_pen_gap, chn_pen_skip, n_seeds, a, &n_regs0, &u, b->km);
	} else {
		a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
						 chn_pen_gap, chn_pen_skip, is_splice, 1, n_seeds, a, &n_regs0, &u, b->km);
	}

	if (opt->bw_long > opt->bw && (opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN)) == 0 && n_regs0 > 1) { // re-chain/long-join for long sequences
		int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
		if (qlen_sum - (en - st) > opt->rmq_rescue_size || en - st > qlen_sum * opt->rmq_rescue_ratio) {
			int32_t i;
			for (i = 0, n_seeds = 0; i < n_regs0; ++i) n_seeds += (int32_t)u[i];
			kfree(b->km, u);
			radix_sort_128x(a, a + n_seeds);
			a = mg_lchain_rmq(opt->max_gap, opt->rmq_inner_dist, opt->bw_long, opt->max_chain_skip, opt->rmq_size_cap, opt->min_cnt, opt->min_chain_score,
							  chn_pen_gap, chn_pen_skip, n_seeds, a, &n_regs0, &u, b->km);
		}
	} else if (opt->max_occ > opt->mid_occ && rep_len > 0 && !(opt->flag & MM_F_RMQ)) { // re-chain, mostly for short reads
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < 1)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, NULL, &mv, qlen_sum, &n_seeds, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, NULL, &mv, qlen_sum, &n_seeds, &rep_len, &n_mini_pos, &mini_pos);
			a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score,
							 chn_pen_gap, chn_pen_skip, is_splice, 1, n_seeds, a, &n_regs0, &u, b->km);
		}
	}
	b->frag_gap = max_chain_gap_ref;
	b->rep_len = rep_len;

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a, !!(opt->flag&MM_F_QSTRAND));
	if (mi->n_alt) {
		mm_mark_alt(mi, n_regs0, regs0);
		mm_hit_sort(b->km, &n_regs0, regs0, opt->alt_drop); // this step can be merged into mm_gen_regs(); will do if this shows up in profile
	}

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, 1, &read_length, &n_regs0, regs0, a);
	if (!is_sr && !(opt->flag&MM_F_QSTRAND)) {

		mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);
		n_regs0 = mm_filter_strand_retained(n_regs0, regs0);
	}

	regs0 = align_regs(opt, mi, b->km, read_length, read_seq, &n_regs0, regs0, a);
	regs0 = (mm_reg1_t*)realloc(regs0, sizeof(*regs0) * n_regs0);
	mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
	n_regs[0] = n_regs0, regs[0] = regs0;
	 
	//kfree(b->km, mv.a);
	kfree(b->km, a); 
	kfree(b->km, u);
	//kfree(b->km, mini_pos);

	if (b->km) {
		km_stat(b->km, &kmst);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > opt->cap_kalloc)) {
			if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
				fprintf(stderr, "[W::%s] reset thread-local memory after read %s\n", __func__, (char *)NULL);
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

void align::mapToTarget(Tree *T, std::unordered_map<std::string, std::vector<int32_t>> &seedToRefPositions, std::vector<std::vector<seed>> &readSeeds, std::vector<std::string> &readSequences, std::vector<std::string> &readQuals, std::vector<std::string> &readNames, std::string &bestMatchSequence, int k, std::string &samFileName, std::string &bamFileName, std::string &mpileupFileName, std::string &vcfFileName) {
    //Collecting reference Seeds
    std::vector<seed> refSeeds;
    for(auto kv : seedToRefPositions) {
        for (int32_t pos : kv.second) {
            seed thisSeed;
            thisSeed.seq = kv.first;
            thisSeed.pos = pos;
            refSeeds.push_back(thisSeed);
        }
    }

    //Sorting ref seeds
    sort(refSeeds.begin(), refSeeds.end(), []( const seed& lhs, const seed& rhs ) { return lhs.seq < rhs.seq; });

    //Sorting read seeds
    for(int i = 0; i < readSeeds.size(); i++){
        sort(readSeeds[i].begin(), readSeeds[i].end(), []( const seed& lhs, const seed& rhs ) { return lhs.seq < rhs.seq; });
    }
    
    //Finding syncmer matches
    for(int r = 0; r < readSequences.size() ; r++){

        int refSeedIndex = 0;
        int readSeedIndex = 0;

        std::vector<seed> matchingSeeds;
        while(readSeedIndex < readSeeds[r].size() && refSeedIndex < refSeeds.size()) {
            if (readSeeds[r][readSeedIndex].seq < refSeeds[refSeedIndex].seq) {
                readSeedIndex++;
            } else if (readSeeds[r][readSeedIndex].seq > refSeeds[refSeedIndex].seq) {
                refSeedIndex++;
            } else {

                matchingSeeds.push_back(readSeeds[r][readSeedIndex]);
                readSeedIndex++;
                refSeedIndex++;
            }
        }
        readSeeds[r] = matchingSeeds;
    }
    
    //Preparing C structures for minimap
    const char *reference = bestMatchSequence.c_str();
    int n_reads = readSequences.size();
    char **read_strings = (char **)malloc(n_reads*sizeof(char *));
    char **qual_strings = (char **)malloc(n_reads*sizeof(char *));
    char **read_names = (char **)malloc(n_reads*sizeof(char *));

    int *r_lens         = (int *)malloc(n_reads*sizeof(int));
    int *seed_counts    = (int *)malloc(n_reads*sizeof(int));

    uint8_t **reversed  = (uint8_t **)malloc(n_reads*sizeof(uint8_t *));
    int **ref_positions = (int **)malloc(n_reads*sizeof(int *));
    int **qry_positions = (int **)malloc(n_reads*sizeof(int *));
    
    for(int i = 0; i < n_reads; i++) {
        int n_seeds = readSeeds[i].size();

        seed_counts[i] = n_seeds;
        read_strings[i] = (char *) readSequences[i].c_str();
        qual_strings[i] = (char *) readQuals[i].c_str();
        read_names[i] =   (char *) readNames[i].c_str();

        r_lens[i] = readSequences[i].length();

        uint8_t *reversed_array = (uint8_t *)malloc(n_seeds*sizeof(uint8_t));
        int *ref_pos_array = (int *)malloc(n_seeds*sizeof(int));
        int *qry_pos_array = (int *)malloc(n_seeds*sizeof(int));


        for(int j = 0; j < n_seeds; j++){
            reversed_array[j] = readSeeds[i][j].reversed;
            qry_pos_array[j] = readSeeds[i][j].pos;
            ref_pos_array[j] = seedToRefPositions[readSeeds[i][j].seq][0] + k - 1;
        }

        reversed[i]      = reversed_array;
        ref_positions[i] = ref_pos_array;
        qry_positions[i] = qry_pos_array;
    }
    
    std::string sam_header = "@SQ\tSN:reference\tLN:";
    sam_header += std::to_string(bestMatchSequence.length());

    char *sam_alignments[n_reads]; //constituants must be freed

    align_reads(reference, n_reads, read_strings, qual_strings, read_names, r_lens, seed_counts, reversed, ref_positions, qry_positions, sam_alignments, k);
    
    /* Sorting the alignments by their reference position */
    std::pair<int, char*> sam_lines[n_reads];
    for(int i = 0; i < n_reads; i++) {
        sam_lines[i] = std::make_pair(r_lens[i],sam_alignments[i]);
    }
    sort(sam_lines, sam_lines + n_reads, [](const std::pair<int, char*>& a, const std::pair<int, char*>& b) {
        return a.first < b.first;
    });

    for(int i = 0; i < n_reads; i++) {
        sam_alignments[i] = sam_lines[i].second;
        if(!sam_alignments[i]){
            n_reads = i; //Some reads failed
        }
    }

    //Print out sam
    if(samFileName.size() > 0){
        std::ofstream outFile{samFileName};
        if (outFile.is_open()) {
            outFile << sam_header << std::endl;
            for(int i = 0; i < n_reads; i++) {
                if(sam_alignments[i]) {
                    outFile << sam_alignments[i] << std::endl;
                }
            }
            std::cout << "Wrote sam data to " << samFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << samFileName << std::endl;
        }
    }

    for(int i = 0; i < n_reads; i++) {
        free(reversed[i]);
        free(ref_positions[i]);
        free(qry_positions[i]);
    }
    free(qry_positions);
    free(ref_positions);
    free(reversed);
    free(seed_counts);
    free(read_strings);
    free(qual_strings);
    free(read_names);
    free(r_lens);

    // Convert to BAM

    sam_hdr_t *header = NULL;

    // Parse SAM header
    header = sam_hdr_parse(sam_header.length(), sam_header.c_str());

    htsFile *bam_file = NULL;
    if (bamFileName.size() > 0) {
        bam_file = hts_open(bamFileName.c_str(), "wb");
        if (!bam_file) {
            fprintf(stderr, "Error: Failed to open output BAM file.\n");
            hts_close(bam_file);
        }
        // Write BAM header
        else if (sam_hdr_write(bam_file, header) < 0) {
            fprintf(stderr, "Error: Failed to write BAM header.\n");
        }
    }

    // Prepare list of bam1_t
    bam1_t **bam_records = (bam1_t **)malloc(sizeof(bam1_t *) * n_reads);
    for (int i = 0; i < n_reads; i++) {
        if(sam_alignments[i]){
            bam_records[i] = bam_init1();
            kstring_t line = KS_INITIALIZE;
            kputs(sam_alignments[i], &line);
            sam_parse1(&line, header, bam_records[i]);
            // Write to bam file
            if (bam_file && bam_write1(bam_file->fp.bgzf, bam_records[i]) < 0) {
                fprintf(stderr, "Error: Failed to write BAM record.\n");
                bam_hdr_destroy(header);
                hts_close(bam_file);
            }
        }
    }
    if(bam_file){
        std::cerr << "Wrote bam files to " << bamFileName << "\n";
    }
    // Converted to Bam
    hts_close(bam_file);
    for(int i = 0; i < n_reads; i++) {
        free(sam_alignments[i]);
    }

    // Convert to mpileup
    // Need to copy reference string because we need a char* rather than const char*
    char* ref_string = new char[bestMatchSequence.length() + 1];
    std::strcpy(ref_string, bestMatchSequence.c_str());

    kstring_t mplp_string = KS_INITIALIZE;
    bam_and_ref_to_mplp(header, bam_records, n_reads, ref_string, bestMatchSequence.size(), &mplp_string);
    
    //Print out mpileup
    if(mpileupFileName.size() > 0){
        std::ofstream outFile{mpileupFileName};
        if (outFile.is_open()) {
            outFile << mplp_string.s;
            std::cout << "Wrote mpileup data to " << mpileupFileName << std::endl;
        } else {
            std::cerr << "Error: failed to write to file " << mpileupFileName << std::endl;
        }
    }

    //Convert to VCF
    //Get mutation matrix.
    mutationMatrices mutMat = mutationMatrices();
    fillMutationMatrices(mutMat, T);

    // Convert c string of mpileup to ifstream
    std::istringstream mpileipStream(mplp_string.s);

    std::ofstream vcfOutFile;
    if(vcfFileName.size() > 0) {
        vcfOutFile.open(vcfFileName);
        if (vcfOutFile.is_open()) {

            genotype::printSamplePlacementVCF(mpileipStream, mutMat, true, 0, vcfOutFile);

            std::cout << "Wrote vcf data to " << vcfFileName << std::endl;
        }else{

            std::cerr << "Error: failed to write to file " << vcfFileName << std::endl;
        }
    }

    free(mplp_string.s);
    
    for(int i = 0; i < n_reads; i++){
        bam_destroy1(bam_records[i]);
    }

    free(bam_records);
}

//Stores sam strings in sam_alignments, and stores reference positions of alignments in r_lens
void align::align_reads(const char *reference, int n_reads, char **reads, char **quality, char **read_names, int *r_lens, int *seed_counts, uint8_t **reversed, int **ref_positions, int **qry_positions, char** sam_alignments, int syncmer_k) {
    mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	int n_threads = 1;

	mm_verbose = 2;             // disable message output to stderr
	mm_set_opt(0, &iopt, &mopt);
	mopt.flag |= MM_F_CIGAR;    // perform alignment

	
	iopt.k = syncmer_k; //setting seed length

	mm_idx_t *mi = mm_idx_str(10, iopt.k, 0, 14, 1, &reference, NULL); //Read reference into structure

	mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence;
	mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread


	
	for(int k = 0; k < n_reads; k++) {
		
		mm_reg1_t *reg;
		int j, i, n_reg;

		align_read_given_seeds(mi,r_lens[k],reads[k], &n_reg, &reg, tbuf, &mopt, iopt.k, seed_counts[k], ref_positions[k], qry_positions[k], reversed[k]);
		
		kstring_t sam = {0,0,0};
		mm_bseq1_t t;
		t.l_seq = r_lens[k];
		t.seq = reads[k];
		t.name = read_names[k];
		t.qual = quality[k];
		const int n_regss[1] = {n_reg};
		const mm_reg1_t *regss = {reg};
		mi->seq[0].name = "reference";

		if(n_reg == 0 || reg->score <= 0 || reg->score > r_lens[k]) { //Maybe remove n_reg==0 check
			sam_alignments[k] = NULL;
			r_lens[k] = -1;
		} else {
			mm_write_sam3( &sam, mi, &t, 0, 0, 1, n_regss, &regss, NULL, 0, 0);
			r_lens[k] = reg->rs+1; //Length of sam line string
			sam_alignments[k] = sam.s;
		}
		


		for (j = 0; j < n_reg; ++j) {
			mm_reg1_t *r = &reg[j];
			assert(r->p); // with MM_F_CIGAR, this should not be NULL

			free(r->p);
		}
		free(reg);

	}
	mm_tbuf_destroy(tbuf);
	mm_idx_destroy(mi);
}