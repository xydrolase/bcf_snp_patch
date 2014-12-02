#include <iostream>
#include <fstream>
#include "bcf_snp_patch.h"

#include <boost/program_options.hpp>

namespace progopt = boost::program_options;

extern "C" {
#include "htslib/vcf.h" // variant calling
#include "htslib/sam.h" // alignment
#include "htslib/bgzf.h"  // BAM
}

void dump_bcf_header(bcf_hdr_t *hdr);
void dump_bcf_hrecs(bcf_hrec_t **hrecs, int nhrec);
void load_snps_bcf(const char *bcf_fpath, snp_map_t &smap, bsp_options &opt);
void output_error_file(const char *bam_fpath, snp_map_t &smap, 
		std::unordered_set<std::string> &fltr, bsp_options &opt);
static inline void print_seq_error(std::ostream &ss, uint32_t is_rev, std::string rid, 
		int rlen, int epos, char ref, char err);

// convert 4-bit repr. of a nucleotide to the "true" literal nucleotide ACTG.
static inline char bam_4bit2char(uint8_t b4) 
{
	// use GCC/GNU builtin to scan the number of trailing zeros.
	// i.e. 1 -> 0, 2 -> 1, 4 -> 2, 8 -> 3
	int idx = __builtin_ctz(b4); 
	return sam_nucleotides[idx & 3];
}

static void parse_options(int argc, char *argv[], bsp_options &opt)
{
	progopt::options_description opts_desc("Allowed options:");
	opts_desc.add_options()
	("help,h", "print usage information")
	("output,o", progopt::value(&opt.out_file), 
	 "redirect output to the file <arg> (instead of stdout)")
	("min_depth,d", progopt::value(&opt.min_read_depth)->default_value(0), 
	 "maximum read depth for a SNP to be considered valid")
	("max_depth,D", progopt::value(&opt.max_read_depth)->default_value(65536), 
	 "maximum read depth for a SNP to be considered valid");

	progopt::options_description pdesc("Positional options:");
	pdesc.add_options()
	("bcf_file", progopt::value(&opt.bcf_file)->required(), 
	 "BCF formatted variant calling file")
	("bam_file", progopt::value(&opt.bam_file)->required(),
	 "BAM formatted alignment results")
	("rid_file", progopt::value(&opt.rid_file)->required(),
	 "a list of read ids serving as read filters");

	progopt::options_description all_desc("All options:");
	all_desc.add(opts_desc).add(pdesc);

	progopt::positional_options_description pos_desc;
	pos_desc.add("bcf_file", 1);
	pos_desc.add("bam_file", 1);
	pos_desc.add("rid_file", 1);

	progopt::variables_map vm;

	// usage information
	std::stringstream usage;
	usage << std::endl<< "Usage: " << argv[0] 
		<< " [options] <bcf_file> <bam_file> <rid_file>"
		<< std::endl << std::endl << opts_desc << std::endl;

	try {
		progopt::store(progopt::command_line_parser(argc, argv)
				.options(all_desc)
				.positional(pos_desc).run(), vm);
		progopt::notify(vm);
	}
	catch(std::exception& e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::cerr << usage.str();
		exit(1);
	}

	if (vm.count("help")) {
		std::cout << usage.str();
		std::cout << opts_desc << std::endl;
		exit(0);
	}
}

int main(int argc, char *argv[])
{
	// bcf_snp_path <bcf_file> <bam_file> <rid_file>
	
	/* PIPELINE:
	 *
	 *   a) Load BCF/VCF file
	 *   b) Load read id filters 
	 *   c) Load BAM/SAM file and apply SNPs
	 */
	
#if 0
	/* ---- FROM vcf.h ---- */
	/* === Dictionary ===
	   The header keeps three dictonaries. The first keeps IDs in the
	   "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
	   in the "contig" lines and the last keeps the sample names. bcf_hdr_t::dict[]
	   is the actual hash table, which is opaque to the end users. In the hash
	   table, the key is the ID or sample name as a C string and the value is a
	   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
	   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
	   size of the hash table or, equivalently, the length of the id[] arrays.
	*/

	/* Useful type definitions */

	typedef struct {
		int type;       // One of the BCF_HL_* type
		char *key;      // The part before '=', i.e. FILTER/INFO/FORMAT/contig/fileformat etc.
		char *value;    // Set only for generic lines, NULL for FILTER/INFO, etc.
		int nkeys;              // Number of structured fields
		char **keys, **vals;    // The key=value pairs
	} bcf_hrec_t;
#endif

	bsp_options opt;
	parse_options(argc, argv, opt);

	std::cerr << "Processing VCF/BCF formatted SNPs..." << std::endl;
	snp_map_t bcf_map;
	load_snps_bcf(opt.bcf_file.c_str(), bcf_map, opt);
	std::cerr << bcf_map.size() << " SNPs loaded." << std::endl;

	std::cerr << "Loading filtering read identifiers..." << std::endl;
	// set of read identifiers (extended) as filter
	std::unordered_set<std::string> rid_fltr_set;
	std::ifstream rid_ifs(opt.rid_file.c_str());
	std::string rid;
	while (std::getline(rid_ifs, rid)) {
		if (rid.back() == ' ')
			rid.pop_back(); // remove the trailing space

		rid_fltr_set.insert(rid);
	}
	std::cerr << rid_fltr_set.size() << " read identifiers loaded." << std::endl;
	
	std::cerr << "Processing SAM/BAM formatted alignments..." << std::endl;
	output_error_file(opt.bam_file.c_str(), bcf_map, rid_fltr_set, opt);
	std::cerr << "Alignments processed." << std::endl;
}


void print_seq_error(std::ostream &ss, uint32_t is_rev, std::string rid, 
		int rlen, int epos, char ref, char err)
{
	if (is_rev) {
		epos = rlen - epos;
		ref = pmr_complements[base_to_bits(ref)];
		err = pmr_complements[base_to_bits(err)];
	}
	else {
		++epos; // convert from 0-base to 1-based coordinate
	}

	ss << rid << " " << epos << " " << ref << " " << err << std::endl;
}

void output_error_file(const char *bam_fpath, snp_map_t &smap, 
		std::unordered_set<std::string> &fltr, bsp_options &opt)
{
	const char aux_md[2] = {'M', 'D'};

	std::ofstream ofs_error;
	std::ostream &ostrm_error = opt.out_file.empty() ? std::cout : ofs_error;
	if (!opt.out_file.empty()) {
		ofs_error.open(opt.out_file.c_str(), std::ofstream::out | std::ofstream::trunc);
	}

	BGZF *bam_fp = bgzf_open(bam_fpath, "r");
	bam_hdr_t *hdr = bam_hdr_read(bam_fp);

	std::cerr << "BAM header loaded. # targets: " << hdr->n_targets << std::endl;

	bam1_t *aln = bam_init1();

	int nbytes = 0;
	while ((nbytes = bam_read1(bam_fp, aln)) >= 0) {
		uint32_t flag = aln->core.flag;

		// check if unmapped, if so, discard immediately
		if (flag & BAM_FUNMAP) continue;

		// expand read identifier according to the flag BAM_FREAD1 / BAM_FREAD2
		std::stringstream erid;
		erid << bam_get_qname(aln) << READ_EXPAND_SEP << 
			((flag & BAM_FREAD1) ? "1" : "2");

		// check if the current read will be selected according to the filter
		auto hit = fltr.find(erid.str());
		if (hit == fltr.end()) continue;

#if DEBUG
		std::cerr << erid.str() << std::endl;
		getchar();
#endif

		// position on the template (reference chromosome)
		int tpos = aln->core.pos;
		// position for the putative error (mismatch)
		uint64_t mism_pos = (uint64_t) tpos;

		// CIGAR is not quite useful in this program (other than computing the 
		// read length) , but might be useful if we'd implement a program 
		// manipulating BAM/SAM file in C++.
		uint32_t *cigar = bam_get_cigar(aln);
		int rlen = bam_cigar2rlen(aln->core.n_cigar, cigar);

		uint8_t *aux = bam_aux_get(aln, aux_md);
		char *MD = bam_aux2Z(aux);

		// pointer to the 4-bit repr. of the read sequence
		const uint8_t *rseq = bam_get_seq(aln);

		// process MD aux tag
		uint32_t epos = 0;
		uint32_t radix = 0;
		for (char *p = MD; *p != '\0'; ++p) {
			char pc = *p;
			if (isalpha(pc)) {
				mism_pos += radix;
				epos += radix;

				//std::cout << radix << ", " << mism_pos << ", " << pc << std::endl;
				uint64_t _loci = aln->core.tid | (mism_pos << 32);
				// check for true SNPs
				auto snp = smap.find(_loci);

				// observed nucleotide on the read:
				// NOTE that bam_seqi returns a 4-bit representation of
				// nucleotide: A = 1, C = 2, G = 4, T = 8, N = 15.
				// We need to convert the 4-bit repr. back to char
				char obs_epos_base = bam_4bit2char(bam_seqi(rseq, epos));

#ifdef DEBUG
				std::cout << "Base @ " << epos << " == " << (int) (bam_seqi(rseq, epos)) <<
					"(4-bit repr.) -> " << obs_epos_base << std::endl;
#endif

				if (snp == smap.end()) {
					// output the error directly
					print_seq_error(ostrm_error, flag & BAM_FREVERSE, erid.str(), 
							rlen, epos, pc, obs_epos_base);
				}
				else {
					// check if this is a valid allele
					uint32_t _allele_flag = (snp->second >> 2);
					if ((_allele_flag & (1 << (base_to_bits(obs_epos_base)))) == 0) {
						// not true SNP, likely a sequencing error
						print_seq_error(ostrm_error, flag & BAM_FREVERSE, erid.str(), 
								rlen, epos, pc, obs_epos_base);
					}
				}

				// acount for this particular nucleotide.
				++mism_pos; 
				++epos;

				radix = 0; // reset radix
			}
			else if (isdigit(pc)) {
				radix = radix * 10 + (pc - '0');
			}
			else {
				// likely a ^ char indicating insertion errors
				// this shouldn't exist if the filter REALLY works.
				std::cerr << "Invalid character: " << pc << std::endl;
				exit(1);
			}
		}
	}

	if (!opt.out_file.empty()) {
		ofs_error.close();
	}
}

void load_snps_bcf(const char *bcf_fpath, snp_map_t &smap, bsp_options &opt)
{
	vcfFile *bcf_fp = bcf_open(bcf_fpath, "r");
	bcf_hdr_t *hdr = bcf_hdr_read(bcf_fp);

#ifdef DEBUG
	dump_bcf_header(hdr);
	dump_bcf_hrecs(hdr->hrec, hdr->nhrec);
#endif

	// extract the INFO index corresponding to "INDEL"
	// BCF_DT_ID corresponds to the first of the three dictionaries encoded by
	// BCF format. 
	int INFO_INDEL_IDX = bcf_hdr_id2int(hdr, BCF_DT_ID, "INDEL");
	if (INFO_INDEL_IDX == -1) {
		std::cerr << "ERROR: Cannot extract index for 'INDEL' tag." << std::endl;
		exit(10);
	}
	int INFO_DP_IDX = bcf_hdr_id2int(hdr, BCF_DT_ID, "DP");
	if (INFO_DP_IDX == -1) {
		std::cerr << "ERROR: Cannot extract index for 'DP' tag." << std::endl;
		exit(10);
	}

	bcf1_t *vcf_line = NULL;
	vcf_line = bcf_init();

#if 0
	/* test for reading a few BCF records */
	for (int i = 0; i < 1000; ++i) {
		bcf_read(bcf_fp, hdr, vcf_line);
		printf("SEQNAME: %s\n", bcf_seqname(hdr, vcf_line));
		printf("[%6d]\n\tREF : %d\n\tPOS : %d\n\tRLEN: %d\n"
				">>> Summary:\n"
				"\tN_SMP: %u\n\tN_ALLELE: %u\n"
				">>> VCF line:\n"
				"\tIndiv : %s\n"
				"\tShared: %s\n",
				i,
				vcf_line->rid, vcf_line->pos, vcf_line->rlen,
				vcf_line->n_sample, vcf_line->n_allele,
				(vcf_line->indiv.s), (vcf_line->shared.s));
		getchar();
	}
#endif

	int retval = 0;
	while ((retval = bcf_read(bcf_fp, hdr, vcf_line)) == 0) {
		/* BCF format packs data into vcf_line->indiv.s and vcf_line->shared.s, 
		 * both are binary-packed "strings". 
		 * Calling bcf_read does not unpack such information (so-called lazy
		 * unpacking).
		 *
		 * BCF_UN_INFO unpacks to the level of "INFO", which encodes
		 * "additional information" like SB (strand bias), DP (depth) etc.
		 * See VCF File format specifications for more details.
		 */
		bcf_unpack(vcf_line, BCF_UN_INFO); 

		/*
		 * vcf_line->n_info is the number of "INFO" fields
		 * To retrieve an "INFO" field, use bcf_get_info(hdr, vcf_line, KEY)
		 * where KEY is a (const char *) typed string.
		 *
		 * Alternatively, one can use bcf_get_info_id(vcf_line, ID)
		 * where ID is the ID associated with the KEY.
		 *
		 * To get the corresponding ID given KEY, 
		 * use bcf_hdr_id2int() for KEY-->ID 
		 * and bcf_hdr_int2id() for ID-->KEY
		 */

		// skip INDELs
		if (bcf_get_info_id(vcf_line, INFO_INDEL_IDX) != NULL) continue;

		// check read depth
		bcf_info_t *bi = bcf_get_info_id(vcf_line, INFO_DP_IDX);
		if (bi == NULL) continue; // no read depth information, skip
		register int vcf_dp = bi->v1.i;
		if (vcf_dp < opt.min_read_depth || vcf_dp > opt.max_read_depth) continue;

		// cache all the SNPs
		bcf_dec_t *d = &(vcf_line->d);
		// first allele is the REF(erence)
		char ref = d->allele[0][0];
		uint32_t _allele_flag = 0;

		int n_alt = vcf_line->n_allele - 1;
		for (int j = 0; j < n_alt; ++j) {
			_allele_flag |= (1 << base_to_bits(d->allele[j+1][0]));
		}

		snp_loci_t _loci = (uint64_t) vcf_line->rid | (
				(uint64_t) vcf_line->pos << 32);
		uint32_t _snp = base_to_bits(ref) | (_allele_flag << 2);

		//smap.insert(std::make_pair<snp_loci_t, uint32_t> (_loci, _snp)); 
		smap.insert(std::pair<snp_loci_t, uint32_t> (_loci, _snp)); 
	}

	if (vcf_line) bcf_destroy(vcf_line);

	bcf_hdr_destroy(hdr);
	bcf_close(bcf_fp);
}

void dump_bcf_header(bcf_hdr_t *hdr)
{
	printf("Numer of:\n\tIDs: %d\n\tSequences: %d\n\tSamples: %d\n"
			"\tHeader records: %d\n"
			"Pointer to hrecs: 0x%x\n",
			hdr->n[0], hdr->n[1], hdr->n[2], hdr->nhrec, hdr->hrec);
}

void dump_bcf_hrecs(bcf_hrec_t **hrecs, int nhrec)
{
	const char *hrec_types[] = {
		"BCF_HL_FLT",
		"BCF_HL_INFO",
		"BCF_HL_FMT", 
		"BCF_HL_CTG", 
		"BCF_HL_STR",
		"BCF_HL_GEN", 
	};

	for (int i = 0; i < nhrec; ++i) {
		bcf_hrec_t *hr = hrecs[i];
		printf("[%4d] %s <%d>\n\t%s = %s\n",
				i, hrec_types[hr->type], 
				hr->nkeys, hr->key, hr->value);

		if (hr->nkeys > 0) {
			printf("\t>>> structured fields:\n");
			for (int j = 0; j < hr->nkeys; ++j) {
				printf("\t\t%s = %s\n",
						hr->keys[j], hr->vals[j]);
			}
		}
		getchar();
	}
}
