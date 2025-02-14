/* bpred.h - branch predictor interfaces */

/* SimpleScalar(TM) Tool Suite
 * Copyright (C) 1994-2003 by Todd M. Austin, Ph.D. and SimpleScalar, LLC.
 * All Rights Reserved. 
 * 
 * THIS IS A LEGAL DOCUMENT, BY USING SIMPLESCALAR,
 * YOU ARE AGREEING TO THESE TERMS AND CONDITIONS.
 * 
 * No portion of this work may be used by any commercial entity, or for any
 * commercial purpose, without the prior, written permission of SimpleScalar,
 * LLC (info@simplescalar.com). Nonprofit and noncommercial use is permitted
 * as described below.
 * 
 * 1. SimpleScalar is provided AS IS, with no warranty of any kind, express
 * or implied. The user of the program accepts full responsibility for the
 * application of the program and the use of any results.
 * 
 * 2. Nonprofit and noncommercial use is encouraged. SimpleScalar may be
 * downloaded, compiled, executed, copied, and modified solely for nonprofit,
 * educational, noncommercial research, and noncommercial scholarship
 * purposes provided that this notice in its entirety accompanies all copies.
 * Copies of the modified software can be delivered to persons who use it
 * solely for nonprofit, educational, noncommercial research, and
 * noncommercial scholarship purposes provided that this notice in its
 * entirety accompanies all copies.
 * 
 * 3. ALL COMMERCIAL USE, AND ALL USE BY FOR PROFIT ENTITIES, IS EXPRESSLY
 * PROHIBITED WITHOUT A LICENSE FROM SIMPLESCALAR, LLC (info@simplescalar.com).
 * 
 * 4. No nonprofit user may place any restrictions on the use of this software,
 * including as modified by the user, by any other authorized user.
 * 
 * 5. Noncommercial and nonprofit users may distribute copies of SimpleScalar
 * in compiled or executable form as set forth in Section 2, provided that
 * either: (A) it is accompanied by the corresponding machine-readable source
 * code, or (B) it is accompanied by a written offer, with no time limit, to
 * give anyone a machine-readable copy of the corresponding source code in
 * return for reimbursement of the cost of distribution. This written offer
 * must permit verbatim duplication by anyone, or (C) it is distributed by
 * someone who received only the executable form, and is accompanied by a
 * copy of the written offer of source code.
 * 
 * 6. SimpleScalar was developed by Todd M. Austin, Ph.D. The tool suite is
 * currently maintained by SimpleScalar LLC (info@simplescalar.com). US Mail:
 * 2395 Timbercrest Court, Ann Arbor, MI 48105.
 * 
 * Copyright (C) 1994-2003 by Todd M. Austin, Ph.D. and SimpleScalar, LLC.
 */


#ifndef BPRED_H
#define BPRED_H

#define dassert(a) assert(a)

#include <stdio.h>

#include "host.h"
#include "misc.h"
#include "machine.h"
#include "stats.h"

/*
 * This module implements a number of branch predictor mechanisms.  The
 * following predictors are supported:
 *
 *	BPred2Level:  two level adaptive branch predictor
 *
 *		It can simulate many prediction mechanisms that have up to
 *		two levels of tables. Parameters are:
 *		     N   # entries in first level (# of shift register(s))
 *		     W   width of shift register(s)
 *		     M   # entries in 2nd level (# of counters, or other FSM)
 *		One BTB entry per level-2 counter.
 *
 *		Configurations:   N, W, M
 *
 *		    counter based: 1, 0, M
 *
 *		    GAg          : 1, W, 2^W
 *		    GAp          : 1, W, M (M > 2^W)
 *		    PAg          : N, W, 2^W
 *		    PAp          : N, W, M (M == 2^(N+W))
 *
 *	BPred2bit:  a simple direct mapped bimodal predictor
 *
 *		This predictor has a table of two bit saturating counters.
 *		Where counter states 0 & 1 are predict not taken and
 *		counter states 2 & 3 are predict taken, the per-branch counters
 *		are incremented on taken branches and decremented on
 *		no taken branches.  One BTB entry per counter.
 *
 *	BPredTaken:  static predict branch taken
 *
 *	BPredNotTaken:  static predict branch not taken
 *
 */

/* branch predictor types */
enum bpred_class {
  BPredComb,                    /* combined predictor (McFarling) */
  BPred2Level,			/* 2-level correlating pred w/2-bit counters */
  BPredTSBP,			/* Temporal Stream Branch Predictor */
  BPredCHBP,			/* Correctness History Branch Predictor (improved TSBP) */
  BPredOB,			/* 2-Lev + OB */
  BPredOHT,			/* 2-Lev + OHT  */
  BPredMBP,			/* Mississippi Branch Predictor (2-lev + OB + OHT)  */
  BPredTSCL,			/* TAGE-SC-L Branch Predictor */
  BPredLLBP,			/* Last-Level Branch Predictor (TSCL + prefetched LLBT) */
  BPred2bit,			/* 2-bit saturating cntr pred (dir mapped) */
  BPredTaken,			/* static predict taken */
  BPredNotTaken,		/* static predict not taken */
  BPred_NUM
};

/* an entry in a BTB */
struct bpred_btb_ent_t {
  md_addr_t addr;		/* address of branch being tracked */
  enum md_opcode op;		/* opcode of branch corresp. to addr */
  md_addr_t target;		/* last destination of branch when taken */
  struct bpred_btb_ent_t *prev, *next; /* lru chaining pointers */
};

/* Future History Buffer (FHB) for each reversible 2-level implementation */
struct bpred_fhb_t {
	int size;	/* number of FHB entries should be L1 width + 1 */
	int top;	/* Top of FHB window */
	int bot;	/* Bottom of FHB Window */
	int fv_out;	/* The buffered output FV from the FHB */
	int rv_out;	/* The buffered output RV from the FHB */
	int o_out;	/* The buffered output O from the FHB */
	int correct_out;	/* The buffered output Correct from the FHB */
	int und_out;	/* The buffered output UND from the FHB */
	int thu_out;	/* The buffered output THU from the FHB */
	md_addr_t addr_out; 	/* The buffered output addr from the FHB */
	md_addr_t key_out; 	/* The buffered output key from the FHB */
	int *fv;    /* fwd valid bits for prediction use in fwd mode */
	int *rv;    /* rev valid bits for prediction use in rev mode */
	int *o;     /* outcome bits for prediction use in fwd & rev modes */
	int *correct;	/* track correctness */
	int *und;		/* Under threshold bit */
	int *thu;		/* threshold update bit */
	md_addr_t *addr;
	md_addr_t *key; // Key component for use with FRMTs 
};

/* direction predictor def */
struct bpred_dir_t {
	enum bpred_class class;	/* type of predictor */
	union {
		struct {
			unsigned int size;	/* number of entries in direct-mapped table */
			unsigned char *table;	/* prediction state table */
		} bimod;
		struct {
			int l1size;		/* level-1 size, number of history regs */
			int l2size;		/* level-2 size, number of pred states */
			int shift_width;		/* amount of history in level-1 shift regs */
			int threshold;			/*	SC/LOOP threshold */
			int xor;			/* history xor address flag */
			int *context;		/* LLBT context bits */
			int *shiftregs;		/* level-1 history table */
			int *hist_length;	/* history lengths for LLBT table */
			unsigned char *l2table;	/* level-2 prediction state table */
			unsigned char *use;	/* useful counter for TAGE tables */
			unsigned char *iter_c;	/* current iter counter for LOOP tables */
			unsigned char *iter_p;	/* past iter counter for LOOP tables */
			md_addr_t *tag;		/* tags TAGE tables */
		} two;
	} config;
};

typedef unsigned long long ts_key_t; /*TS Key data type*/

/* temporal stream predictor def */
struct bpred_ts_t {
  enum bpred_class class;				/* type of predictor */
  struct {
    unsigned int head_table_size;		/* head table size, number of head table entries */
    unsigned int head_table_width;		/* head table width, correlated w/ correctness buffer header size */
    unsigned int correctness_width;		/* correctness circular buffer width, 2^head_table_width */
    bool_t *correctness_buffer;			/* correctness circular buffer */
    unsigned int *head_table;			/* head table containing pointers to circular buffer */
	unsigned int head;					/* head of the circular buffer */
	unsigned int tail;					/* tail of the circular buffer */
	bool_t replay;						/* replay flag */
	bool_t enabled;						/* TSBP Enabled flag */
  } ts;
};

/* temporal stream predictor def */
struct bpred_chbp_t {
  enum bpred_class class;				/* type of predictor */
  struct {
    unsigned int cht_size;				/* correctness history table size, number of table entries */
	bool_t enabled;						/* CHBP Enabled flag */
	md_addr_t *cht_spc;					/* Correctness history table source pc bits*/
	bool_t *cht_replay;					/* Correctness history table replay bits*/
	bool_t *cht_correct;				/* Correctness history table correct bits*/
	bool_t *cht_valid;					/* Correctness history table valid bits*/
	md_addr_t *cht_dpc;					/* Correctness history table destination pc bits*/
  } chbp;
};

/* Outcome History Table (OHT) def */
struct bpred_oht_t {
  enum bpred_class class;				/* type of predictor */
  struct {
	unsigned int size;				/* outcome history table size, number of table entries */
	int *oc;							/* outcome bits */
	int *fwd_valid;							/* FWD valid bits; FWD and REV are needed for FRMT implementation */
	int *rev_valid;							/* REV valid bits */
  } oht;
};

/* Create FWD-to-REV Mapping Table struct */
struct bpred_frmt_t {
	enum bpred_class class;	/* FRMT setup type */
	unsigned int budget;	/* FRMT budget size */
	md_addr_t *addr;		/* FRMT addr entries */
	md_addr_t *tag;		/* FRMT addr entries */
	md_addr_t *context;			/* FRMT addr entries */
};

// Making dirpred a struct accessible outside of bpred_t
struct bpred_dirpred_t {
    struct bpred_dir_t *bimod;	  /* first direction predictor */
    struct bpred_dir_t *twolev;	  /* second direction predictor */
	struct bpred_ts_t  *tsbp;	  /* temporal stream */
	struct bpred_chbp_t *chbp;	  /* correction history branch predictor */
    struct bpred_dir_t *meta;	  /* meta predictor */
	struct bpred_oht_t *oht;	  /* Outcome History Table (OHT) for use with reversible 2-level/MBP */
	struct bpred_frmt_t *frmt;	  /* FWD-to-REV Mapping Tables */
};

/* branch predictor def */
struct bpred_t {
  enum bpred_class class;	/* type of predictor */
  
  int tage_depth;	/* Depth of tables for TSCL/LLBP */
  int num_fhbs;		/* Number of FHBs in use */
  
  /* FWD Structs */
  struct bpred_dirpred_t fwd_dirpred;
  struct bpred_dirpred_t* fwd_tage_dirpred;	//	FWD TAGE Dirpreds
  struct bpred_dirpred_t* fwd_sc_dirpred;	//	FWD SC Dirpreds
  struct bpred_dirpred_t fwd_loop_dirpred;	//	FWD LOOP Dirpreds
  struct bpred_dirpred_t* fwd_pb_dirpred;	//	FWD LLBP Pattern Buffer (PB) Dirpreds
  struct bpred_dirpred_t* fwd_llbt_dirpred;	//	FWD LLBP LLBT Dirpreds
  
  /* An FHB is needed for every set of reversible 2 level preds */
  struct bpred_fhb_t* fhb; 	// Can have an array of FHBs of various lengths 
  
  /* Reversible Structs */
  struct bpred_dirpred_t rev_dirpred;
  struct bpred_dirpred_t* rev_tage_dirpred;	//	REV TAGE Dirpreds
  struct bpred_dirpred_t* rev_sc_dirpred;	//	REV SC Dirpreds
  struct bpred_dirpred_t rev_loop_dirpred;	//	REV LOOP Dirpreds
  struct bpred_dirpred_t* rev_pb_dirpred;	//	REV LLBP Pattern Buffer (PB) Dirpreds
  struct bpred_dirpred_t* rev_llbt_dirpred;	//	REV LLBT Dirpreds
  
  /* Rolling Context Register (RCR) for LLBP */
  struct {
	int size;	/* number of RCR entries */
	int top;	/* Top of RCR */
	int win_top;	/* Top of CCID window */
	int fpre_top;	/* Bottom of CCID Window */
	int rpre_bot;	/* Bottom of CCID Window */
	int win_bot;	/* Bottom of CCID Window */
	int bot;	/* Bottom of RCR */
	int fccid;	/* FWD pre CCID */
	int rccid;	/* REV pre CCID */
	int ccid;		/* RCR Current Context ID (CCID) */
	int *fv;    /* fwd valid bits for update use in rev mode */
	int *rv;    /* rev valid bits for update use in fwd mode */
	md_addr_t *addr;	/* current RCR address window */
  } rcr;
  
  /* Future Context Buffer (FCB) for reversible LLBP */
  struct {
	int size;	/* number of FHB entries should be L1 width + 1 */
	int top;	/* Top of FHB window */
	int bot;	/* Bottom of FHB Window */
	int *fv;    /* fwd valid bits for prediction use in fwd mode from FHB */
	int *rv;    /* rev valid bits for prediction use in rev mode from FHB */
	int *o;     /* outcome bits for prediction use in fwd & rev modes from FHB */
	md_addr_t *context;	/* context from RCR */
	md_addr_t *tag;	/* tag from Addr and GH after leaving FHB */
	int *tv;     /* tag valid bits for determining where to store the tag bits after acquired */
  } fcb;
  
  /* Outcome Buffer (OB) for use with reversible 2-level/MBP */
  struct {
	int width;	/* number of OB entries; independent of any L table dimenssions */
	int end;	/* End of buffer window at width - 1 */
	int beg;	/* Beginning of buffer window at 0 */
	int *fv;    /* fwd valid bits for prediction use in fwd mode */
	int *rv;    /* rev valid bits for prediction use in rev mode */
	int *oc;     /* outcome bits for prediction use in fwd & rev modes */
  } ob;

  struct {
    int sets;			/* num BTB sets */
    int assoc;			/* BTB associativity */
    struct bpred_btb_ent_t *btb_data; /* BTB addr-prediction table */
  } btb;

  struct {
    int size;			/* return-address stack size */
    int tos;			/* top-of-stack */
    struct bpred_btb_ent_t *stack; /* return-address stack */
  } retstack;

  /* stats */
  counter_t addr_hits;		/* num correct addr-predictions */
  counter_t dir_hits;		/* num correct dir-predictions (incl addr) */
  counter_t used_ras;		/* num RAS predictions used */
  counter_t used_bimod;		/* num bimodal predictions used (BPredComb) */
  counter_t used_2lev;		/* num 2-level predictions used (BPredComb) */
  counter_t jr_hits;		/* num correct addr-predictions for JR's */
  counter_t jr_seen;		/* num JR's seen */
  counter_t jr_non_ras_hits;	/* num correct addr-preds for non-RAS JR's */
  counter_t jr_non_ras_seen;	/* num non-RAS JR's seen */
  counter_t misses;		/* num incorrect predictions */
  
  counter_t lookups;		/* num lookups */

  counter_t replays;	        /* num of replays */

  counter_t retstack_pops;	/* number of times a value was popped */
  counter_t retstack_pushes;	/* number of times a value was pushed */
  counter_t ras_hits;		/* num correct return-address predictions */
  
  /* reverse stats */
  counter_t reverse_addr_hits;		/* num correct addr-predictions */
  counter_t reverse_dir_hits;		/* num correct dir-predictions (incl addr) */
  counter_t reverse_used_ras;		/* num RAS predictions used */
  counter_t reverse_used_bimod;		/* num bimodal predictions used (BPredComb) */
  counter_t reverse_used_2lev;		/* num 2-level predictions used (BPredComb) */
  counter_t reverse_jr_hits;		/* num correct addr-predictions for JR's */
  counter_t reverse_jr_seen;		/* num JR's seen */
  counter_t reverse_jr_non_ras_hits;	/* num correct addr-preds for non-RAS JR's */
  counter_t reverse_jr_non_ras_seen;	/* num non-RAS JR's seen */
  counter_t reverse_misses;		/* num incorrect reverse predictions */

  counter_t reverse_lookups;		/* num lookups */

  counter_t reverse_replays;	        /* num of replays */

  counter_t reverse_retstack_pops;	/* number of times a value was popped */
  counter_t reverse_retstack_pushes;	/* number of times a value was pushed */
  counter_t reverse_ras_hits;		/* num correct return-address predictions */
};

/* branch predictor update information */
struct bpred_update_t {
  char *fwd_pdir1;		/* direction-1 predictor counter */
  char *fwd_pdir2;		/* direction-2 predictor counter */
  char *fwd_pmeta;		/* meta predictor counter */
  struct {		/* predicted directions */
    unsigned int ras    : 1;	/* RAS used */
    unsigned int bimod  : 1;    /* bimodal predictor */
    unsigned int twolev : 1;    /* 2-level predictor */
    unsigned int meta   : 1;    /* meta predictor (0..bimod / 1..2lev) */
	int sum;				/* SC tables sum */
  } fwd_dir;
  int fwd_tage_pred;		/*	For TSCL/LLBP, track what the tage pred was*/
  
  /* Reversible Structs */
  char *rev_pdir1;		/* direction-1 predictor counter */
  char *rev_pdir2;		/* direction-2 predictor counter */
  char *rev_pmeta;		/* meta predictor counter */
  struct {		/* predicted directions */
    unsigned int ras    : 1;	/* RAS used */
    unsigned int bimod  : 1;    /* bimodal predictor */
    unsigned int twolev : 1;    /* 2-level predictor */
    unsigned int meta   : 1;    /* meta predictor (0..bimod / 1..2lev) */
	int sum;			/* SC tables sum */
  } rev_dir;
  int rev_tage_pred;		/*	For TSCL/LLBP, track what the tage pred was*/
};

/* create a branch predictor */
struct bpred_t *			/* branch predictory instance */
bpred_create(enum bpred_class class,	/* type of predictor to create */
	     unsigned int bimod_size,	/* bimod table size */
	     unsigned int l1size,	/* level-1 table size */
	     unsigned int l2size,	/* level-2 table size */
	     unsigned int meta_size,	/* meta predictor table size */
	     unsigned int shift_width,	/* history register width */
	     unsigned int xor,		/* history xor address flag */
       	     unsigned int head_table_width, /*TS head table width*/
	     unsigned int cht_size,	 /*CHBP correctness history table size*/
	     unsigned int btb_sets,	/* number of sets in BTB */ 
	     unsigned int btb_assoc,	/* BTB associativity */
	     unsigned int retstack_size, /* num entries in ret-addr stack */
		 unsigned int frmt); /* FRMTs enabled flag */

/* Create FWD-to-REV Mapping Tables */	 
struct bpred_frmt_t *		/* FRMT Inst */
bpred_frmt_create (enum bpred_class class,	/* type of predictor to create */
	unsigned int budget,		// budget (total storage size in entries) 
	unsigned int m_depth);		// m_depth (depth of contexts to track)

/* create a branch direction predictor */
struct bpred_dir_t *		/* branch direction predictor instance */
bpred_dir_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int l1size,		/* level-1 table size */
  unsigned int l2size,		/* level-2 table size (if relevant) */
  unsigned int shift_width,	/* history register width */
  unsigned int xor);	   	/* history xor address flag */


/* create a TS branch direction predictor */
struct bpred_ts_t *		/* temporal stream branch predictor instance */
bpred_ts_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int ts_enabled,              /* TSBP Enabled Flag */
  unsigned int head_table_width,	 	/* head table width */
  unsigned int head_table_size);			/* header table size */

/* create a correction history branch direction predictor */
struct bpred_chbp_t *		/* temporal stream branch predictor instance */
bpred_chbp_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int chbp_enabled,              /* CHBP Enabled Flag */
  unsigned int cht_size);			/* header table size */
  
/* create an outcome history table predictor */
struct bpred_oht_t *		/* outcome history table instance */
bpred_oht_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int oht_size);			/* OHT size */
  
/* Create an instance of the Future History Buffer (FHB) */
void bpred_rcr_create (struct bpred_t *pred, int size);
  
/* create a future history buffer (FHB) */
void bpred_fhb_create (
	struct bpred_fhb_t *pred_fhb,	// bpred pointer 
	int size);				// FHB size in bits

/* Create an outcome buffer (OB) */
void bpred_ob_create (struct bpred_t *pred, int width);
	
/* print branch predictor configuration */
void
bpred_config(struct bpred_t *pred,	/* branch predictor instance */
	     FILE *stream);		/* output stream */

/* print predictor stats */
void
bpred_stats(struct bpred_t *pred,	/* branch predictor instance */
	    FILE *stream);		/* output stream */

/* register branch predictor stats */
void
bpred_reg_stats(struct bpred_t *pred,	/* branch predictor instance */
		struct stat_sdb_t *sdb);/* stats database */

/* reset stats after priming, if appropriate */
void bpred_after_priming(struct bpred_t *bpred);

/* probe a predictor for a next fetch address, the predictor is probed
   with branch address BADDR, the branch target is BTARGET (used for
   static predictors), and OP is the instruction opcode (used to simulate
   predecode bits; a pointer to the predictor state entry (or null for jumps)
   is returned in *DIR_UPDATE_PTR (used for updating predictor state),
   and the non-speculative top-of-stack is returned in stack_recover_idx 
   (used for recovering ret-addr stack after mis-predict).  */
md_addr_t				/* predicted branch target addr */
bpred_lookup(struct bpred_t *pred,	/* branch predictor instance */
	     md_addr_t baddr,		/* branch address */
	     md_addr_t btarget,		/* branch target if taken */
	     enum md_opcode op,		/* opcode of instruction */
	     int is_call,		/* non-zero if inst is fn call */
	     int is_return,		/* non-zero if inst is fn return */
	     struct bpred_update_t *dir_update_ptr, /* pred state pointer */
	     int *stack_recover_idx,	/* Non-speculative top-of-stack; * used on mispredict recovery */
		 int flow_mode,  /* Flow mode (0=FWD; 1=REV) */
		 int frmt);		/* FRMTs Enabled flag*/

/* Speculative execution can corrupt the ret-addr stack.  So for each
 * lookup we return the top-of-stack (TOS) at that point; a mispredicted
 * branch, as part of its recovery, restores the TOS using this value --
 * hopefully this uncorrupts the stack. */
void
bpred_recover(struct bpred_t *pred,	/* branch predictor instance */
	      md_addr_t baddr,		/* branch address */
	      int stack_recover_idx);	/* Non-speculative top-of-stack;
					 * used on mispredict recovery */

/* update the branch predictor, only useful for stateful predictors; updates
   entry for instruction type OP at address BADDR.  BTB only gets updated
   for branches which are taken.  Inst was determined to jump to
   address BTARGET and was taken if TAKEN is non-zero.  Predictor 
   statistics are updated with result of prediction, indicated by CORRECT and 
   PRED_TAKEN, predictor state to be updated is indicated by *DIR_UPDATE_PTR 
   (may be NULL for jumps, which shouldn't modify state bits).  Note if
   bpred_update is done speculatively, branch-prediction may get polluted. */
void
bpred_update(struct bpred_t *pred,	/* branch predictor instance */
	     md_addr_t baddr,		/* branch address */
	     md_addr_t btarget,		/* resolved branch target */
	     int taken,			/* non-zero if branch was taken */
	     int pred_taken,		/* non-zero if branch was pred taken */
	     int correct,		/* was earlier prediction correct? */
	     enum md_opcode op,		/* opcode of instruction */
	     struct bpred_update_t *dir_update_ptr, /* pred state pointer */
		 int flow_mode,  /* Flow mode (0=FWD; 1=REV) */
		 int frmt);			/* FRMTs Enabled flag */


#ifdef foo0
/* OBSOLETE */
/* dump branch predictor state (for debug) */
void
bpred_dump(struct bpred_t *pred,	/* branch predictor instance */
	   FILE *stream);		/* output stream */
#endif

#endif /* BPRED_H */
