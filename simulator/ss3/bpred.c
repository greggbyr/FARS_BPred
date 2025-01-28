/* bpred.c - branch predictor routines */

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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "host.h"
#include "misc.h"
#include "machine.h"
#include "bpred.h"

/* turn this on to enable the SimpleScalar 2.0 RAS bug */
/* #define RAS_BUG_COMPATIBLE */

/* create a branch predictor */
struct bpred_t *			/* branch predictory instance */
bpred_create(enum bpred_class class,	/* type of predictor to create */
	     unsigned int bimod_size,	/* bimod table size */
	     unsigned int l1size,	/* 2lev l1 table size */
	     unsigned int l2size,	/* 2lev l2 table size */
	     unsigned int meta_size,	/* meta table size */
	     unsigned int shift_width,	/* history register width */
	     unsigned int xor,  	/* history xor address flag */
     	 unsigned int head_table_width, /*TS head table width*/
	     unsigned int cht_size, /*CHBP correctness history table width*/
	     unsigned int btb_sets,	/* number of sets in BTB */ 
	     unsigned int btb_assoc,	/* BTB associativity */
	     unsigned int retstack_size, /* num entries in ret-addr stack */
	     unsigned int frmt)	/* FRMTs Enabled Flag */
{
  struct bpred_t *pred;

  if (!(pred = calloc(1, sizeof(struct bpred_t))))
    fatal("out of virtual memory");

  pred->class = class;

  switch (class) {
  case BPredComb:
    /* bimodal component */
    pred->fwd_dirpred.bimod = 
      bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used

	// If FRMTs enabled, only one is needed.
	if (frmt==1)
      pred->fwd_dirpred.frmt =
        bpred_frmt_create(class, l2size, shift_width);

    //Need reverse bimod now too
	if (frmt==0)
      pred->rev_dirpred.bimod =
        bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);
    	
    /* 2-level component */
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	//Need reverse 2-level now too
	if (frmt==0)
		pred->rev_dirpred.twolev = 
			bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	
	if (!(pred->fhb = (struct bpred_fhb_t*)malloc(pred->num_fhbs * sizeof(struct bpred_fhb_t))))
			fatal("cannot allocate FHBs");
	
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(&pred->fhb[0], (shift_width + 1));
	
    /* metapredictor component */
    pred->fwd_dirpred.meta = 
      bpred_dir_create(BPred2bit, meta_size, 0, 0, 0);
	  
	//Need reverse meta too
	if (frmt==0)
		pred->rev_dirpred.meta = 
			bpred_dir_create(BPred2bit, meta_size, 0, 0, 0);
	
    break;

  case BPred2Level:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(class, l1size, l2size, shift_width, xor);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used
	
	// If FRMTs enabled, only one is needed.
	if (frmt==1)
      pred->fwd_dirpred.frmt =
        bpred_frmt_create(class, l2size, shift_width);
	
	//Need reverse 2-level now too
	if (frmt==0)
		pred->rev_dirpred.twolev = 
			bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	//else //With FRMTs, rev tables are the same as fwd
	//	pred->rev_dirpred.twolev = pred->fwd_dirpred.twolev;
	  
	if (!(pred->fhb = (struct bpred_fhb_t*)malloc(pred->num_fhbs * sizeof(struct bpred_fhb_t))))
			fatal("cannot allocate FHBs");
	  
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(&pred->fhb[0], (shift_width + 1));
	
    break;
	
  case BPredOB:
	pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used
	
	// If FRMTs enabled, only one is needed.
	if (frmt==1)
      pred->fwd_dirpred.frmt =
        bpred_frmt_create(class, l2size, shift_width);
	
	//Need reverse 2-level now too
	if (frmt==0)
		pred->rev_dirpred.twolev = 
			bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	
	if (!(pred->fhb = (struct bpred_fhb_t*)malloc(pred->num_fhbs * sizeof(struct bpred_fhb_t))))
			fatal("cannot allocate FHBs");
	
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(&pred->fhb[0], (shift_width + 1));
	
	// Outcome Buffer (OB)
	bpred_ob_create(pred, 16384); /* Hard setting OB width to 16k for now*/
	
    break;
	
  case BPredOHT:
	pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used
	 
	// If FRMTs enabled, only one is needed.
	if (frmt==1)
      pred->fwd_dirpred.frmt =
        bpred_frmt_create(class, l2size, shift_width);
	 
	//Need reverse 2-level now too
	if (frmt==0)
		pred->rev_dirpred.twolev = 
			bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
		
	if (!(pred->fhb = (struct bpred_fhb_t*)malloc(pred->num_fhbs * sizeof(struct bpred_fhb_t))))
			fatal("cannot allocate FHBs");
		
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(&pred->fhb[0], (shift_width + 1));
	
	//Need FWD and REV OHTs
	pred->fwd_dirpred.oht = 
      bpred_oht_create(class, l2size);
	
	// Leaving OHTs separate for now since they would have the same Valid bits
	if (frmt==0)
		pred->rev_dirpred.oht = 
			bpred_oht_create(class, l2size);
			
    break;
	
  case BPredMBP:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used
	
	// If FRMTs enabled, only one is needed.
	if (frmt==1)
      pred->fwd_dirpred.frmt =
        bpred_frmt_create(class, l2size, shift_width);
	
	//Need reverse 2-level now too
	if (frmt==0)
		pred->rev_dirpred.twolev = 
			bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
		
	if (!(pred->fhb = (struct bpred_fhb_t*)malloc(pred->num_fhbs * sizeof(struct bpred_fhb_t))))
			fatal("cannot allocate FHBs");
		
	// Reversible 2-Level also needs one FHB
	bpred_fhb_create(&pred->fhb[0], (shift_width + 1));
	
	//Need FWD and REV OHTs
	pred->fwd_dirpred.oht = 
      bpred_oht_create(class, l2size);
	  
	// Leaving OHTs separate for now since they would have the same Valid bits
	if (frmt==0)
		pred->rev_dirpred.oht = 
			bpred_oht_create(class, l2size);
		
	// Outcome Buffer (OB)
	bpred_ob_create(pred, 16384); /* Hard setting OB width to 16k for now*/
	
    break;	

	case BPredLLBP:
		//Last-Level Branch Predictor consists of a TAGE-SC-L predictor along with a prefetching cache
		//which contains longer histories of hard to predict patterns
		
		// LLBT
		// Allocating "shift_width" wide Pattern Buffer 
		if (!(pred->fwd_pb_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
			fatal("cannot allocate LLBP tables");
		
		for (int i = 0; i < shift_width; i++) {
			pred->fwd_pb_dirpred[i].twolev = 
				bpred_dir_create(class, 1, (l2size/(2<<(shift_width-2))), (int)(pow((double)xor,(double)(shift_width-2))*(double)l1size+0.5), 0);
		}
		
		// Allocating "shift_width" wide LLBT
		if (!(pred->fwd_llbt_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
			fatal("cannot allocate LLBP tables");
		
		for (int i = 0; i < shift_width; i++) {
			pred->fwd_llbt_dirpred[i].twolev = 
				bpred_dir_create(class, 1, (4*(l2size/(2<<(shift_width-2)))), (int)(pow((double)xor,(double)(shift_width-2))*(double)l1size+0.5), 0);
		}
		
		// Need reverse table as well
		if (!frmt) {
			// Allocating "shift_width" wide Pattern Buffer 
			if (!(pred->rev_pb_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
				fatal("cannot allocate LLBP tables");
			
			for (int i = 0; i < shift_width; i++) {
				pred->rev_pb_dirpred[i].twolev = 
					bpred_dir_create(class, 1, (l2size/(2<<(shift_width-2))), (int)(pow((double)xor,(double)(shift_width-2))*(double)l1size+0.5), 0);
			}
			
			// Allocating "shift_width" wide LLBT
			if (!(pred->rev_llbt_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
				fatal("cannot allocate LLBP tables");
			
			for (int i = 0; i < shift_width; i++) {
				pred->rev_llbt_dirpred[i].twolev = 
					bpred_dir_create(class, 1, (4*(l2size/(2<<(shift_width-2)))), (int)(pow((double)xor,(double)(shift_width-2))*(double)l1size+0.5), 0);
			}
		}
		
		// Rolling Context Register (RCR); only one needed for reversible LLBP
		bpred_rcr_create(pred, shift_width);
		
	case BPredTSCL:	
		// TAGE-SC-L
		pred->tage_depth = shift_width;	//Gate TAGE tables depth for later prediction/updates
		pred->num_fhbs = 1;	// Also set the number of FHBs equal to TAGE tables
		
		// Base pred is bimod
		pred->fwd_dirpred.bimod = 
			bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);
		
		// TAGE and SC table allocation based on shift width
		if (!(pred->fwd_tage_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
			fatal("cannot allocate TAGE tables");
		
		if (!(pred->fwd_sc_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
			fatal("cannot allocate SC tables");
		
		for (int i = 1; i < shift_width; i++) {
			pred->fwd_tage_dirpred[i].twolev = 
				bpred_dir_create(class, 1, (l2size/(2<<(i-1))), (int)(pow((double)xor,(double)(i-1))*(double)l1size+0.5), 1);
			
			pred->fwd_sc_dirpred[i].twolev = 
				bpred_dir_create(class, 1, (l2size/(2<<(i-1))), (int)(pow((double)xor,(double)(i-1))*(double)l1size+0.5), 1);
		}
		
		// Loop Tables
		pred->fwd_loop_dirpred.twolev = 
			bpred_dir_create(class, 1, l1size, shift_width, 1);
			
		//Need OHT
		pred->fwd_dirpred.oht =
			bpred_oht_create(class, l2size);
		
		// Need the REV unless FRMTs enabled
		// If FRMTs enabled, only one is needed.
		if (frmt==1) {
			pred->fwd_dirpred.frmt =
				bpred_frmt_create(class, bimod_size, shift_width);
		} else {
			// TAGE-SC-L
			
			// Base pred is bimod
			pred->rev_dirpred.bimod =
				bpred_dir_create(BPred2bit, bimod_size, 0, 0, 0);
			
			// TAGE and SC table allocation based on shift width
			if (!(pred->rev_tage_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
				fatal("cannot allocate TAGE tables");
			
			if (!(pred->rev_sc_dirpred = (struct bpred_dirpred_t*)malloc(shift_width * sizeof(struct bpred_dirpred_t))))
				fatal("cannot allocate SC tables");
			
			for (int i = 1; i < shift_width; i++) {
				pred->rev_tage_dirpred[i].twolev = 
					bpred_dir_create(class, 1, (l2size/(2<<(i-1))), (int)(pow((double)xor,(double)(i-1))*(double)l1size+0.5), 1);
				
				pred->rev_sc_dirpred[i].twolev = 
					bpred_dir_create(class, 1, (l2size/(2<<(i-1))), (int)(pow((double)xor,(double)(i-1))*(double)l1size+0.5), 1);
			}
			
			// Loop Tables
			pred->rev_loop_dirpred.twolev = 
				bpred_dir_create(class, 1, l1size, shift_width, 1);
			
			// Need OHT
			pred->rev_dirpred.oht = 
				bpred_oht_create(class, l2size);
		}
		
		if (!(pred->fhb = (struct bpred_fhb_t*)malloc(pred->num_fhbs * sizeof(struct bpred_fhb_t))))
			fatal("cannot allocate FHBs");
			
		// Reversible TSCL/LLBP also needs one FHB for every TAGE table
		bpred_fhb_create(&pred->fhb[0], ((int)(pow((double)xor,(double)(shift_width-2))*(double)l1size+0.5) + 1));
		
		// Reversible TSCL/LLBP also needs one Outcome Buffer (OB)
		bpred_ob_create(pred, 16384); /* Hard setting OB width to 16k for now*/
    break;

  case BPredTSBP:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used
	
	pred->fwd_dirpred.tsbp =
	  bpred_ts_create(class, 1, head_table_width, ((unsigned int)l2size << 3)); /* md_addr_t is the size of the PC*/
    break;
	
  case BPredCHBP:
    pred->fwd_dirpred.twolev = 
      bpred_dir_create(BPred2Level, l1size, l2size, shift_width, xor);
	  
	pred->num_fhbs = 1; //To work with TAGE, should always be 1 more than actually used
	
	pred->fwd_dirpred.chbp =
	  bpred_chbp_create(class, 1, cht_size);
    break;

  case BPred2bit:
    pred->fwd_dirpred.bimod = 
      bpred_dir_create(class, bimod_size, 0, 0, 0);

	// If FRMTs enabled, only one is needed.
	if (frmt==1)
      pred->fwd_dirpred.frmt =
        bpred_frmt_create(class, bimod_size, shift_width);

    //Need reverse bimod now too
    if (frmt==0)
		pred->rev_dirpred.bimod =
			bpred_dir_create(class, bimod_size, 0, 0, 0);
	
  case BPredTaken:
  case BPredNotTaken:
    /* no other state */
    break;

  default:
    panic("bogus predictor class");
  }

  /* allocate ret-addr stack */
  switch (class) {
  case BPredComb:
  case BPred2Level:
  case BPredTSBP:
  case BPredCHBP:
  case BPred2bit:
  case BPredOB:
  case BPredOHT:
  case BPredMBP:
  case BPredTSCL:
  case BPredLLBP:
    {
      int i;

      /* allocate BTB */
      if (!btb_sets || (btb_sets & (btb_sets-1)) != 0)
	fatal("number of BTB sets must be non-zero and a power of two");
      if (!btb_assoc || (btb_assoc & (btb_assoc-1)) != 0)
	fatal("BTB associativity must be non-zero and a power of two");

      if (!(pred->btb.btb_data = calloc(btb_sets * btb_assoc,
					sizeof(struct bpred_btb_ent_t))))
	fatal("cannot allocate BTB");

      pred->btb.sets = btb_sets;
      pred->btb.assoc = btb_assoc;

      if (pred->btb.assoc > 1)
	for (i=0; i < (pred->btb.assoc*pred->btb.sets); i++)
	  {
	    if (i % pred->btb.assoc != pred->btb.assoc - 1)
	      pred->btb.btb_data[i].next = &pred->btb.btb_data[i+1];
	    else
	      pred->btb.btb_data[i].next = NULL;
	    
	    if (i % pred->btb.assoc != pred->btb.assoc - 1)
	      pred->btb.btb_data[i+1].prev = &pred->btb.btb_data[i];
	  }

      /* allocate retstack */
      if ((retstack_size & (retstack_size-1)) != 0)
	fatal("Return-address-stack size must be zero or a power of two");
      
      pred->retstack.size = retstack_size;
      if (retstack_size)
	if (!(pred->retstack.stack = calloc(retstack_size, 
					    sizeof(struct bpred_btb_ent_t))))
	  fatal("cannot allocate return-address-stack");
      pred->retstack.tos = retstack_size - 1;
      
      break;
    }

  case BPredTaken:
  case BPredNotTaken:
    /* no other state */
    break;

  default:
    panic("bogus predictor class");
  }

  return pred;
}

struct bpred_frmt_t *		/* FRMT Inst */
bpred_frmt_create (enum bpred_class class,	/* type of predictor to create */
	unsigned int budget,		// budget (total storage size in entries) 
	unsigned int m_depth		// m_depth (depth of contexts to track)
) {
	struct bpred_frmt_t *frmts;
	
	if (!(frmts = calloc(1, sizeof(struct bpred_frmt_t))))
		fatal("out of virtual memory");
	
	frmts->class = class;
	frmts->budget = budget;
	
	switch (class) {
		case BPredLLBP:
			//Allocate FRMT context entries which are needed only for LLBP
			if (!(frmts->context = calloc(budget, sizeof(md_addr_t))))
				fatal("cannot allocate FRMT context array"); 
		
		case BPredComb:
		case BPred2Level:
		case BPredTSBP:
		case BPredCHBP:
		case BPredOB:
		case BPredOHT:
		case BPredMBP:
		case BPredTSCL:
			//Allocate FRMT tag entries which are needed for all tag-type preds starting at BPredComb and below 
			if (!(frmts->tag = calloc(budget, sizeof(md_addr_t))))
				fatal("cannot allocate FRMT tag array");
		
		case BPred2bit:
			//Allocate FRMT addr entries which are needed for 2-bit and all following predictors
			if (!(frmts->addr = calloc(budget, sizeof(md_addr_t))))
				fatal("cannot allocate FRMT tag array");
		
		case BPredTaken:
		case BPredNotTaken:
			break;
		
		/* no other state */
		default:
			panic("bogus predictor class");
	}
			
	return frmts;
}

/* Create an instance of the Outcome Buffer (OB) */
void bpred_ob_create (struct bpred_t *pred, int width) {
	// Outcome Buffer (OB)
	pred->ob.width = width; /* hardcoding to 16 kb for now, will make configurable later */
	
	pred->ob.end = width - 1;	/* End of buffer window at width - 1 */
	
	pred->ob.beg = 0;	/* Beginning of buffer window at 0 */
	
	if (!(pred->ob.fv = calloc(width, sizeof(int))))
		fatal("cannot allocate OB FV bits");
	
	if (!(pred->ob.rv = calloc(width, sizeof(int))))
		fatal("cannot allocate OB RV bits");
	
	if (!(pred->ob.oc = calloc(width, sizeof(int))))
		fatal("cannot allocate OB OC bits");
}

/* Create an instance of the Future History Buffer (FHB) */
void bpred_fhb_create (struct bpred_fhb_t *pred_fhb, int size) {
	pred_fhb->size = size;
	
	pred_fhb->top = size - 1;
	
	pred_fhb->bot = 0;
	
	pred_fhb->fv_out = 0;
	pred_fhb->rv_out = 0;
	pred_fhb->o_out = 0;
	pred_fhb->correct_out = 0;
	pred_fhb->und_out = 0;
	pred_fhb->thu_out = 0;
	pred_fhb->addr_out = 0;
	pred_fhb->key_out = 0;
	
	
	if (!(pred_fhb->fv = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB FV bits");
	
	if (!(pred_fhb->rv = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB RV bits");
	
	if (!(pred_fhb->o = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB O bits");
	
	if (!(pred_fhb->correct = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB Correctness bits");
	
	if (!(pred_fhb->und = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB under threshold bits");
	
	if (!(pred_fhb->thu = calloc(size, sizeof(int))))
		fatal("cannot allocate FHB threshold update bits");
	
	if (!(pred_fhb->addr = calloc(size, sizeof(md_addr_t))))
		fatal("cannot allocate FHB ADDR bits");
	
	//For FRMT setups, need to also keep track of key from current direction
	if (!(pred_fhb->key = calloc(size, sizeof(md_addr_t))))
		fatal("cannot allocate FHB key bits");
}

/* Create an instance of the Future History Buffer (FHB) */
void bpred_rcr_create (struct bpred_t *pred, int size) {
	pred->rcr.size = size+4;	//Need extra 4 entries; 2 for prefetch and 2 for post fetch 
	
	pred->rcr.top = pred->rcr.size - 1;
	
	pred->rcr.win_top = pred->rcr.top - 2;
	
	pred->rcr.fpre_top = pred->rcr.win_top - 2;
	
	pred->rcr.rpre_bot = 4;
	
	pred->rcr.win_bot = 2;
	
	pred->rcr.bot = 0;
	
	pred->rcr.fccid = 0;
	pred->rcr.ccid = 0;
	pred->rcr.rccid = 0;
	
	if (!(pred->rcr.fv = calloc(size, sizeof(int))))
		fatal("cannot allocate RCR FV bits");
	
	if (!(pred->rcr.rv = calloc(size, sizeof(int))))
		fatal("cannot allocate RCR RV bits");
	
	if (!(pred->rcr.addr = calloc(size, sizeof(md_addr_t))))
		fatal("cannot allocate RCR ADDR bits");
}

/* create a branch direction predictor */
struct bpred_dir_t *		/* branch direction predictor instance */
bpred_dir_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int l1size,	 	/* level-1 table size */
  unsigned int l2size,	 	/* level-2 table size (if relevant) */
  unsigned int shift_width,	/* history register width */
  unsigned int xor)	    	/* history xor address flag */
{
  struct bpred_dir_t *pred_dir;
  unsigned int cnt;
  int flipflop;

  if (!(pred_dir = calloc(1, sizeof(struct bpred_dir_t))))
    fatal("out of virtual memory");

  pred_dir->class = class;

  cnt = -1;
  switch (class) {
  case BPredLLBP:
  case BPredTSCL:
  case BPred2Level:
    {
	if (!l1size || (l1size & (l1size-1)) != 0)
		fatal("level-1 size, `%d', must be non-zero and a power of two", l1size);
	
	pred_dir->config.two.l1size = l1size;
      
	if (!l2size || (l2size & (l2size-1)) != 0) {
		fatal("`%d' level-2 size, `%d', must be non-zero and a power of two; l1size=%d", class, l2size, l1size);
	}
	
	pred_dir->config.two.l2size = l2size;
      
	if (!shift_width || shift_width > 30)
		fatal("shift register width, `%d', must be non-zero and positive", shift_width);
    
	pred_dir->config.two.shift_width = shift_width;
      
	pred_dir->config.two.xor = xor;
 
	if (!(pred_dir->config.two.shiftregs = calloc(l1size, sizeof(int))))
		fatal("cannot allocate shift register table");
      
	if (!(pred_dir->config.two.l2table = calloc(l2size, sizeof(unsigned char))))
		fatal("cannot allocate second level table");

	/* initialize counters to weakly this-or-that */
	flipflop = 1;
	for (cnt = 0; cnt < l2size; cnt++) {
		pred_dir->config.two.l2table[cnt] = flipflop;
		flipflop = 3 - flipflop;
	}

	pred_dir->config.two.threshold = 6;			/*	SC/LOOP threshold */
	
	if (!(pred_dir->config.two.context = calloc(l2size, sizeof(int))))
		fatal("cannot allocate context entries");
	
	if (!(pred_dir->config.two.hist_length = calloc(l2size, sizeof(int))))
		fatal("cannot allocate history length entries");
	
	if (!(pred_dir->config.two.use = calloc(l2size, sizeof(unsigned char))))
		fatal("cannot allocate useful counter entries");
	
	if (!(pred_dir->config.two.iter_c = calloc(l2size, sizeof(unsigned char))))
		fatal("cannot allocate current iteration counter entries");
	
	if (!(pred_dir->config.two.iter_p = calloc(l2size, sizeof(unsigned char))))
		fatal("cannot allocate past iteration counter entries");
	
	if (!(pred_dir->config.two.tag = calloc(l2size, sizeof(md_addr_t))))
		fatal("cannot allocate tag entries");		/* tags for TAGE tables */
	
      break;
    }

  case BPred2bit:
    if (!l1size || (l1size & (l1size-1)) != 0)
      fatal("2bit table size, `%d', must be non-zero and a power of two", 
	    l1size);
    pred_dir->config.bimod.size = l1size;
    if (!(pred_dir->config.bimod.table =
	  calloc(l1size, sizeof(unsigned char))))
      fatal("cannot allocate 2bit storage");
    /* initialize counters to weakly this-or-that */
    flipflop = 1;
    for (cnt = 0; cnt < l1size; cnt++)
      {
	pred_dir->config.bimod.table[cnt] = flipflop;
	flipflop = 3 - flipflop;
      }

    break;

  case BPredTaken:
  case BPredNotTaken:
    /* no other state */
    break;

  default:
    panic("bogus branch direction predictor class");
  }

  return pred_dir;
}

/* create a branch direction predictor */
struct bpred_ts_t *		/* temporal stream branch predictor instance */
bpred_ts_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int ts_enabled,              /* TSBP Enabled Flag */
  unsigned int head_table_width,	 	/* head table width */
  unsigned int head_table_size)			/* head table size */
{
  struct bpred_ts_t *pred_ts;
  unsigned int key;
  
  if (!(pred_ts = calloc(1, sizeof(struct bpred_ts_t))))
    fatal("out of virtual memory");

  pred_ts->class = class;

  if (!head_table_size || (head_table_size & (head_table_size - 1)) != 0)
	fatal("head table size, `%d', must be non-zero and a power of two", head_table_size);
  
  pred_ts->ts.head_table_size = head_table_size;
  
  if (!head_table_width || head_table_width > 30)
	fatal("head table width, `%d', must be non-zero and positive", head_table_width);
  
  pred_ts->ts.head_table_width = head_table_width;
 
  pred_ts->ts.head_table = calloc(head_table_size, sizeof(unsigned int));
  //if (head_table_width % 8) { 
  //	pred_ts->ts.head_table = calloc(head_table_size, ((head_table_width / 8) + 1));
  //} else {
//	pred_ts->ts.head_table = calloc(head_table_size, (head_table_width / 8));
  //}
  
  if (!pred_ts->ts.head_table)
	fatal("cannot allocate head table");

  /* initializing head table entries to NULL*/
  for (key = 0; key < head_table_size; key++)
	  pred_ts->ts.head_table[key] = NULL;

  pred_ts->ts.correctness_width = 2 << (head_table_width - 1);
  
  pred_ts->ts.correctness_buffer = calloc(pred_ts->ts.correctness_width, sizeof(bool_t));

  if (!pred_ts->ts.correctness_buffer)
	fatal("cannot allocate correctness buffer");

   /* initializing CB bits to 1*/
  for (key = 0; key < pred_ts->ts.correctness_width; key++)
  	pred_ts->ts.correctness_buffer[key] = 1;

  /* initialize current head of the correctness buffer */
  pred_ts->ts.head = 0;

  /* initialize current tail of the correctness buffer */
  pred_ts->ts.tail = 0;
  
  /* initialize replay to false */
  pred_ts->ts.replay = FALSE;

  /* initialize enabled flag */
  if (ts_enabled == 0)
          pred_ts->ts.enabled = FALSE;
  else
          pred_ts->ts.enabled = TRUE;

  return pred_ts;
}

/* create a branch direction predictor */
struct bpred_chbp_t *		/* temporal stream branch predictor instance */
bpred_chbp_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int chbp_enabled,              /* CHBP Enabled Flag */
  unsigned int cht_size)			/* Correctness History Table size */
{
  struct bpred_chbp_t *pred_chbp;
  //unsigned int key;
  
  if (!(pred_chbp = calloc(1, sizeof(struct bpred_chbp_t))))
    fatal("out of virtual memory");

  pred_chbp->class = class;

  if (!cht_size || (cht_size & (cht_size - 1)) != 0)
	fatal("correctness history table size, `%d', must be non-zero and a power of two", cht_size);
  
  pred_chbp->chbp.cht_size = cht_size;
  
  /* Need to allocate for source PC, replay bit, correctness bit, valid bit, and destination PC*/
  pred_chbp->chbp.cht_spc = calloc(cht_size, sizeof(md_addr_t)); 
  
  if (!pred_chbp->chbp.cht_spc)
	fatal("cannot allocate correctness history table source pc bits");

  pred_chbp->chbp.cht_replay = calloc(cht_size, sizeof(bool_t)); 
  
  if (!pred_chbp->chbp.cht_replay)
	fatal("cannot allocate correctness history table replay bits");

  pred_chbp->chbp.cht_correct = calloc(cht_size, sizeof(bool_t)); 
  
  if (!pred_chbp->chbp.cht_correct)
	fatal("cannot allocate correctness history table correct bits");

  pred_chbp->chbp.cht_valid = calloc(cht_size, sizeof(bool_t)); 
  
  if (!pred_chbp->chbp.cht_valid)
	fatal("cannot allocate correctness history table valid bits");

  pred_chbp->chbp.cht_dpc = calloc(cht_size, sizeof(md_addr_t)); 
  
  if (!pred_chbp->chbp.cht_dpc)
	fatal("cannot allocate correctness history table destination pc bits");

  /* initialize enabled flag */
  if (chbp_enabled == 0)
          pred_chbp->chbp.enabled = FALSE;
  else
          pred_chbp->chbp.enabled = TRUE;

  return pred_chbp;
}

/* create an outcome history table */
struct bpred_oht_t *		/* outcome history table instance */
bpred_oht_create (
  enum bpred_class class,	/* type of predictor to create */
  unsigned int oht_size)			/* Outcome History Table size */
{
  struct bpred_oht_t *pred_oht;
  //unsigned int key;
  
  if (!(pred_oht = calloc(1, sizeof(struct bpred_oht_t))))
    fatal("out of virtual memory");

  pred_oht->class = class;

  if (!oht_size || (oht_size & (oht_size - 1)) != 0)
	fatal("outcome history table size, `%d', must be non-zero and a power of two", oht_size);
  
  pred_oht->oht.size = oht_size;
  
  /* Need to allocate for outcome (OC) and valid bits */
  pred_oht->oht.oc = calloc(oht_size, sizeof(int)); 
  
  if (!pred_oht->oht.oc)
	fatal("cannot allocate outcome history table outcome (pc) bits");

  pred_oht->oht.fwd_valid = calloc(oht_size, sizeof(int)); 
  
  if (!pred_oht->oht.fwd_valid)
	fatal("cannot allocate outcome history table fwd valid bits");

  pred_oht->oht.rev_valid = calloc(oht_size, sizeof(int)); 
  
  if (!pred_oht->oht.rev_valid)
	fatal("cannot allocate outcome history table rev valid bits");

  return pred_oht;
}

/* print branch direction predictor configuration */
void
bpred_dir_config(
  struct bpred_dir_t *pred_dir,	/* branch direction predictor instance */
  char name[],			/* predictor name */
  FILE *stream)			/* output stream */
{
  switch (pred_dir->class) {
  case BPred2Level:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredOB:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredOHT:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredMBP:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPredTSBP:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;
	
  case BPredCHBP:
    fprintf(stream,
      "pred_dir: %s: 2-lvl: %d l1-sz, %d bits/ent, %s xor, %d l2-sz, direct-mapped\n",
      name, pred_dir->config.two.l1size, pred_dir->config.two.shift_width,
      pred_dir->config.two.xor ? "" : "no", pred_dir->config.two.l2size);
    break;

  case BPred2bit:
    fprintf(stream, "pred_dir: %s: 2-bit: %d entries, direct-mapped\n",
      name, pred_dir->config.bimod.size);
    break;

  case BPredTaken:
    fprintf(stream, "pred_dir: %s: predict taken\n", name);
    break;

  case BPredNotTaken:
    fprintf(stream, "pred_dir: %s: predict not taken\n", name);
    break;

  default:
    panic("bogus branch direction predictor class");
  }
}

/* print temporal stream predictor configuration */
void
bpred_ts_config(
  struct bpred_ts_t *pred_ts,	/* branch direction predictor instance */
  char name[],			/* predictor name */
  FILE *stream)			/* output stream */
{
  fprintf(stream,
    "pred_ts: %s: tsbp: %d ht-sz, %d ht-wd, %d cr-wd\n",
    name, pred_ts->ts.head_table_size, pred_ts->ts.head_table_width,
	pred_ts->ts.correctness_width);
}

/* print mississippi branch predictor configuration */
void
bpred_chbp_config(
  struct bpred_chbp_t *pred_chbp,	/* branch direction predictor instance */
  char name[],			/* predictor name */
  FILE *stream)			/* output stream */
{
  fprintf(stream,
    "pred_chbp: %s: chbp: %d cht-sz\n",
    name, pred_chbp->chbp.cht_size);
}

/* print branch predictor configuration */
void
bpred_config(struct bpred_t *pred,	/* branch predictor instance */
	     FILE *stream)		/* output stream */
{
  switch (pred->class) {
  case BPredComb:
    bpred_dir_config (pred->fwd_dirpred.bimod, "bimod", stream);
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    bpred_dir_config (pred->fwd_dirpred.meta, "meta", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPred2Level:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;
	
  case BPredOB:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity",
            pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPredOHT:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity",
            pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;
  
  case BPredMBP:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
    fprintf(stream, "btb: %d sets x %d associativity",
            pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPredTSBP:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
	bpred_ts_config (pred->fwd_dirpred.tsbp, "tsbp", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;
	
  case BPredCHBP:
    bpred_dir_config (pred->fwd_dirpred.twolev, "2lev", stream);
	bpred_chbp_config (pred->fwd_dirpred.chbp, "chbp", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPred2bit:
    bpred_dir_config (pred->fwd_dirpred.bimod, "bimod", stream);
    fprintf(stream, "btb: %d sets x %d associativity", 
	    pred->btb.sets, pred->btb.assoc);
    fprintf(stream, "ret_stack: %d entries", pred->retstack.size);
    break;

  case BPredTaken:
    bpred_dir_config (pred->fwd_dirpred.bimod, "taken", stream);
    break;
  case BPredNotTaken:
    bpred_dir_config (pred->fwd_dirpred.bimod, "nottaken", stream);
    break;

  default:
    panic("bogus branch predictor class");
  }
}

/* print predictor stats */
void
bpred_stats(struct bpred_t *pred,	/* branch predictor instance */
	    FILE *stream)		/* output stream */
{
	/* stats */
  fprintf(stream, "pred: addr-prediction rate = %f\n",
	  (double)pred->addr_hits/(double)(pred->addr_hits+pred->misses));
  fprintf(stream, "pred: dir-prediction rate = %f\n",
	  (double)pred->dir_hits/(double)(pred->dir_hits+pred->misses));
  
	/* reverse stats */
  fprintf(stream, "pred: reverse addr-prediction rate = %f\n",
	  (double)pred->reverse_addr_hits/(double)(pred->reverse_addr_hits+pred->reverse_misses));
  fprintf(stream, "pred: reverse dir-prediction rate = %f\n",
	  (double)pred->reverse_dir_hits/(double)(pred->reverse_dir_hits+pred->reverse_misses));
}

/* register branch predictor stats */
void
bpred_reg_stats(struct bpred_t *pred,	/* branch predictor instance */
		struct stat_sdb_t *sdb)	/* stats database */
{
  char buf[512], buf1[512], *name;

  /* get a name for this predictor */
  switch (pred->class)
    {
    case BPredComb:
      name = "bpred_comb";
      break;
    case BPred2Level:
      name = "bpred_2lev";
      break;
     case BPredOB:
      name = "bpred_ob";
      break;
     case BPredOHT:
      name = "bpred_oht";
      break;
     case BPredMBP:
      name = "bpred_mbp";
      break;
	case BPredTSCL:
      name = "bpred_tscl";
      break;
	case BPredLLBP:
      name = "bpred_llbp";
      break;
	case BPredTSBP:
      name = "bpred_tsbp";
      break;
	case BPredCHBP:
      name = "bpred_chbp";
      break;
    case BPred2bit:
      name = "bpred_bimod";
      break;
    case BPredTaken:
      name = "bpred_taken";
      break;
    case BPredNotTaken:
      name = "bpred_nottaken";
      break;
    default:
      panic("bogus branch predictor class");
    }

  sprintf(buf, "%s.lookups", name);
  stat_reg_counter(sdb, buf, "total number of bpred lookups",
		   &pred->lookups, 0, NULL);
  sprintf(buf, "%s.updates", name);
  sprintf(buf1, "%s.dir_hits + %s.misses", name, name);
  stat_reg_formula(sdb, buf, "total number of updates", buf1, "%12.0f");
  sprintf(buf, "%s.addr_hits", name);
  stat_reg_counter(sdb, buf, "total number of address-predicted hits", 
		   &pred->addr_hits, 0, NULL);
  sprintf(buf, "%s.dir_hits", name);
  stat_reg_counter(sdb, buf, 
		   "total number of direction-predicted hits "
		   "(includes addr-hits)", 
		   &pred->dir_hits, 0, NULL);
  if (pred->class == BPredComb)
    {
      sprintf(buf, "%s.used_bimod", name);
      stat_reg_counter(sdb, buf, 
		       "total number of bimodal predictions used", 
		       &pred->used_bimod, 0, NULL);
      sprintf(buf, "%s.used_2lev", name);
      stat_reg_counter(sdb, buf, 
		       "total number of 2-level predictions used", 
		       &pred->used_2lev, 0, NULL);
    }
  if (pred->class == BPredTSBP || pred->class == BPredCHBP) {
       sprintf(buf, "%s.replays", name);
       stat_reg_counter(sdb, buf, "total number of replays", &pred->replays, 0, NULL);
  }
  sprintf(buf, "%s.misses", name);
  stat_reg_counter(sdb, buf, "total number of misses", &pred->misses, 0, NULL);
  sprintf(buf, "%s.jr_hits", name);
  stat_reg_counter(sdb, buf,
		   "total number of address-predicted hits for JR's",
		   &pred->jr_hits, 0, NULL);
  sprintf(buf, "%s.jr_seen", name);
  stat_reg_counter(sdb, buf,
		   "total number of JR's seen",
		   &pred->jr_seen, 0, NULL);
  sprintf(buf, "%s.jr_non_ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of address-predicted hits for non-RAS JR's",
		   &pred->jr_non_ras_hits, 0, NULL);
  sprintf(buf, "%s.jr_non_ras_seen.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of non-RAS JR's seen",
		   &pred->jr_non_ras_seen, 0, NULL);
  sprintf(buf, "%s.bpred_addr_rate", name);
  sprintf(buf1, "%s.addr_hits / %s.updates", name, name);
  stat_reg_formula(sdb, buf,
		   "branch address-prediction rate (i.e., addr-hits/updates)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.bpred_dir_rate", name);
  sprintf(buf1, "%s.dir_hits / %s.updates", name, name);
  stat_reg_formula(sdb, buf,
		  "branch direction-prediction rate (i.e., all-hits/updates)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.bpred_jr_rate", name);
  sprintf(buf1, "%s.jr_hits / %s.jr_seen", name, name);
  stat_reg_formula(sdb, buf,
		  "JR address-prediction rate (i.e., JR addr-hits/JRs seen)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.bpred_jr_non_ras_rate.PP", name);
  sprintf(buf1, "%s.jr_non_ras_hits.PP / %s.jr_non_ras_seen.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "non-RAS JR addr-pred rate (ie, non-RAS JR hits/JRs seen)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.retstack_pushes", name);
  stat_reg_counter(sdb, buf,
		   "total number of address pushed onto ret-addr stack",
		   &pred->retstack_pushes, 0, NULL);
  sprintf(buf, "%s.retstack_pops", name);
  stat_reg_counter(sdb, buf,
		   "total number of address popped off of ret-addr stack",
		   &pred->retstack_pops, 0, NULL);
  sprintf(buf, "%s.used_ras.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of RAS predictions used",
		   &pred->used_ras, 0, NULL);
  sprintf(buf, "%s.ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of RAS hits",
		   &pred->ras_hits, 0, NULL);
  sprintf(buf, "%s.ras_rate.PP", name);
  sprintf(buf1, "%s.ras_hits.PP / %s.used_ras.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "RAS prediction rate (i.e., RAS hits/used RAS)",
		   buf1, "%9.4f");
		   
		   
		  /*  Reverse Data Sets */
		   
  sprintf(buf, "%s.reverse_lookups", name);
  stat_reg_counter(sdb, buf, "total number of reverse bpred lookups",
		   &pred->reverse_lookups, 0, NULL);
  sprintf(buf, "%s.reverse_updates", name);
  sprintf(buf1, "%s.reverse_dir_hits + %s.reverse_misses", name, name);
  stat_reg_formula(sdb, buf, "total number of reverse updates", buf1, "%12.0f");
  sprintf(buf, "%s.reverse_addr_hits", name);
  stat_reg_counter(sdb, buf, "total number of reverse address-predicted hits", 
		   &pred->reverse_addr_hits, 0, NULL);
  sprintf(buf, "%s.reverse_dir_hits", name);
  stat_reg_counter(sdb, buf, 
		   "total number of reverse direction-predicted hits "
		   "(includes addr-hits)", 
		   &pred->reverse_dir_hits, 0, NULL);
  if (pred->class == BPredComb)
    {
      sprintf(buf, "%s.reverse_used_bimod", name);
      stat_reg_counter(sdb, buf, 
		       "total number of reverse bimodal predictions used", 
		       &pred->reverse_used_bimod, 0, NULL);
      sprintf(buf, "%s.reverse_used_2lev", name);
      stat_reg_counter(sdb, buf, 
		       "total number of reverse 2-level predictions used", 
		       &pred->reverse_used_2lev, 0, NULL);
    }
  if (pred->class == BPredTSBP || pred->class == BPredCHBP) {
       sprintf(buf, "%s.reverse_replays", name);
       stat_reg_counter(sdb, buf, "total number of reverse replays", &pred->reverse_replays, 0, NULL);
  }
  sprintf(buf, "%s.reverse_misses", name);
  stat_reg_counter(sdb, buf, "total number of reverse misses", &pred->reverse_misses, 0, NULL);
  sprintf(buf, "%s.reverse_jr_hits", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address-predicted hits for JR's",
		   &pred->reverse_jr_hits, 0, NULL);
  sprintf(buf, "%s.reverse_jr_seen", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse JR's seen",
		   &pred->reverse_jr_seen, 0, NULL);
  sprintf(buf, "%s.reverse_jr_non_ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address-predicted hits for non-RAS JR's",
		   &pred->reverse_jr_non_ras_hits, 0, NULL);
  sprintf(buf, "%s.reverse_jr_non_ras_seen.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse non-RAS JR's seen",
		   &pred->reverse_jr_non_ras_seen, 0, NULL);
  sprintf(buf, "%s.reverse_bpred_addr_rate", name);
  sprintf(buf1, "%s.reverse_addr_hits / %s.reverse_updates", name, name);
  stat_reg_formula(sdb, buf,
		   "branch address-prediction rate (i.e., addr-hits/updates)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.reverse_bpred_dir_rate", name);
  sprintf(buf1, "%s.reverse_dir_hits / %s.reverse_updates", name, name);
  stat_reg_formula(sdb, buf,
		  "branch direction-prediction rate (i.e., all-hits/updates)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.reverse_bpred_jr_rate", name);
  sprintf(buf1, "%s.reverse_jr_hits / %s.reverse_jr_seen", name, name);
  stat_reg_formula(sdb, buf,
		  "JR address-prediction rate (i.e., JR addr-hits/JRs seen)",
		  buf1, "%9.4f");
  sprintf(buf, "%s.reverse_bpred_jr_non_ras_rate.PP", name);
  sprintf(buf1, "%s.reverse_jr_non_ras_hits.PP / %s.reverse_jr_non_ras_seen.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "non-RAS JR addr-pred rate (ie, non-RAS JR hits/JRs seen)",
		   buf1, "%9.4f");
  sprintf(buf, "%s.reverse_retstack_pushes", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address pushed onto ret-addr stack",
		   &pred->reverse_retstack_pushes, 0, NULL);
  sprintf(buf, "%s.reverse_retstack_pops", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse address popped off of ret-addr stack",
		   &pred->reverse_retstack_pops, 0, NULL);
  sprintf(buf, "%s.reverse_used_ras.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse RAS predictions used",
		   &pred->reverse_used_ras, 0, NULL);
  sprintf(buf, "%s.reverse_ras_hits.PP", name);
  stat_reg_counter(sdb, buf,
		   "total number of reverse RAS hits",
		   &pred->reverse_ras_hits, 0, NULL);
  sprintf(buf, "%s.reverse_ras_rate.PP", name);
  sprintf(buf1, "%s.reverse_ras_hits.PP / %s.reverse_used_ras.PP", name, name);
  stat_reg_formula(sdb, buf,
		   "RAS prediction rate (i.e., RAS hits/used RAS)",
		   buf1, "%9.4f");
}

void
bpred_after_priming(struct bpred_t *bpred)
{
  if (bpred == NULL)
    return;

	/* Stats */
  bpred->lookups = 0;
  bpred->addr_hits = 0;
  bpred->dir_hits = 0;
  bpred->used_ras = 0;
  bpred->used_bimod = 0;
  bpred->used_2lev = 0;
  bpred->jr_hits = 0;
  bpred->jr_seen = 0;
  bpred->misses = 0;
  bpred->replays = 0;
  bpred->retstack_pops = 0;
  bpred->retstack_pushes = 0;
  bpred->ras_hits = 0;
  
	/* Reverse Stats */ //Set to 1 for now to avoid runtime issues
  bpred->reverse_lookups = 0;
  bpred->reverse_addr_hits = 0;
  bpred->reverse_dir_hits = 0;
  bpred->reverse_used_ras = 0;
  bpred->reverse_used_bimod = 0;
  bpred->reverse_used_2lev = 0;
  bpred->reverse_jr_hits = 0;
  bpred->reverse_jr_seen = 0;
  bpred->reverse_misses = 0;
  bpred->reverse_replays = 0;
  bpred->reverse_retstack_pops = 0;
  bpred->reverse_retstack_pushes = 0;
  bpred->reverse_ras_hits = 0;
}

#define BIMOD_HASH(PRED, ADDR)						\
  ((((ADDR) >> 19) ^ ((ADDR) >> MD_BR_SHIFT)) & ((PRED)->config.bimod.size-1))
    /* was: ((baddr >> 16) ^ baddr) & (pred->fwd_dirpred.bimod.size-1) */

int key_from_llbp_features (struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
		 md_addr_t baddr,		/* branch address */
		 int ccid)			/* Current Context ID */				
{
	int l1index, key;

    /* traverse 2-level tables */
    l1index = (baddr >> MD_BR_SHIFT) & (pred_dir->config.two.l1size - 1);
    key = pred_dir->config.two.shiftregs[l1index];
        
	if (pred_dir->config.two.xor) {
	    /* this L2 index computation is more "compatible" to McFarling's
	       verison of it, i.e., if the PC xor address component is only
	       part of the index, take the lower order address bits for the
	       other part of the index, rather than the higher order ones */
#if 1  
	    key = (((key ^ (baddr >> MD_BR_SHIFT))
			& ((1 << pred_dir->config.two.hist_length[ccid]) - 1))
		    | ((baddr >> MD_BR_SHIFT)
			<< pred_dir->config.two.hist_length[ccid]));
#else   
		key = key ^ (baddr >> MD_BR_SHIFT);
#endif
	} else {
		key = key | ((baddr >> MD_BR_SHIFT) << pred_dir->config.two.hist_length[ccid]);
	}
	
	return key;
}

/* Used to calculate future 2nd level table indexes/keys for TSCL/LLBP*/
int future_key_from_tage (struct bpred_t *pred,	/* branch pred inst*/
		struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
		 md_addr_t baddr,			/* branch address */
		 int tage_table,			/* TAGE table num */
		 int flow_mode)				/* flow_mode flag for determining key bits*/
{
	int l1index, key;

    /* traverse 2-level tables */
    l1index = (baddr >> MD_BR_SHIFT) & (pred_dir->config.two.l1size - 1);
    key = pred_dir->config.two.shiftregs[l1index];
	
	// Need to hijack the key depending on the flow_mode and tage table num
	if (flow_mode) {	// Need to shift key right bc REV mode the latest history is on the left
		key = key >> (pred_dir->config.two.shift_width - pred->fwd_tage_dirpred[tage_table].twolev->config.two.shift_width);
	}
	
	key = key ^ ((1 << pred->fwd_tage_dirpred[tage_table].twolev->config.two.shift_width) - 1);
        
	if (pred_dir->config.two.xor) {
	    /* this L2 index computation is more "compatible" to McFarling's
	       verison of it, i.e., if the PC xor address component is only
	       part of the index, take the lower order address bits for the
	       other part of the index, rather than the higher order ones */
#if 1  
	    key = (((key ^ (baddr >> MD_BR_SHIFT))
			& ((1 << pred_dir->config.two.shift_width) - 1))
			| ((baddr >> MD_BR_SHIFT)
			<< pred_dir->config.two.shift_width));
#else   
		key = key ^ (baddr >> MD_BR_SHIFT);
#endif
	} else {
		key = key | ((baddr >> MD_BR_SHIFT) << pred_dir->config.two.shift_width);
	}
	
	return key;
}

/* Used to calculate 2nd level table indexes/keys for 2lev, tsbp, and chbp*/
int key_from_features (struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
		 md_addr_t baddr)		/* branch address */
{
	int l1index, key;

    /* traverse 2-level tables */
    l1index = (baddr >> MD_BR_SHIFT) & (pred_dir->config.two.l1size - 1);
    key = pred_dir->config.two.shiftregs[l1index];
        
	if (pred_dir->config.two.xor) {
	    /* this L2 index computation is more "compatible" to McFarling's
	       verison of it, i.e., if the PC xor address component is only
	       part of the index, take the lower order address bits for the
	       other part of the index, rather than the higher order ones */
#if 1  
	    key = (((key ^ (baddr >> MD_BR_SHIFT))
			& ((1 << pred_dir->config.two.shift_width) - 1))
		    | ((baddr >> MD_BR_SHIFT)
			<< pred_dir->config.two.shift_width));
#else   
		key = key ^ (baddr >> MD_BR_SHIFT);
#endif
	} else {
		key = key | ((baddr >> MD_BR_SHIFT) << pred_dir->config.two.shift_width);
	}
	
	return key;
}

/* FRMT Functions for mapping and getting keys, contexts, and addrs */
int unmasked_key_from_frmt(
	struct bpred_t *pred, 		/* branch predictor instance */
	int key						/* branch key/tag to exchange */
) {
	key = key & (pred->fwd_dirpred.frmt->budget - 1); //Mask based on FRMT budget
	return pred->fwd_dirpred.frmt->tag[key];
}

void map_unmasked_key_to_frmt(
	struct bpred_t *pred, 		/* branch predictor instance */
	int fwd_key,						/* fwd branch key/tag to map */
	int rev_key						/* rev branch key/tag to map */
) {
	rev_key = rev_key & (pred->fwd_dirpred.frmt->budget - 1); //Mask based on FRMT budget
	pred->fwd_dirpred.frmt->tag[rev_key] = fwd_key;
}

md_addr_t addr_from_frmt(
	struct bpred_t *pred, 		/* branch predictor instance */
	int baddr						/* branch addr to exchange */
) {
	baddr = baddr & (pred->fwd_dirpred.frmt->budget - 1); //Mask based on FRMT budget
	return pred->fwd_dirpred.frmt->addr[baddr];
}

void map_addr_to_frmt(
	struct bpred_t *pred, 		/* branch predictor instance */
	int fbaddr,						/* fwd branch addr to map */
	int rbaddr						/* rev branch addr to map */
) {
	rbaddr = rbaddr & (pred->fwd_dirpred.frmt->budget - 1); //Mask based on FRMT budget
	pred->fwd_dirpred.frmt->addr[rbaddr] = fbaddr;
}

int context_from_frmt(
	struct bpred_t *pred, 		/* branch predictor instance */
	int cont						/* branch context to exchange */
) {
	cont = cont & (pred->fwd_dirpred.frmt->budget - 1); //Mask based on FRMT budget
	return pred->fwd_dirpred.frmt->context[cont];
}

void map_context_to_frmt(
	struct bpred_t *pred, 		/* branch predictor instance */
	int fcont,						/* fwd branch context to map */
	int rcont						/* rev branch context to map */
) {
	rcont = rcont & (pred->fwd_dirpred.frmt->budget - 1); //Mask based on FRMT budget
	pred->fwd_dirpred.frmt->context[rcont] = fcont;
}

/* predicts a branch direction */
char *						/* pointer to counter */
bpred_dir_lookup(struct bpred_t *pred,	/* branch predictor instance */
	struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
	md_addr_t baddr,		/* branch address */
	int flow_mode,				/* Flow mode (0=FWD; 1=REV) */
	int frmt)  				/* FRMTs enabled flag */
{
  unsigned char *p = NULL;

  /* Except for jumps, get a pointer to direction-prediction bits */
  switch (pred_dir->class) {
    case BPredLLBP:		/*Add LLBP case, should be same as bimod to get base prediction*/
		{
		//if (flow_mode && frmt)	// FRMT was enabled so REV ADDR needs to be exchanged for FWD ADDR
		//	baddr = addr_from_frmt(pred, baddr);
		
		int l2index = key_from_llbp_features (pred_dir, baddr, pred->rcr.ccid); //Get un-masked l2index  
		
		//if (flow_mode && frmt)	// FRMT was enabled so REV key needs to be exchanged for FWD key
		//	l2index = unmasked_key_from_frmt(pred, l2index);
		
        l2index = l2index & (pred_dir->config.two.l2size - 1);
		
        /* get a pointer to prediction state information */
        p = &pred_dir->config.two.l2table[l2index];
		}
		break;
	case BPred2Level:
    case BPredOB:
    case BPredOHT:
    case BPredMBP:  
    case BPredTSBP:         /*Add TSBP case, should be same as 2 level to get base prediction*/
	case BPredCHBP:         /*Add CHBP case, should be same as 2 level to get base prediction*/
    case BPredTSCL:		/*Add TAGE-SC-L case, should be same as bimod to get base prediction*/
		{
		//if (flow_mode && frmt)	// FRMT was enabled so REV ADDR needs to be exchanged for FWD ADDR
		//	baddr = addr_from_frmt(pred, baddr);
		
		int l2index = key_from_features (pred_dir, baddr); //Get un-masked l2index  
		
		//if (flow_mode && frmt)	// FRMT was enabled so REV key needs to be exchanged for FWD key
		//	l2index = unmasked_key_from_frmt(pred, l2index);
		
        l2index = l2index & (pred_dir->config.two.l2size - 1);
		
        /* get a pointer to prediction state information */
        p = &pred_dir->config.two.l2table[l2index];
      }
      break;
    case BPred2bit:
	{
		//if (flow_mode && frmt)	// FRMT was enabled so REV ADDR needs to be exchanged for FWD ADDR
		//	baddr = addr_from_frmt(pred, baddr);
	
		p = &pred_dir->config.bimod.table[BIMOD_HASH(pred_dir, baddr)];
	}
		break;
    case BPredTaken:
    case BPredNotTaken:
      break;
    default:
      panic("bogus branch direction predictor class");
    }

  return (char *)p;
}

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
	     int *stack_recover_idx,	/* Non-speculative top-of-stack; used on mispredict recovery */
		 int flow_mode,				/* Flow mode (0=FWD; 1=REV) */
		 int frmt)  				/* FRMTs enabled flag */
{
	struct bpred_btb_ent_t *pbtb = NULL;
	int index, i;
	bool_t invert = FALSE;
	  
	int fwd_valid_outcome = NULL;
	int rev_valid_outcome = NULL;

	if (!dir_update_ptr)
		panic("no bpred update record");

	/* if this is not a branch, return not-taken */
	if (!(MD_OP_FLAGS(op) & F_CTRL))
		return 0;
   
	// FWD & REV addresses; 
	// These are supposed to be opposites but 
	// without ISA updates and knowledge of addrs and targets, 
	// these need to be the same.
	md_addr_t fbaddr = baddr;		/* FWD Branch ADDR */		
	md_addr_t fbtarget = btarget;	/* FWD Branch Target */
	md_addr_t rbaddr = baddr;		/* REV Branch ADDR; for now, setting to baddr */
	md_addr_t rbtarget = btarget;	/* REV Branch Target; for now, setting to btarget */
	
	if (frmt)	// FRMT was enabled so REV ADDR needs to be exchanged for FWD ADDR
		rbaddr = addr_from_frmt(pred, rbaddr);

	if (!flow_mode) {
		pred->lookups++;
	} else {
		pred->reverse_lookups++;
	}

	// FWD DIR Pointers
	dir_update_ptr->fwd_dir.ras = FALSE;
	dir_update_ptr->fwd_pdir1 = NULL;
	dir_update_ptr->fwd_pdir2 = NULL;
	dir_update_ptr->fwd_pmeta = NULL;
	dir_update_ptr->fwd_tage_pred = NULL;
	
	// REV DIR Pointers
	dir_update_ptr->rev_dir.ras = FALSE;
	dir_update_ptr->rev_pdir1 = NULL;
	dir_update_ptr->rev_pdir2 = NULL;
	dir_update_ptr->rev_pmeta = NULL;
	dir_update_ptr->rev_tage_pred = NULL;
	
	/* Except for jumps, get a pointer to direction-prediction bits */
	switch (pred->class) {
		case BPredComb:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				char *bimod, *twolev, *meta;
			  
				bimod = bpred_dir_lookup(pred, pred->fwd_dirpred.bimod, fbaddr, 0, frmt);
				twolev = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev, fbaddr, 0, frmt);
				meta = bpred_dir_lookup(pred, pred->fwd_dirpred.meta, fbaddr, 0, frmt);
				  
				dir_update_ptr->fwd_pmeta = meta;
				dir_update_ptr->fwd_dir.meta  = (*meta >= 2);
				dir_update_ptr->fwd_dir.bimod = (*bimod >= 2);
				dir_update_ptr->fwd_dir.twolev  = (*twolev >= 2);
				  
				if (*meta >= 2) {
					dir_update_ptr->fwd_pdir1 = twolev;
					dir_update_ptr->fwd_pdir2 = bimod;
				} else {
					dir_update_ptr->fwd_pdir1 = bimod;
					dir_update_ptr->fwd_pdir2 = twolev;
				}
			}
			break;
			
		case BPredTSBP:   
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev, fbaddr, 0, frmt);  //get 2level base outcome prediction
				
				int key = key_from_features (pred->fwd_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
				
				key = key & (pred->fwd_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
				
				/* incr head but prevent from going out of bounds*/
				if (pred->fwd_dirpred.tsbp->ts.head >= pred->fwd_dirpred.tsbp->ts.correctness_width) {
					pred->fwd_dirpred.tsbp->ts.head = 0;
				} else {
					pred->fwd_dirpred.tsbp->ts.head++;
				}

				/*if in replay mode and corretness buffer head indicates base predictor mistake*/
				if(pred->fwd_dirpred.tsbp->ts.replay 
					&& pred->fwd_dirpred.tsbp->ts.enabled 
					&& (pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.head] == 0)) {
					invert = TRUE; 
				}
			}
			break;
			
		case BPredCHBP:   
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				dir_update_ptr->fwd_pdir1 = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev, fbaddr, 0, frmt);  //get 2level base outcome prediction
				
				int key = key_from_features (pred->fwd_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
				
				key = key & (pred->fwd_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size

				/*if enabled, replay bit set, correctness bits is 0, and src_pc matches fbaddr predictor is inverted*/
				if (pred->fwd_dirpred.chbp->chbp.enabled 
					&& pred->fwd_dirpred.chbp->chbp.cht_replay[key] 				// Only perform correction if replay is on
					&& !pred->fwd_dirpred.chbp->chbp.cht_correct[key] 			// Check for past correctness history
					&& (pred->fwd_dirpred.chbp->chbp.cht_spc[key] == fbaddr)) { 	// Check stored source pc is same as fbaddr
					invert = TRUE; 
				}
			}
			break;
		
		case BPred2Level:
		case BPredOB:
		case BPredOHT:
		case BPredMBP:
		case BPredTSCL:
		case BPredLLBP:
		case BPred2bit:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				//Check if any valid past outcome results for FWD (and REV if FRMT set)
				if (((pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||(pred->class==BPredOB)) &&
				(!flow_mode && pred->ob.fv[pred->ob.end])) {
					//Check if OB has valid past outcome results
					fwd_valid_outcome = pred->ob.oc[pred->ob.end];
				}
				
				if ((fwd_valid_outcome==NULL) && (pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||
				(pred->class==BPredOB)||(pred->class==BPredOHT)||(pred->class==BPred2Level)) {
					if (!flow_mode && pred->fhb[0].fv[pred->fhb[0].top]) {
						//Check if FHB has valid past outcome results
						fwd_valid_outcome = pred->fhb[0].o[pred->fhb[0].top];
					}
				}
					
				if ((fwd_valid_outcome==NULL) && (pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||
				(pred->class==BPredOHT)) {
					//Check if OHT has valid past outcome results
					int key;
					
					if ((pred->class==BPredMBP)||(pred->class==BPredOHT))
						key = key_from_features (pred->fwd_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
					else 
						key = key_from_features (pred->fwd_tage_dirpred[pred->tage_depth-1].twolev, fbaddr); // Get unmasked key from GHR and PC
					
					key = key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on predictor table size
					
					if (!flow_mode && pred->fwd_dirpred.oht->oht.fwd_valid[key]) {
						fwd_valid_outcome = pred->fwd_dirpred.oht->oht.oc[key];
					}
				}			
					
				if ((pred->class==BPredMBP)||(pred->class==BPredOB)||(pred->class==BPredOHT)||(pred->class==BPred2Level)) {
					dir_update_ptr->fwd_pdir1 = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev, fbaddr, 0, frmt);
				} else {
					char *bimod, *tage, *loop, *llbt;
					int tage_match = 0;
					int loop_match = 0;
					int llbt_match = 0;
					
					bimod = bpred_dir_lookup(pred, pred->fwd_dirpred.bimod, fbaddr, 0, frmt);
					
					if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
						dir_update_ptr->fwd_dir.sum = 0;
						
						if (pred->class==BPredLLBP) {
							for (int i = pred->tage_depth-1; i >= 1; i--) {
								// Get LLBT tag and compare to current tag
								int llbt_key = key_from_features (pred->fwd_pb_dirpred[i].twolev, fbaddr); // Get unmasked key from GHR and PC
								
								// LLBP is indexed by RCR CCID!
								// If matched unmasked TAGE key, set tage pointer and note it and hist length
								if (llbt_key==pred->fwd_pb_dirpred[i].twolev->config.two.tag[pred->rcr.ccid] && 
								pred->rcr.ccid==pred->fwd_pb_dirpred[i].twolev->config.two.context[pred->rcr.ccid]) {
									llbt = bpred_dir_lookup(pred, pred->fwd_pb_dirpred[i].twolev, fbaddr, 0, frmt);
									llbt_match = pred->fwd_pb_dirpred[i].twolev->config.two.hist_length[pred->rcr.ccid];
									break;	// Break at first (highest) match
								}
							}
						}
						
						for (int i = pred->tage_depth-1; i >= 1; i--) {
							// Get tage and sc tags and compare to current tag
							int tage_key = key_from_features (pred->fwd_tage_dirpred[i].twolev, fbaddr); // Get unmasked key from GHR and PC
							int masked_tage_key = tage_key & (pred->fwd_tage_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
							
							// If matched unmasked TAGE key, set tage pointer
							if (tage_key==pred->fwd_tage_dirpred[i].twolev->config.two.tag[masked_tage_key]) {
								tage = bpred_dir_lookup(pred, pred->fwd_tage_dirpred[i].twolev, fbaddr, 0, frmt);
								tage_match = pred->fwd_tage_dirpred[i].twolev->config.two.shift_width;
								break;	// Break at first (highest) match
							}
						}
						
						for (int i = 1; i < pred->tage_depth; i++) {
							int sc_key = key_from_features (pred->fwd_sc_dirpred[i].twolev, fbaddr); // Get unmasked key from GHR and PC
							int masked_sc_key = sc_key & (pred->fwd_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
							
							if (sc_key==pred->fwd_sc_dirpred[i].twolev->config.two.tag[masked_sc_key]) {
								dir_update_ptr->fwd_dir.sum = dir_update_ptr->fwd_dir.sum + (pred->fwd_sc_dirpred[i].twolev->config.two.l2table[masked_sc_key] >= 2); 
							}
						}
						
						// Get loop tag and compare to current tag
						int loop_key = key_from_features (pred->fwd_loop_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
						int masked_loop_key = loop_key & (pred->fwd_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
						
						if (loop_key==pred->fwd_loop_dirpred.twolev->config.two.tag[masked_loop_key]) {
							loop = bpred_dir_lookup(pred, pred->fwd_loop_dirpred.twolev, fbaddr, 0, frmt);
							
							if (*loop >= pred->fwd_loop_dirpred.twolev->config.two.threshold)
								loop_match = 1;
						}
					}
					
					if (loop_match) {
						dir_update_ptr->fwd_pdir1 = loop;
					} else {
						if (llbt_match>tage_match) {
							dir_update_ptr->fwd_pdir1 = llbt;
						} else if (tage_match) {
							dir_update_ptr->fwd_pdir1 = tage;
						} else {
							dir_update_ptr->fwd_pdir1 = bimod;
						}
						
						// Check for statistical correction
						if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
							dir_update_ptr->fwd_dir.sum = dir_update_ptr->fwd_dir.sum + (*(dir_update_ptr->fwd_pdir1) >= 2);
							
							// invert on overthreshold summation
							if (dir_update_ptr->fwd_dir.sum >= pred->fwd_sc_dirpred[pred->tage_depth-1].twolev->config.two.threshold)
								invert = TRUE;
						}
					}
					
					// Save TAGE prediction as a valid outcome could have been used.
					if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
						dir_update_ptr->fwd_tage_pred = (*(dir_update_ptr->fwd_pdir1) >= 2);
					}
				}
			}
			break;
		
		case BPredTaken:
			return fbtarget;
			
		case BPredNotTaken:
			if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
				return fbaddr + sizeof(md_inst_t);
			} else {
				return fbtarget;
			}
			
		default:
			panic("bogus predictor class");
	}
	
	// REV Predictors
	/* Except for jumps, get a pointer to direction-prediction bits */
	if (frmt) {
		switch (pred->class) {
			case BPredComb:
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					char *bimod, *twolev, *meta;
					
					bimod = bpred_dir_lookup(pred, pred->fwd_dirpred.bimod,  rbaddr, 1, frmt);
					twolev = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev,  rbaddr, 1, frmt);
					meta = bpred_dir_lookup(pred, pred->fwd_dirpred.meta,  rbaddr, 1, frmt);
					
					dir_update_ptr->rev_pmeta = meta;
					dir_update_ptr->rev_dir.meta  = (*meta >= 2);
					dir_update_ptr->rev_dir.bimod = (*bimod >= 2);
					dir_update_ptr->rev_dir.twolev  = (*twolev >= 2);
					
					if (*meta >= 2) {
						dir_update_ptr->rev_pdir1 = twolev;
						dir_update_ptr->rev_pdir2 = bimod;
					} else {
						dir_update_ptr->rev_pdir1 = bimod;
						dir_update_ptr->rev_pdir2 = twolev;
					}
				}
				break;
				
			case BPredTSBP:   
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					dir_update_ptr->rev_pdir1 = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev,  rbaddr, 1, frmt);  //get 2level base outcome prediction
					
					int key = key_from_features (pred->fwd_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
					
					key = key & (pred->fwd_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
					
					/* incr head but prevent from going out of bounds*/
					if (pred->fwd_dirpred.tsbp->ts.head >= pred->fwd_dirpred.tsbp->ts.correctness_width) {
						pred->fwd_dirpred.tsbp->ts.head = 0;
					} else {
						pred->fwd_dirpred.tsbp->ts.head++;
					}

					/*if in replay mode and corretness buffer head indicates base predictor mistake*/
					if(pred->fwd_dirpred.tsbp->ts.replay 
						&& pred->fwd_dirpred.tsbp->ts.enabled 
						&& (pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.head] == 0)) {
						invert = TRUE; 
					}
				}
				break;
				
			case BPredCHBP:   
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					dir_update_ptr->rev_pdir1 = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev,  rbaddr, 1, frmt);  //get 2level base outcome prediction
					
					int key = key_from_features (pred->fwd_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
					
					key = key & (pred->fwd_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size

					/*if enabled, replay bit set, correctness bits is 0, and src_pc matches rbaddr predictor is inverted*/
					if (pred->fwd_dirpred.chbp->chbp.enabled 
						&& pred->fwd_dirpred.chbp->chbp.cht_replay[key] 				// Only perform correction if replay is on
						&& !pred->fwd_dirpred.chbp->chbp.cht_correct[key] 			// Check for past correctness history
						&& (pred->fwd_dirpred.chbp->chbp.cht_spc[key] == rbaddr)) { 	// Check stored source pc is same as rbaddr
						invert = TRUE; 
					}
				}
				break;
				
			case BPred2Level:
			case BPredOB:
			case BPredOHT:
			case BPredMBP:
			case BPredTSCL:
			case BPredLLBP:
			case BPred2bit:
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					//Check if any valid past outcome results
					if (((pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||(pred->class==BPredOB)) &&
					(flow_mode && pred->ob.rv[pred->ob.beg])) {
						//Check if OB has valid past outcome results
						rev_valid_outcome = pred->ob.oc[pred->ob.beg];
					} 
					
					if ((rev_valid_outcome==NULL) && (pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||
					(pred->class==BPredOB)||(pred->class==BPredOHT)||(pred->class==BPred2Level)) {
						if (flow_mode && pred->fhb[0].rv[pred->fhb[0].bot]) {
							//Check if FHB has valid past outcome results
							rev_valid_outcome = pred->fhb[0].o[pred->fhb[0].bot];
						} 
					}
						
					if ((rev_valid_outcome==NULL) && (pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||
					(pred->class==BPredOHT)) {
						//Check if OHT has valid past outcome results
						int key;
						
						if ((pred->class==BPredMBP)||(pred->class==BPredOHT))
							key = key_from_features (pred->fwd_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
						else 
							key = key_from_features (pred->fwd_tage_dirpred[pred->tage_depth-1].twolev, rbaddr); // Get unmasked key from GHR and PC
					
						// FRMT was enabled so REV key needs to be exchanged for FWD key
						//key = unmasked_key_from_frmt(pred, key);
					
						key = key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on predictor table size
						
						if (flow_mode && pred->fwd_dirpred.oht->oht.rev_valid[key]) {
							rev_valid_outcome = pred->fwd_dirpred.oht->oht.oc[key];
						}
					} 
						
					if ((pred->class==BPredMBP)||(pred->class==BPredOB)||(pred->class==BPredOHT)||(pred->class==BPred2Level)) {
						dir_update_ptr->rev_pdir1 = bpred_dir_lookup(pred, pred->fwd_dirpred.twolev, rbaddr, 1, frmt);
					} else {
						char *bimod, *tage, *loop, *llbt;
						int tage_match = 0;
						int loop_match = 0;
						int llbt_match = 0;
						
						bimod = bpred_dir_lookup(pred, pred->fwd_dirpred.bimod, rbaddr, 1, frmt);
						
						if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
							dir_update_ptr->fwd_dir.sum = 0;
							
							if (pred->class==BPredLLBP) {
								for (int i = pred->tage_depth-1; i >= 1; i--) {
									// Get LLBT tag and compare to current tag
									int llbt_key = key_from_features (pred->fwd_pb_dirpred[i].twolev, rbaddr); // Get unmasked key from GHR and PC
									
									// LLBP is indexed by RCR CCID!
									// If matched unmasked TAGE key, set tage pointer and note it and hist length
									if (llbt_key==pred->fwd_pb_dirpred[i].twolev->config.two.tag[pred->rcr.ccid] && 
									pred->rcr.ccid==pred->fwd_pb_dirpred[i].twolev->config.two.context[pred->rcr.ccid]) {
										llbt = bpred_dir_lookup(pred, pred->fwd_pb_dirpred[i].twolev, rbaddr, 1, frmt);
										llbt_match = pred->fwd_pb_dirpred[i].twolev->config.two.hist_length[pred->rcr.ccid];
										break;	// Break at first (highest) match
									}
								}
							}
							
							for (int i = pred->tage_depth-1; i >= 1; i--) {
								// Get tage and sc tags and compare to current tag
								int tage_key = key_from_features (pred->fwd_tage_dirpred[i].twolev, rbaddr); // Get unmasked key from GHR and PC
								int masked_tage_key = tage_key & (pred->fwd_tage_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
								
								// If matched unmasked TAGE key, set tage pointer
								if (tage_key==pred->fwd_tage_dirpred[i].twolev->config.two.tag[masked_tage_key]) {
									tage = bpred_dir_lookup(pred, pred->fwd_tage_dirpred[i].twolev, rbaddr, 1, frmt);
									tage_match = pred->fwd_tage_dirpred[i].twolev->config.two.shift_width;
									break;	// Break at first (highest) match
								}
							}
							
							for (int i = 1; i < pred->tage_depth; i++) {
								int sc_key = key_from_features (pred->fwd_sc_dirpred[i].twolev, rbaddr); // Get unmasked key from GHR and PC
								int masked_sc_key = sc_key & (pred->fwd_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
								
								if (sc_key==pred->fwd_sc_dirpred[i].twolev->config.two.tag[masked_sc_key]) {
									dir_update_ptr->fwd_dir.sum = dir_update_ptr->fwd_dir.sum + (pred->fwd_sc_dirpred[i].twolev->config.two.l2table[masked_sc_key] >= 2); 
								}
							}
							
							// Get loop tag and compare to current tag
							int loop_key = key_from_features (pred->fwd_loop_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
							int masked_loop_key = loop_key & (pred->fwd_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
							
							if (loop_key==pred->fwd_loop_dirpred.twolev->config.two.tag[masked_loop_key]) {
								loop = bpred_dir_lookup(pred, pred->fwd_loop_dirpred.twolev, rbaddr, 1, frmt);
								
								if (*loop >= pred->fwd_loop_dirpred.twolev->config.two.threshold)
									loop_match = 1;
							}
						}
						
						if (loop_match) {
							dir_update_ptr->rev_pdir1 = loop;
						} else {
							if (llbt_match>tage_match) {
								dir_update_ptr->rev_pdir1 = llbt;
							} else if (tage_match) {
								dir_update_ptr->rev_pdir1 = tage;
							} else {
								dir_update_ptr->rev_pdir1 = bimod;
							}
							
							// Check for statistical correction
							if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
								dir_update_ptr->fwd_dir.sum = dir_update_ptr->fwd_dir.sum + (*(dir_update_ptr->rev_pdir1) >= 2);
								
								// invert on overthreshold summation
								if (dir_update_ptr->fwd_dir.sum >= pred->fwd_sc_dirpred[pred->tage_depth-1].twolev->config.two.threshold)
									invert = TRUE;
							}
						}
						
						// Save TAGE prediction as a valid outcome could have been used.
						if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
							dir_update_ptr->fwd_tage_pred = (*(dir_update_ptr->rev_pdir1) >= 2);
						}
					}
				}
				break;
				
			case BPredTaken:
				return rbtarget;
				
			case BPredNotTaken:
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					return rbaddr + sizeof(md_inst_t);
				} else {
					return rbtarget;
				}
				
			default:
				panic("bogus predictor class");
		}
	} else {
		switch (pred->class) {
			case BPredComb:
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					char *bimod, *twolev, *meta;
					
					bimod = bpred_dir_lookup(pred, pred->rev_dirpred.bimod,  rbaddr, 1, frmt);
					twolev = bpred_dir_lookup(pred, pred->rev_dirpred.twolev,  rbaddr, 1, frmt);
					meta = bpred_dir_lookup(pred, pred->rev_dirpred.meta,  rbaddr, 1, frmt);
					
					dir_update_ptr->rev_pmeta = meta;
					dir_update_ptr->rev_dir.meta  = (*meta >= 2);
					dir_update_ptr->rev_dir.bimod = (*bimod >= 2);
					dir_update_ptr->rev_dir.twolev  = (*twolev >= 2);
					
					if (*meta >= 2) {
						dir_update_ptr->rev_pdir1 = twolev;
						dir_update_ptr->rev_pdir2 = bimod;
					} else {
						dir_update_ptr->rev_pdir1 = bimod;
						dir_update_ptr->rev_pdir2 = twolev;
					}
				}
				break;
				
			case BPredTSBP:   
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					dir_update_ptr->rev_pdir1 = bpred_dir_lookup(pred, pred->rev_dirpred.twolev,  rbaddr, 1, frmt);  //get 2level base outcome prediction
					
					int key = key_from_features (pred->rev_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
					
					key = key & (pred->rev_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
					
					/* incr head but prevent from going out of bounds*/
					if (pred->rev_dirpred.tsbp->ts.head >= pred->rev_dirpred.tsbp->ts.correctness_width) {
						pred->rev_dirpred.tsbp->ts.head = 0;
					} else {
						pred->rev_dirpred.tsbp->ts.head++;
					}

					/*if in replay mode and corretness buffer head indicates base predictor mistake*/
					if(pred->rev_dirpred.tsbp->ts.replay 
						&& pred->rev_dirpred.tsbp->ts.enabled 
						&& (pred->rev_dirpred.tsbp->ts.correctness_buffer[pred->rev_dirpred.tsbp->ts.head] == 0)) {
						invert = TRUE; 
					}
				}
				break;
				
			case BPredCHBP:   
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					dir_update_ptr->rev_pdir1 = bpred_dir_lookup(pred, pred->rev_dirpred.twolev,  rbaddr, 1, frmt);  //get 2level base outcome prediction
					
					int key = key_from_features (pred->rev_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
					
					key = key & (pred->rev_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size

					/*if enabled, replay bit set, correctness bits is 0, and src_pc matches rbaddr predictor is inverted*/
					if (pred->rev_dirpred.chbp->chbp.enabled 
						&& pred->rev_dirpred.chbp->chbp.cht_replay[key] 				// Only perform correction if replay is on
						&& !pred->rev_dirpred.chbp->chbp.cht_correct[key] 			// Check for past correctness history
						&& (pred->rev_dirpred.chbp->chbp.cht_spc[key] == rbaddr)) { 	// Check stored source pc is same as rbaddr
						invert = TRUE; 
					}
				}
				break;
				
			case BPred2Level:
			case BPredOB:
			case BPredOHT:
			case BPredMBP:
			case BPredTSCL:
			case BPredLLBP:
			case BPred2bit:
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					//Check if any valid past outcome results
					if (((pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||(pred->class==BPredOB))&&
					(flow_mode && pred->ob.rv[pred->ob.beg])) {
						//Check if OB has valid past outcome results
						rev_valid_outcome = pred->ob.oc[pred->ob.beg];
					} 
					
					if ((rev_valid_outcome==NULL) && ((pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||
					(pred->class==BPredOB)||(pred->class==BPredOHT)||(pred->class==BPred2Level))) {
						if (flow_mode && pred->fhb[0].rv[pred->fhb[0].bot]) {
							//Check if FHB has valid past outcome results
							rev_valid_outcome = pred->fhb[0].o[pred->fhb[0].bot];
						}
					} 
					
					if ((rev_valid_outcome==NULL) && (pred->class==BPredMBP)||(pred->class==BPredTSCL)||(pred->class==BPredLLBP)||
					(pred->class==BPredOHT)) {
						//Check if OHT has valid past outcome results
						int key ;
						
						if ((pred->class==BPredMBP)||(pred->class==BPredOHT))
							key = key_from_features (pred->rev_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
						else 
							key = key_from_features (pred->rev_tage_dirpred[pred->tage_depth-1].twolev, rbaddr); // Get unmasked key from GHR and PC
					
						key = key & (pred->rev_dirpred.oht->oht.size - 1); // mask key based on predictor table size
						
						if (flow_mode && pred->rev_dirpred.oht->oht.rev_valid[key]) {
							rev_valid_outcome = pred->rev_dirpred.oht->oht.oc[key];
						}
					} 
						
					if ((pred->class==BPredMBP)||(pred->class==BPredOB)||(pred->class==BPredOHT)||(pred->class==BPred2Level)) {
						dir_update_ptr->rev_pdir1 = bpred_dir_lookup(pred, pred->rev_dirpred.twolev, rbaddr, 1, frmt);
					} else {
						char *bimod, *tage, *loop, *llbt;
						int tage_match = 0;
						int loop_match = 0;
						int llbt_match = 0;
						
						bimod = bpred_dir_lookup(pred, pred->rev_dirpred.bimod, rbaddr, 1, frmt);
						
						if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
							dir_update_ptr->rev_dir.sum = 0;
							
							if (pred->class==BPredLLBP) {
								for (int i = pred->tage_depth-1; i >= 1; i--) {
									// Get LLBT tag and compare to current tag
									int llbt_key = key_from_features (pred->rev_pb_dirpred[i].twolev, rbaddr); // Get unmasked key from GHR and PC
									
									// LLBP is indexed by RCR CCID!
									// If matched unmasked TAGE key, set tage pointer and note it and hist length
									if (llbt_key==pred->rev_pb_dirpred[i].twolev->config.two.tag[pred->rcr.ccid] && 
									pred->rcr.ccid==pred->rev_pb_dirpred[i].twolev->config.two.context[pred->rcr.ccid]) {
										llbt = bpred_dir_lookup(pred, pred->rev_pb_dirpred[i].twolev, rbaddr, 1, frmt);
										llbt_match = pred->rev_pb_dirpred[i].twolev->config.two.hist_length[pred->rcr.ccid];
										break;	// Break at first (highest) match
									}
								}
							}
							
							for (int i = pred->tage_depth-1; i >= 1; i--) {
								// Get tage and sc tags and compare to current tag
								int tage_key = key_from_features (pred->rev_tage_dirpred[i].twolev, rbaddr); // Get unmasked key from GHR and PC
								int masked_tage_key = tage_key & (pred->rev_tage_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
								
								// If matched unmasked TAGE key, set tage pointer
								if (tage_key==pred->rev_tage_dirpred[i].twolev->config.two.tag[masked_tage_key]) {
									tage = bpred_dir_lookup(pred, pred->rev_tage_dirpred[i].twolev, rbaddr, 1, frmt);
									tage_match = pred->rev_tage_dirpred[i].twolev->config.two.shift_width;
									break;	// Break at first (highest) match
								}
							}
							
							for (int i = 1; i < pred->tage_depth; i++) {
								int sc_key = key_from_features (pred->rev_sc_dirpred[i].twolev, rbaddr); // Get unmasked key from GHR and PC
								int masked_sc_key = sc_key & (pred->rev_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
								
								if (sc_key==pred->rev_sc_dirpred[i].twolev->config.two.tag[masked_sc_key]) {
									dir_update_ptr->rev_dir.sum = dir_update_ptr->rev_dir.sum + (pred->rev_sc_dirpred[i].twolev->config.two.l2table[masked_sc_key] >= 2); 
								}
							}
							
							// Get loop tag and compare to current tag
							int loop_key = key_from_features (pred->rev_loop_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
							int masked_loop_key = loop_key & (pred->rev_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
							
							if (loop_key==pred->rev_loop_dirpred.twolev->config.two.tag[masked_loop_key]) {
								loop = bpred_dir_lookup(pred, pred->rev_loop_dirpred.twolev, rbaddr, 1, frmt);
								
								if (*loop >= pred->rev_loop_dirpred.twolev->config.two.threshold)
									loop_match = 1;
							}
						}
						
						if (loop_match) {
							dir_update_ptr->rev_pdir1 = loop;
						} else {
							if (llbt_match>tage_match) {
								dir_update_ptr->rev_pdir1 = llbt;
							} else if (tage_match) {
								dir_update_ptr->rev_pdir1 = tage;
							} else {
								dir_update_ptr->rev_pdir1 = bimod;
							}
							
							// Check for statistical correction
							if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
								dir_update_ptr->rev_dir.sum = dir_update_ptr->rev_dir.sum + (*(dir_update_ptr->rev_pdir1) >= 2);
								
								// invert on overthreshold summation
								if (dir_update_ptr->rev_dir.sum >= pred->rev_sc_dirpred[pred->tage_depth-1].twolev->config.two.threshold)
									invert = TRUE;
							}
						}
						
						// Save TAGE prediction as a valid outcome could have been used.
						if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
							dir_update_ptr->rev_tage_pred = (*(dir_update_ptr->rev_pdir1) >= 2);
						}
					}
				}
				break;
				
			case BPredTaken:
				return rbtarget;
				
			case BPredNotTaken:
				if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) {
					return rbaddr + sizeof(md_inst_t);
				} else {
					return rbtarget;
				}
				
			default:
				panic("bogus predictor class");
		}
	}
	
   /*
	* We have a stateful predictor, and have gotten a pointer into the
	* direction predictor (except for jumps, for which the ptr is null)
	*
	* record pre-pop TOS; if this branch is executed speculatively
	* and is squashed, we'll restore the TOS and hope the data
	* wasn't corrupted in the meantime.
	*
	* For REV mode, return requires push and call requires pop 
	*
	*/
	
	if (!flow_mode) {
		if (pred->retstack.size)
			*stack_recover_idx = pred->retstack.tos;
		else
			*stack_recover_idx = 0;

		/* if this is a return, pop return-address stack */
		if (is_return && pred->retstack.size) {
			md_addr_t target = pred->retstack.stack[pred->retstack.tos].target;
			
			pred->retstack.tos = (pred->retstack.tos + pred->retstack.size - 1) % pred->retstack.size;
			pred->retstack_pops++;
			
			dir_update_ptr->fwd_dir.ras = TRUE; /* using RAS here */
			
			return target;
		}

#ifndef RAS_BUG_COMPATIBLE
		/* if function call, push return-address onto return-address stack */
		if (is_call && pred->retstack.size) {
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = fbaddr + sizeof(md_inst_t);
			pred->retstack_pushes++;
		}
#endif /* !RAS_BUG_COMPATIBLE */
	} else {
		if (pred->retstack.size)
			*stack_recover_idx = pred->retstack.tos;
		else
			*stack_recover_idx = 0;

		/* Since this is reversed, return constitutes a push and call a pop  */
		if (is_return && pred->retstack.size) {
#ifndef RAS_BUG_COMPATIBLE
			/* if this is a return, push rbtarget onto return-address stack */
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = rbtarget;
			
			pred->reverse_retstack_pushes++;
#endif /* !RAS_BUG_COMPATIBLE */
			// Still need to track this like in FWD mode until ISA changes
			dir_update_ptr->rev_dir.ras = TRUE; /* using RAS here */

			return rbtarget;
		}
		
		/* if function call, pop return-address stack */
		if (is_call && pred->retstack.size) {
			// Don't need target for this as there is no reversible ISA yet implemented
			//md_addr_t target = pred->retstack.stack[pred->retstack.tos].target;
			
			pred->retstack.tos = (pred->retstack.tos + pred->retstack.size - 1) % pred->retstack.size;
			pred->reverse_retstack_pops++;
			
			//dir_update_ptr->rev_dir.ras = TRUE; /* using RAS here */
			
			//return target;
		}
	}
  
	/* Get a pointer into the BTB early as we may need target and it's already set after FWD mode */
	index = (baddr >> MD_BR_SHIFT) & (pred->btb.sets - 1);

	if (pred->btb.assoc > 1) {
		index *= pred->btb.assoc;

		/* Now we know the set; look for a PC match */
		for (i = index; i < (index+pred->btb.assoc) ; i++)
			if (pred->btb.btb_data[i].addr == baddr) {
				/* match */
				pbtb = &pred->btb.btb_data[i];
				break;
			}
	} else {
		pbtb = &pred->btb.btb_data[index];
		
		if (pbtb->addr != baddr)
			pbtb = NULL;
	}

  /*
   * We now also have a pointer into the BTB for a hit, or NULL otherwise
   */

	/* if this is a jump, ignore predicted direction; we know it's taken. */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) == (F_CTRL|F_UNCOND)) {
		return (pbtb ? pbtb->target : 1);
	}
	
	int prediction;
		
	/* otherwise we have a conditional branch */
	if (!flow_mode) {
		// There was a fwd valid outcome
		if (fwd_valid_outcome != NULL) {
			prediction = fwd_valid_outcome;
		} else {
			prediction = (*(dir_update_ptr->fwd_pdir1) >= 2);
			
			if (invert) {
				prediction = !prediction;
			}
		}
	} else {
		// There was a rev valid outcome
		if (rev_valid_outcome != NULL) {
			prediction = rev_valid_outcome;
		} else {
			prediction = (*(dir_update_ptr->rev_pdir1) >= 2);
			
			if (invert) {
				prediction = !prediction;
			}
		}
	}
	  
	md_addr_t ret_target; 
	
	if (pbtb == NULL) {
		/* BTB miss -- just return a predicted direction */
		ret_target = 1;
	} else {
		/* BTB hit, so return target if it's a predicted-taken branch */
		ret_target = pbtb->target;
	}
	
	return (prediction
		? /* taken */ ret_target
		: /* not taken */ 0
	);
}

/* Speculative execution can corrupt the ret-addr stack.  So for each
 * lookup we return the top-of-stack (TOS) at that point; a mispredicted
 * branch, as part of its recovery, restores the TOS using this value --
 * hopefully this uncorrupts the stack. */
void
bpred_recover(struct bpred_t *pred,	/* branch predictor instance */
	      md_addr_t baddr,		/* branch address */
	      int stack_recover_idx)	/* Non-speculative top-of-stack;
					 * used on mispredict recovery */
{
  if (pred == NULL)
    return;

  pred->retstack.tos = stack_recover_idx;
}

// Fetch LLBP row to PB
void bpred_llbp_fetch(struct bpred_t *pred,	/* branch predictor instance */
			struct bpred_dirpred_t *pb_pred_dirpred,	/* branch dir predictor inst */
			struct bpred_dirpred_t *llbt_pred_dirpred,	/* branch dir predictor inst */
			int ccid)				/* CCID */
{
	for (int i = pred->tage_depth-1; i >= 1; i--) {
		// LLBT indexed by RCR CCID!!	
		// Check to make sure contexts match
		if (ccid!=pb_pred_dirpred[i].twolev->config.two.context[ccid]) { // In this case, they don't
			// Row now has to be replaced
			// First evict current row
			int victim_context = pb_pred_dirpred[i].twolev->config.two.context[ccid];
			int victim_tag = pb_pred_dirpred[i].twolev->config.two.tag[ccid];
			int victim_hist_length = pb_pred_dirpred[i].twolev->config.two.hist_length[ccid];
			int victim_count = pb_pred_dirpred[i].twolev->config.two.l2table[ccid];
			int victim_use = pb_pred_dirpred[i].twolev->config.two.use[ccid];
			
			// Next, pull in row from LLBT
			pb_pred_dirpred[i].twolev->config.two.context[ccid] = llbt_pred_dirpred[i].twolev->config.two.context[ccid];
			pb_pred_dirpred[i].twolev->config.two.tag[ccid] = llbt_pred_dirpred[i].twolev->config.two.tag[ccid];
			pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] = llbt_pred_dirpred[i].twolev->config.two.hist_length[ccid];
			pb_pred_dirpred[i].twolev->config.two.l2table[ccid] = llbt_pred_dirpred[i].twolev->config.two.l2table[ccid];
			pb_pred_dirpred[i].twolev->config.two.use[ccid] = llbt_pred_dirpred[i].twolev->config.two.use[ccid];
			
			// Lastly, push victim row to LLBT
			llbt_pred_dirpred[i].twolev->config.two.context[victim_context] = victim_context;
			llbt_pred_dirpred[i].twolev->config.two.tag[victim_context] = victim_tag;
			llbt_pred_dirpred[i].twolev->config.two.hist_length[victim_context] = victim_hist_length;
			llbt_pred_dirpred[i].twolev->config.two.l2table[victim_context] = victim_count;
			llbt_pred_dirpred[i].twolev->config.two.use[victim_context] = victim_use;
		}
	}
}

// Update LLBP based on outcomes and pre-/post-fetch
void bpred_llbp_update(struct bpred_t *pred,	/* branch predictor instance */
			struct bpred_dirpred_t *pb_pred_dirpred,	/* branch dir predictor inst */
			struct bpred_dirpred_t *llbt_pred_dirpred,	/* branch dir predictor inst */
			int ccid,				/* CCID */
			md_addr_t addr,			/* addr to get keys */
			int *replacement_max,
			int* *tage_replacement_counters,				/* replacement pointers */
			int* *tage_replacement_tags,				/* replacement pointers */
			int* *tage_replacement_lengths,				/* replacement pointers */
			int taken,					/* branch outcome */
			int flow_mode,
			int update_history)
{
	// Get LLBT tag(s) and compare to current PB tag(s)
	for (int i = pred->tage_depth-1; i >= 1; i--) {
		// Get unmasked PB tag from GHR and PC
		// Need to use a custom method for LLBP's
		int pb_tag = key_from_llbp_features (pb_pred_dirpred[i].twolev, addr, ccid); 
		
		// LLBT indexed by RCR CCID!!	
		// Check to make sure contexts match
		if (ccid!=pb_pred_dirpred[i].twolev->config.two.context[ccid]) { // In this case, they don't
			// Row now has to be replaced
			// First evict current row
			int victim_context = pb_pred_dirpred[i].twolev->config.two.context[ccid];
			int victim_tag = pb_pred_dirpred[i].twolev->config.two.tag[ccid];
			int victim_hist_length = pb_pred_dirpred[i].twolev->config.two.hist_length[ccid];
			int victim_count = pb_pred_dirpred[i].twolev->config.two.l2table[ccid];
			int victim_use = pb_pred_dirpred[i].twolev->config.two.use[ccid];
			
			// Next, pull in row from LLBT
			pb_pred_dirpred[i].twolev->config.two.context[ccid] = llbt_pred_dirpred[i].twolev->config.two.context[ccid];
			pb_pred_dirpred[i].twolev->config.two.tag[ccid] = llbt_pred_dirpred[i].twolev->config.two.tag[ccid];
			pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] = llbt_pred_dirpred[i].twolev->config.two.hist_length[ccid];
			pb_pred_dirpred[i].twolev->config.two.l2table[ccid] = llbt_pred_dirpred[i].twolev->config.two.l2table[ccid];
			pb_pred_dirpred[i].twolev->config.two.use[ccid] = llbt_pred_dirpred[i].twolev->config.two.use[ccid];
			
			// Lastly, push victim row to LLBT
			llbt_pred_dirpred[i].twolev->config.two.context[victim_context] = victim_context;
			llbt_pred_dirpred[i].twolev->config.two.tag[victim_context] = victim_tag;
			llbt_pred_dirpred[i].twolev->config.two.hist_length[victim_context] = victim_hist_length;
			llbt_pred_dirpred[i].twolev->config.two.l2table[victim_context] = victim_count;
			llbt_pred_dirpred[i].twolev->config.two.use[victim_context] = victim_use;
		} 
		
		if (pb_tag!=pb_pred_dirpred[i].twolev->config.two.tag[ccid]) {
			// If didn't match unmasked LLBT key, decrement useful counter
			// If usefulness counter is not at 0, then decrement
			if (!pb_pred_dirpred[i].twolev->config.two.use[ccid]) {
				pb_pred_dirpred[i].twolev->config.two.use[ccid]--;
			}
			
			// If now a replacement candidate, replace with longer history tage victim if the masked keys match
			if (!pb_pred_dirpred[i].twolev->config.two.use[ccid]) {
				// check if any replacement candidates from TAGE; go with largest first and swap replacements
				for (int j = *replacement_max-1; j >= 0; j--) {
					if ((*tage_replacement_lengths[j] > pb_pred_dirpred[i].twolev->config.two.hist_length[ccid])) {
						// Temp vars for swapping
						int temp_tag = pb_pred_dirpred[i].twolev->config.two.tag[ccid];
						int temp_hist = pb_pred_dirpred[i].twolev->config.two.hist_length[ccid];
						int temp_count = pb_pred_dirpred[i].twolev->config.two.l2table[ccid];
						
						pb_pred_dirpred[i].twolev->config.two.tag[ccid] = *tage_replacement_tags[j];
						pb_pred_dirpred[i].twolev->config.two.use[ccid] = 2; 	//weakly set use
						pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] = *tage_replacement_lengths[j];
						pb_pred_dirpred[i].twolev->config.two.l2table[ccid] = *tage_replacement_counters[j];
						
						*tage_replacement_tags[j] = temp_tag;
						*tage_replacement_lengths[j] = temp_hist;
						*tage_replacement_counters[j] = temp_count;
						
						break;
					}
				}
			}
		} else {
			//increment usefulness unless at 3 already
			if (pb_pred_dirpred[i].twolev->config.two.use[ccid]<3) {
				pb_pred_dirpred[i].twolev->config.two.use[ccid]++;
			}
			
			if (taken) {
				if (pb_pred_dirpred[i].twolev->config.two.l2table[ccid] < 3)
					pb_pred_dirpred[i].twolev->config.two.l2table[ccid]++;
			} else {
				if (pb_pred_dirpred[i].twolev->config.two.l2table[ccid])
					pb_pred_dirpred[i].twolev->config.two.l2table[ccid]--;
			}
		}
		
		if (update_history) {
			if (!flow_mode) { // FWD mode
				shift_history_left(pb_pred_dirpred[i].twolev, addr, taken);
			} else { // REV mode
				shift_history_right(pb_pred_dirpred[i].twolev, addr, taken);
			}
		}
	}
	
	int needs_sorting = 1;
	
	//Sort LLBTs 
	// There aren't a lot so simple bubble sort is fine
	while (needs_sorting) {
		int sorted = 0;
			
		for (int i = 1; i < pred->tage_depth - 1; i++) {
			if (pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] > pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid]) {
				// Temp vars for sorting
				int temp_tag = pb_pred_dirpred[i].twolev->config.two.tag[ccid];
				int temp_hist = pb_pred_dirpred[i].twolev->config.two.hist_length[ccid];
				int temp_use = pb_pred_dirpred[i].twolev->config.two.hist_length[ccid];
				int temp_count = pb_pred_dirpred[i].twolev->config.two.hist_length[ccid];
				
				// start sort
				pb_pred_dirpred[i].twolev->config.two.tag[ccid] = pb_pred_dirpred[i+1].twolev->config.two.tag[ccid];
				pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] = pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid];
				pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] = pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid];
				pb_pred_dirpred[i].twolev->config.two.hist_length[ccid] = pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid];
				
				// complete sort
				pb_pred_dirpred[i+1].twolev->config.two.tag[ccid] = temp_tag;
				pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid] = temp_hist;
				pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid] = temp_use;
				pb_pred_dirpred[i+1].twolev->config.two.hist_length[ccid] = temp_count;
				
				sorted = 1;
			}
		}
			
		if (!sorted)	// No more sorting was done so it doesn't need it anymore
			needs_sorting = 0;
	}
}

// Update Loop table based on outcome
void bpred_loop_update(struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
	int loop_key,
	int unmasked_loop_key,
	int taken)
{
	if (unmasked_loop_key==pred_dir->config.two.tag[loop_key]) {
		if (taken) {	// Branch was taken
			pred_dir->config.two.iter_c[loop_key]++;
		} else {	// Branch was not taken so moving iteration counter to past
			pred_dir->config.two.iter_p[loop_key] = pred_dir->config.two.iter_c[loop_key];
			pred_dir->config.two.iter_c[loop_key] = 0;
		}
		
		// hit the same number of past and current iterations
		if (pred_dir->config.two.iter_c[loop_key]==pred_dir->config.two.iter_p[loop_key]) {
			if (pred_dir->config.two.l2table[loop_key] < 3)
				pred_dir->config.two.l2table[loop_key]++; // incr counter
		}
		
		if (pred_dir->config.two.use[loop_key] < 3)
			pred_dir->config.two.use[loop_key]++;
	} else if (unmasked_loop_key!=pred_dir->config.two.tag[loop_key]) {
		if (pred_dir->config.two.use[loop_key])
			pred_dir->config.two.use[loop_key]--;
	}
	
	if (!pred_dir->config.two.use[loop_key]) {
		pred_dir->config.two.l2table[loop_key] = 0; //strongly unset counter
		pred_dir->config.two.tag[loop_key] = unmasked_loop_key; //overwrite tag
		pred_dir->config.two.use[loop_key] = 0; //set strongly not useful
	}
}

// Update SC counters based on correctness 
void bpred_sc_update(struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
	int sc_key,
	int unmasked_sc_key,
	int correct,		// correctness flag
	int und,			// under threshold flag
	int thu)			// threshold update
{
	//correct = (!correct ^ und);		// Invert correction if over threshold to get tage correction
	
	if (unmasked_sc_key==pred_dir->config.two.tag[sc_key]) {
		if (!correct && und) {
			if (pred_dir->config.two.l2table[sc_key] < 3)
				pred_dir->config.two.l2table[sc_key]++;
		}
		
		if (pred_dir->config.two.use[sc_key] < 3)
			pred_dir->config.two.use[sc_key]++;
	} else if (unmasked_sc_key!=pred_dir->config.two.tag[sc_key]) {
		if (pred_dir->config.two.use[sc_key])
			pred_dir->config.two.use[sc_key]--;
	}
	
	if (!pred_dir->config.two.use[sc_key]) {
		pred_dir->config.two.l2table[sc_key] = 0; //strongly don't set counter
		pred_dir->config.two.tag[sc_key] = unmasked_sc_key; //overwrite tag
		pred_dir->config.two.use[sc_key] = 2; //set weakly useful
	}
}

// Update TAGE counters based on dirpred pointer
void bpred_tage_update (struct bpred_t *pred,	/* branch predictor inst */
	int *replacement_size,
	int *replacement_max,
	int* *tage_replacement_counters,				/* replacement pointers */
	int* *tage_replacement_tags,				/* replacement pointers */
	int* *tage_replacement_lengths,				/* replacement pointers */
	md_addr_t baddr,			/* Branch address */
	int correct,				/* predictor correctness bit */
	int taken,					/* Branch Outcome */
	int flow_mode,			/* flow mode flag; FWD=0,REV=1 */
	int future_mode,			/* future mode enabled; changes addr, correct, and taken to FHB values */
	int frmt)					/* FRMTs enabled flag */
{
	int matching_key = 0;
	int matching_table = 0;
	
	int altpred_key = 0;
	int altpred_table = 0;
	
	char *matching_pred, *altpred, *bimod_pred;
	
	struct bpred_dir_t *pred_dir;	/* branch dir predictor inst */
	int unmasked_tage_key, tage_key;	
	
	int* this_outcome = (int*)malloc(pred->tage_depth * sizeof(int));
	int* this_correct = (int*)malloc(pred->tage_depth * sizeof(int));
	int* this_addr = (int*)malloc(pred->tage_depth * sizeof(int));
	
	// Get bimod pred just in case it is the altpred
	if (flow_mode && !frmt) {	//bimod is rev bimod
		bimod_pred = bpred_dir_lookup(pred, pred->rev_dirpred.bimod, baddr, flow_mode, frmt);
	} else {
		bimod_pred = bpred_dir_lookup(pred, pred->fwd_dirpred.bimod, baddr, flow_mode, frmt);
	}
	
	// Iterate through each table to find a match
	for (int i = 1; i < pred->tage_depth; i++) {
		if (flow_mode && !frmt) {
			pred_dir = pred->rev_tage_dirpred[i].twolev;
		} else {
			pred_dir = pred->fwd_tage_dirpred[i].twolev;
		}
		
		// Get tage keys
		if (future_mode) {	// NEED TO ACCOUNT FOR FUTURE FHB VALUES !!!!!!!
			unmasked_tage_key = future_key_from_tage(pred, pred_dir, baddr, i, flow_mode);
		} else {
			unmasked_tage_key = key_from_features (pred_dir, baddr); // Get unmasked key from GHR and PC
		}
		tage_key = unmasked_tage_key & (pred_dir->config.two.l2size - 1); // mask key based on predictor table size
		
		if (unmasked_tage_key==pred_dir->config.two.tag[tage_key]) {
			if (matching_key) {	// Set altpred if already matched at lower table
				altpred_key = matching_key;
				altpred_table = matching_table;
				altpred = matching_pred;
			} 
			
			matching_key = tage_key;
			matching_table = i;
			matching_pred = bpred_dir_lookup(pred, pred_dir, baddr, flow_mode, frmt);
		}
	}
	
	if (!altpred_key) { // If no altpred, then bimod is the altpred
		altpred = bimod_pred;
	}
	
	if (!matching_key) {	// If not match, will set bimod as predictor
		matching_pred = bimod_pred;
	}
	
	//Update useful counter only on altpred!=pred and matching_pred is not bimod_pred
	if ((*matching_pred >= 2)!=(*altpred >= 2) && matching_key) {
		if (flow_mode && !frmt) {
			pred_dir = pred->rev_tage_dirpred[matching_table].twolev;
		} else {
			pred_dir = pred->fwd_tage_dirpred[matching_table].twolev;
		}
		
		if (correct) {
			if (pred_dir->config.two.use[matching_key] < 3)
				pred_dir->config.two.use[matching_key]++;
		} else {
			if (pred_dir->config.two.use[matching_key])
				pred_dir->config.two.use[matching_key]--;
		}
	}
	
	// Update providing component pred counter
	// Only if not bimod which is updated on it's own.
	if (matching_key) {
		if (taken) {
			if (*matching_pred < 3) {
				++*matching_pred;
			}
		} else {
			if (*matching_pred) {
				--*matching_pred;
			}
		}
	}
	
	// If incorrect, reallocate longer history entry if applicable 
	// Based on usefulness
	// Using policy defined in https://inria.hal.science/hal-03408381/document
	if (!correct && matching_table!=(pred->tage_depth-1)) {
		// viable entries j & k
		int j = 0;
		int k = 0;	
		
		// First find viable entries
		for (int i = matching_table + 1; i < pred->tage_depth; i ++) {
			if (flow_mode && !frmt) {
				pred_dir = pred->rev_tage_dirpred[i].twolev;
			} else {
				pred_dir = pred->fwd_tage_dirpred[i].twolev;
			}
			
			// Get tage keys
			if (future_mode) {	// NEED TO ACCOUNT FOR FUTURE FHB VALUES !!!!!!!
				unmasked_tage_key = future_key_from_tage(pred, pred_dir, baddr, i, flow_mode);
			} else {
				unmasked_tage_key = key_from_features (pred_dir, baddr); // Get unmasked key from GHR and PC
			}
			tage_key = unmasked_tage_key & (pred_dir->config.two.l2size - 1); // mask key based on predictor table size
			
			if (!pred_dir->config.two.use[tage_key]) {
				if (k) {
					j = k;
				}
				
				k = i;
			}
		}
		
		if (j) {	// j is set so there exists some additional option
			// per Seznec and Michaud, j should have twice the probability of being set
			
			int prob = rand() % 100;
			
			if (prob > 33) {
				k = j;
			}
		}
		
		if (k) {
			if (flow_mode && !frmt) {
				pred_dir = pred->rev_tage_dirpred[k].twolev;
			} else {
				pred_dir = pred->fwd_tage_dirpred[k].twolev;
			}
			
			// Get tage keys
			if (future_mode) {	// NEED TO ACCOUNT FOR FUTURE FHB VALUES !!!!!!!
				unmasked_tage_key = future_key_from_tage(pred, pred_dir, baddr, k, flow_mode);
			} else {
				unmasked_tage_key = key_from_features (pred_dir, baddr); // Get unmasked key from GHR and PC
			}
			tage_key = unmasked_tage_key & (pred_dir->config.two.l2size - 1); // mask key based on predictor table size
			
			if (*replacement_size >= *replacement_max) {
				*replacement_max++;
				*tage_replacement_tags = (int*)realloc(*tage_replacement_tags, *replacement_max * sizeof(int));
				*tage_replacement_counters = (int*)realloc(*tage_replacement_counters, *replacement_max * sizeof(int));
				*tage_replacement_lengths = (int*)realloc(*tage_replacement_lengths, *replacement_max * sizeof(int));
			}
			
			*tage_replacement_counters = pred_dir->config.two.l2table[tage_key]; //copy current counter
			*tage_replacement_tags = pred_dir->config.two.tag[tage_key]; //copy current tag
			*tage_replacement_lengths = pred_dir->config.two.shift_width; //copy hist_length
			*replacement_size++;
			
			pred_dir->config.two.l2table[tage_key] = 2; //weakly set counter
			pred_dir->config.two.tag[tage_key] = unmasked_tage_key; //overwrite tag
			pred_dir->config.two.use[tage_key] = 0; //strongly not useful
		} else {	// If no options exist, decrement usefulness of all Tj i < j < M
			for (int i = matching_table + 1; i < pred->tage_depth; i ++) {
				if (flow_mode && !frmt) {
					pred_dir = pred->rev_tage_dirpred[i].twolev;
				} else {
					pred_dir = pred->fwd_tage_dirpred[i].twolev;
				}
				
				// Get tage keys
				if (future_mode) {	// NEED TO ACCOUNT FOR FUTURE FHB VALUES !!!!!!!
					unmasked_tage_key = future_key_from_tage(pred, pred_dir, baddr, i, flow_mode);
				} else {
					unmasked_tage_key = key_from_features (pred_dir, baddr); // Get unmasked key from GHR and PC
				}
				tage_key = unmasked_tage_key & (pred_dir->config.two.l2size - 1); // mask key based on predictor table size
				
				if (pred_dir->config.two.use[tage_key])
					pred_dir->config.two.use[tage_key]--;
			}
		}
	}
}

// Shift L1 History Right
void shift_history_left(struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
			md_addr_t baddr,			/* addr to index into l1 table */
			int taken)				/* outcome */
{
	int l1index, shift_reg;		
	l1index = (baddr >> MD_BR_SHIFT) & (pred_dir->config.two.l1size - 1);
	shift_reg = (pred_dir->config.two.shiftregs[l1index] << 1) | (!!taken);
	pred_dir->config.two.shiftregs[l1index] = shift_reg & ((1 << pred_dir->config.two.shift_width) - 1);
}

// Shift L1 History Left
void shift_history_right(struct bpred_dir_t *pred_dir,	/* branch dir predictor inst */
			md_addr_t baddr,			/* addr to index into l1 table */
			int taken)				/* outcome */
{
	int l1index, shift_reg;		
	l1index = (baddr >> MD_BR_SHIFT) & (pred_dir->config.two.l1size - 1);
	shift_reg = (pred_dir->config.two.shiftregs[l1index] >> 1) | ((!!taken) << (pred_dir->config.two.shift_width - 1));
	pred_dir->config.two.shiftregs[l1index] = shift_reg & ((1 << pred_dir->config.two.shift_width) - 1);
}

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
	     int correct,		/* was earlier addr prediction ok? */
	     enum md_opcode op,		/* opcode of instruction */
	     struct bpred_update_t *dir_update_ptr,  /* FWD pred state pointer */
		 int flow_mode,  /* Flow mode (0=FWD; 1=REV) */
		 int frmt)		 /* FRMTs Enabled flag */
{
	struct bpred_btb_ent_t *pbtb = NULL;
	struct bpred_btb_ent_t *lruhead = NULL;
	struct bpred_btb_ent_t *lruitem = NULL;
	int index, i;
	
	// FWD & REV addresses; 
	// These are supposed to be opposites but 
	// without ISA updates and knowledge of addrs and targets, 
	// these need to be the same.
	md_addr_t fbaddr = baddr;		/* FWD Branch ADDR */		
	md_addr_t fbtarget = btarget;	/* FWD Branch Target */
	md_addr_t rbaddr = baddr;		/* REV Branch ADDR; for now, setting to baddr */
	md_addr_t rbtarget = btarget;	/* REV Branch Target; for now, setting to btarget */
	
	// Need to update FRMT ADDR before overwriting 
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (!flow_mode) && (frmt) && 
	((pred->class == BPred2bit) || (pred->class == BPred2Level) || (pred->class == BPredOB) || 
	(pred->class == BPredOHT) || (pred->class == BPredMBP) || (pred->class == BPredTSCL) || 
	(pred->class == BPredLLBP))) {
		switch (pred->class) {
			case BPredLLBP:
			case BPredComb:
			case BPred2Level:
			case BPredTSBP:
			case BPredCHBP:
			case BPredOB:
			case BPredOHT:
			case BPredMBP:
			case BPredTSCL:
			case BPred2bit:
				//Need to map FWD to REV addr
				map_addr_to_frmt(pred, fbaddr, rbaddr);
			
			case BPredTaken:
			case BPredNotTaken:
				break;
			
			/* no other state */
			default:
				panic("bogus predictor class");
		}
	}
	
	if (frmt)	// FRMT was enabled so REV ADDR needs to be exchanged for FWD ADDR
		rbaddr = addr_from_frmt(pred, rbaddr);	
  
   /* don't change bpred state for non-branch instructions or if this
	* is a stateless predictor*/
	if (!(MD_OP_FLAGS(op) & F_CTRL))
		return;

	/* Have a branch here */

	/* For UNCOND branches using LLBP, need to update the RCR */
	if (((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) == (F_CTRL|F_UNCOND)) && (pred->class == BPredLLBP)){
		//Update RCR and shift up or down depending on flow mode
		if (!flow_mode) {
			// We don't care about saving the address after it's gone off
			//rcr_addr = pred->rcr.addr[pred->rcr.top];
			
			if (pred->rcr.top > 0) {
				pred->rcr.top--;
			} else {
				pred->rcr.top = pred->rcr.size - 1;
			}
			
			if (pred->rcr.win_top > 0) {
				pred->rcr.win_top--;
			} else {
				pred->rcr.win_top = pred->rcr.size - 1;
			}
			
			if (pred->rcr.fpre_top > 0) {
				pred->rcr.fpre_top--;
			} else {
				pred->rcr.fpre_top = pred->rcr.size - 1;
			}
			
			if (pred->rcr.rpre_bot > 0) {
				pred->rcr.rpre_bot--;
			} else {
				pred->rcr.rpre_bot = pred->rcr.size - 1;
			}
			
			if (pred->rcr.win_bot > 0) {
				pred->rcr.win_bot--;
			} else {
				pred->rcr.win_bot = pred->rcr.size - 1;
			}
			
			if (pred->rcr.bot > 0) {
				pred->rcr.bot--;
			} else {
				pred->rcr.bot = pred->rcr.size - 1;
			}
			
			if (pbtb) {
				pred->rcr.addr[pred->rcr.bot] = pbtb->target;
			} else {
				pred->rcr.addr[pred->rcr.bot] = NULL;
			}
			
			pred->rcr.addr[pred->rcr.bot] = fbaddr;
		} else {
			// We don't care about saving the address after it's gone off
			//rcr_addr = pred->rcr.addr[pred->rcr.bot];
			
			if (pred->rcr.top < (pred->rcr.size - 1)) {
				pred->rcr.top++;
			} else {
				pred->rcr.top = 0;
			}
			
			if (pred->rcr.win_top < (pred->rcr.size - 1)) {
				pred->rcr.win_top++;
			} else {
				pred->rcr.win_top = 0;
			}
			
			if (pred->rcr.fpre_top < (pred->rcr.size - 1)) {
				pred->rcr.fpre_top++;
			} else {
				pred->rcr.fpre_top = 0;
			}
			
			if (pred->rcr.rpre_bot < (pred->rcr.size - 1)) {
				pred->rcr.rpre_bot++;
			} else {
				pred->rcr.rpre_bot = 0;
			}
			
			if (pred->rcr.win_bot < (pred->rcr.size - 1)) {
				pred->rcr.win_bot++;
			} else {
				pred->rcr.win_bot = 0;
			}
			
			if (pred->rcr.bot < (pred->rcr.size - 1)) {
				pred->rcr.bot++;
			} else {
				pred->rcr.bot = 0;
			}
			
			if (pbtb) {
				pred->rcr.addr[pred->rcr.top] = pbtb->addr;
			} else {
				pred->rcr.addr[pred->rcr.top] = NULL;
			}
			
			pred->rcr.addr[pred->rcr.top] = rbaddr;
		}
		
		// Update CCIDs
		pred->rcr.fccid = 0;
		pred->rcr.ccid = 0;
		pred->rcr.rccid = 0;
		
		for (int i = 0; i < pred->rcr.size - 4; i++) {
			pred->rcr.fccid = pred->rcr.fccid ^ pred->rcr.addr[pred->rcr.bot + i];
			pred->rcr.ccid = pred->rcr.ccid ^ pred->rcr.addr[pred->rcr.win_bot + i];
			pred->rcr.rccid = pred->rcr.rccid ^ pred->rcr.addr[pred->rcr.rpre_bot + i];
		}
	}

	/* Now accounting for FWD or REV stats */
	if (!flow_mode) {
		if (correct)
			pred->addr_hits++;

		if (!!pred_taken == !!taken)
			pred->dir_hits++;
		else
			pred->misses++;

		if (dir_update_ptr->fwd_dir.ras) {
			pred->used_ras++;
			
			if (correct)
				pred->ras_hits++;
		} else if ((MD_OP_FLAGS(op) & (F_CTRL|F_COND)) == (F_CTRL|F_COND)) {
			if (dir_update_ptr->fwd_dir.meta)
				pred->used_2lev++;
			else
				pred->used_bimod++;
		}
		
		/* keep stats about JR's; also, but don't change any bpred state for JR's
		* which are returns unless there's no retstack */
		if (MD_IS_INDIR(op)) {
			pred->jr_seen++;
			if (correct)
				pred->jr_hits++;
		  
			if (!dir_update_ptr->fwd_dir.ras) {
				pred->jr_non_ras_seen++;
				
				if (correct)
					pred->jr_non_ras_hits++;
			} else {
				/* return that used the ret-addr stack; no further work to do */
				return;
			}
		}
	} else {
		if (correct)
			pred->reverse_addr_hits++;

		if (!!pred_taken == !!taken)
			pred->reverse_dir_hits++;
		else
			pred->reverse_misses++;

		if (dir_update_ptr->rev_dir.ras) {
			pred->reverse_used_ras++;
			
			if (correct)
				pred->reverse_ras_hits++;
		} else if ((MD_OP_FLAGS(op) & (F_CTRL|F_COND)) == (F_CTRL|F_COND)) {
			if (dir_update_ptr->rev_dir.meta)
				pred->reverse_used_2lev++;
			else
				pred->reverse_used_bimod++;
		}
		
		/* keep stats about rev JR's (should be CALLs but will update with ISA changes later); also, but don't change any bpred state for JR's
		* which are returns unless there's no retstack */
		if (MD_IS_INDIR(op)) {
			pred->reverse_jr_seen++;
			
			if (correct)
				pred->reverse_jr_hits++;
		  
			if (!dir_update_ptr->rev_dir.ras) {
				pred->reverse_jr_non_ras_seen++;
				
				if (correct)
					pred->reverse_jr_non_ras_hits++;
			} else {
				/* return that used the ret-addr stack; no further work to do */
				return;
			}
		}
	}

	/* Can exit now if this is a stateless predictor */
	if (pred->class == BPredNotTaken || pred->class == BPredTaken)
		return;

   /* 
	* Now we know the branch didn't use the ret-addr stack, and that this
	* is a stateful predictor 
	*
	* Added handling for REV where return constitutes a push
	*/
 
#ifdef RAS_BUG_COMPATIBLE
	if (!flow_mode) {
		/* if function call, push return-address onto return-address stack */
		if (MD_IS_CALL(op) && pred->retstack.size) {
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = fbaddr + sizeof(md_inst_t);
		
			pred->retstack_pushes++;
		}
	} else {
		if (MD_IS_INDIR(op) && pred->retstack.size) {
			pred->retstack.tos = (pred->retstack.tos + 1)% pred->retstack.size;
			pred->retstack.stack[pred->retstack.tos].target = rbtarget;
		
			pred->reverse_retstack_pushes++;
		}
	}
#endif /* RAS_BUG_COMPATIBLE */

	// Updating BTB earlier to get target for rev pred updates
	
	/* find BTB entry if it's a taken branch (don't allocate for non-taken) */
	if (taken) {
		index = (baddr >> MD_BR_SHIFT) & (pred->btb.sets - 1);
      
		if (pred->btb.assoc > 1) {
			index *= pred->btb.assoc;
	  
		   /* Now we know the set; look for a PC match; also identify
			* MRU and LRU items */
			for (i = index; i < (index+pred->btb.assoc) ; i++) {
				if (pred->btb.btb_data[i].addr == baddr) {
					/* match */
					assert(!pbtb);
					pbtb = &pred->btb.btb_data[i];
				}
	      
				dassert(pred->btb.btb_data[i].prev != pred->btb.btb_data[i].next);
				if (pred->btb.btb_data[i].prev == NULL) {
					/* this is the head of the lru list, ie current MRU item */
					dassert(lruhead == NULL);
					lruhead = &pred->btb.btb_data[i];
				}
				
				if (pred->btb.btb_data[i].next == NULL) {
					/* this is the tail of the lru list, ie the LRU item */
					dassert(lruitem == NULL);
					lruitem = &pred->btb.btb_data[i];
				}
			}
			dassert(lruhead && lruitem);
	  
			if (!pbtb)
				/* missed in BTB; choose the LRU item in this set as the victim */
				pbtb = lruitem;	
				/* else hit, and pbtb points to matching BTB entry */
			  
		   /* Update LRU state: selected item, whether selected because it
			* matched or because it was LRU and selected as a victim, becomes 
			* MRU */
			if (pbtb != lruhead) {
				/* this splices out the matched entry... */
				if (pbtb->prev)
					pbtb->prev->next = pbtb->next;
				if (pbtb->next)
					pbtb->next->prev = pbtb->prev;
				/* ...and this puts the matched entry at the head of the list */
				pbtb->next = lruhead;
				pbtb->prev = NULL;
				lruhead->prev = pbtb;
				dassert(pbtb->prev || pbtb->next);
				dassert(pbtb->prev != pbtb->next);
			}
			/* else pbtb is already MRU item; do nothing */
		} else {
			pbtb = &pred->btb.btb_data[index];
		}
      
	   /* 
		* Now 'p' is a possibly null pointer into the direction prediction table, 
		* and 'pbtb' is a possibly null pointer into the BTB (either to a 
		* matched-on entry or a victim which was LRU in its set)
		*/
		   
		/* update BTB (but only for taken branches) */
		if (pbtb) {
			/* update current information */
			dassert(taken);

			if (pbtb->addr == baddr) {
				if (!correct)
					pbtb->target = btarget;
			} else {
				/* enter a new branch in the table */
				pbtb->addr = baddr;
				pbtb->op = op;
				pbtb->target = btarget;
			}
		}
	}

	int fwd_key, rev_key;
	
	/* Get keys before updating L1 table*/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && ((pred->class == BPredTSBP) || (pred->class == BPredCHBP) || (pred->class == BPredOHT) || (pred->class == BPredMBP))) {
		// Get keys before updating GHR!!!!
		fwd_key = key_from_features (pred->fwd_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
		
		// Get REV key from correct GHR based on FRMT setting
		if (!frmt)
			rev_key = key_from_features (pred->rev_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
		else 
			rev_key = key_from_features (pred->fwd_dirpred.twolev, rbaddr); // Get unmasked key from GHR and PC
	}
	
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && 
	((pred->class == BPredComb) || (pred->class == BPred2Level) || (pred->class == BPredTSBP) || 
	(pred->class == BPredCHBP) || (pred->class == BPredOB) || (pred->class == BPredOHT) || 
	(pred->class == BPredMBP) || (pred->class == BPredTSCL) || (pred->class == BPredLLBP))) {
		//Update FHB and shift up or down depending on flow mode
		if (!flow_mode) {
			pred->fhb[0].addr_out = pred->fhb[0].addr[pred->fhb[0].top];
			pred->fhb[0].fv_out = pred->fhb[0].fv[pred->fhb[0].top];
			pred->fhb[0].rv_out = pred->fhb[0].rv[pred->fhb[0].top];
			pred->fhb[0].o_out = pred->fhb[0].o[pred->fhb[0].top];
			pred->fhb[0].correct_out = pred->fhb[0].correct[pred->fhb[0].top];
			pred->fhb[0].und_out = pred->fhb[0].und[pred->fhb[0].top];
			pred->fhb[0].thu_out = pred->fhb[0].thu[pred->fhb[0].top];
			pred->fhb[0].key_out = pred->fhb[0].key[pred->fhb[0].top];
			
			if (pred->fhb[0].top > 0) {
				pred->fhb[0].top--;
			} else {
				pred->fhb[0].top = pred->fhb[0].size - 1;
			}
			
			if (pred->fhb[0].bot > 0) {
				pred->fhb[0].bot--;
			} else {
				pred->fhb[0].bot = pred->fhb[0].size - 1;
			}
			
			if (pbtb) {
				pred->fhb[0].addr[pred->fhb[0].bot] = pbtb->target;
			} else {
				pred->fhb[0].addr[pred->fhb[0].bot] = NULL;
			}
			
			pred->fhb[0].addr[pred->fhb[0].bot] = rbaddr;
			pred->fhb[0].fv[pred->fhb[0].bot] = 0;
			pred->fhb[0].rv[pred->fhb[0].bot] = 1;
			pred->fhb[0].o[pred->fhb[0].bot] = taken;
			// For TSCL/LLBP correctness is determined by TAGE prediction since a valid outcome could've been used
			if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
				pred->fhb[0].correct[pred->fhb[0].bot] = (dir_update_ptr->fwd_tage_pred==taken);
			} else {
				pred->fhb[0].correct[pred->fhb[0].bot] = correct;
			}
			if (pred->class==BPredLLBP||pred->class==BPredTSCL)
				pred->fhb[0].und[pred->fhb[0].bot] = (dir_update_ptr->fwd_dir.sum < pred->fwd_sc_dirpred[1].twolev->config.two.threshold);
			pred->fhb[0].thu[pred->fhb[0].bot] = 0;
			pred->fhb[0].key[pred->fhb[0].bot] = fwd_key;
		} else {
			pred->fhb[0].addr_out = pred->fhb[0].addr[pred->fhb[0].bot];
			pred->fhb[0].fv_out = pred->fhb[0].fv[pred->fhb[0].bot];
			pred->fhb[0].rv_out = pred->fhb[0].rv[pred->fhb[0].bot];
			pred->fhb[0].o_out = pred->fhb[0].o[pred->fhb[0].bot];
			pred->fhb[0].correct_out = pred->fhb[0].correct[pred->fhb[0].bot];
			pred->fhb[0].und_out = pred->fhb[0].und[pred->fhb[0].bot];
			pred->fhb[0].thu_out = pred->fhb[0].thu[pred->fhb[0].bot];
			pred->fhb[0].key_out = pred->fhb[0].key[pred->fhb[0].bot];
			
			if (pred->fhb[0].top < (pred->fhb[0].size - 1)) {
				pred->fhb[0].top++;
			} else {
				pred->fhb[0].top = 0;
			}
			
			if (pred->fhb[0].bot < (pred->fhb[0].size - 1)) {
				pred->fhb[0].bot++;
			} else {
				pred->fhb[0].bot = 0;
			}
			
			if (pbtb) {
				pred->fhb[0].addr[pred->fhb[0].top] = pbtb->addr;
			} else {
				pred->fhb[0].addr[pred->fhb[0].top] = NULL;
			}
			
			pred->fhb[0].addr[pred->fhb[0].top] = fbaddr;
			pred->fhb[0].fv[pred->fhb[0].top] = 1;
			pred->fhb[0].rv[pred->fhb[0].top] = 0;
			pred->fhb[0].o[pred->fhb[0].top] = taken;
			// For TSCL/LLBP correctness is determined by TAGE prediction since a valid outcome could've been used
			if ((pred->class==BPredTSCL)||(pred->class==BPredLLBP)) {
				if (frmt) {
					pred->fhb[0].correct[pred->fhb[0].bot] = (dir_update_ptr->fwd_tage_pred==taken);
				} else {
					pred->fhb[0].correct[pred->fhb[0].bot] = (dir_update_ptr->rev_tage_pred==taken);
				}
			} else {
				pred->fhb[0].correct[pred->fhb[0].top] = correct;
			}
			if (pred->class==BPredLLBP||pred->class==BPredTSCL){
				if (frmt)
					pred->fhb[0].und[pred->fhb[0].top] = (dir_update_ptr->rev_dir.sum < pred->fwd_sc_dirpred[1].twolev->config.two.threshold);
				else 
					pred->fhb[0].und[pred->fhb[0].top] = (dir_update_ptr->rev_dir.sum < pred->rev_sc_dirpred[1].twolev->config.two.threshold);
			}
			pred->fhb[0].thu[pred->fhb[0].top] = 0;
			// FHB deals with unmasked keys!!!
			pred->fhb[0].key[pred->fhb[0].top] = rev_key;
		}
		
		//if (flow_mode && frmt)	// FRMT was enabled so REV ADDR needs to be exchanged for FWD ADDR
		//	pred->fhb[0].addr_out = addr_from_frmt(pred, pred->fhb[0].addr_out);
    }

	/***********************IF TS, also update correctness buffer*****************************/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (pred->class == BPredTSBP)) {
		//int l1index, shift_reg;
		
		// Get key before updating GHR!!!!
		//int key = key_from_features (pred->fwd_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
    
		/* update L1 table, same as 2lev/comb predictors above this */
		//l1index = (fbaddr >> MD_BR_SHIFT) & (pred->fwd_dirpred.twolev->config.two.l1size - 1);
		//shift_reg = (pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] << 1) | (!!taken);
		//pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->fwd_dirpred.twolev->config.two.shift_width) - 1);

		int base_outcome;
		int ts_outcome;
		//unsigned int key;  /*added typedef in tsbp.h file*/
      	
		fwd_key = fwd_key & (pred->fwd_dirpred.tsbp->ts.head_table_size - 1); // mask key based on predictor table size
		
		/*Set the base outcome; Also if in replay mode and ts_outcome is incorrect, turn off replay mode*/
		if (pred->fwd_dirpred.tsbp->ts.replay) {
			pred->replays++;
			ts_outcome = pred_taken;

			if (pred->fwd_dirpred.tsbp->ts.enabled && (pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.head] == 0)) {
				base_outcome = !pred_taken;
			} else {
				base_outcome = pred_taken;
			}

			if (!!ts_outcome != !!taken) {
				pred->fwd_dirpred.tsbp->ts.replay = FALSE;
			}
		} else {
			base_outcome = pred_taken;
		}
		
		/* incr tail but prevent from going out of bounds*/
		if (pred->fwd_dirpred.tsbp->ts.tail >= pred->fwd_dirpred.tsbp->ts.correctness_width) {
			pred->fwd_dirpred.tsbp->ts.tail = 0;
		} else {
			pred->fwd_dirpred.tsbp->ts.tail++;
		}
		
		/*determine if actual outcome (taken) of predicted direction is correct and update correctness buffer*/
		pred->fwd_dirpred.tsbp->ts.correctness_buffer[pred->fwd_dirpred.tsbp->ts.tail] = (!!base_outcome == !!taken);  /*1 = base predictor correct. 0 = prediction incorrect*/
      
		/*if incorrect base prediction, update head table*/
		if (!!base_outcome != !!taken) {
			if(!pred->fwd_dirpred.tsbp->ts.replay && pred->fwd_dirpred.tsbp->ts.head_table[fwd_key] != NULL) { /*if not in replay mode, update head and set replay flag*/
				pred->fwd_dirpred.tsbp->ts.head = pred->fwd_dirpred.tsbp->ts.head_table[fwd_key];
				pred->fwd_dirpred.tsbp->ts.replay = TRUE;
			}
		}

		pred->fwd_dirpred.tsbp->ts.head_table[fwd_key] = pred->fwd_dirpred.tsbp->ts.tail;   //else update head table
	}
  
	/***********************IF CHBP, also update correctness history table*****************************/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (pred->class == BPredCHBP)) {
		//int l1index, shift_reg;
		
		// Get key before updating GHR!!!!
		//int key = key_from_features (pred->fwd_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
		
		/* update L1 table, same as 2lev/comb predictors above this */
		//l1index = (fbaddr >> MD_BR_SHIFT) & (pred->fwd_dirpred.twolev->config.two.l1size - 1);
		//shift_reg = (pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] << 1) | (!!taken);
		//pred->fwd_dirpred.twolev->config.two.shiftregs[l1index] = shift_reg & ((1 << pred->fwd_dirpred.twolev->config.two.shift_width) - 1);

		int base_outcome;
		int chbp_outcome;
		
		fwd_key = fwd_key & (pred->fwd_dirpred.chbp->chbp.cht_size - 1); // mask key based on predictor table size
		
		bool_t just_disabled_replay = FALSE;

		/*Set the base outcome; Also if in replay mode and chbp_outcome is incorrect, turn off replay mode*/
		if (pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] && (pred->fwd_dirpred.chbp->chbp.cht_spc[fwd_key] == fbaddr)) {
			pred->replays++;
			chbp_outcome = pred_taken;

			if (pred->fwd_dirpred.chbp->chbp.enabled && !pred->fwd_dirpred.chbp->chbp.cht_correct[fwd_key]) {
				base_outcome = !pred_taken;
			} else {
				base_outcome = pred_taken;
			}

			if (!!chbp_outcome != !!taken) {
				pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] = FALSE;
				just_disabled_replay = TRUE;
			}
		} else {
			base_outcome = pred_taken;
		}
		
		/*determine if actual outcome (taken) of predicted direction is correct and update correctness bit*/
		pred->fwd_dirpred.chbp->chbp.cht_correct[fwd_key] = (!!base_outcome == !!taken);  /*1 = base predictor correct. 0 = prediction incorrect*/
		  
		/*if not in replay mode and base outcome incorrect, set replay flag as long as it wasn't just disabled*/
		if ((!!base_outcome != !!taken) && !pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] && !just_disabled_replay) {
			pred->fwd_dirpred.chbp->chbp.cht_replay[fwd_key] = TRUE;
		}
		
		/* Update other correctness history table bits */
		pred->fwd_dirpred.chbp->chbp.cht_dpc[fwd_key] = fbtarget;
		pred->fwd_dirpred.chbp->chbp.cht_valid[fwd_key] = TRUE;
		pred->fwd_dirpred.chbp->chbp.cht_spc[fwd_key] = fbaddr;
		
		/* make sure no other matching spc keys are valid*/
		//int valid_key;
		
		//for (valid_key = 0; valid_key < pred->fwd_dirpred.chbp->chbp.cht_size; valid_key++) {
		//	if ((pred->fwd_dirpred.chbp->chbp.cht_spc[valid_key] == fbaddr) && pred->fwd_dirpred.chbp->chbp.cht_valid[valid_key] && (valid_key != key)) {
		//		pred->fwd_dirpred.chbp->chbp.cht_valid[valid_key] = FALSE;
		//	}
		//}
		
	}
	
	/* Update OB if OB or MBP or TSCL or LLBP */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && ((pred->class == BPredOB) || (pred->class == BPredMBP) || (pred->class == BPredTSCL) || (pred->class == BPredLLBP))) {
		//Update OB and shift up (right) or down (left) depending on flow mode
		if (!flow_mode) {
			if (pred->ob.end > 0) {
				pred->ob.end--;
			} else {
				pred->ob.end = pred->ob.width - 1;
			}
			
			if (pred->ob.beg > 0) {
				pred->ob.beg--;
			} else {
				pred->ob.beg = pred->ob.width - 1;
			}
			
			pred->ob.fv[pred->ob.beg] = 0;
			pred->ob.rv[pred->ob.beg] = 1;
			pred->ob.oc[pred->ob.beg] = taken;
		} else {
			if (pred->ob.end < (pred->ob.width - 1)) {
				pred->ob.end++;
			} else {
				pred->ob.end = 0;
			}
			
			if (pred->ob.beg < (pred->ob.width - 1)) {
				pred->ob.beg++;
			} else {
				pred->ob.beg = 0;
			}
			
			pred->ob.fv[pred->ob.end] = 1;
			pred->ob.rv[pred->ob.end] = 0;
			pred->ob.oc[pred->ob.end] = taken;
		}
	}
	
	/* Update FRMTs keys and contexts if enabled */
	/* This needs to be updated before the OHT(s) to ensure proper mapping for the FOHT*/
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (!flow_mode) && (frmt) && 
	((pred->class == BPred2bit) || (pred->class == BPred2Level) || (pred->class == BPredOB) || 
	(pred->class == BPredOHT) || (pred->class == BPredMBP) || (pred->class == BPredTSCL) || 
	(pred->class == BPredLLBP))) {
		switch (pred->class) {
			case BPredLLBP:
				// Map FRMT context entries; skipping until after FCB update
			case BPredComb:
			case BPred2Level:
			case BPredTSBP:
			case BPredCHBP:
			case BPredOB:
			case BPredOHT:
			case BPredMBP:
			case BPredTSCL:
				// Only update if a valid future history is available
				if (pred->fhb[0].rv_out) {
					/* Get future key from FHB */
					int future_key;
					
					if ((pred->class == BPredTSCL) || (pred->class == BPredLLBP))
						future_key = key_from_features (pred->fwd_tage_dirpred[pred->tage_depth-1].twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
					else
						future_key = key_from_features (pred->fwd_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
					
					//future_key = future_key & (pred->fwd_dirpred.twolev->config.two.l2size - 1); // mask key based on L2 table size
					//pred->fhb[0].key_out = pred->fhb[0].key_out & (pred->fwd_dirpred.twolev->config.two.l2size - 1); // mask key based on L2 table size
					
					// Need to map FWD to REV keys based on mode
					map_unmasked_key_to_frmt(pred, pred->fhb[0].key_out, future_key);
				}
			
			case BPred2bit:
			case BPredTaken:
			case BPredNotTaken:
				break;
			
			/* no other state */
			default:
				panic("bogus predictor class");
		}
	}
	
	/* Update OHT if OHT or MBP or TSCL or LLBP */
	if (((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND)) && ((pred->class == BPredOHT) || (pred->class == BPredMBP) || 
	(pred->class == BPredTSCL) || (pred->class == BPredLLBP))) {
		// Get future and current addresses
		/* Get FWD keys */
		int fwd_future_key;
		
		if ((pred->class == BPredOHT) || (pred->class == BPredMBP)) {
			fwd_future_key = key_from_features (pred->fwd_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
		} else {
			fwd_future_key = key_from_features (pred->fwd_tage_dirpred[pred->tage_depth-1].twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
		}
		
		fwd_future_key = fwd_future_key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
		fwd_key = fwd_key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
		
		/* Get REV keys; dpending on FRMT, getting from FWD table instead */
		int rev_future_key;
		if (!frmt) {
			if ((pred->class == BPredOHT) || (pred->class == BPredMBP)) {
				rev_future_key = key_from_features (pred->rev_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
			} else {
				rev_future_key = key_from_features (pred->rev_tage_dirpred[pred->tage_depth-1].twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
			}
			
			rev_future_key = rev_future_key & (pred->rev_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
			rev_key = rev_key & (pred->rev_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
		} else {
			if ((pred->class == BPredOHT) || (pred->class == BPredMBP)) {
				rev_future_key = key_from_features (pred->fwd_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
			} else {
				rev_future_key = key_from_features (pred->fwd_tage_dirpred[pred->tage_depth-1].twolev, pred->fhb[0].addr_out); // Get unmasked key from GHR and PC
			}
			
			rev_future_key = rev_future_key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
			rev_key = rev_key & (pred->fwd_dirpred.oht->oht.size - 1); // mask key based on outcome history table size
		}

		
		// Update only one OHT depending on the flow and 
		if (flow_mode) {
			// updating to future fwd if valid fhb value
			if(pred->fhb[0].fv_out) {
				pred->fwd_dirpred.oht->oht.oc[fwd_future_key] = pred->fhb[0].o_out;
				pred->fwd_dirpred.oht->oht.fwd_valid[fwd_future_key] = pred->fhb[0].fv_out;
				pred->fwd_dirpred.oht->oht.rev_valid[fwd_future_key] = pred->fhb[0].rv_out;
			}
			
			// Unset REV valid bit
			if (!frmt) {	// If not FRMT, using both sets of tables
				pred->rev_dirpred.oht->oht.rev_valid[rev_key] = 0;          //Even if not used for predict, need to invalidate REV mode OHT Val at this addr.
			} else {
				pred->fwd_dirpred.oht->oht.rev_valid[rev_key] = 0;          //Even if not used for predict, need to invalidate REV mode OHT Val at this addr.
			}
		} else {
			// Unset FWD valid bit
			pred->fwd_dirpred.oht->oht.fwd_valid[fwd_key] = 0;		//Even if not used for predict, need to invalidate FWD mode OHT Val at this addr.
			
			// updating to future rev if valid fhb value
			if(pred->fhb[0].rv_out) {
				if (!frmt) {	// If not FRMT, using both sets of tables
					pred->rev_dirpred.oht->oht.oc[rev_future_key] = pred->fhb[0].o_out; //Update with pred->fhb[0].o_out since it is the past outcome that was waiting in the FHB
					pred->rev_dirpred.oht->oht.fwd_valid[rev_future_key] = pred->fhb[0].fv_out;
					pred->rev_dirpred.oht->oht.rev_valid[rev_future_key] = pred->fhb[0].rv_out;
				} else {
					pred->fwd_dirpred.oht->oht.oc[rev_future_key] = pred->fhb[0].o_out; //Update with pred->fhb[0].o_out since it is the past outcome that was waiting in the FHB
					pred->fwd_dirpred.oht->oht.fwd_valid[rev_future_key] = pred->fhb[0].fv_out;
					pred->fwd_dirpred.oht->oht.rev_valid[rev_future_key] = pred->fhb[0].rv_out;
				}
			}
		}
	}
	
	int replacement_max = 1;
	int replacement_size = 0;
	int* tage_replacement_tags = (int*)malloc(replacement_max * sizeof(int));
	int* tage_replacement_counters = (int*)malloc(replacement_max * sizeof(int));
	int* tage_replacement_lengths = (int*)malloc(replacement_max * sizeof(int));
	
	/* Update TAGE, SC, and Loop Preds if TSCL/LLBP; If LLBP, save TAGE victims as replacement candidates for LLBT */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (pred->class == BPredLLBP||pred->class == BPredTSCL)) {
		// tage correct value is determined based on tage_pred which is saved pre-inversion; all we care about is mode & frmt
		int tage_correct;
		if (!flow_mode || frmt)	{
			tage_correct = (dir_update_ptr->fwd_tage_pred==taken);	 
		} else {
			tage_correct = (dir_update_ptr->rev_tage_pred==taken);
		}
		
		// depending on flow mode & frmt, different entries may be updated with different keys
		if (flow_mode) { // The full REV mode need to update both sides
			// Checking against future fwd if valid fhb value and not FRMTs
			if(pred->fhb[0].fv_out && !frmt) {		// don't update if FRMT set as this is now shared and update already occured
				bpred_tage_update(
					pred,
					&replacement_size,
					&replacement_max,
					&tage_replacement_counters,			/* replacement pointers */
					&tage_replacement_tags,				/* replacement pointers */
					&tage_replacement_lengths,			/* replacement pointers */
					pred->fhb[0].addr_out,				/* Branch address */
					pred->fhb[0].correct_out,			/* predictor correctness bit */
					pred->fhb[0].o_out,					/* Branch Outcome */
					0,									/* flow_mode flag; for updating means update function mode not actual flow_mode */
					1,									/* future mode */
					frmt								/* FRMTs enabled flag*/
				);
			}
			
			// Checking against current rev
			bpred_tage_update(
				pred,
				&replacement_size,
				&replacement_max,
				&tage_replacement_counters,			/* replacement pointers */
				&tage_replacement_tags,				/* replacement pointers */
				&tage_replacement_lengths,			/* replacement pointers */
				rbaddr,								/* Branch address */
				tage_correct,							/* predictor correctness bit */
				taken,								/* Branch Outcome */
				1,									/* flow_mode flag; for updating means update function mode not actual flow_mode */
				0,									/* future mode */
				frmt								/* FRMTs enabled flag*/
			);
		} else {
			// Checking against current fwd
			bpred_tage_update(
				pred,
				&replacement_size,
				&replacement_max,
				&tage_replacement_counters,			/* replacement pointers */
				&tage_replacement_tags,				/* replacement pointers */
				&tage_replacement_lengths,			/* replacement pointers */
				fbaddr,								/* Branch address */
				tage_correct,							/* predictor correctness bit */
				taken,								/* Branch Outcome */
				0,									/* flow_mode flag; for updating means update function mode not actual flow_mode */
				0,									/* future mode */
				frmt								/* FRMTs enabled flag*/
			);
			
			// Checking against future rev if valid fhb value and not FRMTs
			if(pred->fhb[0].rv_out && !frmt) {		// don't update if FRMT set as this is now shared and update already occured
				bpred_tage_update(
					pred,
					&replacement_size,
					&replacement_max,
					&tage_replacement_counters,			/* replacement pointers */
					&tage_replacement_tags,				/* replacement pointers */
					&tage_replacement_lengths,			/* replacement pointers */
					pred->fhb[0].addr_out,				/* Branch address */
					pred->fhb[0].correct_out,			/* predictor correctness bit */
					pred->fhb[0].o_out,					/* Branch Outcome */
					1,									/* flow_mode flag; for updating means update function mode not actual flow_mode */
					1,									/* future mode */
					frmt								/* FRMTs enabled flag*/
				);
			}
		}
		
		for (int i = 1; i < pred->tage_depth; i++) {
			// Update history
			shift_history_left(pred->fwd_tage_dirpred[i].twolev, fbaddr, taken);
			
			if (!frmt) {
				// Update history
				shift_history_left(pred->rev_tage_dirpred[i].twolev, rbaddr, taken);
			}
		}
		
		for (int i = 1; i < pred->tage_depth; i++) {
			int unmasked_fwd_sc_key, 
				unmasked_fwd_future_sc_key, 
				fwd_sc_key, fwd_future_sc_key, 
				unmasked_rev_sc_key, 
				unmasked_rev_future_sc_key, 
				rev_sc_key, 
				rev_future_sc_key;
			
			// Get sc tags and compare to current tag
			unmasked_fwd_sc_key = key_from_features (pred->fwd_sc_dirpred[i].twolev, fbaddr); // Get unmasked key from GHR and PC
			unmasked_fwd_future_sc_key = key_from_features (pred->fwd_sc_dirpred[i].twolev, pred->fhb[0].addr_out); // Get unmasked future key from GHR and PC
			
			fwd_sc_key = unmasked_fwd_sc_key & (pred->fwd_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
			fwd_future_sc_key = unmasked_fwd_future_sc_key & (pred->fwd_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
			
			if (frmt) {
				unmasked_rev_sc_key = key_from_features (pred->fwd_sc_dirpred[i].twolev, fbaddr); // Get unmasked key from GHR and PC
				unmasked_rev_future_sc_key = key_from_features (pred->fwd_sc_dirpred[i].twolev, pred->fhb[0].addr_out); // Get unmasked future key from GHR and PC
				
				//unmasked_rev_sc_key = unmasked_key_from_frmt(pred, unmasked_rev_sc_key); // FRMT was enabled so REV key needs to be exchanged for FWD key
				//unmasked_rev_future_sc_key = unmasked_key_from_frmt(pred, unmasked_rev_future_sc_key); // FRMT was enabled so REV key needs to be exchanged for FWD key
				
				rev_sc_key = unmasked_fwd_sc_key & (pred->fwd_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
				rev_future_sc_key = unmasked_fwd_future_sc_key & (pred->fwd_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
			} else {
				unmasked_rev_sc_key = key_from_features (pred->rev_sc_dirpred[i].twolev, fbaddr); // Get unmasked key from GHR and PC
				unmasked_rev_future_sc_key = key_from_features (pred->rev_sc_dirpred[i].twolev, pred->fhb[0].addr_out); // Get unmasked future key from GHR and PC
				
				rev_sc_key = unmasked_fwd_sc_key & (pred->rev_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
				rev_future_sc_key = unmasked_fwd_future_sc_key & (pred->rev_sc_dirpred[i].twolev->config.two.l2size - 1); // mask key based on predictor table size
			}
			
			// If matched unmasked SC key, increment useful counter
			// depending on flow mode & frmt, different entries may be updated with different keys
			if (flow_mode) { // The full REV mode need to update both sides
				// Checking against future fwd if valid fhb value
				if(pred->fhb[0].fv_out && !frmt) {	// don't update if FRMT set as this is now shared and update already occured
					bpred_sc_update(
						pred->fwd_sc_dirpred[i].twolev,
						fwd_future_sc_key,
						unmasked_fwd_future_sc_key,
						pred->fhb[0].correct_out,
						pred->fhb[0].und_out,
						pred->fhb[0].thu_out
					);
				}
				
				// Checking against current rev
				if (!frmt) {	// If not FRMT, using both sets of tables
					bpred_sc_update(
						pred->rev_sc_dirpred[i].twolev,
						rev_sc_key,
						unmasked_rev_sc_key,
						tage_correct,
						(dir_update_ptr->rev_dir.sum < pred->rev_sc_dirpred[i].twolev->config.two.threshold),
						0
					);
					
					// Update history
					shift_history_right(pred->rev_sc_dirpred[i].twolev, rbaddr, taken);
				} else {
					bpred_sc_update(
						pred->fwd_sc_dirpred[i].twolev,
						rev_sc_key,
						unmasked_rev_sc_key,
						tage_correct,
						(dir_update_ptr->rev_dir.sum < pred->fwd_sc_dirpred[i].twolev->config.two.threshold),
						0
					);
				}

				// Update history
				shift_history_right(pred->fwd_sc_dirpred[i].twolev, fbaddr, taken);
			} else {
				// Checking against current fwd
				bpred_sc_update(
					pred->fwd_sc_dirpred[i].twolev,
					fwd_sc_key,
					unmasked_fwd_sc_key,
					tage_correct,
					(dir_update_ptr->fwd_dir.sum < pred->fwd_sc_dirpred[i].twolev->config.two.threshold),
					0
				);
				
				// Checking against future rev if valid fhb value
				if(pred->fhb[0].rv_out && !frmt) {		// don't update if FRMT set as this is now shared and update already occured
					// If not FRMT, using both sets of tables
					bpred_sc_update(
						pred->rev_sc_dirpred[i].twolev,
						rev_future_sc_key,
						unmasked_rev_future_sc_key,
						pred->fhb[0].correct_out,
						pred->fhb[0].und_out,
						pred->fhb[0].thu_out
					);
					
					// Update history
					shift_history_left(pred->rev_sc_dirpred[i].twolev, rbaddr, taken);
				}
				
				// Update history
				shift_history_left(pred->fwd_sc_dirpred[i].twolev, fbaddr, taken);
			}
		}
			
		// Upate loop table and compare to current 
		int unmasked_fwd_loop_key, 
			unmasked_fwd_future_loop_key, 
			fwd_loop_key, fwd_future_loop_key, 
			unmasked_rev_loop_key, 
			unmasked_rev_future_loop_key, 
			rev_loop_key, 
			rev_future_loop_key;
		
		// Get loop tags and compare to current tag
		unmasked_fwd_loop_key = key_from_features (pred->fwd_loop_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
		unmasked_fwd_future_loop_key = key_from_features (pred->fwd_loop_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked future key from GHR and PC
		
		fwd_loop_key = unmasked_fwd_loop_key & (pred->fwd_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
		fwd_future_loop_key = unmasked_fwd_future_loop_key & (pred->fwd_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
		
		if (frmt) {
			unmasked_rev_loop_key = key_from_features (pred->fwd_loop_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
			unmasked_rev_future_loop_key = key_from_features (pred->fwd_loop_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked future key from GHR and PC
			
			//unmasked_rev_loop_key = unmasked_key_from_frmt(pred, unmasked_rev_loop_key); // FRMT was enabled so REV key needs to be exchanged for FWD key
			//unmasked_rev_future_loop_key = unmasked_key_from_frmt(pred, unmasked_rev_future_loop_key); // FRMT was enabled so REV key needs to be exchanged for FWD key
			
			rev_loop_key = unmasked_fwd_loop_key & (pred->fwd_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
			rev_future_loop_key = unmasked_fwd_future_loop_key & (pred->fwd_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
		} else {
			unmasked_rev_loop_key = key_from_features (pred->rev_loop_dirpred.twolev, fbaddr); // Get unmasked key from GHR and PC
			unmasked_rev_future_loop_key = key_from_features (pred->rev_loop_dirpred.twolev, pred->fhb[0].addr_out); // Get unmasked future key from GHR and PC
			
			rev_loop_key = unmasked_fwd_loop_key & (pred->rev_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
			rev_future_loop_key = unmasked_fwd_future_loop_key & (pred->rev_loop_dirpred.twolev->config.two.l1size - 1); // mask key based on predictor table size
		}
			
		if (flow_mode) { // The full REV mode need to update both sides
			// Checking against future fwd if valid fhb value
			if(pred->fhb[0].fv_out && !frmt) {	// don't update if FRMT set as this is now shared and update already occured
				bpred_loop_update(
					pred->fwd_loop_dirpred.twolev,
					fwd_future_loop_key,
					unmasked_fwd_future_loop_key,
					pred->fhb[0].o_out
				);
			}
			
			// Checking against current rev
			if (!frmt) {	// If not FRMT, using both sets of tables
				bpred_loop_update(
					pred->rev_loop_dirpred.twolev,
					rev_loop_key,
					unmasked_rev_loop_key,
					taken
				);
				
				// Update history
				shift_history_right(pred->rev_loop_dirpred.twolev, rbaddr, taken);
			} else {
				bpred_loop_update(
					pred->fwd_loop_dirpred.twolev,
					rev_loop_key,
					unmasked_rev_loop_key,
					taken
				);
			}
			
			// Update history
			shift_history_right(pred->fwd_loop_dirpred.twolev, fbaddr, taken);
		} else {
			// Checking against current fwd
			bpred_loop_update(
				pred->fwd_loop_dirpred.twolev,
				fwd_loop_key,
				unmasked_fwd_loop_key,
				taken
			);
			
			// Checking against future rev if valid fhb value
			if(pred->fhb[0].rv_out && !frmt) {		// don't update if FRMT set as this is now shared and update already occured
				// If not FRMT, using both sets of tables
				bpred_loop_update(
					pred->rev_loop_dirpred.twolev,
					rev_future_loop_key,
					unmasked_rev_future_loop_key,
					pred->fhb[0].o_out
				);
				
				// Update history
				shift_history_left(pred->rev_loop_dirpred.twolev, rbaddr, taken);
			}
			
			// Update history
			shift_history_left(pred->fwd_loop_dirpred.twolev, fbaddr, taken);
		}
	}
	
	/* Update LLBT/PB if LLBP; opposite flow updates depend on FCB to have a */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (pred->class == BPredLLBP)) {
		if (flow_mode) { // The full REV mode need to update both sides
			// "Post-Fetch" fccid to FWD PB
			bpred_llbp_fetch(
				pred,	/* branch predictor instance */
				pred->fwd_pb_dirpred,	/* branch dir predictor inst */
				pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
				pred->rcr.fccid					/* CCID */
			);				
			
			// Checking against future fwd if valid fhb value
			if(pred->fhb[0].fv_out) {
				bpred_llbp_update(
					pred,	/* branch predictor instance */
					pred->fwd_pb_dirpred,	/* branch dir predictor inst */
					pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.ccid,				/* CCID */
					pred->fhb[0].addr_out,			/* addr to get keys */
					&replacement_max,
					&tage_replacement_counters,				/* replacement pointers */
					&tage_replacement_tags,				/* replacement pointers */
					&tage_replacement_lengths,				/* replacement pointers */
					pred->fhb[0].o_out,
					flow_mode,
					1
				);
			}
			
			// Checking against current rev
			if (!frmt) {	// If not FRMT, using both sets of tables
				bpred_llbp_update(
					pred,	/* branch predictor instance */
					pred->rev_pb_dirpred,	/* branch dir predictor inst */
					pred->rev_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.ccid,				/* CCID */
					rbaddr,			/* addr to get keys */
					&replacement_max,
					&tage_replacement_counters,				/* replacement pointers */
					&tage_replacement_tags,				/* replacement pointers */
					&tage_replacement_lengths,				/* replacement pointers */
					taken,
					flow_mode,
					1
				);
				
				// pre-fetch rccid to REV PB
				bpred_llbp_fetch(
					pred,	/* branch predictor instance */
					pred->rev_pb_dirpred,	/* branch dir predictor inst */
					pred->rev_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.rccid					/* CCID */
				);
			} else {
				bpred_llbp_update(
					pred,	/* branch predictor instance */
					pred->fwd_pb_dirpred,	/* branch dir predictor inst */
					pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.ccid,				/* CCID */
					&replacement_max,
					&tage_replacement_counters,				/* replacement pointers */
					&tage_replacement_tags,				/* replacement pointers */
					&tage_replacement_lengths,				/* replacement pointers */
					rbaddr,			/* addr to get keys */
					taken,
					flow_mode,
					0
				);
				
				// pre-fetch rccid to FWD PB
				bpred_llbp_fetch(
					pred,	/* branch predictor instance */
					pred->fwd_pb_dirpred,	/* branch dir predictor inst */
					pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.rccid					/* CCID */
				);
			}	
		} else {
			// "Post-Fetch" rccid to REV PB
			bpred_llbp_fetch(
				pred,	/* branch predictor instance */
				pred->rev_pb_dirpred,	/* branch dir predictor inst */
				pred->rev_llbt_dirpred,	/* branch dir predictor inst */
				pred->rcr.rccid					/* CCID */
			);				
			
			// Checking against current fwd
			bpred_llbp_update(
				pred,	/* branch predictor instance */
				pred->fwd_pb_dirpred,	/* branch dir predictor inst */
				pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
				pred->rcr.ccid,				/* CCID */
				fbaddr,			/* addr to get keys */
				&replacement_max,
				&tage_replacement_counters,				/* replacement pointers */
				&tage_replacement_tags,				/* replacement pointers */
				&tage_replacement_lengths,				/* replacement pointers */
				taken,
				flow_mode,
				1
			);
			
			// Checking against future rev
			if (!frmt) {	// If not FRMT, using both sets of tables
				// Only if valid fhb value
				if(pred->fhb[0].rv_out) {
					bpred_llbp_update(
						pred,	/* branch predictor instance */
						pred->rev_pb_dirpred,	/* branch dir predictor inst */
						pred->rev_llbt_dirpred,	/* branch dir predictor inst */
						pred->rcr.ccid,				/* CCID */
						pred->fhb[0].addr_out,			/* addr to get keys */
						&replacement_max,
						&tage_replacement_counters,				/* replacement pointers */
						&tage_replacement_tags,				/* replacement pointers */
						&tage_replacement_lengths,				/* replacement pointers */
						pred->fhb[0].o_out,
						flow_mode,
						1
					);
				}
				
				// pre-fetch fccid to REV PB
				bpred_llbp_fetch(
					pred,	/* branch predictor instance */
					pred->rev_pb_dirpred,	/* branch dir predictor inst */
					pred->rev_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.fccid					/* CCID */
				);
			} else {
				// Only if valid fhb value
				if(pred->fhb[0].rv_out) {
					bpred_llbp_update(
						pred,	/* branch predictor instance */
						pred->fwd_pb_dirpred,	/* branch dir predictor inst */
						pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
						pred->rcr.ccid,				/* CCID */
						pred->fhb[0].addr_out,			/* addr to get keys */
						&replacement_max,
						&tage_replacement_counters,				/* replacement pointers */
						&tage_replacement_tags,				/* replacement pointers */
						&tage_replacement_lengths,				/* replacement pointers */
						pred->fhb[0].o_out,
						flow_mode,
						0
					);
				}
				
				// pre-fetch fccid to FWD PB
				bpred_llbp_fetch(
					pred,	/* branch predictor instance */
					pred->fwd_pb_dirpred,	/* branch dir predictor inst */
					pred->fwd_llbt_dirpred,	/* branch dir predictor inst */
					pred->rcr.fccid					/* CCID */
				);
			}	
		}
	}
	
	// Update future l2table info if Combo, 2Level, OHT, OB, MBP, CHBP or TSBP ??and not FRMTs??
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) && (!frmt) && ((pred->class == BPredComb) || 
	(pred->class == BPred2Level) || (pred->class == BPredTSBP) || (pred->class == BPredCHBP) || 
	(pred->class == BPredOB) || (pred->class == BPredOHT) || (pred->class == BPredMBP))) {
		int fwd_future_key, rev_future_key; 
		
		fwd_future_key = key_from_features(pred->fwd_dirpred.twolev, pred->fhb[0].addr_out);
		fwd_future_key = fwd_future_key & (pred->fwd_dirpred.twolev->config.two.l2size-1);
		
		if (frmt) {
			rev_future_key = key_from_features(pred->fwd_dirpred.twolev, pred->fhb[0].addr_out);
			rev_future_key = rev_future_key & (pred->fwd_dirpred.twolev->config.two.l2size-1);
		} else {
			rev_future_key = key_from_features(pred->rev_dirpred.twolev, pred->fhb[0].addr_out);
			rev_future_key = rev_future_key & (pred->rev_dirpred.twolev->config.two.l2size-1);
		}
		
		if (flow_mode) {
			// Need to update future FWD table if valid FHB data from REV 
			if (pred->fhb[0].fv_out) {
				if (pred->fhb[0].o_out) {
					if (pred->fwd_dirpred.twolev->config.two.l2table[fwd_future_key] < 3)
						pred->fwd_dirpred.twolev->config.two.l2table[fwd_future_key]++;
				} else { /* not taken */
					if (pred->fwd_dirpred.twolev->config.two.l2table[fwd_future_key])
						pred->fwd_dirpred.twolev->config.two.l2table[fwd_future_key]--;
				}
			}
		} else {
			// Need to update future REV table if valid FHB data from FWD 
			if (pred->fhb[0].rv_out) {
				if (frmt) {
					if (0) {
					if (pred->fhb[0].o_out) {
						if (pred->fwd_dirpred.twolev->config.two.l2table[rev_future_key] < 3)
							pred->fwd_dirpred.twolev->config.two.l2table[rev_future_key]++;
					} else { /* not taken */
						if (pred->fwd_dirpred.twolev->config.two.l2table[rev_future_key])
							pred->fwd_dirpred.twolev->config.two.l2table[rev_future_key]--;
					}
					}
				} else {
					if (pred->fhb[0].o_out) {
						if (pred->rev_dirpred.twolev->config.two.l2table[rev_future_key] < 3)
							pred->rev_dirpred.twolev->config.two.l2table[rev_future_key]++;
					} else { /* not taken */
						if (pred->rev_dirpred.twolev->config.two.l2table[rev_future_key])
							pred->rev_dirpred.twolev->config.two.l2table[rev_future_key]--;
					}
				}
			}
		}
	}
	
	/* update dir_update_ptr state(s) (but not for jumps) */
	if ((pred->class == BPredTSCL) || (pred->class == BPredLLBP)) {
		// Need to update the base bimod when using TSCL or LLBP
		if (frmt) {
			if (flow_mode) {
				// Only update REV table in REV mode with FRMT set
				if (taken) {
					if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, rbaddr)] < 3)
						pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, rbaddr)]++;
				} else { /* not taken */
					if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, rbaddr)] > 0)
						pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, rbaddr)]--;
				}
			} else {
				// Only update FWD table in FWD mode with FRMT set
				if (taken) {
					if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)] < 3)
						pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)]++;
				} else { /* not taken */
					if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)] > 0)
						pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)]--;
				}
			}
		} else {
			// Update both REV and FWD but depending on mode, addr is determined by FHB
			if (flow_mode) {	// FWD is future
				if (taken) {
					if (pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, rbaddr)] < 3)
						pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, rbaddr)]++;
				} else { /* not taken */
					if (pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, rbaddr)] > 0)
						pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, rbaddr)]--;
				}
				
				if (pred->fhb[0].fv_out) {
					if (pred->fhb[0].o_out) {
						if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, pred->fhb[0].addr_out)] < 3)
							pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, pred->fhb[0].addr_out)]++;
					} else { /* not taken */
						if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, pred->fhb[0].addr_out)] > 0)
							pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, pred->fhb[0].addr_out)]--;
					}
				}
			} else {	// REV is future
				if (taken) {
					if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)] < 3)
						pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)]++;
				} else { /* not taken */
					if (pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)] > 0)
						pred->fwd_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->fwd_dirpred.bimod, fbaddr)]--;
				}
				
				if (pred->fhb[0].rv_out) {
					if (pred->fhb[0].o_out) {
						if (pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, pred->fhb[0].addr_out)] < 3)
							pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, pred->fhb[0].addr_out)]++;
					} else { /* not taken */
						if (pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, pred->fhb[0].addr_out)] > 0)
							pred->rev_dirpred.bimod->config.bimod.table[BIMOD_HASH(pred->rev_dirpred.bimod, pred->fhb[0].addr_out)]--;
					}
				}
			}
		}
	} else {	// Do things the old way		
		if ((!flow_mode) || (pred->class == BPred2bit && !frmt)) { // don't update in REV mode if FRMTs set and it's Bimod
			if (dir_update_ptr->fwd_pdir1) {  
				if (taken) {
					if (*dir_update_ptr->fwd_pdir1 < 3)
						++*dir_update_ptr->fwd_pdir1;
				} else { /* not taken */
					if (*dir_update_ptr->fwd_pdir1 > 0)
						--*dir_update_ptr->fwd_pdir1;
				}
			}
			
			/* combining predictor also updates second predictor and meta predictor */
			/* second direction predictor */
			if (dir_update_ptr->fwd_pdir2) {
				if (taken) {
					if (*dir_update_ptr->fwd_pdir2 < 3)
						++*dir_update_ptr->fwd_pdir2;
				} else { /* not taken */
					if (*dir_update_ptr->fwd_pdir2 > 0)
						--*dir_update_ptr->fwd_pdir2;
				}
			}

			/* meta predictor */
			if (dir_update_ptr->fwd_pmeta) {
				if (dir_update_ptr->fwd_dir.bimod != dir_update_ptr->fwd_dir.twolev) {
					/* we only update meta predictor if directions were different */
					if (dir_update_ptr->fwd_dir.twolev == (unsigned int)taken) {
						/* 2-level predictor was correct */
						if (*dir_update_ptr->fwd_pmeta < 3)
							++*dir_update_ptr->fwd_pmeta;
					} else {
						/* bimodal predictor was correct */
						if (*dir_update_ptr->fwd_pmeta > 0)
							--*dir_update_ptr->fwd_pmeta;
					}
				}
			}
		}
		
		/* Acounting for reverse dir pointer now */
		if ((flow_mode) || (pred->class == BPred2bit && !frmt)) {	// don't update in FWD mode if FRMTs set and it's Bimod
			if (dir_update_ptr->rev_pdir1) {  
				if (taken) {
					if (*dir_update_ptr->rev_pdir1 < 3)
						++*dir_update_ptr->rev_pdir1;
				} else { /* not taken */
					if (*dir_update_ptr->rev_pdir1 > 0)
						--*dir_update_ptr->rev_pdir1;
				}
			}
			
			if (dir_update_ptr->rev_pdir2) {  
				if (taken) {
					if (*dir_update_ptr->rev_pdir2 < 3)
						++*dir_update_ptr->rev_pdir2;
				} else { /* not taken */
					if (*dir_update_ptr->rev_pdir2 > 0)
						--*dir_update_ptr->rev_pdir2;
				}
			}
			
			if (dir_update_ptr->rev_pmeta) {
				if (dir_update_ptr->rev_dir.bimod != dir_update_ptr->rev_dir.twolev) {
					/* we only update meta predictor if directions were different */
					if (dir_update_ptr->rev_dir.twolev == (unsigned int)taken) {
						/* 2-level predictor was correct */
						if (*dir_update_ptr->rev_pmeta < 3)
							++*dir_update_ptr->rev_pmeta;
					} else {
						/* bimodal predictor was correct */
						if (*dir_update_ptr->rev_pmeta > 0)
							--*dir_update_ptr->rev_pmeta;
					}
				}
			}
		}
	}
	
	/* update L1 table if appropriate */
	/* L1 table is updated unconditionally for combining predictor too */
	if ((MD_OP_FLAGS(op) & (F_CTRL|F_UNCOND)) != (F_CTRL|F_UNCOND) &&
	((pred->class == BPredComb) || (pred->class == BPred2Level) || 
	(pred->class == BPredTSBP) || (pred->class == BPredCHBP) || (pred->class == BPredOB) || 
	(pred->class == BPredOHT) || (pred->class == BPredMBP))) {
		/* also update appropriate L1 history registers; Shift R or L depending on flow mode */
		/* if FRMTs enabled, only one shift is needed per mode */
		if (!flow_mode) { // FWD mode
			shift_history_left(pred->fwd_dirpred.twolev, fbaddr, taken);
			
			if (frmt==0) {
				shift_history_left(pred->rev_dirpred.twolev, rbaddr, taken);
			}
		} else { // REV mode
			shift_history_right(pred->fwd_dirpred.twolev, fbaddr, taken);
			
			if (frmt==0) {
				shift_history_right(pred->rev_dirpred.twolev, rbaddr, taken);
			}
		}
	}
}