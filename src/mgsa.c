/**
 * @file mgsa.c
 *
 * @author Sebastian Bauer
 *
 * @note if you want to compile as standalone, invoke
 *
 *  "gcc mgsa.c -DSTANDALONE `R CMD config --cppflags` `R CMD config --ldflags` -o mgsa"
 *
 *  If you want to test using R something like
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);mgsa(list(c(1,2),c(3)),3,o=c(1,2),steps=100)" | R --vanilla)"
 *  or (for invoking the function directly)
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_mcmc\",list(c(1,2),c(3)),3,o=c(1,2),4,5,6,steps=100)" | R --vanilla)"
 *  or (for the test function)
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_test\")" | R --vanilla)"
 */
#include <stdio.h>
#include <stdint.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Enable debugging */
//#define DEBUG

/* Define if file should be compiled as standalone (mainly for testing purposes) */
/* #define STANDALONE */

struct summary_for_cont_var
{
	double min;
	double max;
	int num_of_discrete_values;
	int values[0];
};

/**
 * Creates the structure needed to summarize a continuous variable.
 *
 * @param min
 * @param max
 * @param number_of_discrete_values
 * @return
 */
struct summary_for_cont_var *new_summary_for_cont_var(double min, double max, int number_of_discrete_values)
{
	struct summary_for_cont_var *sum;

	sum = (struct summary_for_cont_var *)R_alloc(1,sizeof(*sum) + number_of_discrete_values * sizeof(int));
	if (sum)
	{
		sum->min = min;
		sum->max = max;
		sum->num_of_discrete_values = number_of_discrete_values;
		memset(sum->values,0,number_of_discrete_values * sizeof(int));
	}

	return sum;
}

/**
 * Adds a new value to a summary.
 *
 * @param sum
 * @param val
 */
void add_to_summary(struct summary_for_cont_var *sum, double val)
{
	int slot;

	val -= sum->min;
	slot = val * sum->num_of_discrete_values / sum->max;
	if (slot < 0) slot=0;
	if (slot >= sum->num_of_discrete_values) slot = sum->num_of_discrete_values - 1;

	sum->values[slot]++;
}

/**
 * Returns a R representation of a summary.
 *
 * @param sum
 *
 * @note you must UNPROTECT one variable, after calling this function.
 */
static SEXP create_R_representation_of_summary(struct summary_for_cont_var *sum)
{
	int i;
	double w;

	SEXP l;
	SEXP l_names;
	SEXP breaks;
	SEXP counts;

	w = (sum->max - sum->min)/sum->num_of_discrete_values;

	PROTECT(l = allocVector(VECSXP,2));
	PROTECT(l_names = allocVector(STRSXP,2));
	PROTECT(counts = allocVector(REALSXP,sum->num_of_discrete_values));
	PROTECT(breaks = allocVector(REALSXP,sum->num_of_discrete_values));

	for (i=0;i<sum->num_of_discrete_values;i++)
		REAL(breaks)[i] = sum->min + w*i;
	for (i=0;i<sum->num_of_discrete_values;i++)
		REAL(counts)[i] = sum->values[i];

	SET_STRING_ELT(l_names,0,mkChar("breaks"));
	SET_STRING_ELT(l_names,1,mkChar("counts"));
	SET_VECTOR_ELT(l,0,breaks);
	SET_VECTOR_ELT(l,1,counts);
	setAttrib(l,R_NamesSymbol,l_names);


	UNPROTECT(3);

	return l;
}

struct context
{
	/** @brief Number of sets */
	int number_of_sets;

	/** @brief Sizes of sets */
	int *sizes_of_sets;

	/** @brief The definitions of the sets */
	int **sets;

	/** @brief Array of size "number_of_sets" indicating the state of the sets (active or not active) */
	int *sets_active;

	/** @brief number of active sets */
	int number_of_inactive_sets;

	/**
	 *  @brief contains indices of sets that are inactive and active in form of a partition.
	 *
	 *  The pivot element is given by number_of_inactive_sets.
	 **/
	int *set_partition;

	/** @brief array that reflects the position of elements in set_inactive_list */
	int *position_of_set_in_partition;

	/** @brief Number of observable */
	int number_of_observable;

	/** @brief Array indicating whether a observable is active or not */
	int *observable;

	/** @brief Array, indicating how much sets lead to the activation of a hidden entry */
	int *hidden_count;

	double alpha_fixed;
	double beta_fixed;
	double p_fixed;

	int n00;
	int n01;
	int n10;
	int n11;
	double alpha;
	double beta;
	double p;

	/* Proposal related */
	int proposal_toggle;
	int proposal_s1;
	int proposal_s2;
	double old_alpha;
	double old_beta;
	double old_p;

	/* Summary related */
	uint64_t *sets_activity_count;
	struct summary_for_cont_var *alpha_summary;
	struct summary_for_cont_var *beta_summary;
	struct summary_for_cont_var *p_summary;
};

/**
 * Initialize the context.
 *
 * @param context defines the context which is to be initalized
 * @param sets
 * @param sizes_of_sets
 * @param number_of_sets
 * @param n
 * @param o
 * @param lo
 * @return
 */
static int init_context(struct context *cn, int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo)
{
	int i;

	cn->number_of_sets = number_of_sets;
	cn->sets = sets;
	cn->sizes_of_sets = sizes_of_sets;
	cn->number_of_observable = n;

	/* Do the lazy stuff, i.e., allocate memory and intialize with meaningful defaults */
	if (!(cn->sets_active = (int*)R_alloc(number_of_sets,sizeof(cn->sets_active[0]))))
		goto bailout;
	memset(cn->sets_active,0,number_of_sets * sizeof(cn->sets_active[0]));
	if (!(cn->set_partition = (int*)R_alloc(number_of_sets,sizeof(cn->set_partition[0]))))
		goto bailout;
	if (!(cn->position_of_set_in_partition = (int*)R_alloc(number_of_sets,sizeof(cn->position_of_set_in_partition[0]))))
		goto bailout;
	for (i=0;i<number_of_sets;i++)
	{
		cn->set_partition[i] = i;
		cn->position_of_set_in_partition[i] = i;
	}
	cn->number_of_inactive_sets = number_of_sets;

	if (!(cn->hidden_count = (int*)R_alloc(n,sizeof(cn->hidden_count[0]))))
		goto bailout;
	memset(cn->hidden_count,0,n * sizeof(cn->hidden_count[0]));

	if (!(cn->observable = (int*)R_alloc(n,sizeof(cn->observable[0]))))
		goto bailout;
	memset(cn->observable,0,n * sizeof(cn->observable[0]));
	for (i=0;i<lo;i++)
		cn->observable[o[i]] = 1;

	/* Summary related */
	if (!(cn->sets_activity_count = (uint64_t *)R_alloc(n,sizeof(cn->sets_activity_count[0]))))
		goto bailout;
	memset(cn->sets_activity_count,0,n * sizeof(cn->sets_activity_count[0]));
	if (!(cn->alpha_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
	if (!(cn->beta_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
	if (!(cn->p_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;


	/* Initially, no set is active, hence all observations that are true are false positive... */
	cn->n10 = lo;
	/* while the rest are true negative */
	cn->n00 = n - lo;
	cn->n01 = cn->n11 = 0;

	cn->alpha = 0.10;
	cn->beta = 0.25;
	cn->p = 1.0 / number_of_sets;

	return 1;

bailout:
	return 0;
}

/**
 * @brief A hidden member has been newly activated.
 *
 * @param cn
 * @param member
 */
static void hidden_member_activated(struct context *cn, int member)
{
	if (cn->observable[member])
	{
		/* observation is true positive rather than false positive
		 * (first digit represents the observation, the second the hidden state) */
		cn->n11++;
		cn->n10--;
	} else
	{
		/* observation is false negative rather than true negative */
		cn->n01++;
		cn->n00--;
	}
}

/**
 * @brief A hidden member has been newly deactivated.
 *
 * @param cn
 * @param member
 */
static void hidden_member_deactivated(struct context *cn, int member)
{
	if (cn->observable[member])
	{
		/* just think about it... */
		cn->n11--;
		cn->n10++;
	} else
	{
		cn->n01--;
		cn->n00++;
	}

}

/**
 * @brief add a new set.
 *
 * @param cn
 * @param to_add
 */
static void add_set(struct context *cn, int to_add)
{
	int i; /* the usual dummy */

	if (cn->sets_active[to_add]) return;
	cn->sets_active[to_add] = 1;

	/* Go through all associations of the term and increase the activation count */
	for (i=0;i<cn->sizes_of_sets[to_add];i++)
	{
		int member = cn->sets[to_add][i];
		if (!cn->hidden_count[member])
		{
			/* A not yet active member is about to be activated */
			hidden_member_activated(cn,member);
		}
		cn->hidden_count[member]++;
	}

	/* Move the added set from the 0 partition to the 1 partition (it essentially becomes
	 * new first element of the 1 element, while the last 0 element gets its original position) */
	cn->number_of_inactive_sets--;
	if (cn->number_of_inactive_sets != 0)
	{
		int pos = cn->position_of_set_in_partition[to_add];
		int e0 = cn->set_partition[cn->number_of_inactive_sets];
		/* Move last element in the partition to left */
		cn->set_partition[pos] = e0;
		cn->position_of_set_in_partition[e0] = pos;
		/* Let be the newly added term the first in the partition */
		cn->set_partition[cn->number_of_inactive_sets] = to_add;
		cn->position_of_set_in_partition[to_add] = cn->number_of_inactive_sets;
	}

#ifdef DEBUG
	{
		printf("Add of %d produced following sets to be inactive:\n",to_add);
		if (cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=0;i<cn->number_of_inactive_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}
		printf("Add of %d produced following sets to be active:\n",to_add);
		if (cn->number_of_sets - cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}

	}
#endif
}

/**
 * @brief remove a set.
 *
 * @param cn
 * @param to_remove
 */
static void remove_set(struct context *cn, int to_remove)
{
	int i;

	if (!cn->sets_active[to_remove]) return;
	cn->sets_active[to_remove] = 0;

	/* Go through all associations of the term and decrease the activation count */
	for (i=0;i<cn->sizes_of_sets[to_remove];i++)
	{
		int member = cn->sets[to_remove][i];
		if (cn->hidden_count[member] == 1)
		{
			/* A active member is about to be deactivated  */
			hidden_member_deactivated(cn,member);
		}
		cn->hidden_count[member]--;
	}

	/* Converse of above. Here the removed set, which belonged to the 1 partition,
	 * is moved at the end of the 0 partition while the element at that place is
	 * pushed to the original position of the removed element. */
	if (cn->number_of_inactive_sets != (cn->number_of_sets - 1))
	{
		int pos = cn->position_of_set_in_partition[to_remove];
		int e1 = cn->set_partition[cn->number_of_inactive_sets];
		cn->set_partition[pos] = e1;
		cn->position_of_set_in_partition[e1] = pos;
		cn->set_partition[cn->number_of_inactive_sets] = to_remove;
		cn->position_of_set_in_partition[to_remove] = cn->number_of_inactive_sets;
	}
	cn->number_of_inactive_sets++;

#ifdef DEBUG
	{
		printf("Remove of %d produced following sets to be inactive:\n",to_remove);
		if (cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=0;i<cn->number_of_inactive_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}
		printf("Remove of %d produced following sets to be active:\n",to_remove);
		if (cn->number_of_sets - cn->number_of_inactive_sets == 0)
			printf(" empty list\n");
		for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		{
			int elem = cn->set_partition[i];
			printf(" %d at pos %d (should be %d)\n",elem,cn->position_of_set_in_partition[elem],i);
		}

	}
#endif

}

/**
 * @brief Switches the state of the given set.
 *
 * @param cn defines the context on which to operate
 * @param to_switch defines the set which should be switched.
 */
static void toggle_state(struct context *cn, int to_switch)
{
	int new_state;

	new_state = !cn->sets_active[to_switch];

	if (new_state)
		add_set(cn,to_switch);
	else
		remove_set(cn,to_switch);
}

/**
 * Returns the current alpha.
 *
 * @param cn
 * @return
 */
static double get_alpha(struct context *cn)
{
	if (cn->alpha_fixed > 0.0 && cn->alpha_fixed < 1)
		return cn->alpha_fixed;
	return cn->alpha;
}

/**
 * Returns the current beta.
 *
 * @param cn
 * @return
 */
static double get_beta(struct context *cn)
{
	if (cn->beta_fixed > 0.0 && cn->beta_fixed < 1)
		return cn->beta_fixed;
	return cn->beta;
}

/**
 * Returns the current p.
 *
 * @param cn
 * @return
 */
static double get_p(struct context *cn)
{
	if (cn->p_fixed > 0.0 && cn->p_fixed < 1)
		return cn->p_fixed;
	return cn->p;
}

/**
 * Returns the size of the on/off neighborhood.
 * @param cn
 * @return
 */
static uint64_t get_neighborhood_size(struct context *cn)
{
	/* First part accounts for the toggles, second part for the exchanges */
	return cn->number_of_sets + cn->number_of_inactive_sets * (cn->number_of_sets - cn->number_of_inactive_sets);

}

/**
 * @brief return the score of the current context.
 *
 * @param cn
 * @return
 */
static double get_score(struct context *cn)
{
	double alpha = get_alpha(cn);
	double beta = get_beta(cn);
	double p = get_p(cn);

	double score = log(alpha) * cn->n10 + log(1.0-alpha) * cn->n00 + log(1-beta)* cn->n11 + log(beta)*cn->n01;

	/* apply prior */
	score += p * (cn->number_of_sets - cn->number_of_inactive_sets) + (1.0-p) * cn->number_of_inactive_sets;

	return score;

}

/**
 * @brief Proposes a new state which can be undone via undo_proposal()
 * @param cn
 */
static void propose_state(struct context *cn)
{
	uint64_t possibilities = get_neighborhood_size(cn);

	cn->proposal_toggle = -1;
	cn->proposal_s1 = -1;
	cn->proposal_s2 = -1;
	cn->old_alpha = -1;
	cn->old_beta = -1;
	cn->old_p = -1;

	if (unif_rand() < 0.5)
	{
		uint32_t proposal = (double)(unif_rand() * possibilities);

		if (proposal < cn->number_of_sets)
		{
			/* on/off */
			cn->proposal_toggle = proposal;
			toggle_state(cn,proposal);
		}	else
		{
			int active_term_pos;
			int inactive_term_pos;

			proposal -= cn->number_of_sets;

			active_term_pos = (int)(proposal / cn->number_of_inactive_sets) +  cn->number_of_inactive_sets;
			inactive_term_pos = (int)(proposal % cn->number_of_inactive_sets);

			cn->proposal_s1 = cn->set_partition[active_term_pos];
			cn->proposal_s2 = cn->set_partition[inactive_term_pos];

			toggle_state(cn,cn->proposal_s1);
			toggle_state(cn,cn->proposal_s2);
		}
	} else
	{
		double which_param = unif_rand();
		if (which_param < (1.0/3.0))
		{
			cn->old_alpha = cn->alpha;
			cn->alpha = unif_rand();
		} else if (which_param < (2.0/3.0))
		{
			cn->old_beta = cn->beta;
			cn->beta = unif_rand();
		} else
		{
			cn->old_p = cn->p;
			cn->p = unif_rand() / 8; /* FIXME: */
		}
	}
}

/**
 * Undoes last proposal of propse_state().
 *
 * @param cn
 */
static void undo_proposal(struct context *cn)
{
	if (cn->proposal_toggle != -1) toggle_state(cn,cn->proposal_toggle);
	if (cn->proposal_s1 != -1)
	{
		toggle_state(cn,cn->proposal_s1);
		toggle_state(cn,cn->proposal_s2);
	}
	if (cn->old_alpha >= 0.0) cn->alpha = cn->old_alpha;
	if (cn->old_beta >= 0.0) cn->beta = cn->old_beta;
	if (cn->old_p >= 0.0) cn->p = cn->old_p;
}

/**
 * Records the current activity of the sets.
 *
 * @param cn
 */
static void record_activity(struct context *cn)
{
	int i;

	/* Remember that sets that are active are stored in the  partition */
	for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		cn->sets_activity_count[cn->set_partition[i]]++;

	add_to_summary(cn->alpha_summary,get_alpha(cn));
	add_to_summary(cn->beta_summary,get_beta(cn));
	add_to_summary(cn->p_summary,get_p(cn));
}

struct result
{
	double *marg_set_activity;
	struct summary_for_cont_var *alpha_summary;
	struct summary_for_cont_var *beta_summary;
	struct summary_for_cont_var *p_summary;
};

/**
 * The work horse.
 *
 * @param sets pointer to the sets. Sets a made of observable entities.
 * @param sizes_of_sets specifies the length of each set
 * @param number_of_sets number of sets (length of sets)
 * @param n the number of observable entities.
 * @param o indices of entities which are "on" (0 based).
 * @param lo length of o
 * @param number_of_steps defines the number of mcmc steps to be performed
 */
static struct result do_mgsa_mcmc(int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo, int64_t number_of_steps, double alpha, double beta, double p)
{
	int i;
	int64_t step;
	uint64_t neighborhood_size;
	struct context cn;
	struct result res;
	double score;

#ifdef DEBUG

	for (i=0;i<number_of_sets;i++)
	{
		int j;

		printf(" %d: ", sizes_of_sets[i]);
		for (j=0;j<sizes_of_sets[i];j++)
			printf("%d ", sets[i][j]);
		printf("\n");
	}

	printf(" n=%d\n",n);

	for (i=0;i<lo;i++)
	{
		printf(" o[%d]=%d\n",i,o[i]);
	}

	printf(" alpha=%g beta=%g p=%g\n",alpha,beta,p);
#endif
	memset(&res,0,sizeof(res));

	GetRNGstate();

	if (!init_context(&cn,sets,sizes_of_sets,number_of_sets,n,o,lo))
		goto bailout;

	cn.alpha_fixed = alpha;
	cn.beta_fixed = beta;
	cn.p_fixed = p;

	score = get_score(&cn);
	neighborhood_size = get_neighborhood_size(&cn);

#ifdef DEBUG
	printf("score=%g\n",score);
#endif

	for (step=0;step<number_of_steps;step++)
	{
		double new_score;
		double accept_probability;
		double u;
		uint64_t new_neighborhood_size;

		propose_state(&cn);
		new_score = get_score(&cn);
		new_neighborhood_size = get_neighborhood_size(&cn);

		accept_probability = exp(new_score - score) * (double)neighborhood_size / (double)new_neighborhood_size; /* last quotient is the hasting ratio */
#ifdef DEBUG
		printf("score = %g new_score=%g alpha=%g beta=%g p=%g\n",score,new_score,get_alpha(&cn),get_beta(&cn),get_p(&cn));
#endif

		u = unif_rand();
		if (u >= accept_probability)
		{
			undo_proposal(&cn);
		} else
		{
			score = new_score;
			neighborhood_size = new_neighborhood_size;
		}

		record_activity(&cn);
	}

#ifdef DEBUG
	{
		int i;

		for (i=0;i<cn.number_of_sets;i++)
		{
			printf(" Set %d: %g\n",i, (double)cn.sets_activity_count[i] / (double)number_of_steps);
		}
	}
#endif

	PutRNGstate();

	if (!(res.marg_set_activity = (double*)R_alloc(number_of_sets,sizeof(res.marg_set_activity[0]))))
		goto bailout;

	for (i=0;i<cn.number_of_sets;i++)
		res.marg_set_activity[i] = (double)cn.sets_activity_count[i] / (double)number_of_steps;

	res.alpha_summary = cn.alpha_summary;
	res.beta_summary = cn.beta_summary;
	res.p_summary = cn.p_summary;

bailout:
	return res;
}

/**
 * @param n the number of observable entities.
 * @param o (1 based)
 *
 * @note TODO: Check whether sets are real sets.
 */
SEXP mgsa_mcmc(SEXP sets, SEXP n, SEXP o, SEXP alpha, SEXP beta, SEXP p, SEXP steps)
{
	int *xo,*no,lo;
	int **nsets, *lset, lsets;
	int i,j;

	SEXP res = NULL_USER_OBJECT;

	if (LENGTH(n) != 1)
		error("Parameter 'n' needs to be atomic!");

	if (LENGTH(steps) != 1)
		error("Parameter 'steps' needs to be atomic!");

	PROTECT(n = AS_INTEGER(n));
	PROTECT(o = AS_INTEGER(o));
	PROTECT(sets = AS_LIST(sets));
	PROTECT(steps = AS_INTEGER(steps));
	PROTECT(alpha = AS_NUMERIC(alpha));
	PROTECT(beta = AS_NUMERIC(beta));
	PROTECT(p = AS_NUMERIC(p));

	/* Observations */
	xo = INTEGER_POINTER(o);
	lo = LENGTH(o);
	if (!(no =  (int*)R_alloc(lo,sizeof(no[0])))) /* R's garbage collection takes care of it */
		goto bailout;

	/* Turn 1 based observation indices into 0 based ones */
	for (i=0;i<lo;i++)
	{
		no[i] = xo[i] - 1; /* 1 based to 0 based */

		if (no[i] >= INTEGER_VALUE(n))
		{
			error("Observation index to high!");
			goto bailout;
		}
	}

	/* Set associations. Turn 1 based observations indices into 0 based ones */
	lsets = LENGTH(sets);
	if (!(lset = (int*)R_alloc(lsets,sizeof(lset[0]))))
		goto bailout;
	if (!(nsets = (int**)R_alloc(lsets,sizeof(nsets[0]))))
		goto bailout;
	for (i=0;i<lsets;i++)
	{
		int *xset;

		SEXP set = VECTOR_ELT(sets, i);

		PROTECT(set = AS_INTEGER(set));
		lset[i] = LENGTH(set);
		xset = INTEGER_POINTER(set);
		if (!(nsets[i] = (int*)R_alloc(lset[i],sizeof(nsets[i][0]))))
		{
			UNPROTECT(1);
			goto bailout;
		}

		xset = INTEGER_POINTER(set);
		for (j=0;j<lset[i];j++)
		{
			nsets[i][j] = xset[j] - 1; /* 1 based to 0 based */
			if (nsets[i][j] >= INTEGER_VALUE(n))
			{
				error("Set index to high (must not exceed 'n')");
				goto bailout;
			}
		}
		UNPROTECT(1);
	}

	/* Create the result. TODO Use a list */
	{
		struct result r;
		SEXP names;

		PROTECT(res = allocVector(VECSXP,4));
		PROTECT(names = allocVector(STRSXP,4));

		r = do_mgsa_mcmc(nsets, lset, lsets, INTEGER_VALUE(n),no,lo,INTEGER_VALUE(steps),DOUBLE_DATA(alpha)[0],DOUBLE_DATA(beta)[0],DOUBLE_DATA(p)[0]);

		if (r.marg_set_activity)
		{
			SEXP marg;
			PROTECT(marg = allocVector(REALSXP,lsets));

			for (i=0;i<lsets;i++)
				REAL(marg)[i] = r.marg_set_activity[i];

			SET_VECTOR_ELT(res,0,marg);
			SET_STRING_ELT(names,0,mkChar("marg"));

			UNPROTECT(1);
		}

		{
			SEXP alpha = create_R_representation_of_summary(r.alpha_summary);

			SET_VECTOR_ELT(res,1,alpha);
			SET_STRING_ELT(names,1,mkChar("alpha"));

			UNPROTECT(1);
		}

		{
			SEXP beta = create_R_representation_of_summary(r.beta_summary);

			SET_VECTOR_ELT(res,2,beta);
			SET_STRING_ELT(names,2,mkChar("beta"));

			UNPROTECT(1);
		}

		{
			SEXP p = create_R_representation_of_summary(r.p_summary);

			SET_VECTOR_ELT(res,3,p);
			SET_STRING_ELT(names,3,mkChar("p"));

			UNPROTECT(1);
		}


		setAttrib(res,R_NamesSymbol,names);
		UNPROTECT(2);
	}

bailout:
	UNPROTECT(7);
	return res;
}


static void print_context(struct context *cn)
{
	printf("n00=%d n01=%d n10=%d n11=%d num_active=%d\n",cn->n00,cn->n01,cn->n10,cn->n11,cn->number_of_sets - cn->number_of_inactive_sets);
}

/**
 * A simple test function for lowlevel functions.
 *
 * @return
 */
SEXP mgsa_test(void)
{
	struct context cn;
	int t1[] = {0,1};
	int t2[] = {1,2};
	int o[] = {0,1};

	int *sets[] = {t1,t2};
	int sizes_of_sets[] = {sizeof(t1)/sizeof(t1[0]),sizeof(t2)/sizeof(t2[0])};

	init_context(&cn,sets,sizes_of_sets,sizeof(sets)/sizeof(sets[0]),3,o,sizeof(o)/sizeof(o[0]));
	printf("no active term: ");print_context(&cn);
	add_set(&cn,0); printf("t1 is active: ");print_context(&cn);
	remove_set(&cn,0);add_set(&cn,1); printf("t2 is active: ");print_context(&cn);
	add_set(&cn,0); printf("t1,t2 is active: ");print_context(&cn);

	return NULL_USER_OBJECT;
}


R_CallMethodDef callMethods[] = {
   {"mgsa_mcmc", (DL_FUNC)&mgsa_mcmc, 7},
   {"mgsa_test", (DL_FUNC)&mgsa_test, 0},
   {NULL, NULL, 0}
};

/**
 * Called by R as the initialization routine.
 *
 * @param info
 */
void R_init_myLib(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

/**
 * Called by R when the object is unloaded.
 *
 * @param info
 */
void R_unload_mylib(DllInfo *info)
{
  /* Release resources. */
}

#ifdef STANDALONE

int main(void)
{
	int t1[] = {0,1};
	int t2[] = {1,2};
	int o[] = {0,1};

	int *sets[] = {t1,t2};
	int sizes_of_sets[] = {sizeof(t1)/sizeof(t1[0]),sizeof(t2)/sizeof(t2[0])};

	do_mgsa_mcmc(sets, sizes_of_sets, sizeof(sets)/sizeof(sets[0]), 3, o, sizeof(o)/sizeof(o[0]));

}

#endif
