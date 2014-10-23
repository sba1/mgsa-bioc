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
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);sets<-list(c(1,2),c(3));names(sets)<-c(\"set 1\", \"set 2\");mgsa:::mgsa.trampoline(c(1,2),sets,3,steps=100,restarts=4)" | R --vanilla)"
 *  or (for invoking the function directly)
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_mcmc\",list(c(1,2),c(3)),3,o=c(1,2),4,5,6,steps=100)" | R --vanilla)"
 *  or (for the test function)
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_test\")" | R --vanilla)"
 */

/*#define DEBUG*/

#include <signal.h>
#include <stdio.h>
#include <stdint.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

#include "mt19937p/mt19937p.c"

/* Enable debugging */
//#define DEBUG

/* Define if file should be compiled as standalone (mainly for testing purposes) */
/* #define STANDALONE */

#define MIN(a,b) ((a)<(b)?(a):(b))

/** Global variable which is set to 1, if thread has been interrupted */
static int is_interrupted;

/**
 * Thread-safe variant of R_alloc().
 *
 * @return memory allocated by R_alloc().
 * @note do not allocate memory too often.
 */
void *ts_R_alloc(size_t n, int size)
{
	void *mem;

#pragma omp critical
	{
		mem = R_alloc(n,size);
	}
	return mem;
}

/** Redirect R_alloc. I know, this is ugly... */
#define R_alloc ts_R_alloc


/*************************************************************/

/**
 * Represents a prior for parameter. Very simple at the moment.
 */
struct parameter_prior
{
	/** @brief is this prior uniform continuous? (otherwise it is uniform discrete) */
	int uniform_continuous;
	double uniform_continous_lower;
	double uniform_continous_upper;

	/** @brief discrete values */
	double *values;

	/** @brief length of values */
	int number_of_states;
};


/**
 * Create a C representation of the parameter range from the
 * R specification.
 *
 * Uses R_alloc(), so R frees memory after C code has been returned.
 *
 * @param sexp
 * @param discrete defines whether the parameter should be handled as discrete.
 * @param param_name defines the name of the parameter (used for error messages only)
 * @return
 */
static struct parameter_prior *create_parameter_prior_from_R(SEXP sexp, int discrete, char *param_name)
{
	int i;
	struct parameter_prior *p;

	if (!(p = R_alloc(1,sizeof(*p))))
		return NULL;

	PROTECT(sexp = AS_NUMERIC(sexp));
	p->number_of_states = LENGTH(sexp);

	if (discrete)
	{
		p->uniform_continuous = 0;

		if (p->number_of_states == 0)
			error("Parameter '%s' has been requested to be discrete but no values were specified");

		if (!(p->values = R_alloc(p->number_of_states,sizeof(p->values[0]))))
			return NULL;

		for (i=0;i<p->number_of_states;i++)
		{
			p->values[i] = REAL(sexp)[i];

			if (isnan(p->values[i]))
				error("NAs are not supported for parameter '%s'",param_name);
		}
	} else
	{
		p->uniform_continuous = 1;

		if (p->number_of_states == 0)
		{
			/* Uniform across the full range */
			p->uniform_continous_lower = 0.0;
			p->uniform_continous_upper = 1.0;
		} else
		{
			if (p->number_of_states == 1)
			{
				/* The fixed case */
				p->uniform_continous_lower =  REAL(sexp)[0];
				p->uniform_continous_upper =  REAL(sexp)[0];

				/* Special case */
				if (isnan(p->uniform_continous_lower))
				{
					p->uniform_continous_lower = 0.0;
					p->uniform_continous_upper = 1.0;
				}
			} else if (p->number_of_states == 2)
			{
				p->uniform_continous_lower =  REAL(sexp)[0];
				p->uniform_continous_upper =  REAL(sexp)[1];
			} else
			{
				error("Only one continuous range is supported at the moment!");
			}

			if (isnan(p->uniform_continous_lower) || isnan(p->uniform_continous_upper))
				error("NAs are not supported for parameter '%s'",param_name);

			/* Check for validity */
			if (p->uniform_continous_lower < 0 || p->uniform_continous_upper > 1)
				error("Range values for '%s' have to be between 0 and 1 (are %f and %f)",param_name,p->uniform_continous_lower,p->uniform_continous_upper);

			/* It is fine if the range is specified in the wrong order */
			if (p->uniform_continous_lower > p->uniform_continous_upper)
			{
				double help;
				help = p->uniform_continous_upper;
				p->uniform_continous_upper = p->uniform_continous_lower;
				p->uniform_continous_lower = help;
			}
		}
	}

	UNPROTECT(1);

	return p;
}

/**
 * Represents a basic sample. Only useful in conjunction with a parameter_prior structure.
 * Can be assigned.
 */
struct prior_sample
{
	union
	{
		int discrete_index;
		double continuous_value;
	} u;
};

/**
 * Samples from the given parameter prior.
 *
 * @param sample holds the sample
 * @param prior from which prior to sample
 * @param mt source for randomness
 *
 * @returns an opaque sample
 */
static void parameter_prior_sample(struct prior_sample *sample, struct parameter_prior *prior, struct mt19937p *mt)
{
	double rnd = genrand(mt);

	if (prior->uniform_continuous)
	{
		sample->u.continuous_value = prior->uniform_continous_lower + rnd * (prior->uniform_continous_upper - prior->uniform_continous_lower);
	} else
	{
		sample->u.discrete_index = rnd * prior->number_of_states;
		if ( sample->u.discrete_index >= prior->number_of_states )
			sample->u.discrete_index = prior->number_of_states - 1;
	}
}

/**
 * Returns the real value represented by the sample.
 *
 * @param sample
 * @param prior
 * @return
 */
static double get_prior_sample_value(struct prior_sample *sample, struct parameter_prior *prior)
{
	if (prior->uniform_continuous)
		return sample->u.continuous_value;
	return	prior->values[sample->u.discrete_index];
}


/*************************************************************/

struct summary
{
	/** The parameter for which this is a summary */
	struct parameter_prior *pdsc;

	/** number of values */
	int num_of_discrete_values;

	/** The breaks, length as specified as above */
	double *breaks;

	/** The actual value counts with same length as breaks */
	int *values;

	/**
	 * Maps indices of the values field of pdsc to the corresponding
	 * index in this values. Invalid, if the described parameter
	 * is continuous.
	 */
	int *dmap;
};

/**
 * Allocates memory in sum to hold number_of_breaks breaks. Element
 * values is initialized while breaks is left uninitialized.
 *
 * @param sum
 * @param number_of_breaks
 * @return
 */
static int init_summary_for_breaks(struct summary *sum, int number_of_breaks)
{
	sum->num_of_discrete_values = number_of_breaks;
	if (!(sum->values = R_alloc(number_of_breaks,sizeof(sum->values[0]))))
		return 0;
	memset(sum->values,0,number_of_breaks * sizeof(sum->values[0]));
	if (!(sum->breaks = R_alloc(number_of_breaks,sizeof(sum->breaks[0]))))
		return 0;
	return 1;
}

/**
 * This call creates a new summary description for the given parameter description and the
 * given breaks.
 *
 * @param pdsc
 * @param breaks
 * @return
 */
static struct summary *create_summary_for_param_from_R(struct parameter_prior *pdsc, SEXP breaks)
{
	int i;
	int default_range;
	int number_of_discrete_values;
	struct summary *sum;

	if (!(sum = (struct summary *)R_alloc(1,sizeof(*sum))))
		error("Couldn't allocate memory for summary statistics");
	memset(sum,0,sizeof(*sum));
	sum->pdsc = pdsc;

	PROTECT((breaks = AS_NUMERIC(breaks)));

	number_of_discrete_values = LENGTH(breaks);
	if (number_of_discrete_values == 0 || (number_of_discrete_values == 1 && isnan(REAL(breaks)[0])))
		default_range = 1;
	else default_range = 0;

	if (pdsc->uniform_continuous)
	{
		double min = pdsc->uniform_continous_lower;
		double max = pdsc->uniform_continous_upper;

		if (default_range)
			number_of_discrete_values = 11;

		if (!(init_summary_for_breaks(sum,number_of_discrete_values)))
			error("Couldn't allocate memory!");

		for (i=0;i<number_of_discrete_values;i++)
			sum->breaks[i] = min + (max - min) / (double)(number_of_discrete_values - 1) * i;
	} else
	{
		number_of_discrete_values = pdsc->number_of_states;

		if (!(init_summary_for_breaks(sum,number_of_discrete_values)))
			error("Couldn't allocate memory!");

		if (!(sum->dmap = R_alloc(number_of_discrete_values,sizeof(sum->dmap[0]))))
			error("Couldn't allocate memory!");

		if (!default_range)
		{
			/* A custom range was specified.
			 *
			 * FIXME: Currently, we don't support values that differ from the real state space
			 */
			if (number_of_discrete_values != pdsc->number_of_states)
				error("Number of breaks (%d) must equals the number of discrete states (%d)!",number_of_discrete_values,pdsc->number_of_states);

			for (i=0;i<number_of_discrete_values;i++)
			{
				if (REAL(breaks)[i] != sum->pdsc->values[i])
					error("Breaks must match states of of discrete values!");

				/* Transfer custom range. Indices map in 1:1 fashion (see FIXME above). */
				sum->dmap[i] = i;
				sum->breaks[i] = REAL(breaks)[i];
			}
		} else
		{
			for (i=0;i<number_of_discrete_values;i++)
			{
				sum->dmap[i] = i;
				sum->breaks[i] = (1.0) / (double)(number_of_discrete_values - 1) * i;
			}
		}
	}

	UNPROTECT(1);
	return sum;
}

/**
 * Duplicates a given summary. The pdsc is not duplicated.
 *
 * @param sum
 * @return
 */
static struct summary *duplicate_summary(struct summary *sum)
{
	struct summary *new_sum = R_alloc(1,sizeof(*sum));
	int num_of_discrete_values = sum->num_of_discrete_values;

	if (!new_sum)
		return NULL;

	*new_sum = *sum;

	if (sum->breaks)
	{
		if (!(new_sum->breaks = (double*)R_alloc(1,sizeof(sum->breaks[0])*num_of_discrete_values)))
			return NULL;
		memcpy(new_sum->breaks,sum->breaks,sizeof(sum->breaks[0])*num_of_discrete_values);
	}

	if (sum->values)
	{
		if (!(new_sum->values = (int*)R_alloc(1,sizeof(sum->values[0])*num_of_discrete_values)))
			return NULL;
		memcpy(new_sum->values,sum->values,sizeof(sum->values[0])*num_of_discrete_values);
	}

	if (sum->dmap)
	{
		if (!(new_sum->dmap = R_alloc(num_of_discrete_values,sizeof(sum->dmap[0]))))
			return NULL;
		memcpy(new_sum->dmap,sum->dmap,sizeof(sum->dmap[0])*num_of_discrete_values);
	}
	return new_sum;
}

/**
 * Adds a new sample for the given summary.
 *
 * @param sum the summary container
 * @param value the actual value
 * @param param the description
 */
void add_to_summary(struct summary *sum, struct prior_sample *value)
{
	int slot;

	if (sum->pdsc->uniform_continuous)
	{
		double val = get_prior_sample_value(value,sum->pdsc);

		val -= sum->pdsc->uniform_continous_lower;
		/* Remember that the last slot is reserved */
		slot = val * (sum->num_of_discrete_values -1) / (sum->pdsc->uniform_continous_upper - sum->pdsc->uniform_continous_lower);
		if (slot < 0) slot=0;
		if (slot >= sum->num_of_discrete_values)
			slot = sum->num_of_discrete_values - 1;
	} else
	{
		slot = sum->dmap[value->u.discrete_index];
	}
	sum->values[slot]++;
}


/**
 * Returns a R representation of a summary.
 *
 * @param sum
 * @param number_of_summaries
 *
 * @note you must UNPROTECT one variable, after calling this function.
 */
static SEXP create_R_representation_of_summary(struct summary **sum, int number_of_summaries)
{
	int num_of_discrete_values;
	int i,j;

	SEXP l;
	SEXP l_names;
	SEXP breaks;
	SEXP counts;

	/* At the moment, we take only the first summary */
	num_of_discrete_values = sum[0]->num_of_discrete_values;
	PROTECT(l = allocVector(VECSXP,2));

	PROTECT(l_names = allocVector(STRSXP,2));
	PROTECT(breaks = allocVector(REALSXP,num_of_discrete_values));
	PROTECT(counts = allocMatrix(REALSXP,num_of_discrete_values, number_of_summaries));

	for (i=0;i<num_of_discrete_values;i++)
		REAL(breaks)[i] = sum[0]->breaks[i];

	for (j=0;j<number_of_summaries;j++)
		for (i=0;i<num_of_discrete_values;i++)
			REAL(counts)[i+j*num_of_discrete_values] = sum[j]->values[i];

	SET_STRING_ELT(l_names,0,mkChar("breaks"));
	SET_STRING_ELT(l_names,1,mkChar("counts"));
	SET_VECTOR_ELT(l,0,breaks);
	SET_VECTOR_ELT(l,1,counts);
	setAttrib(l,R_NamesSymbol,l_names);

	UNPROTECT(3);
	return l;
}

/** @brief MGSA MCMC parameters */
struct mcmc_params
{
	/** @brief number of steps in MCMC */
	int64_t nsteps;
	/** @brief number of burn-in MCMC steps */
	int64_t nsteps_burnin;
	/** @brief nsteps_thin sample collecting period */
	int nsteps_thin;
	/** @brief the frequency of set active state flipping Gibbs step */
	double flip_freq;
};

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

	/** @brief description of the prior for alpha */
	struct parameter_prior *alpha_prior;

	/** @brief prior for beta */
	struct parameter_prior *beta_prior;

	/** @brief prior for p */
	struct parameter_prior *p_prior;

	/* first digit is the observation state, second the hidden state */
	int n00;
	int n01;
	int n10;
	int n11;

	/* The actual values of the parameters */
	struct prior_sample alpha;
	struct prior_sample beta;
	struct prior_sample p;

	/* Proposal related, used to undo proposal */
	int proposal_toggle;
	int proposal_s1;
	int proposal_s2;

	int proposal_which_parameter_changed; /**! @brief defines which parameter (0-2) has been changed */
	struct prior_sample old_alpha; /* In theory, we would need only one old parameter sample (at the moment we can only change one parameter) */
	struct prior_sample old_beta;
	struct prior_sample old_p;

	/* The remainder is summary related and updated in record_activity() */
	uint64_t nsamples; /**! @brief number of samples collected */
	uint64_t *sets_activity_count;
	struct summary *alpha_summary;
	struct summary *beta_summary;
	struct summary *p_summary;

	/** @brief represents the max log likelihood encountered so far */
	double max_score;

	/** @brief in which step max_score was recorded */
	int64_t max_score_step;

	double max_score_alpha;
	double max_score_beta;
	double max_score_p;
	int *max_score_sets_active;
	int max_score_sets_active_length;
};

/**
 * Initialize the context.
 *
 * @param context defines the context which is to be initialized
 * @param sets
 * @param sizes_of_sets
 * @param number_of_sets
 * @param n number of items
 * @param o
 * @param lo number of items that are observed as true (<=n)
 * @return
 */
static int init_context(struct context *cn, int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo)
{
	int i;
#ifdef DEBUG
	printf("init_context(number_of_sets=%d,n=%d,lo=%d)\n",number_of_sets,n,lo);
#endif

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

	if (!(cn->max_score_sets_active = (int*)R_alloc(number_of_sets,sizeof(cn->max_score_sets_active[0]))))
		goto bailout;

	/* Summary related, default values, some will be overwritten later */
	cn->nsamples = 0;
	if (!(cn->sets_activity_count = (uint64_t *)R_alloc(number_of_sets,sizeof(cn->sets_activity_count[0]))))
		goto bailout;
	memset(cn->sets_activity_count,0,number_of_sets * sizeof(cn->sets_activity_count[0]));
#if 0
	if (!(cn->alpha_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
	if (!(cn->beta_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
	if (!(cn->p_summary = new_summary_for_cont_var(0,1,10)))
		goto bailout;
#endif
	/* Initially, no set is active, hence all observations that are true are false positives... */
	cn->n10 = lo;
	/* ...while the rest are true negatives */
	cn->n00 = n - lo;
	cn->n01 = cn->n11 = 0;

	cn->max_score = -DBL_MAX;
	cn->max_score_sets_active_length = 0;

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

#ifdef DEBUG_VERBOSE
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
#ifdef DEBUG_VERBOSE
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
	return get_prior_sample_value(&cn->alpha,cn->alpha_prior);
}

/**
 * Returns the current beta.
 *
 * @param cn
 * @return
 */
static double get_beta(struct context *cn)
{
	return get_prior_sample_value(&cn->beta,cn->beta_prior);
}

/**
 * Returns the current p.
 *
 * @param cn
 * @return
 */
static double get_p(struct context *cn)
{
	return get_prior_sample_value(&cn->p,cn->p_prior);
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

	double score = log(alpha) * cn->n10 + log1p(-alpha) * cn->n00 + log1p(-beta)* cn->n11 + log(beta)*cn->n01;

	/* apply prior */
	score += log(p) * (cn->number_of_sets - cn->number_of_inactive_sets) + log1p(-p) * cn->number_of_inactive_sets;
	return score;
}

/**
 * @brief Proposes a new state which can be undone via undo_proposal()
 * @param cn
 * @param step current MCMC step number
 */
static void propose_state(struct context *cn, struct mcmc_params *params, struct mt19937p *mt, int64_t step)
{
	uint64_t possibilities = get_neighborhood_size(cn);

	cn->proposal_toggle = -1;
	cn->proposal_s1 = -1;
	cn->proposal_s2 = -1;

	if (step >= 0.5 * params->nsteps_burnin || genrand(mt) <= params->flip_freq)
	{
		/* toggle inactive/active states */
		uint32_t proposal = (double)(genrand(mt) * possibilities);

		if (proposal < cn->number_of_sets)
		{
			/* on/off for a single set */
			cn->proposal_toggle = proposal;
			toggle_state(cn,proposal);
		}	else
		{
			/* on/off for a pair of sets */
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
		/* FIXME: If a parameter is fixed it still is considered here (this wastes steps) */
		double which_param = genrand(mt);
		if (which_param < (1.0/3.0))
		{
			cn->old_alpha = cn->alpha;
			parameter_prior_sample(&cn->alpha,cn->alpha_prior, mt);
			cn->proposal_which_parameter_changed = 0;
		} else if (which_param < (2.0/3.0))
		{
			cn->old_beta = cn->beta;
			parameter_prior_sample(&cn->beta,cn->beta_prior, mt);
			cn->proposal_which_parameter_changed = 1;
		} else
		{
			cn->old_p = cn->p;
			parameter_prior_sample(&cn->p,cn->p_prior, mt);
			cn->proposal_which_parameter_changed = 2;
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
	else if (cn->proposal_s1 != -1)
	{
		toggle_state(cn,cn->proposal_s1);
		toggle_state(cn,cn->proposal_s2);
	} else
	{
		switch (cn->proposal_which_parameter_changed)
		{
		case	0: cn->alpha = cn->old_alpha; break;
		case	1: cn->beta = cn->old_beta; break;
		default: cn->p = cn->old_p; break;
		}
	}
}

/**
 * Records the current activity of the sets. Should be called on each step.
 *
 * @param cn
 * @param score specifies the current score (log)
 */
static void record_activity(struct context *cn, int64_t step, double score)
{
	int i,j;

	cn->nsamples++;

	/* Remember that sets that are active are stored in the second partition */
	for (i=cn->number_of_inactive_sets;i<cn->number_of_sets;i++)
		cn->sets_activity_count[cn->set_partition[i]]++;

	add_to_summary(cn->alpha_summary,&cn->alpha);
	add_to_summary(cn->beta_summary,&cn->beta);
	add_to_summary(cn->p_summary,&cn->p);

	if (score > cn->max_score)
	{
		/* Remember the score and the configuration */
		cn->max_score = score;
		cn->max_score_step = step;
		cn->max_score_alpha = get_alpha(cn);
		cn->max_score_beta = get_beta(cn);
		cn->max_score_p = get_p(cn);

		for (i=cn->number_of_inactive_sets,j=0;i<cn->number_of_sets;i++,j++)
			cn->max_score_sets_active[j] = cn->set_partition[i];
		cn->max_score_sets_active_length = cn->number_of_sets - cn->number_of_inactive_sets;
	}
}

/* This in essence replicates the summary stuff in context. TODO: unify this */
struct result
{
	uint64_t nsamples;

	double *marg_set_activity;
	struct summary *alpha_summary;
	struct summary *beta_summary;
	struct summary *p_summary;

	double max_score;
	double max_score_alpha;
	double max_score_beta;
	double max_score_p;
	int *max_score_sets_active;
	int max_score_sets_active_length;
};

/**
 * The work horse, the single MCMC run.
 *
 * @param sets pointer to the sets. Sets a made of observable entities.
 * @param sizes_of_sets specifies the length of each set
 * @param number_of_sets number of sets (length of sets)
 * @param n the number of observable entities.
 * @param o indices of entities which are "on" (0 based).
 * @param lo length of o
 * @param params MGSA MCMC params
 * @param alpha
 * @param beta
 * @param p
 * @param mt pre-seeded structure used for random number generation.
 */
static struct result do_mgsa_mcmc(int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo,
		struct mcmc_params *params,
		struct summary *alpha_summary, struct summary *beta_summary, struct summary *p_summary,
		struct mt19937p *mt)
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

	printf(" p=%d p=%d p=%d\n",alpha_prior->uniform_continuous,beta_prior->uniform_continuous,p_prior->uniform_continuous);
#endif
	memset(&res,0,sizeof(res));

	if (!init_context(&cn,sets,sizes_of_sets,number_of_sets,n,o,lo))
		goto bailout;

	if (!(alpha_summary = duplicate_summary(alpha_summary)))
		goto bailout;
	if (!(beta_summary = duplicate_summary(beta_summary)))
		goto bailout;
	if (!(p_summary = duplicate_summary(p_summary)))
		goto bailout;

	cn.alpha_summary = alpha_summary;
	cn.beta_summary = beta_summary;
	cn.p_summary = p_summary;

	cn.alpha_prior = alpha_summary->pdsc;
	cn.beta_prior = beta_summary->pdsc;
	cn.p_prior = p_summary->pdsc;

	/* Sample initial state */
	parameter_prior_sample(&cn.alpha,	cn.alpha_prior, mt);
	parameter_prior_sample(&cn.beta,	cn.beta_prior, mt);
	parameter_prior_sample(&cn.p,		cn.p_prior, mt);
	for (i=0;i<number_of_sets;i++) {
		if (genrand(mt) < 0.5) toggle_state(&cn,i);
	}

	score = get_score(&cn);
	neighborhood_size = get_neighborhood_size(&cn);

#ifdef DEBUG
	printf("score=%g\n",score);
#endif

	for (step=0;step<params->nsteps;step++)
	{
		double new_score;
		double accept_probability;
		double u;
		uint64_t new_neighborhood_size;

		/* Break, if we are interrupted (we cannot use R_CheckUserInterrupt() here,
		 * as this is not thread-safe */
		if (is_interrupted) break;

		propose_state(&cn, params, mt, step);
		new_score = get_score(&cn);
		new_neighborhood_size = get_neighborhood_size(&cn);

		accept_probability = exp(new_score - score) * (double)neighborhood_size / (double)new_neighborhood_size; /* last quotient is the hasting ratio */
#ifdef DEBUG
		printf("score = %g new_score=%g alpha=%g beta=%g p=%g\n",score,new_score,get_alpha(&cn),get_beta(&cn),get_p(&cn));
#endif

		u = genrand(mt);
		if (u >= accept_probability)
		{
			undo_proposal(&cn);
		} else
		{
			score = new_score;
			neighborhood_size = new_neighborhood_size;
		}

		if ( step >= params->nsteps_burnin && (step - params->nsteps_burnin) % params->nsteps_thin == 0 ) {
			record_activity(&cn, step, score);
		}
	}

#ifdef DEBUG
	{
		int i;

		for (i=0;i<cn.number_of_sets;i++)
		{
			printf(" Set %d: %g\n",i, (double)cn.sets_activity_count[i] / (double)cn.nsamples);
		}
	}
#endif

	res.nsamples = cn.nsamples;

	if (!(res.marg_set_activity = (double*)R_alloc(number_of_sets,sizeof(res.marg_set_activity[0]))))
		goto bailout;

	for (i=0;i<cn.number_of_sets;i++)
		res.marg_set_activity[i] = (double)cn.sets_activity_count[i] / (double)cn.nsamples;

	/* Note, we are using here stuff that is allocated in init_context() */
	res.alpha_summary = cn.alpha_summary;
	res.beta_summary = cn.beta_summary;
	res.p_summary = cn.p_summary;
	res.max_score = cn.max_score;
	res.max_score_alpha = cn.max_score_alpha;
	res.max_score_beta = cn.max_score_beta;
	res.max_score_p = cn.max_score_p;
	res.max_score_sets_active = cn.max_score_sets_active;
	res.max_score_sets_active_length = cn.max_score_sets_active_length;
	bailout:
	return res;
}

/**
 * Interrupt handler. Simply sets global is_interrupted flag.
 *
 * @param sig
 */
static void signal_handler(int sig)
{
	/* Set interruption flag */
	is_interrupted = 1;
}

/**
 *
 * @param sets (1 based)
 * @param n the number of observable entities.
 * @param o (1 based, just the indices)
 * @param alpha domain of the alpha parameter
 * @param beta domain of the beta parameter
 * @param p domain of the p parameter
 * @param discrete a logical (Boolean) vector which specifies the whether alpha,beta and p are discrete.
 * @param alpha_breaks breaks used for summary. Depending on type of alpha, certain restrictions apply.
 * @param beta_breaks breaks used for summary. Depending on type of beta, certain restrictions apply.
 * @param p_breaks  breaks used for summary.  Depending on type of p, certain restrictions apply.
 * @param steps number of steps in each Monte-Carlo Markov chain
 * @param burnin number of burn-in MCMC steps
 * @param thin sampling period, every thin-th step after the chain burn-in the parameters sample is collected
 * @param flip_freq the frequency of set flipping Gibbs step
 * @param restarts number of MCMC chains in a single thread
 * @param threads number of parallel MCMC threads
 * @param as (1 based, just indices as o, not that the shape for the return will change in that case).
 * @return
 */
SEXP mgsa_mcmc(SEXP sets, SEXP n, SEXP o,
		SEXP alpha, SEXP beta, SEXP p, SEXP discrete,
		SEXP alpha_breaks, SEXP beta_breaks, SEXP p_breaks,
		SEXP steps, SEXP burnin, SEXP thin, SEXP flip_freq,
		SEXP restarts, SEXP threads, SEXP as)
{
	struct parameter_prior *alpha_prior, *beta_prior, *p_prior;
	struct summary *alpha_summary, *beta_summary, *p_summary;
	int *xo,*no,lo;
	int *nas = NULL, las;
	int **nsets, *lset, lsets;
	int *ndiscrete;
	int i,j;

	/* Clear interruption flag */
	is_interrupted = 0;

	/* Call this function, when interruption occurs */
	signal(SIGINT,signal_handler);

	SEXP res = NULL_USER_OBJECT;

	if (LENGTH(n) != 1)
		error("Parameter 'n' needs to be atomic!");

	if (LENGTH(steps) != 1)
		error("Parameter 'steps' needs to be atomic!");

	if (LENGTH(burnin) != 1)
		error("Parameter 'burnin' needs to be atomic!");

	if (LENGTH(thin) != 1)
		error("Parameter 'thin' needs to be atomic!");

	if (LENGTH(flip_freq) != 1)
		error("Parameter 'flip_freq' needs to be atomic!");

	if (LENGTH(restarts) != 1)
		error("Parameter 'restarts' needs to be atomic!");

	if (LENGTH(threads) != 1)
		error("Parameter 'threads' needs to be atomic!");

	if (LENGTH(discrete) != 3)
		error("Parameter 'discrete' needs to be a vector of length 3!");

	PROTECT(n = AS_INTEGER(n));
	PROTECT(o = AS_INTEGER(o));
	PROTECT(sets = AS_LIST(sets));
	PROTECT(restarts = AS_INTEGER(restarts));
	PROTECT(threads = AS_INTEGER(threads));
	PROTECT(as = AS_INTEGER(as));
	PROTECT(discrete = AS_LOGICAL(discrete));
	PROTECT(alpha_breaks = AS_NUMERIC(alpha_breaks));
	PROTECT(beta_breaks = AS_NUMERIC(beta_breaks));
	PROTECT(p_breaks = AS_NUMERIC(p_breaks));

	if (INTEGER_VALUE(restarts) < 1)
	{
		error("Parameter 'restarts' needs to be larger than 0!");
		goto bailout;
	}

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
			error("Observation index given in o %d exceeds the total number of observations (n=%d)!",no[i],INTEGER_VALUE(n));
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
				UNPROTECT(1);
				goto bailout;
			}
		}
		UNPROTECT(1);
	}

	/* Convert active sets */
	las = LENGTH(as);
	if (las > 0)
	{
		nas = R_alloc(las,sizeof(nas[0]));

		for (i=0;i<las;i++)
		{
			nas[i] =  INTEGER_POINTER(as)[i] - 1;
			if (nas[i] >= lsets)
			{
				error("Parameter 'as' refers to a non-existing set %d.",nas[i]);
				goto bailout;
			}
		}
	}

	ndiscrete = LOGICAL(discrete);

	/* Make prior stuff */
	if (!(alpha_prior = create_parameter_prior_from_R(alpha,ndiscrete[0],"alpha")))
		goto bailout;
	if (!(beta_prior = create_parameter_prior_from_R(beta,ndiscrete[1],"beta")))
		goto bailout;
	if (!(p_prior = create_parameter_prior_from_R(p,ndiscrete[2],"p")))
		goto bailout;

	if (!(alpha_summary = create_summary_for_param_from_R(alpha_prior,alpha_breaks)))
		goto bailout;
	if (!(beta_summary = create_summary_for_param_from_R(beta_prior,beta_breaks)))
		goto bailout;
	if (!(p_summary = create_summary_for_param_from_R(p_prior,p_breaks)))
		goto bailout;

	if (nas)
	{
		/* Mode 1: Just calculate the score */
		struct context cn;
		double score;
		int i;

		if (!alpha_prior->uniform_continuous && alpha_prior->number_of_states != 1)
		{
			error("Parameter 'alpha' needs to be atomic if 'as' is given.");
			goto bailout;
		}

		if (!beta_prior->uniform_continuous && beta_prior->number_of_states != 1)
		{
			error("Parameter 'beta' needs to be atomic if 'as' is given.");
			goto bailout;
		}

		if (!p_prior->uniform_continuous && p_prior->number_of_states != 1)
		{
			error("Parameter 'p' needs to be atomic if 'as' is given.");
			goto bailout;
		}

		if (!init_context(&cn,nsets,lset,lsets,INTEGER_VALUE(n),no,lo))
			goto bailout;

		cn.alpha_summary = alpha_summary;
		cn.beta_summary = beta_summary;
		cn.p_summary = p_summary;

		cn.alpha_prior = alpha_prior;
		cn.beta_prior = beta_prior;
		cn.p_prior = p_prior;

		/* This simply assigns the user-desired parameter */
		cn.alpha.u.discrete_index = 0;
		cn.beta.u.discrete_index = 0;
		cn.p.u.discrete_index = 0;

		if (alpha_prior->uniform_continuous)
		{
			error("Parameter 'alpha' needs to be specified!");
			goto bailout;
		}

		if (beta_prior->uniform_continuous)
		{
			error("Parameter 'beta' needs to be specified!");
			goto bailout;
		}

		if (p_prior->uniform_continuous)
		{
			error("Parameter 'p' needs to be specified!");
			goto bailout;
		}

		for (i=0;i<las;i++)
			add_set(&cn,nas[i]);

		score = get_score(&cn);

		PROTECT(res = allocVector(REALSXP,1));
		REAL(res)[0] = score;
		UNPROTECT(1);
	} else
	{
		/* Mode 2: Perform MCMC. Here, we organize the creation of the result frame and the threading */
		struct mcmc_params params;
		struct result *r;
		int run;
		int have_margs, have_alphas, have_betas, have_ps;
		int irestarts = INTEGER_VALUE(restarts);

		struct summary *summaries[irestarts];

		SEXP names;

		if (!(r = (struct result*)R_alloc(irestarts,sizeof(*r))))
			goto bailout;

		have_margs = have_alphas = have_betas = have_ps = 0;

		/* TODO: The size is not really fixed */
		PROTECT(res = allocVector(VECSXP,6));
		PROTECT(names = allocVector(STRSXP,6));

		GetRNGstate();

		PROTECT(steps = AS_INTEGER(steps));
		PROTECT(burnin = AS_INTEGER(burnin));
		PROTECT(thin = AS_INTEGER(thin));
		PROTECT(flip_freq = AS_NUMERIC(flip_freq));
		params.nsteps = INTEGER_VALUE(steps);
		params.nsteps_burnin = INTEGER_VALUE(burnin);
		params.nsteps_thin = INTEGER_VALUE(thin);
		params.flip_freq = NUMERIC_VALUE(flip_freq);
		UNPROTECT(4);

#ifdef SUPPORT_OPENMP
		{
			int ithreads = INTEGER_VALUE(threads);
			if (ithreads != 0)
			{
				omp_set_num_threads(MIN(ithreads,omp_get_num_procs()));
			}
		}
#endif

#pragma omp parallel for
		for (run=0;run<irestarts;run++)
		{
			struct mt19937p mt;
			unsigned long seed;

#pragma omp critical
			{
				seed = (unif_rand() * (double)(INT32_MAX));
			}

			sgenrand(seed, &mt);

			/* Run! */
			r[run] = do_mgsa_mcmc(nsets, lset, lsets, INTEGER_VALUE(n),no,lo,
							&params,
							alpha_summary,beta_summary,p_summary, &mt);

			if (r[run].marg_set_activity) have_margs = 1;
			if (r[run].alpha_summary) have_alphas = 1;
			if (r[run].beta_summary) have_betas = 1;
			if (r[run].p_summary) have_ps = 1;
		}

		if (is_interrupted)
		{
			error("***Break");
			goto bailout;
		}
		PutRNGstate();

		/* Store set marginals */
		if (have_margs)
		{
			SEXP marg;
			PROTECT(marg = allocMatrix(REALSXP,lsets,irestarts));
			double *rmarg = REAL(marg);

			for (run=0;run<irestarts;run++)
				for (i=0;i<lsets;i++)
					rmarg[run * lsets + i] = r[run].marg_set_activity[i];

			SET_VECTOR_ELT(res,0,marg);
			SET_STRING_ELT(names,0,mkChar("marg"));

			/* Provide rownames, if any */
			if (getAttrib(sets,R_NamesSymbol) != R_NilValue)
			{
				SEXP dimnames;

				PROTECT(dimnames = allocVector(VECSXP,2));
				SET_VECTOR_ELT(dimnames, 0, getAttrib(sets,R_NamesSymbol));
				setAttrib(marg, R_DimNamesSymbol, dimnames);

				UNPROTECT(1);
			}
			UNPROTECT(1);
		}

		/* Store alpha marginals */
		if (have_alphas)
		{
			/* Build summaries array */
			for (run=0;run<irestarts;run++)
				summaries[run] = r[run].alpha_summary;

			SEXP alpha = create_R_representation_of_summary(summaries,irestarts);

			SET_VECTOR_ELT(res,1,alpha);
			SET_STRING_ELT(names,1,mkChar("alpha"));

			UNPROTECT(1);
		}

		/* Store beta marginals */
		if (have_betas)
		{
			/* Build summaries array */
			for (run=0;run<irestarts;run++)
				summaries[run] = r[run].beta_summary;

			SEXP beta = create_R_representation_of_summary(summaries,irestarts);

			SET_VECTOR_ELT(res,2,beta);
			SET_STRING_ELT(names,2,mkChar("beta"));

			UNPROTECT(1);
		}

		/* Store p marginals */
		if (have_ps)
		{
			/* Build summaries array */
			for (run=0;run<irestarts;run++)
				summaries[run] = r[run].p_summary;

			SEXP p = create_R_representation_of_summary(summaries,irestarts);

			SET_VECTOR_ELT(res,3,p);
			SET_STRING_ELT(names,3,mkChar("p"));

			UNPROTECT(1);
		}

		/* Store max stuff for each run (TODO: Add names) */
		{
			SEXP max;

			PROTECT(max = allocVector(VECSXP,irestarts));

			for (run=0;run<irestarts;run++)
			{
				SEXP max_run, max_run_names;
				SEXP el;

				int o;

				PROTECT(max_run = allocVector(VECSXP, 5));
				PROTECT(max_run_names = allocVector(STRSXP,5));

				PROTECT(el = allocVector(REALSXP,1));
				REAL(el)[0] = r[run].max_score;
				SET_VECTOR_ELT(max_run,0,el);
				SET_STRING_ELT(max_run_names,0,mkChar("score"));
				UNPROTECT(1);

				PROTECT(el = allocVector(REALSXP,1));
				REAL(el)[0] = r[run].max_score_alpha;
				SET_VECTOR_ELT(max_run,1,el);
				SET_STRING_ELT(max_run_names,1,mkChar("alpha"));
				UNPROTECT(1);

				PROTECT(el = allocVector(REALSXP,1));
				REAL(el)[0] = r[run].max_score_beta;
				SET_VECTOR_ELT(max_run,2,el);
				SET_STRING_ELT(max_run_names,2,mkChar("beta"));
				UNPROTECT(1);

				PROTECT(el = allocVector(REALSXP,1));
				REAL(el)[0] = r[run].max_score_p;
				SET_VECTOR_ELT(max_run,3,el);
				SET_STRING_ELT(max_run_names,3,mkChar("p"));
				UNPROTECT(1);

				PROTECT(el = allocVector(REALSXP,r[run].max_score_sets_active_length));
				for (o=0;o<r[run].max_score_sets_active_length;o++)
					REAL(el)[o] = r[run].max_score_sets_active[o] + 1; /* Don't forget to add one, as R starts from one */
				SET_VECTOR_ELT(max_run,4,el);
				SET_STRING_ELT(max_run_names,4,mkChar("sets"));
				UNPROTECT(1);

				setAttrib(max_run,R_NamesSymbol,max_run_names);
				SET_VECTOR_ELT(max,run,max_run);

				UNPROTECT(2); /* max_run and max_run_names */
			}

			SET_VECTOR_ELT(res,4,max);
			SET_STRING_ELT(names,4,mkChar("max"));

			UNPROTECT(1);
		}
		SEXP nsamples;
		PROTECT(nsamples = NEW_INTEGER(1));
		INTEGER_POINTER(nsamples)[0] = r[0].nsamples;
		SET_VECTOR_ELT(res,5,nsamples);
		SET_STRING_ELT(names,5,mkChar("nsamples"));
		UNPROTECT(1);

		setAttrib(res,R_NamesSymbol,names);
		UNPROTECT(2);
	}

	bailout:
	UNPROTECT(10);
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
		{"mgsa_mcmc", (DL_FUNC)&mgsa_mcmc, 14},
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
