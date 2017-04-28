/**
 * @file mgsa.c
 *
 * @author Sebastian Bauer
 *
 * @note
 *
 *  If you want to test using R something like
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);sets<-list(c(1,2),c(3));names(sets)<-c(\"set 1\", \"set 2\");mgsa:::mgsa.trampoline(c(1,2),sets,3,steps=100,restarts=4)" | R --vanilla)"
 *  or (for invoking the function directly)
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_mcmc\",list(c(1,2),c(3)),3,o=c(1,2),4,5,6,steps=100)" | R --vanilla)"
 *  or (for the test function)
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_test\")" | R --vanilla)"
 */

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

#include "mgsa-core.h"
#include "mt.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

/** Global variable which is set to 1, if thread has been interrupted */
static int is_interrupted;

/*************************************************************/

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

	if (!(p = (struct parameter_prior*)R_alloc(1,sizeof(*p))))
		return NULL;

	PROTECT(sexp = AS_NUMERIC(sexp));
	p->number_of_states = LENGTH(sexp);

	if (discrete)
	{
		p->uniform_continuous = 0;

		if (p->number_of_states == 0)
			error("Parameter '%s' has been requested to be discrete but no values were specified");

		if (!(p->values = (double*)R_alloc(p->number_of_states,sizeof(p->values[0]))))
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

		if (!(sum->dmap = (int*)R_alloc(number_of_discrete_values,sizeof(sum->dmap[0]))))
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
		nas = (int*)R_alloc(las,sizeof(nas[0]));

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
							alpha_summary,beta_summary,p_summary, &mt, &is_interrupted);

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
