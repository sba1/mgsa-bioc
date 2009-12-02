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
 *
 *   "R CMD INSTALL ../workspace/mgsa/ && (echo "library(mgsa);.Call(\"mgsa_mcmc\",list(c(2,1,3),c(1,2,3)),10,o=c(3,1),4,5,6)" | R --vanilla)"
 */
#include <stdio.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Enable debugging */
#define DEBUG

/* Define if file should be compiled as standalone (mainly for testing purposes) */
/* #define STANDALONE */

struct context
{
	/** @brief Number of sets */
	int number_of_sets;

	/** @brief Sizes of sets */
	int *sizes_of_sets;

	/** The definitions of the sets */
	int **sets;

	/** @brief Array of size "number_of_sets" indicating the state of the sets (active or not active) */
	int *sets_active;

	/** @brief Number of observable */
	int number_of_observable;

	/** @brief Array indicating whether a observable is active or not */
	int *observable;

	/** @brief Array, indicating how much sets lead to the activation of a hidden entry */
	int *hidden_count;

	/** @brief The number hidden elements that are active */
	int number_of_hidden_active;

	/**
	 * @brief Array that contains the hidden that are active.
	 *
	 * Only the first number_of_hidden_active elements are interesting.
	 **/
	int *hidden_active;

	int n00;
	int n01;
	int n10;
	int n11;
	double alpha;
	double beta;
	double p;

	/** @brief the current score */
	double score;
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

	if (!(cn->hidden_active = (int*)R_alloc(n,sizeof(cn->hidden_active[0]))))
		goto bailout;
	memset(cn->hidden_active,0,n * sizeof(cn->hidden_active[0]));
	cn->number_of_hidden_active = 0;

	if (!(cn->hidden_count = (int*)R_alloc(n,sizeof(cn->hidden_count[0]))))
		goto bailout;
	memset(cn->hidden_count,0,n * sizeof(cn->hidden_count[0]));

	if (!(cn->observable = (int*)R_alloc(n,sizeof(cn->observable[0]))))
		goto bailout;
	memset(cn->observable,0,n * sizeof(cn->observable[0]));
	for (i=0;i<lo;i++)
		cn->observable[o[i]] = 1;

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
		/* observation is true positive rather than false postive
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
			/* A not yet active member is about to be activated, append it at the hidden_active list */
//			cn->hidden_active[cn->number_of_hidden_active] = member;
			cn->number_of_hidden_active++;
			hidden_member_activated(cn,member);
		}
		cn->hidden_count[member]++;
	}

	/* TODO: Update the inactive list */
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
//			cn->hidden_active[cn->number_of_hidden_active] = member;
			cn->number_of_hidden_active--;
			hidden_member_deactivated(cn,member);
		}
		cn->hidden_count[member]--;
	}

	/* TODO: Update the active list */
}

/**
 * @brief Switches the state of the given set.
 *
 * @param cn defines the context on which to operate
 * @param to_switch defines the set which should be switched.
 */
static void switch_state(struct context *cn, int to_switch)
{
	int new_state;

	new_state = cn->sets_active[to_switch];

	if (new_state)
		add_set(cn,to_switch);
	else
		remove_set(cn,to_switch);
}

/**
 * The work horse.
 *
 * @param sets pointer to the sets. Sets a made of observable entities.
 * @param sizes_of_sets specifies the length of each set
 * @param number_of_sets number of sets (length of sets)
 * @param n the number of observable entities.
 * @param o indices of entities which are "on" (0 based).
 * @param lo length of o
 */
static void do_mgsa_mcmc(int **sets, int *sizes_of_sets, int number_of_sets, int n, int *o, int lo)
{
	int i,j;
	struct context cn;

#ifdef DEBUG

	for (i=0;i<number_of_sets;i++)
	{
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
#endif

	if (!init_context(&cn,sets,sizes_of_sets,number_of_sets,n,o,lo))
		goto bailout;


bailout:
	printf("huhu");
}

/**
 * @param n the number of observable entities.
 * @param o (1 based)
 *
 * @note TODO: Check whether sets are real sets.
 */
SEXP mgsa_mcmc(SEXP sets, SEXP n, SEXP o, SEXP alpha, SEXP beta, SEXP p)
{
	int *xo,*no,lo;
	int **nsets, *lset, lsets;
	int i,j;

	if (LENGTH(n) != 1)
		error("Parameter 'n' needs to be atomic!");

	PROTECT(n = AS_INTEGER(n));
	PROTECT(o = AS_INTEGER(o));
	PROTECT(sets = AS_LIST(sets));

	/* Observations */
	xo = INTEGER_POINTER(o);
	lo = LENGTH(o);
	if (!(no =  (int*)R_alloc(lo,sizeof(no[0])))) /* R's garbage collection takes care of it */
		goto bailout;

	/* Turn 1 based observation indices into 0 based ones */
	for (i=0;i<lo;i++)
		no[i] = xo[i] - 1;

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
			nsets[i][j] = xset[j];
		UNPROTECT(1);
	}

	do_mgsa_mcmc(nsets, lset, lsets, INTEGER_VALUE(n),no,lo);

bailout:
	UNPROTECT(3);
	return NULL_USER_OBJECT;
}

R_CallMethodDef callMethods[] = {
   {"mgsa_mcmc", (DL_FUNC)&mgsa_mcmc, 6},
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
