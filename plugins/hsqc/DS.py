"""
Evidencial step
"""
import numpy as np
from copy import copy

def open01(x, limit=1.e-6):
    """Constrain numbers to (0,1) interval"""
    try:
        return np.array([min(max(y, limit), 1.-limit) for y in x])
    except TypeError:
        return min(max(x, limit), 1.-limit)

def diagnostic(f):
    """
    This decorator allows  various types to be passed to
    the diagnostic functions. 

    """

    def wrapper(pymc_obj, *args, **kwargs):

        # Figure out what type of object it is
        try:
            values = {}
            # First try Model type
            for variable in pymc_obj._variables_to_tally:
                data = variable.trace()
                name = variable.__name__
                if kwargs.get('verbose'):
                    print "\nDiagnostic for %s ..." % name
                values[name] = f(data, *args, **kwargs)
            return values
        except AttributeError:
            pass

        try:
            # Then try Node type
            data = pymc_obj.trace()
            name = pymc_obj.__name__
            return f(data, *args, **kwargs)
        except (AttributeError,ValueError):
            pass

        # If others fail, assume that raw data is passed
        return f(pymc_obj, *args, **kwargs)

    return wrapper


def validate(sampler, replicates=20, iterations=10000, burn=5000, thin=1, deterministic=False, db='ram', plot=True, verbose=0):
    """
    

    Generates posterior samples based on 'true' parameter values and data simulated
    from the priors. The quantiles of the parameter values are calculated, based on
    the samples. If the model is valid, the quantiles should be uniformly distributed
    over [0,1].

 
    """
    import scipy as sp
    # Set verbosity for models to zero
    sampler.verbose = 0

    # Specify parameters to be evaluated
    parameters = sampler.stochastics
    if deterministic:
        # Add deterministics to the mix, if requested
        parameters = parameters.union(sampler.deterministics)

    # Assign database backend
    original_backend = sampler.db.__name__
    sampler._assign_database_backend(db)

    # Empty lists for quantiles
    quantiles = {}

    if verbose:
        print "\nExecuting Cook et al. (2006) validation procedure ...\n"

    # Loop over replicates
    for i in range(replicates):

        # Sample from priors
        for p in sampler.stochastics:
            if not p.extended_parents:
                p.random()

        # Sample "true" data values
        for o in sampler.observed_stochastics:
            # Generate simuated data for data stochastic
            o.set_value(o.random(), force=True)
            if verbose:
                print "Data for %s is %s" % (o.__name__, o.value)

        param_values = {}
        # Record data-generating parameter values
        for s in parameters:
            param_values[s] = s.value

        try:
            # Fit models given parameter values
            sampler.sample(iterations, burn=burn, thin=thin)

            for s in param_values:

                if not i:
                    # Initialize dict
                    quantiles[s.__name__] = []
                trace = s.trace()
                q = sum(trace<param_values[s], 0)/float(len(trace))
                quantiles[s.__name__].append(open01(q))

            # Replace data values
            for o in sampler.observed_stochastics:
                o.revert()

        finally:
            # Replace data values
            for o in sampler.observed_stochastics:
                o.revert()

            # Replace backend
            sampler._assign_database_backend(original_backend)

        if not i % 10 and i and verbose:
            print "\tCompleted validation replicate", i


    # Replace backend
    sampler._assign_database_backend(original_backend)

    stats = {}
    # Calculate chi-square statistics
    for param in quantiles:
        q = quantiles[param]
        # Calculate chi-square statistics
        X2 = sum(sp.special.ndtri(q)**2)
        # Calculate p-value
        p = sp.special.chdtrc(replicates, X2)

        stats[param] = (X2, p)

    if plot:
        # Convert p-values to z-scores
        p = copy(stats)
        for i in p:
            p[i] = p[i][1]
        pymc.Matplot.zplot(p, verbose=verbose)

    return stats


@diagnostic
def geweke(x, first=.1, last=.5, intervals=20):
    """Return z-scores for convergence diagnostics.

 



    """
    # Filter out invalid intervals
    if first + last >= 1:
        raise "Invalid intervals for Geweke convergence analysis",(first,last)

    # Initialize list of z-scores
    zscores = []

    # Last index value
    end = len(x) - 1
