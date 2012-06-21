"""
Base classes Model and Sampler are defined here.
"""
from numpy import zeros, floor
import sys,os
from copy import copy
from threading import Thread
import thread
from time import sleep
import pdb
import warnings, exceptions, traceback

GuiInterrupt = 'Computation halt'
Paused = 'Computation paused'

class EndofSampling(Exception):
    pass


class Model():
    """
    The base class for all objects that fit probability models. Model is initialized with:

      >>> A = Model(input, verbose=0)

    """
    def __init__(self, input=None, name=None, verbose=0):
        """Initialize a Model instance.

        :Parameters:
          - input : module, list, tuple, dictionary, set, object or nothing.
              Model definition, in terms of Stochastics, Deterministics, Potentials and Containers.
              If nothing, all nodes are collected from the base namespace.
        """

        # Get stochastics, deterministics, etc.
        if input is None:
            import warnings
            warnings.warn('The MCMC() syntax is deprecated. Please pass in nodes explicitly via M = MCMC(input).')
            import __main__
            __main__.__dict__.update(self.__class__.__dict__)
            input = __main__

        ObjectContainer.__init__(self, input)

        if name is not None:
            self.__name__ = name
        self.verbose = verbose

    def _get_generations(self):
        if not hasattr(self, '_generations'):
            self._generations = utils.find_generations(self)
        return self._generations
    generations = property(_get_generations)

    def draw_from_prior(self):
        """
        Sets all variables to random values drawn from joint 'prior', meaning contributions
        of data and potentials to the joint distribution are not considered.
        """

        for generation in self.generations:
            for s in generation:
                s.random()

    def seed(self):
        """
        Seed new initial values for the stochastics.
        """

        for generation in self.generations:
            for s in generation:
                try:
                    if s.rseed is not None:
                        value = s.random(**s.parents.value)
                except:
                    pass




class Sampler(Model):
    """
    The base class for all objects that fit probability models using Monte Carlo methods.
    Sampler is initialized with:

      >>> A = Sampler(input, db, output_path=None, verbose=0)


    """
    def __init__(self, input=None, db='ram', name='Sampler', reinit_model=True, calc_deviance=False, verbose=0, **kwds):
        """Initialize a Sampler instance.


        """


        # Initialize superclass
        if reinit_model:
            Model.__init__(self, input, name, verbose)

        # Initialize deviance, if asked
        if calc_deviance:
            self._funs_to_tally = {'deviance': self._sum_deviance}
        else:
            self._funs_to_tally = {}

        # Specify database backend and save its keywords
        self._db_args = kwds
        self._assign_database_backend(db)

        # Flag for model state
        self.status = 'ready'

        self._current_iter = None
        self._iter = None

        self._state = ['status', '_current_iter', '_iter']

    def _sum_deviance(self):
        # Sum deviance from all stochastics

        return -2*sum([v.get_logp() for v in self.observed_stochastics])

    def sample(self, iter, length=None, verbose=0):
        """
        Draws iter samples from the posterior.
        """
        self._cur_trace_index=0
        self.max_trace_length = iter
        self._iter = iter

        if verbose>0:
            self.verbose = verbose
        self.seed()

        # Initialize database -> initialize traces.
        if length is None:
            length = iter
        self.db._initialize(self._funs_to_tally, length)

        # Put traces on objects
        for v in self._variables_to_tally:
            v.trace = self.db._traces[v.__name__]

        # Loop
        self._current_iter = 0
        self._loop()
        self._finalize()

    def _finalize(self):
        """Reset the status and tell the database to finalize the traces."""
        if self.status in ['running', 'halt']:
            if self.verbose > 0:
                print 'Sampling finished normally.'
            self.status = 'ready'

        self.save_state()
        self.db._finalize()

    def _loop(self):
        """

        """
        self.status='running'

        try:
            while self._current_iter < self._iter and not self.status == 'halt':
                if self.status == 'paused':
                    break

                i = self._current_iter

                self.draw()
                self.tally()

                if not i % 10000 and self.verbose > 0:
                    print 'Iteration ', i, ' of ', self._iter
                    sys.stdout.flush()

                self._current_iter += 1

        except KeyboardInterrupt:
            self.status='halt'

        if self.status == 'halt':
            self._halt()

    def draw(self):
        """
        Either draw() or _loop() must be overridden in subclasses of Sampler.
        """
        pass

    def stats(self, alpha=0.05, start=0, variable=None):
        """
        Statistical output for variables.
        """
        if variable is not None:
            return variable.stats(alpha=alpha, start=start)

        else:
            stat_dict = {}

            # Loop over nodes
            for variable in self._variables_to_tally:
                # Plot object
                stat_dict[variable.__name__] = variable.stats(alpha=alpha, start=start)

            return stat_dict

    # Property --- status : the sampler state.
    def status():
        doc = \
        """Status of sampler. May be one of running, paused, halt or ready.
        """
        def fget(self):
            return self.__status
        def fset(self, value):
            if value in ['running', 'paused', 'halt', 'ready']:
                self.__status=value
            else:
                raise AttributeError, value
        return locals()
    status = property(**status())


    def _assign_database_backend(self, db):
        """Assign Trace instance to stochastics and deterministics and Database instance
        to self.

        :Parameters:
          - `db` : string, Database instance

        """
        # Objects that are not to be tallied are assigned a no_trace.Trace
        # Tallyable objects are listed in the _nodes_to_tally set.

        no_trace = getattr(database, 'no_trace')
        self._variables_to_tally = set()
        for object in self.stochastics | self.deterministics :

            if object.trace:
                self._variables_to_tally.add(object)
                self._funs_to_tally[object.__name__] = object.get_value
            else:
                object.trace = no_trace.Trace(object.__name__)

        check_valid_object_name(self._variables_to_tally)

        # If not already done, load the trace backend from the database
        # module, and assign a database instance to Model.
        if type(db) is str:
            if db in dir(database):
                module = getattr(database, db)

                # Assign a default name for the database output file.
                if self._db_args.get('dbname') is None:
                    self._db_args['dbname'] = self.__name__ + '.' + db

                self.db = module.Database(**self._db_args)
            elif db in database.__modules__:
                raise ImportError, \
                    'Database backend `%s` is not properly installed. Please see the documentation for instructions.' % db
            else:
                raise AttributeError, \
                    'Database backend `%s` is not defined in pymc.database.'%db
        elif isinstance(db, database.base.Database):
            self.db = db
            self.restore_sampler_state()
        else:   # What is this for? DH. If it's a user defined backend, it doesn't initialize a Database.
            self.db = db.Database(**self._db_args)

        # Assign Trace instances to tallyable objects.
        self.db.connect_model(self)


    def pause(self):
        """Pause the sampler. Sampling can be resumed by calling `icontinue`.
        """
        self.status = 'paused'
        # The _loop method will react to 'paused' status and stop looping.
        if hasattr(self, '_sampling_thread') and self._sampling_thread.isAlive():
            print 'Waiting for current iteration to finish...'
            while self._sampling_thread.isAlive():
                sleep(.1)

    def halt(self):
        """Halt a sampling running in another thread."""
        self.status = 'halt'
        # The _halt method is called by _loop.
        if hasattr(self, '_sampling_thread') and self._sampling_thread.isAlive():
            print 'Waiting for current iteration to finish...'
            while self._sampling_thread.isAlive():
                sleep(.1)

    def _halt(self):
        print 'Halting at iteration ', self._current_iter, ' of ', self._iter
        self.db.truncate(self._cur_trace_index)
        self._finalize()

    #
    # Tally
    #
    def tally(self):
       """
       tally()

       Records the value of all tracing variables.
       """
       if self.verbose > 2:
           print self.__name__ + ' tallying.'
       if self._cur_trace_index < self.max_trace_length:
           self.db.tally()

       self._cur_trace_index += 1
       if self.verbose > 2:
           print self.__name__ + ' done tallying.'

    def commit(self):
        """
        Tell backend database to commit.
        """

        self.db.commit()


    def isample(self, *args, **kwds):
        """
        Samples in interactive mode. Main thread of control stays in this function.
        """
        self._exc_info = None
        out = kwds.pop('out',  sys.stdout)
        def samp_targ(*args, **kwds):
            try:
                self.sample(*args, **kwds)
            except:
                self._exc_info = sys.exc_info()

        self._sampling_thread = Thread(target=samp_targ, args=args, kwargs=kwds)
        self.status = 'running'
        self._sampling_thread.start()
        self.iprompt(out=out)

    def icontinue(self):
        """
        Restarts thread in interactive mode
        """
        if self.status != 'paused':
            print "No sampling to continue. Please initiate sampling with isample."
            return

        def sample_and_finalize():
            self._loop()
            self._finalize()

        self._sampling_thread = Thread(target=sample_and_finalize)
        self.status = 'running'
        self._sampling_thread.start()
        self.iprompt()

    def iprompt(self, out=sys.stdout):
        """Start a prompt listening to user input."""

        cmds = """
        Commands:
          i -- index: print current iteration index
          p -- pause: interrupt sampling and return to the main console.
                      Sampling can be resumed later with icontinue().
          h -- halt:  stop sampling and truncate trace. Sampling cannot be
                      resumed for this chain.
          b -- bg:    return to the main console. The sampling will still
                      run in a background thread. There is a possibility of
                      malfunction if you interfere with the Sampler's
                      state or the database during sampling. Use this at your
                      own risk.
        """

        print >> out, """==============
 
==============

        
        """
        print >> out, cmds

        prompt = True
        try:
            while self.status in ['running', 'paused']:
                    # sys.stdout.write('pymc> ')
                    if prompt:
                        out.write('pymc > ')
                        out.flush()

                    if self._exc_info is not None:
                        a,b,c = self._exc_info
                        raise a, b, c

                    cmd = utils.getInput().strip()
                    if cmd == 'i':
                        print >> out,  'Current iteration: ', self._current_iter
                        prompt = True
                    elif cmd == 'p':
                        self.status = 'paused'
                        break
                    elif cmd == 'h':
                        self.status = 'halt'
                        break
                    elif cmd == 'b':
                        return
                    elif cmd == '\n':
                        prompt = True
                        pass
                    elif cmd == '':
                        prompt = False
                    else:
                        print >> out, 'Unknown command: ', cmd
                        print >> out, cmds
                        prompt = True

        except KeyboardInterrupt:
            if not self.status == 'ready':
                self.status = 'halt'


        if self.status == 'ready':
            print >> out, "Sampling terminated successfully."
        else:
            print >> out, 'Waiting for current iteration to finish...'
            while self._sampling_thread.isAlive():
                sleep(.1)
            print >> out, 'Exiting interactive prompt...'
            if self.status == 'paused':
                print >> out, 'Call icontinue method to continue, or call halt method to truncate traces and stop.'


    def get_state(self):
        """
        Return the sampler's current state in order to
        restart sampling at a later time.
        """
        state = dict(sampler={}, stochastics={})
        # The state of the sampler itself.
        for s in self._state:
            state['sampler'][s] = getattr(self, s)

        # The state of each stochastic parameter
        for s in self.stochastics:
            state['stochastics'][s.__name__] = s.value
        return state

    def save_state(self):
        """
        Tell the database to save the current state of the sampler.
        """
        try:
            self.db.savestate(self.get_state())
        except:
            print 'Warning, unable to save state.'
            print 'Error message:'
            traceback.print_exc()

    def restore_sampler_state(self):
        """
        Restore the state of the sampler and to
        the state stored in the database.
        """

        state = self.db.getstate() or {}

        # Restore sampler's state
        sampler_state = state.get('sampler', {})
        self.__dict__.update(sampler_state)

        # Restore stochastic parameters state
        stoch_state = state.get('stochastics', {})
        for sm in self.stochastics:
            try:
                sm.value = stoch_state[sm.__name__]
            except:
                warnings.warn(\
    'Failed to restore state of stochastic %s from %s backend'%(sm.__name__, self.db.__name__), exceptions.UserWarning)
                #print 'Error message:'
                #traceback.print_exc()


    def remember(self, chain=-1, trace_index = None):
        """
        remember(trace_index = randint(trace length to date))

        Sets the value of all tracing variables to a value recorded in
        their traces.
        """
        if trace_index is None:
            trace_index = randint(self.cur_trace_index)

        for variable in self._variables_to_tally:
            if isinstance(variable, Stochastic):
                variable.value = self.trace(variable.__name__, chain=chain)[trace_index]

    def trace(self, name, chain=-1):
        """Return the trace of a tallyable object stored in the database.

        :Parameters:
        name : string
          The name of the tallyable object.
        chain : int
          The trace index. Setting `chain=i` will return the trace created by
          the ith call to `sample`.
        """
        if type(name) is str:
            return self.db.trace(name, chain)
        elif isinstance(name, Variable):
            return self.db.trace(name.__name__, chain)
        else:
            raise ValueError, 'Name argument must be string or Variable, got %s.'%name

    def _get_deviance(self):
        return self._sum_deviance()
    deviance = property(_get_deviance)

def check_valid_object_name(sequence):
    """Check that the names of the objects are all different."""
    names = []
    for o in sequence:
        if o.__name__ in names:
            raise ValueError, 'A tallyable PyMC object called %s already exists. This will cause problems for some database backends.'%o.__name__
        else:
            names.append(o.__name__)
