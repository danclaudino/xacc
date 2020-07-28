import xacc
from pelix.ipopo.decorators import ComponentFactory, Property, Requires, Provides, \
    Validate, Invalidate, Instantiate

@ComponentFactory("scipy_optimizer_factory")
@Provides("optimizer")
@Property("_optimizer", "optimizer", "scipy")
@Property("_name", "name", "scipy")
@Instantiate("scipy_instance")
class scipyOptimizer(xacc.Optimizer):
    def __init__(self):
        xacc.Optimizer.__init__(self)
        self.options = {}
        self.sigma = .1
        self.hetOpts = None

    def name(self):
        return 'scipy'

    def isGradientBased(self):
        if self.hetOpts['scipy-optimizer'] == 'BFGS' or self.hetOpts['scipy-optimizer'] == 'L-BFGS-B' or self.hetOpts['scipy-optimizer'] == 'CG' or self.hetOpts['scipy-optimizer'] == 'Newton-CG':
            return True
        else:
            return False

    def get_algorithm(self):
        return self.hetOpts['scipy-optimizer']

    def setOptions(self, opts):
        self.hetOpts = opts
        CURRENT_SCIPY_OPT_ALGO = ['Nelder-Mead', 'Powell', 'BFGS', 'L-BFGS-B', 'CG', 'Newton-CG']

        if 'scipy-optimizer' in opts and opts['scipy-optimizer'] in CURRENT_SCIPY_OPT_ALGO:
            self.scipy_opt_name = opts['scipy-optimizer']
        else:
            xacc.error("Please provide an optimization algorithm among: " % CURRENT_SCIPY_OPT_ALGO)

    def optimize(self, function):
        import scipy.optimize
        params = [] 
        if 'initial-parameters' in self.hetOpts:
            params = self.hetOpts['initial-parameters']
        else:
            params = list(function.dimensions() * [0.])
        r = scipy.optimize.minimize(function, params, method = self.scipy_opt_name)

        return (r['fun'], r['x'])