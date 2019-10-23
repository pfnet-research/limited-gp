from cocoex import default_observers
from cocoex import Observer
from cocoex import Suite
from cocoex.utilities import ObserverOptions
from tqdm import tqdm
from typing import Callable  # NOQA
from typing import Optional  # NOQA


class Experiment(object):
    def __init__(self,
                 solver,
                 suite_name="bbob",
                 suite_instance="",
                 suite_options="dimensions: 2,3",
                 algorithm_name=None):
        # type: (Callable, str, str, str, Optional[str]) -> None

        self._solver = solver
        self._suite_name = suite_name
        self._suite_instance = suite_instance
        self._suite_options = suite_options

        default_algorithm_name = '{}({})'.format(solver.__name__, solver.__module__)
        self._algorithm_name = algorithm_name or default_algorithm_name

    def run(self, budget=100, current_batch=1, number_of_batches=15):
        # type: (int, int, int) -> None

        suite = Suite(self._suite_name, self._suite_instance,
                      self._suite_options)

        observer_name = default_observers()[self._suite_name]

        observer_options = self._build_observer_options(budget)
        observer = Observer(observer_name, observer_options.as_string)

        for problem_index, problem in enumerate(tqdm(suite)):
            if (problem_index % number_of_batches) != current_batch - 1:
                continue

            observer.observe(problem)

            max_evals = budget * problem.dimension
            self._solver(problem,
                         problem.lower_bounds,
                         problem.upper_bounds,
                         max_evals - problem.evaluations_constraints,
                         verbose=False)

    def _build_observer_options(self, budget):
        # type: (int) -> ObserverOptions

        opts = {
            'result_folder':
            '"%s/on_%s_budget%04dxD"' %
            (self._algorithm_name, self._suite_name, budget),
            'algorithm_name': self._algorithm_name
        }
        return ObserverOptions(opts)
