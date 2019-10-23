import argparse
from experiment import Experiment
import Limited_GP

parser = argparse.ArgumentParser()
parser.add_argument('--number-of-batches', type=int, default=1)
parser.add_argument('--current-batch', type=int, default=1)
parser.add_argument('--budget', type=int, default=10000)
parser.add_argument('--suite-name', default='bbob')
args = parser.parse_args()

e = Experiment(suite_name=args.suite_name, solver=Limited_GP.solve, algorithm_name='Limited-GP')
e.run(budget=args.budget, current_batch=args.current_batch, number_of_batches=args.number_of_batches)
