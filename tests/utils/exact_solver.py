import sys
import argparse
import dimod
from dimod.serialization import coo

__copyright__ = "hyperQ – Ewa Hendzel"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Dev"

parser = argparse.ArgumentParser(description='Solve QUBO problem.')
parser.add_argument('--input', type=str, help='input qubo file')
parser.add_argument('--output', type=str, help='output csv file')
args = parser.parse_args()

if not args.input or not args.output:
    parser.print_help()
    exit(-1)

with open(args.input) as f:
    qubo = coo.load(f, vartype=dimod.BINARY)
    result = dimod.ExactSolver().sample(qubo)
    result.to_pandas_dataframe().sort_values("energy").to_csv(args.output)