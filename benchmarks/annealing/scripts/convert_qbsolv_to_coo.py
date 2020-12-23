
__copyright__ = "hyperQ – Ewa Hendzel"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Dev"

"""Convert QUBO instance written in Qbsolv format to COO format."""
from dimod import BQM
from argparse import ArgumentParser, FileType


if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument("input_file", help="File to be converted.", type=FileType("r"))
    parser.add_argument("output_file", help="Path to the output file", type=FileType("w"))

    args = parser.parse_args()

    print(args.input_file)
    print(args.output_file)

    ising = BQM.from_coo(args.input_file, vartype="SPIN")

    qubo_bqm = ising.change_vartype("BINARY", inplace=False)
    qubo_bqm.relabel_variables({i: i - 1 for i in qubo_bqm.variables}, inplace=True)

    num_variables = len(qubo_bqm.variables)
    num_couplers = len(qubo_bqm.quadratic)

    header = f"p qubo 0 {num_variables} {num_variables} {num_couplers}\n"
    qubo_dict, _ = qubo_bqm.to_qubo()

    args.output_file.write(header)

    for (i, j), coef in sorted(qubo_dict.items()):
        args.output_file.write(f"{i} {j} {coef}\n")

    args.output_file.close()


