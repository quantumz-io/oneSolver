
__copyright__ = "hyperQ – Ewa Hendzel"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Dev"

PATTERN = (
    r"tmp/benchmark_(?P<size>\d+)_(?P<size>\d+)_"
    r"betamin_(?P<beta_min>.+)_num_iter_(?P<num_iter>.+)_"
    r"num_tries_(?P<num_tries>\d+)_schedule_(?P<schedule>[a-z]+).csv"
)
if __name__ == '__main__':
