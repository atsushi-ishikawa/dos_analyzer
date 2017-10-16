import sys

argvs = sys.argv
from tools import json_to_csv
json_to_csv(argvs[1],argvs[2])

