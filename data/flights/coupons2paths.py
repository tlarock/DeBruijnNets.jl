import pathpy as pp
import pandas as pd
import argparse
from random import random
import progressbar

# turns dictionary with flight segments into tuple
def tickets_to_tuple(path):
    if 1 in path:
        t = (path[1][0],)
    else: 
        print('Skipped inconsistent itinerary')
        return ()

    for i in range(2, len(path)+1):
        if i in path and path[i][0] == path[i-1][1]:
            t += (path[i][0],)
        else:
            print('Skipped inconsistent itinerary')
            return ()
    t += (path[len(path)][1],)
    return t

parser = argparse.ArgumentParser(description='Extracts flight paths from DB1DB Coupon data files.')
parser.add_argument('coupon_data', help='path to coupon csv data file.', type=str)
parser.add_argument('outfile', help='path to ngram file storing results.', type=str)
parser.add_argument('--fraction', help='fraction of itineraries to take.', type=float)
args = parser.parse_args()

fraction = 1.0
if args.fraction:
    fraction = args.fraction


print('Reading coupon data ...')
df = pd.read_csv(args.coupon_data)
itins = df.groupby('ITIN_ID')
print('done.')

# Paths object
p = pp.Paths()

# collect a `fraction` of all paths
est_num = int(1.03*fraction*itins.ngroups)
with progressbar.ProgressBar(max_value=est_num) as bar:
    i=0
    for itname, seq in itins:
        # only process "p_to_take" percent of itineraries
        if random() > fraction:
            continue

        # dictionary indexed by SeqNum, supports out-of-order flight segments
        current_itin = {}    
        for _,leg in seq.iterrows():
            current_itin[leg['SEQ_NUM']] = (leg['ORIGIN'], leg['DEST'])

        if current_itin:
            x = tickets_to_tuple(current_itin)
            if x:
                p.add_path(x)

        if i < est_num:
            bar.update(i)
        i += 1

print('Writing target file ...')
p.write_file(args.outfile)
print('done.')
