#!/usr/bin/env python

import os
import shutil
import sys
import subprocess
import argparse
import time
import random
from multiprocessing import Pool

# import math, datetime
# import rpy2.robjects.lib.ggplot2 as ggplot2
# from rpy2.robjects import *
# from rpy2.robjects.packages import importr
# base = importr('base')
# grdevices = importr('grDevices')

parser = argparse.ArgumentParser()
parser.add_argument("bp_file", help="BED file of breakpoint regions")
parser.add_argument("bp_name", help="Name of this set of breakpoint regions")
parser.add_argument("feature_file", help="BED file of features to test against")
parser.add_argument("feature_name", help="name of the features being tested against")
parser.add_argument("genome", help="FAI file of the genome")
parser.add_argument("--gaps", help="BED file of gaps in the genome")
parser.add_argument("--process_iteration_chunk", help="one of bp_hits|feature_hits|bases_overlapped")
parser.add_argument("--iteration_number", type=int, help="starting iteration number")
parser.add_argument("--iteration_chunk_size", type=int, default=1, help="number of iterations in this chunk")
parser.add_argument("--permutations", type=int, default=10000, help="iteration number")
parser.add_argument("--cores", type=int, default=100, help="number of worker cores on the cluster to use")
parser.add_argument("--replot", action='store_true', dest="replot", help="only replot a directory, don't redo the permutations")
parser.add_argument("--shift_only", action='store_true', dest="shift_only", help="only do shift analysis")

args = parser.parse_args()

bp_file = args.bp_file
bp_name = args.bp_name
feature_file = args.feature_file
feature_name = args.feature_name
genome = args.genome
gaps = args.gaps
process_iteration_cmd = args.process_iteration_chunk
iteration_number = args.iteration_number
iteration_chunk_size = args.iteration_chunk_size
permutations = args.permutations
cores = args.cores
replot = args.replot
shift_only = args.shift_only

script_path = os.path.dirname( os.path.realpath( __file__ ) )

def process_iteration(i, fun, gaps, bp_file, genome, feature_file):
    if gaps != None:
        shufflep = subprocess.Popen(['shuffleBed', '-excl', gaps, '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    else:
        shufflep = subprocess.Popen(['shuffleBed', '-i', bp_file, '-g', genome, '-chrom'], stdout=subprocess.PIPE)
    sortp = subprocess.Popen('sort -k1,1 -k2,2n', stdin=shufflep.stdout, stdout=subprocess.PIPE, shell=True)
    shuffled = sortp.communicate()[0]
    hits = fun(feature_file, '-', input=shuffled)
    return (str(i) + "\t" + str(hits))

# doesn't support passing gaps yet
def wrap_process_iteration(arg_tuple):
    myargs = arg_tuple
    time.sleep(random.random() * 60)
    iteration_command = ['condor_run', 'python', script_path + '/overlap_bps_with_feature.py', myargs[3], "none", myargs[5], "none", myargs[4], "--process_iteration", myargs[1], "--iteration_number", str(myargs[0]), "--iteration_chunk_size", str(myargs[6])]
    result = subprocess.Popen(iteration_command, stdout=subprocess.PIPE).communicate()[0]
    return result

def compute_bases_overlapped(feature_file, bp_file, input=None):
    if bp_file == '-':
#        random_bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wo', '-sorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
        random_bp_hits = subprocess.Popen(['bedmap', '--bases-uniq', bp_file, feature_file], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
    else:
        random_bp_hits = subprocess.Popen(['bedmap', '--bases-uniq', bp_file, feature_file], stdout=subprocess.PIPE).communicate()[0]
    bases_overlapped = 0
    for line in random_bp_hits.split("\n"):        
        if line != "":
            bases_overlapped += int(line)
    return bases_overlapped

# def compute_bp_hits(feature_file, bp_file, input=None):
#     if bp_file == 'stdin':
#         bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wb', '-sorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
#     else:
#         bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wb', '-sorted'], stdout=subprocess.PIPE).communicate()[0]
#     hits = 0
#     current_bp = ""
#     bps = set()
#     for line in bp_hits.split("\n"):
#         fields = line.split("\t")
#         if len(fields) > 6:
#             sys.exit("please make sure your feature file is in BED3 format")
#         if len(fields) == 6:
#             bp = "\t".join(fields[3:6])
#             bps.add(bp)
#     return len(bps)

def compute_bp_hits(feature_file, bp_file, input=None):
    if bp_file == '-':
#        bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wa', '-u', '-sorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
        bp_hits = subprocess.Popen(['bedops', '--element-of', '-1',  bp_file,  feature_file], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
    else:
#        bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wa', '-u', '-sorted'], stdout=subprocess.PIPE).communicate()[0]
        bp_hits = subprocess.Popen(['bedops', '--element-of', '-1', bp_file,  feature_file], stdout=subprocess.PIPE).communicate()[0]
    hits = len(bp_hits.split("\n")) - 1
    return hits


def compute_feature_hits(feature_file, bp_file, input=None):
    if bp_file == '-':
#        bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wa', '-u', '-sorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
        bp_hits = subprocess.Popen(['bedops', '--element-of', '-1',  feature_file,  bp_file], stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate(input)[0]
    else:
#        bp_hits = subprocess.Popen(['intersectBed', '-a', feature_file, '-b', bp_file, '-wa', '-u', '-sorted'], stdout=subprocess.PIPE).communicate()[0]
        bp_hits = subprocess.Popen(['bedops', '--element-of', '-1', feature_file,  bp_file], stdout=subprocess.PIPE).communicate()[0]
    hits = len(bp_hits.split("\n")) - 1
    return hits

if process_iteration_cmd == None:
    if not replot and not shift_only:
        iterations = permutations
        p = Pool(cores)
    # compute the real number of bps that overlap a feature
    num_real_bp_hits = compute_bp_hits(feature_file, bp_file)
    print "num bps that overlap features: {0}".format(num_real_bp_hits)

    # compute the real number of features that overlap a bp region
    num_real_feature_hits = compute_feature_hits(feature_file, bp_file)
    print "num features that overlap bps: {0}".format(num_real_feature_hits)

    # compute the real number of bases overlapped
    real_bases_overlapped = compute_bases_overlapped(feature_file, bp_file)
    print "num bases in bps that overlap features: {0}".format(real_bases_overlapped)

    if not replot:
        if not shift_only:
            # use shuffleBed to permute bp locations
            print "Computing permutations for BP hits"
            bps_with_hits_file = open('bps_with_hits.txt', 'w')
            num_random_bp_hits = p.map(wrap_process_iteration, [(i, "bp_hits", gaps, bp_file, genome, feature_file, iteration_chunk_size) for i in xrange(0, iterations, iteration_chunk_size)])
            for ip in ("".join(num_random_bp_hits)).rstrip().split("\n"):
                fields = ip.split("\t")
                bps_with_hits_file.write(str(fields[0]) + "\t" + str(fields[1]) + "\n")
            bps_with_hits_file.close()                             

            print "Computing permutations for feature hits"
            feature_hits_file = open('feature_hits.txt', 'w')
            num_random_feature_hits = p.map(wrap_process_iteration, [(i, "feature_hits", gaps, bp_file, genome, feature_file, iteration_chunk_size) for i in xrange(0, iterations, iteration_chunk_size)])
            for ip in ("".join(num_random_feature_hits)).rstrip().split("\n"):
                fields = ip.split("\t")
                feature_hits_file.write(str(fields[0]) + "\t" + str(fields[1]) + "\n")
            feature_hits_file.close()                            

            print "Computing permutations for bases overlapped"
            bases_overlap_file = open('bases_overlapped.txt', 'w')
            num_bases_overlapped = p.map(wrap_process_iteration, [(i, "bases_overlapped", gaps, bp_file, genome, feature_file, iteration_chunk_size) for i in xrange(0, iterations, iteration_chunk_size)])
            for ip in ("".join(num_bases_overlapped)).rstrip().split("\n"):
                fields = ip.split("\t")
                bases_overlap_file.write(str(fields[0]) + "\t" + str(fields[1]) + "\n")
            bases_overlap_file.close()                            

        print "Computing shifts"
        # shift the regions over and compute overlaps
        shifted_bp_hits_file = open('shifted_bp_hits.txt', 'w')
        shifted_feature_hits_file = open('shifted_feature_hits.txt', 'w')
        shifted_bases_overlapped_file = open('shifted_bases_overlapped.txt', 'w')
        for i in xrange(-1000000,1000000,25000):
            #sys.stderr.write("shift: " + str(i) + "\n")
            #sys.stderr.write("\t".join(['slopBed', '-i', bp_file, '-g', genome, '-l', str(i), '-r', str(-1 * i)]))
            #sys.stderr.write("\n")
            shifted_regions = subprocess.Popen(['slopBed', '-i', bp_file, '-g', genome, '-l', str(i), '-r', str(-1 * i)], stdout=subprocess.PIPE).communicate()[0]

            # regions very close to the ends of chromosomes can end up with starts greater than ends; need to filter them out
            filtered_shifted_regions = []
            for line in shifted_regions.split("\n"):
                if line == "":
                    continue
                fields = line.split("\t")
                start = int(fields[1])
                end = int(fields[2])
                if start < end:
                    filtered_shifted_regions.append(line)
                else:
                    sys.stderr.write("filtering out shifted region: " + line + "\n")
            
            #sys.stderr.write("filtered_shifted_regions: " + "\n".join(filtered_shifted_regions))
            #sys.stderr.write("left with {0} shifted regions\n".format(len(filtered_shifted_regions)))
            bp_hits = float(compute_bp_hits(feature_file, '-', "\n".join(filtered_shifted_regions))) / len(filtered_shifted_regions)
            shifted_bp_hits_file.write(str(i) + "\t" + str(bp_hits) + "\n")
            feature_hits = float(compute_feature_hits(feature_file, '-', "\n".join(filtered_shifted_regions))) / len(filtered_shifted_regions)
            shifted_feature_hits_file.write(str(i) + "\t" + str(feature_hits) + "\n")
            bases_overlapped = float(compute_bases_overlapped(feature_file, '-', "\n".join(filtered_shifted_regions))) / len(filtered_shifted_regions)
            shifted_bases_overlapped_file.write(str(i) + "\t" + str(bases_overlapped) + "\n")
        shifted_bp_hits_file.close()
        shifted_feature_hits_file.close()
        shifted_bases_overlapped_file.close()

    if not shift_only:
        results_dirname = bp_name.replace(' ', '_') + "_to_" + feature_name.replace(' ', '_')
        # plot the results
        if replot:
            subprocess.call(map(str, ['/g/whelanch/software/bin/Rscript', script_path + '/plot_results.R', feature_name, bp_name, num_real_bp_hits, num_real_feature_hits, real_bases_overlapped, results_dirname]))
        else:
            subprocess.call(map(str, ['/g/whelanch/software/bin/Rscript', script_path + '/plot_results.R', feature_name, bp_name, num_real_bp_hits, num_real_feature_hits, real_bases_overlapped, '.']))

            os.mkdir(results_dirname)
            shutil.move('bases_overlapped.txt', results_dirname + "/")
            shutil.move('bps_with_hits.txt', results_dirname + "/")
            shutil.move('feature_hits.txt', results_dirname + "/")
            shutil.move('shifted_bp_hits.txt', results_dirname + "/")
            shutil.move('shifted_feature_hits.txt', results_dirname + "/")
            shutil.move('shifted_bases_overlapped.txt', results_dirname + "/")
            shutil.move('results.txt', results_dirname + "/")
            shutil.move(results_dirname + ".pdf", results_dirname + "/")

elif process_iteration_cmd == "feature_hits":
    for i in xrange(iteration_number, iteration_number + iteration_chunk_size):
        print process_iteration(i, compute_feature_hits, gaps, bp_file, genome, feature_file)
elif process_iteration_cmd == "bp_hits":
    for i in xrange(iteration_number, iteration_number + iteration_chunk_size):
        print process_iteration(i, compute_bp_hits, gaps, bp_file, genome, feature_file)
elif process_iteration_cmd == "bases_overlapped":
    for i in xrange(iteration_number, iteration_number + iteration_chunk_size):
        print process_iteration(i, compute_bases_overlapped, gaps, bp_file, genome, feature_file)
else:
    sys.stderr.write("Bad process_iteration commmand!\n")
#sys.stderr.write("all done\n")
