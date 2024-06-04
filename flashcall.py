#!/usr/bin/env python3
import pyBigWig
import os
import concurrent.futures
import statistics
import argparse
from pybedtools import BedTool
import scipy.stats as stats
import random
import numpy as np
import time
import random
from concurrent.futures import ThreadPoolExecutor



# Specify a custom directory for temporary files
custom_temp_dir = '/path/to/define'

# Set the TEMP environment variable to the custom directory
os.environ['TEMP'] = custom_temp_dir

window_size = 10000  # Set your window size here


def create_background_regions_parallel(chrom_lengths, valid_chromosomes, num_regions=10000, region_size=10000):
    background_regions = []
    
    def create_regions_for_chromosome(chrom):
        chrom_length = chrom_lengths[chrom]
        regions = []
        while len(regions) < num_regions_per_chrom:
            if chrom_length > region_size:
                start = random.randint(0, chrom_length - region_size)
                end = start + region_size
                regions.append((chrom, start, end))
        return regions

    num_regions_per_chrom = num_regions // len(valid_chromosomes)
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(create_regions_for_chromosome, chrom) for chrom in valid_chromosomes]
        for future in concurrent.futures.as_completed(futures):
            background_regions.extend(future.result())

    return background_regions[:num_regions]
    

def calculate_statistics_for_region(bigwig_file, region):
    chrom, start, end = region
    with pyBigWig.open(bigwig_file) as bw:
        values = bw.values(chrom, start, end, numpy=True)
    return np.mean(values), np.std(values)

def calculate_background_statistics(bigwig_file, background_regions):
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(calculate_statistics_for_region, bigwig_file, region) for region in background_regions]
        results = [future.result() for future in futures]

    mean_values = [res[0] for res in results]
    std_values = [res[1] for res in results]

    final_mean = np.mean(mean_values)
    final_std = np.mean(std_values)

    return final_mean, final_std

def is_in_blacklist(chrom, start, end, blacklist):
    for blk_chrom, blk_start, blk_end in blacklist:
        if chrom == blk_chrom and blk_start < end and blk_end > start:
            return True
    return False

def calculate_p_value_for_interval(interval, background_mean, background_std):
    z_score = (interval[3] - background_mean) / background_std
    return interval, stats.norm.sf(abs(z_score))  # Two-tailed p-value

def calculate_p_values_parallel(intervals, background_mean, background_std):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_interval = {executor.submit(calculate_p_value_for_interval, interval, background_mean, background_std): interval for interval in intervals}
        results = {}
        for future in concurrent.futures.as_completed(future_to_interval):
            interval = future_to_interval[future]
            try:
                results[interval] = future.result()
            except Exception as exc:
                print(f'{interval} generated an exception: {exc}')
        return results

def calculate_auc(bigwig_file, interval):
    chrom, start, end = interval[:3]
    with pyBigWig.open(bigwig_file) as bw:
        values = bw.values(chrom, start, end)
    auc = sum(values)  # or any other method to calculate AUC
    return auc

def write_auc_values(bigwig_file, intervals, output_file):
    with open(output_file, 'w') as f:
        for interval in intervals:
            auc = calculate_auc(bigwig_file, interval)
            f.write(f'{interval[0]}\t{interval[1]}\t{interval[2]}\t{auc}\n')

def process_chromosome(args):
    chromosome, length, valid_chromosomes, bigwig_file, threshold, blacklist = args

    if chromosome not in valid_chromosomes:
        return []
    intervals = []
    with pyBigWig.open(bigwig_file) as bw:
        start = None
        max_val = None
        max_pos = None
        last_val = None

        for i in range(0, length, window_size):
            values = bw.values(chromosome, i, min(i + window_size, length))
            for j, value in enumerate(values):
                if value > threshold:  # Start interval only if value is greater than threshold
                    if start is None:
                        start = i + j
                        max_val = value
                        max_pos = j
                    elif value > max_val:
                        max_val = value
                        max_pos = i + j - start
                    last_val = value
                elif start is not None and last_val > threshold:
                    end = i + j
                    if end - start > 1:  # Check for intervals longer than a single point
                        if not is_in_blacklist(chromosome, start, end, blacklist):
                            intervals.append((chromosome, start, end, max_val, max_pos))
                    # Reset for next interval
                    start = None
                    max_val = None
                    max_pos = None
                    last_val = None

        # Check for the last interval in the chromosome
        if start is not None and last_val > threshold:
            end = length
            if end - start > 1:
                if not is_in_blacklist(chromosome, start, end, blacklist):
                    intervals.append((chromosome, start, end, max_val, max_pos))

    return intervals


    
    
def fuse_intervals(intervals, distance_threshold):
    fused_intervals = []
    current_interval = None

    for interval in sorted(intervals, key=lambda x: (x[0], x[1])):
        chrom, start, end, max_val, max_pos = interval

        if current_interval and chrom == current_interval[0] and start - current_interval[2] <= distance_threshold:
            # Fuse intervals
            new_end = max(end, current_interval[2])
            new_max_val = max(max_val, current_interval[3])
            new_max_pos = max_pos if new_max_val == max_val else current_interval[4]
            current_interval = (chrom, current_interval[1], new_end, new_max_val, new_max_pos)
        else:
            # Save the previous interval and start a new one
            if current_interval:
                fused_intervals.append(current_interval)
            current_interval = interval

    # Add the last interval
    if current_interval:
        fused_intervals.append(current_interval)

    return fused_intervals

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Process BigWig files.')
    parser.add_argument('-b', '--bigwig', required=True, help='Path to the BigWig file to analyze.')
    parser.add_argument('-o', '--output', required=True, help='Path and name for the output file.')
    parser.add_argument('-g', '--genome', required=True, choices=['human', 'mouse'], help='Genome type (human or mouse).')
    parser.add_argument('-p', '--pvalue_threshold', type=float, default=0.05, help='P-value threshold for significance.')
    parser.add_argument('-t', '--threshold', type=float, default=1.0, help='Threshold value for processing (default: 1)')
    parser.add_argument('-d', '--distance', type=int, default=50, help='Distance threshold for fusing intervals (default: 50 bp)')
    parser.add_argument('-a', '--auc_threshold', type=float, default=5, help='AUC threshold for filtering (default: 5.0)')


    args = parser.parse_args()

    bigwig_file = args.bigwig
    output_file = args.output

    # Determine the path to the blacklist and set valid chromosomes
    if args.genome == 'human':
        blacklist_path = 'path'
    elif args.genome == 'mouse':
        blacklist_path = 'path'

    if args.genome == 'human':
        valid_chromosomes = ['chr' + str(i) for i in range(1, 22)] + ['chrX']
    elif args.genome == 'mouse':
        valid_chromosomes = ['chr' + str(i) for i in range(1, 20)] + ['chrX']
    # Load blacklist regions from the determined path
    blacklist_bed = BedTool(blacklist_path)
    blacklist_regions = [(interval.chrom, int(interval.start), int(interval.end)) for interval in blacklist_bed]    
    # Open the BigWig file and get chromosomes information
    chrom_lengths = pyBigWig.open(bigwig_file).chroms()

    # Prepare arguments for the process_chromosome function
    process_args = [(chrom, length, valid_chromosomes, bigwig_file, args.threshold, blacklist_regions)
                for chrom, length in chrom_lengths.items() 
                if chrom in valid_chromosomes]

    start_time = time.time()
    # Process each chromosome in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        all_results = executor.map(process_chromosome, process_args)
    end_time = time.time()
    print(f"Time taken for process_chromosome: {end_time - start_time} seconds")

    start_time = time.time()
    # Flatten the list of results
    all_intervals = [interval for sublist in all_results for interval in sublist]
    end_time = time.time()
    print(f"Time taken for flatten_chromosome: {end_time - start_time} seconds")    

 
    start_time = time.time()    
    # Fuse intervals based on distance threshold
    fused_intervals = fuse_intervals(all_intervals, args.distance)
    end_time = time.time()
    print(f"Time taken for fuse_intervals: {end_time - start_time} seconds") 
    
    # Calculate background statistics
    start_time = time.time()       
    background_regions = create_background_regions_parallel(chrom_lengths, valid_chromosomes, 10000, 10000)


    end_time = time.time()
    print(f"Time taken for background_regions: {end_time - start_time} seconds")    

    start_time = time.time()    
    background_mean, background_std = calculate_background_statistics(bigwig_file, background_regions)
    end_time = time.time()
    print(f"Time taken for background_mean: {end_time - start_time} seconds")    



    # Calculate p-values for each interval
    start_time = time.time()
    interval_pvalues = calculate_p_values_parallel(fused_intervals, background_mean, background_std)
    end_time = time.time()
    print(f"Time taken for interval_pvalues: {end_time - start_time} seconds")
    # Filter intervals based on p-values and blacklist
    start_time = time.time()
    significantly_enriched_intervals = [interval for interval, (interval_inner, p_value) in interval_pvalues.items() if p_value < args.pvalue_threshold]
    end_time = time.time()
    print(f"Time taken for significantly_enriched_intervals: {end_time - start_time} seconds")


    start_time = time.time()
    # Split intervals by chromosome
    intervals_by_chromosome = {}
    for interval in fused_intervals:
        chrom = interval[0]
        if chrom not in intervals_by_chromosome:
            intervals_by_chromosome[chrom] = []
        intervals_by_chromosome[chrom].append(interval)


    # Calculate AUC for each interval
    interval_aucs = {interval: calculate_auc(bigwig_file, interval) for interval in fused_intervals}
    end_time = time.time()
    print(f"Time taken for filtered_blacklist: {end_time - start_time} seconds")
    
    start_time = time.time()
    # Filter intervals based on p-values, AUC and blacklist
    filtered_intervals = []
    for interval in fused_intervals:
        p_value = interval_pvalues[interval][1]
        auc = interval_aucs[interval]
        if p_value < args.pvalue_threshold and auc > args.auc_threshold:
            filtered_intervals.append(interval)
    end_time = time.time()
    print(f"Time taken for filtered_interval_pval_auc: {end_time - start_time} seconds")


    start_time = time.time()
    # Write filtered intervals and their p-values to files
    enriched_output_file = output_file.replace('.bed', '_enriched.bed')
    pvalue_output_file = output_file.replace('.bed', '_pvalues.txt')

    with open(enriched_output_file, 'w') as f:
        for interval in filtered_intervals:
            f.write(f'{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\t{interval[4]}\n')

    with open(pvalue_output_file, 'w') as f:
        for interval in filtered_intervals:
            p_value = interval_pvalues[interval][1]
            auc = interval_aucs[interval]
            f.write(f'{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\t{interval[4]}\t{p_value}\t{auc}\n')

    end_time = time.time()
    print(f"Time taken for write_intervals: {end_time - start_time} seconds")

if __name__ == '__main__':
    main()
