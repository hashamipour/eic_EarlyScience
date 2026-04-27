#!/usr/bin/env python3
"""
XRootD file downloader with parallel downloads and random selection.
Usage: python xrootd_downloader.py <links_file> <num_files> [--workers N] [--output-dir DIR]
"""

import argparse
import random
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def download_file(url, output_dir):
    """Download a single file using xrdcp."""
    filename = url.split('/')[-1]
    output_path = os.path.join(output_dir, filename)
    
    print(f"Starting: {filename}")
    
    try:
        result = subprocess.run(
            ['xrdcp', '-f', url, output_path],
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            print(f"Done: {filename}")
            return (url, True, None)
        else:
            print(f"Failed: {filename} - {result.stderr.strip()}")
            return (url, False, result.stderr)
            
    except Exception as e:
        print(f"Error: {filename} - {str(e)}")
        return (url, False, str(e))


def main():
    parser = argparse.ArgumentParser(description='Download xrootd files in parallel')
    parser.add_argument('links_file', help='Text file with xrootd links (one per line)')
    parser.add_argument('num_files', type=int, help='Number of files to download')
    parser.add_argument('--workers', '-w', type=int, default=4, 
                        help='Number of parallel downloads (default: 4)')
    parser.add_argument('--output-dir', '-o', default='.', 
                        help='Output directory (default: current dir)')
    parser.add_argument('--seed', '-s', type=int, default=None,
                        help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # read links
    with open(args.links_file, 'r') as f:
        links = [line.strip() for line in f if line.strip()]
    
    if len(links) == 0:
        print("No links found in file!")
        return
    
    # check if we have enough links
    if args.num_files > len(links):
        print(f"Warning: requested {args.num_files} but only {len(links)} available")
        args.num_files = len(links)
    
    # random selection
    if args.seed is not None:
        random.seed(args.seed)
    selected = random.sample(links, args.num_files)
    
    print(f"Selected {len(selected)} files randomly from {len(links)} total")
    print(f"Using {args.workers} parallel workers")
    print(f"Output directory: {args.output_dir}\n")
    
    # create output dir
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # download in parallel
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_file, url, args.output_dir): url 
            for url in selected
        }
        
        for future in as_completed(futures):
            results.append(future.result())
    
    # summary
    success = sum(1 for _, ok, _ in results if ok)
    failed = len(results) - success
    
    print(f"\n{'='*40}")
    print(f"Downloaded: {success}/{len(selected)}")
    if failed > 0:
        print(f"Failed: {failed}")
        for url, ok, err in results:
            if not ok:
                print(f"  - {url.split('/')[-1]}")


if __name__ == '__main__':
    main()
