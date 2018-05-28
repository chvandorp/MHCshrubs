"""@package progressbar
A progress bar for a list of jobs that is run in parallel.

Credits go to Dan Shiebler (http://danshiebler.com/2016-09-14-parallel-progress-bar/)
"""
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def parallel_map(fun, jobs, max_workers=1):
    """
    A parallel version of the map function with a progress bar.

    Args:
        fun (function): A python function to apply to the elements of the array jobs
        jobs (array-like): An array to iterate over.

    Kwargs:
        max_workers (int, default=1): The number of cores to use

    Returns:
        a list with fun mapped over jobs
    """
    ## If we set num_workers to 1, just run a list comprehension.
    ## This is useful for benchmarking and debugging.
    if max_workers == 1:
        return [fun(job) for job in tqdm(jobs)]
    ## Assemble the workers
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        ## Pass the elements of array into function
        futures = [pool.submit(fun, job) for job in jobs]
        ## Print out the progress as tasks complete
        for f in tqdm(as_completed(futures), total=len(futures)):
            pass
    results = []
    ## Get the results from the futures.
    for i, future in enumerate(futures):
        results.append(future.result())
    return results
