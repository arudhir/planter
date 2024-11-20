import requests
import logging
import xml.etree.ElementTree as ET
from typing import Dict, Optional
import time
from functools import wraps
import random

logger = logging.getLogger(__name__)

def rate_limit(calls: int = 1, period: int = 1):
    """
    Rate limiting decorator.
    
    Args:
        calls: Number of calls allowed per period
        period: Time period in seconds
    """
    min_time_between_calls = period / calls
    last_call_time = [0.0]  # List to make it mutable in closure
    
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            elapsed = time.time() - last_call_time[0]
            if elapsed < min_time_between_calls:
                time.sleep(min_time_between_calls - elapsed)
            result = func(*args, **kwargs)
            last_call_time[0] = time.time()
            return result
        return wrapper
    return decorator

def retry_with_backoff(retries=3, backoff_in_seconds=1):
    """
    Retry decorator with exponential backoff.
    
    Args:
        retries: Number of retries
        backoff_in_seconds: Initial backoff time
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Initialize variables
            num_tries = 0
            sleep_time = backoff_in_seconds
            
            while True:
                try:
                    return func(*args, **kwargs)
                except requests.exceptions.RequestException as e:
                    num_tries += 1
                    if num_tries >= retries:
                        raise
                    
                    # If we get a 429, use the Retry-After header if available
                    if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 429:
                        sleep_time = int(e.response.headers.get('Retry-After', sleep_time))
                    
                    # Add some random jitter to prevent thundering herd
                    jitter = random.uniform(0, 0.1 * sleep_time)
                    total_sleep = sleep_time + jitter
                    
                    logger.warning(
                        f"Attempt {num_tries} failed. Retrying in {total_sleep:.1f} seconds... "
                        f"Error: {str(e)}"
                    )
                    
                    time.sleep(total_sleep)
                    # Exponential backoff
                    sleep_time *= 2
        return wrapper
    return decorator

class NCBIRateLimiter:
    """Rate limits for NCBI API"""
    def __init__(self, calls_per_second=3):
        self.calls_per_second = calls_per_second
        self.last_call = 0
    
    def wait(self):
        """Wait if necessary before making next call"""
        elapsed = time.time() - self.last_call
        min_time_between_calls = 1.0 / self.calls_per_second
        
        if elapsed < min_time_between_calls:
            time.sleep(min_time_between_calls - elapsed)
        
        self.last_call = time.time()

# Global rate limiter
ncbi_limiter = NCBIRateLimiter(calls_per_second=2)  # 2 calls per second to be safe

@retry_with_backoff(retries=3, backoff_in_seconds=1)
@rate_limit(calls=2, period=1)  # 2 calls per second
def fetch_ncbi_data(url: str) -> str:
    """Fetch data from NCBI with rate limiting"""
    ncbi_limiter.wait()
    response = requests.get(url)
    response.raise_for_status()
    return response.content

def get_sra_info(srr_id: str) -> Optional[Dict]:
    """
    Retrieve comprehensive information for a given SRR ID from NCBI's SRA database.
    
    Args:
        srr_id: The SRR ID to look up
        
    Returns:
        Dictionary containing metadata or None if fetch fails
    """
    logger.info(f"Fetching metadata for {srr_id}")
    
    try:
        # Step 1: Use esearch to find the SRA record
        esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={srr_id}"
        content = fetch_ncbi_data(esearch_url)
        root = ET.fromstring(content)
        id_list = root.find('IdList')
        
        if id_list is None or len(id_list) == 0:
            logger.warning(f"No record found for {srr_id}")
            return None
            
        sra_id = id_list.find('Id').text

        # Step 2: Use efetch to retrieve the full record
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={sra_id}"
        content = fetch_ncbi_data(efetch_url)
        root = ET.fromstring(content)

        # Helper function to safely extract text from XML elements
        def safe_find_text(element, path, default="Not found"):
            found = element.find(path)
            return found.text if found is not None else default

        info = {}

        # Study information
        study = root.find(".//STUDY")
        if study is not None:
            info['study_title'] = safe_find_text(study, ".//STUDY_TITLE")
            info['study_abstract'] = safe_find_text(study, ".//STUDY_ABSTRACT")
            info['bioproject'] = safe_find_text(study, ".//EXTERNAL_ID[@namespace='BioProject']")

        # Sample information
        sample = root.find(".//SAMPLE")
        if sample is not None:
            info['organism'] = safe_find_text(sample, ".//SCIENTIFIC_NAME")
            info['biosample'] = safe_find_text(sample, ".//EXTERNAL_ID[@namespace='BioSample']")

        # Experiment information
        experiment = root.find(".//EXPERIMENT")
        if experiment is not None:
            library = experiment.find(".//LIBRARY_DESCRIPTOR")
            if library is not None:
                layout = library.find("LIBRARY_LAYOUT")
                layout_type = next(iter(layout)) if layout is not None and len(layout) > 0 else None
                info['library'] = {
                    'strategy': safe_find_text(library, "LIBRARY_STRATEGY"),
                    'source': safe_find_text(library, "LIBRARY_SOURCE"),
                    'selection': safe_find_text(library, "LIBRARY_SELECTION"),
                    'layout': layout_type.tag if layout_type is not None else "Not found"
                }

            platform = experiment.find(".//PLATFORM")
            if platform is not None:
                instrument = next(iter(platform), None)
                info['instrument'] = instrument.tag if instrument is not None else "Not found"

        # Run information
        run = root.find(".//RUN")
        if run is not None:
            info['run'] = {
                'spots': run.get('total_spots'),
                'bases': run.get('total_bases'),
                'published': run.get('published')
            }

        return info

    except Exception as e:
        logger.error(f"Error fetching data for {srr_id}: {str(e)}")
        return None