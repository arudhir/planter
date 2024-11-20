import requests
import logging
import xml.etree.ElementTree as ET
from typing import Dict, Optional
import time

def get_sra_info(srr_id: str) -> Optional[Dict]:
    """
    Retrieve comprehensive information for a given SRR ID from NCBI's SRA database.
    
    Args:
        srr_id: The SRR ID to look up
        
    Returns:
        Dictionary containing the retrieved information or None if error occurs
    """
    logger = logging.getLogger(__name__)
    
    try:
        # Step 1: Use esearch to find the SRA record
        esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={srr_id}"
        response = requests.get(esearch_url)
        response.raise_for_status()
        root = ET.fromstring(response.content)
        id_list = root.find('IdList')
        
        if id_list is None or len(id_list) == 0:
            logger.error(f"No record found for {srr_id}")
            return None
            
        sra_id = id_list.find('Id').text

        # Step 2: Use efetch to retrieve the full record
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={sra_id}"
        response = requests.get(efetch_url)
        response.raise_for_status()
        root = ET.fromstring(response.content)

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

    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching data for {srr_id}: {str(e)}")
        return None
    except ET.ParseError as e:
        logger.error(f"Error parsing XML data for {srr_id}: {str(e)}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error occurred for {srr_id}: {str(e)}")
        return None