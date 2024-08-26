#!/usr/bin/env python
import requests
import xml.etree.ElementTree as ET
import argparse
from pprint import pprint
import time

def get_sra_info(srr_id):
    """
    Retrieve comprehensive information for a given SRR ID from NCBI's SRA database.

    Args:
    srr_id (str): The SRR ID to look up.

    Returns:
    dict: A dictionary containing the retrieved information, or
    str: An error message if the record is not found or an error occurs.
    """
    try:
        # Step 1: Use esearch to find the SRA record
        esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={srr_id}"
        response = requests.get(esearch_url)
        response.raise_for_status()
        root = ET.fromstring(response.content)
        id_list = root.find('IdList')
        if id_list is None or len(id_list) == 0:
            return f"No record found for {srr_id}"
        sra_id = id_list.find('Id').text

        # Step 2: Use efetch to retrieve the full record
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={sra_id}"
        response = requests.get(efetch_url)
        response.raise_for_status()
        root = ET.fromstring(response.content)

        # Step 3: Parse the XML to find all required information
        info = {}

        # Helper function to safely extract text from XML elements
        def safe_find_text(element, path, default="Not found"):
            found = element.find(path)
            return found.text if found is not None else default

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
                    'name': safe_find_text(library, "LIBRARY_NAME"),
                    'strategy': safe_find_text(library, "LIBRARY_STRATEGY"),
                    'source': safe_find_text(library, "LIBRARY_SOURCE"),
                    'selection': safe_find_text(library, "LIBRARY_SELECTION"),
                    'layout': layout_type.tag if layout_type is not None else "Not found"
                }

            platform = experiment.find(".//PLATFORM")
            if platform is not None:
                # The instrument is the first child of the PLATFORM element
                instrument = next(iter(platform), None)
                info['instrument'] = instrument.tag if instrument is not None else "Not found"

        # Run information
        run = root.find(".//RUN")
        if run is not None:
            info['run'] = {
                'accession': run.get('accession'),
                'spots': run.get('total_spots'),
                'bases': run.get('total_bases'),
                'size': run.get('size'),
                'published': run.get('published')
            }

        return info

    except requests.exceptions.RequestException as e:
        return f"Error fetching data for {srr_id}: {str(e)}"
    except ET.ParseError as e:
        return f"Error parsing XML data for {srr_id}: {str(e)}"
    except Exception as e:
        return f"Unexpected error occurred for {srr_id}: {str(e)}"

def main():
    """
    Main function to parse command-line arguments and display SRA information for one or multiple SRR IDs.
    """
    parser = argparse.ArgumentParser(description="Retrieve comprehensive SRA information for one or multiple SRR IDs")
    parser.add_argument("--srr", nargs='+', required=True, help="One or more SRR IDs to look up")
    args = parser.parse_args()

    for srr_id in args.srr:
        try:
            info = get_sra_info(srr_id)
            print(f"\nInformation for {srr_id}:")
            if isinstance(info, str):
                print(info)
            else:
                pprint(info, width=100, sort_dicts=False)
        except Exception as e:
            print(f"\nError processing {srr_id}: {str(e)}")
        finally:
            print("-" * 80)  # Separator between multiple SRR outputs

        # Add a delay between processing each SRR ID
        time.sleep(1)

if __name__ == "__main__":
    main()
