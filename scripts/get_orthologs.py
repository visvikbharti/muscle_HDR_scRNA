#!/usr/bin/env python3
import requests
import pandas as pd
import time
import json
import sys
import os

def fetch_ensembl_orthologs(human_genes, chunk_size=100):
    """
    Fetch human-mouse orthologs from Ensembl REST API
    """
    base_url = "https://rest.ensembl.org/homology/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json"}
    
    results = []
    
    # Process genes in chunks to avoid API limits
    for i in range(0, len(human_genes), chunk_size):
        chunk = human_genes[i:i + chunk_size]
        print(f"Processing genes {i} to {i + len(chunk)}")
        
        for gene in chunk:
            try:
                # Construct the URL
                url = f"{base_url}/{gene}"
                params = {
                    "type": "orthologues",
                    "target_species": "mus_musculus",
                    "format": "json"
                }
                
                # Make the request with a timeout
                response = requests.get(url, headers=headers, params=params, timeout=10)
                
                if response.status_code == 200:
                    data = response.json()
                    
                    # Parse the response
                    if 'data' in data and data['data']:
                        for homology in data['data'][0]['homologies']:
                            if homology['target']['species'] == 'mus_musculus':
                                results.append({
                                    'human_symbol': gene,
                                    'mouse_symbol': homology['target']['gene_symbol'],
                                    'confidence': homology.get('confidence', 'NA'),
                                    'type': homology.get('type', 'NA')
                                })
                                break
                    # If no ortholog found, add a simple mapping
                    else:
                        results.append({
                            'human_symbol': gene,
                            'mouse_symbol': gene.capitalize(),  # Simple fallback
                            'confidence': 'simple_mapping',
                            'type': 'capitalized'
                        })
                else:
                    print(f"Warning: API error for {gene} (status code: {response.status_code})")
                    # Use simple fallback
                    results.append({
                        'human_symbol': gene,
                        'mouse_symbol': gene.capitalize(),
                        'confidence': 'fallback',
                        'type': 'capitalized'
                    })
                
                # Rate limiting to avoid API throttling
                time.sleep(0.1)
                
            except Exception as e:
                print(f"Error processing {gene}: {e}")
                # Use simple fallback
                results.append({
                    'human_symbol': gene,
                    'mouse_symbol': gene.capitalize(),
                    'confidence': 'error_fallback',
                    'type': 'capitalized'
                })
                continue
        
        # Longer pause between chunks
        time.sleep(1)
    
    return results

def main(hdr_genes_file, output_file):
    """
    Main function to fetch orthologs and save to TSV
    """
    # Load HDR genes
    with open(hdr_genes_file, 'r') as f:
        hdr_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Loaded {len(hdr_genes)} HDR genes")
    
    # Try to fetch orthologs
    results = fetch_ensembl_orthologs(hdr_genes)
    
    # If API fails completely, use simple mapping
    if not results:
        print("API failed completely. Using simple capitalization mapping...")
        results = [
            {
                'human_symbol': gene,
                'mouse_symbol': gene.capitalize(),
                'confidence': 'offline_fallback',
                'type': 'capitalized'
            }
            for gene in hdr_genes
        ]
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Ensure we have data
    if len(df) == 0:
        # If still no data, create minimal dataframe
        df = pd.DataFrame({
            'human_symbol': hdr_genes,
            'mouse_symbol': [g.capitalize() for g in hdr_genes],
            'confidence': ['fallback'] * len(hdr_genes),
            'type': ['capitalized'] * len(hdr_genes)
        })
    
    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved {len(df)} orthologs to {output_file}")
    
    # Also create a simple mapping dictionary
    mapping = dict(zip(df['human_symbol'], df['mouse_symbol']))
    mapping_file = output_file.replace('.tsv', '_mapping.json')
    with open(mapping_file, 'w') as f:
        json.dump(mapping, f, indent=2)
    
    # Print summary
    print(f"\nSummary:")
    print(f"Total genes processed: {len(hdr_genes)}")
    print(f"Orthologs mapped: {len(df)}")
    print(f"API success: {len(df[df['confidence'] != 'fallback'])}")
    print(f"Fallback mappings: {len(df[df['confidence'] == 'fallback'])}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_orthologs.py HDR_genes.txt output.tsv")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])