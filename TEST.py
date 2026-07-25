import sqlite3
import itertools
import json
from datetime import datetime
# Adjust import based on your exact project structure
from dna_matcher.compare import compare_species_matrix 

def precompute_comparisons(db_path="dna_cache.sqlite3", max_len=600):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Retrieve all species that successfully cached a COI marker[cite: 1]
    cursor.execute("SELECT species_name, sequence FROM marker_sequences WHERE gene = 'COI'")
    records = cursor.fetchall()
    
    species_list = [row[0] for row in records]
    sequences = {row[0]: row[1] for row in records}

    print(f"Found {len(species_list)} species with COI markers.")
    
    # Generate all unique pairs to avoid redundant A->B and B->A calculations
    pairs = list(itertools.combinations(species_list, 2))
    print(f"Total comparisons to compute: {len(pairs)}")

    for species_a, species_b in pairs:
        seq_a = sequences[species_a]
        seq_b = sequences[species_b]
        
        try:
            # Execute your existing comparison logic
            result = compare_species_matrix(seq_a, seq_b, max_len=max_len)
            
            # Serialize the results dictionary to match the 'details' TEXT column[cite: 1]
            details = json.dumps(result)
            computed_at = datetime.utcnow().isoformat()
            
            # Insert or update the comparison_results table[cite: 1]
            cursor.execute("""
                INSERT INTO comparison_results 
                (species_a, species_b, marker, max_len, details, computed_at)
                VALUES (?, ?, ?, ?, ?, ?)
                ON CONFLICT(species_a, species_b, marker, max_len) DO UPDATE SET
                    details=excluded.details,
                    computed_at=excluded.computed_at
            """, (species_a, species_b, "COI", max_len, details, computed_at))
            
        except Exception as e:
            print(f"Error computing {species_a} vs {species_b}: {e}")

    # Commit all transactions and close
    conn.commit()
    conn.close()
    print("All comparisons computed and cached successfully.")

if __name__ == "__main__":
    precompute_comparisons()