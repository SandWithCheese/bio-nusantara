import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
from Bio import Entrez
import time
import os # Pastikan os diimpor di atas
import json # Untuk menyimpan dan memuat progres

# ==============================================================================
# WAJIB: Ganti ini dengan alamat email aktif Anda sebelum menjalankan skrip!
# NCBI memerlukan ini untuk melacak penggunaan dan menghubungi Anda jika ada masalah.
Entrez.email = "irvansyahshafiqs@gmail.com"
# ==============================================================================

class DataProcessor:
    def __init__(self):
        self.species_data = None
        self.parks_data = None

    def _fetch_ncbi_sequence(self, species_name):
        """
        Fetches a nucleotide sequence for a given species name from NCBI Entrez.
        Tries to fetch one relevant sequence (e.g., mitochondrial or COI if common).
        """
        if pd.isna(species_name) or not species_name.strip():
            return None
        
        search_terms = [
            f"{species_name}[Organism] AND (COI[Gene Name] OR COX1[Gene Name] OR CO1[Gene Name] OR 'cytochrome c oxidase subunit I'[Gene Name]) AND mitochondrion[Filter]",
            f"{species_name}[Organism] AND (CytB[Gene Name] OR 'cytochrome b'[Gene Name]) AND mitochondrion[Filter]",
            f"{species_name}[Organism] AND mitochondrion[Filter] AND complete[Title]",
            f"{species_name}[Organism] AND complete genome[Title]",
            f"{species_name}[Organism]"
        ]

        seq_id = None
        for term_idx, term in enumerate(search_terms):
            try:
                handle_search = Entrez.esearch(db="nucleotide", term=term, retmax=1, idtype="acc")
                record_search = Entrez.read(handle_search)
                handle_search.close()
                ids = record_search["IdList"]
                if ids:
                    seq_id = ids[0]
                    break
            except Exception as e_search:
                print(f"      Error during esearch for {species_name} with term '{term[:30]}...': {e_search}")
                if "timeout" in str(e_search).lower() or "NameResolutionError" in str(e_search).lower() or "http" in str(e_search).lower():
                    print("      Retrying esearch after a short delay...")
                    time.sleep(5) # Jeda lebih lama untuk error jaringan
                    try:
                        handle_search = Entrez.esearch(db="nucleotide", term=term, retmax=1, idtype="acc")
                        record_search = Entrez.read(handle_search)
                        handle_search.close()
                        ids = record_search["IdList"]
                        if ids: seq_id = ids[0]; break
                    except Exception as e_retry:
                        print(f"      Retry failed for {species_name}: {e_retry}")
                continue

        if not seq_id:
            # print(f"      No sequence ID found for {species_name} after trying all search terms.") # Dihilangkan agar tidak terlalu verbose
            return None

        try:
            handle_fetch = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            fasta_sequence = handle_fetch.read()
            handle_fetch.close()
            if fasta_sequence and fasta_sequence.startswith(">"):
                return fasta_sequence.split('\n', 1)[1].replace('\n', '') if len(fasta_sequence.split('\n', 1)) > 1 else fasta_sequence
            return fasta_sequence
        except Exception as e_fetch:
            print(f"      Error fetching FASTA for {species_name} (ID: {seq_id}): {e_fetch}")
            return None

    def _load_progress(self, progress_file_path):
        """Memuat progres pengambilan sekuens dari file JSON."""
        if os.path.exists(progress_file_path):
            try:
                with open(progress_file_path, 'r') as f:
                    # Pastikan file tidak kosong sebelum memuat JSON
                    content = f.read()
                    if not content:
                        print(f"Warning: File progres {progress_file_path} kosong, memulai dari awal.")
                        return {}
                    return json.loads(content)
            except json.JSONDecodeError:
                print(f"Warning: File progres {progress_file_path} rusak atau bukan JSON valid, memulai dari awal.")
                # Anda bisa mencoba mem-backup file yang rusak di sini jika perlu
                # os.rename(progress_file_path, progress_file_path + ".corrupted")
                return {}
            except Exception as e:
                print(f"Error memuat file progres {progress_file_path}: {e}. Memulai dari awal.")
                return {}
        return {}

    def _save_progress(self, progress_data, progress_file_path):
        """Menyimpan progres pengambilan sekuens ke file JSON."""
        try:
            # Pastikan direktori untuk file progres ada
            progress_dir = os.path.dirname(progress_file_path)
            if progress_dir and not os.path.exists(progress_dir):
                os.makedirs(progress_dir)
                print(f"Created directory for progress file: {progress_dir}")

            with open(progress_file_path, 'w') as f:
                json.dump(progress_data, f, indent=4)
        except Exception as e:
            print(f"Error menyimpan file progres {progress_file_path}: {e}")


    def load_and_clean_gbif_data(self, file_paths, fetch_sequences=False, progress_file_path_override=None):
        """Load and clean GBIF occurrence data from multiple CSV files"""
        all_data = []
        
        for file_path in file_paths:
            print(f"Loading {file_path}...")
            try:
                df = pd.read_csv(file_path, low_memory=False)
                print(f"Loaded {len(df)} records from {file_path}")
                all_data.append(df)
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
        
        if not all_data:
            print("No data files could be loaded. Returning empty GeoDataFrame.")
            return gpd.GeoDataFrame(columns=['geometry'], geometry='geometry', crs='EPSG:4326')

        combined_df = pd.concat(all_data, ignore_index=True)
        print(f"Combined dataset: {len(combined_df)} total records")
        
        cleaned_df = self.clean_gbif_data(combined_df, fetch_sequences=fetch_sequences, progress_file_path_override=progress_file_path_override)
        
        self.species_data = self.create_geodataframe(cleaned_df)
        
        return self.species_data
    
    def clean_gbif_data(self, df, fetch_sequences=False, progress_file_path_override=None):
        """Clean GBIF data based on the actual format"""
        print("Cleaning GBIF data...")
        
        if 'key' in df.columns:
            initial_count_dup = len(df)
            df = df.drop_duplicates(subset=['key'], keep='first')
            print(f"Removed {initial_count_dup - len(df)} duplicate records based on 'key'")

        initial_count = len(df)
        df = df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
        print(f"Removed {initial_count - len(df)} records without coordinates")
        
        df['decimalLatitude'] = pd.to_numeric(df['decimalLatitude'], errors='coerce')
        df['decimalLongitude'] = pd.to_numeric(df['decimalLongitude'], errors='coerce')
        df = df.dropna(subset=['decimalLatitude', 'decimalLongitude'])

        if 'coordinateUncertaintyInMeters' in df.columns:
            df['coordinateUncertaintyInMeters'] = pd.to_numeric(df['coordinateUncertaintyInMeters'], errors='coerce')
            high_quality = df['coordinateUncertaintyInMeters'].fillna(9999) < 1000 
            df = df[high_quality]
            print(f"Kept {len(df)} records with coordinate uncertainty < 1000m (NaNs treated as >1000m)")
        
        valid_coords = (
            (df['decimalLatitude'].between(-90, 90)) & 
            (df['decimalLongitude'].between(-180, 180))
        )
        df = df[valid_coords]
        print(f"Kept {len(df)} records with valid coordinates after numeric conversion")
        
        indonesia_bounds = (
            (df['decimalLatitude'].between(-11, 6)) &
            (df['decimalLongitude'].between(95, 141))
        )
        df = df[indonesia_bounds]
        print(f"Kept {len(df)} records within Indonesia bounds")
        
        essential_columns = [
            'key', 'scientificName', 'acceptedScientificName', 'taxonKey',
            'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
            'decimalLatitude', 'decimalLongitude', 'coordinateUncertaintyInMeters',
            'year', 'month', 'day', 'eventDate',
            'stateProvince', 'locality', 'country', 'countryCode',
            'basisOfRecord', 'occurrenceStatus', 'iucnRedListCategory',
            'datasetName', 'recordedBy', 'identifiedBy'
        ]
        
        available_columns = [col for col in essential_columns if col in df.columns]
        df = df[available_columns].copy()
        
        df['taxonomic_group'] = df['class'].apply(self.classify_taxonomic_group)
        
        if fetch_sequences:
            print("\nFetching sequences from NCBI Entrez...")
            if 'scientificName' not in df.columns:
                print("Warning: 'scientificName' column not found. Cannot fetch sequences.")
                df['gen_seq'] = "Sequence Fetching Skipped"
            else:
                unique_species_names = df['scientificName'].dropna().unique()
                total_unique_species = len(unique_species_names)
                
                # Path untuk file progres
                if progress_file_path_override:
                    progress_file = progress_file_path_override
                else:
                    # Default path jika tidak di-override
                    # Menggunakan direktori data/processed relatif terhadap lokasi eksekusi skrip
                    # Ini mungkin perlu disesuaikan jika skrip tidak dijalankan dari direktori 'scripts'
                    current_script_path = os.path.dirname(os.path.abspath(__file__))
                    default_progress_dir = os.path.join(current_script_path, '..', 'data', 'processed')
                    progress_file = os.path.join(default_progress_dir, 'sequence_fetch_progress.json')

                species_sequence_map = self._load_progress(progress_file)
                
                processed_count_this_session = 0
                newly_fetched_count = 0
                print(f"Ditemukan {total_unique_species} nama spesies unik.")
                print(f"{len(species_sequence_map)} spesies sudah ada di file progres ({progress_file}).")
                print("PENTING: Proses ini bisa sangat lambat...")
                print(f"Pastikan Entrez.email ({Entrez.email}) sudah diatur dengan benar.")

                for i, species_name in enumerate(unique_species_names):
                    if species_name in species_sequence_map and species_sequence_map[species_name] != "Not Found or Error": # Hanya lewati jika sudah ada dan valid
                        # print(f"  Melewati '{species_name}', sudah ada di progres dengan sekuens valid.")
                        continue

                    print(f"  Mengambil sekuens untuk: '{species_name}' ({i+1}/{total_unique_species})...")
                    cleaned_species_name = ' '.join(species_name.split()[:2]) if isinstance(species_name, str) else species_name
                    sequence = self._fetch_ncbi_sequence(cleaned_species_name)
                    
                    species_sequence_map[species_name] = sequence if sequence else "Not Found or Error"
                    if sequence: # Hanya hitung jika sekuens benar-benar didapatkan (bukan None dari _fetch_ncbi_sequence)
                        newly_fetched_count +=1

                    processed_count_this_session += 1
                    
                    if processed_count_this_session % 1 == 0: # Simpan setiap 1 item baru diproses
                        self._save_progress(species_sequence_map, progress_file)
                        # print(f"    Progres disimpan untuk '{species_name}'.")

                    if i < total_unique_species - 1:
                        time.sleep(0.4) 
                
                self._save_progress(species_sequence_map, progress_file)
                print(f"Pengambilan sekuens selesai untuk sesi ini. {newly_fetched_count} sekuens baru diambil.")
                
                df['gen_seq'] = df['scientificName'].map(species_sequence_map)
        else:
            df['gen_seq'] = "Sequence Fetching Skipped"

        print(f"\nFinal cleaned dataset: {len(df)} records")
        return df
    
    def classify_taxonomic_group(self, class_name):
        if pd.isna(class_name): return 'Unknown'
        class_name_str = str(class_name).lower()
        if class_name_str == 'mammalia': return 'Mammals'
        elif class_name_str == 'aves': return 'Birds'
        elif class_name_str == 'amphibia': return 'Amphibians'
        elif class_name_str in ['reptilia', 'crocodylia', 'squamata', 'testudines', 'sphenodontia']: return 'Reptiles'
        elif class_name_str in ['actinopterygii', 'chondrichthyes', 'sarcopterygii', 'agnatha']: return 'Fish'
        elif class_name_str in ['insecta', 'arachnida', 'malacostraca', 'chilopoda', 'diplopoda', 'merostomata']: return 'Arthropods (Non-Insect Invertebrates)'
        elif class_name_str in ['gastropoda', 'bivalvia', 'cephalopoda']: return 'Molluscs'
        else: return 'Other Invertebrates' if class_name_str not in ['magnoliopsida', 'liliopsida', 'pinopsida', 'polypodiopsida'] else 'Plants'

    def create_geodataframe(self, df):
        if df.empty:
            print("Input DataFrame is empty. Returning empty GeoDataFrame.")
            return gpd.GeoDataFrame(columns=['geometry', 'gen_seq'] if 'gen_seq' in df.columns else ['geometry'], 
                                    geometry='geometry', crs='EPSG:4326')
        print("Creating GeoDataFrame...")
        geometry = [Point(lon, lat) for lon, lat in zip(df['decimalLongitude'], df['decimalLatitude'])]
        gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
        print(f"Created GeoDataFrame with {len(gdf)} records")
        return gdf
    
    def load_parks_data(self, shapefile_path):
        print(f"Loading national parks from {shapefile_path}")
        try:
            self.parks_data = gpd.read_file(shapefile_path)
            if self.parks_data.crs != 'EPSG:4326':
                self.parks_data = self.parks_data.to_crs('EPSG:4326')
            print(f"Loaded {len(self.parks_data)} national parks")
            return self.parks_data
        except Exception as e:
            print(f"Error loading parks data: {e}")
            return None
    
    def spatial_join_species_parks(self):
        if self.species_data is None or self.species_data.empty or \
           self.parks_data is None or self.parks_data.empty:
            print("Species data or parks data is missing or empty. Cannot perform spatial join.")
            return gpd.GeoDataFrame(), pd.DataFrame() 
        print("Performing spatial join...")
        species_in_parks = gpd.sjoin(self.species_data, self.parks_data, how='inner', predicate='intersects')
        print(f"Found {len(species_in_parks)} species occurrences within national parks")
        if species_in_parks.empty:
            print("No species occurrences found within any national parks.")
            return species_in_parks, pd.DataFrame()
        park_stats = self.calculate_park_statistics(species_in_parks)
        return species_in_parks, park_stats
    
    def calculate_park_statistics(self, species_in_parks):
        if species_in_parks.empty or self.parks_data is None or self.parks_data.empty:
            print("No species in parks or parks data available for statistics.")
            return pd.DataFrame()
        print("Calculating park statistics...")
        stats = []
        if 'index_right' not in species_in_parks.columns:
            print("Error: 'index_right' column not found in species_in_parks. Cannot calculate park statistics.")
            return pd.DataFrame()
        for park_idx_from_sjoin in species_in_parks['index_right'].unique():
            try:
                park_detail = self.parks_data.loc[park_idx_from_sjoin]
                park_name = park_detail.get('NAME', f'Park_Index_{park_idx_from_sjoin}')
            except KeyError:
                print(f"Warning: Park with index {park_idx_from_sjoin} not found in self.parks_data. Skipping.")
                continue
            park_species_subset = species_in_parks[species_in_parks['index_right'] == park_idx_from_sjoin]
            if not park_species_subset.empty:
                park_stat = {
                    'park_name': park_name,
                    'park_id_original': park_idx_from_sjoin,
                    'total_records': len(park_species_subset),
                    'unique_species': park_species_subset['scientificName'].nunique(),
                    'taxonomic_groups': park_species_subset['taxonomic_group'].value_counts().to_dict(),
                    'families': park_species_subset['family'].nunique() if 'family' in park_species_subset.columns else 0,
                    'genera': park_species_subset['genus'].nunique() if 'genus' in park_species_subset.columns else 0,
                    'threatened_species': len(park_species_subset[
                        park_species_subset['iucnRedListCategory'].isin(['CR', 'EN', 'VU'])
                    ]) if 'iucnRedListCategory' in park_species_subset.columns else 0,
                    'endemic_candidates': self.identify_endemic_candidates(park_species_subset, species_in_parks)
                }
                stats.append(park_stat)
        stats_df = pd.DataFrame(stats)
        if not stats_df.empty: print(f"Calculated statistics for {len(stats_df)} parks with species data")
        else: print("No park statistics could be calculated.")
        return stats_df
    
    def identify_endemic_candidates(self, park_specific_species, all_species_in_all_parks):
        species_in_this_park = set(park_specific_species['scientificName'].unique())
        current_park_id = park_specific_species['index_right'].iloc[0]
        other_parks_species = all_species_in_all_parks[all_species_in_all_parks['index_right'] != current_park_id]
        species_in_other_parks = set()
        if not other_parks_species.empty:
            species_in_other_parks = set(other_parks_species['scientificName'].unique())
        endemic_candidates_set = species_in_this_park - species_in_other_parks
        return len(endemic_candidates_set)

# Usage example
def process_gbif_data(fetch_sequences_flag=False, progress_file=None): # Tambahkan argumen progress_file
    """Example of how to process GBIF data"""
    
    print(f"PENTING: Pastikan Entrez.email diatur ke email Anda yang valid di awal skrip data_processing.py!")
    print(f"Saat ini diatur ke: {Entrez.email}")
    if Entrez.email == "your_actual_email@example.com": # Periksa placeholder email
        print("PERINGATAN: Email Entrez belum diubah dari placeholder. Pengambilan sekuens mungkin gagal atau dibatasi.")
        # fetch_sequences_flag = False # Atau paksa false jika email tidak diatur

    base_data_path = '../data/raw/'
    gbif_files = [
        os.path.join(base_data_path, 'gbif_mammalia_indonesia.csv'),
        os.path.join(base_data_path, 'gbif_aves_indonesia.csv'),
        os.path.join(base_data_path, 'gbif_amphibia_indonesia.csv'),
        os.path.join(base_data_path, 'gbif_reptilia_combined_indonesia.csv'),
    ]
    
    existing_gbif_files = [f for f in gbif_files if os.path.exists(f)]
    if not existing_gbif_files:
        print(f"Error: Tidak ada file data GBIF yang ditemukan di path yang ditentukan. Cek path: {base_data_path}")
        return None
    
    print(f"File GBIF yang akan diproses: {existing_gbif_files}")

    processor = DataProcessor()
    
    # Tentukan path file progres di sini atau gunakan default di dalam kelas
    # Jika Anda ingin mengontrolnya dari sini:
    # script_dir_main = os.path.dirname(os.path.abspath(__file__)) if '__file__' in locals() else os.getcwd()
    # default_progress_file_path = os.path.join(script_dir_main, '..', 'data', 'processed', 'sequence_fetch_progress.json')
    # progress_file_to_use = progress_file if progress_file else default_progress_file_path
    # print(f"File progres akan disimpan/dimuat dari: {progress_file_to_use}")

    species_data = processor.load_and_clean_gbif_data(existing_gbif_files, 
                                                      fetch_sequences=fetch_sequences_flag,
                                                      progress_file_path_override=progress_file) # Teruskan path progres
    
    if species_data is None or species_data.empty:
        print("Gagal memproses data spesies atau data kosong.")
        return None

    output_geojson_path = '../data/processed/cleaned_species_data_with_seq.geojson' if fetch_sequences_flag else '../data/processed/cleaned_species_data.geojson'
    
    output_dir = os.path.dirname(output_geojson_path)
    if output_dir and not os.path.exists(output_dir): # Pastikan output_dir tidak kosong
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    try:
        species_data.to_file(output_geojson_path, driver='GeoJSON')
        print(f"Data berhasil disimpan ke: {output_geojson_path}")
    except Exception as e:
        print(f"Error menyimpan GeoJSON: {e}")
        if 'gen_seq' in species_data.columns:
            print("Mencoba menyimpan lagi tanpa kolom 'gen_seq'...")
            try:
                # Pastikan path untuk fallback juga valid
                fallback_path = output_geojson_path.replace("_with_seq","").replace(".geojson", "_no_seq.geojson")
                species_data.drop(columns=['gen_seq']).to_file(fallback_path, driver='GeoJSON')
                print(f"Data berhasil disimpan (tanpa gen_seq) ke: {fallback_path}")
            except Exception as e2:
                print(f"Gagal menyimpan GeoJSON bahkan tanpa gen_seq: {e2}")
    
    return species_data

if __name__ == "__main__":
    FETCH_SEQUENCES_FROM_NCBI = True 
    
    # Tentukan path file progres di sini jika ingin mengontrol dari luar kelas
    # Jika None, path default di dalam kelas akan digunakan
    # (dengan asumsi __file__ terdefinisi atau skrip dijalankan dari direktori yang benar)
    custom_progress_file = None 
    # Contoh path kustom jika skrip dijalankan dari root proyek:
    # custom_progress_file = os.path.join('data', 'processed', 'sequence_fetch_progress.json')


    print(f"Proses data GBIF dimulai. Pengambilan sekuens dari NCBI diatur ke: {FETCH_SEQUENCES_FROM_NCBI}")
    # Pastikan email Entrez sudah diubah dari placeholder sebelum menjalankan!
    if Entrez.email == "your_actual_email@example.com":
        print("KRUSIAL: Harap ubah 'Entrez.email' di bagian atas skrip dengan email Anda yang valid!")
    
    species_data_result = process_gbif_data(fetch_sequences_flag=FETCH_SEQUENCES_FROM_NCBI, 
                                            progress_file=custom_progress_file)
    
    if species_data_result is not None and not species_data_result.empty:
        print(f"\nProcessing complete! Final dataset has {len(species_data_result)} records.")
        if FETCH_SEQUENCES_FROM_NCBI and 'gen_seq' in species_data_result.columns:
            valid_seq_count = species_data_result[~species_data_result['gen_seq'].isin(["Not Found or Error", "Sequence Fetching Skipped", None])]['gen_seq'].count()
            print(f"  Berhasil mengambil {valid_seq_count} sekuens unik dari NCBI (mungkin termasuk yang dimuat dari progres).")
    else:
        print("\nProcessing ended with no data or an error.")

