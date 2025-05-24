import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np

class DataProcessor:
    def __init__(self):
        self.species_data = None
        self.parks_data = None
        
    def load_and_clean_gbif_data(self, file_paths):
        """Load and clean GBIF occurrence data from multiple CSV files"""
        all_data = []
        
        for file_path in file_paths:
            print(f"Loading {file_path}...")
            try:
                df = pd.read_csv(file_path)
                print(f"Loaded {len(df)} records from {file_path}")
                all_data.append(df)
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
        
        if not all_data:
            raise ValueError("No data files could be loaded")
        
        # Combine all datasets
        combined_df = pd.concat(all_data, ignore_index=True)
        print(f"Combined dataset: {len(combined_df)} total records")
        
        # Clean the data based on GBIF format
        cleaned_df = self.clean_gbif_data(combined_df)
        
        # Convert to GeoDataFrame
        self.species_data = self.create_geodataframe(cleaned_df)
        
        return self.species_data
    
    def clean_gbif_data(self, df):
        """Clean GBIF data based on the actual format"""
        print("Cleaning GBIF data...")
        
        # Remove records without coordinates
        initial_count = len(df)
        df = df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
        print(f"Removed {initial_count - len(df)} records without coordinates")
        
        # Filter for high-quality coordinates (< 1000m uncertainty)
        if 'coordinateUncertaintyInMeters' in df.columns:
            high_quality = df['coordinateUncertaintyInMeters'].fillna(0) < 1000
            df = df[high_quality]
            print(f"Kept {len(df)} records with coordinate uncertainty < 1000m")
        
        # Remove invalid coordinates
        valid_coords = (
            (df['decimalLatitude'].between(-90, 90)) & 
            (df['decimalLongitude'].between(-180, 180))
        )
        df = df[valid_coords]
        print(f"Kept {len(df)} records with valid coordinates")
        
        # Filter for Indonesia (rough bounding box)
        indonesia_bounds = (
            (df['decimalLatitude'].between(-11, 6)) &  # Indonesia latitude range
            (df['decimalLongitude'].between(95, 141))   # Indonesia longitude range
        )
        df = df[indonesia_bounds]
        print(f"Kept {len(df)} records within Indonesia bounds")
        
        # Keep only essential columns for analysis
        essential_columns = [
            'key', 'scientificName', 'acceptedScientificName', 'taxonKey',
            'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
            'decimalLatitude', 'decimalLongitude', 'coordinateUncertaintyInMeters',
            'year', 'month', 'day', 'eventDate',
            'stateProvince', 'locality', 'country', 'countryCode',
            'basisOfRecord', 'occurrenceStatus', 'iucnRedListCategory',
            'datasetName', 'recordedBy', 'identifiedBy'
        ]
        
        # Keep only columns that exist in the dataset
        available_columns = [col for col in essential_columns if col in df.columns]
        df = df[available_columns].copy()
        
        # Add taxonomic group identifier for easier analysis
        df['taxonomic_group'] = df['class'].apply(self.classify_taxonomic_group)
        
        print(f"Final cleaned dataset: {len(df)} records")
        return df
    
    def classify_taxonomic_group(self, class_name):
        """Classify records into major taxonomic groups"""
        if pd.isna(class_name):
            return 'Unknown'
        
        class_name = str(class_name).lower()
        
        if class_name == 'mammalia':
            return 'Mammals'
        elif class_name == 'aves':
            return 'Birds'
        elif class_name == 'amphibia':
            return 'Amphibians'
        elif class_name in ['reptilia', 'crocodylia', 'squamata', 'testudines', 'sphenodontia']:
            return 'Reptiles'
        else:
            return 'Other'
    
    def create_geodataframe(self, df):
        """Convert pandas DataFrame to GeoDataFrame with Point geometries"""
        print("Creating GeoDataFrame...")
        
        # Create Point geometries
        geometry = [
            Point(lon, lat) for lon, lat in 
            zip(df['decimalLongitude'], df['decimalLatitude'])
        ]
        
        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
        
        print(f"Created GeoDataFrame with {len(gdf)} records")
        return gdf
    
    def load_parks_data(self, shapefile_path):
        """Load national parks shapefile"""
        print(f"Loading national parks from {shapefile_path}")
        try:
            self.parks_data = gpd.read_file(shapefile_path)
            # Ensure CRS is WGS84
            if self.parks_data.crs != 'EPSG:4326':
                self.parks_data = self.parks_data.to_crs('EPSG:4326')
            print(f"Loaded {len(self.parks_data)} national parks")
            return self.parks_data
        except Exception as e:
            print(f"Error loading parks data: {e}")
            return None
    
    def spatial_join_species_parks(self):
        """Associate species occurrences with national parks"""
        if self.species_data is None or self.parks_data is None:
            raise ValueError("Load species and parks data first")
        
        print("Performing spatial join...")
        
        # Spatial join to find species within parks
        species_in_parks = gpd.sjoin(
            self.species_data, 
            self.parks_data, 
            how='inner', 
            predicate='intersects'
        )
        
        print(f"Found {len(species_in_parks)} species occurrences within national parks")
        
        # Calculate statistics per park
        park_stats = self.calculate_park_statistics(species_in_parks)
        
        return species_in_parks, park_stats
    
    def calculate_park_statistics(self, species_in_parks):
        """Calculate biodiversity statistics for each park"""
        print("Calculating park statistics...")
        
        stats = []
        
        for park_idx, park in self.parks_data.iterrows():
            park_species = species_in_parks[
                species_in_parks.index_right == park_idx
            ]
            
            if len(park_species) > 0:
                park_stat = {
                    'park_name': park.get('NAME', f'Park_{park_idx}'),
                    'park_id': park_idx,
                    'total_records': len(park_species),
                    'unique_species': park_species['scientificName'].nunique(),
                    'taxonomic_groups': park_species['taxonomic_group'].value_counts().to_dict(),
                    'families': park_species['family'].nunique() if 'family' in park_species.columns else 0,
                    'genera': park_species['genus'].nunique() if 'genus' in park_species.columns else 0,
                    'threatened_species': len(park_species[
                        park_species['iucnRedListCategory'].isin(['CR', 'EN', 'VU'])
                    ]) if 'iucnRedListCategory' in park_species.columns else 0,
                    'endemic_candidates': self.identify_endemic_candidates(park_species)
                }
                stats.append(park_stat)
        
        stats_df = pd.DataFrame(stats)
        print(f"Calculated statistics for {len(stats_df)} parks with species data")
        
        return stats_df
    
    def identify_endemic_candidates(self, park_species):
        """Identify potential endemic species (species found only in this park)"""
        # This is a simplified approach - in reality, you'd need more comprehensive data
        species_counts = park_species['scientificName'].value_counts()
        # Species with only 1-2 records might be endemic candidates
        rare_species = species_counts[species_counts <= 2]
        return len(rare_species)

# Usage example
def process_gbif_data():
    """Example of how to process GBIF data"""
    
    # File paths for downloaded GBIF data
    gbif_files = [
        '../data/raw/gbif_mammalia_indonesia.csv',
        '../data/raw/gbif_aves_indonesia.csv',
        '../data/raw/gbif_amphibia_indonesia.csv',
        '../data/raw/gbif_crocodylia_indonesia.csv',
        '../data/raw/gbif_squamata_indonesia.csv',
        '../data/raw/gbif_testudines_indonesia.csv',
        # '../data/raw/gbif_sphenodontia_indonesia.csv'
    ]
    
    # Initialize processor
    processor = DataProcessor()
    
    # Load and clean data
    species_data = processor.load_and_clean_gbif_data(gbif_files)
    
    # Save cleaned data
    species_data.to_file('../data/processed/cleaned_species_data.geojson', driver='GeoJSON')
    
    # Load parks data (you'll need to get Indonesian national parks shapefile)
    # parks_data = processor.load_parks_data('../data/shapefiles/indonesia_national_parks.shp')
    
    # Perform spatial join (when parks data is available)
    # species_in_parks, park_stats = processor.spatial_join_species_parks()
    
    return species_data

if __name__ == "__main__":
    species_data = process_gbif_data()
    print(f"Processing complete! Final dataset has {len(species_data)} records")