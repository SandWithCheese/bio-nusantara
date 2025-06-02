import folium
from folium import plugins
import geopandas as gpd
from shapely.geometry import LineString
import pandas as pd
import re # Modul untuk regular expression (membersihkan nama file)
import os # Modul untuk path operations (opsional, untuk memeriksa keberadaan file jika diperlukan di sisi Python)

class BiogeographicAnalyzer:
    def __init__(self):
        self.wallace_line = None
        self.weber_line = None
        self.map = None

    def create_wallace_weber_lines(self):
        """Create Wallace and Weber lines from coordinates"""
        # Wallace Line coordinates (simplified)
        wallace_coords = [
            (117.0, -8.5),  # Bali-Lombok
            (119.0, -2.0),  # Sulawesi
            (127.0, 2.0),  # Philippines
        ]

        # Weber Line coordinates (simplified)
        weber_coords = [(129.0, -8.0), (130.0, -2.0), (132.0, 2.0)]

        self.wallace_line = gpd.GeoDataFrame(
            {"name": ["Wallace Line"]},
            geometry=[LineString(wallace_coords)],
            crs="EPSG:4326",
        )

        self.weber_line = gpd.GeoDataFrame(
            {"name": ["Weber Line"]},
            geometry=[LineString(weber_coords)],
            crs="EPSG:4326",
        )

    def classify_biogeographic_zones(self, species_gdf):
        """Classify species points into Western, Central, Eastern zones"""

        def get_zone(point_geom):
            if hasattr(point_geom, "x"):
                longitude = point_geom.x
            else:
                # For LineString or Polygon, take the longitude of the first coordinate of the centroid
                # This might need adjustment if geometries are complex
                if point_geom.geom_type == 'Point':
                    longitude = point_geom.coords[0][0]
                else: # For other types, use centroid
                    longitude = point_geom.centroid.x


            # Simplified classification based on longitude
            if longitude < 119:  # West of Wallace Line
                return "Western"
            elif longitude < 129:  # Between lines
                return "Central"
            else:  # East of Weber Line
                return "Eastern"

        species_gdf["biogeographic_zone"] = species_gdf["geometry"].apply(get_zone)
        return species_gdf

    def _generate_phylo_image_path(self, scientific_name):
        """Generate the expected path for the phylogenetic tree image."""
        if pd.isna(scientific_name) or scientific_name == "":
            return ""
        # Sanitize the name: replace spaces with underscores, remove non-alphanumeric except underscore
        sanitized_name = re.sub(r'[^\w_]', '', scientific_name.replace(' ', '_'))
        image_filename = f"{sanitized_name}_phylo.png"
        # Path relative to the HTML file in the output directory
        return f"phylo_trees/{image_filename}"

    def _create_species_popup(self, species):
        """Create detailed popup content for species, including a button for phylogenetic tree."""
        scientific_name = species.get('scientificName', 'Unknown')
        
        # Basic info part of the popup
        popup_html_parts = [
            f'<div style="width: 350px;">',
            f'<h4 style="color: darkgreen; margin-bottom: 10px;"><i>{scientific_name}</i></h4>',
            '<table style="width: 100%; font-size: 12px;">',
            f"<tr><td><b>Common Name:</b></td><td>{species.get('vernacularName', 'N/A')}</td></tr>",
            f"<tr><td><b>Family:</b></td><td>{species.get('family', 'N/A')}</td></tr>",
            f"<tr><td><b>Order:</b></td><td>{species.get('order', 'N/A')}</td></tr>",
            f"<tr><td><b>Class:</b></td><td>{species.get('class', 'N/A')}</td></tr>",
            f"<tr><td><b>IUCN Status:</b></td><td>{self._format_iucn_status(species.get('iucnRedListCategory'))}</td></tr>",
            f"<tr><td><b>Biogeographic Zone:</b></td><td>{species.get('biogeographic_zone', 'N/A')}</td></tr>",
            f"<tr><td><b>Location:</b></td><td>{species.get('stateProvince', 'N/A')}</td></tr>",
            f"<tr><td><b>Date:</b></td><td>{species.get('eventDate', 'N/A')}</td></tr>",
            f"<tr><td><b>Recorded by:</b></td><td>{species.get('recordedBy', 'N/A')}</td></tr>",
            f"<tr><td><b>Data Source:</b></td><td>{species.get('datasetName', 'N/A')}</td></tr>",
            '</table>'
        ]

        # Add phylogenetic tree button
        phylo_image_path = self._generate_phylo_image_path(scientific_name)
        if scientific_name != 'Unknown' and phylo_image_path:
            # Escape single quotes in scientific_name for JavaScript string
            js_safe_scientific_name = scientific_name.replace("'", "\\'")
            popup_html_parts.append('<hr style="margin-top:10px; margin-bottom:10px;">')
            popup_html_parts.append(
                f'<button onclick="showPhyloModal(\'{js_safe_scientific_name}\', \'{phylo_image_path}\')" '
                'style="background-color: #4CAF50; color: white; padding: 8px 12px; border: none; '
                'border-radius: 4px; cursor: pointer; font-size: 12px;">'
                'View Phylogenetic Tree</button>'
            )
        
        popup_html_parts.append('</div>')
        return "\n".join(popup_html_parts)

    def _format_iucn_status(self, status):
        """Format IUCN Red List status with colors"""
        if pd.isna(status) or status == "":
            return "Not Assessed"

        status_colors = {
            "CR": '<span style="color: #d73027; font-weight: bold;">Critically Endangered</span>',
            "EN": '<span style="color: #fc8d59; font-weight: bold;">Endangered</span>',
            "VU": '<span style="color: #fee08b; font-weight: bold;">Vulnerable</span>', # Adjusted color for better visibility
            "NT": '<span style="color: #9ECAE1;">Near Threatened</span>', # Adjusted color
            "LC": '<span style="color: #4292C6;">Least Concern</span>', # Adjusted color
            "DD": '<span style="color: #999999;">Data Deficient</span>',
            "NE": '<span style="color: #cccccc;">Not Evaluated</span>',
        }
        return status_colors.get(status, status)

    def _add_modal_html_css_js(self):
        """Adds HTML, CSS, and JavaScript for the phylogenetic tree modal to the map."""
        
        modal_html = """
        <div id="phyloModal" class="phylo-modal">
            <div class="phylo-modal-content">
                <span class="phylo-modal-close" onclick="closePhyloModal()">&times;</span>
                <h3 id="phyloModalTitle" style="margin-top: 0; color: darkgreen;">Phylogenetic Tree</h3>
                <img id="phyloModalImage" src="" alt="Phylogenetic Tree" style="width: 100%; height: auto; max-height: 70vh; object-fit: contain; border: 1px solid #ddd; border-radius: 4px; margin-bottom:10px;">
                <p id="phyloModalError" style="color:red; display:none; font-size:12px;"></p>
            </div>
        </div>
        """

        modal_css = """
        <style>
            .phylo-modal {
                display: none; /* Hidden by default */
                position: fixed; /* Stay in place */
                z-index: 10000; /* Sit on top - higher than folium controls */
                left: 0;
                top: 0;
                width: 100%; /* Full width */
                height: 100%; /* Full height */
                overflow: auto; /* Enable scroll if needed */
                background-color: rgba(0,0,0,0.7); /* Black w/ opacity */
                padding-top: 30px; /* Space from top */
            }
            .phylo-modal-content {
                background-color: #fefefe;
                margin: 2% auto; /* Centered */
                padding: 25px;
                border: 1px solid #888;
                width: 85%; /* Responsive width */
                max-width: 800px; /* Max width for large screens */
                border-radius: 10px;
                position: relative;
                box-shadow: 0 5px 15px rgba(0,0,0,0.3);
            }
            .phylo-modal-close {
                color: #aaa;
                float: right;
                font-size: 32px;
                font-weight: bold;
                position: absolute;
                top: 10px;
                right: 20px;
                line-height: 1;
            }
            .phylo-modal-close:hover,
            .phylo-modal-close:focus {
                color: black;
                text-decoration: none;
                cursor: pointer;
            }
        </style>
        """

        modal_js = """
        <script type="text/javascript">
            function showPhyloModal(speciesName, imagePath) {
                var modal = document.getElementById('phyloModal');
                var modalImg = document.getElementById('phyloModalImage');
                var modalTitle = document.getElementById('phyloModalTitle');
                var modalError = document.getElementById('phyloModalError');

                modalTitle.innerHTML = 'Phylogenetic Tree for <i>' + speciesName + '</i>';
                modalImg.src = imagePath;
                modalImg.style.display = 'block'; // Ensure image is visible
                modalError.style.display = 'none'; // Hide error by default

                modalImg.onerror = function() {
                    modalImg.style.display = 'none'; // Hide broken image icon
                    // Fallback to a more generic placeholder if the specific one also fails
                    var placeholderUrl = 'https://placehold.co/600x400/DDDDDD/AAAAAA?text=Tree+Not+Available';
                    modalImg.src = placeholderUrl; // Attempt to load placeholder
                    modalImg.style.display = 'block';

                    modalError.textContent = 'Phylogenetic tree image for ' + speciesName + ' not found at expected path: ' + imagePath + '. Displaying a placeholder.';
                    modalError.style.display = 'block';
                    
                    // Second level onerror to handle placeholder failure
                    modalImg.onerror = function() {
                        modalImg.style.display = 'none';
                        modalError.textContent = 'Phylogenetic tree image and placeholder could not be loaded.';
                        modalError.style.display = 'block';
                        modalImg.onerror = null; // Prevent infinite loop
                    }
                };
                
                modal.style.display = "block";
            }

            function closePhyloModal() {
                var modal = document.getElementById('phyloModal');
                modal.style.display = "none";
                // Clear image src to prevent it from showing briefly next time and free memory
                var modalImg = document.getElementById('phyloModalImage');
                modalImg.src = ""; 
                modalImg.onerror = null; // Reset error handler
                document.getElementById('phyloModalError').style.display = 'none';
            }

            // Close modal if user clicks outside of the modal content
            window.addEventListener('click', function(event) {
                var modal = document.getElementById('phyloModal');
                if (event.target == modal) {
                    closePhyloModal();
                }
            });

            // Close modal if user presses Escape key
            window.addEventListener('keydown', function(event) {
                var modal = document.getElementById('phyloModal');
                if (event.key === 'Escape' && modal.style.display === "block") {
                    closePhyloModal();
                }
            });
        </script>
        """
        self.map.get_root().html.add_child(folium.Element(modal_html))
        self.map.get_root().html.add_child(folium.Element(modal_css))
        self.map.get_root().html.add_child(folium.Element(modal_js))


    def create_biodiversity_map(
        self, species_gdf, title="BioNusantara: Indonesian Biodiversity Map"
    ):
        """Create comprehensive biodiversity visualization using Folium"""

        # Create base map centered on Indonesia
        indonesia_center = [-2.5, 118.0]
        self.map = folium.Map(
            location=indonesia_center,
            zoom_start=5,
            tiles="OpenStreetMap", # Using a more standard tile
            control_scale=True,
        )

        # Add title
        title_html = f"""
        <h3 align="center" style="font-size:20px; color:darkgreen; margin-top:10px; margin-bottom:5px;">
        <b>{title}</b>
        </h3>
        """
        self.map.get_root().html.add_child(folium.Element(title_html))

        # 1. Add biogeographic lines first (background)
        self._add_biogeographic_lines()

        # 2. Add species occurrence points
        self._add_species_points(species_gdf)

        # 3. Add species density heatmap (optional layer)
        self._add_species_heatmap(species_gdf)

        # 4. Add biodiversity statistics
        self._add_biodiversity_statistics(species_gdf)

        # 5. Add taxonomic group layers
        self._add_taxonomic_group_layers(species_gdf)
        
        # 6. Add HTML, CSS, JS for Modal
        self._add_modal_html_css_js() # Crucial step

        # 7. Add layer control (add after all layers are defined)
        folium.LayerControl(collapsed=False).add_to(self.map)

        # 8. Add custom legend
        self._add_custom_legend()

        # 9. Add search functionality
        # self._add_search_functionality(species_gdf) # Pass species_gdf if needed for search layer

        return self.map

    def _add_biogeographic_lines(self):
        """Add Wallace and Weber lines to the map"""
        if self.wallace_line is not None:
            folium.GeoJson(
                self.wallace_line.__geo_interface__,
                style_function=lambda feature: {
                    "color": "red",
                    "weight": 4,
                    "opacity": 0.8,
                    "dashArray": "10, 5",
                },
                tooltip=folium.Tooltip(
                    "Wallace Line - Separates Asian and Transitional fauna"
                ),
                popup=folium.Popup(
                    "<b>Wallace Line</b><br>Biogeographic boundary separating Asian fauna from Australasian species",
                    max_width=300,
                ),
                name="Wallace Line" # Name for LayerControl
            ).add_to(self.map)

        if self.weber_line is not None:
            folium.GeoJson(
                self.weber_line.__geo_interface__,
                style_function=lambda feature: {
                    "color": "blue",
                    "weight": 4,
                    "opacity": 0.8,
                    "dashArray": "10, 5",
                },
                tooltip=folium.Tooltip(
                    "Weber Line - Separates Transitional and Australasian fauna"
                ),
                popup=folium.Popup(
                    "<b>Weber Line</b><br>Biogeographic boundary separating Transitional fauna from Australasian species",
                    max_width=300,
                ),
                name="Weber Line" # Name for LayerControl
            ).add_to(self.map)

    def _add_species_points(self, species_gdf):
        """Add species occurrence points with clustering"""
        marker_cluster = plugins.MarkerCluster(
            name="All Species Occurrences", # Layer name for LayerControl
            show_coverage_on_hover=False,
            options={"maxClusterRadius": 40, "spiderfyDistanceMultiplier": 1.5}, # Adjusted cluster radius
        ).add_to(self.map) # Add to map immediately to be part of LayerControl

        color_map = {
            "Mammals": "#E41A1C", # Red
            "Birds": "#377EB8",   # Blue
            "Reptiles": "#4DAF4A",# Green
            "Amphibians": "#FF7F00", # Orange
            "Other": "#984EA3",   # Purple
            "Unknown": "#AAAAAA" # Grey
        }

        for idx, species in species_gdf.iterrows():
            if pd.notna(species.geometry) and species.geometry.geom_type == 'Point':
                lat, lon = species.geometry.y, species.geometry.x

                taxonomic_group = species.get("taxonomic_group", "Unknown")
                color = color_map.get(taxonomic_group, color_map["Unknown"])

                popup_content = self._create_species_popup(species)

                folium.CircleMarker(
                    location=[lat, lon],
                    radius=5, # Slightly smaller for less clutter
                    popup=folium.Popup(popup_content, max_width=400, sticky=False), # Non-sticky popup
                    tooltip=f"{species.get('scientificName', 'Unknown species')} ({taxonomic_group})",
                    color="black",
                    weight=0.5, # Thinner border
                    fillColor=color,
                    fillOpacity=0.8,
                ).add_to(marker_cluster)
        # marker_cluster.add_to(self.map) # Already added when created

    def _add_species_heatmap(self, species_gdf):
        """Add species density heatmap layer"""
        heat_data = []
        for idx, species in species_gdf.iterrows():
            if pd.notna(species.geometry) and species.geometry.geom_type == 'Point':
                lat, lon = species.geometry.y, species.geometry.x
                heat_data.append([lat, lon])
        
        if heat_data: # Only add heatmap if there's data
            plugins.HeatMap(
                heat_data,
                name="Species Density Heatmap",
                min_opacity=0.3,
                max_zoom=18, # Show up to zoom level 18
                radius=12, # Adjusted radius
                blur=8,   # Adjusted blur
                gradient={0.2: 'blue', 0.4: 'cyan', 0.6: 'lime', 0.8: 'yellow', 1.0: 'red'},
                show=False,  # Start hidden
            ).add_to(self.map)

    def _add_biodiversity_statistics(self, species_gdf):
        """Add biodiversity statistics panel"""
        total_species = species_gdf["scientificName"].nunique()
        total_records = len(species_gdf)

        taxonomic_counts = species_gdf["taxonomic_group"].value_counts().to_dict()
        zone_counts = species_gdf["biogeographic_zone"].value_counts().to_dict()

        threatened_count = 0
        if "iucnRedListCategory" in species_gdf.columns:
            threatened_count = len(
                species_gdf[species_gdf["iucnRedListCategory"].isin(["CR", "EN", "VU"])]
            )

        stats_html = f"""
        <div style="position: fixed; 
                    top: 70px; left: 10px; width: 260px; max-height: calc(100vh - 140px); 
                    overflow-y: auto;
                    background-color: rgba(255,255,255,0.9); border: 1px solid #bbb; z-index:9998; 
                    font-size:12px; padding: 12px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1);">
            
            <h4 style="color: darkgreen; margin-top: 0; margin-bottom: 10px; font-size: 14px; border-bottom: 1px solid #eee; padding-bottom: 5px;">📊 Biodiversity Statistics</h4>
            
            <div style="margin-bottom: 8px;">
                <b>🌍 Overall Diversity:</b><br>
                &bull; Total Species: {total_species:,}<br>
                &bull; Total Records: {total_records:,}<br>
                &bull; Threatened (CR,EN,VU): {threatened_count:,}<br>
            </div>
            
            <div style="margin-bottom: 8px;">
                <b>🏛️ Taxonomic Groups:</b><br>
                {'<br>'.join([f"&bull; {group}: {count:,}" for group, count in taxonomic_counts.items()])}
            </div>
            
            <div>
                <b>🗺️ Biogeographic Zones:</b><br>
                {'<br>'.join([f"&bull; {zone}: {count:,}" for zone, count in zone_counts.items()])}
            </div>
            
            <div style="font-size: 10px; color: #777; margin-top: 12px; border-top: 1px solid #eee; padding-top: 8px;">
                Data from GBIF &bull; BioNusantara Project
            </div>
        </div>
        """
        self.map.get_root().html.add_child(folium.Element(stats_html))

    def _add_taxonomic_group_layers(self, species_gdf):
        """Add separate layers for each taxonomic group"""
        taxonomic_groups = sorted([tg for tg in species_gdf["taxonomic_group"].unique() if pd.notna(tg)])

        color_map = {
            "Mammals": "#E41A1C", "Birds": "#377EB8", "Reptiles": "#4DAF4A",
            "Amphibians": "#FF7F00", "Other": "#984EA3", "Unknown": "#AAAAAA"
        }

        for group in taxonomic_groups:
            group_data = species_gdf[species_gdf["taxonomic_group"] == group]
            if group_data.empty:
                continue

            feature_group = folium.FeatureGroup(
                name=f"{group} ({len(group_data):,} records)",
                show=False,  # Start hidden
            )

            for idx, species in group_data.iterrows():
                if pd.notna(species.geometry) and species.geometry.geom_type == 'Point':
                    lat, lon = species.geometry.y, species.geometry.x
                    folium.CircleMarker(
                        location=[lat, lon],
                        radius=5,
                        popup=folium.Popup(self._create_species_popup(species), max_width=400),
                        tooltip=f"{species.get('scientificName', 'Unknown')} ({group})",
                        color="black", weight=0.5,
                        fillColor=color_map.get(group, color_map["Unknown"]),
                        fillOpacity=0.8,
                    ).add_to(feature_group)
            feature_group.add_to(self.map)

    def _add_custom_legend(self):
        """Add custom legend to the map"""
        legend_html = """
        <div style="position: fixed; 
                    bottom: 20px; left: 10px; width: 220px; 
                    background-color: rgba(255,255,255,0.9); border:1px solid #bbb; z-index:9998; 
                    font-size:11px; padding: 10px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1);">
            
            <h4 style="margin-top: 0; margin-bottom:8px; color: darkgreen; font-size: 13px; border-bottom: 1px solid #eee; padding-bottom: 4px;">🗺️ Map Legend</h4>
            
            <div style="margin-bottom: 6px;">
                <b>Biogeographic Lines:</b><br>
                <svg width="20" height="10" style="vertical-align: middle;"><line x1="0" y1="5" x2="20" y2="5" style="stroke:red; stroke-width:3; stroke-dasharray:4,2;"/></svg> Wallace Line<br>
                <svg width="20" height="10" style="vertical-align: middle;"><line x1="0" y1="5" x2="20" y2="5" style="stroke:blue; stroke-width:3; stroke-dasharray:4,2;"/></svg> Weber Line
            </div>
            
            <div>
                <b>Taxonomic Groups (Example):</b><br>
                <span style="color: #E41A1C; font-size: 18px; vertical-align: middle;">●</span> Mammals<br>
                <span style="color: #377EB8; font-size: 18px; vertical-align: middle;">●</span> Birds<br>
                <span style="color: #4DAF4A; font-size: 18px; vertical-align: middle;">●</span> Reptiles<br>
                <span style="color: #FF7F00; font-size: 18px; vertical-align: middle;">●</span> Amphibians
            </div>
            <div style="font-size: 9px; color: #777; margin-top: 8px; border-top: 1px solid #eee; padding-top: 6px;">
                Toggle layers via control panel.<br>Click points for species details & tree.
            </div>
        </div>
        """
        self.map.get_root().html.add_child(folium.Element(legend_html))

    def _add_search_functionality(self, species_gdf_for_search):
        """Add search functionality to the map.
           This requires creating a GeoJson layer from points for the search plugin.
        """
        # Prepare GeoJSON data for search
        features = []
        for _, row in species_gdf_for_search.iterrows():
            if pd.notna(row.geometry) and row.geometry.geom_type == 'Point':
                features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "Point",
                        "coordinates": [row.geometry.x, row.geometry.y]
                    },
                    "properties": {
                        "scientificName": row.get('scientificName', 'N/A'),
                        "popup": self._create_species_popup(row) # Search results can also have popups
                    }
                })
        
        search_geojson = {"type": "FeatureCollection", "features": features}

        search_layer = folium.GeoJson(
            search_geojson,
            name="Searchable Species Layer",
            popup=folium.GeoJsonPopup(fields=["scientificName"], aliases=["Species:"]),
            show=False # Typically, this layer is not shown directly, only used by search
        ).add_to(self.map)

        try:
            # Search plugin for specific layer
            plugins.Search(
                layer=search_layer, # Search within this GeoJSON layer
                search_label="scientificName", # Property to search in GeoJSON
                placeholder="Search species by scientific name...",
                collapsed=True, # Start collapsed
                position='topright'
            ).add_to(self.map)
        except Exception as e:
            print(f"Could not add Search plugin: {e}. Ensure folium.plugins.Search is available.")


    def save_map(self, filename="bionusantara_map.html"):
        """Save the map to HTML file"""
        if self.map:
            # Ensure output directory exists
            output_dir = os.path.dirname(filename)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
                
            self.map.save(filename)
            print(f"Map saved as {filename}")
        else:
            print("No map to save. Create map first using create_biodiversity_map()")

    def display_map(self):
        """Display the map (for Jupyter notebooks)"""
        if self.map:
            return self.map
        else:
            print("No map to display. Create map first using create_biodiversity_map()")


# Usage example
if __name__ == "__main__":
    analyzer = BiogeographicAnalyzer()
    analyzer.create_wallace_weber_lines()

    # Define paths (relative to the script's location)
    script_dir = os.path.dirname(__file__) # Gets the directory where the script is located
    species_data_path = os.path.join(script_dir, "../data/processed/cleaned_species_data.geojson")
    output_map_path = os.path.join(script_dir, "../output/bionusantara_biodiversity_map.html")
    
    # Create dummy phylo_trees directory if it doesn't exist for testing
    phylo_tree_dir = os.path.join(script_dir, "../output/phylo_trees")
    if not os.path.exists(phylo_tree_dir):
        os.makedirs(phylo_tree_dir)
        print(f"Created directory: {phylo_tree_dir} (for placeholder tree images)")
        # Optional: Create a dummy image for testing if needed
        # try:
        #     from PIL import Image, ImageDraw, ImageFont
        #     img = Image.new('RGB', (400, 300), color = (220, 220, 220))
        #     d = ImageDraw.Draw(img)
        #     try:
        #         font = ImageFont.truetype("arial.ttf", 20)
        #     except IOError:
        #         font = ImageFont.load_default()
        #     d.text((10,10), "Placeholder Phylo Tree", fill=(50,50,50), font=font)
        #     # Save a common placeholder that might be picked up by a common species name
        #     # For example, if your data has 'Homo sapiens'
        #     # img.save(os.path.join(phylo_tree_dir, "Homo_sapiens_phylo.png"))
        #     print("Created a dummy placeholder image for testing (if a common species name is used).")
        # except ImportError:
        #     print("Pillow not installed, cannot create dummy image.")
        # except Exception as e_img:
        #     print(f"Error creating dummy image: {e_img}")


    try:
        species_gdf = gpd.read_file(species_data_path)
        print(f"Loaded {len(species_gdf)} species records from {species_data_path}")

        # Ensure 'taxonomic_group' exists, if not, create a default one
        if 'taxonomic_group' not in species_gdf.columns:
            print("'taxonomic_group' column not found. Adding a default 'Unknown' group.")
            species_gdf['taxonomic_group'] = 'Unknown'
            # You might want to call a classification method here if available from DataProcessor
            # For example:
            # from data_processing import DataProcessor
            # processor = DataProcessor()
            # species_gdf['taxonomic_group'] = species_gdf['class'].apply(processor.classify_taxonomic_group)


        classified_species = analyzer.classify_biogeographic_zones(species_gdf.copy()) # Use .copy() to avoid SettingWithCopyWarning

        biodiversity_map = analyzer.create_biodiversity_map(
            classified_species, title="BioNusantara: Indonesian Biodiversity Explorer"
        )
        
        # Add search functionality (optional, can be resource-intensive for large datasets)
        # analyzer._add_search_functionality(classified_species)


        analyzer.save_map(output_map_path)

        print("✅ Biodiversity map created successfully!")
        print(f"📁 Saved as '{output_map_path}'")
        print("🌐 Open the HTML file in your browser to explore the map.")
        print("❗ Remember to place your phylogenetic tree images in '../output/phylo_trees/'")
        print("   with filenames like 'Genus_species_phylo.png' (e.g., 'Panthera_tigris_phylo.png').")

    except FileNotFoundError:
        print(f"❌ Species data file not found at: {species_data_path}")
        print("💡 Please ensure 'cleaned_species_data.geojson' exists in the correct directory.")
        print("   You might need to run the data processing script first.")

    except Exception as e:
        import traceback
        print(f"❌ An error occurred while creating the map: {e}")
        print(traceback.format_exc())
        print("💡 Check data integrity, file paths, and dependencies.")

