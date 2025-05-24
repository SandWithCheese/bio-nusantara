import folium
from folium import plugins
import geopandas as gpd
from shapely.geometry import LineString
import pandas as pd


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
                longitude = point_geom.coords[0][0]

            # Simplified classification based on longitude
            if longitude < 119:  # West of Wallace Line
                return "Western"
            elif longitude < 129:  # Between lines
                return "Central"
            else:  # East of Weber Line
                return "Eastern"

        species_gdf["biogeographic_zone"] = species_gdf["geometry"].apply(get_zone)
        return species_gdf

    def create_biodiversity_map(
        self, species_gdf, title="BioNusantara: Indonesian Biodiversity Map"
    ):
        """Create comprehensive biodiversity visualization using Folium"""

        # Create base map centered on Indonesia
        indonesia_center = [-2.5, 118.0]
        self.map = folium.Map(
            location=indonesia_center,
            zoom_start=5,
            tiles="OpenStreetMap",
            control_scale=True,
        )

        # Add title
        title_html = f"""
        <h3 align="center" style="font-size:20px; color:darkgreen; margin-top:10px;">
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

        # 6. Add layer control
        folium.LayerControl(collapsed=False).add_to(self.map)

        # 7. Add custom legend
        self._add_custom_legend()

        # 8. Add search functionality
        self._add_search_functionality()

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
            ).add_to(self.map)

    def _add_species_points(self, species_gdf):
        """Add species occurrence points with clustering"""
        # Create marker cluster
        marker_cluster = plugins.MarkerCluster(
            name="Species Occurrences",
            show_coverage_on_hover=False,
            options={"maxClusterRadius": 50, "spiderfyDistanceMultiplier": 1.5},
        )

        # Color mapping for taxonomic groups
        color_map = {
            "Mammals": "red",
            "Birds": "blue",
            "Reptiles": "green",
            "Amphibians": "orange",
            "Other": "gray",
        }

        # Add points to cluster
        for idx, species in species_gdf.iterrows():
            if pd.notna(species.geometry):
                lat, lon = species.geometry.y, species.geometry.x

                # Get taxonomic group color
                taxonomic_group = species.get("taxonomic_group", "Other")
                color = color_map.get(taxonomic_group, "gray")

                # Create popup content
                popup_content = self._create_species_popup(species)

                # Create marker
                folium.CircleMarker(
                    location=[lat, lon],
                    radius=6,
                    popup=folium.Popup(popup_content, max_width=400),
                    tooltip=f"{species.get('scientificName', 'Unknown species')}",
                    color="black",
                    weight=1,
                    fillColor=color,
                    fillOpacity=0.7,
                ).add_to(marker_cluster)

        marker_cluster.add_to(self.map)

    def _create_species_popup(self, species):
        """Create detailed popup content for species"""
        popup_html = f"""
        <div style="width: 350px;">
            <h4 style="color: darkgreen; margin-bottom: 10px;">
                <i>{species.get('scientificName', 'Unknown')}</i>
            </h4>
            
            <table style="width: 100%; font-size: 12px;">
                <tr><td><b>Common Name:</b></td><td>{species.get('vernacularName', 'N/A')}</td></tr>
                <tr><td><b>Family:</b></td><td>{species.get('family', 'N/A')}</td></tr>
                <tr><td><b>Order:</b></td><td>{species.get('order', 'N/A')}</td></tr>
                <tr><td><b>Class:</b></td><td>{species.get('class', 'N/A')}</td></tr>
                <tr><td><b>IUCN Status:</b></td><td>{self._format_iucn_status(species.get('iucnRedListCategory'))}</td></tr>
                <tr><td><b>Biogeographic Zone:</b></td><td>{species.get('biogeographic_zone', 'N/A')}</td></tr>
                <tr><td><b>Location:</b></td><td>{species.get('stateProvince', 'N/A')}</td></tr>
                <tr><td><b>Date:</b></td><td>{species.get('eventDate', 'N/A')}</td></tr>
                <tr><td><b>Recorded by:</b></td><td>{species.get('recordedBy', 'N/A')}</td></tr>
                <tr><td><b>Data Source:</b></td><td>{species.get('datasetName', 'N/A')}</td></tr>
            </table>
        </div>
        """
        return popup_html

    def _format_iucn_status(self, status):
        """Format IUCN Red List status with colors"""
        if pd.isna(status) or status == "":
            return "Not Assessed"

        status_colors = {
            "CR": '<span style="color: #d73027; font-weight: bold;">Critically Endangered</span>',
            "EN": '<span style="color: #fc8d59; font-weight: bold;">Endangered</span>',
            "VU": '<span style="color: #fee08b; font-weight: bold;">Vulnerable</span>',
            "NT": '<span style="color: #d9ef8b;">Near Threatened</span>',
            "LC": '<span style="color: #91bfdb;">Least Concern</span>',
            "DD": '<span style="color: #999999;">Data Deficient</span>',
            "NE": '<span style="color: #cccccc;">Not Evaluated</span>',
        }

        return status_colors.get(status, status)

    def _add_species_heatmap(self, species_gdf):
        """Add species density heatmap layer"""
        # Prepare heatmap data
        heat_data = []
        for idx, species in species_gdf.iterrows():
            if pd.notna(species.geometry):
                lat, lon = species.geometry.y, species.geometry.x
                heat_data.append([lat, lon])

        # Create heatmap
        heatmap = plugins.HeatMap(
            heat_data,
            name="Species Density Heatmap",
            min_opacity=0.4,
            max_zoom=18,
            radius=15,
            blur=10,
            gradient={0.4: "blue", 0.6: "cyan", 0.7: "lime", 0.8: "yellow", 1.0: "red"},
            show=False,  # Start hidden
        )

        heatmap.add_to(self.map)

    def _add_biodiversity_statistics(self, species_gdf):
        """Add biodiversity statistics panel"""
        # Calculate statistics
        total_species = species_gdf["scientificName"].nunique()
        total_records = len(species_gdf)

        # Taxonomic group counts
        taxonomic_counts = species_gdf["taxonomic_group"].value_counts().to_dict()

        # Zone distribution
        zone_counts = species_gdf["biogeographic_zone"].value_counts().to_dict()

        # IUCN threatened species
        threatened = (
            species_gdf[species_gdf["iucnRedListCategory"].isin(["CR", "EN", "VU"])]
            if "iucnRedListCategory" in species_gdf.columns
            else pd.DataFrame()
        )

        # Create statistics HTML
        stats_html = f"""
        <div style="position: fixed; 
                    top: 80px; left: 10px; width: 250px; height: auto;
                    background-color: white; border: 2px solid grey; z-index:9999; 
                    font-size:12px; padding: 10px; border-radius: 5px;">
            
            <h4 style="color: darkgreen; margin-top: 0;">📊 Biodiversity Statistics</h4>
            
            <div style="margin-bottom: 10px;">
                <b>🌍 Overall Diversity:</b><br>
                • Total Species: {total_species:,}<br>
                • Total Records: {total_records:,}<br>
                • Threatened Species: {len(threatened)}<br>
            </div>
            
            <div style="margin-bottom: 10px;">
                <b>🏛️ Taxonomic Groups:</b><br>
                {'<br>'.join([f"• {group}: {count}" for group, count in taxonomic_counts.items()])}
            </div>
            
            <div style="margin-bottom: 10px;">
                <b>🗺️ Biogeographic Zones:</b><br>
                {'<br>'.join([f"• {zone}: {count}" for zone, count in zone_counts.items()])}
            </div>
            
            <div style="font-size: 10px; color: #666; margin-top: 10px;">
                Data from GBIF • BioNusantara Project
            </div>
        </div>
        """

        self.map.get_root().html.add_child(folium.Element(stats_html))

    def _add_taxonomic_group_layers(self, species_gdf):
        """Add separate layers for each taxonomic group"""
        taxonomic_groups = species_gdf["taxonomic_group"].unique()

        color_map = {
            "Mammals": "red",
            "Birds": "blue",
            "Reptiles": "green",
            "Amphibians": "orange",
            "Other": "gray",
        }

        for group in taxonomic_groups:
            if pd.notna(group):
                group_data = species_gdf[species_gdf["taxonomic_group"] == group]

                # Create feature group for this taxonomic group
                feature_group = folium.FeatureGroup(
                    name=f"{group} ({len(group_data)} records)",
                    show=False,  # Start hidden
                )

                # Add points for this group
                for idx, species in group_data.iterrows():
                    if pd.notna(species.geometry):
                        lat, lon = species.geometry.y, species.geometry.x

                        folium.CircleMarker(
                            location=[lat, lon],
                            radius=5,
                            popup=folium.Popup(
                                self._create_species_popup(species), max_width=400
                            ),
                            tooltip=f"{species.get('scientificName', 'Unknown')}",
                            color="black",
                            weight=1,
                            fillColor=color_map.get(group, "gray"),
                            fillOpacity=0.8,
                        ).add_to(feature_group)

                feature_group.add_to(self.map)

    def _add_custom_legend(self):
        """Add custom legend to the map"""
        legend_html = """
        <div style="position: fixed; 
                    bottom: 50px; left: 10px; width: 200px; height: auto; 
                    background-color: white; border:2px solid grey; z-index:9999; 
                    font-size:12px; padding: 10px; border-radius: 5px;">
            
            <h4 style="margin-top: 0; color: darkgreen;">🗺️ Map Legend</h4>
            
            <div style="margin-bottom: 8px;">
                <b>Biogeographic Lines:</b><br>
                <span style="color: red;">━━━</span> Wallace Line<br>
                <span style="color: blue;">━━━</span> Weber Line
            </div>
            
            <div style="margin-bottom: 8px;">
                <b>Taxonomic Groups:</b><br>
                <span style="color: red;">●</span> Mammals<br>
                <span style="color: blue;">●</span> Birds<br>
                <span style="color: green;">●</span> Reptiles<br>
                <span style="color: orange;">●</span> Amphibians
            </div>
            
            <div style="font-size: 10px; color: #666;">
                Click on points for details<br>
                Use layer control to toggle views
            </div>
        </div>
        """

        self.map.get_root().html.add_child(folium.Element(legend_html))

    def _add_search_functionality(self):
        """Add search functionality to the map"""
        # Add search plugin (requires folium-plugins)
        try:
            search = plugins.Search(
                layer=None,
                search_label="scientificName",
                placeholder="Search for species...",
                collapsed=True,
            )
            search.add_to(self.map)
        except:
            # If search plugin not available, skip
            pass

    def save_map(self, filename="bionusantara_map.html"):
        """Save the map to HTML file"""
        if self.map:
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
    # Initialize analyzer
    analyzer = BiogeographicAnalyzer()

    # Create biogeographic lines
    analyzer.create_wallace_weber_lines()

    # Load species data
    try:
        species_gdf = gpd.read_file("../data/processed/cleaned_species_data.geojson")
        print(f"Loaded {len(species_gdf)} species records")

        # Classify species by biogeographic zones
        classified_species = analyzer.classify_biogeographic_zones(species_gdf)

        # Create comprehensive biodiversity map
        biodiversity_map = analyzer.create_biodiversity_map(
            classified_species, title="BioNusantara: Indonesian Biodiversity Explorer"
        )

        # Save map
        analyzer.save_map("../output/bionusantara_biodiversity_map.html")

        print("✅ Biodiversity map created successfully!")
        print("📁 Saved as 'bionusantara_biodiversity_map.html'")
        print("🌐 Open the HTML file in your browser to explore the map")

    except FileNotFoundError:
        print("❌ Species data file not found. Please run data processing first.")
        print("💡 Expected file: '../data/processed/cleaned_species_data.geojson'")

    except Exception as e:
        print(f"❌ Error creating map: {e}")
        print("💡 Make sure you have the required dependencies installed:")
