import folium
from folium import plugins
import geopandas as gpd
from shapely.geometry import LineString
import pandas as pd
import re 
import os
import json 
import math 

class BiogeographicAnalyzer:
    def __init__(self):
        self.wallace_line = None
        self.weber_line = None
        self.map = None

    def _standardize_species_name_for_phylo_file(self, original_name):
        """
        Convert nama spesies menjadi format "Genus_spesies".
        """
        if pd.isna(original_name) or not str(original_name).strip():
            return "unknown_species_id" 

        name = str(original_name)
        name = re.sub(r'\s*\(.*\)\s*', '', name).strip()
        name = re.sub(r'[, ]*\b\d{4}\b$', '', name).strip()
        name = re.sub(r'\s*\b(spp?|cf|aff)\.?\b.*', '', name, flags=re.IGNORECASE).strip()
        name = re.sub(r'[^a-zA-Z0-9\s]+$', '', name).strip()
        
        name_parts = name.split()
        
        if len(name_parts) >= 2:
            genus = name_parts[0].capitalize()
            species = name_parts[1].lower()
            base_name = f"{genus}_{species}"
            return re.sub(r'[^\w-]', '', base_name) 
        elif len(name_parts) == 1:
            base_name = name_parts[0].capitalize()
            return re.sub(r'[^\w-]', '', base_name)
        else:
            return "unknown_species_id"

    def create_wallace_weber_lines(self):
        # Plot wallace
        wallace_coords = [
                (114.770503, -11.054354),
                (115.880123, -8.390789),
                (116.517330, -6.298842),  
                (119.329829, 0.510935),  
                (119.505611, 1.576416),   
                (119.989009, 2.334026 ),  
                (123.481628, 4.309107),   
                (128.230147, 5.248134 )
        ]

        # Plot weber
        weber_coords = [
            (131.651688, 3.938487),
            (127.974819, 3.127062),
            (126.201054, 1.317781),
            (126.127147  , -0.418912),
            (126.810786, -1.675094),
            (126.496681, -2.727476),
            (125.406554, -3.151880),
            (125.498938, -4.110730),
            (128.178063, -5.804167),
            (130.450701, -7.089332),
            (130.284410, -8.407480),
            (126.755355, -9.520776),
            (123.281731, -11.482735)
        ]
        self.wallace_line = gpd.GeoDataFrame(
            {"name": ["Wallace Line (Detailed)"]}, geometry=[LineString(wallace_coords)], crs="EPSG:4326",
        )
        self.weber_line = gpd.GeoDataFrame(
            {"name": ["Weber Line (Detailed)"]}, geometry=[LineString(weber_coords)], crs="EPSG:4326",
        )

    def classify_biogeographic_zones(self, species_gdf):
        def get_zone(point_geom):
            if hasattr(point_geom, "x"): longitude = point_geom.x
            else: longitude = point_geom.centroid.x 
            if longitude < 119: return "Western (Asian)" 
            elif longitude < 129: return "Central (Wallacea)"
            else: return "Eastern (Australasian)"
        species_gdf["biogeographic_zone"] = species_gdf["geometry"].apply(get_zone)
        return species_gdf

    def _generate_phylo_image_path(self, scientific_name_original):
        """
        Menghasilkan path ke file gambar filogenetik berdasarkan nama ilmiah .
        """
        if pd.isna(scientific_name_original) or not str(scientific_name_original).strip():
            base_name = self._standardize_species_name_for_phylo_file(scientific_name_original) 
        else:
            base_name = self._standardize_species_name_for_phylo_file(scientific_name_original)
        
        if not base_name or base_name == "unknown_species_id": 
            return "" 

        image_filename = f"{base_name}_phylo.png"
        return f"phylo_trees/{image_filename}"


    def _create_species_popup(self, species):
        """
        Modal Informasi Spesies
        """
        scientific_name = species.get('scientificName', 'Unknown')
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
        phylo_image_path = self._generate_phylo_image_path(scientific_name) 
        if scientific_name != 'Unknown' and phylo_image_path: 
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
        if pd.isna(status) or status == "": return "Not Assessed"
        status_colors = {
            "CR": '<span style="color: #d73027; font-weight: bold;">Critically Endangered</span>',
            "EN": '<span style="color: #fc8d59; font-weight: bold;">Endangered</span>',
            "VU": '<span style="color: #fee08b; font-weight: bold;">Vulnerable</span>',
            "NT": '<span style="color: #9ECAE1;">Near Threatened</span>',
            "LC": '<span style="color: #4292C6;">Least Concern</span>',
            "DD": '<span style="color: #999999;">Data Deficient</span>',
            "NE": '<span style="color: #cccccc;">Not Evaluated</span>',
        }
        return status_colors.get(status, status)

    def _add_modal_html_css_js(self):
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
            .phylo-modal { display: none; position: fixed; z-index: 10000; left: 0; top: 0; width: 100%; height: 100%; overflow: auto; background-color: rgba(0,0,0,0.7); padding-top: 30px; }
            .phylo-modal-content { background-color: #fefefe; margin: 2% auto; padding: 25px; border: 1px solid #888; width: 85%; max-width: 800px; border-radius: 10px; position: relative; box-shadow: 0 5px 15px rgba(0,0,0,0.3); }
            .phylo-modal-close { color: #aaa; float: right; font-size: 32px; font-weight: bold; position: absolute; top: 10px; right: 20px; line-height: 1; }
            .phylo-modal-close:hover, .phylo-modal-close:focus { color: black; text-decoration: none; cursor: pointer; }
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
                modalImg.src = imagePath; modalImg.style.display = 'block'; modalError.style.display = 'none';
                modalImg.onerror = function() {
                    modalImg.style.display = 'none'; 
                    var placeholderUrl = 'https://placehold.co/600x400/DDDDDD/AAAAAA?text=Tree+Not+Available';
                    modalImg.src = placeholderUrl; modalImg.style.display = 'block';
                    modalError.textContent = 'Phylogenetic tree image for ' + speciesName + ' not found: ' + imagePath + '. Placeholder shown.';
                    modalError.style.display = 'block';
                    modalImg.onerror = function() {
                        modalImg.style.display = 'none';
                        modalError.textContent = 'Phylogenetic tree image and placeholder could not be loaded.';
                        modalError.style.display = 'block'; modalImg.onerror = null; 
                    }
                };
                modal.style.display = "block";
            }
            function closePhyloModal() {
                var modal = document.getElementById('phyloModal'); modal.style.display = "none";
                var modalImg = document.getElementById('phyloModalImage'); modalImg.src = ""; 
                modalImg.onerror = null; document.getElementById('phyloModalError').style.display = 'none';
            }
            window.addEventListener('click', function(event) {
                var modal = document.getElementById('phyloModal'); if (event.target == modal) { closePhyloModal(); }
            });
            window.addEventListener('keydown', function(event) {
                var modal = document.getElementById('phyloModal');
                if (event.key === 'Escape' && modal.style.display === "block") { closePhyloModal(); }
            });
        </script>
        """
        self.map.get_root().html.add_child(folium.Element(modal_html))
        self.map.get_root().html.add_child(folium.Element(modal_css))
        self.map.get_root().html.add_child(folium.Element(modal_js))

    def _add_pulsing_dot_css(self):
        pulsing_dot_css = """
        <style>
            .pulsing-dot-container {}
            .pulsing-dot { width: 24px; height: 24px; background-color: #28a745; border-radius: 50%; border: 3px solid white; box-shadow: 0 0 5px rgba(0,0,0,0.6); position: relative; }
            .pulsing-dot::before { content: ''; position: absolute; display: block; width: 200%; height: 200%; top: -50%; left: -50%; background-color: #28a745; border-radius: 50%; opacity: 0.7; animation: pulse-animation 2s infinite; z-index: -1; }
            @keyframes pulse-animation { 0% { transform: scale(0.5); opacity: 0.7; } 50% { opacity: 0.2; } 100% { transform: scale(1.5); opacity: 0; } }
        </style>
        """
        self.map.get_root().html.add_child(folium.Element(pulsing_dot_css))

    def _add_national_parks_layer(self, national_parks_data_path):
        try:
            with open(national_parks_data_path, 'r', encoding='utf-8') as f:
                parks_data = json.load(f)
        except FileNotFoundError:
            print(f"Peringatan: File data taman nasional tidak ditemukan di {national_parks_data_path}")
            return
        except json.JSONDecodeError:
            print(f"Peringatan: File data taman nasional {national_parks_data_path} bukan JSON valid.")
            return
        except Exception as e:
            print(f"Error memuat data taman nasional: {e}")
            return
        if not parks_data:
            print("Tidak ada data taman nasional untuk ditampilkan.")
            return
        parks_fg = folium.FeatureGroup(name="Taman Nasional (Titik Berdenyut)", show=True)
        for park in parks_data:
            name = park.get("name", "Nama Tidak Diketahui")
            lat = park.get("latitude")
            lon = park.get("longitude")
            size = park.get("size", 0) 
            if lat is None or lon is None:
                print(f"Data koordinat tidak lengkap untuk taman: {name}, dilewati.")
                continue
            tooltip_text = f"<b>{name}</b><br>Luas: {size:,.2f} km²"
            icon_html = '<div class="pulsing-dot"></div>'
            folium.Marker(
                location=[lat, lon], tooltip=tooltip_text, popup=tooltip_text,
                icon=folium.DivIcon(html=icon_html, icon_size=(24,24), icon_anchor=(12,12))
            ).add_to(parks_fg)
        parks_fg.add_to(self.map)
        print(f"Lapisan Taman Nasional ditambahkan dengan {len(parks_data)} taman (efek berdenyut).")

    def create_biodiversity_map(
        self, species_gdf, title="BioNusantara: Indonesian Biodiversity Map"
    ):
        indonesia_center = [-2.5, 118.0]
        self.map = folium.Map(
            location=indonesia_center, zoom_start=5, tiles="OpenStreetMap", control_scale=True,
        )
        title_html = f"""<h3 align="center" style="font-size:20px; color:darkgreen; margin-top:10px; margin-bottom:5px;"><b>{title}</b></h3>"""
        self.map.get_root().html.add_child(folium.Element(title_html))
        self._add_biogeographic_lines()
        self._add_species_points(species_gdf)
        script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in locals() else os.getcwd()
        national_parks_file = os.path.join(script_dir, "..", "data", "processed", "national_parks.json")
        self._add_national_parks_layer(national_parks_file)
        self._add_species_heatmap(species_gdf)
        self._add_biodiversity_statistics(species_gdf)
        self._add_taxonomic_group_layers(species_gdf)
        self._add_modal_html_css_js()
        self._add_pulsing_dot_css() 
        folium.LayerControl(collapsed=False).add_to(self.map)
        self._add_custom_legend()
        return self.map

    def _add_biogeographic_lines(self):
        if self.wallace_line is not None:
            folium.GeoJson(
                self.wallace_line, 
                style_function=lambda feature: {"color": "red", "weight": 4, "opacity": 0.8, "dashArray": "10, 5"},
                tooltip=folium.Tooltip("Wallace Line - Separates Asian and Transitional fauna"),
                popup=folium.Popup("<b>Wallace Line</b><br>Biogeographic boundary separating Asian fauna from Australasian species", max_width=300),
                name="Wallace Line (Detailed)" 
            ).add_to(self.map)
        if self.weber_line is not None:
            folium.GeoJson(
                self.weber_line, 
                style_function=lambda feature: {"color": "blue", "weight": 4, "opacity": 0.8, "dashArray": "10, 5"},
                tooltip=folium.Tooltip("Weber Line - Separates Transitional and Australasian fauna"),
                popup=folium.Popup("<b>Weber Line</b><br>Biogeographic boundary separating Transitional fauna from Australasian species", max_width=300),
                name="Weber Line (Detailed)" 
            ).add_to(self.map)

    def _add_species_points(self, species_gdf):
        marker_cluster = plugins.MarkerCluster(name="All Species Occurrences", show_coverage_on_hover=False, options={"maxClusterRadius": 40, "spiderfyDistanceMultiplier": 1.5}).add_to(self.map)
        color_map = {"Mammals": "#E41A1C", "Birds": "#377EB8", "Reptiles": "#4DAF4A", "Amphibians": "#FF7F00", "Other": "#984EA3", "Unknown": "#AAAAAA"}
        for idx, species in species_gdf.iterrows():
            if pd.notna(species.geometry) and species.geometry.geom_type == 'Point':
                lat, lon = species.geometry.y, species.geometry.x
                taxonomic_group = species.get("taxonomic_group", "Unknown")
                color = color_map.get(taxonomic_group, color_map["Unknown"])
                popup_content = self._create_species_popup(species)
                folium.CircleMarker(
                    location=[lat, lon], radius=5, popup=folium.Popup(popup_content, max_width=400, sticky=False),
                    tooltip=f"{species.get('scientificName', 'Unknown species')} ({taxonomic_group})",
                    color="black", weight=0.5, fillColor=color, fillOpacity=0.8,
                ).add_to(marker_cluster)

    def _add_species_heatmap(self, species_gdf):
        heat_data = []
        for idx, species in species_gdf.iterrows():
            if pd.notna(species.geometry) and species.geometry.geom_type == 'Point':
                lat, lon = species.geometry.y, species.geometry.x
                heat_data.append([lat, lon])
        if heat_data:
            plugins.HeatMap(
                heat_data, name="Species Density Heatmap", min_opacity=0.3, max_zoom=18, radius=12, blur=8,
                gradient={0.2: 'blue', 0.4: 'cyan', 0.6: 'lime', 0.8: 'yellow', 1.0: 'red'}, show=False,
            ).add_to(self.map)

    def _add_biodiversity_statistics(self, species_gdf):
        total_species = species_gdf["scientificName"].nunique()
        total_records = len(species_gdf)
        taxonomic_counts = species_gdf["taxonomic_group"].value_counts().to_dict()
        zone_counts = species_gdf["biogeographic_zone"].value_counts().to_dict()
        threatened_count = 0
        if "iucnRedListCategory" in species_gdf.columns:
            threatened_count = len(species_gdf[species_gdf["iucnRedListCategory"].isin(["CR", "EN", "VU"])])
        stats_html = f"""<div style="position: fixed; top: 70px; left: 10px; width: 260px; max-height: calc(100vh - 140px); overflow-y: auto; background-color: rgba(255,255,255,0.9); border: 1px solid #bbb; z-index:9998; font-size:12px; padding: 12px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1);">
            <h4 style="color: darkgreen; margin-top: 0; margin-bottom: 10px; font-size: 14px; border-bottom: 1px solid #eee; padding-bottom: 5px;">📊 Biodiversity Statistics</h4>
            <div style="margin-bottom: 8px;"><b>🌍 Overall Diversity:</b><br>&bull; Total Species: {total_species:,}<br>&bull; Total Records: {total_records:,}<br>&bull; Threatened (CR,EN,VU): {threatened_count:,}<br></div>
            <div style="margin-bottom: 8px;"><b>🏛️ Taxonomic Groups:</b><br>{'<br>'.join([f"&bull; {group}: {count:,}" for group, count in taxonomic_counts.items()])}</div>
            <div><b>🗺️ Biogeographic Zones:</b><br>{'<br>'.join([f"&bull; {zone}: {count:,}" for zone, count in zone_counts.items()])}</div>
            <div style="font-size: 10px; color: #777; margin-top: 12px; border-top: 1px solid #eee; padding-top: 8px;">Data from GBIF &bull; BioNusantara Project</div></div>"""
        self.map.get_root().html.add_child(folium.Element(stats_html))

    def _add_taxonomic_group_layers(self, species_gdf):
        taxonomic_groups = sorted([tg for tg in species_gdf["taxonomic_group"].unique() if pd.notna(tg)])
        color_map = {"Mammals": "#E41A1C", "Birds": "#377EB8", "Reptiles": "#4DAF4A", "Amphibians": "#FF7F00", "Other": "#984EA3", "Unknown": "#AAAAAA"}
        for group in taxonomic_groups:
            group_data = species_gdf[species_gdf["taxonomic_group"] == group]
            if group_data.empty: continue
            feature_group = folium.FeatureGroup(name=f"{group} ({len(group_data):,} records)", show=False)
            for idx, species in group_data.iterrows():
                if pd.notna(species.geometry) and species.geometry.geom_type == 'Point':
                    lat, lon = species.geometry.y, species.geometry.x
                    folium.CircleMarker(
                        location=[lat, lon], radius=5, popup=folium.Popup(self._create_species_popup(species), max_width=400),
                        tooltip=f"{species.get('scientificName', 'Unknown')} ({group})", color="black", weight=0.5,
                        fillColor=color_map.get(group, color_map["Unknown"]), fillOpacity=0.8,
                    ).add_to(feature_group)
            feature_group.add_to(self.map)

    def _add_custom_legend(self):
        legend_html = """<div style="position: fixed; bottom: 20px; left: 10px; width: 220px; background-color: rgba(255,255,255,0.9); border:1px solid #bbb; z-index:9998; font-size:11px; padding: 10px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1);">
            <h4 style="margin-top: 0; margin-bottom:8px; color: darkgreen; font-size: 13px; border-bottom: 1px solid #eee; padding-bottom: 4px;">🗺️ Map Legend</h4>
            <div style="margin-bottom: 6px;"><b>Biogeographic Lines:</b><br><svg width="20" height="10" style="vertical-align: middle;"><line x1="0" y1="5" x2="20" y2="5" style="stroke:red; stroke-width:3; stroke-dasharray:4,2;"/></svg> Wallace Line<br><svg width="20" height="10" style="vertical-align: middle;"><line x1="0" y1="5" x2="20" y2="5" style="stroke:blue; stroke-width:3; stroke-dasharray:4,2;"/></svg> Weber Line</div>
            <div><b>Taxonomic Groups (Example):</b><br><span style="color: #E41A1C; font-size: 18px; vertical-align: middle;">●</span> Mammals<br><span style="color: #377EB8; font-size: 18px; vertical-align: middle;">●</span> Birds<br><span style="color: #4DAF4A; font-size: 18px; vertical-align: middle;">●</span> Reptiles<br><span style="color: #FF7F00; font-size: 18px; vertical-align: middle;">●</span> Amphibians</div>
            <div style="margin-bottom: 6px;"><b>Lainnya:</b><br><div style="display:inline-block; vertical-align: middle; width:24px; height:24px; background-color: #28a745; border-radius:50%; border:3px solid white; position:relative; box-shadow: 0 0 5px rgba(0,0,0,0.6);"><div style="position:absolute; width:200%; height:200%; top:-50%; left:-50%; background-color:#28a745; border-radius:50%; opacity:0.7; animation:pulse-animation 2s infinite; z-index:-1;"></div></div> Taman Nasional</div>
            <div style="font-size: 9px; color: #777; margin-top: 8px; border-top: 1px solid #eee; padding-top: 6px;">Toggle layers via control panel.<br>Click points for species details & tree.</div></div>"""
        self.map.get_root().html.add_child(folium.Element(legend_html))
    
    def _add_search_functionality(self, species_gdf_for_search):
        pass

    def save_map(self, filename="bionusantara_map.html"):
        if self.map:
            output_dir = os.path.dirname(filename)
            if output_dir and not os.path.exists(output_dir): os.makedirs(output_dir)
            self.map.save(filename)
            print(f"Map saved as {filename}")
        else: print("No map to save. Create map first using create_biodiversity_map()")

    def display_map(self):
        if self.map: return self.map
        else: print("No map to display. Create map first using create_biodiversity_map()")

if __name__ == "__main__":
    analyzer = BiogeographicAnalyzer()
    analyzer.create_wallace_weber_lines() 
    script_dir = os.path.dirname(__file__) if '__file__' in locals() else os.getcwd()
    species_data_path = os.path.join(script_dir, "../data/processed/cleaned_species_data.geojson")
    output_map_path = os.path.join(script_dir, "../output/bionusantara_biodiversity_map.html")
    phylo_tree_dir = os.path.join(script_dir, "../output/phylo_trees")
    if not os.path.exists(phylo_tree_dir): os.makedirs(phylo_tree_dir)
    national_parks_json_path = os.path.join(script_dir, "../data/processed/national_parks.json")
    if not os.path.exists(national_parks_json_path):
        print(f"Membuat file dummy {national_parks_json_path} untuk pengujian...")
        dummy_parks_data = [
            {"name": "Taman Nasional Bali Barat", "latitude": -8.1333, "longitude": 114.4833, "size": 190},
            {"name": "Taman Nasional Gunung Rinjani", "latitude": -8.40806, "longitude": 116.44944, "size": 413.3},
            {"name": "Taman Nasional Komodo", "latitude": -8.5513, "longitude": 119.4832, "size": 1733}
        ]
        processed_dir = os.path.dirname(national_parks_json_path)
        if not os.path.exists(processed_dir): os.makedirs(processed_dir)
        with open(national_parks_json_path, 'w', encoding='utf-8') as f_np:
            json.dump(dummy_parks_data, f_np, indent=4, ensure_ascii=False)
    try:
        species_gdf = gpd.read_file(species_data_path)
        print(f"Loaded {len(species_gdf)} species records from {species_data_path}")
        if 'taxonomic_group' not in species_gdf.columns:
            species_gdf['taxonomic_group'] = 'Unknown'
        classified_species = analyzer.classify_biogeographic_zones(species_gdf.copy())
        biodiversity_map = analyzer.create_biodiversity_map(
            classified_species, title="BioNusantara: Indonesian Biodiversity Explorer"
        )
        analyzer.save_map(output_map_path)
        print("\033[92mBiodiversity map created successfully!\033[0m")
        print(f"\033[93mSaved in '{output_map_path}'\033[0m")
    except FileNotFoundError:
        print("\033[91mSpecies data file not found at: {species_data_path}\033[0m")
    except Exception as e:
        import traceback
        print(f"\033[91mAn error occurred while creating the map: {e}\033[0m")
        print(traceback.format_exc())
