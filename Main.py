import streamlit as st
from streamlit_folium import st_folium
from streamlit.components.v1 import html
import folium
import time
from folium.plugins import Draw, MousePosition
from geopy.geocoders import Nominatim
from geopy.exc import GeocoderTimedOut, GeocoderServiceError
from pyproj import CRS, Transformer

# streamlit run C:\Users\leech\Desktop\Quellcode\Main.py

# Solar panel specifications (in meters)
solar_panel_width = 1.96
solar_panel_height = 3.66 #Für 4 module
row_spacing = 3.5

# Wind turbines specifications (in meters)
min_spacing_x = 1270
min_spacing_y = 762
hub_height = 135

def get_local_crs(lon, lat):
    return CRS.from_proj4(
        f"+proj=tmerc +lat_0={lat} +lon_0={lon} +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    )

def find_circle_markers(obj):
    markers = []
    if isinstance(obj, folium.CircleMarker):
        markers.append(obj)
    if hasattr(obj, '_children'):
        for child in obj._children.values():
            markers.extend(find_circle_markers(child))
    return markers

@st.dialog("Anleitung und Hinweise", width="large") 
def show_instructions():
    st.markdown(
        """
        **Anleitung:** 
        1. Wählen Sie die Art der Anlage, die Sie planen möchten (FFPV, WEA oder Hybrid).  
        2. Geben Sie einen Ort oder eine Adresse ein, um die Karte zu zentrieren. Alternativ können Sie die Karte manuell verschieben und zoomen.
        3. Klicken Sie auf das Kreissymbol oben rechts auf der Karte und ziehen Sie den Kreis auf der Karte auf, indem Sie auf eine Stelle auf der Karte klicken und nach außen ziehen. Lassen Sie die Maus los und das Programm wird automatisch durchgeführt. Der Kreis stellt das zu beplanende Gebiet dar.
        4. Scrollen Sie nach unten, um die Ergebnisse der Simulation zu sehen.
        5. Um einen neuen Kreis zu zeichnen, klicken Sie zunächst auf das Mülleimersymbol auf der rechten Seite der Karte unter dem Kreissymbol. Klicken Sie dann auf „Alles löschen“, um den aktuellen Kreis zu löschen.
        6. Zeichnen Sie einen neuen Kreis, indem Sie erneut auf das Kreissymbol klicken und den Kreis auf der Karte ziehen.

        **Hinweise:**  
        *Dieses Programm funktioniert derzeit nur für Regionen innerhalb Deutschlands.  
        *Eine Infobox unterhalb der Karte zeigt den Fortschritt der Simulation an. Alternativ bedeutet ein „RUNNING...“,-Symbol oben rechts im Browser, dass das Programm gerade läuft.   
        *In der Infobox wird „Erfolgreich“ angezeigt, wenn der Vorgang abgeschlossen ist.  
        *Sie können im Diagramm verschieben und zoomen. Klicken Sie auf das Symbol oben rechts im Diagramm, um es als Vollbild anzuzeigen. Klicken Sie auf die drei Punkte oben rechts, um das Diagramm zu exportieren. 
        *Sie können ebenfalls die Karte verschieben und zoomen. Oben rechts auf der Karte können Sie die Layers ein- oder ausblenden.  
        *Klicken Sie auf die Schaltfläche „Karte als HTML herunterladen“, um die Karte als HTML-Datei herunterzuladen.  

        Die in diesem Tool erstellten Darstellungen und Simulationen stellen keine vollständige Planung realer 
        Solar- oder Windparks dar. Vielmehr handelt es sich um eine theoretische Abschätzung des Flächenpotenzials 
        und der potenziellen Energieerzeugung in einem definierten Gebiet. Die Positionierung der Solarmodule und 
        Windturbinen basiert ausschließlich auf öffentlich zugänglichen Geodaten (z.B. OpenStreetMap) und standardisierten 
        Ausschlusskriterien, ohne Prüfung standortbezogener Genehmigungsbedingungen, Netzanschlussoptionen oder detaillierter 
        Umweltverträglichkeitsprüfungen. Dieses Werkzeug dient daher in erster Linie der automatisierten Potenzialanalyse 
        und nicht der Erstellung genehmigungsfähiger Planungsunterlagen.
        """
    )

# Selector at the top
st.title("Programm zur Planung von Solar- und Windkraftanlagen")

if st.button("Anleitung und Hinweise"):
    show_instructions()  

option = st.radio(
    "Welche Art von Anlage möchten Sie planen?",
    ("FFPV", "WEA", "Hybrid (FFPV + WEA)")
)

# Search functionality
geolocator = Nominatim(user_agent="solar-farm-planner")
location_input = st.text_input("Suchen Sie nach einem Ort oder einer Adresse:", "")

# Step 2: Initialize session state to avoid repeated queries
if 'geocoded_results' not in st.session_state:
    st.session_state['geocoded_results'] = None
if 'last_query' not in st.session_state:
    st.session_state['last_query'] = ""

# Default center
map_center = [51.1657, 10.4515]
zoom_level = 6
location = None

# Step 4: Perform geocoding only when query is new
if location_input and location_input != st.session_state['last_query']:
    max_retries = 3
    for attempt in range(max_retries):
        try:
            locations = geolocator.geocode(location_input, exactly_one=False, addressdetails=True, limit=5)
            if locations:
                st.session_state['geocoded_results'] = locations
                st.session_state['last_query'] = location_input
            else:
                st.warning("Keine Adresse gefunden")
                st.session_state['geocoded_results'] = None
        except (GeocoderTimedOut, GeocoderServiceError) as e:
            if attempt < max_retries - 1:
                time.sleep(1)
            else:
                st.error(f"Geocoding error: {e}")
                st.session_state['geocoded_results'] = None

# Step 5: Show dropdown and update map if results exist
if st.session_state['geocoded_results']:
    locations = st.session_state['geocoded_results']
    location_options = [f"{loc.address} ({loc.latitude:.4f}, {loc.longitude:.4f})" for loc in locations]
    selection = st.selectbox("Ausgewählte Adresse:", location_options)

    selected_index = location_options.index(selection)
    selected_location = locations[selected_index]
    st.session_state["selected_location"] = selected_location
    map_center = [selected_location.latitude, selected_location.longitude]
    zoom_level = 15

m = folium.Map(location=map_center, zoom_start=zoom_level)
mouse_position = MousePosition(position='bottomright', separator=' | ', prefix="Coordinates:", num_digits=6)
m.add_child(mouse_position)

# drop marker at the searched location
if "selected_location" in st.session_state:
    folium.Marker(
        location=[st.session_state["selected_location"].latitude, st.session_state["selected_location"].longitude],
        icon=folium.Icon(color="blue", icon="search")
    ).add_to(m)

# Add a draw tool to the map (for drawing circles only)
draw = Draw(
    export=False,  # Enable exporting to GeoJSON
    draw_options={
        'polyline': False,    # Disable lines
        'polygon': False,     # Disable polygons
        'rectangle': False,   # Disable rectangles
        'circle': True,       # Enable circles
        'marker': False,      # Disable markers
        'circlemarker': False # Disable circle markers
    },
    edit_options={
        'edit': False,        # Disable editing
        'remove': True        # Enable deleting
    }
)
draw.add_to(m)

# Display the map and get draw data
with st.container():
    output = st_folium(m, width=700, height=500, key="map_draw")

status_box = st.empty()

if "last_circle_radius" not in st.session_state:
    st.session_state["last_circle_radius"] = -1

if "last_circle_coordinates" not in st.session_state:
    st.session_state["last_circle_coordinates"] = [0, 0]

if "num_panels" not in st.session_state:
    st.session_state["num_panels"] = None

if "num_turbines" not in st.session_state:
    st.session_state["num_turbines"] = None

if "second_map" not in st.session_state:
    st.session_state["second_map"] = None

if "results_df" not in st.session_state:
    st.session_state["results_df"] = None

if "total_energy" not in st.session_state:
    st.session_state["total_energy"] = None

if "results_ac" not in st.session_state:
    st.session_state["results_ac"] = None

if "rated_power_solar" not in st.session_state:
    st.session_state["rated_power_solar"] = None

if "rated_power_wind" not in st.session_state:
    st.session_state["rated_power_wind"] = None

# Check if the user has drawn a circle
if output['last_active_drawing']:

    current_circle_radius = output['last_active_drawing']['properties']['radius']
    current_circle_coordinates = output['last_active_drawing']['geometry']['coordinates']

    if current_circle_radius != st.session_state["last_circle_radius"] or \
    st.session_state["last_circle_coordinates"][0] != current_circle_coordinates[0] or \
    st.session_state["last_circle_coordinates"][1] != current_circle_coordinates[1]:

        st.session_state["last_circle_radius"] = current_circle_radius
        st.session_state["last_circle_coordinates"] = current_circle_coordinates

        center = output['last_active_drawing']['geometry']['coordinates']
        radius = output['last_active_drawing']['properties']['radius']  # Radius is in meters

        lat_center, lon_center = center[1], center[0]

        # Initalize m2
        m2 = folium.Map(location=[lat_center, lon_center], zoom_start=20)
        mouse_position = MousePosition(position='bottomright', separator=' | ', prefix="Coordinates:", num_digits=6)

        if option == "FFPV":
            from SolarFunctions import fetch_obstacles_solar, packing_solar, simulate_solarfarm_output
            
            obstacles = fetch_obstacles_solar(output, center, status_box)

            m2_solar, num_panels = packing_solar(lat_center, lon_center, radius, solar_panel_width, 
                                                 solar_panel_height, row_spacing, obstacles, 
                                                 status_box, m2)

            st.session_state["num_panels"] = num_panels
            st.session_state["second_map"] = m2_solar
            
            st.session_state["second_map_html"] = m2_solar._repr_html_()
            
            status_box.info("🧮 Solaranlagen simulieren...")
            results_ac, rated_power_solar = simulate_solarfarm_output(lat_center, lon_center, num_panels)

            st.session_state["results_ac"] = results_ac
            st.session_state["rated_power_solar"] = rated_power_solar

        elif option == "WEA":
            from WindFunctions import fetch_obstacles_wind, packing_wind, simulate_windfarm_output, get_weather_for_windpowerlib
            
            obstacles = fetch_obstacles_wind(output, center, status_box, min_spacing_x, min_spacing_y)
        
            status_box.info("☁️ Wetterdaten abrufen...")

            weather_df, main_dir = get_weather_for_windpowerlib(st.session_state["last_circle_coordinates"][1], st.session_state["last_circle_coordinates"][0], year=2024)
            m2_wind, num_turbines = packing_wind(lat_center, lon_center, radius, 
                                                 min_spacing_x, min_spacing_y, obstacles, 
                                                 main_dir, status_box, m2, option)

            st.session_state["num_turbines"] = num_turbines
            st.session_state["second_map"] = m2_wind

            st.session_state["second_map_html"] = m2_wind._repr_html_()

            status_box.info("🧮 Windturbinen simulieren...")
            results_df, total_energy, rated_power_wind = simulate_windfarm_output(weather_df, num_turbines, hub_height)

            st.session_state["results_df"] = results_df
            st.session_state["total_energy"] = total_energy
            st.session_state["rated_power_wind"] = rated_power_wind
        
        elif option == "Hybrid (FFPV + WEA)":
            from WindFunctions import fetch_obstacles_wind, packing_wind, simulate_windfarm_output, get_weather_for_windpowerlib
            from SolarFunctions import fetch_obstacles_solar, packing_solar, simulate_solarfarm_output
            import geopandas as gpd
            import pandas as pd
            from shapely.geometry import Point

            obstacles_wind = fetch_obstacles_wind(output, center, status_box, min_spacing_x, min_spacing_y)

            status_box.info("☁️ Wetterdaten abrufen...")

            weather_df, main_dir = get_weather_for_windpowerlib(st.session_state["last_circle_coordinates"][1], st.session_state["last_circle_coordinates"][0], year=2024)
            m2_wind, num_turbines = packing_wind(lat_center, lon_center, radius, 
                                                 min_spacing_x, min_spacing_y, obstacles_wind, 
                                                 main_dir, status_box, m2, option)

            st.session_state["num_turbines"] = num_turbines

            status_box.info("🧮 Windturbinen simulieren...")
            results_df, total_energy, rated_power_wind = simulate_windfarm_output(weather_df, num_turbines, hub_height)

            st.session_state["results_df"] = results_df
            st.session_state["total_energy"] = total_energy
            st.session_state["rated_power_wind"] = rated_power_wind
            
            # Project turbine points and buffer them in meters
            crs_local = get_local_crs(lon_center, lat_center)
            to_local = Transformer.from_crs("EPSG:4326", crs_local, always_xy=True)
            to_wgs = Transformer.from_crs(crs_local, "EPSG:4326", always_xy=True)       

            circle_markers = find_circle_markers(m2_wind)

            # 3. Create buffers
            turbine_polygons = []
            shadow_polygons = []

            for marker in circle_markers:
                lat, lon = marker.location
                x, y = to_local.transform(lon, lat)
                pt = Point(x, y)

                # --- Turbinenfundament (20m radius circle) ---
                fundament = pt.buffer(20)
                if fundament.is_valid and not fundament.is_empty:
                    turbine_polygons.append(fundament)
                        
            # Convert to GeoDataFrames
            obstacle_frames = []

            if turbine_polygons:
                fundament_gdf = gpd.GeoDataFrame(
                    {"buffered_geometry": turbine_polygons, "obstacle_type": "Turbinenfundament"},
                    geometry="buffered_geometry",
                    crs=crs_local
                ).to_crs("EPSG:4326")
                obstacle_frames.append(fundament_gdf)

            # Merge all turbine obstacles
            if obstacle_frames:
                tower_gdf = pd.concat(obstacle_frames, ignore_index=True)
            else:
                tower_gdf = gpd.GeoDataFrame(
                    {"buffered_geometry": [], "obstacle_type": []},
                    geometry="buffered_geometry",
                    crs="EPSG:4326"
                )

            # Place solar panels between turbines
            obstacles_solar = fetch_obstacles_solar(output, center, status_box)
            combined_obstacles = pd.concat([obstacles_solar, tower_gdf], ignore_index=True)

            m2_wind_solar, num_panels = packing_solar(lat_center, lon_center, radius, solar_panel_width, solar_panel_height, row_spacing, combined_obstacles, status_box, m2_wind)
            results_ac, rated_power_solar = simulate_solarfarm_output(lat_center, lon_center, num_panels)
            
            st.session_state["num_panels"] = num_panels
            st.session_state["second_map"] = m2_wind_solar
            st.session_state["second_map_html"] = m2_wind_solar._repr_html_()
            st.session_state["results_ac"] = results_ac
            st.session_state["rated_power_solar"] = rated_power_solar

if st.session_state['last_circle_coordinates'][0] != 0 and st.session_state['last_circle_coordinates'][1] != 0:
    st.subheader("Kreisinformationen")
    st.write(f"Kreiszentrum (Längengrad, Breitengrad): {st.session_state['last_circle_coordinates']}")

if st.session_state['last_circle_radius'] != -1:
    st.write(f"Kreisradius: {st.session_state['last_circle_radius']:.2f} meter")

if option == "FFPV":
    if st.session_state['num_panels'] or st.session_state['num_panels'] == 0:
        st.subheader("Simulationsergebnisse")
        st.write(f"Anzahl von Solarmodulen: {st.session_state['num_panels']}")

    if st.session_state["results_ac"] is not None and not st.session_state["results_ac"].empty and st.session_state["rated_power_solar"] is not None:
        total_gwh = st.session_state["results_ac"].sum() / 1_000_000_000
        total_kwh = st.session_state["results_ac"].sum() / 1_000
        total_kwp = st.session_state["rated_power_solar"] / 1_000
        kWh_kWp = total_kwh / total_kwp if total_kwp > 0 else 0

        st.write(f"Gesamtenergieerzeugung pro Jahr: {total_gwh:.2f} GWh")
        st.write(f"kWh/kWp: {kWh_kWp:.2f} kWh/kWp/a")
        st.write("Stromproduktion der Solaranlagen im Jahresverlauf (W)")
        st.line_chart(st.session_state["results_ac"])
        
elif option == "WEA":
    if st.session_state['num_turbines'] or st.session_state['num_turbines'] == 0:
        st.subheader("Simulationsergebnisse")
        st.write(f"Anzahl von Windturbinen: {st.session_state['num_turbines']}")

    if st.session_state["results_df"] is not None and not st.session_state["results_df"].empty and st.session_state["total_energy"] is not None:
        total_gwh = st.session_state["total_energy"] / 1_000
        total_kwh = st.session_state["total_energy"] * 1_000
        total_kwp = st.session_state["rated_power_wind"]
        vollaststunden = total_kwh / total_kwp if total_kwp > 0 else 0

        st.write(f"Gesamtenergieerzeugung pro Jahr: {total_gwh:.2f} GWh")
        st.write(f"Volllaststunden: {vollaststunden:.2f} h")
        st.write("Stromproduktion der Windturbinen im Jahresverlauf (MW)")
        daily_df = st.session_state["results_df"]["power_output_MW"].resample("D").mean()
        st.line_chart(daily_df)
        
elif option == "Hybrid (FFPV + WEA)":
    if st.session_state['num_panels'] or st.session_state['num_panels'] == 0:
        st.subheader("Solar Simulationsergebnisse")
        st.write(f"Anzahl von Solarmodulen: {st.session_state['num_panels']}")

    if st.session_state["results_ac"] is not None and not st.session_state["results_ac"].empty and st.session_state["rated_power_solar"] is not None:
        total_gwh = st.session_state["results_ac"].sum() / 1_000_000_000
        total_kwh = st.session_state["results_ac"].sum() / 1_000
        total_kwp = st.session_state["rated_power_solar"] / 1_000
        kWh_kWp = total_kwh / total_kwp if total_kwp > 0 else 0

        st.write(f"Gesamtenergieerzeugung pro Jahr: {total_gwh:.2f} GWh")
        st.write(f"kWh/kWp: {kWh_kWp:.2f} kWh/kWp/a")
        st.write("Stromproduktion der Solaranlagen im Jahresverlauf (W)")
        st.line_chart(st.session_state["results_ac"])

    if st.session_state['num_turbines'] or st.session_state['num_turbines'] == 0:
        st.subheader("Wind Simulationsergebnisse")
        st.write(f"Anzahl von Windturbinen: {st.session_state['num_turbines']}")

    if st.session_state["results_df"] is not None and not st.session_state["results_df"].empty and st.session_state["total_energy"] is not None:
        total_gwh = st.session_state["total_energy"] / 1_000
        total_kwh = st.session_state["total_energy"] * 1_000
        total_kwp = st.session_state["rated_power_wind"]
        vollaststunden = total_kwh / total_kwp if total_kwp > 0 else 0

        st.write(f"Gesamtenergieerzeugung pro Jahr: {total_gwh:.2f} GWh")
        st.write(f"Volllaststunden: {vollaststunden:.2f} h")
        st.write("Stromproduktion der Windturbinen im Jahresverlauf (MW)")
        daily_df = st.session_state["results_df"]["power_output_MW"].resample("D").mean()
        st.line_chart(daily_df)

if st.session_state["second_map"] and "second_map_html" in st.session_state:

    if option == "FFPV":
        st.subheader("Karte mit Solarmodulen")
        status_box.info("🗺️ Karte mit Solarmodulen erstellen...")
        with st.container():
            legend_template = """
                <style>
                    @import url('https://fonts.googleapis.com/css2?family=Source+Sans+Pro&display=swap');

                    .legend-container {
                        position: absolute;
                        top: 520px;
                        left: 0px;
                        width: 700px;
                        padding: 10px;
                    }
                        
                    .legend-box {
                        font-family: 'Source Sans Pro', sans-serif;
                        font-size: 16px;
                        color: rgb(250, 250, 250);
                    }
                    .legend-grid {
                        display: grid;
                        grid-template-columns: repeat(3, 1fr);
                        gap: 8px 15px;
                        margin-top: 10px;
                    }

                    .legend-icon {
                        width: 16px;
                        height: 16px;
                        display: inline-block;
                        border: 1px solid white;
                        margin-right: 6px;
                        vertical-align: middle;
                        border-radius: 2px;
                    }
                </style>

                <div class="legend-container">
                    <div class="legend-box">
                        <strong>Hindernis-Legende</strong>
                        <div class="legend-grid">
                            <div><i style="background: gray;" class="legend-icon"></i> Verkehr</div>
                            <div><i style="background: green;" class="legend-icon"></i> Landnutzung</div>
                            <div><i style="background: red;" class="legend-icon"></i> Gebäude</div>
                            <div><i style="background: blue;" class="legend-icon"></i> Gewässer</div>
                            <div><i style="background: purple;" class="legend-icon"></i> Schutzgebiet</div>
                            <div><i style="background: orange;" class="legend-icon"></i> Stromleitung</div>
                        </div>
                    </div>
                </div>
            """
            map_with_legend = f"""
                <style>
                    /* Remove all margin/padding from the outer container */
                    .map-wrapper {{
                        margin: 0;
                        padding: 0;
                        width: 700px;
                        height: 500px;
                    }}

                    /* Force the iframe Folium uses to render map */
                    .map-wrapper iframe {{
                        width: 700px !important;
                        height: 500px !important;
                        margin: 0 !important;
                        padding: 0 !important;
                        border: none !important;
                        display: block;
                        position: relative;
                    }}
                </style>

                <div class="map-wrapper">
                    {st.session_state["second_map_html"]}
                    {legend_template}
                </div>
            """
            html(map_with_legend, height=700)

    elif option == "WEA":
        st.subheader("Karte mit Windturbinen")
        status_box.info("🗺️ Karte mit Windturbinen erstellen...")
        with st.container():
            legend_template = """
                <style>
                    @import url('https://fonts.googleapis.com/css2?family=Source+Sans+Pro&display=swap');

                    .legend-container {
                        position: absolute;
                        top: 520px;
                        left: 0px;
                        width: 700px;
                        padding: 10px;
                    }
                        
                    .legend-box {
                        font-family: 'Source Sans Pro', sans-serif;
                        font-size: 16px;
                        color: rgb(250, 250, 250);
                    }
                    .legend-grid {
                        display: grid;
                        grid-template-columns: repeat(3, 1fr);
                        gap: 8px 15px;
                        margin-top: 10px;
                    }

                    .legend-icon {
                        width: 16px;
                        height: 16px;
                        display: inline-block;
                        border: 1px solid white;
                        margin-right: 6px;
                        vertical-align: middle;
                        border-radius: 2px;
                    }
                </style>

                <div class="legend-container">
                    <div class="legend-box">
                        <strong>Hindernis-Legende</strong>
                        <div class="legend-grid">
                            <div><i style="background: rgb(128, 128, 128);" class="legend-icon"></i> Verkehr</div>
                            <div><i style="background: rgb(34, 139, 34);" class="legend-icon"></i> Landnutzung</div>
                            <div><i style="background: rgb(220, 20, 60);" class="legend-icon"></i> Infrastruktur</div>
                            <div><i style="background: rgb(255, 140, 0);" class="legend-icon"></i> Militär</div>
                            <div><i style="background: rgb(128, 0, 128);" class="legend-icon"></i> Artenschutz</div>
                            <div><i style="background: rgb(30, 144, 255);" class="legend-icon"></i> Natur & Landschaft</div>
                            <div><i style="background: rgb(107, 142, 35);" class="legend-icon"></i> Wald</div>
                            <div><i style="background: rgb(32, 178, 170);" class="legend-icon"></i> Gewässer</div>
                            <div><i style="background: rgb(139, 69, 19);" class="legend-icon"></i> Sonstiges</div>
                        </div>
                    </div>
                </div>
            """
            map_with_legend = f"""
                <style>
                    /* Remove all margin/padding from the outer container */
                    .map-wrapper {{
                        margin: 0;
                        padding: 0;
                        width: 700px;
                        height: 500px;
                    }}

                    /* Force the iframe Folium uses to render map */
                    .map-wrapper iframe {{
                        width: 700px !important;
                        height: 500px !important;
                        margin: 0 !important;
                        padding: 0 !important;
                        border: none !important;
                        display: block;
                        position: relative;
                    }}
                </style>

                <div class="map-wrapper">
                    {st.session_state["second_map_html"]}
                    {legend_template}
                </div>
            """
            html(map_with_legend, height=700)
    
    elif option == "Hybrid (FFPV + WEA)":
        st.subheader("Karte mit Solarmodulen und Windturbinen")
        status_box.info("🗺️ Karte mit Windturbinen erstellen...")
        with st.container():
            legend_template = """
                <style>
                    @import url('https://fonts.googleapis.com/css2?family=Source+Sans+Pro&display=swap');

                    .legend-container {
                        position: relative;
                        margin-top: 100px;
                        width: 700px;
                        padding: 10px;
                        background: none;
                        color: white;
                    }
                        
                    .legend-box {
                        font-family: 'Source Sans Pro', sans-serif;
                        font-size: 16px;
                        color: rgb(250, 250, 250);
                        margin-bottom: 30px;
                    }
                    .legend-grid {
                        display: grid;
                        grid-template-columns: repeat(3, 1fr);
                        gap: 8px 15px;
                        margin-top: 10px;
                    }

                    .legend-icon {
                        width: 16px;
                        height: 16px;
                        display: inline-block;
                        border: 1px solid white;
                        margin-right: 6px;
                        vertical-align: middle;
                        border-radius: 2px;
                    }
                </style>

                <div class="legend-container">
                    <div class="legend-box">
                        <div class="legend-title">Hindernisse-Legende für Solar ☀️</div>
                        <div class="legend-grid">
                            <div><i style="background: gray;" class="legend-icon"></i> Verkehr</div>
                            <div><i style="background: green;" class="legend-icon"></i> Landnutzung</div>
                            <div><i style="background: red;" class="legend-icon"></i> Gebäude</div>
                            <div><i style="background: blue;" class="legend-icon"></i> Gewässer</div>
                            <div><i style="background: purple;" class="legend-icon"></i> Schutzgebiet</div>
                            <div><i style="background: orange;" class="legend-icon"></i> Stromleitung</div>
                        </div>
                    </div>

                    <div class="legend-box">
                        <div class="legend-title">Hindernisse-Legende für Wind 💨</div>
                        <div class="legend-grid">
                            <div><i style="background: rgb(128, 128, 128);" class="legend-icon"></i> Verkehr</div>
                            <div><i style="background: rgb(34, 139, 34);" class="legend-icon"></i> Landnutzung</div>
                            <div><i style="background: rgb(220, 20, 60);" class="legend-icon"></i> Infrastruktur</div>
                            <div><i style="background: rgb(255, 140, 0);" class="legend-icon"></i> Militär</div>
                            <div><i style="background: rgb(128, 0, 128);" class="legend-icon"></i> Artenschutz</div>
                            <div><i style="background: rgb(30, 144, 255);" class="legend-icon"></i> Natur & Landschaft</div>
                            <div><i style="background: rgb(107, 142, 35);" class="legend-icon"></i> Wald</div>
                            <div><i style="background: rgb(32, 178, 170);" class="legend-icon"></i> Gewässer</div>
                            <div><i style="background: rgb(139, 69, 19);" class="legend-icon"></i> Sonstiges</div>
                        </div>
                    </div>
                </div>
            """
            map_with_legend = f"""
                <style>
                    .map-wrapper {{
                        margin: 0;
                        padding: 0;
                        width: 700px;
                    }}

                    .map-wrapper iframe {{
                        width: 700px !important;
                        height: 500px !important;
                        margin: 0 !important;
                        padding: 0 !important;
                        border: none !important;
                        display: block;
                        position: relative;
                    }}
                </style>

                <div class="map-wrapper">
                    {st.session_state["second_map_html"]}
                </div>

                {legend_template}
            """
            html(map_with_legend, height=1000)

    st.download_button(
        "Karte als HTML herunterladen",
        data=st.session_state["second_map_html"],
        file_name="map_fragment.html",
        mime="text/html",
    )

    status_box.info("✅ Erfolgreich!")