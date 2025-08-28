import streamlit as st
import folium
from shapely.geometry import Polygon, mapping 
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
import osmnx as ox
import geopandas as gpd
from shapely.ops import unary_union
import pvlib
from pvlib.modelchain import ModelChain
from pvlib.location import Location
from pvlib.pvsystem import PVSystem
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS

def get_local_crs(lon, lat):
    return CRS.from_proj4(
        f"+proj=tmerc +lat_0={lat} +lon_0={lon} +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    )

def fetch_obstacles_solar(output, center, status_box):
    status_box.info("🔍 Hindernisse aus OSM abrufen...")
    # Get the bounding box coordinates for the ox.features_from_bbox function
    padding = 0.0004
    north = output['last_circle_polygon']['coordinates'][0][0][1] + padding
    south = output['last_circle_polygon']['coordinates'][0][16][1] - padding
    east = output['last_circle_polygon']['coordinates'][0][8][0] + padding
    west = output['last_circle_polygon']['coordinates'][0][24][0] - padding

    # Alle erforderliche Schlüssel (Tags)
    tags = {
        "highway": True,    "railway": True,    "landuse": True,
        "natural": True,    "waterway": True,   "leisure": True,
        "boundary": True,   "power": True,      "protect_class": True,
        "aeroway": True,    "building": True,   "flood_prone": True
    }

    # Fetch features from OSM using osmnx and return empty GeoDataFrame if no features found
    try:
        features = ox.features_from_bbox((west, south, east, north), tags=tags)
    except ox._errors.InsufficientResponseError:
        return gpd.GeoDataFrame({'buffered_geometry': []}, geometry='buffered_geometry', crs="EPSG:4326")

    status_box.info("⚙️ Hindernisse verarbeiten...")
    # Ensure all relevant columns exist
    for col in ["natural", "waterway", "landuse", "highway", 
                "railway", "leisure", "boundary", "power", 
                "building", "aeroway", "protect_class", "flood_prone"]:
        if col not in features.columns:
            features[col] = np.nan

    # Get local crs
    lat_center, lon_center = center[1], center[0]
    crs_local = get_local_crs(lon_center, lat_center)

    obstacle_geoms = []

    # Define traffic features
    traffic = features[
    (
        (features["highway"].notna()) |
        (features["railway"].notna()) |
        (features["aeroway"].isin(["Navigation aid", "aerodrome", "airstrip"]))
        ) &
        (features.geometry.type.isin(["LineString", "Polygon", "MultiPolygon", "Point"]))
    ].copy()
        # f"+proj=tmerc +lat_0={lat} +lon_0={lon} +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"


    # Process roads if any
    if not traffic.empty:
        traffic = traffic.to_crs(crs_local)

        # Puffer für Kategorie Verkehr
        def traffic_buffer(row):
            if row.get("aeroway") == "aerodrome":
                return 0
            elif row.get("aeroway") == "airstrip": 
                return 0
            elif row.get("railway") == "rail":
                return 15
            elif row.get("highway") == "motorway":
                return 40
            elif row.get("highway") == "primary":
                return 20
            elif pd.notna(row.get("highway")):
                return 5
            return 0

        # Apply buffer to traffic
        traffic["buffer_distance"] = traffic.apply(traffic_buffer, axis=1)
        traffic["buffered_geometry"] = traffic.apply(
            lambda row: row.geometry.buffer(row.buffer_distance), axis=1
        )

        # Append to GeoDataFrame
        gdf = gpd.GeoDataFrame(traffic[['buffered_geometry']], geometry="buffered_geometry")
        gdf.set_crs(crs_local, inplace=True)
        gdf['obstacle_type'] = "Verkehr" #traffic
        obstacle_geoms.append(gdf)

    # Define protected features
    protected_features = features[
        (
            (features["leisure"] == "nature_reserve") |
            (features["boundary"] == "protected_area") |
            (features["protect_class"].isin(["12", "2", "4", "97", "1", "7"]))
        ) &
        (features.geometry.type == "Polygon")
    ].copy()

    # Process protected features if any
    if not protected_features.empty:
        protected_features = protected_features.to_crs(crs_local)

        # Apply buffer to protected areas
        protected_features["buffered_geometry"] = protected_features.geometry.buffer(0)

        # Append to GeoDataFrame
        gdf = gpd.GeoDataFrame(protected_features[["buffered_geometry"]], geometry="buffered_geometry")
        gdf.set_crs(crs_local, inplace=True)
        gdf['obstacle_type'] = "Schutzgebiet" #protected
        obstacle_geoms.append(gdf)

    # Define water features
    water = features[
    (
        ((features["natural"] == "water") & (features.geometry.type.isin(["Polygon", "MultiPolygon"]))) |
        (features["waterway"].isin(["river", "canal", "dock"])) |
        (features["flood_prone"] == "yes") |
        (features["leisure"] == "swimming_pool")
    )
    ].copy()

    # Process water if any
    if not water.empty:
        water = water.to_crs(crs_local)

        # Apply buffer to water bodies
        water["buffered_geometry"] = water.geometry.buffer(0)

        # Append to GeoDataFrame
        gdf = gpd.GeoDataFrame(water[['buffered_geometry']], geometry="buffered_geometry")
        gdf.set_crs(crs_local, inplace=True)
        gdf['obstacle_type'] = "Gewässer"  #water
        obstacle_geoms.append(gdf)

    # Define power lines features
    power_lines = features[
        (features["power"].notna())
        & (features["power"] != "cable")  
    ].copy()

    # Process power lines if any
    if not power_lines.empty:
        power_lines = power_lines.to_crs(crs_local)

        # Define buffer distances by line type
        def power_buffer(row):
            if row.get("power") == "line":
                return 25
            elif row.get("power") == "minor_line":
                return 10
            elif pd.notna(row.get("power")):
                return 0
            return 0  # fallback

        # Apply buffer to power lines
        power_lines["buffer_distance"] = power_lines.apply(power_buffer, axis=1)
        power_lines["buffered_geometry"] = power_lines.apply(
            lambda row: row.geometry.buffer(row.buffer_distance), axis=1
        )

        # Append to GeoDataFrame
        gdf = gpd.GeoDataFrame(power_lines[["buffered_geometry"]], geometry="buffered_geometry")
        gdf.set_crs(crs_local, inplace=True)
        gdf['obstacle_type'] = "Stromleitung" #power lines
        obstacle_geoms.append(gdf)

    # Define landuse features
    landuse_all = features[(features["landuse"].notna()) & (features.geometry.type == "Polygon")].copy()
    allowed_landuse_types = ["farmland", "meadow", "grass", "greenfield"]
    landuse = landuse_all[~landuse_all["landuse"].isin(allowed_landuse_types)].copy()

    # Process landuse if any
    if not landuse.empty:
        landuse = landuse.to_crs(crs_local)

        # Apply buffer to landuse
        landuse['buffered_geometry'] = landuse.geometry.buffer(0)

        # Append to GeoDataFrame
        gdf = gpd.GeoDataFrame(landuse[['buffered_geometry']], geometry="buffered_geometry")
        gdf.set_crs(crs_local, inplace=True)
        gdf['obstacle_type'] = "Landnutzung" #landuse
        obstacle_geoms.append(gdf)

    # Define building features
    buildings = features[(features["building"].notna()) & (features.geometry.type == "Polygon")].copy()

    # Process buildings if any
    if not buildings.empty:
        buildings = buildings.to_crs(crs_local)

        # Define buffer distances by buildings type
        def buildings_buffer(row):
            if row.get("building") == "commercial":
                return 2
            elif row.get("building") == "industrial":
                return 2
            elif pd.notna(row.get("building")):
                return 2
            return 0  # fallback

        # Apply buffer to buildings
        buildings["buffer_distance"] = buildings.apply(buildings_buffer, axis=1)
        buildings["buffered_geometry"] = buildings.apply(
            lambda row: row.geometry.buffer(row.buffer_distance), axis=1
        )

        # Append to GeoDataFrame
        gdf = gpd.GeoDataFrame(buildings[['buffered_geometry']], geometry="buffered_geometry")
        gdf.set_crs(crs_local, inplace=True)
        gdf['obstacle_type'] = "Gebäude" #building
        obstacle_geoms.append(gdf)

    # Combine all obstacles (or create empty if none)
    if obstacle_geoms:
        combined = pd.concat(obstacle_geoms, ignore_index=True)
        combined = gpd.GeoDataFrame(combined, geometry='buffered_geometry', crs=crs_local).to_crs("EPSG:4326")
    else:
        # Return empty GeoDataFrame with correct structure if nothing is found
        combined = gpd.GeoDataFrame({'buffered_geometry': []}, geometry='buffered_geometry', crs="EPSG:4326")

    return combined
    

def packing_solar(lat_center, lon_center, radius_meters, solar_panel_width, solar_panel_height, row_spacing, obstacles, status_box, m2):
    # Convert lat/lon to UTM (meters)

    crs_local = get_local_crs(lon_center, lat_center)
    to_local = Transformer.from_crs("EPSG:4326", crs_local, always_xy=True)
    to_wgs = Transformer.from_crs(crs_local, "EPSG:4326", always_xy=True)

    x_center, y_center = to_local.transform(lon_center, lat_center)

    # Create a grid inside the circle
    x_min, x_max = x_center - radius_meters, x_center + radius_meters
    y_min, y_max = y_center - radius_meters, y_center + radius_meters

    # Generate grid cells
    x_values = np.arange(
        x_min + solar_panel_width / 2,
        x_max - solar_panel_width / 2 + 1e-9,
        solar_panel_width
    )
    y_values = np.arange(
        y_min + solar_panel_height / 2,
        y_max - solar_panel_height / 2 + row_spacing + 1e-9,
        solar_panel_height + row_spacing
    )

    # Initialize layers
    panel_layer = folium.FeatureGroup(name="Solaranlagen")
    solar_obstacle_layer = folium.FeatureGroup(name="Solar_Hindernisse")

    rect_polys = []
    num_panels = 0

    progress_bar = st.progress(0)
    total_rows = len(y_values)
    
    status_box.info("☀️ Solarmodule innerhalb des Kreises platzieren...")
    # Greedy placement: row by row
    for i, y in enumerate(y_values):
        for x in x_values:
            # Check if the rectangle is within the circle
            corners_within_circle = all(
                (corner_x - x_center) ** 2 + (corner_y - y_center) ** 2 <= radius_meters ** 2
                for corner_x, corner_y in [
                    (x - solar_panel_width / 2, y - solar_panel_height / 2),  # Bottom-left corner
                    (x + solar_panel_width / 2, y - solar_panel_height / 2),  # Bottom-right corner
                    (x + solar_panel_width / 2, y + solar_panel_height / 2),  # Top-right corner
                    (x - solar_panel_width / 2, y + solar_panel_height / 2),  # Top-left corner
                ]
            )
            if corners_within_circle:

                 # Define only Bottom-Left and Top-Right
                bottom_left = (x - solar_panel_width / 2, y - solar_panel_height / 2)
                top_right = (x + solar_panel_width / 2, y + solar_panel_height / 2)

                # Convert to lat/lon
                bottom_left_latlon = to_wgs.transform(*bottom_left)
                top_right_latlon = to_wgs.transform(*top_right)
                
                panel_polygon = Polygon([
                    (bottom_left_latlon[0], bottom_left_latlon[1]),
                    (top_right_latlon[0], bottom_left_latlon[1]),
                    (top_right_latlon[0], top_right_latlon[1]),
                    (bottom_left_latlon[0], top_right_latlon[1])
                ])

        #"+proj=tmerc +lat_0={lat} +lon_0={lon} +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"


                # Prüfung auf Überschneidung mit Hindernissen
                intersecting_obstacles = obstacles[
                    obstacles['buffered_geometry'].intersects(panel_polygon) |
                    obstacles['buffered_geometry'].contains(panel_polygon)
                ]
                if not intersecting_obstacles.empty:
                    continue
                else:
                    rect_polys.append(panel_polygon)
                    num_panels += 1

        # Update progress bar after each row
        progress_bar.progress((i + 1) / total_rows)
    
    feature_collection = {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "geometry": mapping(poly), "properties": {}}
            for poly in rect_polys
        ],
    }

    folium.GeoJson(
        feature_collection,
        name="Solaranlagen",
        style_function=lambda _: {
            "color": "orange",      # thin outline keeps the grid visible
            "weight": 1,
            "fillColor": "orange",
            "fillOpacity": 0.5,
        },
    ).add_to(panel_layer)

    # Define obstacle colors
    color_map = {
        "Verkehr": "gray",
        "Landnutzung": "green",
        "Gebäude": "red",
        "Gewässer": "blue",
        "Schutzgebiet": "purple",
        "Stromleitung": "orange",
    }

    merged_obstacles = []

    for obstacle_type in obstacles['obstacle_type'].unique():
        category_gdf = obstacles[obstacles['obstacle_type'] == obstacle_type]
        unioned = unary_union(category_gdf['buffered_geometry'])

        # Ensure union result is iterable
        geoms_to_draw = []

        if unioned.geom_type in ["Polygon", "MultiPolygon"]:
            geoms_to_draw = list(unioned.geoms) if unioned.geom_type == "MultiPolygon" else [unioned]

        elif unioned.geom_type == "GeometryCollection":
            geoms_to_draw = [g for g in unioned.geoms if g.geom_type in ["Polygon", "MultiPolygon"]]
        
        for geom in geoms_to_draw:
            merged_obstacles.append({
                "geometry": geom,
                "obstacle_type": obstacle_type
            })

    # Convert to GeoDataFrame
    merged_gdf = gpd.GeoDataFrame(merged_obstacles, geometry="geometry", crs=obstacles.crs)

    # Draw shaded obstacles
    for _, row in merged_gdf.iterrows():
        obstacle_type = row["obstacle_type"]
        color = color_map.get(obstacle_type, "red")
        geom = row["geometry"]

        if geom.geom_type == 'Polygon':
            if geom.is_empty or geom.exterior.is_empty:
                continue
            exterior_coords = [(lat, lon) for lon, lat in geom.exterior.coords]
            hole_coords = [[(lat, lon) for lon, lat in hole.coords] for hole in geom.interiors]
            folium.Polygon(
                locations=[exterior_coords] + hole_coords,
                color=color,
                weight=0.5,
                fill=True,
                fill_color=color,
                fill_opacity=0.5,
                tooltip=obstacle_type
            ).add_to(solar_obstacle_layer)
        
        elif geom.geom_type == 'MultiPolygon':
            for poly in geom.geoms:
                if poly.is_empty or poly.exterior.is_empty:
                    continue
                exterior_coords = [(lat, lon) for lon, lat in poly.exterior.coords]
                hole_coords = [[(lat, lon) for lon, lat in hole.coords] for hole in poly.interiors]

                folium.Polygon(
                    locations=[exterior_coords] + hole_coords,
                    color=color,
                    weight=0.5,
                    fill=True,
                    fill_color=color,
                    fill_opacity=0.5,
                    tooltip=obstacle_type
                ).add_to(solar_obstacle_layer)

    folium.Circle(
        location=[lat_center, lon_center],
        radius=radius_meters,
        color='rgb(255, 225, 25)',
        weight=2,
        fill=False,
        fill_opacity=0.2,
    ).add_to(m2)

    progress_bar.empty()
    m2.add_child(panel_layer)
    m2.add_child(solar_obstacle_layer)
    m2.add_child(folium.LayerControl())

    return m2, num_panels * 4 # Assuming 4 panels per rectangle arranged horizontally above each other


def simulate_solarfarm_output(lat_center, lon_center, num_panels):
    location = Location(lat_center, lon_center, tz='Europe/Berlin', altitude=0)
    cec_modules = pvlib.pvsystem.retrieve_sam('CECMod')
    cec_inverters = pvlib.pvsystem.retrieve_sam('CECInverter')
    
    module_key = next(
        k for k in cec_modules if "CS6U_330P" in k and "Canadian_Solar" in k
    )
    module = cec_modules[module_key]

    inverter_key = next(k for k in cec_inverters if "IQ7PLUS" in k)
    inverter = cec_inverters[inverter_key]

    rated_power_watts = module['I_mp_ref'] * module['V_mp_ref'] * num_panels * 4

    temperature_parameters = TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']

    system = PVSystem(surface_tilt=25, surface_azimuth=180,
                      module_parameters=module, inverter_parameters=inverter,
                      temperature_model_parameters=temperature_parameters,
                      modules_per_string=1, strings_per_inverter=1)
    
    modelchain = ModelChain(system, location, aoi_model="no_loss")
    
    times = pd.date_range(start='2024-01-01', end='2024-12-31', freq='H', tz=location.tz)
    clear_sky = location.get_clearsky(times)

    # pvlib's Perez transposition model (default) expects additional columns.
    # Provide minimal yet reasonable defaults to avoid KeyErrors.
    clear_sky["temp_air"] = 20.0            # °C – nominal operating temperature
    clear_sky["wind_speed"] = 1.0           # m/s – light breeze
    clear_sky["precipitable_water"] = 1.5   # cm – typical mid‑latitude value

    modelchain.run_model(clear_sky)
    ac_power = modelchain.results.ac * num_panels * 4
    
    return ac_power, rated_power_watts
