import osmnx as ox
ox.settings.overpass_rate_limit = False  # skip /status pre-check that hung last time
ox.settings.requests_timeout = 60
ox.settings.http_user_agent = "ffpv-wea-planungstool (Bachelorarbeit, contact: your-email@example.com)"
ox.settings.http_referer = "https://planungstool.streamlit.app/"

OVERPASS_URLS = [
    "https://overpass-api.de/api",
    "https://ovp.dresel.systems/api",
    "https://overpass.private.coffee/api",
    "https://maps.mail.ru/osm/tools/overpass/api",
]

def features_from_bbox_with_fallback(bbox, tags):
    last_error = None
    for url in OVERPASS_URLS:
        ox.settings.overpass_url = url
        try:
            features = ox.features_from_bbox(bbox, tags=tags)
            print(f"[Overpass] OK: {url}")
            return features
        except ox._errors.InsufficientResponseError:
            print(f"[Overpass] OK (no features in area): {url}")
            raise
        except Exception as e:
            print(f"[Overpass] FAILED: {url} -> {type(e).__name__}: {e}")
            last_error = e
    raise last_error
