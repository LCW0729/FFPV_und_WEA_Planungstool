import osmnx as ox
ox.settings.overpass_rate_limit = False  # skip the /status pre-check that hangs
ox.settings.requests_timeout = 60        # default 180 s makes each failed mirror cost 3 minutes

OVERPASS_URLS = [
    "https://overpass-api.de/api",
    "https://overpass.private.coffee/api",
    "https://overpass.kumi.systems/api",
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
