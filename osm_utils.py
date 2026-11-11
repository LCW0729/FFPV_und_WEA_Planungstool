import osmnx as ox

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
            return ox.features_from_bbox(bbox, tags=tags)
        except ox._errors.InsufficientResponseError:
            raise  # no data in the area, handled by the callers
        except Exception as e:  # connection refused, timeout, HTTP 429/504, ...
            last_error = e
    raise last_error
