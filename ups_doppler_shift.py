from datetime import datetime, timezone, timedelta
import numpy as np
import math
import matplotlib.pyplot as plt

from skyfield.api import load, wgs84, EarthSatellite

GS_LAT_DEG = 37.2299444
GS_LON_DEG = -80.4396389
GS_ALT_M = 634 

F_TX_HZ = 401_000_000   # Hz
C_M_S = 299_792_458.0    # m/s

ts = load.timescale()  # leap seconds aware
gs = wgs84.latlon(GS_LAT_DEG, GS_LON_DEG, elevation_m=GS_ALT_M)

#Orbital Elements
BEST_PASS_INDEX = 1      # assume it exists
SMA_KM       = 6853      # semi-major axis a [km]
ECC          = 0.01 
INC_DEG      = 45.01
RAAN_DEG     = 0
ARGP_DEG     = 0
TA_DEG       = 0         # true anomaly at epoch
EPOCH_STR = "2026-03-01T00:00:00.000000"

MU_EARTH_KM3_S2 = 398600.4418  # Earth's gravitational parameter [km^3/s^2]

def true_to_mean_anomaly_deg(e, ta_deg):
    """Convert true anomaly (deg) → mean anomaly (deg)."""
    ta = math.radians(ta_deg)
    cosE = (e + math.cos(ta)) / (1.0 + e * math.cos(ta))
    sinE = (math.sqrt(1.0 - e*e) * math.sin(ta)) / (1.0 + e * math.cos(ta))
    E = math.atan2(sinE, cosE)  # radians in [-pi, pi]
    M = E - e * math.sin(E)     # radians
    return math.degrees(M) % 360.0

def sma_km_to_mean_motion_rev_per_day(a_km, mu_km3_s2=MU_EARTH_KM3_S2):
    """Convert semi-major axis to mean motion (rev/day)."""
    n_rad_s = math.sqrt(mu_km3_s2 / (a_km**3))
    return n_rad_s / (2.0 * math.pi) * 86400.0

def parse_iso_utc(s: str) -> datetime:
    """
    Parse ISO-8601 strings like 'YYYY-MM-DDTHH:MM:SS[.ffffff]Z' or with offsets.
    Returns a timezone-aware UTC datetime.
    """
    if s.endswith('Z'):
        s = s[:-1] + '+00:00'
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)

def group_triplets(events):
    out = []
    i = 0
    while i + 2 < len(events):
        if events[i] == 0 and events[i+1] == 1 and events[i+2] == 2:
            out.append((i, i+1, i+2))
            i += 3
        else:
            i += 1
    return out

M_deg = true_to_mean_anomaly_deg(ECC, TA_DEG)
n_rev_per_day = sma_km_to_mean_motion_rev_per_day(SMA_KM)

fields = {
    "OBJECT_NAME": "UtProSat-1",
    "OBJECT_ID": "2025-001A",
    "EPOCH": EPOCH_STR,  # microseconds, no 'Z'
    "MEAN_MOTION": n_rev_per_day,           # rev/day
    "MEAN_MOTION_DOT": 0.0,                 # 1/day^2 (set if known)
    "MEAN_MOTION_DDOT": 0.0,                # 1/day^3 (set if known)
    "ECCENTRICITY": ECC,
    "INCLINATION": INC_DEG,
    "RA_OF_ASC_NODE": RAAN_DEG,
    "ARG_OF_PERICENTER": ARGP_DEG,
    "MEAN_ANOMALY": M_deg,
    "EPHEMERIS_TYPE": 0,
    "CLASSIFICATION_TYPE": "U",
    "BSTAR": 0.0,
    "NORAD_CAT_ID": 0,
    "ELEMENT_SET_NO": 1,
    "REV_AT_EPOCH": 0,
}

sat = EarthSatellite.from_omm(ts, fields)

WINDOW_START_STR = "2026-03-01T00:00:00Z"
WINDOW_END_STR   = "2026-03-01T03:00:00Z"

#Pass window
t0 = ts.from_datetime(parse_iso_utc(WINDOW_START_STR))
t1 = ts.from_datetime(parse_iso_utc(WINDOW_END_STR))

#Pass
t_events, events = sat.find_events(gs, t0, t1, altitude_degrees=0.0)
triplets = group_triplets(np.array(events, dtype=int))
i_r, i_c, i_s = triplets[BEST_PASS_INDEX]
trise, tculm, tset = t_events[i_r], t_events[i_c], t_events[i_s]

start_dt = trise.utc_datetime().replace(tzinfo=timezone.utc)
end_dt   = tset.utc_datetime().replace(tzinfo=timezone.utc)
n_sec = max(1, int((end_dt - start_dt).total_seconds()))
dts  = [start_dt + timedelta(seconds=s) for s in np.arange(0, n_sec + 1, 1.0)] #1s step
t    = ts.from_datetimes(dts)

#Range and range rate
pos = (sat - gs).at(t)
_, _, rng, _, _, rng_rate = pos.frame_latlon_and_rates(gs)
range_km = rng.km
range_rate_m_s = rng_rate.m_per_s

doppler_hz = - (range_rate_m_s / C_M_S) * F_TX_HZ

#Text file for GNU Radio
ts_unix = [dt.timestamp() for dt in dts]
with open("doppler_curve.txt","w") as f:
    f.writelines(f"{t:.6f} {df:.6f}\n" for t, df in zip(ts_unix, doppler_hz))

# Print info
print(f"Pass (UTC)  Rise {trise.utc_strftime('%Y-%m-%d %H:%M:%S')}  "
      f"Culm {tculm.utc_strftime('%H:%M:%S')}  Set {tset.utc_strftime('%H:%M:%S')}")
print(f"Max elevation @ culm: {(sat - gs).at(tculm).altaz()[0].degrees:.1f}°")
print(f"Doppler min/max (Hz) @ {F_TX_HZ/1e6:.3f} MHz: {doppler_hz.min():.1f} / {doppler_hz.max():.1f}")
print(f"Minimum range rate (m/s): {range_rate_m_s.min():.1f} / Maximum range rate (m/s): {range_rate_m_s.max():.1f}")

#Plot
t_rel_min = (np.array([dt.timestamp() for dt in dts]) - dts[0].timestamp()) / 60.0
# Find index where doppler_hz is closest to zero
closest_idx = np.argmin(np.abs(doppler_hz))
# Shift t_rel_min so that closest approach (doppler_hz ~ 0) is at zero
t_rel_min_centered = t_rel_min - t_rel_min[closest_idx]

plt.figure(figsize=(9, 4.5))
plt.plot(t_rel_min_centered, doppler_hz, linewidth=2)
plt.axhline(0.0, linestyle='--', linewidth=1)
plt.title("Doppler S-curve (Selected Pass)")
plt.xlabel("Time relative to closest approach (min)")
plt.ylabel("Doppler shift (Hz)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()