#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
# ==============================================================================
#  Copyright (c) 2026 JL7KHN/栃技研
#
#  [Terms of Use & Disclaimer]
#  - Free for personal use and modification.
#  - THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
#    JL7KHN / Tochigi Tech Lab SHALL NOT BE LIABLE FOR ANY CLAIMS OR
#    DAMAGES ARISING FROM THE USE OR RESULTS OF THIS PROGRAM.
#
#  [利用・免責条件]
#  ・個人的利用における本ソースコードの改変・利用は自由です。
#  ・本プログラムの結果および使用に伴ういかなる不利益・損害についても、
#    作成者（JL7KHN/栃技研）は一切の責任を負いません。
# ==============================================================================
"""
ESDUCT.py  Ver1.0  All-in-One (単一ファイル版)
[ESDUCT: Es & Duct Propagation Forecast Program]
  Es(スポラディックE層)予測 + DUCT(対流圏ダクト)予測の統合プログラム。
  対応環境: Windows / Raspberry Pi / Android(Pydroid3) / Chromebook / Termux

  F2層(電離層F2層)伝搬予測機能は、NICT側からの安定したデータ取得が実運用上
  不可能なため廃止した(fetch_nict_foF2/compute_f2_layer_forecastは常に空を
  返すのみとし、関連パネルは自動的に非表示となる)。

  本ファイルの詳しい変更履歴・依存パッケージ・起動オプション・理論的根拠の
  解説は、ファイル末尾の付録コメント(「付録A: 変更履歴」)にまとめて記載している。

  ------------------------------------------------------------------------
  目次 (本ファイル内を検索する際は "■ PART" でジャンプすると分かりやすい)
  ------------------------------------------------------------------------
    PART 1  : ヘッダー・基本設定 / 定数
    PART 2  : 共通ユーティリティ (プラットフォーム判定・汎用関数)
    PART 3  : Es予測エンジン - NICTデータ取得・解析・履歴・自動補正
    PART 4  : Es予測エンジン - 物理特徴量・宇宙天気・FT8/PSKReporter伝播
    PART 5  : Es予測エンジン - エリア別学習・成熟度管理・複数端末統合
    PART 6  : Es予測エンジン - 判定ロジック・コンソール表示・起動処理(main)
    PART 7  : Webダッシュボード (GUI、Es/F2タブ・DUCTタブ共通)
    PART 8  : DUCT予測エンジン (対流圏ダクト伝搬)
    PART 9  : 付録A - 変更履歴
  ------------------------------------------------------------------------
"""


import argparse
import json
import math
import os
import platform
import re
import sys
import time
import select
try:
    import msvcrt
    _MSVCRT_AVAILABLE = True
except ImportError:
    _MSVCRT_AVAILABLE = False
import warnings
import xml.etree.ElementTree as ET
import io
import html as _html_mod
import threading
import http.server
import socketserver
import functools
import socket
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from PIL import Image, ImageFilter

try:
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    try:
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score as _sk_silhouette
        from sklearn.metrics import adjusted_rand_score as _sk_ari
        _SKLEARN_CLUSTER_OK = True
    except ImportError:
        _SKLEARN_CLUSTER_OK = False
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False



JST = timezone(timedelta(hours=9))

SOURCE_VERSION = "1.0"
SYSTEM_NAME    = "ESDUCT"
VERSION = SOURCE_VERSION
VERSION_LABEL = f"Ver{SOURCE_VERSION}"

def display_version() -> str:
    """システム名＋バージョンの固定表示文字列を返す。
    成熟度表示はprint_edfs_status()のLevel(0-5)で行うため、
    バージョン文字列はモデル状態に連動しない。
    """
    return f"{SYSTEM_NAME} {VERSION_LABEL}"

TITLE = f"{SYSTEM_NAME} {VERSION_LABEL}  (Es & Duct Propagation Forecast Program)"

NICT_FXES_URL_LEGACY = "https://wdc.nict.go.jp/Ionosphere/realtime/fxEs/latest-fxEs.html"
NICT_FXES_URL_DEFAULT = NICT_FXES_URL_LEGACY
NICT_IONOS_URL_CANDIDATES = [
    "https://wdc.nict.go.jp/Ionosphere/realtime/ISDJ/ionospheric-alert.html",
    "https://wdc.nict.go.jp/Ionosphere/realtime/ISDJ/ionospheric-signal.html",
    "https://wdc.nict.go.jp/Ionosphere/realtime/ionos/latest-ionos.html",
]
NICT_IONOS_URL = NICT_IONOS_URL_CANDIDATES[0]

PSKREPORTER_QUERY_URL = "https://retrieve.pskreporter.info/query"
DEFAULT_PSK_CALLSIGN = "JL7KHN/P"
PSK_MODE = "FT8"
PSK_FREQ_RANGE_HZ = (28000000, 28100000)
FT8_STALE_MIN_THRESHOLD = 60.0

NOAA_SWPC_MAG_URL      = "https://services.swpc.noaa.gov/products/solar-wind/mag-1-day.json"
NOAA_SWPC_PLASMA_URL   = "https://services.swpc.noaa.gov/products/solar-wind/plasma-1-day.json"
NOAA_SWPC_PROTON_URL   = "https://services.swpc.noaa.gov/json/goes/primary/xrays-1-day.json"
NOAA_SWPC_KP_URL       = "https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json"
NOAA_SWPC_ELECTRON_URL = "https://services.swpc.noaa.gov/json/goes/primary/electrons-1-day.json"

NOAA_SWPC_F107_SUMMARY_URL = "https://services.swpc.noaa.gov/products/summary/10cm-flux.json"
NOAA_SWPC_F107_URL         = "https://services.swpc.noaa.gov/json/f107_cm_flux.json"
SW_FETCH_TIMEOUT       = 10
SW_FETCH_TIMEOUT_BARGRAPH = 5
SW_HISTORY_MAXLEN      = 48
SW_CACHE_TTL_HOURS     = 24

NICT_IONOGRAM_URLS = {
    "okinawa":   "https://wdc.nict.go.jp/Ionosphere/data/ionogram/latest_oki.gif",
    "yamagawa":  "https://wdc.nict.go.jp/Ionosphere/data/ionogram/latest_yam.gif",
    "kokubunji": "https://wdc.nict.go.jp/Ionosphere/data/ionogram/latest_kok.gif",
    "wakkanai":  "https://wdc.nict.go.jp/Ionosphere/data/ionogram/latest_wak.gif",
}

EQI_W_FXES       = 0.30
EQI_W_FOES       = 0.15
EQI_W_LAYER_CNT  = 0.10
EQI_W_REFLECT    = 0.15
EQI_W_SPREAD     = 0.10
EQI_W_TRACE      = 0.20

ES_PHASE_FORMATION_FXES = 5.0
ES_PHASE_GROWTH_FXES    = 8.0
ES_PHASE_PEAK_FXES      = 12.0
ES_PHASE_DECAY_DEQI     = -0.05
ES_PHASE_COLLAPSE_FXES  = 5.0
ES_PHASE_MIN_SAMPLES    = 5

PERM_IMP_MIN_RATIO      = 0.05

LAG_STEPS = [0, 1, 2, 3, 4]

INTERACTION_PAIRS = [
    ("eqi_norm",     "fxEs_norm"),
    ("eqi_norm",     "mstid_approach_index"),
    ("sw_score",     "fxEs_norm"),
    ("sw_score",     "mstid_approach_index"),
    ("eqi_norm",     "sw_score"),
]

FORECAST_CYCLES       = 2
FXES_THRESHOLD        = 8.0
ES_SAMPLE_THR_RPI     = 0.5
ES_SAMPLE_THR_MODE    = 0.5
NICT_FXES_PLAUSIBLE_MIN = 0.0
NICT_FXES_PLAUSIBLE_MAX = 30.0
NICT_GRADIENT_MAX_AGE_MIN = 60.0

AREA_PRED_LEARN_START_TOTAL  = 100
AREA_PRED_LEARN_START_ES     = 30
AREA_PRED_LOW_TOTAL          = 100; AREA_PRED_LOW_ES   = 30;  AREA_PRED_LOW_R2   = 0.20
AREA_PRED_MED_TOTAL          = 200; AREA_PRED_MED_ES   = 60;  AREA_PRED_MED_R2   = 0.35
AREA_PRED_HIGH_TOTAL         = 500; AREA_PRED_HIGH_ES  = 150; AREA_PRED_HIGH_R2  = 0.50
AREA_PRED_MIN_HOLDOUT_SKILL  = 0.05

SKILL_RECENT_WEIGHT_ALPHA = 0.6
SKILL_REFERENCE = 0.3
GEOMETRY_CONSISTENCY_RPI_THRESHOLD = 0.4
GEOMETRY_CONFLICT_WEIGHT_DAMPING = 0.5
FORECAST_HORIZON_MIN_APPROX = FORECAST_CYCLES * (900 / 60.0)
ENSEMBLE_SHADOW_EVAL_MAXLEN = 200

MAX_BIAS_CORRECTION   = 15.0
BIAS_STD_HEALTH_THR   = 20.0
BIAS_HISTORY_MAXLEN   = 100

VER113_SAMPLE_MIN     = 300
VER113_ES_MIN         = 100
VER113_SILHOUETTE_MIN = 0.45
VER113_R2GAIN_MIN     = 0.10
VER113_ARI_MIN        = 0.70
VER113_SEASON_DAYS    = 30

MODEL_SCORE_W_SAMPLE  = 0.30
MODEL_SCORE_W_ES      = 0.20
MODEL_SCORE_W_SIL     = 0.20
MODEL_SCORE_W_ARI     = 0.20
MODEL_SCORE_W_R2GAIN  = 0.10

BUILTIN_FIXED_MONITOR = [
    ("JA1OVD",35.45,139.63,"1"),("JH1IFS",35.68,139.69,"1"),
    ("JG1GCO",35.60,140.12,"1"),("JA1RUJ",35.85,139.64,"1"),
    ("JQ1BWT",36.34,140.44,"1"),("JE1SGH",35.53,139.46,"1"),
    ("JF1EEA",35.70,139.94,"1"),("JH1AZO",35.51,134.12,"1"),
    ("JI1IWI",35.31,139.55,"1"),("JK1GPD",36.56,139.88,"1"),
    ("JA2GRC",34.97,138.38,"2"),("JE2SJH",35.18,136.90,"2"),
    ("JF2IWL",35.41,136.75,"2"),("JG2VTD",34.73,136.51,"2"),
    ("JK2XXK",35.18,136.90,"2"),("JA2ACX",34.96,137.16,"2"),
    ("JE2FJI",34.70,137.73,"2"),("JH2AMH",34.97,138.38,"2"),
    ("JI2XWB",35.39,136.72,"2"),("JR2AWS",35.10,136.81,"2"),
    ("JA3ENE",34.69,135.50,"3"),("JH3ECA",34.69,135.19,"3"),
    ("JI3GGO",35.01,135.75,"3"),("JO3RUL",35.00,135.86,"3"),
    ("JA3UWB",34.68,135.83,"3"),("JE3KBP",34.76,135.41,"3"),
    ("JF3KNW",34.54,135.51,"3"),("JH3HGI",34.81,135.54,"3"),
    ("JI3BNB",34.40,135.30,"3"),("JN3ANO",34.73,135.33,"3"),
    ("JA4GXS",34.66,133.93,"4"),("JH4MGU",34.39,132.46,"4"),
    ("JR4OZR",34.18,131.47,"4"),("JN4MMO",35.47,133.05,"4"),
    ("JA4LCI",35.50,134.23,"4"),("JA4BUA",34.39,132.46,"4"),
    ("JH4UTP",34.66,133.92,"4"),("JI4RGO",34.18,131.47,"4"),
    ("JO4EFC",34.61,133.64,"4"),("JR4DHK",34.80,132.71,"4"),
    ("JA5AQC",34.34,134.04,"5"),("JH5GEN",34.06,134.55,"5"),
    ("JR5JAQ",33.84,132.76,"5"),("JA5RE",33.56,133.53,"5"),
    ("JR5PPN",34.34,134.04,"5"),("JA5FBZ",33.84,132.75,"5"),
    ("JH5HDA",34.06,134.56,"5"),("JI5XSL",33.56,133.53,"5"),
    ("JJ5AKK",34.13,133.11,"5"),("JR5EFI",34.14,134.54,"5"),
    ("JA6BZI",33.60,130.41,"6"),("JH6BPG",32.79,130.74,"6"),
    ("JE6IYN",32.74,129.87,"6"),("JM6EBU",31.56,130.55,"6"),
    ("JA6SRB",33.23,131.61,"6"),("JF6RIM",33.59,130.40,"6"),
    ("JH6ARA",32.51,130.60,"6"),("JI6DUE",31.56,130.55,"6"),
    ("JH6QIL",33.23,131.61,"6"),("JA6FFK",33.88,130.88,"6"),
    ("JA7BJS",38.26,140.87,"7"),("JH7RTQ",37.75,140.46,"7"),
    ("JA7QOU",38.25,140.33,"7"),("JR7IWL",39.70,141.15,"7"),
    ("JE7JZK",39.71,140.10,"7"),("JA7FKF",40.82,140.74,"7"),
    ("JH7VVR",38.26,140.88,"7"),("JI7FBM",37.49,140.38,"7"),
    ("JK7LXU",39.70,141.15,"7"),("JR7MAZ",37.75,140.47,"7"),
    ("JA8EAT",43.06,141.35,"8"),("JH8DBJ",43.77,142.35,"8"),
    ("JR8OFE",41.77,140.74,"8"),("JE8KKX",42.98,144.38,"8"),
    ("JA8DJY",42.92,143.20,"8"),("JA8CDG",43.06,141.35,"8"),
    ("JH8FIH",43.77,142.36,"8"),("JI8FXA",44.02,144.27,"8"),
    ("JK8PBO",43.80,143.89,"8"),("JR8SWR",42.63,141.60,"8"),
    ("JA9FFZ",36.06,136.22,"9"),("JH9FCP",36.59,136.62,"9"),
    ("JA9MAT",36.69,137.21,"9"),("JE9WWB",36.59,136.62,"9"),
    ("JF9JTS",36.06,136.22,"9"),("JA9AMR",36.70,137.22,"9"),
    ("JH9DRA",36.59,136.63,"9"),("JI9IJP",36.06,136.21,"9"),
    ("JJ9ENO",35.74,136.19,"9"),("JR9DFO",36.39,136.45,"9"),
    ("JA0BDJ",36.65,138.18,"0"),("JH0INX",37.90,139.02,"0"),
    ("JR0DVM",36.65,138.18,"0"),("JA0FVU",37.90,139.02,"0"),
    ("JH0BQF",36.65,138.18,"0"),("JA0AVS",37.90,139.03,"0"),
    ("JH0CEO",36.31,137.97,"0"),("JI0VWL",36.65,138.19,"0"),
    ("JJ0FSM",37.40,138.92,"0"),("JR0HYT",36.14,137.96,"0"),
    ("JR6CSY",26.21,127.68,"JR6"),("JS6TMW",26.20,127.65,"JR6"),
    ("JS6RRR",26.33,127.80,"JR6"),("JS6QNM",26.21,127.68,"JR6"),
    ("JS6PSV",26.40,127.85,"JR6"),
    ("BV2CH",25.03,121.56,"BV"),("BU2FP",25.05,121.50,"BV"),
    ("BV5OC",24.14,120.67,"BV"),("BX4AD",23.00,120.20,"BV"),
    ("BV3CE",24.95,121.22,"BV"),
]


try:
    _SCRIPT_DIR = Path(__file__).resolve().parent
except NameError:
    _SCRIPT_DIR = Path.cwd()


def _init_cache_dir_with_fallback() -> Tuple[Path, str]:
    """学習データキャッシュフォルダを多段フォールバック付きで確保する。
    戻り値: (実際に使用するPath, 採用した経路を示す短い説明文字列)
    """
    candidates: List[Tuple[Path, str]] = [
        (_SCRIPT_DIR / "cache", "スクリプト直下"),
    ]
    is_probably_android = bool(
        os.environ.get("ANDROID_ROOT") or os.environ.get("ANDROID_ARGUMENT")
        or "pydroid" in sys.executable.lower()
    )
    if is_probably_android:
        candidates.append(
            (Path("/storage/emulated/0/ESDUCT/cache"), "Android共有ストレージ(Files/USB経由で参照可)")
        )
    try:
        candidates.append((Path.home() / "ESDUCT" / "cache", "ホームディレクトリ配下"))
    except Exception:
        pass
    import tempfile as _tempfile_mod
    candidates.append(
        (Path(_tempfile_mod.gettempdir()) / "ESDUCT_cache",
         "OS標準の一時フォルダ(最終手段・再起動で消える場合あり)")
    )

    for path, label in candidates:
        try:
            path.mkdir(parents=True, exist_ok=True)
            probe = path / ".write_test"
            with open(probe, "w") as f:
                f.write("ok")
            probe.unlink(missing_ok=True)
            return path, label
        except Exception:
            continue
    fallback = Path.cwd() / "cache"
    try:
        fallback.mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    return fallback, "カレントディレクトリ(全候補失敗のため最終フォールバック)"


DEFAULT_DATA_DIR, _CACHE_DIR_SOURCE_LABEL = _init_cache_dir_with_fallback()
print(f"[起動] 学習データキャッシュフォルダ: {DEFAULT_DATA_DIR.resolve()} ({_CACHE_DIR_SOURCE_LABEL})")
DEFAULT_WIND_PROFILE_CSV = DEFAULT_DATA_DIR / "wind_profile.csv"
DEFAULT_MSTID_DIR = DEFAULT_DATA_DIR / "mstid"
try:
    SW_CACHE_DIR = DEFAULT_DATA_DIR / "ionogram"
    SW_CACHE_DIR.mkdir(parents=True, exist_ok=True)
except Exception as _e_sw_cache:
    print(f"[起動] 警告: イオノグラム画像キャッシュフォルダの作成に失敗しました "
          f"({_e_sw_cache})。画像キャッシュ機能は無効化されますが、"
          f"予測本体の動作には影響しません。")
    SW_CACHE_DIR = DEFAULT_DATA_DIR / "ionogram"
MSTID_NICT_BASE_URL = "https://www2.nict.go.jp/spe/gps/REALTIME_GEONET"
MSTID_NICT_KIND_TABLE = {
    'A':   ('AMAP',   'Amap'),
    'H60': ('MAP',    'Hmap'),
    'H15': ('MAP15',  'Hmap15'),
    'R':   ('RMAP',   'Rmap'),
    'L':   ('LMAP',   'Lmap'),
}
MSTID_NICT_KIND = 'H15'
MSTID_CROP_FRAC = (0.08, 0.10, 0.82, 0.91)
MSTID_NICT_AUTO_ENABLED = True
MSTID_NICT_PUBLISH_LAG_MIN = 15
MSTID_AUTO_FETCH_URL = ""
MSTID_KEEP_LATEST_N = 6
DEFAULT_PROFILE_JSON = DEFAULT_DATA_DIR / "calibration_profile.json"
DEFAULT_DEBUG_HTML = DEFAULT_DATA_DIR / "debug_last_nict.html"
DEFAULT_DEBUG_TEXT = DEFAULT_DATA_DIR / "debug_last_nict_text.txt"
DEFAULT_DEBUG_TOKENS = DEFAULT_DATA_DIR / "debug_last_nict_tokens.txt"
DEFAULT_FIXED_MONITOR_CSV = DEFAULT_DATA_DIR / "fixed_monitor.csv"

TARGET_18_MHZ = 18.0
TARGET_27_MHZ = 27.0
HTTP_TIMEOUT_SEC = 20
DEFAULT_INTERVAL_SEC = 900
DASHBOARD_BIND_HOST = "127.0.0.1"
MAX_HISTORY_ROWS = 1000
MIN_ANALYSIS_POINTS = 8
MIN_AUTO_CAL_ROWS = 20
MIN_SEGMENT_DURATION_MIN = 30
MAX_ABS_TREND_MHZ_PER_HOUR = 12.0
MIN_SPEED_KMH = 20.0
MAX_SPEED_KMH = 1500.0
PROB_MAX_PERCENT = 95.0

ES_HOP_DIST_MIN_KM = 600.0
ES_HOP_DIST_MAX_KM = 2400.0
ES_MODE_SCORE_GATE = 0.15
AREA_REACH_TREND_EPS = 0.02
MIN_FIXED_STATIONS_FOR_CONFIDENCE = 5

AREA_IDW_POWER = 2.0
AREA_IDW_MIN_SOURCES = 2
AREA_IDW_MAX_NEIGHBOR_KM = 1200.0

GEO_ETA_MIN_PROJECTED_SPEED_KMH = 5.0
GEO_ETA_MAX_MINUTES = 240.0

HOUR_BIN_SIZE = 3

SELF_CLIMATOLOGY_MIN_SAMPLES = 8
SELF_CLIMATOLOGY_REFRESH_MIN = 360

NICT_CLIMATOLOGY_STATION_CODES = {
    'wakkanai':  'WK546',
    'kokubunji': 'TO536',
    'yamagawa':  'YG431',
    'okinawa':   'OK426',
}
NICT_STATION_NAME_JP = {
    'wakkanai': '稚内', 'kokubunji': '国分寺', 'yamagawa': '山川', 'okinawa': '沖縄',
}
NICT_CLIMATOLOGY_BASE_URL = "https://wdc.nict.go.jp/Ionosphere/archive/observation-history/factor-manual-{code}-{year}H.sjis.txt"
NICT_CLIMATOLOGY_YEARS_BACK = 3
NICT_CLIMATOLOGY_TIMEOUT_SEC = 60
NICT_CLIMATOLOGY_CACHE_FILE = "nict_foes_climatology.csv"
NICT_CLIMATOLOGY_RAW_CACHE_FILE = "nict_foes_raw_cache.csv"
NICT_CLIMATOLOGY_REFRESH_DAYS = 30
NICT_CLIMATOLOGY_FOES_THRESHOLD = FXES_THRESHOLD * 0.8
NICT_FOES_RAW_SCALE_DIVISOR = 10.0
NICT_CLIMATOLOGY_MIN_SOURCES = 20
NICT_CLIMATOLOGY_AUTO_ENABLED_DEFAULT = True
NICT_CLIMATOLOGY_AUTO_ENABLED = NICT_CLIMATOLOGY_AUTO_ENABLED_DEFAULT

BASE_PROFILE = {
    "prob_intercept": -2.00,
    "w_fxes": 0.32,
    "w_trend": 0.12,
    "w_cluster": 0.75,
    "w_ns": 0.25,
    "w_shear": 0.45,
    "w_layer": 0.25,
    "w_descent": 0.20,
    "w_gw_damp": 0.20,
    "w_ef_penalty": 0.15,
    "w_rising": 0.20,
    "w_mstid_approach": 0.00,
    "w_fof2": 0.00,
    "w_dfof2": 0.00,
    "w_hmf2": 0.00,
    "w_dhmf2": 0.00,
    "w_m3000f2": 0.00,
    "sign_fof2": 0.0,
    "sign_dfof2": 0.0,
    "sign_hmf2": 0.0,
    "sign_dhmf2": 0.0,
    "sign_m3000f2": 0.0,
    "lag_stability_score": None,
    "ridge_lag_stable": None,
    "f_scaler_mean": {},
    "f_scaler_scale": {},
    "w_rpi": 0.00,
    "sign_rpi": 0.0,
    "w_es_mode": 0.00,
    "sign_es_mode": 0.0,
    "w_eqi":         0.00,
    "sign_eqi":      0.0,
    "w_sw_score":    0.00,
    "sign_sw_score": 0.0,
    "w_es_phase":    0.00,
    "sign_es_phase": 0.0,
    "w_ia_eqi_fxes":    0.00,
    "w_ia_eqi_mstid":   0.00,
    "w_ia_sw_fxes":     0.00,
    "w_ia_sw_mstid":    0.00,
    "w_ia_eqi_sw":      0.00,
    "eqi_confidence":  1.0,
    "sw_quality":      1.0,
    "trace_quality":   1.0,
    "strong_open_thr": 75.0,
    "open_thr": 55.0,
    "watch_thr": 35.0,
    "weak_thr": 20.0,
    "source": "default",
}

STATIONS = {
    "okinawa": {"jp": "沖縄", "lat_order": 0, "lat": 26.68},
    "yamagawa": {"jp": "山川", "lat_order": 1, "lat": 31.20},
    "kokubunji": {"jp": "国分寺", "lat_order": 2, "lat": 35.71},
    "wakkanai": {"jp": "稚内", "lat_order": 3, "lat": 45.16},
}
DISPLAY_ORDER = ["wakkanai", "kokubunji", "yamagawa", "okinawa"]

NICT_STATION_LATLON = {
    "okinawa":   (26.68, 127.68),
    "yamagawa":  (31.20, 130.62),
    "kokubunji": (35.71, 139.49),
    "wakkanai":  (45.16, 141.80),
}

HF_VOACAP_BANDS = [
    (1.9,  "1.9MHz(160m)"),
    (3.5,  "3.5MHz(80m)"),
    (7.0,  "7MHz(40m)"),
    (10.1, "10MHz(30m)"),
    (14.0, "14MHz(20m)"),
    (18.1, "18MHz(17m)"),
    (21.0, "21MHz(15m)"),
    (24.9, "24MHz(12m)"),
    (28.5, "28MHz(10m)"),
    (50.0, "50MHz(6m)"),
]
ES_BAND_FREQS = {28.5, 50.0}

AREA_NAME_JP = {
    "0": "信越", "1": "関東", "2": "東海", "3": "関西", "4": "中国",
    "5": "四国", "6": "九州", "7": "東北", "8": "北海道", "9": "北陸",
    "JR6": "沖縄", "BV": "台湾",
}
AREA_DISPLAY_ORDER = ["8", "7", "9", "0", "1", "2", "3", "4", "5", "6", "JR6", "BV"]


def format_area_label(area_code: str) -> str:
    """Rev16.2: エリアコードをJCC/JARL方式の「Nエリア(地域名)」表示ラベルに変換する。

    【修正内容】HTML/CLI双方の表示で使うエリア名称のソースが不整合だった
    問題を修正した。具体的には、ダッシュボード側専用の _v13_AREA_JP が
    AREA_NAME_JP(こちらが正しい定義)とは異なる独自のマッピングを持っており、
    "0":"関東","1":"東京" のようにエリア番号が1つずれた誤った対応表に
    なっていた(正しくは "0":"信越","1":"関東"…であり、そもそも"東京"は
    JCC/JARL方式のエリア番号として単独では存在せず、東京は1エリア(関東)の
    一部)。本関数はAREA_NAME_JPを単一の正式なソースとして参照し、
    数字エリア(0-9)は「Nエリア(地域名)」、沖縄(JR6)/台湾(BV)のような
    非数字コードは「地域名(コード)」の形式で統一する。
    """
    code = str(area_code)
    jp = AREA_NAME_JP.get(code, code)
    if code.isdigit():
        return f"{code}エリア({jp})"
    return f"{jp}({code})"

C_RESET = "\033[0m"
C_RED = "\033[91m"
C_YELLOW = "\033[93m"
C_GREEN = "\033[92m"
C_GRAY = "\033[90m"

DISP_W_NORMAL  = 78
DISP_W_ANDROID = 48


UI_MODE_NORMAL = "NORMAL"
UI_MODE_DETAIL = "DETAIL"
UI_MODE_FACTOR = "FACTOR"
UI_MODE_EDFS   = "EDFS"
UI_MODE_DEBUG  = "DEBUG"
_edfs_stop_event = threading.Event()


def request_edfs_stop() -> None:
    """常駐中のEDFSサイクルループを安全に停止させる。

    Android/APK版main.pyの「常駐停止」ボタン、およびHTMLダッシュボードの
    「常駐停止」ボタン(/api/stop経由)の両方から呼び出される想定。
    CLI版(python EDFS_Ver15_9.py)でもimport後に呼べば同様に機能する。
    実行中のHTTPダッシュボードサーバーが動いていれば、それも停止する。
    """
    _edfs_stop_event.set()
    try:
        bridge = globals().get('_active_html_bridge')
        if bridge is not None:
            bridge.stop()
    except Exception:
        pass


def is_edfs_stop_requested() -> bool:
    return _edfs_stop_event.is_set()


_ui_display_mode = UI_MODE_NORMAL
_ui_last_results = None
_ui_last_meta = None
_ui_last_android = False
_ui_last_area_df = pd.DataFrame()
_ui_last_area_fvecs = {}
_ui_last_area_geo_eta = {}
_ui_last_area_climatology = {}
_ui_last_f2_forecast: List[Dict] = []
_ui_manual_location = None
_ui_gps_enabled = True
_ui_last_location = None
_ui_location_dialog_done = False
_ui_gps_status = "未試行"
_ui_gps_reason = ""
_ui_location_source = "未設定"
_ui_termux_gps_file = "/storage/emulated/0/EDFS/gps_location.json"
_ui_termux_gps_max_age_sec = 600



def is_windows() -> bool:
    return platform.system().lower().startswith("win")


def is_android() -> bool:
    return bool(os.environ.get("ANDROID_ARGUMENT") or os.environ.get("ANDROID_ROOT") or "pydroid" in sys.executable.lower())


def is_chromebook_arc() -> bool:
    """Rev15.9: Chromebook上のARC(Androidアプリ実行環境)で動作しているかを判定する。
    判定できない場合は安全側としてFalseを返す(=これまでの挙動を変えない)。
    ChromeOSのARCコンテナはBuild.DEVICE/Build.BRANDが"cheets"系になることが多く、
    機種によってはBuild.MODELに"Chromebook"が含まれることもあるため、両方を見る。"""
    if not is_android():
        return False
    try:
        from jnius import autoclass
        Build = autoclass('android.os.Build')
        device = str(Build.DEVICE or '').lower()
        brand = str(Build.BRAND or '').lower()
        model = str(Build.MODEL or '').lower()
        if 'cheets' in device or 'cheets' in brand:
            return True
        if 'chromebook' in model:
            return True
    except Exception:
        pass
    return False


def supports_ansi_color() -> bool:
    """Ver10.8J: ANSIカラー対応の簡易ヒューリスティック判定。
    Pydroid3/Androidのターミナルは機種・アプリ更新によって対応状況が
    バラつくため、実行時に「描画が崩れたかどうか」を確実に検知する手段は
    ない（エラーは発生せず、見た目が崩れるだけ）。そのため本関数は
    「対応していなさそうな場合のみFalse」というブラックリスト方式の
    安全側ヒューリスティックとし、最終判断は --no-color による手動復帰に委ねる。
    """
    term = os.environ.get("TERM", "").lower()
    if term == "dumb":
        return False
    return True


def enable_windows_vt100():
    if is_windows():
        try:
            os.system("")
        except Exception:
            pass
        try:
            os.system("chcp 65001 > nul 2>&1")
        except Exception:
            pass
    for stream_name in ('stdout', 'stderr'):
        try:
            getattr(sys, stream_name).reconfigure(encoding='utf-8', errors='replace')
        except AttributeError:
            try:
                stream = getattr(sys, stream_name)
                setattr(sys, stream_name, io.TextIOWrapper(
                    stream.buffer, encoding='utf-8', errors='replace',
                    line_buffering=True))
            except Exception:
                pass
        except Exception:
            try:
                getattr(sys, stream_name).reconfigure(errors='replace')
            except Exception:
                pass

enable_windows_vt100()


def now_jst() -> datetime:
    return datetime.now(JST)


def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + np.exp(-float(x)))


def clamp(x, lo, hi):
    if pd.isna(x):
        return np.nan
    return max(lo, min(hi, float(x)))


F_STD_Z_CLIP = 3.0


def std_z(raw_val, feat_key: str, profile: Dict, clip: float = F_STD_Z_CLIP) -> float:
    """Ver10.8J: _ridge_calibrate_f_features()が保存したStandardScalerの
    平均(f_scaler_mean)・標準偏差(f_scaler_scale)を使って生値を標準化する。
    学習時と推論時で同じ正規化を使うことで、Ridge係数(標準化スケール上の重み)
    と整合させる。スケーラー未学習（起動直後で自動補正が一度も成功していない等）
    の場合はNaNを返し、呼び出し側はその特徴量の寄与を0として扱う。
    """
    if pd.isna(raw_val):
        return np.nan
    mean_ = (profile.get("f_scaler_mean") or {}).get(feat_key)
    scale_ = (profile.get("f_scaler_scale") or {}).get(feat_key)
    if mean_ is None or scale_ is None or scale_ == 0:
        return np.nan
    return clamp((float(raw_val) - float(mean_)) / float(scale_), -clip, clip)


def fmt(x, nd=1, na="--") -> str:
    if x is None or pd.isna(x):
        return na
    return f"{float(x):.{nd}f}"


def fmt_mhz(x, nd=1, na="--MHz") -> str:
    if x is None or pd.isna(x):
        return na
    return f"{float(x):.{nd}f}MHz"


def safe_float(v):
    if v is None:
        return np.nan
    s = str(v).strip()
    if s in ("", "---", "----", "--", "Undefined", "undefined", "nan", "NaN", "None", "*"):
        return np.nan
    m = re.search(r"[-+]?\d+(?:\.\d+)?", s)
    return float(m.group(0)) if m else np.nan


def _color_code_for_state(text: str) -> str:
    """状態テキスト（フルラベル）から対応するANSIカラーコードを返す。空文字=無色。"""
    if "強い開通" in text or text == "強開":
        return C_RED
    if "開通" in text:
        return C_GREEN
    if "注意" in text or "監視" in text:
        return C_YELLOW
    if "弱い" in text or text == "弱兆":
        return C_GRAY
    return ""


def color_by_state(text: str, use_color: bool = True) -> str:
    if not use_color:
        return text
    code = _color_code_for_state(text)
    return code + text + C_RESET if code else text


def short_state(text: str) -> str:
    return {"強い開通": "強開", "注意監視": "注意", "弱い兆候": "弱兆"}.get(text, text)


def state_symbol(text: str) -> str:
    """Ver10.7J: カラー不使用時の状態記号ラベル（Pydroid3向け）"""
    return {
        "強い開通": "[★強開★]",
        "開通":     "[●開通 ]",
        "注意監視": "[▲注意 ]",
        "弱い兆候": "[△弱兆 ]",
        "低調":     "[×低調 ]",
        "暫定":     "[?暫定 ]",
        "判定不可": "[--不明]",
    }.get(text, f"[{text[:4]}]")


def color_or_symbol_state(text: str, use_color: bool, android: bool) -> str:
    """Ver10.7J: android=True時は記号、use_color=True時はANSIカラー付き記号、両立対応
    Ver10.8J修正: 旧実装はandroid=True時にuse_colorの値を一切見ておらず、
    --color を強制指定してもPydroid3では絶対に色が付かないバグがあった。
    android時もuse_color=Trueなら記号にANSIカラーを重ねて返すよう修正。
    色判定は記号(例: [★強開★])ではなく元のフルテキスト(例: 強い開通)で行う
    （記号への部分一致では「強い開通」「弱い兆候」が拾えず無色になるため）。
    """
    sym = state_symbol(text)
    if not android:
        return color_by_state(short_state(text), use_color=use_color)
    if not use_color:
        return sym
    code = _color_code_for_state(text)
    return f"{code}{sym}{C_RESET}" if code else sym


def bar_graph(pct: float, width: int = 10) -> str:
    """Ver10.7J: 実効寄与をブロックバーで表示 (0-100% → width文字)
    █░ Block Elements (U+2588/U+2591) を使用。Pydroid3を含むAndroid環境でも
    デフォルトのNoto Sansフォントで対応済み（Ver11.0J: ASCII代替を廃止）。
    """
    filled = int(round(pct / 100.0 * width))
    filled = max(0, min(width, filled))
    return "█" * filled + "░" * (width - filled)


def analysis_state_label(stage: int) -> str:
    return {1: "簡易解析", 2: "標準解析", 3: "高度解析"}.get(stage, "標準解析")


def phase_label(cycles: int) -> str:
    if cycles < 3:
        return "初期観測"
    if cycles < MIN_ANALYSIS_POINTS:
        return "参考予測"
    return "標準解析"


def ratio_label(r) -> str:
    if r is None or pd.isna(r):
        return "不明"
    r = float(r)
    if r >= 0.70:
        return "強"
    if r >= 0.55:
        return "中"
    if r >= 0.40:
        return "弱"
    if r >= 0.25:
        return "微"
    return "低"


def trend_label(v) -> str:
    if pd.isna(v):
        return "不明"
    if v >= 2.0:
        return "急上昇"
    if v >= 0.5:
        return "上昇"
    if v <= -2.0:
        return "急低下"
    if v <= -0.5:
        return "低下"
    return "横ばい"


def f_region_phrase(masking_risk) -> str:
    if masking_risk is None or pd.isna(masking_risk):
        return "不明"
    x = float(masking_risk)
    if x < 0.2:
        return "小"
    if x < 0.5:
        return "中"
    return "大"


def layer_thickness_phrase(thickness) -> str:
    if thickness is None or pd.isna(thickness):
        return "不明"
    x = float(thickness)
    if x <= 1.0:
        return "薄め"
    if x <= 2.5:
        return "普通"
    return "厚め"


def blanketing_phrase(blanketing_index) -> str:
    if blanketing_index is None or pd.isna(blanketing_index):
        return "不明"
    return "あり" if float(blanketing_index) >= 0.5 else "なし"



def fetch_nict_once(url: str) -> str:
    headers = {"User-Agent": "CBM2.py Es predictor; one GET per cycle"}
    r = requests.get(url, timeout=HTTP_TIMEOUT_SEC, headers=headers)
    r.raise_for_status()
    r.encoding = r.apparent_encoding or "utf-8"
    return r.text


def parse_datetime_jst(date_s: str, time_s: str):
    s = f"{date_s} {time_s}".replace("年", "/").replace("月", "/").replace("日", " ")
    s = s.replace("JST", "").strip()
    s = re.sub(r"\s+", " ", s)
    for fmt_s in ["%Y/%m/%d %H:%M", "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M:%S", "%Y-%m-%d %H:%M:%S"]:
        try:
            return datetime.strptime(s, fmt_s).replace(tzinfo=JST)
        except ValueError:
            pass
    return pd.NaT


def html_to_text(html: str) -> str:
    soup = BeautifulSoup(html, "html.parser")
    return soup.get_text("\n", strip=True)


def normalize_text_for_tokens(text: str) -> str:
    t = text.replace("\r", "\n")
    t = t.replace("年", "/").replace("月", "/").replace("日", " ")
    t = t.replace("JST", " ")
    t = t.replace("\u3000", " ")
    return t


def tokenize_nict_text(text: str):
    t = normalize_text_for_tokens(text)
    return re.findall(
        r"20\d{2}[/-]\d{1,2}[/-]\d{1,2}|"
        r"\d{1,2}:\d{2}(?::\d{2})?|"
        r"[-+]?\d+(?:\.\d+)?|"
        r"----|---|--|Undefined|undefined|"
        r"Okinawa|Yamagawa|Kokubunji|Wakkanai|Ogimi|"
        r"OKINAWA|YAMAGAWA|KOKUBUNJI|WAKKANAI|OGIMI|"
        r"沖縄|山川|国分寺|稚内|大宜味|"
        r"Date|Time|fxEs|foEs|fbEs|hEs|h'Es|foF2|hmF2|MUF",
        t,
    )


def build_records_from_date_time_values(date_s: str, time_s: str, values: List[str]):
    obs = parse_datetime_jst(date_s, time_s)
    if pd.isna(obs):
        return []
    keys = ["okinawa", "yamagawa", "kokubunji", "wakkanai"]
    records = []
    for key, val in zip(keys, values[:4]):
        fx = safe_float(val)
        if not pd.isna(fx):
            records.append({"obs_jst": obs.isoformat(), "station": key, "fxEs": fx})
    return records


def parse_fxes_from_tokens(tokens: List[str]) -> pd.DataFrame:
    records = []
    i = 0
    date_re = re.compile(r"^20\d{2}[/-]\d{1,2}[/-]\d{1,2}$")
    time_re = re.compile(r"^\d{1,2}:\d{2}(?::\d{2})?$")
    while i < len(tokens) - 5:
        if date_re.match(tokens[i]) and time_re.match(tokens[i + 1]):
            vals = []
            j = i + 2
            while j < len(tokens) and len(vals) < 4:
                if date_re.match(tokens[j]):
                    break
                if re.fullmatch(r"[-+]?\d+(?:\.\d+)?|----|---|--|Undefined|undefined|\*", str(tokens[j]).strip()):
                    vals.append(tokens[j])
                j += 1
            if len(vals) >= 4:
                records.extend(build_records_from_date_time_values(tokens[i], tokens[i + 1], vals))
                i = j
                continue
        i += 1
    out = pd.DataFrame(records)
    if out.empty:
        return out
    return out.drop_duplicates(["obs_jst", "station"]).sort_values(["obs_jst", "station"])


def parse_fxes_from_line_regex(text: str) -> pd.DataFrame:
    records = []
    normalized = normalize_text_for_tokens(text)
    lines = [re.sub(r"\s+", " ", ln).strip() for ln in normalized.split("\n")]
    row_re = re.compile(
        r"(20\d{2}[/-]\d{1,2}[/-]\d{1,2})\s+"
        r"(\d{1,2}:\d{2}(?::\d{2})?)\s+"
        r"([+\-]?\d+(?:\.\d+)?|----|---|--|Undefined|undefined|\*)\s+"
        r"([+\-]?\d+(?:\.\d+)?|----|---|--|Undefined|undefined|\*)\s+"
        r"([+\-]?\d+(?:\.\d+)?|----|---|--|Undefined|undefined|\*)\s+"
        r"([+\-]?\d+(?:\.\d+)?|----|---|--|Undefined|undefined|\*)"
    )
    for line in lines:
        m = row_re.search(line)
        if m:
            records.extend(build_records_from_date_time_values(m.group(1), m.group(2), [m.group(3), m.group(4), m.group(5), m.group(6)]))
    out = pd.DataFrame(records)
    if out.empty:
        return out
    return out.drop_duplicates(["obs_jst", "station"]).sort_values(["obs_jst", "station"])


def parse_fxes_from_tables(html: str) -> pd.DataFrame:
    def cell_text(cell):
        return re.sub(r"\s+", " ", cell.get_text(" ", strip=True)).strip()
    soup = BeautifulSoup(html, "html.parser")
    tables = soup.find_all("table")
    if not tables:
        return pd.DataFrame()
    station_keywords = ["Okinawa", "Yamagawa", "Kokubunji", "Wakkanai", "沖縄", "山川", "国分寺", "稚内"]
    best_rows = None
    best_score = -1
    for table in tables:
        rows = []
        for tr in table.find_all("tr"):
            cells = tr.find_all(["th", "td"])
            if cells:
                rows.append([cell_text(c) for c in cells])
        if not rows:
            continue
        joined = " ".join(" ".join(r) for r in rows[:10])
        score = sum(1 for k in station_keywords if k in joined)
        if "Date" in joined or "日" in joined:
            score += 1
        if "Time" in joined or "時" in joined:
            score += 1
        if score > best_score:
            best_score = score
            best_rows = rows
    if best_rows is None or best_score <= 0:
        return pd.DataFrame()
    text = "\n".join("\n".join(r) for r in best_rows)
    return parse_fxes_from_tokens(tokenize_nict_text(text))


def parse_fxes_html(html: str, debug: bool = False,
                    debug_html: Path = DEFAULT_DEBUG_HTML,
                    debug_text: Path = DEFAULT_DEBUG_TEXT,
                    debug_tokens: Path = DEFAULT_DEBUG_TOKENS) -> pd.DataFrame:
    text = html_to_text(html)
    if debug:
        debug_html.write_text(html, encoding="utf-8", errors="ignore")
        debug_text.write_text(text, encoding="utf-8", errors="ignore")
    out = parse_fxes_from_tables(html)
    if not out.empty:
        return out
    out = parse_fxes_from_line_regex(text)
    if not out.empty:
        return out
    tokens = tokenize_nict_text(text)
    if debug:
        debug_tokens.write_text("\n".join(tokens[:3000]), encoding="utf-8", errors="ignore")
    out = parse_fxes_from_tokens(tokens)
    if not out.empty:
        return out
    rough = re.sub(r"<[^>]+>", " ", html)
    rough = re.sub(r"\s+", " ", rough)
    tokens = tokenize_nict_text(rough)
    out = parse_fxes_from_tokens(tokens)
    if not out.empty:
        return out
    debug_html.write_text(html, encoding="utf-8", errors="ignore")
    debug_text.write_text(text, encoding="utf-8", errors="ignore")
    debug_tokens.write_text("\n".join(tokens[:3000]), encoding="utf-8", errors="ignore")
    raise RuntimeError("fxEs値が抽出できませんでした。NICTの構造変更の可能性があります。")


def latest_observation_only(parsed_df: pd.DataFrame) -> pd.DataFrame:
    if parsed_df is None or parsed_df.empty:
        return pd.DataFrame(columns=["obs_jst", "station", "fxEs"])
    df = parsed_df.copy()
    df["obs_dt_tmp"] = pd.to_datetime(df["obs_jst"], errors="coerce")
    latest_time = df["obs_dt_tmp"].max()
    latest = df[df["obs_dt_tmp"] == latest_time].drop(columns=["obs_dt_tmp"])
    return latest.drop_duplicates(["obs_jst", "station"]).sort_values(["obs_jst", "station"])


def empty_history() -> pd.DataFrame:
    return pd.DataFrame(columns=["obs_jst", "station", "fxEs"])


def merge_history_in_memory(history: pd.DataFrame, latest_sample: pd.DataFrame):
    if history is None or history.empty:
        before_keys = set()
        merged = latest_sample.copy()
    else:
        before_keys = set(zip(history["obs_jst"], history["station"]))
        merged = pd.concat([history, latest_sample], ignore_index=True)
    merged = merged.drop_duplicates(["obs_jst", "station"]).sort_values(["obs_jst", "station"])
    merged = merged.groupby("station", group_keys=False).tail(MAX_HISTORY_ROWS).reset_index(drop=True)
    after_keys = set(zip(merged["obs_jst"], merged["station"]))
    return merged, len(after_keys - before_keys)


SESSION_ARCHIVE_RETENTION_HOURS = 6


def archive_old_session_logs(log_dir: Path) -> Tuple[int, int]:
    """Ver13.15: 起動のたびに増え続けるセッション別CSV
    (prediction_YYYYMMDD_HHMMSS.csv / prediction_calibration_YYYYMMDD_HHMMSS.csv)を、
    累積アーカイブファイル(prediction_archive.csv / prediction_calibration_archive.csv)
    へ安全に統合してから元ファイルを削除する。

    【過去の学習データは失われない】
    auto_calibrate_profile() は list_session_csvs()/load_all_session_logs() で
    log_dir 内の "prediction_*.csv" 系ファイルを全てglobして学習データとして使う。
    アーカイブファイル名(prediction_archive.csv 等)もこのglobパターンに一致するため、
    元ファイルを削除してもアーカイブ側に同じ行データがそのまま残り、
    学習に使えるデータ量は変化しない。ファイル「個数」だけが減る。

    【安全策】
      1) 直近 SESSION_ARCHIVE_RETENTION_HOURS 時間以内に更新されたセッションCSVは
         対象外にする（実行中プロセス・他端末との同時書き込み競合を避けるため）。
      2) 読み込みに失敗した（壊れている）CSVはアーカイブせず、削除もしない。
      3) アーカイブへの追記(_csv_append、ロックファイルによる排他制御つき原子追記)が
         成功したことを確認してから、初めて元ファイルを削除する。
         追記に失敗した場合は元ファイルを一切消さない。
    分析専用レポート(*_corr.csv, *_summary.txt, partial_corr_*.csv)は
    学習の重みに使われない解析結果（コード内コメント「分析専用」参照）なので、
    アーカイブはせず、対応するセッションCSVがアーカイブされたタイミングで
    一緒に削除する。

    戻り値: (アーカイブに統合したセッションCSV数, 削除した分析専用レポート数)
    """
    archived = 0
    reports_removed = 0
    threshold = time.time() - SESSION_ARCHIVE_RETENTION_HOURS * 3600.0

    for prefix in ("prediction", "prediction_calibration"):
        archive_path = log_dir / f"{prefix}_archive.csv"
        pattern = f"{prefix}_[0-9]*.csv"
        for p in sorted(log_dir.glob(pattern)):
            if not p.is_file() or p.name == archive_path.name or p.name.endswith("_corr.csv"):
                continue
            try:
                if p.stat().st_mtime >= threshold:
                    continue
            except OSError:
                continue

            try:
                df = pd.read_csv(p)
            except Exception:
                continue
            if df.empty:
                continue

            try:
                _csv_append(archive_path, df)
            except Exception:
                continue

            try:
                p.unlink()
                archived += 1
            except Exception:
                continue

            session_id = p.stem[len(prefix) + 1:] if p.stem.startswith(prefix + "_") else None
            if session_id:
                for report_name in (f"{prefix}_{session_id}_corr.csv",
                                     f"{prefix}_{session_id}_summary.txt",
                                     f"partial_corr_{session_id}.csv"):
                    rp = log_dir / report_name
                    if rp.exists():
                        try:
                            rp.unlink()
                            reports_removed += 1
                        except Exception:
                            pass

    return archived, reports_removed


def list_session_csvs(log_dir: Path) -> List[Path]:
    paths = set(log_dir.glob("prediction_*.csv")) | set(log_dir.glob("prediction_calibration_*.csv"))
    return sorted(p for p in paths if p.is_file() and not p.name.endswith("_corr.csv"))


def load_all_session_logs(log_dir: Path) -> pd.DataFrame:
    frames = []
    for p in list_session_csvs(log_dir):
        try:
            frames.append(pd.read_csv(p))
        except Exception:
            pass
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def split_into_segments(ds: pd.DataFrame) -> List[pd.DataFrame]:
    """obs_dt昇順のDataFrameを、Δt(分)が10〜20の範囲を外れた時点で
    新規セグメントに分割する。各セグメントはobs_dt連続・等間隔観測のみを含む。
    """
    ds = ds.sort_values("obs_dt").reset_index(drop=True)
    if ds.empty:
        return []
    dt_min = ds["obs_dt"].diff().dt.total_seconds() / 60.0
    segments: List[pd.DataFrame] = []
    start = 0
    for i in range(1, len(ds)):
        d = dt_min.iloc[i]
        if pd.isna(d) or d > 20.0 or d < 10.0:
            segments.append(ds.iloc[start:i].reset_index(drop=True))
            start = i
    segments.append(ds.iloc[start:].reset_index(drop=True))
    return segments


def compute_dt_quality(df: pd.DataFrame) -> Dict[str, float]:
    """観測局ごとのΔt(分)から等間隔率を算出する。
    equal_ratio = 有効Δt数(10〜20分) / 全Δt数
    """
    valid_count = 0
    total_count = 0
    for _, ds in df.groupby("station"):
        ds = ds.sort_values("obs_dt")
        dt_min = ds["obs_dt"].diff().dt.total_seconds() / 60.0
        dt_min = dt_min.dropna()
        total_count += int(len(dt_min))
        valid_count += int(((dt_min >= 10.0) & (dt_min <= 20.0)).sum())
    equal_ratio = (valid_count / total_count) if total_count > 0 else 0.0
    return {"valid_dt": valid_count, "total_dt": total_count, "equal_ratio": equal_ratio}


def observation_quality_label(equal_ratio: float) -> str:
    if equal_ratio >= 0.8:
        return "良好"
    if equal_ratio >= 0.5:
        return "普通"
    return "不十分"


def estimate_cycles_to_good_quality(
    valid_dt: int, total_dt: int, equal_ratio: float, interval_sec: int = DEFAULT_INTERVAL_SEC,
) -> Tuple[int, float]:
    """等間隔率(equal_ratio)を「良好」(>=0.8)まで引き上げるのに必要な
    追加の等間隔観測回数の目安を返す。

    以後の観測が全て等間隔(10〜20分)で継続すると仮定した場合の下限値:
      (valid_dt + X) / (total_dt + X) >= 0.8  を解くと
      X >= 4 * (0.8*total_dt - valid_dt)
    戻り値: (必要回数X, 目安時間(分))
    """
    if equal_ratio >= 0.8 or total_dt <= 0:
        return 0, 0.0
    needed = 4.0 * (0.8 * total_dt - valid_dt)
    x = max(0, int(np.ceil(needed)))
    minutes = x * (interval_sec / 60.0)
    return x, minutes


def format_quality_progress_hint(label: str, equal_ratio: float, x_cycles: int, minutes: float) -> str:
    """画面表示用の「良好まで」小型ヒント文字列を作る（小さく1行で収まる形式）。"""
    if label == "良好":
        return "品質: 良好 到達済み"
    if x_cycles <= 0:
        return f"品質: {label}({equal_ratio*100:.0f}%)"
    if minutes >= 120:
        time_str = f"約{minutes/60:.1f}時間"
    else:
        time_str = f"約{minutes:.0f}分"
    return f"品質: {label}({equal_ratio*100:.0f}%) 良好まであと{x_cycles}回({time_str})"


FUTURE_LAG_COLUMNS = {"future15": 1, "future30": 2, "future45": 3, "future60": 4, "future90": 6}


def add_future_deltas_segmented(segments: List[pd.DataFrame], col: str = "fxEs") -> pd.DataFrame:
    """有効セグメントのリストを受け取り、各セグメント内のみでshiftベースの
    将来変化量(future15〜future90)を計算して結合したDataFrameを返す。
    セグメント境界を越えた未来値参照は発生しない（shiftはセグメント内のみ）。
    """
    out_frames = []
    for seg in segments:
        s = seg.reset_index(drop=True).copy()
        fx = pd.to_numeric(s[col], errors="coerce")
        for name, lag in FUTURE_LAG_COLUMNS.items():
            s[name] = fx.shift(-lag) - fx
        out_frames.append(s)
    if not out_frames:
        return pd.DataFrame()
    return pd.concat(out_frames, ignore_index=True)


def build_valid_segment_dataset(
    df: pd.DataFrame,
    min_segment_len: int = 3,
    min_duration_min: float = MIN_SEGMENT_DURATION_MIN,
) -> Tuple[pd.DataFrame, int, int]:
    """観測局ごとにセグメント分割し、サンプル数 >= min_segment_len かつ
    経過時間 >= min_duration_min分 のセグメントのみを採用して
    将来変化量を付与した結合データセットを返す。

    Ver11.5: 点数条件だけでは --interval を短く変更した環境で
    「点数は足りるが実際は30分未満」のスポット起動が学習対象に
    混入し得たため、経過時間(obs_dt最終行-先頭行)による判定を独立に追加した。
    どちらか一方でも満たさないセグメントは除外し、除外件数を返す。

    戻り値: (有効セグメント結合データ, 有効セグメント数, 除外セグメント数)
    """
    frames = []
    seg_count = 0
    rejected_count = 0
    for _, ds in df.groupby("station"):
        ds = ds.sort_values("obs_dt").reset_index(drop=True)
        segments = split_into_segments(ds)
        valid_segs = []
        for s in segments:
            if len(s) < min_segment_len:
                rejected_count += 1
                continue
            span_min = (s["obs_dt"].iloc[-1] - s["obs_dt"].iloc[0]).total_seconds() / 60.0
            if span_min < min_duration_min:
                rejected_count += 1
                continue
            valid_segs.append(s)
        seg_count += len(valid_segs)
        if valid_segs:
            frames.append(add_future_deltas_segmented(valid_segs))
    if not frames:
        return pd.DataFrame(), 0, rejected_count
    return pd.concat(frames, ignore_index=True), seg_count, rejected_count


def compute_vif(df_feat: pd.DataFrame) -> Dict[str, float]:
    """Ver10.6J: VIF(分散膨張係数)を算出する。
    各特徴量を他の特徴量で線形回帰し、決定係数R²からVIF=1/(1-R²)を計算。
    返り値: {特徴量名: VIF値}。行数不足または演算失敗時はNaNを返す。
    VIF目安: <5=問題なし / 5〜10=注意 / >10=多重共線性の可能性
    注: 多重共線性は予測精度を必ずしも下げないが重みの不安定化を招く。
    """
    result: Dict[str, float] = {}
    cols = [c for c in df_feat.columns if df_feat[c].notna().sum() >= 10]
    if len(cols) < 2:
        return {c: np.nan for c in df_feat.columns}
    data = df_feat[cols].dropna()
    if len(data) < 10:
        return {c: np.nan for c in df_feat.columns}
    for col in cols:
        y = data[col].values
        X = data[[c for c in cols if c != col]].values
        try:
            coef, res, rank, sv = np.linalg.lstsq(
                np.column_stack([np.ones(len(X)), X]), y, rcond=None
            )
            y_pred = np.column_stack([np.ones(len(X)), X]) @ coef
            ss_res = float(np.sum((y - y_pred) ** 2))
            ss_tot = float(np.sum((y - y.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
            r2 = float(np.clip(r2, 0.0, 0.9999))
            result[col] = float(1.0 / (1.0 - r2))
        except Exception:
            result[col] = np.nan
    for col in df_feat.columns:
        if col not in result:
            result[col] = np.nan
    return result


_F_FEAT_RIDGE = [
    ("foF2_es_ratio", "w_fof2"),
    ("dfoF2",         "w_dfof2"),
    ("hmF2",          "w_hmf2"),
    ("dhmF2",         "w_dhmf2"),
    ("m3000f2_norm",  "w_m3000f2"),
    ("rpi_norm",      "w_rpi"),
    ("es_mode_score",  "w_es_mode"),
    ("eqi_norm",      "w_eqi"),
    ("sw_score",      "w_sw_score"),
    ("es_phase_num",  "w_es_phase"),
    ("ia_eqi_fxes",   "w_ia_eqi_fxes"),
    ("ia_eqi_mstid",  "w_ia_eqi_mstid"),
    ("ia_sw_fxes",    "w_ia_sw_fxes"),
    ("ia_sw_mstid",   "w_ia_sw_mstid"),
    ("ia_eqi_sw",     "w_ia_eqi_sw"),
]
_F_FEAT_FALLBACK = "foF2_ratio"
RIDGE_ALPHA = 1.0
RIDGE_WEIGHT_CAP = 0.20
LAG_STABILITY_MIN_THRESHOLD = 3.0
VIF_EXCLUDE_THR = 20.0

OPEN_RPI_THRESHOLD      = 20.0

EDFS_LV1_TOTAL          = 100;   EDFS_LV1_ES       = 30
EDFS_LV2_AREA_TEACHERS  = 500
EDFS_LV3_CENTROID_HIST  = 150
EDFS_LV4_LIFETIME       = 300
EDFS_LV5_AREA_TEACHERS  = 3000;  EDFS_LV5_CENTROID = 2000;  EDFS_LV5_LIFETIME = 1000

FEAT_AREA_MIN           = 500
FEAT_MOTION_MIN         = 500
FEAT_ADV_MOTION_MIN     = 1000
FEAT_DYNAMICS_AREA      = 3000;  FEAT_DYNAMICS_CENT = 2000;  FEAT_DYNAMICS_LIFE = 1000

EDFS_NORMAL_TRAIN_INTERVAL_CYCLES = 10
EDFS_MATURE_TRAIN_INTERVAL_CYCLES = 40

EDFS_DRIFT_CHECK_INTERVAL_CYCLES  = 20
EDFS_DRIFT_LOCATION_SIGMA         = 3.0
EDFS_DRIFT_SPREAD_RATIO           = 2.2
EDFS_DRIFT_SFI_RATIO              = 1.4
EDFS_DRIFT_SFI_ABS_MIN            = 15.0
EDFS_DRIFT_SFI_SMOOTH_N           = 3
EDFS_DRIFT_SNAPSHOT_WINDOW        = 200
EDFS_DRIFT_MIN_SAMPLES            = 50
EDFS_DRIFT_CONFIRM_COUNT          = 2
EDFS_DRIFT_COOLDOWN_CYCLES        = 100
EDFS_TRAINING_STATE_FILE          = "edfs_training_state.json"

FUTURE_LAGS = {15: 1, 30: 2, 45: 3, 60: 4, 90: 6}

EDFS_AREA_HISTORY_FILE      = "edfs_area_history.csv"
EDFS_CENTROID_HISTORY_FILE  = "edfs_centroid_history.csv"
EDFS_PREDICTION_ERROR_FILE  = "edfs_prediction_error.csv"
EDFS_RIDGE_TRAINING_FILE    = "edfs_ridge_training.csv"

EDFS_MERGE_INTERVAL_SEC   = 900
EDFS_LOCK_TIMEOUT_SEC     = 5.0
EDFS_LOCK_POLL_SEC        = 0.05
EDFS_V12_MIGRATION_MARKER = ".edfs_v12_migrated"
MERGE_ACCELERATOR_FILENAME = "es_data_merge_accelerator.py"


def _ridge_calibrate_f_features(
    valid_df: pd.DataFrame,
    profile: Dict,
    future_cols: List[str],
) -> Tuple[List[str], bool]:
    """Ver10.6J: Ridge回帰でF層特徴量の重みを学習する。
    全future lagについてRidgeを実行し、最良R²のlagの係数を採用。
    係数は正規化スケール上の値をそのまま用い、RIDGE_WEIGHT_CAPでクリップする。

    戻り値: (cal_log行リスト, ridge_used_flag)
    sklearn未インストール時はFalseを返し、呼び出し元でフォールバックする。
    """
    if not _SKLEARN_AVAILABLE:
        return ["Ridge: sklearn未インストール → 単純相関にフォールバック"], False

    feat_pairs = []
    for feat_col, wkey in _F_FEAT_RIDGE:
        if feat_col in valid_df.columns:
            feat_pairs.append((feat_col, wkey))
        elif feat_col == "foF2_es_ratio" and _F_FEAT_FALLBACK in valid_df.columns:
            feat_pairs.append((_F_FEAT_FALLBACK, wkey))

    if not feat_pairs:
        return ["Ridge: F層特徴量列なし → 重みすべて0"], True

    feat_names = [fp[0] for fp in feat_pairs]
    wkeys      = [fp[1] for fp in feat_pairs]

    vif_excluded: List[str] = []
    cal_log_vif: List[str] = []
    vif_check_df = valid_df[feat_names].apply(pd.to_numeric, errors="coerce")
    vif_for_ridge = compute_vif(vif_check_df)
    high_vif = {k: v for k, v in vif_for_ridge.items() if not pd.isna(v) and v > VIF_EXCLUDE_THR}
    if high_vif:
        cal_log_vif.append(
            f"Ridge: VIF>{VIF_EXCLUDE_THR:.0f}で動的排除: "
            + ", ".join(f"{k}(VIF={v:.1f})" for k, v in high_vif.items())
        )
        feat_pairs = [(f, w) for f, w in feat_pairs if f not in high_vif]
        vif_excluded = list(high_vif.keys())
        feat_names = [fp[0] for fp in feat_pairs]
        wkeys      = [fp[1] for fp in feat_pairs]
        if not feat_pairs:
            cal_log_vif.append("Ridge: 全特徴量がVIF排除 → 重みすべて0")
            return cal_log_vif, True

    cal_log = cal_log_vif + [f"Ridge(alpha={RIDGE_ALPHA}): 説明変数={feat_names}"]
    if vif_excluded:
        cal_log.append(f"  VIF排除済(重み=0維持): {vif_excluded}")

    best_r2    = -np.inf
    best_coefs: np.ndarray = np.zeros(len(feat_names))
    best_lag   = None
    best_scaler_mean: Optional[np.ndarray] = None
    best_scaler_scale: Optional[np.ndarray] = None
    lag_r2_log = []
    r2_by_lag: Dict[str, float] = {}

    for fc in future_cols:
        if fc not in valid_df.columns:
            continue
        sub = valid_df[feat_names + [fc]].copy()
        sub = sub.apply(pd.to_numeric, errors="coerce").dropna()
        if len(sub) < 20:
            lag_r2_log.append(f"{fc}:n={len(sub)}(不足)")
            continue

        X_raw = sub[feat_names].values
        y_raw = sub[fc].values

        try:
            scaler = StandardScaler()
            X_sc = scaler.fit_transform(X_raw)
            reg = Ridge(alpha=RIDGE_ALPHA, fit_intercept=True)
            reg.fit(X_sc, y_raw)
            y_pred = reg.predict(X_sc)
            ss_res = float(np.sum((y_raw - y_pred) ** 2))
            ss_tot = float(np.sum((y_raw - y_raw.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0
            lag_r2_log.append(f"{fc}:R²={r2:.3f}(n={len(sub)})")
            r2_by_lag[fc] = r2
            if r2 > best_r2:
                best_r2  = r2
                best_coefs = reg.coef_.copy()
                best_lag = fc
                best_scaler_mean = scaler.mean_.copy()
                best_scaler_scale = scaler.scale_.copy()
        except Exception as e:
            lag_r2_log.append(f"{fc}:error({e})")

    cal_log.append(f"  lag別R²: {', '.join(lag_r2_log)}")

    r2_vals = list(r2_by_lag.values())
    if len(r2_vals) >= 2:
        r2_std = float(np.std(r2_vals))
        lag_stability_score = float(1.0 / r2_std) if r2_std > 1e-6 else 999.0
        cal_log.append(f"  lag_stability_score={lag_stability_score:.2f} (n_lag={len(r2_vals)})")
    else:
        lag_stability_score = None
        cal_log.append("  lag_stability_score=算出不可(有効lag<2)")
    profile["lag_stability_score"] = lag_stability_score

    if best_lag is None or best_r2 < 0.0:
        for wkey in wkeys:
            profile[wkey] = 0.0
        for _, wkey in _F_FEAT_RIDGE:
            profile["sign_" + wkey[2:]] = 0.0
        profile["ridge_lag_stable"] = None
        cal_log.append(f"  最良R²={best_r2:.3f} < 0 → 全重み=0")
        return cal_log, True

    if lag_stability_score is not None and lag_stability_score < LAG_STABILITY_MIN_THRESHOLD:
        for wkey in wkeys:
            profile[wkey] = 0.0
        for _, wkey in _F_FEAT_RIDGE:
            profile["sign_" + wkey[2:]] = 0.0
        profile["ridge_lag_stable"] = False
        cal_log.append(
            f"  採用候補lag={best_lag} R²={best_r2:.3f} だが "
            f"lag_stability_score={lag_stability_score:.2f} < {LAG_STABILITY_MIN_THRESHOLD} "
            f"→ 偶然相関とみなし全重み=0"
        )
        return cal_log, True
    profile["ridge_lag_stable"] = True

    cal_log.append(f"  採用lag={best_lag} R²={best_r2:.3f}")

    f_scaler_mean = dict(profile.get("f_scaler_mean") or {})
    f_scaler_scale = dict(profile.get("f_scaler_scale") or {})
    if best_scaler_mean is not None and best_scaler_scale is not None:
        for i, (feat_col, _) in enumerate(feat_pairs):
            f_scaler_mean[feat_col] = float(best_scaler_mean[i])
            f_scaler_scale[feat_col] = float(best_scaler_scale[i]) if best_scaler_scale[i] > 1e-9 else 1.0
        cal_log.append(
            "  scaler統計: " + ", ".join(
                f"{fc}(mean={f_scaler_mean[fc]:.3f},scale={f_scaler_scale[fc]:.3f})"
                for fc, _ in feat_pairs
            )
        )
    profile["f_scaler_mean"] = f_scaler_mean
    profile["f_scaler_scale"] = f_scaler_scale

    for i, (feat_col, wkey) in enumerate(feat_pairs):
        raw_coef = float(best_coefs[i])
        w = float(np.clip(abs(raw_coef), 0.0, RIDGE_WEIGHT_CAP))
        sign = float(np.sign(raw_coef))
        sign_key = "sign_" + wkey[2:]
        profile[wkey] = w
        profile[sign_key] = sign
        cal_log.append(f"  {feat_col}: coef={raw_coef:.4f}(sign={sign:+.0f}) → {wkey}={w:.4f}")

    used_wkeys = set(wkeys)
    for feat_col, wkey in _F_FEAT_RIDGE:
        if wkey not in used_wkeys:
            profile[wkey] = 0.0
            profile["sign_" + wkey[2:]] = 0.0
            if feat_col in vif_excluded:
                cal_log.append(f"  {feat_col}: VIF排除 → {wkey}=0.0")

    return cal_log, True


def auto_calibrate_profile(log_dir: Path, profile_path: Path, min_rows: int = MIN_AUTO_CAL_ROWS) -> Tuple[Dict, str]:
    """Ver10.3J: 不等間隔ログ自動補正除外機能。
    Pydroid3スポット運用やRaspberryPi常時監視のログが混在しても、
    10〜20分の等間隔で連続するセグメント（サンプル数3以上）のみを
    補正学習対象とする。
    """
    profile = dict(BASE_PROFILE)
    if profile_path.exists():
        try:
            profile.update(json.loads(profile_path.read_text(encoding='utf-8')))
        except Exception:
            pass

    df = load_all_session_logs(log_dir)
    if df.empty or len(df) < min_rows:
        profile["source"] = "default_or_previous"
        shortfall = max(0, min_rows - len(df))
        msg_lines = [
            "自動補正: 未実施",
            "理由: ログデータ不足",
            f"有効データ: {len(df)}件",
            f"必要データ: {min_rows}件",
            f"不足: あと{shortfall}件",
        ]
        return profile, "\n".join(msg_lines)

    required_cols = {"datetime_jst", "station", "fxEs"}
    if not required_cols.issubset(set(df.columns)):
        profile["source"] = "default_or_previous"
        return profile, "自動補正: 未実施\n理由: 必要列が見つかりません"

    df = df.copy()
    df["obs_dt"] = pd.to_datetime(df["datetime_jst"], errors="coerce")
    df = df.dropna(subset=["obs_dt", "station"])
    if df.empty:
        profile["source"] = "default_or_previous"
        return profile, "自動補正: 未実施\n理由: 観測時刻を解釈できませんでした"

    quality = compute_dt_quality(df)
    equal_ratio = quality["equal_ratio"]
    quality_label = observation_quality_label(equal_ratio)
    _good_x, _good_min = estimate_cycles_to_good_quality(
        quality["valid_dt"], quality["total_dt"], equal_ratio,
    )
    quality_hint = format_quality_progress_hint(quality_label, equal_ratio, _good_x, _good_min)
    profile["_quality_hint"] = quality_hint

    valid_df, seg_count, rejected_seg_count = build_valid_segment_dataset(
        df, min_segment_len=3, min_duration_min=MIN_SEGMENT_DURATION_MIN,
    )
    valid_data_count = int(len(valid_df))

    quality_lines = [
        f"観測品質: {quality_label}",
        f"等間隔率: {equal_ratio*100:.0f}%",
        f"有効データ: {valid_data_count}件",
        f"除外セグメント: {rejected_seg_count}件"
        f"（{MIN_SEGMENT_DURATION_MIN:.0f}分未満のスポット起動など）",
        quality_hint,
    ]

    if not (seg_count >= 1 and valid_data_count >= MIN_AUTO_CAL_ROWS and equal_ratio >= 0.5):
        profile["source"] = "default_or_previous"
        shortfall = max(0, MIN_AUTO_CAL_ROWS - valid_data_count)
        msg_lines = [
            "自動補正: 未実施",
            "理由: 等間隔観測不足、または最低起動時間未達",
            f"有効データ: {valid_data_count}件",
            f"必要データ: {MIN_AUTO_CAL_ROWS}件",
            f"不足: あと{shortfall}件",
            f"等間隔率: {equal_ratio*100:.0f}%",
            f"除外セグメント: {rejected_seg_count}件"
            f"（1回の起動は{MIN_SEGMENT_DURATION_MIN:.0f}分以上継続してください）",
        ]
        if shortfall > 0:
            msg_lines.append(
                f"目安: {MIN_SEGMENT_DURATION_MIN:.0f}分起動をあと{shortfall}回"
                f"（1回あたり最短でも1件分。長め起動ならより早く到達）"
            )
        msg_lines.append(quality_hint)
        return profile, "\n".join(msg_lines)

    corr_rows = []
    if "mstid_approach_index" in valid_df.columns:
        x_all = pd.to_numeric(valid_df["mstid_approach_index"], errors="coerce")
        for future_col, lag in FUTURE_LAG_COLUMNS.items():
            y_all = pd.to_numeric(valid_df.get(future_col), errors="coerce") if future_col in valid_df.columns else pd.Series(dtype=float)
            valid_mask = x_all.notna() & y_all.notna()
            if int(valid_mask.sum()) >= 20:
                corr = float(np.corrcoef(x_all[valid_mask], y_all[valid_mask])[0, 1])
                corr_rows.append((lag * 15, corr, int(valid_mask.sum())))

    if corr_rows:
        best = sorted(corr_rows, key=lambda r: abs(r[1]) if not pd.isna(r[1]) else -1, reverse=True)[0]
        lag_min, best_corr, n_used = best
        if (not pd.isna(best_corr)) and best_corr > 0.10:
            profile["w_mstid_approach"] = float(min(0.15, max(0.0, best_corr) * 0.20))
        else:
            profile["w_mstid_approach"] = 0.0
        profile["mstid_best_lag_min"] = int(lag_min)
        profile["mstid_best_corr"] = None if pd.isna(best_corr) else float(best_corr)
        profile["mstid_best_n"] = int(n_used)
    else:
        profile["w_mstid_approach"] = 0.0

    vif_feat_cols = [c for c in [
        "fxEs", "trend", "cluster", "foF2_es_ratio", "foF2_ratio",
        "dfoF2", "hmF2", "dhmF2", "m3000f2_norm", "mstid_approach_index",
        "fxes_std_60", "d2_eff",
        "rpi_norm", "es_mode_score",
    ] if c in valid_df.columns]
    vif_df_input = valid_df[vif_feat_cols].apply(pd.to_numeric, errors="coerce")
    vif_result = compute_vif(vif_df_input)
    vif_warnings = [f"{k}:VIF={v:.1f}(>10)" for k, v in vif_result.items()
                    if not pd.isna(v) and v > 10.0]
    profile["vif"] = {k: (None if pd.isna(v) else round(float(v), 2))
                      for k, v in vif_result.items()}
    if vif_warnings:
        profile["vif_warnings"] = vif_warnings

    future_cols_list = list(FUTURE_LAG_COLUMNS.keys())
    f_cal_log, ridge_used = _ridge_calibrate_f_features(valid_df, profile, future_cols_list)

    if not ridge_used:
        f_feat_map = {
            "foF2_es_ratio": "w_fof2",
            "dfoF2":         "w_dfof2",
            "hmF2":          "w_hmf2",
            "dhmF2":         "w_dhmf2",
            "m3000f2_norm":  "w_m3000f2",
            "rpi_norm":      "w_rpi",
            "es_mode_score": "w_es_mode",
        }
        LAG_CORR_THR = 0.10
        LAG_CONSISTENCY_MIN = 2
        for feat_col, weight_key in f_feat_map.items():
            if feat_col not in valid_df.columns:
                if feat_col == "foF2_es_ratio" and "foF2_ratio" in valid_df.columns:
                    feat_col_use = "foF2_ratio"
                else:
                    profile[weight_key] = 0.0
                    f_cal_log.append(f"{feat_col}: 列なし → 重み=0")
                    continue
            else:
                feat_col_use = feat_col
            x_f = pd.to_numeric(valid_df[feat_col_use], errors="coerce")
            lag_corrs = []
            lag_results = []
            for fc in future_cols_list:
                if fc not in valid_df.columns:
                    continue
                y_f = pd.to_numeric(valid_df[fc], errors="coerce")
                vm = x_f.notna() & y_f.notna()
                n_valid = int(vm.sum())
                if n_valid < 20:
                    lag_results.append(f"{fc}:n={n_valid}(不足)")
                    continue
                c = float(np.corrcoef(x_f[vm], y_f[vm])[0, 1])
                lag_corrs.append((fc, c, n_valid))
                lag_results.append(f"{fc}:corr={c:.3f}(n={n_valid})")
            above_thr = [(fc, c, n) for fc, c, n in lag_corrs if abs(c) > LAG_CORR_THR]
            if len(above_thr) < LAG_CONSISTENCY_MIN:
                profile[weight_key] = 0.0
                f_cal_log.append(
                    f"{feat_col_use}: 閾値超えlag={len(above_thr)}本<{LAG_CONSISTENCY_MIN}本"
                    f" → ノイズ判定 重み=0  [{', '.join(lag_results)}]"
                )
                continue
            best_fc, best_corr_f, best_n = max(above_thr, key=lambda r: abs(r[1]))
            w = float(min(0.20, abs(best_corr_f) * 0.20))
            profile[weight_key] = w
            f_cal_log.append(
                f"{feat_col_use}: best={best_fc} corr={best_corr_f:.3f} n={best_n}"
                f" (閾値超え{len(above_thr)}本) → {weight_key}={w:.4f}  [{', '.join(lag_results)}]"
            )

    profile["f_layer_cal_log"] = f_cal_log
    profile["ridge_used"] = ridge_used

    if ridge_used and len(valid_df) >= 30:
        perm_exclude = permutation_importance_check(profile, valid_df)
        if perm_exclude:
            f_cal_log.append(
                f"Ver11.4自己学習: PermImportance低寄与除外: {perm_exclude}"
            )
            for feat in perm_exclude:
                wkey = next((w for f, w in _F_FEAT_RIDGE if f == feat), None)
                if wkey:
                    profile[wkey] = 0.0
                    profile["sign_" + wkey[2:]] = 0.0
                    f_cal_log.append(f"  {feat} → {wkey}=0.0 (自己学習除外)")
        profile["perm_exclude_features"] = perm_exclude
        profile["f_layer_cal_log"] = f_cal_log

    if "guide_percent" in valid_df.columns and "future30" in valid_df.columns:
        gp = pd.to_numeric(valid_df["guide_percent"], errors="coerce")
        future30 = pd.to_numeric(valid_df["future30"], errors="coerce")
        valid_mask = gp.notna() & future30.notna()
        if int(valid_mask.sum()) >= 20:
            rising = future30[valid_mask] > 1.5
            if rising.any() and (~rising).any():
                pos = gp[valid_mask][rising]
                neg = gp[valid_mask][~rising]
                new_open = float(np.clip((pos.median() + neg.median()) / 2.0, 45.0, 70.0))
                profile["open_thr"] = 0.8 * float(profile.get("open_thr", 55.0)) + 0.2 * new_open
                profile["strong_open_thr"] = max(profile["open_thr"] + 15.0, float(profile.get("strong_open_thr", 75.0)))

    profile["source"] = "auto-cal"
    profile["updated_jst"] = now_jst().strftime("%Y-%m-%d %H:%M:%S JST")
    profile["rows_used"] = int(len(df))
    profile["valid_rows_used"] = valid_data_count
    profile["valid_segments"] = int(seg_count)
    profile["equal_ratio"] = float(equal_ratio)

    _learnable_weights = {
        "dfoF2":            profile.get("w_dfof2", 0.0),
        "mstid_approach":   profile.get("w_mstid_approach", 0.0),
        "fxEs_trend":       profile.get("w_trend", 0.12),
        "cluster":          profile.get("w_cluster", 0.75),
        "foF2_es_ratio":    profile.get("w_fof2", 0.0),
        "hmF2":             profile.get("w_hmf2", 0.0),
        "dhmF2":            profile.get("w_dhmf2", 0.0),
        "m3000f2":          profile.get("w_m3000f2", 0.0),
        "ns_wave":          profile.get("w_ns", 0.25),
        "shear":            profile.get("w_shear", 0.45),
        "gw_damp":          profile.get("w_gw_damp", 0.20),
        "rpi":              profile.get("w_rpi", 0.0),
        "es_mode":          profile.get("w_es_mode", 0.0),
    }
    total_w = sum(abs(v) for v in _learnable_weights.values() if not pd.isna(v))
    if total_w > 0:
        ranking = sorted(
            [(k, abs(float(v))) for k, v in _learnable_weights.items() if not pd.isna(v)],
            key=lambda x: x[1], reverse=True
        )
        ranking_pct = [(k, round(w / total_w * 100, 1)) for k, w in ranking]
    else:
        ranking_pct = []

    profile["feature_ranking"] = ranking_pct
    profile["feature_ranking_date"] = now_jst().strftime("%Y-%m-%d")

    EFF_RANK_LOOKBACK_ROWS = 96
    eff_hist_csv = log_dir / "effective_contribution_history.csv"
    effective_ranking_pct: List[Tuple[str, float]] = []
    if eff_hist_csv.exists():
        try:
            eff_hist_df = pd.read_csv(eff_hist_csv)
            eff_cols = [c for c in eff_hist_df.columns if c.startswith("eff_")]
            if eff_cols:
                recent_eff = eff_hist_df[eff_cols].tail(EFF_RANK_LOOKBACK_ROWS).apply(pd.to_numeric, errors="coerce")
                means = recent_eff.mean(numeric_only=True)
                total_m = float(means.abs().sum())
                if total_m > 0:
                    effective_ranking_pct = sorted(
                        [(c[len("eff_"):], round(float(v) / total_m * 100, 1))
                         for c, v in means.items() if not pd.isna(v)],
                        key=lambda x: x[1], reverse=True
                    )
        except Exception:
            effective_ranking_pct = []

    profile["effective_feature_ranking"] = effective_ranking_pct

    ranking_csv = log_dir / "feature_ranking_history.csv"
    try:
        rank_row = {"date": profile["feature_ranking_date"],
                    "equal_ratio": round(float(equal_ratio), 3),
                    "valid_rows": valid_data_count,
                    "ridge_used": profile.get("ridge_used", False)}
        for k, pct in ranking_pct:
            rank_row[f"pct_{k}"] = pct
        for k, pct in effective_ranking_pct:
            rank_row[f"effrank_{k}"] = pct
        rank_df_new = pd.DataFrame([rank_row])
        if ranking_csv.exists():
            rank_df_old = pd.read_csv(ranking_csv)
            pd.concat([rank_df_old, rank_df_new], ignore_index=True).to_csv(ranking_csv, index=False)
        else:
            rank_df_new.to_csv(ranking_csv, index=False)
    except Exception:
        pass

    profile_path.write_text(json.dumps(profile, ensure_ascii=False, indent=2), encoding='utf-8')

    rank_lines = []
    for i, (k, pct) in enumerate(ranking_pct[:5], 1):
        rank_lines.append(f"  {i}位 {k:<18} {pct:.1f}%")
    eff_rank_lines = []
    for i, (k, pct) in enumerate(effective_ranking_pct[:5], 1):
        eff_rank_lines.append(f"  {i}位 {k:<18} {pct:.1f}%")

    msg_lines = quality_lines + [
        f"自動補正: 実施（有効セグメント={seg_count} / MSTID重み={fmt(profile.get('w_mstid_approach',0),3)}"
        f" / open閾={fmt(profile.get('open_thr',55),1)}）",
    ]
    if rank_lines:
        msg_lines += ["特徴量寄与・固定重み込み(上位5):"] + rank_lines
    if eff_rank_lines:
        msg_lines += ["実効寄与ランキング(上位5):"] + eff_rank_lines
    return profile, "\n".join(msg_lines)


class SessionLogger:
    def __init__(self, enabled: bool, calibration_mode: bool,
                 log_interval: int, log_dir: Path, use_tmpfs: bool):
        self.enabled = bool(enabled)
        self.calibration_mode = bool(calibration_mode)
        self.log_interval = max(1, int(log_interval))
        self.log_dir = Path('/tmp') if use_tmpfs else Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.buffer: List[Dict] = []
        self.flush_counter = 0
        self.session_id = now_jst().strftime('%Y%m%d_%H%M%S')
        prefix = 'prediction_calibration' if calibration_mode else 'prediction'
        self.file_path = self.log_dir / f'{prefix}_{self.session_id}.csv'
        self.corr_csv = self.log_dir / f'{prefix}_{self.session_id}_corr.csv'
        self.summary_txt = self.log_dir / f'{prefix}_{self.session_id}_summary.txt'

    def add_rows(self, rows: List[Dict]):
        if not self.enabled or not rows:
            return
        self.buffer.extend(rows)
        self.flush_counter += 1
        if self.flush_counter >= self.log_interval:
            self.flush()

    def flush(self):
        """Ver10.8J: 終了時一括保存から「サイクル毎・上書き保存」に変更。
        Ver11.0J(改): mode='w'直接書き込みを廃止し、一時ファイル(.tmp)へ先に書いて
        os.replace()で原子的に置き換えるアトミック書き込みに変更。
        Pydroid3のOSキル・ラズパイの電源断・OOM killerでプロセスが書き込み途中に
        強制終了された場合、.tmpファイルのみが不完全になり、本番CSV(.file_path)は
        直前のflush時点の完全な状態で残る。mode='w'直接書き込みでは書き込み中断が
        即座にCSV破損となり、セッション全ログを失うリスクがあった。
        os.replace()はPOSIX/Windows両対応の原子的ファイル置換（Python 3.3+）。
        """
        if not self.enabled or not self.buffer:
            return
        tmp_path = self.file_path.with_suffix('.tmp')
        df = pd.DataFrame(self.buffer)
        try:
            df.to_csv(tmp_path, mode='w', header=True, index=False)
            os.replace(tmp_path, self.file_path)
        except Exception:
            try:
                tmp_path.unlink(missing_ok=True)
            except Exception:
                pass
        self.flush_counter = 0

    def close(self):
        if self.enabled:
            self.flush()
            self.generate_correlation_report()

    def generate_correlation_report(self):
        if not self.file_path.exists():
            return
        try:
            df = pd.read_csv(self.file_path)
        except Exception:
            return
        if df.empty or 'fxEs' not in df.columns:
            return
        rows = []
        lines = [f'session_id={self.session_id}', f'rows={len(df)}']
        for st in sorted(df['station'].dropna().unique()):
            ds = df[df['station']==st].copy().reset_index(drop=True)
            for lag in [1,2,3,4,6]:
                if len(ds) <= lag:
                    continue
                future_delta = pd.to_numeric(ds['fxEs'], errors='coerce').shift(-lag) - pd.to_numeric(ds['fxEs'], errors='coerce')
                for col in ['mstid_approach_index','mstid_strength','gw_modulation_index','trend','cluster',
                            'foF2_ratio','foF2_es_ratio','dfoF2','hmF2','dhmF2','m3000f2_norm',
                            'rpi_norm','es_mode_score','es_hop_ratio','distance_std','bearing_std','dist_in_window_count',
                            'bearing_axial_r','bearing_orientation_deg','bearing_ns_alignment']:
                    if col not in ds.columns:
                        continue
                    x = pd.to_numeric(ds[col], errors='coerce')
                    y = pd.to_numeric(future_delta, errors='coerce')
                    valid = x.notna() & y.notna()
                    corr = float(np.corrcoef(x[valid], y[valid])[0,1]) if valid.sum() >= 3 else np.nan
                    rows.append({'station': st, 'lag_min': lag*15, 'feature': col,
                                 'corr_with_future_fxEs_delta': corr, 'n': int(valid.sum())})
        if rows:
            out = pd.DataFrame(rows)
            out.to_csv(self.corr_csv, index=False)
            tmp = out.dropna(subset=['corr_with_future_fxEs_delta'])
            if not tmp.empty:
                best = tmp.iloc[tmp['corr_with_future_fxEs_delta'].abs().argmax()]
                lines.append(f"best_corr feature={best['feature']} station={best['station']} lag_min={int(best['lag_min'])} corr={best['corr_with_future_fxEs_delta']:.3f} n={int(best['n'])}")

        pcorr_rows = []
        pcorr_feat_cols = [c for c in [
            'fxEs', 'trend', 'cluster', 'foF2_es_ratio', 'dfoF2',
            'hmF2', 'dhmF2', 'm3000f2_norm', 'mstid_approach_index',
            'rpi_norm', 'es_mode_score', 'distance_std', 'bearing_std',
            'bearing_axial_r', 'bearing_orientation_deg', 'bearing_ns_alignment',
        ] if c in df.columns]

        for st in sorted(df['station'].dropna().unique()):
            ds = df[df['station'] == st].copy().reset_index(drop=True)
            if len(ds) < 10:
                continue
            fx = pd.to_numeric(ds['fxEs'], errors='coerce')
            future30 = fx.shift(-2) - fx

            avail = [c for c in pcorr_feat_cols if c in ds.columns]
            if len(avail) < 2:
                continue

            sub = ds[avail].apply(pd.to_numeric, errors='coerce').copy()
            sub['__y__'] = future30

            sub = sub.dropna()
            if len(sub) < 10:
                continue

            for target in avail:
                other = [c for c in avail if c != target]
                if not other:
                    continue
                try:
                    X_oth = np.column_stack([np.ones(len(sub))] + [sub[c].values for c in other])
                    coef_x, *_ = np.linalg.lstsq(X_oth, sub[target].values, rcond=None)
                    resid_x = sub[target].values - X_oth @ coef_x
                    coef_y, *_ = np.linalg.lstsq(X_oth, sub['__y__'].values, rcond=None)
                    resid_y = sub['__y__'].values - X_oth @ coef_y
                    if resid_x.std() < 1e-12 or resid_y.std() < 1e-12:
                        pcorr = np.nan
                    else:
                        pcorr = float(np.corrcoef(resid_x, resid_y)[0, 1])
                except Exception:
                    pcorr = np.nan
                pcorr_rows.append({
                    'station': st, 'feature': target,
                    'partial_corr_future30': pcorr, 'n': int(len(sub)),
                })

        if pcorr_rows:
            pcorr_df = pd.DataFrame(pcorr_rows)
            pcorr_path = self.log_dir / f"partial_corr_{self.session_id}.csv"
            pcorr_df.to_csv(pcorr_path, index=False)
            tmp_pc = pcorr_df.dropna(subset=['partial_corr_future30'])
            if not tmp_pc.empty:
                top3 = tmp_pc.reindex(
                    tmp_pc['partial_corr_future30'].abs().nlargest(3).index
                )
                lines.append("偏相関TOP3(分析専用):")
                for _, r in top3.iterrows():
                    lines.append(
                        f"  {r['feature']}@{r['station']}: pcorr={r['partial_corr_future30']:.3f} n={int(r['n'])}"
                    )

        self.summary_txt.write_text('\n'.join(lines), encoding='utf-8')



def add_dynamics(history: pd.DataFrame) -> pd.DataFrame:
    df = history.copy()
    if df.empty:
        return df
    df['obs_dt'] = pd.to_datetime(df['obs_jst'], errors='coerce')
    df = df.sort_values(['station','obs_dt'])
    out = []
    for _, g in df.groupby('station'):
        g = g.copy().sort_values('obs_dt')
        dt_h = g['obs_dt'].diff().dt.total_seconds()/3600.0
        g['d_inst'] = (g['fxEs'].diff()/dt_h).clip(-MAX_ABS_TREND_MHZ_PER_HOUR, MAX_ABS_TREND_MHZ_PER_HOUR)
        def slope(y):
            yy = pd.Series(y).dropna().values
            if len(yy) < 3:
                return np.nan
            return np.polyfit(np.arange(len(yy))*0.25, yy, 1)[0]
        g['d_reg'] = g['fxEs'].rolling(5, min_periods=3).apply(slope, raw=False)
        g['d_reg'] = g['d_reg'].clip(-MAX_ABS_TREND_MHZ_PER_HOUR, MAX_ABS_TREND_MHZ_PER_HOUR)
        g['d_eff'] = 0.7*g['d_reg'] + 0.3*g['d_inst']
        conflict = np.sign(g['d_inst']) != np.sign(g['d_reg'])
        g.loc[conflict, 'd_eff'] *= 0.5
        g['fxes_std_60'] = g.rolling('60min', on='obs_dt', min_periods=3)['fxEs'].std()
        g['d2_eff'] = g['d_eff'].diff()
        out.append(g)
    return pd.concat(out, ignore_index=True) if out else df


def latest_network_sample(dyn: pd.DataFrame) -> pd.DataFrame:
    if dyn.empty:
        return dyn
    latest_time = dyn['obs_dt'].max()
    s = dyn[(dyn['obs_dt'] >= latest_time - pd.Timedelta(minutes=8)) & (dyn['obs_dt'] <= latest_time + pd.Timedelta(minutes=8))]
    return s.sort_values('obs_dt').drop_duplicates('station', keep='last')


def cluster_score(sample: pd.DataFrame, station: str):
    if sample.empty:
        return np.nan
    active = sample.assign(active=lambda x: (x['fxEs'] >= 7.0).astype(float))
    vals = active.set_index('station')['active'].to_dict()
    if station == 'wakkanai':
        peers = ['wakkanai','kokubunji']
    elif station == 'okinawa':
        peers = ['okinawa','yamagawa']
    elif station == 'kokubunji':
        peers = ['wakkanai','kokubunji','yamagawa']
    else:
        peers = ['kokubunji','yamagawa','okinawa']
    present = [vals[p] for p in peers if p in vals]
    return float(np.mean(present)) if present else np.nan


def cluster_phrase(score, station: str) -> str:
    if pd.isna(score):
        return '不明'
    total = 2 if station in ('wakkanai','okinawa') else 3
    active = int(round(float(score)*total))
    if active <= 0:
        label='孤立的'
    elif active == total:
        label='広域的'
    elif total == 2:
        label='近傍的'
    elif active == 1:
        label='局所的'
    else:
        label='拡大的'
    return f'{label}({total}局中{active}局)'


def mstid_sign_change_score(sample: pd.DataFrame):
    if sample.empty or len(sample) < 3:
        return np.nan, 0
    s = sample.copy(); s['ord'] = s['station'].map(lambda k: STATIONS[k]['lat_order']); s = s.sort_values('ord')
    grads = np.diff(s['fxEs'].values)
    signs = [np.sign(g) for g in grads if abs(g) >= 0.3]
    if len(signs) < 2:
        return 0.0, 0
    changes = sum(1 for a,b in zip(signs[:-1], signs[1:]) if a != b)
    return clamp(changes/2.0,0.0,1.0), changes


def north_south_phrase(changes) -> str:
    if changes is None or pd.isna(changes):
        return '不明'
    if changes <= 0:
        return '傾斜型'
    if changes == 1:
        return '波状傾向'
    return '波状強め'



LOCAL_SHAPE_ORDER = ["Cluster", "Shear", "Wave", "Diffuse"]
LOCAL_SHAPE_JP = {"Cluster": "Cluster", "Shear": "Shear", "Wave": "Wave", "Diffuse": "Diffuse"}


def compute_local_shape_scores(
    cscore: float,
    shear_feats: Dict[str, float],
    gw_feats: Dict[str, float],
    fxes_std60: float,
    mstid_wave_clear: bool,
) -> Dict[str, float]:
    """局ごとに独立したEs形状スコア(Cluster/Shear/Wave/Diffuse)を算出する。
    仕様書2章の局別特徴量のうち、既存パイプラインで既に自局値として
    得られているものを流用する（cscore=局所cluster_score、
    shear_feats/gw_feats/fxes_std60は全て当該局gから計算済み）。

    Wave成分は仕様書5章のMSTID整合性チェックに従い、
    MSTID波面が明瞭な場合のみ積極採用し、不明瞭・クラスター優勢の
    場合は減衰係数を掛ける。
    """
    c = 0.0 if pd.isna(cscore) else clamp(float(cscore), 0.0, 1.0)

    conv = shear_feats.get('convergence_index', np.nan)
    s = 0.0 if pd.isna(conv) else clamp(float(conv) / 1.5, 0.0, 1.0)

    gwm = gw_feats.get('gw_modulation_index', np.nan)
    w_raw = 0.0 if pd.isna(gwm) else clamp(abs(float(gwm)), 0.0, 1.0)
    w = clamp(w_raw * 1.2, 0.0, 1.0) if mstid_wave_clear else w_raw * 0.4

    d = 0.0 if pd.isna(fxes_std60) else clamp(float(fxes_std60) / 3.0, 0.0, 1.0)

    return {"Cluster": c, "Shear": s, "Wave": w, "Diffuse": d}


def local_shape_label(scores: Dict[str, float]) -> str:
    if not scores or all(pd.isna(v) or v <= 0.0 for v in scores.values()):
        return "不明"
    return max(scores, key=lambda k: scores[k])


def spatial_continuity_penalty(shapes_north_to_south: List[str]) -> float:
    """北→南に並べた形状列の空間連続性を評価する（仕様書6章）。
    Cluster→Shear→Wave→Diffuseのような滑らかな順序変化ならペナルティ0に近く、
    Cluster→Wave→Cluster→Waveのような交互反転が多いほどペナルティが大きくなる。
    形状ラベル自体は変更せず、戻り値はConfidence補正にのみ用いる。
    """
    idx = {name: i for i, name in enumerate(LOCAL_SHAPE_ORDER)}
    steps = [idx[s] for s in shapes_north_to_south if s in idx]
    if len(steps) < 3:
        return 0.0
    diffs = np.diff(steps)
    nonzero = [d for d in diffs if d != 0]
    if len(nonzero) < 2:
        return 0.0
    reversals = sum(1 for a, b in zip(nonzero[:-1], nonzero[1:]) if np.sign(a) != np.sign(b))
    return clamp(reversals / (len(nonzero) - 1), 0.0, 1.0)


def local_shape_confidence(
    scores: Dict[str, float],
    analysis_rate: float,
    eqi: float,
    mstid_wave_clear: bool,
    continuity_penalty: float,
) -> float:
    """局別Es形状のConfidence(0-100%)を算出する（仕様書7章）。
    要素: イオノグラム品質(analysis_rate代替)・EQI・MSTID整合性・空間連続性。
    """
    if not scores:
        return 0.0
    vals = sorted(scores.values(), reverse=True)
    top = vals[0]
    second = vals[1] if len(vals) > 1 else 0.0
    margin = clamp(top - second, 0.0, 1.0)
    q = 0.0 if pd.isna(analysis_rate) else clamp(float(analysis_rate) / 100.0, 0.0, 1.0)
    e = 0.0 if pd.isna(eqi) else clamp(float(eqi) / 100.0, 0.0, 1.0)
    m = 1.0 if mstid_wave_clear else 0.7
    base = 0.40 * margin + 0.25 * q + 0.20 * e + 0.15 * m
    base *= (1.0 - 0.5 * continuity_penalty)
    return round(clamp(base, 0.0, 1.0) * 100.0, 0)


def compute_vertical_structure_features(g: pd.DataFrame) -> Dict[str,float]:
    return {'layer_height': np.nan, 'layer_thickness': np.nan, 'blanketing_index': np.nan, 'height_descent_rate': np.nan}


def load_wind_profile_nearest(obs_dt, station: str, wind_file: Path) -> pd.DataFrame:
    wind_path = Path(wind_file)
    if not wind_path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(wind_path)
    except Exception:
        return pd.DataFrame()
    req = {'datetime_jst','station','z_km','u_ms','v_ms'}
    if not req.issubset(set(df.columns)):
        return pd.DataFrame()
    df = df[df['station']==station].copy()
    if df.empty:
        return pd.DataFrame()
    df['obs_dt'] = pd.to_datetime(df['datetime_jst'], errors='coerce')
    if df['obs_dt'].isna().all():
        return pd.DataFrame()
    nearest_time = min(df['obs_dt'].dropna().unique(), key=lambda x: abs(pd.Timestamp(x)-pd.Timestamp(obs_dt)))
    prof = df[df['obs_dt']==nearest_time].copy().sort_values('z_km')
    return prof[['z_km','u_ms','v_ms']]


def compute_shear_features(profile: pd.DataFrame, vstruct: Dict[str,float], g: pd.DataFrame, force_no_wind: bool=False) -> Dict[str,float]:
    feats = {'shear_peak_strength': np.nan, 'shear_peak_height': np.nan, 'convergence_index': np.nan, 'shear_mode': 'none'}
    if (not force_no_wind) and profile is not None and not profile.empty and len(profile) >= 3:
        z = profile['z_km'].astype(float).values; u = profile['u_ms'].astype(float).values
        try:
            du_dz = np.gradient(u, z)
            neg = -du_dz
            idx = int(np.nanargmax(neg))
            feats['shear_peak_strength'] = float(neg[idx])
            feats['shear_peak_height'] = float(z[idx])
            feats['convergence_index'] = float(max(0.0, neg[idx]))
            feats['shear_mode'] = 'profile'
            return feats
        except Exception:
            pass
    trend = np.nan if g.empty else g.sort_values('obs_dt').iloc[-1].get('d_eff', np.nan)
    rising = 0.0 if pd.isna(trend) else max(0.0, float(trend))
    pseudo = 0.40 * rising
    feats['convergence_index'] = pseudo
    feats['shear_peak_strength'] = pseudo
    feats['shear_mode'] = 'pseudo'
    return feats


def compute_lowfreq_baseline(g: pd.DataFrame):
    if g.empty:
        return np.nan
    latest = g.sort_values('obs_dt').iloc[-1]
    trend = latest.get('d_eff', np.nan)
    trend = 0.0 if pd.isna(trend) else float(trend)
    return float(latest['fxEs']) + trend*0.25


def compute_gw_residual_features(g: pd.DataFrame) -> Dict[str,float]:
    feats = {'gw_residual': np.nan, 'gw_modulation_index': 0.0, 'forecast_damping': 0.0}
    if g.empty:
        return feats
    gg = g.sort_values('obs_dt').copy()
    if len(gg) < 2:
        return feats
    baseline_prev = compute_lowfreq_baseline(gg.iloc[:-1]) if len(gg) >= 2 else np.nan
    if not pd.isna(baseline_prev):
        feats['gw_residual'] = float(gg.iloc[-1]['fxEs'] - baseline_prev)
    if len(gg) >= 3:
        x0 = float(gg.iloc[-3]['fxEs']); x1 = float(gg.iloc[-2]['fxEs']); x2 = float(gg.iloc[-1]['fxEs'])
        feats['gw_modulation_index'] = abs(x2 - 2*x1 + x0)
    else:
        feats['gw_modulation_index'] = abs(float(gg.iloc[-1]['fxEs'] - gg.iloc[-2]['fxEs']))
    feats['forecast_damping'] = 0.5 * feats['gw_modulation_index']
    return feats




def fetch_nict_foF2(url: Optional[str] = None, timeout: int = HTTP_TIMEOUT_SEC,
                     debug: bool = False) -> pd.DataFrame:
    """ITPFS(整理版): F2層(電離層F2層)伝搬予測機能の廃止に伴い無効化。
    NICT側の該当ページ群(斜め伝搬可能周波数/電離層概況詳細版/旧foF2ページ)から
    安定してfoF2/hmF2/M3000F2を取得することが実運用上不可能と判断したため、
    本関数は常に空DataFrameを返す(ネットワーク通信自体を行わない)。
    唯一の呼び出し元compute_f2_layer_forecast()は空リストを返す縮退実装のため、
    本関数を無効化するだけで安全に縮退動作する。
    ESDUCT(整理版・不具合修正): 以前は本関数のfoF2取得失敗時に、HFバンド
    モニタ側でfoF2=6.5MHz/M3000F2=2.9固定のフォールバック代表値を使い続ける
    残骸コードが存在していたが、これは削除した(F2予測機能の復活ではなく、
    常に取得不能な値を前提とした紛らわしい表示だったため)。現在のHFバンド
    モニタ(print_hf_band_bargraph_monitor_silent()/_simple())はF2データに
    一切依存せず、実測fxEsのみに基づくEsバンド(28.5/50MHz)判定のみを表示する。
    columns(常に空): station, obs_jst, foF2, hmF2, M3000F2
    """
    return pd.DataFrame()




_sw_history: List[Dict] = []
_eqi_history: Dict[str, List[Dict]] = {}
_es_phase_history: Dict[str, List[Dict]] = {}

_latest_sw_features: Dict = {}
_latest_eqi_by_station: Dict[str, Dict] = {}
_latest_es_phase_by_station: Dict[str, str] = {}
_latest_es_phase_num_by_station: Dict[str, float] = {}

ES_PHASE_LABELS = ["FORMATION", "GROWTH", "PEAK", "DECAY", "COLLAPSE"]
_ES_PHASE_NUM = {p: float(i) for i, p in enumerate(ES_PHASE_LABELS)}



def purge_old_cache(cache_dir: Path = SW_CACHE_DIR,
                    ttl_hours: float = SW_CACHE_TTL_HOURS) -> int:
    """cacheディレクトリ内の古いファイルを削除する。
    戻り値: 削除したファイル数。
    """
    if not cache_dir.exists():
        return 0
    threshold = time.time() - ttl_hours * 3600.0
    removed = 0
    for f in cache_dir.iterdir():
        if f.is_file() and f.stat().st_mtime < threshold:
            try:
                f.unlink()
                removed += 1
            except Exception:
                pass
    return removed


def cache_path_for(url: str, ext: str = ".bin") -> Path:
    """URLから一意なキャッシュパスを生成する。"""
    import hashlib
    h = hashlib.md5(url.encode()).hexdigest()[:12]
    SW_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return SW_CACHE_DIR / f"{h}{ext}"



def _fetch_json_noaa(url: str, timeout: int = SW_FETCH_TIMEOUT) -> Optional[List]:
    """NOAA SWPC JSONエンドポイントから配列データを取得する。
    失敗時はNoneを返す（呼び出し側でグレースフルデグレード）。
    """
    try:
        r = requests.get(url, timeout=timeout,
                         headers={"User-Agent": "CBM2.py EDFS Ver11.4; JL7KHN"})
        r.raise_for_status()
        data = r.json()
        if not isinstance(data, list) or len(data) < 2:
            return None
        return data
    except Exception:
        return None


def parse_noaa_mag(data: Optional[List]) -> pd.DataFrame:
    """IMF磁場データ(mag-1-day.json)をDataFrameに変換する。
    columns: time_tag, bz, bt, bx, by
    """
    if not data:
        return pd.DataFrame()
    try:
        header = data[0]
        rows = data[1:]
        df = pd.DataFrame(rows, columns=header)
        df["time_tag"] = pd.to_datetime(df["time_tag"], errors="coerce", utc=True)
        for col in ["bz_gsm", "bt", "bx_gsm", "by_gsm"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df.dropna(subset=["time_tag"]).sort_values("time_tag").tail(96)
    except Exception:
        return pd.DataFrame()


def parse_noaa_plasma(data: Optional[List]) -> pd.DataFrame:
    """太陽風プラズマデータをDataFrameに変換する。
    columns: time_tag, speed, density, temperature
    """
    if not data:
        return pd.DataFrame()
    try:
        header = data[0]
        rows = data[1:]
        df = pd.DataFrame(rows, columns=header)
        df["time_tag"] = pd.to_datetime(df["time_tag"], errors="coerce", utc=True)
        for col in ["speed", "density", "temperature"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df.dropna(subset=["time_tag"]).sort_values("time_tag").tail(96)
    except Exception:
        return pd.DataFrame()


def parse_noaa_kp(data: Optional[List]) -> pd.DataFrame:
    """Kp指数データをDataFrameに変換する。"""
    if not data:
        return pd.DataFrame()
    try:
        header = data[0]
        rows = data[1:]
        df = pd.DataFrame(rows, columns=header)
        df["time_tag"] = pd.to_datetime(df["time_tag"], errors="coerce", utc=True)
        kp_col = next((c for c in df.columns if "kp" in c.lower()), None)
        if kp_col:
            df["kp"] = pd.to_numeric(df[kp_col], errors="coerce")
        return df.dropna(subset=["time_tag"]).sort_values("time_tag").tail(96)
    except Exception:
        return pd.DataFrame()


def compute_sw_time_features(df: pd.DataFrame, col: str,
                              windows_min: List[int] = [15, 30, 60]) -> Dict[str, float]:
    """時系列DataFrameの指定列から時間変化特徴量を生成する。
    生成: delta_15m, delta_30m, mean_60m, slope_60m, std_60m
    仕様書: 「取得値そのものではなく時間変化を特徴量として学習器へ入力する」
    """
    out: Dict[str, float] = {}
    if df.empty or col not in df.columns:
        return out
    df2 = df.copy()
    if "time_tag" not in df2.columns:
        return out
    df2 = df2.dropna(subset=["time_tag", col]).sort_values("time_tag")
    series = df2.set_index("time_tag")[col].apply(pd.to_numeric, errors="coerce").dropna()
    if series.empty:
        return out
    now_utc = pd.Timestamp.utcnow()
    for wmin in windows_min:
        cut = now_utc - pd.Timedelta(minutes=wmin)
        seg = series[series.index >= cut]
        if len(seg) >= 2:
            delta = float(seg.iloc[-1] - seg.iloc[0])
            out[f"{col}_delta_{wmin}m"] = delta
        elif len(seg) == 1:
            out[f"{col}_delta_{wmin}m"] = 0.0
    cut60 = now_utc - pd.Timedelta(minutes=60)
    seg60 = series[series.index >= cut60]
    if len(seg60) >= 3:
        out[f"{col}_mean_60m"] = float(seg60.mean())
        out[f"{col}_std_60m"]  = float(seg60.std())
        t_idx = np.arange(len(seg60), dtype=float)
        try:
            slope = float(np.polyfit(t_idx, seg60.values, 1)[0])
        except Exception:
            slope = 0.0
        out[f"{col}_slope_60m"] = slope
    elif len(seg60) >= 1:
        out[f"{col}_mean_60m"] = float(seg60.mean())
        out[f"{col}_std_60m"]  = 0.0
        out[f"{col}_slope_60m"] = 0.0
    return out


def fetch_space_weather_features() -> Dict[str, float]:
    """NOAA SWPCから宇宙天気特徴量を取得・生成する。
    仕様書: 取得後に時間変化特徴量へ変換する（単純保存ではない）。
    グレースフルデグレード: 各エンドポイント失敗時は空辞書とマージ。
    """
    feats: Dict[str, float] = {}

    mag_data = _fetch_json_noaa(NOAA_SWPC_MAG_URL)
    mag_df = parse_noaa_mag(mag_data)
    if not mag_df.empty:
        bz_col = "bz_gsm" if "bz_gsm" in mag_df.columns else None
        if bz_col:
            feats.update(compute_sw_time_features(mag_df, bz_col))
            bz_series = mag_df.set_index("time_tag")[bz_col].dropna().sort_index()
            if not bz_series.empty:
                neg_streak = 0
                for v in reversed(bz_series.values):
                    if v < 0:
                        neg_streak += 1
                    else:
                        break
                feats["bz_southward_streak_min"] = float(neg_streak) * 1.0
            bt_col = "bt" if "bt" in mag_df.columns else None
            if bt_col:
                feats.update(compute_sw_time_features(mag_df, bt_col))

    plasma_data = _fetch_json_noaa(NOAA_SWPC_PLASMA_URL)
    plasma_df = parse_noaa_plasma(plasma_data)
    if not plasma_df.empty:
        for col in ["speed", "density"]:
            if col in plasma_df.columns:
                feats.update(compute_sw_time_features(plasma_df, col))
        if "speed" in plasma_df.columns and "density" in plasma_df.columns:
            plasma_df["dyn_pressure"] = (
                pd.to_numeric(plasma_df["speed"], errors="coerce") ** 2
                * pd.to_numeric(plasma_df["density"], errors="coerce")
            )
            feats.update(compute_sw_time_features(plasma_df, "dyn_pressure"))

    kp_data = _fetch_json_noaa(NOAA_SWPC_KP_URL)
    kp_df = parse_noaa_kp(kp_data)
    if not kp_df.empty and "kp" in kp_df.columns:
        feats.update(compute_sw_time_features(kp_df, "kp"))
        feats["kp_latest"] = float(kp_df["kp"].dropna().iloc[-1]) if kp_df["kp"].notna().any() else np.nan
        feats["kp_delta_60m_rate"] = feats.get("kp_delta_60m", np.nan)

    el_data = _fetch_json_noaa(NOAA_SWPC_ELECTRON_URL)
    if el_data:
        try:
            el_header = el_data[0]; el_rows = el_data[1:]
            el_df = pd.DataFrame(el_rows, columns=el_header)
            el_df["time_tag"] = pd.to_datetime(el_df["time_tag"], errors="coerce", utc=True)
            flux_col = next((c for c in el_df.columns if "flux" in c.lower()), None)
            if flux_col:
                el_df[flux_col] = pd.to_numeric(el_df[flux_col], errors="coerce")
                feats.update(compute_sw_time_features(el_df, flux_col))
                feats["electron_flux_latest"] = float(
                    el_df[flux_col].dropna().iloc[-1]) if el_df[flux_col].notna().any() else np.nan
        except Exception:
            pass

    return feats


def compute_sw_score(feats: Dict[str, float]) -> float:
    """Space Weather特徴量から総合スコア(0-1)を生成する。
    低スコア=宇宙天気穏やか(Es有利)、高スコア=擾乱大(Es不安定化リスク)。
    仕様書因果モデル: Space Weather → EQI悪化 → FT8減少 の入力となる。
    """
    score = 0.0
    wt = 0.0

    def _add(val, w, direction=1.0, scale=1.0):
        nonlocal score, wt
        if not pd.isna(val):
            score += direction * clamp(float(val) / scale, -1.0, 1.0) * w
            wt += w

    _add(feats.get("bz_gsm_delta_60m", np.nan), 0.30, direction=-1.0, scale=10.0)
    _add(feats.get("bz_southward_streak_min", np.nan), 0.10, direction=1.0, scale=60.0)
    _add(feats.get("kp_latest", np.nan), 0.25, direction=1.0, scale=5.0)
    _add(feats.get("kp_delta_60m", np.nan), 0.10, direction=1.0, scale=3.0)
    _add(feats.get("speed_delta_30m", np.nan), 0.10, direction=1.0, scale=100.0)
    _add(feats.get("electron_flux_latest", np.nan), 0.05, direction=1.0, scale=1e5)
    _add(feats.get("dyn_pressure_delta_30m", np.nan), 0.10, direction=1.0, scale=5e9)

    if wt <= 0:
        return np.nan
    raw = score / wt
    return float(clamp((raw + 1.0) / 2.0, 0.0, 1.0))


def update_sw_history(feats: Dict[str, float], sw_score: float) -> None:
    """Space Weather特徴量をRAM履歴に追加する。"""
    _sw_history.append({
        "obs_jst": now_jst().isoformat(),
        "sw_score": sw_score,
        **feats,
    })
    if len(_sw_history) > SW_HISTORY_MAXLEN:
        _sw_history[:] = _sw_history[-SW_HISTORY_MAXLEN:]



_KP_TO_AP_TABLE = {
    0.000: 0,   0.333: 2,   0.667: 3,   1.000: 4,   1.333: 5,
    1.667: 6,   2.000: 7,   2.333: 9,   2.667: 12,  3.000: 15,
    3.333: 18,  3.667: 22,  4.000: 27,  4.333: 32,  4.667: 39,
    5.000: 48,  5.333: 56,  5.667: 67,  6.000: 80,  6.333: 94,
    6.667: 111, 7.000: 132, 7.333: 154, 7.667: 179, 8.000: 207,
    8.333: 236, 8.667: 300, 9.000: 400,
}


def kp_to_ap(kp: float) -> float:
    """3時間Kp値を ap(等価地磁気活動度)へ変換する。日次Apが未発表の時間帯でも
    直近Kpからリアルタイムに擾乱度合いを推定するための代理値。"""
    if kp is None or pd.isna(kp):
        return np.nan
    keys = sorted(_KP_TO_AP_TABLE.keys())
    nearest = min(keys, key=lambda k: abs(k - float(kp)))
    return float(_KP_TO_AP_TABLE[nearest])


def xray_flux_to_class(flux_wm2: float) -> str:
    """GOES X線強度[W/m^2]をA/B/C/M/Xクラス表記へ変換する。"""
    if flux_wm2 is None or pd.isna(flux_wm2) or flux_wm2 <= 0:
        return "--"
    if flux_wm2 >= 1e-4:
        return f"X{flux_wm2/1e-4:.1f}"
    if flux_wm2 >= 1e-5:
        return f"M{flux_wm2/1e-5:.1f}"
    if flux_wm2 >= 1e-6:
        return f"C{flux_wm2/1e-6:.1f}"
    if flux_wm2 >= 1e-7:
        return f"B{flux_wm2/1e-7:.1f}"
    return f"A{flux_wm2/1e-8:.1f}"


def fetch_realtime_solar_indices() -> Dict[str, float]:
    """SFI(F10.7実測)・Kp・Ap換算・X線強度をNOAA SWPCから1回だけリアルタイム取得する。
    Ver12.1までの15分自動マージ・Space Weather RAM履歴とは独立した単発取得。
    取得失敗時は該当キーがNaN/"--"のまま返る（グレースフルデグレード）。
    """
    out: Dict[str, float] = {"sfi": np.nan, "kp": np.nan, "ap": np.nan,
                              "xray_flux": np.nan, "xray_class": "--"}
    try:
        r = requests.get(NOAA_SWPC_F107_SUMMARY_URL, timeout=SW_FETCH_TIMEOUT,
                          headers={"User-Agent": "CBM2.py EDFS Ver12.1; JL7KHN"})
        r.raise_for_status()
        d = r.json()
        v = d.get("solar_radio_flux", d.get("flux", d.get("Flux")))
        out["sfi"] = safe_float(v)
    except Exception:
        pass
    if pd.isna(out["sfi"]):
        try:
            data = _fetch_json_noaa(NOAA_SWPC_F107_URL)
            if data:
                header = data[0]; rows = data[1:]
                df = pd.DataFrame(rows, columns=header)
                flux_col = next((c for c in df.columns if "flux" in c.lower()), None)
                if flux_col:
                    df[flux_col] = pd.to_numeric(df[flux_col], errors="coerce")
                    s = df[flux_col].dropna()
                    if not s.empty:
                        out["sfi"] = float(s.iloc[-1])
        except Exception:
            pass
    try:
        kp_data = _fetch_json_noaa(NOAA_SWPC_KP_URL)
        kp_df = parse_noaa_kp(kp_data)
        if not kp_df.empty and "kp" in kp_df.columns and kp_df["kp"].notna().any():
            out["kp"] = float(kp_df["kp"].dropna().iloc[-1])
            out["ap"] = kp_to_ap(out["kp"])
    except Exception:
        pass
    try:
        xr_data = _fetch_json_noaa(NOAA_SWPC_PROTON_URL)
        if xr_data:
            header = xr_data[0]; rows = xr_data[1:]
            xdf = pd.DataFrame(rows, columns=header)
            flux_col = next((c for c in xdf.columns if "flux" in c.lower()), None)
            if flux_col:
                xdf[flux_col] = pd.to_numeric(xdf[flux_col], errors="coerce")
                s = xdf[flux_col].dropna()
                if not s.empty:
                    fx = float(s.iloc[-1])
                    out["xray_flux"] = fx
                    out["xray_class"] = xray_flux_to_class(fx)
    except Exception:
        pass
    return out


def solar_zenith_deg(lat: float, lon: float, dt_utc: datetime) -> float:
    """太陽天頂角の実用近似計算（均時差は無視、時刻誤差±10分程度相当）。
    D層吸収モデルの昼夜判定・吸収量推定に使用する。高精度な天体暦計算は行わない。
    """
    try:
        n = dt_utc.timetuple().tm_yday
        decl = -23.44 * math.cos(math.radians(360.0 / 365.0 * (n + 10)))
        frac_hour = dt_utc.hour + dt_utc.minute / 60.0 + dt_utc.second / 3600.0
        solar_time = frac_hour + lon / 15.0
        hour_angle = 15.0 * (solar_time - 12.0)
        lat_r = math.radians(lat); decl_r = math.radians(decl); ha_r = math.radians(hour_angle)
        cos_z = (math.sin(lat_r) * math.sin(decl_r)
                 + math.cos(lat_r) * math.cos(decl_r) * math.cos(ha_r))
        cos_z = max(-1.0, min(1.0, cos_z))
        return math.degrees(math.acos(cos_z))
    except Exception:
        return np.nan


def nearest_nict_station(lat: float, lon: float) -> str:
    """指定緯度経度に最も近いNICT電離層観測局を返す。"""
    best, best_d = "kokubunji", 1e9
    for st, (slat, slon) in NICT_STATION_LATLON.items():
        d = haversine_km(lat, lon, slat, slon)
        if d < best_d:
            best_d = d; best = st
    return best


def sfi_to_ssn(sfi: float) -> float:
    """SFI→SSN概算変換（標準的経験式 SFI=63.75+0.53*SSN の逆算）。"""
    if pd.isna(sfi):
        return np.nan
    return max(0.0, (float(sfi) - 63.75) / 0.53)






def compute_f2_layer_forecast() -> List[Dict]:
    """ITPFS(整理版): F2層(電離層F2層)伝搬予測機能は廃止した。
    NICTからのfoF2/hmF2データ取得が実運用上不可能と判断したため、Es予測とは
    独立していたこの機能(遠隔方面へのMUF/LUF評価・開通期待度表示)を無効化する。
    fetch_nict_foF2()が常に空DataFrameを返すため以下の計算自体は成立しないが、
    呼び出し側(print_unified_display/ダッシュボード)との互換性のため、
    空リストを返す形はそのまま維持する(F2関連パネルは自動的に非表示になる)。
    """
    return []


def empirical_probability_snapshot(fx_for_prob: float, ef_feats: Dict[str, float],
                                    profile: Dict) -> float:
    """起動直後（セッション履歴0件）でも算出できる範囲で empirical_probability() の
    z式を再現するスナップショット版。トレンド/クラスタ/シア/GW減衰/MSTID等、
    複数サイクルの観測蓄積が前提の特徴量は使わず(寄与0固定)、fxEs実測値と、
    その端末の過去ログから auto_calibrate_profile() が自動較正した
    prob_intercept / w_fxes / F層学習重み(w_fof2・w_hmf2・w_m3000f2等)のみで
    確率を算出する。empirical_probability()側にあるサイクル数キャップ
    （暫定/上限50%等）は適用しない。トレンドを使わない以上、サイクル数に
    依存する確信度キャップの前提がそもそも当てはまらないため。
    """
    if pd.isna(fx_for_prob):
        return np.nan
    z = float(profile.get('prob_intercept', -2.0))
    z += float(profile.get('w_fxes', 0.32)) * (float(fx_for_prob) - 6.0)
    w_fof2 = float(profile.get('w_fof2', 0.0))
    if w_fof2 != 0.0:
        fes = ef_feats.get('foF2_es_ratio', np.nan)
        z_fes = std_z(fes, 'foF2_es_ratio', profile)
        if not pd.isna(z_fes):
            sign_fof2 = float(profile.get('sign_fof2', 0.0))
            z += w_fof2 * sign_fof2 * z_fes
    w_hmf2 = float(profile.get('w_hmf2', 0.0))
    if w_hmf2 != 0.0:
        z_hmf2 = std_z(ef_feats.get('hmF2', np.nan), 'hmF2', profile)
        if not pd.isna(z_hmf2):
            sign_hmf2 = float(profile.get('sign_hmf2', 0.0))
            z += w_hmf2 * sign_hmf2 * z_hmf2
    w_m3000f2 = float(profile.get('w_m3000f2', 0.0))
    if w_m3000f2 != 0.0:
        z_m3000f2 = std_z(ef_feats.get('m3000f2_norm', np.nan), 'm3000f2_norm', profile)
        if not pd.isna(z_m3000f2):
            sign_m3000f2 = float(profile.get('sign_m3000f2', 0.0))
            z += w_m3000f2 * sign_m3000f2 * z_m3000f2
    p = sigmoid(z) * 100.0
    return clamp(p, 1.0, PROB_MAX_PERCENT)


def state_from_prob_snapshot(p: float, profile: Dict) -> str:
    """state_from_prob() からサイクル数ゲート(暫定判定/確信度キャップ)を除いた版。
    閾値(strong_open_thr/open_thr/watch_thr/weak_thr)は、その端末のログ履歴から
    auto_calibrate_profile()が自動較正した値をそのまま使う（未較正時はBASE_PROFILE既定値）。
    """
    if pd.isna(p):
        return '判定不可'
    if p >= float(profile.get('strong_open_thr', 75.0)):
        return '強い開通'
    if p >= float(profile.get('open_thr', 55.0)):
        return '開通'
    if p >= float(profile.get('watch_thr', 35.0)):
        return '注意監視'
    if p >= float(profile.get('weak_thr', 20.0)):
        return '弱い兆候'
    return '低調'


def cbm2_es_band_state(freq_mhz: float, fxEs_repr: float, profile: Dict,
                        ef_feats: Optional[Dict[str, float]] = None) -> Tuple[str, float, float]:
    """CBM2 Es確率モデルを使ってバンド判定する。固定閾値ではなく、その端末に
    蓄積されているセッションログから auto_calibrate_profile() が学習した
    プロファイル（prob_intercept・F層学習重み・open_thr等の較正閾値）を使う。
    fxEsはMUF∝foEsの関係(周波数比)で対象バンドへ換算してから確率評価する。
    戻り値: (記号付き状態ラベル, 換算fxEs相当値, 確率%)
    """
    if pd.isna(fxEs_repr):
        return "?不明", np.nan, np.nan
    scale = max(freq_mhz, 0.1) / 28.5
    fx_equiv = float(fxEs_repr) / scale
    p = empirical_probability_snapshot(fx_equiv, ef_feats or {}, profile)
    state_full = state_from_prob_snapshot(p, profile)
    return state_symbol(state_full), fx_equiv, p




def _print_es_band_monitor_common(profile: Dict, use_color: bool, android: bool, header_label: str) -> None:
    """ESDUCT(整理版・不具合修正): print_hf_band_bargraph_monitor_silent()/_simple()の
    共通実装。
    【修正内容】旧実装は毎回NICT foF2ページへ問い合わせ、常に取得失敗する
    (fetch_nict_foF2()は仕様上常に空DataFrameを返す)ため、フォールバック代表値
    (foF2=6.5MHz固定, M3000F2=2.9固定)を使ってVOACAP簡易MUF/LUFモデルを動かし、
    全10バンド(1.9〜50MHz)の「開/中/弱」判定を表示していた。この判定は実質的に
    「常に同じ固定入力から算出される固定結果」であり、あたかもリアルタイム判定の
    ように見えるのに実体は無意味だった(F2予測機能自体は既に廃止済みのため)。
    さらに、Esバンド(28.5/50MHz)についてはcbm2_es_band_state()で実測fxEsに基づく
    正しい判定を別途計算していたが、その結果(es_state)を一切表示に使わず、
    同じ固定値ベースのVOACAP判定を表示していた(計算だけして結果を捨てる死んだ
    コードだった)。
    本関数では、F2に依存する処理(NICT foF2取得・VOACAP MUF/LUFモデル・
    フォールバック代表値)を完全に削除し、実測fxEs(NICT Es観測)のみに基づく
    Esバンド(28.5MHz/50MHz)のCBM2判定のみを表示する。F2データが無ければ
    判定できない他の8バンド(1.9〜24MHz)は、もはや判定材料が無いため表示しない。
    """
    W = DISP_W_ANDROID if android else DISP_W_NORMAL
    try:
        loc = ui_resolve_location()
        if loc:
            lat, lon = loc[0], loc[1]
        else:
            lat, lon = NICT_STATION_LATLON.get("kokubunji", (35.6903, 139.4914))
    except Exception:
        lat, lon = NICT_STATION_LATLON.get("kokubunji", (35.6903, 139.4914))
    st = nearest_nict_station(lat, lon)

    fxEs_repr = np.nan
    fxEs_fail_reason = None
    try:
        html = fetch_nict_once(NICT_FXES_URL_DEFAULT)
        fx_df = latest_observation_only(parse_fxes_html(html))
        if fx_df.empty:
            fxEs_fail_reason = "NICT fxEsページ取得/解析0件"
        else:
            row = fx_df[fx_df["station"] == st]
            if row.empty:
                fxEs_repr = safe_float(fx_df["fxEs"].iloc[0])
            else:
                fxEs_repr = safe_float(row["fxEs"].iloc[0])
    except Exception as e:
        fxEs_fail_reason = f"例外:{type(e).__name__}"

    print("-" * W)
    print(f" {header_label} (Esバンドのみ・実測fxEs基準)")
    print("-" * W)
    if pd.isna(fxEs_repr):
        print(f" ※NICT実測fxEsを取得できず判定できません [理由:{fxEs_fail_reason or '不明'}]")
        print("-" * W)
        return
    for freq, band_label in HF_VOACAP_BANDS:
        if freq not in ES_BAND_FREQS:
            continue
        short_label = band_label.split('(')[0].strip()
        try:
            es_state, _fx_equiv, _p = cbm2_es_band_state(freq, fxEs_repr, profile)
            if use_color:
                if "★" in es_state:
                    colored_status = f"{C_RED}{es_state}{C_RESET}"
                elif "●" in es_state:
                    colored_status = f"{C_GREEN}{es_state}{C_RESET}"
                elif "▲" in es_state:
                    colored_status = f"{C_YELLOW}{es_state}{C_RESET}"
                else:
                    colored_status = f"{C_GRAY}{es_state}{C_RESET}"
            else:
                colored_status = es_state
            print(f" {short_label:<7} {colored_status}")
        except Exception:
            print(f" {short_label:<7} ?")
    print("-" * W)


def print_hf_band_bargraph_monitor_silent(profile: Dict, use_color: bool = True, android: bool = False) -> None:
    """Ver12.1: GPS入力を避けた無音バーグラフモニタ（高速版）
    メニュー表示の直前に、データ取得なしでバンド状態を簡潔表示。
    Pydroid3でのユーザー入力をスキップ。
    ESDUCT(整理版): F2関連の残骸削除に伴い、実装は_print_es_band_monitor_common()
    に統一した(詳細はそちらのdocstring参照)。
    """
    _print_es_band_monitor_common(profile, use_color, android, "Esバンド開通状態")


def print_hf_band_bargraph_monitor_simple(profile: Dict, use_color: bool = True, android: bool = False) -> None:
    """Ver12.1: 簡潔版バンド別バーグラフ（Pydroid3用フォールバック）
    ESDUCT(整理版): F2関連の残骸削除に伴い、実装は_print_es_band_monitor_common()
    に統一した(詳細はそちらのdocstring参照)。
    """
    _print_es_band_monitor_common(profile, use_color, android, "Esバンド開通状態")


def fetch_ionogram_image(station: str, timeout: int = SW_FETCH_TIMEOUT) -> Optional[np.ndarray]:
    """NICTイオノグラム最新画像を取得してグレースケール配列で返す。
    キャッシュディレクトリへの保存は行わない（メモリのみ）。
    取得失敗時はNoneを返す。
    """
    url = NICT_IONOGRAM_URLS.get(station)
    if not url:
        return None
    try:
        r = requests.get(url, timeout=timeout,
                         headers={"User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                                                   "AppleWebKit/537.36 (KHTML, like Gecko) "
                                                   "Chrome/120.0 Safari/537.36")})
        r.raise_for_status()
        from io import BytesIO
        img = Image.open(BytesIO(r.content)).convert("L")
        img = img.resize((256, 256))
        return np.asarray(img, dtype=np.float32) / 255.0
    except Exception:
        return None


def analyze_ionogram_image(arr: Optional[np.ndarray],
                            fxEs: float = np.nan,
                            foEs: float = np.nan) -> Dict[str, float]:
    """イオノグラム画像から特徴量を抽出する。
    仕様書:
      Layer Count / Spread Width / Trace Continuity / Reflection Height
      Layer Thickness / Trace Intensity / Diffusion / Saturation / Trace State
    Pillow+numpy実装（OpenCV不要）。画像不可時はfxEs/foEsのみで代替算出。
    """
    result = {
        "layer_count":      np.nan,
        "spread_width":     np.nan,
        "trace_continuity": np.nan,
        "reflection_height":np.nan,
        "layer_thickness":  np.nan,
        "trace_intensity":  np.nan,
        "diffusion":        np.nan,
        "saturation":       np.nan,
        "trace_state":      np.nan,
        "ionogram_quality": 0.0,
    }
    if arr is None:
        if not pd.isna(fxEs):
            result["trace_intensity"] = clamp(fxEs / 15.0, 0.0, 1.0)
            result["trace_state"] = 1.0 if fxEs >= 8.0 else (0.5 if fxEs >= 5.0 else 0.0)
        if not pd.isna(foEs):
            result["saturation"] = clamp(foEs / 10.0, 0.0, 1.0)
        result["ionogram_quality"] = 0.3
        return result

    H, W = arr.shape
    gy = np.gradient(arr, axis=0)
    gx = np.gradient(arr, axis=1)
    edge = np.sqrt(gx**2 + gy**2)

    thr = float(np.quantile(edge, 0.95))
    result["trace_intensity"] = float(thr)

    result["saturation"] = float((arr > 0.85).mean())

    col_edge = edge.mean(axis=0)
    if col_edge.sum() > 0:
        result["spread_width"] = float(col_edge.std() / (W / 4.0))

    row_edge = edge.mean(axis=1)
    if row_edge.sum() > 0:
        cy = float(np.average(np.arange(H), weights=row_edge))
        result["reflection_height"] = float(1.0 - cy / H)

    row_bright = arr.mean(axis=1)
    smooth = np.convolve(row_bright, np.ones(8) / 8, mode="same")
    peaks = 0
    for i in range(1, len(smooth) - 1):
        if smooth[i] > smooth[i-1] and smooth[i] > smooth[i+1] and smooth[i] > 0.4:
            peaks += 1
    result["layer_count"] = float(min(peaks, 5))

    strong_rows = row_edge > float(np.quantile(row_edge, 0.80))
    if strong_rows.sum() > 0:
        rows_idx = np.where(strong_rows)[0]
        result["layer_thickness"] = float((rows_idx[-1] - rows_idx[0]) / H)

    thr_edge = float(np.quantile(edge, 0.80))
    strong_edge_rows = (edge > thr_edge).any(axis=1)
    if len(strong_edge_rows) > 0:
        result["trace_continuity"] = float(strong_edge_rows.mean())

    weak_mask = edge < float(np.quantile(edge, 0.50))
    if weak_mask.sum() > 0:
        result["diffusion"] = float(arr[weak_mask].std())

    cont = result.get("trace_continuity", 0.0) or 0.0
    intn = result.get("trace_intensity", 0.0) or 0.0
    if cont > 0.6 and intn > 0.15:
        result["trace_state"] = 1.0
    elif cont > 0.3 or intn > 0.08:
        result["trace_state"] = 0.5
    else:
        result["trace_state"] = 0.0

    if not pd.isna(fxEs):
        if fxEs >= 10.0:
            result["trace_state"] = max(result["trace_state"], 1.0)
        elif fxEs >= 7.0:
            result["trace_state"] = max(result["trace_state"], 0.5)

    result["ionogram_quality"] = 1.0
    return result


def compute_eqi(fxEs: float, foEs: float,
                ionogram_feats: Dict[str, float]) -> Tuple[float, float]:
    """EQI (Es Quality Index, 0-100) を生成する。
    仕様書生成式: fxEs / foEs / Layer Count / Reflection Quality / Spread / Trace State / Peak Height / Thickness
    戻り値: (eqi, eqi_confidence)
    """
    if pd.isna(fxEs):
        return np.nan, 0.0

    contributions: Dict[str, float] = {}

    contributions["fxEs"] = clamp(fxEs / 15.0, 0.0, 1.0)

    if not pd.isna(foEs):
        contributions["foEs"] = clamp(foEs / 10.0, 0.0, 1.0)
    else:
        contributions["foEs"] = contributions["fxEs"] * 0.8

    contributions["layer_count"]    = clamp((ionogram_feats.get("layer_count", 0.0) or 0.0) / 3.0, 0.0, 1.0)
    contributions["reflect_quality"]= clamp(ionogram_feats.get("trace_continuity", 0.5) or 0.5, 0.0, 1.0)
    contributions["spread"]         = clamp(1.0 - (ionogram_feats.get("spread_width", 0.5) or 0.5), 0.0, 1.0)
    contributions["trace_state"]    = clamp(ionogram_feats.get("trace_state", 0.5) or 0.5, 0.0, 1.0)
    contributions["peak_height"]    = clamp(ionogram_feats.get("reflection_height", 0.5) or 0.5, 0.0, 1.0)
    contributions["thickness"]      = clamp(1.0 - (ionogram_feats.get("layer_thickness", 0.3) or 0.3), 0.0, 1.0)

    w = {
        "fxEs":           EQI_W_FXES,
        "foEs":           EQI_W_FOES,
        "layer_count":    EQI_W_LAYER_CNT,
        "reflect_quality":EQI_W_REFLECT,
        "spread":         EQI_W_SPREAD,
        "trace_state":    EQI_W_TRACE,
        "peak_height":    0.0,
        "thickness":      0.0,
    }
    eqi_raw = sum(contributions[k] * w[k] for k in w if k in contributions)
    eqi = float(clamp(eqi_raw * 100.0, 0.0, 100.0))

    iq = float(ionogram_feats.get("ionogram_quality", 0.3) or 0.3)
    eqi_confidence = iq

    return eqi, eqi_confidence


def estimate_es_phase(fxEs: float, eqi: float, prev_eqi: float,
                      prev_fxEs: float, sw_score: float,
                      d_fxEs: float) -> str:
    """Es Phase を推定する。
    Phase: FORMATION → GROWTH → PEAK → DECAY → COLLAPSE
    仕様書: ΔfxEs / ΔEQI / ΔReflection / ΔSpaceWeather を入力とする。
    """
    if pd.isna(fxEs):
        return "FORMATION"

    d_eqi = (eqi - prev_eqi) if (not pd.isna(eqi) and not pd.isna(prev_eqi)) else 0.0
    d_fx  = d_fxEs if not pd.isna(d_fxEs) else 0.0
    sw    = sw_score if not pd.isna(sw_score) else 0.5

    if fxEs < ES_PHASE_COLLAPSE_FXES and d_fx < -2.0:
        return "COLLAPSE"

    if d_eqi < ES_PHASE_DECAY_DEQI * 100 or (fxEs < ES_PHASE_GROWTH_FXES and d_fx < -0.5):
        return "DECAY"

    if fxEs >= ES_PHASE_PEAK_FXES and abs(d_eqi) < 5.0:
        return "PEAK"

    if fxEs >= ES_PHASE_GROWTH_FXES or (d_fx > 0.5 and d_eqi > 2.0):
        return "GROWTH"

    if fxEs >= ES_PHASE_FORMATION_FXES:
        return "FORMATION"

    return "FORMATION"


def fetch_and_update_eqi(station: str, fxEs: float, foEs: float = np.nan,
                          fetch_image: bool = True) -> Tuple[float, float, Dict]:
    """1観測局のEQIを更新する（イオノグラム画像取得→解析→EQI生成）。
    fetch_image=Falseの場合は画像取得をスキップし既存データで代理算出。
    戻り値: (eqi, eqi_confidence, ionogram_feats)
    """
    if fetch_image:
        arr = fetch_ionogram_image(station)
    else:
        arr = None
    ionogram_feats = analyze_ionogram_image(arr, fxEs=fxEs, foEs=foEs)
    eqi, eqi_conf = compute_eqi(fxEs, foEs, ionogram_feats)
    _eqi_history.setdefault(station, [])
    _eqi_history[station].append({
        "obs_jst": now_jst().isoformat(),
        "eqi": eqi, "eqi_confidence": eqi_conf,
        **ionogram_feats,
    })
    if len(_eqi_history[station]) > SW_HISTORY_MAXLEN:
        _eqi_history[station] = _eqi_history[station][-SW_HISTORY_MAXLEN:]
    return eqi, eqi_conf, ionogram_feats


def get_prev_eqi(station: str) -> float:
    """1サイクル前のEQIをRAM履歴から取得する。"""
    hist = _eqi_history.get(station, [])
    if len(hist) >= 2:
        return float(hist[-2].get("eqi", np.nan))
    return np.nan


def update_es_phase_all(sample: pd.DataFrame, sw_score: float,
                         dyn: pd.DataFrame) -> None:
    """全観測局のEs Phaseを推定してグローバル状態を更新する。"""
    for _, row in sample.iterrows():
        st = row["station"]
        fxEs = float(row["fxEs"]) if not pd.isna(row.get("fxEs")) else np.nan
        eqi = _latest_eqi_by_station.get(st, {}).get("eqi", np.nan)
        prev_eqi = get_prev_eqi(st)
        d_fx = np.nan
        if not dyn.empty:
            g = dyn[dyn["station"] == st]
            if not g.empty:
                d_fx = float(g.sort_values("obs_dt").iloc[-1].get("d_eff", np.nan))
        phase = estimate_es_phase(fxEs, eqi, prev_eqi, np.nan, sw_score, d_fx)
        if len(_es_phase_history.get(st, [])) < ES_PHASE_MIN_SAMPLES:
            phase = "FORMATION"
        _latest_es_phase_by_station[st] = phase
        _latest_es_phase_num_by_station[st] = _ES_PHASE_NUM.get(phase, 0.0)
        _es_phase_history.setdefault(st, [])
        _es_phase_history[st].append({
            "obs_jst": now_jst().isoformat(),
            "phase": phase,
            "phase_num": _ES_PHASE_NUM.get(phase, 0.0),
        })
        if len(_es_phase_history[st]) > SW_HISTORY_MAXLEN:
            _es_phase_history[st] = _es_phase_history[st][-SW_HISTORY_MAXLEN:]



def compute_interaction_terms(feat_dict: Dict[str, float]) -> Dict[str, float]:
    """仕様書定義の交互作用項を生成する。
    交互作用項 = A × B (両者をclamp済みの0-1正規化値で計算)。
    """
    pairs = INTERACTION_PAIRS
    out: Dict[str, float] = {}
    for a, b in pairs:
        va = feat_dict.get(a, np.nan)
        vb = feat_dict.get(b, np.nan)
        key = f"ia_{a.replace('_norm','')}_{b.replace('_norm','').replace('_approach_index','mstid')}"
        key_map = {
            "ia_eqi_fxes":  ("eqi_norm",     "fxEs_norm"),
            "ia_eqi_mstid": ("eqi_norm",     "mstid_approach_index"),
            "ia_sw_fxes":   ("sw_score",     "fxEs_norm"),
            "ia_sw_mstid":  ("sw_score",     "mstid_approach_index"),
            "ia_eqi_sw":    ("eqi_norm",     "sw_score"),
        }
        for ckey, (ca, cb) in key_map.items():
            if a == ca and b == cb:
                va2 = feat_dict.get(ca, np.nan)
                vb2 = feat_dict.get(cb, np.nan)
                if pd.isna(va2) or pd.isna(vb2):
                    out[ckey] = 0.0
                else:
                    out[ckey] = float(va2) * float(vb2)
                break
    return out


def build_lag_features(history_dict: Dict[str, List], key: str,
                        lag_steps: List[int] = LAG_STEPS) -> Dict[str, float]:
    """指定キーの時系列履歴からラグ特徴量(t, t-15, t-30, t-45, t-60)を生成する。
    仕様書: ラグ特徴量 t/t-15/t-30/t-45/t-60 を追加。
    """
    out: Dict[str, float] = {}
    hist = history_dict if isinstance(history_dict, list) else []
    for step in lag_steps:
        idx = -(1 + step)
        try:
            val = hist[idx].get(key, np.nan) if len(hist) > step else np.nan
        except (IndexError, AttributeError):
            val = np.nan
        suffix = f"_t" if step == 0 else f"_t{step*15}m"
        out[f"{key}{suffix}"] = float(val) if not pd.isna(val) else np.nan
    return out



def permutation_importance_check(profile: Dict,
                                   valid_df: pd.DataFrame,
                                   future_col: str = "future30") -> List[str]:
    """Permutation Importanceで寄与率が低い特徴量を特定して返す。
    仕様書: 寄与率が一定以下なら学習対象から自動除外する。
    対象: Permutation Importance / Ridge係数 / VIF
    戻り値: 除外推奨特徴量名リスト
    """
    exclude: List[str] = []
    if not _SKLEARN_AVAILABLE:
        return exclude
    if valid_df.empty or future_col not in valid_df.columns:
        return exclude

    feat_names = [fp[0] for fp in _F_FEAT_RIDGE
                  if fp[0] in valid_df.columns]
    if not feat_names:
        return exclude

    sub = valid_df[feat_names + [future_col]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 30:
        return exclude

    X = sub[feat_names].values
    y = sub[future_col].values
    try:
        scaler = StandardScaler()
        Xs = scaler.fit_transform(X)
        ridge = Ridge(alpha=RIDGE_ALPHA)
        ridge.fit(Xs, y)
        base_score = float(np.corrcoef(ridge.predict(Xs), y)[0, 1] ** 2)
        if base_score < 0.01:
            return exclude

        rng = np.random.default_rng(seed=42)
        importances: Dict[str, float] = {}
        for i, feat in enumerate(feat_names):
            Xs_perm = Xs.copy()
            Xs_perm[:, i] = rng.permutation(Xs_perm[:, i])
            perm_score = float(np.corrcoef(ridge.predict(Xs_perm), y)[0, 1] ** 2)
            importances[feat] = max(0.0, base_score - perm_score)

        total_imp = sum(importances.values())
        if total_imp <= 0:
            return exclude

        for feat, imp in importances.items():
            ratio = imp / total_imp
            if ratio < PERM_IMP_MIN_RATIO:
                exclude.append(feat)
    except Exception:
        pass
    return exclude



PRED_FEATURE_CSV = "prediction_feature.csv"
TEACHER_CSV      = "teacher.csv"
PREDICTION_CSV   = "prediction.csv"

_PRED_FEATURE_EXTRA_COLS = [
    "space_weather_score", "eqi", "reflection_quality",
    "es_phase", "prediction_confidence",
]
_TEACHER_EXTRA_COLS = [
    "eqi", "es_phase", "space_weather_score", "reflection_quality",
]
_PREDICTION_EXTRA_COLS = [
    "prediction_confidence", "eqi", "space_weather_score",
]


def compute_prediction_confidence(profile: Dict, ft8_feats: Dict,
                                   eqi: float, eqi_confidence: float,
                                   sw_score: float, trace_quality: float) -> float:
    """Prediction Confidence (0-100) を算出する。
    仕様書: 既存のモデル成熟度・取得成功率に加えてEQI信頼度・SW品質・Trace品質を反映。
    """
    base = float(_model_score) * 0.50

    eqi_c = clamp(float(eqi_confidence) if not pd.isna(eqi_confidence) else 0.3, 0.0, 1.0)
    base += eqi_c * 20.0

    sw_ok = not pd.isna(sw_score)
    base += (15.0 if sw_ok else 5.0)

    tq = clamp(float(trace_quality) if not pd.isna(trace_quality) else 0.3, 0.0, 1.0)
    base += tq * 15.0

    conf = ft8_feats.get("confidence_label", "LOW")
    base += {"HIGH": 10.0, "MEDIUM": 5.0, "LOW": 0.0}.get(conf, 0.0)

    profile["eqi_confidence"]  = float(eqi_c)
    profile["sw_quality"]      = 1.0 if sw_ok else 0.0
    profile["trace_quality"]   = float(tq)

    return clamp(base, 0.0, 100.0)


def append_to_prediction_feature_csv(log_dir: Path, row: Dict) -> None:
    """prediction_feature.csv に1行を追記する（新規CSV禁止・既存統合）。"""
    p = log_dir / PRED_FEATURE_CSV
    _csv_append(p, pd.DataFrame([row]))


def append_to_teacher_csv(log_dir: Path, row: Dict) -> None:
    """teacher.csv に1行を追記する（新規CSV禁止・既存統合）。"""
    p = log_dir / TEACHER_CSV
    _csv_append(p, pd.DataFrame([row]))


def append_to_prediction_csv(log_dir: Path, row: Dict) -> None:
    """prediction.csv に1行を追記する（新規CSV禁止・既存統合）。"""
    p = log_dir / PREDICTION_CSV
    _csv_append(p, pd.DataFrame([row]))


_fof2_history: Dict[str, List[Dict]] = {}
_run_cycle_mock_mode: bool = False
_FOF2_HISTORY_MAXLEN = 20


def update_fof2_history(foF2_df: pd.DataFrame) -> None:
    """fetch_nict_foF2()の結果をセッション内RAM履歴に追記する。"""
    global _fof2_history
    if foF2_df.empty:
        return
    for _, row in foF2_df.iterrows():
        st = row.get("station")
        if st is None:
            continue
        if st not in _fof2_history:
            _fof2_history[st] = []
        entry = {
            "obs_jst": row.get("obs_jst"),
            "foF2": safe_float(row.get("foF2")),
            "hmF2": safe_float(row.get("hmF2")),
            "M3000F2": safe_float(row.get("M3000F2")),
        }
        if any(e["obs_jst"] == entry["obs_jst"] for e in _fof2_history[st]):
            continue
        _fof2_history[st].append(entry)
        if len(_fof2_history[st]) > _FOF2_HISTORY_MAXLEN:
            _fof2_history[st] = _fof2_history[st][-_FOF2_HISTORY_MAXLEN:]


def compute_f_region_background(obs_latest: pd.DataFrame, station: str = "", fxEs: float = np.nan) -> Dict[str, float]:
    """Ver10.5J: F層背景特徴量。
    foF2, hmF2, M3000F2 の現在値と時間微分(dfoF2/dt, dhmF2/dt)を返す。

    Ver10.5J変更点:
      - masking_risk_index を foF2/27MHz から foF2/fxEs に変更。
        foF2>fxEs(>1.0) でF層がEsを上回る = 物理的マスキング状態。
      - foF2_es_ratio = foF2/fxEs を特徴量として追加。
      - M3000F2を正規化して特徴量として追加(m3000f2_norm = M3000F2/3.0)。
      - fxEs=NaN時は旧互換で foF2/27MHz にフォールバック。

    返り値キー:
        foF2               : 最新foF2 [MHz]
        foF2_ratio         : foF2 / 27.0 (旧互換用、以降は foF2_es_ratio を推奨)
        foF2_es_ratio      : foF2 / fxEs (物理的マスキング指標、>1.0=マスキング状態)
        dfoF2              : dfoF2/dt [MHz/h]
        hmF2               : 最新hmF2 [km]
        dhmF2              : dhmF2/dt [km/h] (負=下降=Es予兆)
        M3000F2            : 最新M3000F2 [-]
        m3000f2_norm       : M3000F2 / 3.0 (正規化済み、学習特徴量用)
        masking_risk_index : foF2_es_ratio (foF2/fxEs)。fxEs不明時は foF2_ratio
        f_region_background_index : dfoF2
        ef_coupling_flag   : dhmF2 < -5 km/h のとき1.0
    """
    nan_result = {
        'foF2': np.nan, 'foF2_ratio': np.nan, 'foF2_es_ratio': np.nan, 'dfoF2': np.nan,
        'hmF2': np.nan, 'dhmF2': np.nan, 'M3000F2': np.nan, 'm3000f2_norm': np.nan,
        'f_region_background_index': np.nan,
        'ef_coupling_flag': 0.0,
        'masking_risk_index': np.nan,
    }
    hist = _fof2_history.get(station, [])
    if not hist:
        return nan_result

    try:
        hist_sorted = sorted(hist, key=lambda x: x["obs_jst"])
    except Exception:
        return nan_result

    latest = hist_sorted[-1]
    foF2    = latest.get("foF2",    np.nan)
    hmF2    = latest.get("hmF2",    np.nan)
    M3000F2 = latest.get("M3000F2", np.nan)

    dfoF2 = np.nan
    dhmF2 = np.nan
    if len(hist_sorted) >= 2:
        prev = hist_sorted[-2]
        try:
            t1 = datetime.fromisoformat(str(latest["obs_jst"]))
            t0 = datetime.fromisoformat(str(prev["obs_jst"]))
            dt_h = (t1 - t0).total_seconds() / 3600.0
            if dt_h > 0:
                pf = prev.get("foF2", np.nan)
                ph = prev.get("hmF2", np.nan)
                if not (pd.isna(foF2) or pd.isna(pf)):
                    dfoF2 = float(np.clip((foF2 - pf) / dt_h, -20.0, 20.0))
                if not (pd.isna(hmF2) or pd.isna(ph)):
                    dhmF2 = float(np.clip((hmF2 - ph) / dt_h, -500.0, 500.0))
        except Exception:
            pass

    foF2_ratio = (float(foF2) / TARGET_27_MHZ) if not pd.isna(foF2) else np.nan

    FXES_FLOOR_MHZ = 5.0
    if not (pd.isna(foF2) or pd.isna(fxEs)):
        fxEs_safe = max(float(fxEs), FXES_FLOOR_MHZ)
        foF2_es_ratio = float(np.clip(float(foF2) / fxEs_safe, 0.0, 5.0))
        masking_risk_index = foF2_es_ratio
    else:
        foF2_es_ratio = np.nan
        masking_risk_index = foF2_ratio

    m3000f2_norm = float(float(M3000F2) / 3.0) if not pd.isna(M3000F2) else np.nan

    ef_coupling_flag = 1.0 if (not pd.isna(dhmF2) and dhmF2 < -5.0) else 0.0

    return {
        'foF2':                    foF2 if not pd.isna(foF2) else np.nan,
        'foF2_ratio':              foF2_ratio,
        'foF2_es_ratio':           foF2_es_ratio,
        'dfoF2':                   dfoF2,
        'hmF2':                    hmF2 if not pd.isna(hmF2) else np.nan,
        'dhmF2':                   dhmF2,
        'M3000F2':                 M3000F2 if not pd.isna(M3000F2) else np.nan,
        'm3000f2_norm':            m3000f2_norm,
        'f_region_background_index': dfoF2,
        'ef_coupling_flag':        ef_coupling_flag,
        'masking_risk_index':      masking_risk_index,
    }


def classify_area_by_callsign(callsign: Optional[str]) -> Optional[str]:
    """受信局コールサインから日本のアマチュア無線エリア(0-9/JR6/BV)を判定する。
    lat/lon逆ジオコーディングではなく、コールサインに埋め込まれたエリア数字を
    直接利用する（PSKReporterのreceiverLocatorはグリッド中心精度しかなく、
    JCC/エリア境界判定には不向きなため）。
    判定不能（海外局・特殊プレフィックス等）の場合はNoneを返す。
    """
    if not callsign:
        return None
    cs = str(callsign).strip().upper().split('/')[0]
    if not cs:
        return None
    if re.match(r'^(BV|BX|BM|BU|BO)', cs):
        return "BV"
    if re.match(r'^(JR6|JS6)', cs):
        return "JR6"
    m = re.match(r'^(J[A-S]|7[J-N]|8[J-N])([0-9])[A-Z]', cs)
    if m and m.group(2) in AREA_NAME_JP:
        return m.group(2)
    return None


def maidenhead_to_latlon(locator: Optional[str]) -> Tuple[float, float]:
    """Maidenhead Grid Locator(4桁/6桁)を緯度経度[度]に変換する。
    6桁未満の場合はその升目の中心を代表点として返す。解釈不能時はNaN。
    """
    if not locator or not isinstance(locator, str) or len(locator) < 4:
        return np.nan, np.nan
    loc = locator.strip().upper()
    try:
        A = ord('A')
        lon = (ord(loc[0]) - A) * 20.0 - 180.0
        lat = (ord(loc[1]) - A) * 10.0 - 90.0
        lon += int(loc[2]) * 2.0
        lat += int(loc[3]) * 1.0
        if len(loc) >= 6 and loc[4].isalpha() and loc[5].isalpha():
            lon += (ord(loc[4]) - A) * (2.0 / 24.0)
            lat += (ord(loc[5]) - A) * (1.0 / 24.0)
            lon += (2.0 / 24.0) / 2.0
            lat += (1.0 / 24.0) / 2.0
        else:
            lon += 1.0
            lat += 0.5
        return lat, lon
    except Exception:
        return np.nan, np.nan


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    if any(pd.isna(x) for x in (lat1, lon1, lat2, lon2)):
        return np.nan
    R = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return R * 2.0 * math.atan2(math.sqrt(a), math.sqrt(max(0.0, 1.0 - a)))


def bearing_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    if any(pd.isna(x) for x in (lat1, lon1, lat2, lon2)):
        return np.nan
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dlambda = math.radians(lon2 - lon1)
    x = math.sin(dlambda) * math.cos(p2)
    y = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dlambda)
    return (math.degrees(math.atan2(x, y)) + 360.0) % 360.0


def fetch_pskreporter_spots(callsign: str, lookback_sec: int,
                             freq_range: Tuple[int, int] = PSK_FREQ_RANGE_HZ,
                             mode: str = PSK_MODE,
                             timeout: int = HTTP_TIMEOUT_SEC
                             ) -> Tuple[pd.DataFrame, bool, str]:
    """PSKReporterから直近lookback_sec秒の受信レポートを取得する。
    戻り値: (DataFrame, success:bool, error_reason:str)
    - success=True, error_reason=''  → 通信成功・パース成功（スポット0件も含む）
    - success=False, error_reason=...→ HTTP障害・タイムアウト・XMLパース失敗等
    呼び出し側でsuccessフラグを見ることで「通信障害」と「本当に0局受信」を区別し、
    整合性判定・lag記録への偽陰性混入を防ぐ（Ver11.0J改）。
    columns(success時): receiver_callsign, receiver_locator, sender_locator,
                        freq_hz, snr, obs_jst
    """
    params = {
        "senderCallsign": callsign,
        "flowStartSeconds": -abs(int(lookback_sec)),
        "frange": f"{freq_range[0]}-{freq_range[1]}",
        "mode": mode,
        "rronly": 1,
    }
    headers = {"User-Agent": "CBM2.py Es predictor FT8 fusion; JL7KHN one GET per cycle"}
    try:
        r = requests.get(PSKREPORTER_QUERY_URL, params=params, timeout=timeout, headers=headers)
        r.raise_for_status()
        xml_text = r.text
    except requests.exceptions.Timeout:
        return pd.DataFrame(), False, f"timeout({timeout}s)"
    except requests.exceptions.HTTPError as e:
        return pd.DataFrame(), False, f"HTTP {e.response.status_code}"
    except requests.exceptions.ConnectionError as e:
        return pd.DataFrame(), False, f"connection_error: {str(e)[:60]}"
    except Exception as e:
        return pd.DataFrame(), False, f"fetch_error: {str(e)[:60]}"
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        return pd.DataFrame(), False, f"xml_parse_error: {str(e)[:60]}"
    records = []
    for rep in root.findall(".//receptionReport"):
        recv_cs = rep.get("receiverCallsign")
        if not recv_cs:
            continue
        flow = rep.get("flowStartSeconds")
        try:
            obs_dt = datetime.fromtimestamp(int(flow), tz=timezone.utc).astimezone(JST)
        except Exception:
            obs_dt = now_jst()
        records.append({
            "receiver_callsign": recv_cs,
            "receiver_locator": rep.get("receiverLocator"),
            "sender_locator": rep.get("senderLocator"),
            "freq_hz": safe_float(rep.get("frequency")),
            "snr": safe_float(rep.get("sNR")),
            "obs_jst": obs_dt.isoformat(),
        })
    return pd.DataFrame(records), True, ''


def enrich_psk_spots_with_geometry(spots: pd.DataFrame) -> pd.DataFrame:
    """受信局の方位・距離（送信地点=JL7KHN/Pの現在地基準）とエリアを付与する。
    送信地点はPSKReporterが返すsenderLocator（その時点でのポータブル運用地）
    から動的に決定する。複数の値が混在する場合は最も出現頻度の高いものを採用。
    """
    if spots is None or spots.empty:
        return spots
    s = spots.copy()
    s["area"] = s["receiver_callsign"].map(classify_area_by_callsign)
    sender_locators = s["sender_locator"].dropna()
    if sender_locators.empty:
        s["distance_km"] = np.nan
        s["bearing_deg"] = np.nan
        return s
    home_loc = sender_locators.mode().iloc[0]
    home_lat, home_lon = maidenhead_to_latlon(home_loc)
    if pd.isna(home_lat):
        s["distance_km"] = np.nan
        s["bearing_deg"] = np.nan
        return s
    latlon = s["receiver_locator"].map(maidenhead_to_latlon)
    s["recv_lat"] = latlon.map(lambda t: t[0])
    s["recv_lon"] = latlon.map(lambda t: t[1])
    s["distance_km"] = s.apply(lambda row: haversine_km(home_lat, home_lon, row["recv_lat"], row["recv_lon"]), axis=1)
    s["bearing_deg"] = s.apply(lambda row: bearing_deg(home_lat, home_lon, row["recv_lat"], row["recv_lon"]), axis=1)
    return s


_AREA_CODE_NORMALIZE = {
    "okinawa": "JR6", "沖縄": "JR6", "jr6": "JR6",
    "taiwan":  "BV",  "台湾": "BV",  "bv":  "BV",
}

def _builtin_fixed_monitor_df() -> pd.DataFrame:
    """Ver11.1J: 組み込み固定観測局DBをDataFrameに変換して返す。"""
    rows = [{"callsign": cs, "lat": lat, "lon": lon, "area": area}
            for cs, lat, lon, area in BUILTIN_FIXED_MONITOR]
    return pd.DataFrame(rows)


def load_fixed_monitor_db(path: Path) -> pd.DataFrame:
    """固定観測局DBを読み込む。優先順:
    1. pathのCSVが存在しフォーマット合格 → CSVを使用
    2. それ以外 → BUILTIN_FIXED_MONITOR（組み込みDB）を使用（Ver11.1J）
    組み込みDBにより、CSVなしで到達率/RPIが最初から有効になる。
    エリアコード正規化: "Okinawa"/"Taiwan"等を"JR6"/"BV"に変換する。
    """
    def _normalize_area(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["area"] = df["area"].astype(str).map(
            lambda x: _AREA_CODE_NORMALIZE.get(x.strip().lower(),
                       _AREA_CODE_NORMALIZE.get(x.strip(), x.strip()))
        )
        return df

    if path.exists():
        try:
            df = pd.read_csv(path)
            req = {"callsign", "lat", "lon", "area"}
            if req.issubset(set(df.columns)):
                return _normalize_area(df)
        except Exception:
            pass
    return _builtin_fixed_monitor_df()


def fixed_monitor_counts_by_area(db: pd.DataFrame) -> Dict[str, int]:
    if db is None or db.empty:
        return {}
    return db.groupby("area")["callsign"].nunique().to_dict()


def compute_area_reachability(spots_enriched: pd.DataFrame, fixed_db: pd.DataFrame) -> pd.DataFrame:
    """エリア毎の到達率（受信局数/固定観測局数）とSNR加重和を算出する（仕様書5章）。
    固定局DB未整備（該当エリアの登録局数0）の場合、reach_ratioはNaNとなり、
    呼び出し側は絶対受信局数(heard_count)のみで表示する。
    area_confidence列（Ver11.0J追加）:
      "NODATA" = 固定局0（到達率算出不可）
      "LOW"    = 固定局1〜MIN_FIXED_STATIONS_FOR_CONFIDENCE-1（到達率は出るが
                 N数が小さく統計的にブレやすい。例: 固定局1局で受信1局なら
                 機械的には100%になるが、これは「その1局がたまたま受信できた」
                 以上の意味を持たない可能性が高い）
      "OK"     = 固定局がMIN_FIXED_STATIONS_FOR_CONFIDENCE以上登録済み
    戻り値columns: area, area_name, heard_count, fixed_count, reach_ratio, weight_sum, area_confidence
    """
    fixed_counts = fixed_monitor_counts_by_area(fixed_db)
    heard_by_area: Dict[str, set] = {}
    weight_by_area: Dict[str, float] = {}
    if spots_enriched is not None and not spots_enriched.empty and "area" in spots_enriched.columns:
        s = spots_enriched.dropna(subset=["area"]).drop_duplicates(["area", "receiver_callsign"])
        if not s.empty:
            heard_by_area = {a: set(g["receiver_callsign"]) for a, g in s.groupby("area")}
            w = s["snr"].map(lambda v: sigmoid(v) if not pd.isna(v) else 0.5)
            weight_by_area = s.assign(__w__=w).groupby("area")["__w__"].sum().to_dict()
    rows = []
    for area in AREA_DISPLAY_ORDER:
        heard = len(heard_by_area.get(area, set()))
        fixed_n = int(fixed_counts.get(area, 0))
        ratio = clamp(heard / fixed_n, 0.0, 1.0) if fixed_n > 0 else np.nan
        if fixed_n <= 0:
            confidence = "NODATA"
        elif fixed_n < MIN_FIXED_STATIONS_FOR_CONFIDENCE:
            confidence = "LOW"
        else:
            confidence = "OK"
        rows.append({
            "area": area, "area_name": format_area_label(area),
            "heard_count": heard, "fixed_count": fixed_n,
            "reach_ratio": ratio, "weight_sum": float(weight_by_area.get(area, 0.0)),
            "area_confidence": confidence,
        })
    return pd.DataFrame(rows)


def interpolate_area_reachability_idw(area_df: pd.DataFrame) -> pd.DataFrame:
    """Rev14.0: 固定観測局が無い/少ないエリアの到達率を、周辺エリアの実測値
    から逆距離加重(IDW)で補完する(提案②(a): 観測点が無い場所の空間推定)。

    設計方針(学習の健全性を守るための境界線):
      - ここで生成する reach_ratio_compensated は「現在状態の表示・入力
        特徴量」としてのみ使用する想定。resolve_pending_samples() が読む
        生の reach_ratio 列はこの関数では一切変更しない。つまりRidge
        モデルの学習ラベル(未来の実測到達率)は常に本物の観測値のみを使い、
        IDW推定値を教師データとして学習させる循環参照は起きない
        (呼び出し側でも reach_ratio 列自体は変更せず、新規列のみ追加する)。
      - 補間元は area_confidence == 'OK' のエリアの実測値のみを使う。
        LOW/NODATAのエリアを補間元として連鎖的に使うことはしない
        (誤差の増幅・伝播を防ぐため)。
      - AREA_IDW_MAX_NEIGHBOR_KM を超える遠方エリアは、Es伝播が地域的な
        現象であるという物理的前提から補間元に含めない。
      - NODATA(固定局0局、reach_ratioがNaN)は全面的にIDW推定値へ置換する。
      - LOW(固定局1局以上だがMIN_FIXED_STATIONS_FOR_CONFIDENCE未満)は、
        実測値とIDW推定値を局数に応じて重み付けブレンドし、小サンプル
        特有の「0%/100%への振れ」を緩和する。
      - 信頼できる観測源(OK確度のエリア)自体が全国的に不足している場合
        (AREA_IDW_MIN_SOURCES未満)は、無理に補間せず「補間不可」のまま
        とする(存在しないデータから確信度高く推定した体裁を作らない)。

    戻り値: area_df に以下の列を追加したコピー
      reach_ratio_compensated : 表示・特徴量入力用の補完済み到達率(0-1)
      is_compensated          : 補間/ブレンドを行ったか(bool)
      compensation_basis      : 'raw' / 'idw_interp' / 'idw_blend' / 'insufficient'
      compensation_n_sources  : 補間に使った近傍エリア数
    """
    if area_df is None or area_df.empty:
        return area_df

    df = area_df.copy()
    df["reach_ratio_compensated"] = df["reach_ratio"]
    df["is_compensated"] = False
    df["compensation_basis"] = "raw"
    df["compensation_n_sources"] = 0

    source_rows = []
    for _, row in df.iterrows():
        if row.get("area_confidence") == "OK" and not pd.isna(row.get("reach_ratio")):
            center = AREA_CENTERS.get(row["area"])
            if center is not None:
                source_rows.append((row["area"], center[0], center[1], float(row["reach_ratio"])))

    if len(source_rows) < AREA_IDW_MIN_SOURCES:
        return df

    for idx, row in df.iterrows():
        area = row["area"]
        conf = row.get("area_confidence")
        raw_ratio = row.get("reach_ratio")
        center = AREA_CENTERS.get(area)
        if center is None:
            continue
        alat, alon = center

        weights = []
        values = []
        for src_area, slat, slon, sval in source_rows:
            if src_area == area:
                continue
            dist = haversine_km(alat, alon, slat, slon)
            if pd.isna(dist) or dist > AREA_IDW_MAX_NEIGHBOR_KM:
                continue
            dist = max(dist, 1.0)
            weights.append(1.0 / (dist ** AREA_IDW_POWER))
            values.append(sval)

        n_sources = len(values)
        idw_estimate = float(np.average(values, weights=weights)) if n_sources >= AREA_IDW_MIN_SOURCES else np.nan

        if conf == "NODATA":
            if not pd.isna(idw_estimate):
                df.at[idx, "reach_ratio_compensated"] = idw_estimate
                df.at[idx, "is_compensated"] = True
                df.at[idx, "compensation_basis"] = "idw_interp"
                df.at[idx, "compensation_n_sources"] = n_sources
            else:
                df.at[idx, "compensation_basis"] = "insufficient"
        elif conf == "LOW" and not pd.isna(raw_ratio) and not pd.isna(idw_estimate):
            fixed_n = float(row.get("fixed_count", 0) or 0)
            w_local = clamp(fixed_n / MIN_FIXED_STATIONS_FOR_CONFIDENCE, 0.0, 1.0)
            blended = w_local * float(raw_ratio) + (1.0 - w_local) * idw_estimate
            df.at[idx, "reach_ratio_compensated"] = clamp(blended, 0.0, 1.0)
            df.at[idx, "is_compensated"] = True
            df.at[idx, "compensation_basis"] = "idw_blend"
            df.at[idx, "compensation_n_sources"] = n_sources

    return df


def compute_rpi(area_df: pd.DataFrame) -> float:
    """Realtime Propagation Index (0-100、仕様書5章)。
    固定局DBに基づき到達率が算出できたエリアのみを対象に、SNR加重和を
    重みとした加重平均を取る。固定局DB未整備で対象エリアが無ければNaN。
    """
    if area_df is None or area_df.empty:
        return np.nan
    valid = area_df.dropna(subset=["reach_ratio"])
    if valid.empty:
        return np.nan
    w = valid["weight_sum"].clip(lower=0.1)
    rpi = float(np.average(valid["reach_ratio"], weights=w) * 100.0)
    return clamp(rpi, 0.0, 100.0)


def judge_propagation_mode(masking_risk: float, es_hop_ratio: float, accel: float) -> Tuple[str, str, str, float]:
    """伝播モード判定（Es優勢/F層優勢/混合、仕様書6章）。
    暫定的なルールベース判定であり、ES_MODE_SCORE_GATE等のしきい値は
    実データ蓄積後の再チューニングが前提（Ver10.8Jのlag_stability_score
    導入時と同様の立ち位置）。
    戻り値: (mode_code['E'/'F'/'M'], mode_label, sub_label, es_mode_score[-1..+1])
    """
    if pd.isna(masking_risk) and pd.isna(es_hop_ratio):
        return 'M', '混合', '判定材料不足', 0.0
    masking = 0.0 if pd.isna(masking_risk) else clamp(float(masking_risk), 0.0, 3.0)
    es_signal = 0.0 if pd.isna(es_hop_ratio) else float(es_hop_ratio)
    f_score = min(1.0, masking)
    e_score = es_signal * (1.0 - min(masking, 1.0))
    es_mode_score = clamp(e_score - f_score, -1.0, 1.0)
    if es_mode_score > ES_MODE_SCORE_GATE:
        sub = '減衰中' if (not pd.isna(accel) and accel < -0.1) else '成長中'
        return 'E', 'Es優勢', sub, es_mode_score
    if -es_mode_score > ES_MODE_SCORE_GATE:
        return 'F', 'F層優勢', '安定', es_mode_score
    return 'M', '混合', '遷移中', es_mode_score


def strength_label_es(x_norm: float) -> str:
    """ratio_label()の5段階(強/中/弱/微/低)を流用し「強Es」等の表記にする（仕様書10章）。"""
    if x_norm is None or pd.isna(x_norm):
        return "不明"
    return f"{ratio_label(x_norm)}Es"


def compute_consistency(fxEs_repr: float, rpi: float) -> Tuple[str, str, float, bool]:
    """NICT観測値(fxEs)とFT8実伝播(RPI)の整合性判定（仕様書10章）。
    両者を0-1正規化しratio_label()の5段階で比較する。15MHzをNICT側の
    「強Es」目安として正規化に用いる（暫定値、F層特徴量と同様の性質）。
    戻り値: (nict_label, ft8_label, consistency[0-1], match)
    """
    if pd.isna(fxEs_repr) or pd.isna(rpi):
        return "不明", "不明", np.nan, False
    nict_norm = clamp(float(fxEs_repr) / 15.0, 0.0, 1.0)
    ft8_norm = clamp(float(rpi) / 100.0, 0.0, 1.0)
    nict_label = strength_label_es(nict_norm)
    ft8_label = strength_label_es(ft8_norm)
    consistency = clamp(1.0 - abs(nict_norm - ft8_norm), 0.0, 1.0)
    return nict_label, ft8_label, consistency, (nict_label == ft8_label)


def compute_distance_distribution_features(spots_enriched: pd.DataFrame) -> Dict[str, float]:
    """FT8受信スポットの距離・方位分布特徴量を算出する（仕様書将来項目、タカユキ指摘）。
    現在のes_hop_ratio（窓内比率の二値判定）は
      [900,1200,1500km] vs [900,1200,1500,5000,7000km]
    をほぼ同じEs強と判定してしまう問題がある。
    分布の形状（σ・ピーク数・方位集中度）を追加することで
    「帯状分布=Es単独」「広域散乱=F寄り」「双峰=Es+F複合」の
    将来的な判別に使えるデータを蓄積する。
    現時点では学習器入力にはせず、corr_csvへのR²追記用途のみ（Ver11.0J）。

    返却キー:
      es_hop_ratio         : 600-2400km窓内比率（既存、ここで正式算出）
      dist_in_window_count : 窓内スポット絶対数
      distance_std         : 距離のσ[km]（大＝広域散乱、小＝集中）
      bearing_std          : 方位のσ[度]（大＝全方位、小＝一方向）
      dist_peak_count      : 距離ヒストグラムの局所ピーク数の概算（1=単純、2+=複合）

    Ver13.20追加（Es表面構造の粗い異方性指標、タカユキ指摘）:
      bearing_axial_r          : 0(等方)〜1(単一軸に強く集中)。帯構造は向きはあるが
                                  方向を持たない(北向きの帯=南向きの帯)ため、方位角を
                                  2倍してから平均する軸性統計(axial statistics)で算出。
      bearing_orientation_deg  : 帯の推定方向 (0〜180度、0=南北方向、90=東西方向)
      bearing_ns_alignment     : +1に近いほど南北方向に整合、-1に近いほど東西方向に整合
      bearing_anisotropy_label : "データ不足"/"等方"/"弱い異方性"/"異方性あり(帯状構造の疑い)"
      bearing_hist16           : 方位16分割(22.5度刻み、真北=0起点、時計回り)の
                                  伝播レポート数リスト。ローズダイアグラム表示用。
      あくまで簡易指標であり、実際の雲形状の再構成はできない前提の粗い判定。
    """
    empty = {
        'es_hop_ratio': np.nan, 'dist_in_window_count': np.nan,
        'distance_std': np.nan, 'bearing_std': np.nan, 'dist_peak_count': np.nan,
        'bearing_axial_r': np.nan, 'bearing_orientation_deg': np.nan,
        'bearing_ns_alignment': np.nan, 'bearing_anisotropy_label': 'データ不足',
        'bearing_hist16': [0] * 16,
    }
    if spots_enriched is None or spots_enriched.empty:
        return empty
    d = spots_enriched["distance_km"].dropna()
    b = spots_enriched["bearing_deg"].dropna()
    if len(d) < 1:
        return empty
    in_window = (d >= ES_HOP_DIST_MIN_KM) & (d <= ES_HOP_DIST_MAX_KM)
    es_hop_ratio = float(in_window.mean())
    dist_in_window_count = int(in_window.sum())
    distance_std = float(d.std()) if len(d) >= 2 else np.nan
    bearing_std = np.nan
    if len(b) >= 2:
        rad = np.deg2rad(b.values)
        bearing_std = float(np.degrees(np.sqrt(-2.0 * np.log(np.clip(
            np.abs(np.mean(np.exp(1j * rad))), 1e-9, 1.0)))))
    dist_peak_count = np.nan
    if len(d) >= 4:
        bins = np.arange(0, d.max() + 200, 200)
        counts, _ = np.histogram(d.values, bins=bins)
        peaks = 0
        for i in range(1, len(counts) - 1):
            if counts[i] > counts[i-1] and counts[i] > counts[i+1] and counts[i] >= 1:
                peaks += 1
        dist_peak_count = float(max(1, peaks))

    bearing_axial_r = np.nan
    bearing_orientation_deg = np.nan
    bearing_ns_alignment = np.nan
    bearing_anisotropy_label = 'データ不足'
    bearing_hist16 = [0] * 16
    if len(b) >= 5:
        rad2 = np.deg2rad(2.0 * (b.values % 360.0))
        c_mean = float(np.mean(np.cos(rad2)))
        s_mean = float(np.mean(np.sin(rad2)))
        bearing_axial_r = float(math.hypot(c_mean, s_mean))
        theta2_deg = math.degrees(math.atan2(s_mean, c_mean)) % 360.0
        bearing_orientation_deg = (theta2_deg / 2.0) % 180.0
        bearing_ns_alignment = float(math.cos(math.radians(2.0 * bearing_orientation_deg)))
        if bearing_axial_r < 0.15:
            bearing_anisotropy_label = '等方(明確な向きなし)'
        elif bearing_axial_r < 0.40:
            bearing_anisotropy_label = '弱い異方性'
        else:
            bearing_anisotropy_label = '異方性あり(帯状構造の疑い)'
        idxs = (np.floor((b.values % 360.0) / 22.5).astype(int)) % 16
        for i in idxs:
            bearing_hist16[int(i)] += 1

    return {
        'es_hop_ratio': es_hop_ratio, 'dist_in_window_count': dist_in_window_count,
        'distance_std': distance_std, 'bearing_std': bearing_std,
        'dist_peak_count': dist_peak_count,
        'bearing_axial_r': bearing_axial_r,
        'bearing_orientation_deg': bearing_orientation_deg,
        'bearing_ns_alignment': bearing_ns_alignment,
        'bearing_anisotropy_label': bearing_anisotropy_label,
        'bearing_hist16': bearing_hist16,
    }


_area_reach_history: Dict[str, List[Tuple[datetime, float]]] = {}
_AREA_HIST_MAXLEN = 12

_last_psk_success_dt: Optional[datetime] = None

_ft8_bg_thread: Optional[threading.Thread] = None
_ft8_bg_lock = threading.Lock()
_ft8_bg_latest_spots: pd.DataFrame = pd.DataFrame()
_ft8_bg_latest_ok: bool = False
_ft8_bg_latest_err: str = ''
_ft8_bg_latest_ts: float = 0.0
_ft8_bg_stop_flag: bool = False
FT8_BG_DEFAULT_INTERVAL_SEC = 180
FT8_BG_STALE_SEC = 600


def _ft8_bg_worker(callsign: str, poll_interval_sec: int) -> None:
    """PSKReporter を独立スレッドでポーリングして最新スポットをキャッシュ。
    レート制限に触れないよう poll_interval_sec は最低 120 秒を強制する。

    ESDUCT(整理版・バグ修正): 引数のcallsignはスレッド起動時点の初期値としてのみ
    使い、実際に問い合わせる対象は毎ポーリングごとにget_es_settings()から
    読み直す(未設定/読込失敗時は起動時引数へフォールバック)。これにより、
    ダッシュボードの設定ボタンでコールサインを変更した場合、スレッドを
    再起動しなくても次回ポーリングから新しいコールサインが使われる。
    """
    global _ft8_bg_latest_spots, _ft8_bg_latest_ok, _ft8_bg_latest_err
    global _ft8_bg_latest_ts, _last_psk_success_dt
    poll_interval_sec = max(120, int(poll_interval_sec))
    while not _ft8_bg_stop_flag:
        try:
            try:
                live_cs = str(get_es_settings().get("psk_callsign") or "").strip()
            except Exception:
                live_cs = ""
            effective_callsign = live_cs or callsign
            lookback = poll_interval_sec + 120
            spots, ok, err = fetch_pskreporter_spots(effective_callsign, lookback_sec=lookback)
            with _ft8_bg_lock:
                _ft8_bg_latest_ok = bool(ok)
                _ft8_bg_latest_err = str(err or '')
                _ft8_bg_latest_ts = time.time()
                if ok:
                    _ft8_bg_latest_spots = spots
                    _last_psk_success_dt = now_jst()
        except Exception as e:
            with _ft8_bg_lock:
                _ft8_bg_latest_ok = False
                _ft8_bg_latest_err = f"bg_worker: {str(e)[:60]}"
        for _ in range(poll_interval_sec):
            if _ft8_bg_stop_flag:
                break
            time.sleep(1)


def start_ft8_background_poller(callsign: str,
                                 poll_interval_sec: int = FT8_BG_DEFAULT_INTERVAL_SEC
                                 ) -> None:
    """FT8 バックグラウンドポーラーを開始する（多重起動なし・冪等）。"""
    global _ft8_bg_thread
    if _ft8_bg_thread is not None and _ft8_bg_thread.is_alive():
        return
    _ft8_bg_thread = threading.Thread(
        target=_ft8_bg_worker,
        args=(callsign, poll_interval_sec),
        name="EDFS-FT8-BG",
        daemon=True,
    )
    _ft8_bg_thread.start()


def get_ft8_bg_latest() -> Tuple[pd.DataFrame, bool, str, float]:
    """バックグラウンドスレッドが取得した最新 FT8 スポットを返す。
    戻り値: (spots_df, ok, err, age_sec)
    """
    with _ft8_bg_lock:
        age = time.time() - _ft8_bg_latest_ts if _ft8_bg_latest_ts > 0 else float('inf')
        return _ft8_bg_latest_spots.copy(), _ft8_bg_latest_ok, _ft8_bg_latest_err, age


def is_ft8_bg_active() -> bool:
    return _ft8_bg_thread is not None and _ft8_bg_thread.is_alive()


ES_SETTINGS_FILE = str(DEFAULT_DATA_DIR / "es_settings.json")
_ES_SETTINGS_LOCK = threading.Lock()
_ES_SETTINGS_DEFAULT = {
    "psk_callsign": "",
    "psk_enabled": True,
    "ft8_bg": False,
    "ft8_bg_interval": FT8_BG_DEFAULT_INTERVAL_SEC,
    "setup_done": False,
}
_ES_SETTINGS = dict(_ES_SETTINGS_DEFAULT)
_ES_SETTINGS_LOADED = False
ES_SETTINGS_APPLY_NOW = threading.Event()


def load_es_settings() -> Dict[str, Any]:
    global _ES_SETTINGS, _ES_SETTINGS_LOADED
    with _ES_SETTINGS_LOCK:
        if os.path.exists(ES_SETTINGS_FILE):
            try:
                with open(ES_SETTINGS_FILE, "r", encoding="utf-8") as f:
                    data = json.load(f)
                for k in _ES_SETTINGS_DEFAULT:
                    if k in data:
                        _ES_SETTINGS[k] = data[k]
            except Exception as e:
                print(f"[settings] es_settings.json 読込失敗: {e}")
        _ES_SETTINGS_LOADED = True
        return dict(_ES_SETTINGS)


def save_es_settings() -> None:
    with _ES_SETTINGS_LOCK:
        snap = dict(_ES_SETTINGS)
    try:
        tmp = ES_SETTINGS_FILE + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(snap, f, ensure_ascii=False, indent=1)
        os.replace(tmp, ES_SETTINGS_FILE)
    except Exception as e:
        print(f"[settings] es_settings.json 保存失敗: {e}")


def get_es_settings() -> Dict[str, Any]:
    """現在有効なEs側設定を返す。初回呼び出し時にファイルからロードする。"""
    if not _ES_SETTINGS_LOADED:
        return load_es_settings()
    with _ES_SETTINGS_LOCK:
        return dict(_ES_SETTINGS)


def update_es_settings(patch: Dict[str, Any]) -> Dict[str, Any]:
    """ダッシュボードの設定ボタンから呼ばれる更新関数。
    保存後、ES_SETTINGS_APPLY_NOWをセットして「次回サイクルを待たずに
    できるだけ早く反映したい」という意図を待機ループ側へ伝える。
    """
    applied = {}
    with _ES_SETTINGS_LOCK:
        if "psk_callsign" in patch:
            cs = str(patch["psk_callsign"] or "").strip()[:32]
            _ES_SETTINGS["psk_callsign"] = cs
            applied["psk_callsign"] = cs
        if "psk_enabled" in patch:
            _ES_SETTINGS["psk_enabled"] = bool(patch["psk_enabled"])
            applied["psk_enabled"] = _ES_SETTINGS["psk_enabled"]
        if "ft8_bg" in patch:
            _ES_SETTINGS["ft8_bg"] = bool(patch["ft8_bg"])
            applied["ft8_bg"] = _ES_SETTINGS["ft8_bg"]
        if "ft8_bg_interval" in patch:
            try:
                iv = max(120, int(patch["ft8_bg_interval"]))
                _ES_SETTINGS["ft8_bg_interval"] = iv
                applied["ft8_bg_interval"] = iv
            except (TypeError, ValueError):
                pass
        if patch.get("mark_setup_done"):
            _ES_SETTINGS["setup_done"] = True
            applied["setup_done"] = True
    save_es_settings()
    ES_SETTINGS_APPLY_NOW.set()
    return applied


_self_climatology_cache: Dict[Tuple[str, int, int], Dict[str, float]] = {}
_self_climatology_built_at: float = 0.0


def _hour_bin(hour: int) -> int:
    return (int(hour) // HOUR_BIN_SIZE) * HOUR_BIN_SIZE


def build_self_climatology(log_dir: Path) -> Dict[Tuple[str, int, int], Dict[str, float]]:
    """edfs_area_history.csv(既存の永続化ログ、学習ラベルとは無関係の観測記録)
    から、エリア×月×時間帯(3時間刻み)ごとの到達率の平均値とサンプル数を集計する。
    外部データ取得不要、EDFS自身の蓄積ログのみで完結する「自己気候値」。

    reach_rate列(生の実測到達率、compute_area_reachability()の出力そのもの)を
    使用しており、Rev14.0のIDW補間値(reach_ratio_compensated)は含まれない
    (save_edfs_area_history()がrow['reach_ratio']=生値のみを保存しているため)。
    """
    path = log_dir / EDFS_AREA_HISTORY_FILE
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path, usecols=['timestamp', 'area', 'reach_rate'], parse_dates=['timestamp'])
    except Exception:
        return {}
    if df.empty:
        return {}
    df['reach_rate'] = pd.to_numeric(df['reach_rate'], errors='coerce')
    df = df.dropna(subset=['reach_rate', 'timestamp'])
    if df.empty:
        return {}
    df['month'] = df['timestamp'].dt.month
    df['hour_bin'] = df['timestamp'].dt.hour.apply(_hour_bin)
    grouped = df.groupby(['area', 'month', 'hour_bin'])['reach_rate'].agg(['mean', 'count'])
    result: Dict[Tuple[str, int, int], Dict[str, float]] = {}
    for (area, month, hbin), row in grouped.iterrows():
        result[(str(area), int(month), int(hbin))] = {'mean': float(row['mean']), 'n': int(row['count'])}
    return result


def get_self_climatology(log_dir: Path, refresh_min: int = SELF_CLIMATOLOGY_REFRESH_MIN,
                          force: bool = False) -> Dict[Tuple[str, int, int], Dict[str, float]]:
    """自己ログ季節統計をキャッシュ付きで取得する。起動時と、以後は
    refresh_min分ごとにのみ再集計する(edfs_area_history.csvは日々肥大化するため
    毎サイクル読み直すのは避ける)。
    """
    global _self_climatology_cache, _self_climatology_built_at
    age_min = (time.time() - _self_climatology_built_at) / 60.0
    if force or (not _self_climatology_cache) or age_min >= refresh_min:
        _self_climatology_cache = build_self_climatology(log_dir)
        _self_climatology_built_at = time.time()
    return _self_climatology_cache


_nict_clim_bg_thread: Optional[threading.Thread] = None
_nict_clim_bg_lock = threading.Lock()
_nict_clim_df_cache: pd.DataFrame = pd.DataFrame()
_nict_clim_bg_status: Dict = {'ok': False, 'err': '', 'in_progress': False, 'last_success_ts': 0.0}
_nict_clim_bg_stop_flag: bool = False

_NICT_MANUAL_LINE_RE = re.compile(r'^([A-Za-z0-9]+),(\d{14}):(.*)$')


def fetch_nict_manual_year_raw(station_key: str, year: int,
                                timeout: int = NICT_CLIMATOLOGY_TIMEOUT_SEC) -> Optional[str]:
    """factor-manual-{code}-{year}H.sjis.txt(手動読取1時間値、Shift-JIS)を
    1年分まとめて取得する。取得失敗(404・タイムアウト等)時はNoneを返し、
    呼び出し側はその年を単純にスキップする(致命的エラーにしない)。
    """
    code = NICT_CLIMATOLOGY_STATION_CODES.get(station_key)
    if not code:
        return None
    url = NICT_CLIMATOLOGY_BASE_URL.format(code=code, year=year)
    try:
        headers = {"User-Agent": "EDFS climatology fetch (Rev14.1); low-frequency bulk archive access"}
        r = requests.get(url, timeout=timeout, headers=headers)
        r.raise_for_status()
        return r.content.decode('shift_jis', errors='replace')
    except Exception:
        return None


def normalize_foes(raw_value: str) -> Optional[float]:
    """NICT手動読取1時間値の foEs 生値をMHzへ正規化する。

    Rev14.1修正: 従来は生値をそのままMHzとして扱いNICT_CLIMATOLOGY_FOES_THRESHOLD
    (6.4MHz)と直接比較していたが、生データの単位(MHzか0.1MHz単位か)を検証する層が
    存在しなかった。単位変換をここに一元化する。
    NICT_FOES_RAW_SCALE_DIVISOR(=10.0、0.1MHz単位)は実アーカイブで検証済み
    (定数定義箇所のコメント参照)。
    """
    try:
        v = float(raw_value)
    except (TypeError, ValueError):
        return None
    return v / NICT_FOES_RAW_SCALE_DIVISOR


def parse_nict_manual_foEs(raw_text: str) -> pd.DataFrame:
    """factor-manual-*.sjis.txt形式(手動読取1時間値)からfoEs列を抽出する。

    フォーマット例:
      #  fmin ,foE ,h'E ,foEs ,h'Es ,TYPES,fbEs ,foF1 , ...
      TO536,20250101120000:     ,369  ,---  , 35  ,---  ,     ,---  , ...

    局コード,YYYYMMDDHHMMSS: の後にカンマ区切りで値が続き、4番目(0-indexで3番目)
    がfoEs。欠測は '---' または空白。'L(54)'のようなフラグ付き値が混在する
    可能性があるため、先頭の数値部分のみを正規表現で抽出する。
    """
    rows = []
    numeric_re = re.compile(r'^-?\d+\.?\d*')
    for line in raw_text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = _NICT_MANUAL_LINE_RE.match(line)
        if not m:
            continue
        ts_str, values_str = m.group(2), m.group(3)
        fields = values_str.split(',')
        if len(fields) < 4:
            continue
        foes_raw = fields[3].strip()
        foes_val = np.nan
        if foes_raw and foes_raw != '---':
            m2 = numeric_re.match(foes_raw)
            if m2:
                normalized = normalize_foes(m2.group())
                foes_val = normalized if normalized is not None else np.nan
        try:
            ts = datetime.strptime(ts_str, '%Y%m%d%H%M%S')
        except ValueError:
            continue
        rows.append({'timestamp': ts, 'foEs': foes_val})
    return pd.DataFrame(rows)


def build_nict_station_climatology(station_key: str, years: List[int]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """1局分、複数年のfoEsデータを取得・結合し、月×時間帯(3時間刻み)ごとに
    平均foEsと『Es開通相当』(foEs>=NICT_CLIMATOLOGY_FOES_THRESHOLD)の割合を集計する。

    戻り値: (集計後df(month, hour_bin, mean_foEs, count, open_frac),
             生の時刻別df(timestamp, foEs))
    Rev14.2追加(#再⑤): 生データも合わせて返すようにした。3時間binの集計方式
    自体は変更しないが、将来1時間bin等との比較評価を行う際に再取得が不要になる
    よう、生データをbuild_nict_climatology_all_stations()側で保存する。
    """
    frames = []
    for yr in years:
        raw = fetch_nict_manual_year_raw(station_key, yr)
        if raw is None:
            continue
        df = parse_nict_manual_foEs(raw)
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame(), pd.DataFrame()
    all_df = pd.concat(frames, ignore_index=True)
    all_df['month'] = all_df['timestamp'].dt.month
    all_df['hour_bin'] = all_df['timestamp'].dt.hour.apply(_hour_bin)
    raw_out = all_df[['timestamp', 'foEs']].copy()

    def _agg(g: pd.DataFrame) -> pd.Series:
        valid = g['foEs'].dropna()
        n = len(valid)
        if n == 0:
            return pd.Series({'mean_foEs': np.nan, 'count': 0, 'open_frac': np.nan})
        return pd.Series({
            'mean_foEs': float(valid.mean()),
            'count': int(n),
            'open_frac': float((valid >= NICT_CLIMATOLOGY_FOES_THRESHOLD).sum()) / n,
        })

    try:
        grouped = all_df.groupby(['month', 'hour_bin']).apply(_agg, include_groups=False).reset_index()
    except TypeError:
        grouped = all_df.groupby(['month', 'hour_bin']).apply(_agg).reset_index()
    return grouped, raw_out


def build_nict_climatology_all_stations(years: Optional[List[int]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """観測4局すべてについて季節統計を構築する(③の自動化本体)。
    ネットワーク取得を伴うため、呼び出し側(バックグラウンドスレッド)で
    メインサイクルをブロックしないよう配慮すること。

    戻り値: (集計後df, 生の時刻別df(station列付き))
    Rev14.2追加(#再⑤): 生データも合わせて返す(理由はbuild_nict_station_climatology
    のdocstring参照)。
    """
    if years is None:
        cur_year = now_jst().year
        years = list(range(cur_year - NICT_CLIMATOLOGY_YEARS_BACK + 1, cur_year + 1))
    frames = []
    raw_frames = []
    for station_key in NICT_CLIMATOLOGY_STATION_CODES.keys():
        stat_clim, stat_raw = build_nict_station_climatology(station_key, years)
        if not stat_clim.empty:
            stat_clim['station'] = station_key
            frames.append(stat_clim)
        if not stat_raw.empty:
            stat_raw = stat_raw.copy()
            stat_raw['station'] = station_key
            raw_frames.append(stat_raw)
    if not frames:
        return pd.DataFrame(), pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out['built_at_epoch'] = time.time()
    raw_out = pd.concat(raw_frames, ignore_index=True) if raw_frames else pd.DataFrame()
    return out[['station', 'month', 'hour_bin', 'mean_foEs', 'count', 'open_frac', 'built_at_epoch']], raw_out


def _refresh_nict_climatology_once(log_dir: Path, refresh_days: int) -> None:
    """キャッシュの鮮度を確認し、必要なら再取得する(1回分の処理本体)。"""
    global _nict_clim_df_cache, _nict_clim_bg_status
    path = log_dir / NICT_CLIMATOLOGY_CACHE_FILE
    raw_path = log_dir / NICT_CLIMATOLOGY_RAW_CACHE_FILE
    if path.exists():
        try:
            cached = pd.read_csv(path)
            if not cached.empty and 'built_at_epoch' in cached.columns:
                age_days = (time.time() - float(cached['built_at_epoch'].max())) / 86400.0
                if age_days < refresh_days:
                    with _nict_clim_bg_lock:
                        _nict_clim_df_cache = cached
                        _nict_clim_bg_status['ok'] = True
                    return
        except Exception:
            pass
    with _nict_clim_bg_lock:
        _nict_clim_bg_status['in_progress'] = True
    df_new, raw_new = build_nict_climatology_all_stations()
    with _nict_clim_bg_lock:
        _nict_clim_bg_status['in_progress'] = False
    if not df_new.empty:
        try:
            df_new.to_csv(path, index=False)
        except Exception:
            pass
        if not raw_new.empty:
            try:
                raw_new.to_csv(raw_path, index=False)
            except Exception:
                pass
        with _nict_clim_bg_lock:
            _nict_clim_df_cache = df_new
            _nict_clim_bg_status['ok'] = True
            _nict_clim_bg_status['err'] = ''
            _nict_clim_bg_status['last_success_ts'] = time.time()
    else:
        if path.exists():
            try:
                stale = pd.read_csv(path)
                with _nict_clim_bg_lock:
                    _nict_clim_df_cache = stale
            except Exception:
                pass
        with _nict_clim_bg_lock:
            _nict_clim_bg_status['ok'] = False
            _nict_clim_bg_status['err'] = 'fetch_failed_or_empty'


def _nict_climatology_bg_worker(log_dir: Path, refresh_days: int) -> None:
    while not _nict_clim_bg_stop_flag:
        try:
            _refresh_nict_climatology_once(log_dir, refresh_days)
        except Exception as e:
            with _nict_clim_bg_lock:
                _nict_clim_bg_status['err'] = f"bg_worker: {str(e)[:80]}"
        for _ in range(24 * 3600):
            if _nict_clim_bg_stop_flag:
                return
            time.sleep(1)


def start_nict_climatology_background(log_dir: Path,
                                       refresh_days: int = NICT_CLIMATOLOGY_REFRESH_DAYS) -> None:
    """NICT気候値バックグラウンド取得を開始する(多重起動なし・冪等)。"""
    global _nict_clim_bg_thread
    if _nict_clim_bg_thread is not None and _nict_clim_bg_thread.is_alive():
        return
    _nict_clim_bg_thread = threading.Thread(
        target=_nict_climatology_bg_worker,
        args=(log_dir, refresh_days),
        name="EDFS-NICT-Climatology-BG",
        daemon=True,
    )
    _nict_clim_bg_thread.start()


def get_nict_climatology_df() -> pd.DataFrame:
    with _nict_clim_bg_lock:
        return _nict_clim_df_cache.copy()


def is_nict_climatology_bg_active() -> bool:
    return _nict_clim_bg_thread is not None and _nict_clim_bg_thread.is_alive()


def nearest_nict_station_for_area(area: str) -> Optional[str]:
    """指定エリアの代表座標(AREA_CENTERS)から最も近いNICT観測局を返す。
    Rev14.1修正: 既存の nearest_nict_station(lat, lon) との名称衝突を避けるため改名。"""
    center = AREA_CENTERS.get(area)
    if center is None:
        return None
    alat, alon = center
    best_station, best_dist = None, float('inf')
    for station_key, (slat, slon) in NICT_STATION_LATLON.items():
        d = haversine_km(alat, alon, slat, slon)
        if d < best_dist:
            best_dist = d
            best_station = station_key
    return best_station


def get_climatology_prior(area: str, dt: datetime,
                           self_clim: Dict[Tuple[str, int, int], Dict[str, float]],
                           nict_clim_df: pd.DataFrame) -> Optional[Dict]:
    """②(自己ログ)を優先し、サンプル不足時は③(NICT統計)にフォールバックする
    コールドスタート用の季節気候値取得。

    優先順位:
      1. 自己ログ(同エリア・同月・同時間帯)がSELF_CLIMATOLOGY_MIN_SAMPLES以上
         → 最優先(自身の運用環境に最も忠実な統計のため)
      2. NICT統計(最寄り観測局・同月・同時間帯)がNICT_CLIMATOLOGY_MIN_SOURCES以上
         → 次点(汎用の電離層物理統計)
      3. 自己ログはあるがサンプル不足 → 参考程度として返す(その旨を明示)
      4. いずれも無い → None(呼び出し側は表示自体を省略する)

    戻り値: {'prob_pct': float, 'source': str, 'n': int, 'detail': str,
             'metric': str, 'interpretation': str, 'role': str, 'is_forecast': bool,
             ...(sourceが'nict'の場合のみ 'cache_age_days', 'updated_at', 'stale')}
             または None

    Rev14.1修正:
      - 'prob_pct' はJSON互換のため維持しつつ、意味を明確化する 'metric' /
        'interpretation' を追加した(6.4MHzしきい値の出現率は27MHz開通確率そのもの
        ではなく、あくまでEs状態の気候学的プロキシであるため)。
      - 本関数の戻り値はいずれもRidge予測を置き換えるものではなく、コールドスタート
        時のフォールバック参考値(fallback_prior)である旨を 'role'/'is_forecast' で明示する。
    """
    month = dt.month
    hbin = _hour_bin(dt.hour)
    self_entry = self_clim.get((str(area), month, hbin))

    if self_entry is not None and self_entry.get('n', 0) >= SELF_CLIMATOLOGY_MIN_SAMPLES:
        return {
            'prob_pct': clamp(self_entry['mean'], 0.0, 1.0) * 100.0,
            'source': 'self_log',
            'n': int(self_entry['n']),
            'detail': f"自己ログ実績(同月同時間帯 N={int(self_entry['n'])})",
            'metric': 'self_log_open_rate',
            'interpretation': 'operational_climatology_proxy',
            'role': 'fallback_prior',
            'is_forecast': False,
        }

    station = nearest_nict_station_for_area(area)
    if station is not None and nict_clim_df is not None and not nict_clim_df.empty:
        match = nict_clim_df[(nict_clim_df['station'] == station) &
                              (nict_clim_df['month'] == month) &
                              (nict_clim_df['hour_bin'] == hbin)]
        if not match.empty:
            r = match.iloc[0]
            n_src = int(r.get('count', 0) or 0)
            if not pd.isna(r.get('open_frac')) and n_src >= NICT_CLIMATOLOGY_MIN_SOURCES:
                station_jp = NICT_STATION_NAME_JP.get(station, station)
                built_at = r.get('built_at_epoch')
                cache_age_days = None
                stale = False
                updated_at_iso = None
                if built_at is not None and not pd.isna(built_at):
                    cache_age_days = (time.time() - float(built_at)) / 86400.0
                    stale = cache_age_days > 45
                    updated_at_iso = datetime.fromtimestamp(float(built_at)).isoformat()
                return {
                    'prob_pct': float(r['open_frac']) * 100.0,
                    'source': 'nict',
                    'n': n_src,
                    'detail': f"NICT統計({station_jp}局 同月同時間帯 N={n_src})",
                    'metric': 'foEs_ge_6.4MHz_rate',
                    'interpretation': 'Es_condition_climatology_proxy',
                    'role': 'fallback_prior',
                    'is_forecast': False,
                    'cache_age_days': cache_age_days,
                    'updated_at': updated_at_iso,
                    'stale': stale,
                }

    if self_entry is not None and self_entry.get('n', 0) > 0:
        return {
            'prob_pct': clamp(self_entry['mean'], 0.0, 1.0) * 100.0,
            'source': 'self_log_thin',
            'n': int(self_entry['n']),
            'detail': f"自己ログ実績(サンプル少 N={int(self_entry['n'])})",
            'metric': 'self_log_open_rate',
            'interpretation': 'operational_climatology_proxy_thin_sample',
            'role': 'fallback_prior',
            'is_forecast': False,
        }
    return None



from collections import deque as _deque

_global_cycle_idx: int = 0

_pending_area_samples: List[Dict] = []

_area_forecast_history: Dict[str, List[Tuple[List[float], float, bool]]] = {}

_area_prediction_error: Dict[str, _deque] = {}

_area_ridge_models: Dict[str, Tuple] = {}
_area_holdout_diag: Dict[str, Dict[str, float]] = {}
_area_state_baseline: Dict[str, Dict] = {}

_area_ensemble_shadow_eval: Dict[str, _deque] = {}

_area_total_samples: int = 0
_area_es_samples: int = 0

_es_spatial_history: List[Dict] = []

_model_version: str = "11.0"
_model_score: float = 0.0
_ver113_unlock: bool = False

_ver113_label_history: List[List[int]] = []
_ver113_silhouette: float = 0.0
_ver113_ari: float = 0.0
_ver113_r2_gain: float = 0.0
_ver113_season_start: Optional[datetime] = None

_edfs_level: int = 0
_edfs_level_label: str = "Learning"
_edfs_area_teachers: int = 0
_edfs_centroid_history_count: int = 0
_edfs_lifetime_samples: int = 0
_edfs_predictor_active: Dict[str, bool] = {
    'area': False, 'motion': False, 'lifetime': False, 'dynamics': False
}
_edfs_feature_set: str = "Basic"
_edfs_prev_centroid: Optional[Dict] = None

_edfs_training_mode: str = "normal"
_edfs_training_baseline: Optional[Dict] = None
_edfs_last_drift_check_cycle: int = 0
_edfs_drift_cooldown_until_cycle: int = 0
_edfs_drift_events: List[Dict] = []
_edfs_training_state_loaded: bool = False
_edfs_drift_pending: Dict[str, int] = {'location': 0, 'spread': 0, 'sfi': 0}
_edfs_sfi_recent: List[float] = []

def _build_area_centers() -> Dict[str, Tuple[float, float]]:
    from collections import defaultdict as _dd
    lats: Dict[str, List[float]] = _dd(list)
    lons: Dict[str, List[float]] = _dd(list)
    for _, lat, lon, area in BUILTIN_FIXED_MONITOR:
        lats[area].append(lat); lons[area].append(lon)
    return {a: (float(np.mean(lats[a])), float(np.mean(lons[a]))) for a in lats}

AREA_CENTERS: Dict[str, Tuple[float, float]] = _build_area_centers()


def update_area_reach_history(area_df: pd.DataFrame, obs_dt: datetime) -> None:
    if area_df is None or area_df.empty:
        return
    for _, row in area_df.iterrows():
        area = row["area"]; ratio = row["reach_ratio"]
        if pd.isna(ratio):
            continue
        _area_reach_history.setdefault(area, [])
        _area_reach_history[area].append((obs_dt, float(ratio)))
        if len(_area_reach_history[area]) > _AREA_HIST_MAXLEN:
            _area_reach_history[area] = _area_reach_history[area][-_AREA_HIST_MAXLEN:]



COMPASS_DIRECTIONS_8 = [('N', 0.0), ('NE', 45.0), ('E', 90.0), ('SE', 135.0),
                         ('S', 180.0), ('SW', 225.0), ('W', 270.0), ('NW', 315.0)]


def compute_weighted_ellipse(lats: np.ndarray, lons: np.ndarray, weights: np.ndarray,
                              min_points: int = 3) -> Optional[Dict]:
    """Rev14.5: 点群(緯度,経度,重み)から重み付き2次元共分散楕円を推定する。

    重心周りのローカル接平面(km、東西=x/南北=y)へ経緯度を投影してから
    重み付き共分散行列を固有値分解する。日本全域程度の範囲であれば
    十分な近似(厳密な測地学的正確性は求めない、方位・比の把握が目的)。

    戻り値: {'centroid_lat','centroid_lon','major_km','minor_km',
             'principal_axis_bearing_deg'(0-180、軸なので双方向),
             'anisotropy_ratio'(=major_km/minor_km、>=1、上限999),
             'n_points'} または None(点が足りない/縮退)。

    【重要な注意】この関数はどんな点群にも数値上は適用できるが、
    「点の数が少ない」「点が概ね1本の線上に近い」場合、anisotropy_ratioは
    実際のデータの異方性ではなく点配置の幾何学的形状そのものを強く反映する。
    呼び出し側で、点群の性質に応じてanisotropy_ratioを信頼してよいかを
    判断すること(例: NICT4局は概ねNNE-SSWの弧状配置であり、比の値は
    局配置の形状に支配されがちなので使わない方針。compute_nict_es_gradient()
    のdocstring参照)。
    """
    valid_mask = ~(pd.isna(lats) | pd.isna(lons) | pd.isna(weights))
    lats, lons, weights = lats[valid_mask], lons[valid_mask], weights[valid_mask]
    if len(lats) < min_points:
        return None
    w = np.clip(weights, 1e-9, None)
    clat = float(np.average(lats, weights=w))
    clon = float(np.average(lons, weights=w))
    xs = np.array([haversine_km(clat, clon, clat, lo) * (1.0 if lo >= clon else -1.0) for lo in lons])
    ys = np.array([haversine_km(clat, clon, la, clon) * (1.0 if la >= clat else -1.0) for la in lats])
    mx = np.average(xs, weights=w)
    my = np.average(ys, weights=w)
    cxx = np.average((xs - mx) ** 2, weights=w)
    cyy = np.average((ys - my) ** 2, weights=w)
    cxy = np.average((xs - mx) * (ys - my), weights=w)
    cov = np.array([[cxx, cxy], [cxy, cyy]])
    try:
        eigvals, eigvecs = np.linalg.eigh(cov)
    except Exception:
        return None
    minor_var, major_var = float(eigvals[0]), float(eigvals[1])
    if major_var <= 1e-9:
        return None
    major_vec = eigvecs[:, 1]
    major_km = math.sqrt(max(major_var, 0.0))
    minor_km = math.sqrt(max(minor_var, 0.0))
    bearing = (math.degrees(math.atan2(major_vec[0], major_vec[1])) + 360.0) % 180.0
    anisotropy_ratio = min(major_km / minor_km, 999.0) if minor_km > 1e-6 else 999.0
    return {
        'centroid_lat': clat, 'centroid_lon': clon,
        'major_km': major_km, 'minor_km': minor_km,
        'principal_axis_bearing_deg': bearing,
        'anisotropy_ratio': anisotropy_ratio,
        'n_points': int(len(lats)),
    }


def directional_expectancy_8way(primary_bearing_deg: Optional[float], magnitude: float,
                                 axis_mode: bool = False) -> Dict[str, float]:
    """Rev14.5: 「主方位+強さ」を8方位(N/NE/E/SE/S/SW/W/NW)の期待度(0-1)へ展開する。

    axis_mode=False(ベクトルモード、FT8観測の重心移動方向など、意味のある
    「片側」がある場合): score(dir) = 0.5 + 0.5*magnitude*cos(Δ)
      Δ=0(同方向)で最大、Δ=180(逆方向)で最小になる。

    axis_mode=True(軸モード、楕円の主軸のように0-180で「両端が等価」な場合):
    score(dir) = 0.5 + 0.5*magnitude*cos(2Δ) を使う(0°と180°を同一視する
    ための2倍角)。ただし実際にどちらの端が有利かはこの関数だけでは
    決められないため、呼び出し側でprimary_bearing_degに「有利な側」の
    方位(例: 主軸のうち重心が偏っている側)を渡すこと。

    magnitude(0-1)は「この情報をどれだけ信用するか」も兼ねる。0なら
    全方位0.5(完全にニュートラル=情報なし)になる。
    """
    out = {name: 0.5 for name, _ in COMPASS_DIRECTIONS_8}
    if primary_bearing_deg is None or (isinstance(primary_bearing_deg, float) and pd.isna(primary_bearing_deg)):
        return out
    mag = float(np.clip(magnitude, 0.0, 1.0))
    for name, deg in COMPASS_DIRECTIONS_8:
        diff = math.radians(angular_diff_deg(deg, primary_bearing_deg) or 0.0)
        if axis_mode:
            out[name] = float(np.clip(0.5 + 0.5 * mag * math.cos(2.0 * diff), 0.0, 1.0))
        else:
            out[name] = float(np.clip(0.5 + 0.5 * mag * math.cos(diff), 0.0, 1.0))
    return out


def compute_nict_es_gradient(station_fxes: Dict[str, float]) -> Optional[Dict]:
    """Rev14.5: FT8(PSKReporter)が無い場合でも使える、NICT4局のfxEsから
    Es場の空間勾配(1次元)を推定する。

    【設計上の重要な判断】当初検討した「NICT4局から2次元共分散楕円を
    フィットする」案は採用しなかった。理由は実データ検証済みの構造的問題:
    NICT4局(稚内/国分寺/山川/沖縄)は日本列島に沿ってNNE-SSWの弧状に
    並んでおり、東西方向にはほとんど情報を持たない。4点の重み付け次第で
    共分散楕円の主軸は多少回転する(実測でおよそ26°〜43°の範囲)ものの、
    常に「東西軸」がほぼ選ばれることはなく、かつanisotropy_ratio(長軸/短軸)
    は現在のEs状態ではなく局配置そのものの形状(南北に長い)に支配される
    (どんな重み付けでも5〜10倍程度になる)。これを「Es空間異方性」として
    提示すると、実際には毎回ほぼ同じ値になる指標を「動的な観測」であるかの
    ように見せてしまう。

    そこで、4局が並ぶ軸(稚内-沖縄を結ぶ固定の地理軸)に沿った「1次元の勾配」
    としてのみ扱う。これは、「北ほど強いか南ほど強いか」という情報を素直に
    表現でき、かつ4局の配置から原理的に得られる情報量とも整合する。
    東西方向の情報は原理的に持たないため、後段のdirectional_expectancy_8way()
    では東西寄りの方位は自動的にニュートラル(0.5)寄りになる
    (cos(90°)=0のため)。

    戻り値: {
        'chain_axis_bearing_deg': 局配置の固定軸方位(稚内→沖縄方向、常に一定),
        'favored_bearing_deg': 現在Esが強い側の方位(北寄りが強ければ
                                chain_axis_bearing_degの逆向き、南寄りが
                                強ければchain_axis_bearing_degそのもの),
        'gradient_strength': 0-1。勾配の明瞭さ(南北で偏りが無ければ0に近い)。
        'weighted_position_frac': 0(稚内側に集中)〜1(沖縄側に集中)。
        'n_stations_used': 有効なfxEsを持つ局数(最低2局必要)。
        'source': 'nict_stations', 'role': 'physical_model_prior',
        'is_directly_observed': False,
    } または None(有効局が2未満)。

    Rev14.8(#4、外部データ異常値耐性): NICT側が欠測・パース異常・非現実的な
    値を返す可能性を考慮し、NICT_FXES_PLAUSIBLE_MIN/MAXの範囲外の値は
    「有効なfxEs」として扱わない(station_fxesに含まれていても無視する)。
    セキュリティ攻撃対策ではなく、EDFS自身の内部状態(directional/Evidence)を
    異常値で汚染しないための入力防御。
    """
    stas = ['wakkanai', 'kokubunji', 'yamagawa', 'okinawa']
    vals = []
    for st in stas:
        v = station_fxes.get(st)
        latlon = NICT_STATION_LATLON.get(st)
        if v is None or latlon is None or (isinstance(v, float) and pd.isna(v)):
            continue
        v = float(v)
        if not (math.isfinite(v) and NICT_FXES_PLAUSIBLE_MIN <= v <= NICT_FXES_PLAUSIBLE_MAX):
            continue
        vals.append((st, v, latlon))
    if len(vals) < 2:
        return None

    lat0, lon0 = NICT_STATION_LATLON['wakkanai']
    lat1, lon1 = NICT_STATION_LATLON['okinawa']
    chain_axis_bearing = bearing_deg(lat0, lon0, lat1, lon1)

    dists = np.array([haversine_km(lat0, lon0, la, lo) for _, _, (la, lo) in vals])
    weights = np.clip(np.array([v for _, v, _ in vals]), 1e-6, None)
    weighted_pos = float(np.average(dists, weights=weights))
    pos_frac = float(np.clip(weighted_pos / max(dists.max(), 1e-6), 0.0, 1.0)) if dists.max() > 0 else 0.5

    gradient_strength = float(np.clip(abs(pos_frac - 0.5) * 2.0, 0.0, 1.0))
    favored_bearing = chain_axis_bearing if pos_frac >= 0.5 else ((chain_axis_bearing + 180.0) % 360.0)

    return {
        'chain_axis_bearing_deg': chain_axis_bearing,
        'favored_bearing_deg': favored_bearing,
        'gradient_strength': gradient_strength,
        'weighted_position_frac': pos_frac,
        'n_stations_used': len(vals),
        'source': 'nict_stations',
        'role': 'physical_model_prior',
        'is_directly_observed': False,
    }


NICT_GRADIENT_TREND_WINDOW = 6
NICT_GRADIENT_TREND_MIN_SAMPLES = 4
NICT_GRADIENT_TREND_MIN_DELTA = 0.08
_nict_gradient_history: _deque = _deque(maxlen=NICT_GRADIENT_TREND_WINDOW)


def update_nict_gradient_trend(weighted_position_frac: float) -> Dict:
    """Rev14.6(#2): 今回のweighted_position_fracを履歴へ積み、直近の移動トレンドを
    算出する(グローバル履歴_nict_gradient_historyを更新する副作用あり。
    1サイクルにつき1回だけ呼ぶこと)。

    戻り値: {
        'trend': 'SOUTHWARD'|'NORTHWARD'|'STABLE'|'INSUFFICIENT_DATA',
        'delta': float|None,                    # 後半平均-前半平均(正=南寄りへ推移)
        'favored_moving_bearing_deg': float|None,  # トレンドが指す方位(SOUTHWARD/NORTHWARDのみ)
        'n_samples': int, 'window': int,
    }
    """
    _nict_gradient_history.append((time.time(), float(weighted_position_frac)))
    n = len(_nict_gradient_history)
    if n < NICT_GRADIENT_TREND_MIN_SAMPLES:
        return {'trend': 'INSUFFICIENT_DATA', 'delta': None,
                'favored_moving_bearing_deg': None, 'n_samples': n,
                'window': NICT_GRADIENT_TREND_WINDOW}
    vals = [v for _, v in _nict_gradient_history]
    mid = n // 2
    older_avg = float(np.mean(vals[:mid]))
    recent_avg = float(np.mean(vals[mid:]))
    delta = recent_avg - older_avg
    if abs(delta) < NICT_GRADIENT_TREND_MIN_DELTA:
        return {'trend': 'STABLE', 'delta': delta, 'favored_moving_bearing_deg': None,
                'n_samples': n, 'window': NICT_GRADIENT_TREND_WINDOW}
    lat0, lon0 = NICT_STATION_LATLON['wakkanai']
    lat1, lon1 = NICT_STATION_LATLON['okinawa']
    chain_axis_bearing = bearing_deg(lat0, lon0, lat1, lon1)
    if delta > 0:
        return {'trend': 'SOUTHWARD', 'delta': delta,
                'favored_moving_bearing_deg': chain_axis_bearing,
                'n_samples': n, 'window': NICT_GRADIENT_TREND_WINDOW}
    else:
        return {'trend': 'NORTHWARD', 'delta': delta,
                'favored_moving_bearing_deg': (chain_axis_bearing + 180.0) % 360.0,
                'n_samples': n, 'window': NICT_GRADIENT_TREND_WINDOW}


def combine_directional_evidence(nict_gradient: Optional[Dict],
                                  ft8_anisotropy: Optional[Dict]) -> Dict:
    """Rev14.5: NICTベースの物理モデルprior(常時利用可)とFT8観測(利用可能な
    ときだけ)を統合し、Evidence用の方位場を1つにまとめる。

    命名方針(検討時の結論を反映): FT8が無い場合の結果を「伝播方位」とは
    呼ばない。あくまで「方位別開通期待度(directional_opening_expectancy)」
    という推定値であることを明示する。FT8観測が加わったときだけ、
    より強い確度を持つ「propagation_directional_confidence」として扱う。

    戻り値: {
        'directional_opening_expectancy': {8方位: 0-1},  # 常時出力(NICTのみでも可)
        'propagation_directional_confidence': {8方位: 0-1} または None,  # FT8有りのみ
        'primary_evidence_source': 'ft8_observed' | 'nict_model_only' | 'none',
        'nict_gradient': nict_gradient,      # 元データ(NoneならNICT側情報なし)
        'ft8_anisotropy': ft8_anisotropy,    # 元データ(NoneならFT8側情報なし)
    }
    """
    out = {
        'directional_opening_expectancy': {name: 0.5 for name, _ in COMPASS_DIRECTIONS_8},
        'propagation_directional_confidence': None,
        'primary_evidence_source': 'none',
        'nict_gradient': nict_gradient,
        'ft8_anisotropy': ft8_anisotropy,
    }
    nict_field = None
    if nict_gradient is not None:
        nict_field = directional_expectancy_8way(
            nict_gradient.get('favored_bearing_deg'),
            nict_gradient.get('gradient_strength', 0.0),
            axis_mode=False,
        )
        out['directional_opening_expectancy'] = nict_field
        out['primary_evidence_source'] = 'nict_model_only'

    if ft8_anisotropy is not None and ft8_anisotropy.get('n_points', 0) >= 3:
        ratio = ft8_anisotropy.get('anisotropy_ratio', 1.0)
        ft8_mag = float(np.clip(math.tanh(max(ratio - 1.0, 0.0) / 3.0), 0.0, 1.0))
        ft8_field = directional_expectancy_8way(
            ft8_anisotropy.get('principal_axis_bearing_deg'), ft8_mag, axis_mode=True,
        )
        if nict_field is not None:
            blended = {name: float(np.clip(0.5 * nict_field[name] + 0.5 * ft8_field[name], 0.0, 1.0))
                       for name, _ in COMPASS_DIRECTIONS_8}
        else:
            blended = ft8_field
        out['propagation_directional_confidence'] = blended
        out['directional_opening_expectancy'] = blended
        out['primary_evidence_source'] = 'ft8_observed'

    return out


def compute_es_spatial_features(spots_enriched: pd.DataFrame,
                                  prev_spatial: Optional[Dict] = None,
                                  dt_min: float = 15.0) -> Dict[str, float]:
    """受信スポットの空間分布からEs空間特徴量を算出する。
    es_area_center_lat/lon: Es電波の重心（SNR重み付き）
    es_area_spread        : 到達局の空間分散[km] = 広がりの指標
    es_area_speed         : 前サイクルからの重心移動速度[km/h]
    es_area_direction     : 前サイクルからの重心移動方向[deg]
    es_area_count         : Es窓内(600-2400km)受信局数
    """
    empty = {k: np.nan for k in ['es_area_center_lat','es_area_center_lon',
                                   'es_area_spread','es_area_speed',
                                   'es_area_direction','es_area_count',
                                   'es_major_km','es_minor_km',
                                   'es_principal_axis_bearing_deg','es_anisotropy_ratio']}
    if spots_enriched is None or spots_enriched.empty:
        return empty
    s = spots_enriched.dropna(subset=['recv_lat','recv_lon'])
    if s.empty:
        return empty
    in_win = (s['distance_km'] >= ES_HOP_DIST_MIN_KM) & (s['distance_km'] <= ES_HOP_DIST_MAX_KM)
    es = s[in_win]
    es_count = float(len(es))
    if len(es) == 0:
        return {**empty, 'es_area_count': 0.0}
    w = es['snr'].map(lambda v: sigmoid(float(v)) if not pd.isna(v) else 0.5).values
    w = np.clip(w, 1e-6, None)
    lats = es['recv_lat'].values
    lons = es['recv_lon'].values
    clat = float(np.average(lats, weights=w))
    clon = float(np.average(lons, weights=w))
    dist_from_center = np.array([haversine_km(clat, clon, la, lo)
                                  for la, lo in zip(lats, lons)])
    spread = float(np.sqrt(np.average(dist_from_center**2, weights=w)))
    speed = np.nan; direction = np.nan
    if prev_spatial and not pd.isna(prev_spatial.get('es_area_center_lat', np.nan)):
        plat = prev_spatial['es_area_center_lat']
        plon = prev_spatial['es_area_center_lon']
        dist_moved = haversine_km(plat, plon, clat, clon)
        if not pd.isna(dist_moved) and dt_min > 0:
            speed = dist_moved / (dt_min / 60.0)
            direction = bearing_deg(plat, plon, clat, clon)
    ellipse = compute_weighted_ellipse(lats, lons, w, min_points=3)
    major_km = ellipse['major_km'] if ellipse else np.nan
    minor_km = ellipse['minor_km'] if ellipse else np.nan
    axis_bearing = ellipse['principal_axis_bearing_deg'] if ellipse else np.nan
    aniso_ratio = ellipse['anisotropy_ratio'] if ellipse else np.nan
    return {'es_area_center_lat': clat, 'es_area_center_lon': clon,
            'es_area_spread': spread, 'es_area_speed': speed,
            'es_area_direction': direction, 'es_area_count': es_count,
            'es_major_km': major_km, 'es_minor_km': minor_km,
            'es_principal_axis_bearing_deg': axis_bearing,
            'es_anisotropy_ratio': aniso_ratio}



def _area_feature_vector(area: str, reach_ratio: float,
                          reach_history: List[Tuple[datetime, float]],
                          ft8_feats: Dict, spatial: Dict,
                          fxEs_repr: float,
                          centroid_feats: Optional[Dict] = None) -> List[float]:
    """エリア別予測の説明変数ベクトルを組み立てる（仕様書3章）。

    Rev14.0: centroid_dist_norm/centroid_approach_cos の2次元を追加
    (提案②(b))。観測点が無い/少ないエリアでも、Es重心との幾何関係を
    入力特徴量としてRidgeモデルに与えることで、既存の空間特徴量
    (center_lat/center_lon等、Es重心そのものの絶対座標)だけでは表現
    できていなかった「このエリアから見てEs重心がどれだけ近く、
    どれだけ自分に向かっているか」という相対関係を補う。
    centroid_feats省略時(呼び出し互換用)は0.0埋めとし、旧12次元相当の
    値に2つの0特徴量が単純に追加されるだけの後方互換動作にする。
    """
    rate_change = 0.0; acceleration = 0.0
    hist = reach_history
    if len(hist) >= 2:
        t0 = hist[0][0]
        times = np.array([(t-t0).total_seconds()/60.0 for t,_ in hist])
        vals  = np.array([v for _,v in hist])
        if len(hist) >= 3 and times[-1] > times[0]:
            c = np.polyfit(times, vals, 2)
            rate_change  = float(c[1])
            acceleration = float(2.0 * c[0])
        else:
            dt = max((times[-1]-times[-2]), 1.0)
            rate_change = float((vals[-1]-vals[-2])/dt)
    rpi      = float(ft8_feats.get('rpi_norm', 0.0) or 0.0)
    es_mode  = float(ft8_feats.get('es_mode_score', 0.0) or 0.0)
    fxes_n   = clamp(float(fxEs_repr)/15.0, 0.0, 1.0) if not pd.isna(fxEs_repr) else 0.0
    clat     = float(spatial.get('es_area_center_lat', 0.0) or 0.0)
    clon     = float(spatial.get('es_area_center_lon', 0.0) or 0.0)
    spread   = float(spatial.get('es_area_spread', 0.0) or 0.0) / 2000.0
    speed    = float(spatial.get('es_area_speed', 0.0) or 0.0) / 500.0
    direction= math.sin(math.radians(spatial.get('es_area_direction', 0.0) or 0.0))
    count    = float(spatial.get('es_area_count', 0.0) or 0.0) / 20.0

    centroid_feats = centroid_feats or {}
    c_dist = centroid_feats.get('distance_to_centroid')
    c_bearing_from = centroid_feats.get('bearing_from_centroid')
    centroid_dist_norm = (float(c_dist) / 2000.0) if (c_dist is not None and not pd.isna(c_dist)) else 0.0
    centroid_approach_cos = 0.0
    es_dir_raw = spatial.get('es_area_direction', np.nan)
    if (c_bearing_from is not None and not pd.isna(c_bearing_from) and not pd.isna(es_dir_raw)):
        _diff = angular_diff_deg(es_dir_raw, c_bearing_from)
        if _diff is not None and not pd.isna(_diff):
            centroid_approach_cos = float(math.cos(math.radians(_diff)))

    return [reach_ratio, rate_change, acceleration,
            rpi, es_mode, fxes_n,
            clat/90.0, clon/180.0, spread, speed, direction, count,
            centroid_dist_norm, centroid_approach_cos]

AREA_FEATURE_NAMES = [
    'reach_ratio','rate_change','acceleration',
    'rpi','es_mode_score','fxEs_norm',
    'center_lat','center_lon','spread','speed','direction','es_count',
    'centroid_dist_norm','centroid_approach_cos',
]


def is_es_sample(ft8_feats: Dict, fxEs_repr: float) -> bool:
    """Es有効サンプル判定（仕様書7章）。"""
    rpi    = float(ft8_feats.get('rpi_norm', 0.0) or 0.0)
    mode   = float(ft8_feats.get('es_mode_score', 0.0) or 0.0)
    fxes   = float(fxEs_repr) if not pd.isna(fxEs_repr) else 0.0
    return (rpi >= ES_SAMPLE_THR_RPI or
            mode >= ES_SAMPLE_THR_MODE or
            fxes >= FXES_THRESHOLD)


def register_pending_sample(area_df: pd.DataFrame, ft8_feats: Dict,
                              spatial: Dict, fxEs_repr: float, cycle_idx: int,
                              centroid_feats_by_area: Optional[Dict] = None) -> None:
    """現サイクルの特徴量を保留キューに登録する（仕様書5章）。
    2サイクル後に実績到達率と突合せて学習サンプルに昇格する。

    Rev14.0: 特徴量ベクトルの入力(reach_ratio)には reach_ratio_compensated
    (IDW補間/ブレンド後の値、無ければ生値)を使う。これによりNODATA/LOW
    エリアも「入力特徴量が揃っていれば」pendingサンプルとして登録できる
    ようになる。ただし resolve_pending_samples() 側の学習ラベル(2サイクル
    後の実績)は引き続き生の reach_ratio のみを参照する(そちらは変更して
    いない)ため、学習ラベルとして推定値を使ってしまう循環参照は起きない。
    """
    global _area_total_samples
    if area_df is None or area_df.empty:
        return
    centroid_feats_by_area = centroid_feats_by_area or {}
    areas_snap = {}
    for _, row in area_df.iterrows():
        area = row['area']
        ratio = row.get('reach_ratio_compensated', row['reach_ratio'])
        if pd.isna(ratio):
            continue
        hist = _area_reach_history.get(area, [])
        fvec = _area_feature_vector(area, float(ratio), hist, ft8_feats, spatial, fxEs_repr,
                                    centroid_feats=centroid_feats_by_area.get(area))
        areas_snap[area] = {
            'ratio': float(ratio), 'fvec': fvec,
            'is_compensated': bool(row.get('is_compensated', False)),
        }
    if areas_snap:
        _pending_area_samples.append({
            'cycle_idx': cycle_idx,
            'obs_dt': now_jst().isoformat(),
            'areas': areas_snap,
            'is_es': is_es_sample(ft8_feats, fxEs_repr),
        })
    while len(_pending_area_samples) > 20:
        _pending_area_samples.pop(0)


def resolve_pending_samples(area_df: pd.DataFrame, cycle_idx: int,
                             log_dir: Optional[Path] = None) -> int:
    """現サイクルの実績到達率で2サイクル前の保留を昇格する（仕様書5章）。
    サンプル数はエリア単位で計上する（Ver11.2J修正）。
    「1サイクル=1サンプル」ではなく「1サイクル×有効エリア数=サンプル数」。
    10エリア有効なら1サイクルで10サンプル増加し、100サンプル到達が
    100サイクル(25時間)ではなく10サイクル(2.5時間)程度になる。
    機械学習的にエリア間の観測は独立した教師データとして扱うのが自然。
    戻り値: 今回昇格したエリアサンプル総数

    Rev14.4追加: log_dirを渡すと、昇格した各サンプルの「実際にRidgeへ入力された
    特徴量ベクトルそのもの」をedfs_ridge_training.csvへ永続化する
    (bootstrap_ridge_training_from_csv()参照)。省略時(None)は従来通り
    永続化をスキップする(呼び出し側でlog_dirが未確定なケースへの後方互換)。
    """
    global _area_total_samples, _area_es_samples
    resolved = 0
    current_ratios = {}
    if area_df is not None and not area_df.empty:
        for _, row in area_df.iterrows():
            if not pd.isna(row['reach_ratio']):
                current_ratios[row['area']] = float(row['reach_ratio'])
    still_pending = []
    ridge_training_rows = []
    schema_names_json = json.dumps(AREA_FEATURE_NAMES, ensure_ascii=False)
    for pend in _pending_area_samples:
        age = cycle_idx - pend['cycle_idx']
        if age >= FORECAST_CYCLES:
            valid_area_count = 0
            for area, snap in pend['areas'].items():
                future_ratio = current_ratios.get(area)
                if future_ratio is None:
                    continue
                _area_forecast_history.setdefault(area, [])
                _area_forecast_history[area].append(
                    (snap['fvec'], future_ratio, bool(snap.get('is_compensated', False))))
                if len(_area_forecast_history[area]) > 2000:
                    _area_forecast_history[area] = _area_forecast_history[area][-2000:]
                if log_dir is not None:
                    ridge_training_rows.append({
                        'timestamp': now_jst().strftime('%Y-%m-%d %H:%M:%S'),
                        'area': area,
                        'schema_names': schema_names_json,
                        'feature_json': json.dumps(list(snap['fvec']), ensure_ascii=False),
                        'target': future_ratio,
                        'is_compensated': bool(snap.get('is_compensated', False)),
                    })
                err_q = _area_prediction_error.setdefault(
                    area, _deque(maxlen=BIAS_HISTORY_MAXLEN))
                if area in _area_ridge_models:
                    pred_raw = _predict_area_raw(area, snap['fvec'])
                    if not pd.isna(pred_raw):
                        err_q.append(future_ratio - pred_raw)
                        try:
                            persistence_rpi = snap.get('ratio')
                            ens = compute_shadow_ensemble(area, snap['fvec'], float(pred_raw),
                                                           persistence_rpi, None)
                            pred_ens = ens.get('pred_rpi_ensemble')
                            if pred_ens is not None:
                                eval_q = _area_ensemble_shadow_eval.setdefault(
                                    area, _deque(maxlen=ENSEMBLE_SHADOW_EVAL_MAXLEN))
                                entry = {
                                    'ridge_err': float(future_ratio - pred_raw),
                                    'ensemble_err': float(future_ratio - pred_ens),
                                }
                                if persistence_rpi is not None and not pd.isna(persistence_rpi):
                                    entry['persistence_err'] = float(future_ratio - float(persistence_rpi))
                                eval_q.append(entry)
                        except Exception:
                            pass
                valid_area_count += 1
            _area_total_samples += valid_area_count
            if pend.get('is_es', False):
                _area_es_samples += valid_area_count
            resolved += valid_area_count
        else:
            still_pending.append(pend)
    _pending_area_samples[:] = still_pending
    if ridge_training_rows and log_dir is not None:
        try:
            _csv_append(Path(log_dir) / EDFS_RIDGE_TRAINING_FILE, pd.DataFrame(ridge_training_rows))
        except Exception:
            pass
    return resolved


def _predict_area_raw(area: str, fvec: List[float]) -> float:
    """学習済みRidgeモデルで生の予測値を返す。モデルがなければNaN。"""
    if area not in _area_ridge_models:
        return np.nan
    ridge, scaler, _, _ = _area_ridge_models[area]
    try:
        X = scaler.transform([fvec])
        return float(np.clip(ridge.predict(X)[0], 0.0, 1.0))
    except Exception:
        return np.nan


def predict_area_reach(area: str, fvec: List[float]) -> Tuple[float, float, float]:
    """バイアス補正済みの到達率予測を返す（仕様書9章）。
    戻り値: (pred_corrected, pred_raw, bias)。学習未達時はすべてNaN。
    """
    if _area_total_samples < AREA_PRED_LEARN_START_TOTAL or \
       _area_es_samples < AREA_PRED_LEARN_START_ES:
        return np.nan, np.nan, np.nan
    pred_raw = _predict_area_raw(area, fvec)
    if pd.isna(pred_raw):
        return np.nan, np.nan, np.nan
    bias = 0.0
    if _area_es_samples >= AREA_PRED_LEARN_START_ES:
        err_q = _area_prediction_error.get(area, _deque())
        if len(err_q) >= 5:
            errs = np.array(list(err_q))
            b = float(np.mean(errs))
            std_b = float(np.std(errs))
            if std_b > BIAS_STD_HEALTH_THR / 100.0:
                b = 0.0
            bias = float(np.clip(b, -MAX_BIAS_CORRECTION/100.0, MAX_BIAS_CORRECTION/100.0))
    pred_corr = float(np.clip(pred_raw + bias, 0.0, 1.0))
    return pred_corr, pred_raw, bias


WALK_FORWARD_N_FOLDS = 3
WALK_FORWARD_MIN_TRAIN = 5
WALK_FORWARD_MIN_TEST_TOTAL = 3


def _v13_walk_forward_skill_score(X_raw: np.ndarray, y: np.ndarray,
                                   comp_flags: Optional[np.ndarray] = None,
                                   n_folds: int = WALK_FORWARD_N_FOLDS) -> Dict[str, float]:
    """時系列順のholdout診断。train_area_predictors()の運用中モデル選択
    (r2/ゲーティング)には一切影響しない。
    「インサンプルR2が persistence(直前値)モデルに勝っているか」だけを別枠で可視化する。

    Rev14.2修正(#再②): 従来は80/20の1回分割のみだった。1回だけの分割だと、
    たまたまEsが弱い期間 / 強い期間がテスト区間に当たるかで skill_score が
    大きく変動し、area_pred_confidence()のゲート(AREA_PRED_MIN_HOLDOUT_SKILL)
    の判定自体が不安定になる問題があった。
    ここではrolling-origin(expanding window)のwalk-forwardに変更する:
        train[0:t1] → test[t1:t2]
        train[0:t2] → test[t2:t3]
        ...
    を n_folds 回繰り返し、各foldのOOS残差をプールしてから
    skill_score/holdout_mae/resid_std/idw・observed別Skill/RPI帯別Skillを
    算出する。foldをまたいでOOS残差をプールすることで、IDW使用/非使用の
    ような小さなサブグループでもサンプル数不足(n_sub<3)に陥りにくくなる。
    データが少ない場合は自動的にfold数を減らし、最悪1foldの単純ホールドアウト
    (従来のVer13.20相当)まで縮退する。O(n_folds)回のRidge再学習のみで、
    K-Fold交差検証のような全組合せ評価は行わない（軽量環境向け）。

    Ver13.21由来:
      resid_std   : 全fold OOS残差の標準偏差（予測の確からしさ = ±レンジに使用）
      state_skill : RPI帯(低<0.4 / 高>=0.4、AREA_FEATURE_NAMESの'rpi'列を使用)別の
                    Skill Score。全期間の平均だけでは見えない「どの状態に強いか」の
                    粗い切り口。プール後サンプルが3未満の帯はNoneのまま返す。

    Rev14.1由来(#4/#再④): comp_flags(各サンプルがIDW補間/ブレンド由来の入力
    だったか)を渡すと、全fold OOS区間をIDW使用/非使用で分割してSkillを別々に
    算出する(idw_skill_score/observed_skill_score)。「IDW補間が予測性能を
    本当に改善しているのか」を実データで検証するための参考値であり、
    モデル選択には使わない。

    Rev14.3追加(#①): 'skill_score'(全fold OOSプール、従来通り。
    area_pred_confidence()のゲートはこの値を使い続ける=本番動作は不変)とは別に、
    'global_skill'(='skill_score'のエイリアス)、'recent_skill'(直近foldのみ)、
    'effective_skill'(=SKILL_RECENT_WEIGHT_ALPHA×recent + (1-α)×global)を追加する。
    これらはRev14.3のShadow Ensemble専用の入力であり、本番のarea_pred_confidence
    には使わない(「古い期間は良かったが直近は悪化している」ような状態変化を
    Shadow Ensembleの重み付けにだけ反映させるため)。
    """
    n = len(y)
    diag = {'n_test': 0, 'n_folds': 0, 'method': 'rolling_origin_walk_forward',
            'skill_score': None, 'holdout_mae': None,
            'resid_std': None, 'state_skill': {},
            'idw_skill_score': None, 'observed_skill_score': None,
            'idw_n_test': 0, 'observed_n_test': 0,
            'global_skill': None, 'recent_skill': None, 'effective_skill': None}
    if n < WALK_FORWARD_MIN_TRAIN + WALK_FORWARD_MIN_TEST_TOTAL:
        return diag

    test_pool = n - max(WALK_FORWARD_MIN_TRAIN, int(n * 0.5))
    test_size = max(WALK_FORWARD_MIN_TEST_TOTAL, test_pool // max(1, n_folds))
    folds = max(1, min(n_folds, test_pool // test_size)) if test_size > 0 else 1
    train_end0 = n - folds * test_size
    if train_end0 < WALK_FORWARD_MIN_TRAIN:
        train_end0 = WALK_FORWARD_MIN_TRAIN
        folds = 1
        test_size = n - train_end0
        if test_size < WALK_FORWARD_MIN_TEST_TOTAL:
            return diag

    resid_parts, ytest_parts, ypersist_parts, comp_parts, rpi_parts = [], [], [], [], []
    try:
        rpi_idx = AREA_FEATURE_NAMES.index('rpi')
    except ValueError:
        rpi_idx = None

    for f in range(folds):
        test_start = train_end0 + f * test_size
        test_end = n if f == folds - 1 else test_start + test_size
        if test_start >= n or test_end <= test_start or test_start < WALK_FORWARD_MIN_TRAIN:
            continue
        X_train, y_train = X_raw[:test_start], y[:test_start]
        X_test, y_test = X_raw[test_start:test_end], y[test_start:test_end]
        try:
            scaler_h = StandardScaler()
            Xh = scaler_h.fit_transform(X_train)
            ridge_h = Ridge(alpha=RIDGE_ALPHA)
            ridge_h.fit(Xh, y_train)
            y_pred = ridge_h.predict(scaler_h.transform(X_test))
        except Exception:
            continue
        y_persist = np.concatenate([[y_train[-1]], y_test[:-1]])
        resid_parts.append(y_test - y_pred)
        ytest_parts.append(y_test)
        ypersist_parts.append(y_persist)
        if comp_flags is not None and len(comp_flags) == n:
            comp_parts.append(comp_flags[test_start:test_end])
        if rpi_idx is not None:
            rpi_parts.append(X_test[:, rpi_idx])

    if not resid_parts:
        return diag
    try:
        resid = np.concatenate(resid_parts)
        y_test_all = np.concatenate(ytest_parts)
        y_persist_all = np.concatenate(ypersist_parts)
        mse_model = float(np.mean(resid ** 2))
        mse_persist = float(np.mean((y_test_all - y_persist_all) ** 2))
        diag['n_test'] = int(len(resid))
        diag['n_folds'] = int(len(resid_parts))
        diag['skill_score'] = 1.0 - (mse_model / (mse_persist + 1e-6))
        diag['global_skill'] = diag['skill_score']
        diag['holdout_mae'] = float(np.mean(np.abs(resid)))
        diag['resid_std'] = float(np.std(resid))

        last_resid = resid_parts[-1]
        last_ytest = ytest_parts[-1]
        last_ypersist = ypersist_parts[-1]
        if len(last_resid) >= 3:
            mse_m_recent = float(np.mean(last_resid ** 2))
            mse_p_recent = float(np.mean((last_ytest - last_ypersist) ** 2))
            diag['recent_skill'] = 1.0 - (mse_m_recent / (mse_p_recent + 1e-6))
            diag['effective_skill'] = (SKILL_RECENT_WEIGHT_ALPHA * diag['recent_skill'] +
                                        (1.0 - SKILL_RECENT_WEIGHT_ALPHA) * diag['global_skill'])
        else:
            diag['recent_skill'] = None
            diag['effective_skill'] = diag['global_skill']

        if comp_parts:
            comp_all = np.concatenate(comp_parts)
            for key, mask in (('idw_skill_score', comp_all == True),
                               ('observed_skill_score', comp_all == False)):
                n_sub = int(np.sum(mask))
                n_key = 'idw_n_test' if key == 'idw_skill_score' else 'observed_n_test'
                diag[n_key] = n_sub
                if n_sub >= 3:
                    mse_m = float(np.mean(resid[mask] ** 2))
                    mse_p = float(np.mean((y_test_all[mask] - y_persist_all[mask]) ** 2))
                    diag[key] = 1.0 - (mse_m / (mse_p + 1e-6))

        if rpi_parts:
            rpi_all = np.concatenate(rpi_parts)
            for band_name, mask in (('rpi_low', rpi_all < 0.4), ('rpi_high', rpi_all >= 0.4)):
                if int(np.sum(mask)) >= 3:
                    mse_m = float(np.mean(resid[mask] ** 2))
                    mse_p = float(np.mean((y_test_all[mask] - y_persist_all[mask]) ** 2))
                    diag['state_skill'][band_name] = {
                        'n': int(np.sum(mask)),
                        'skill_score': 1.0 - (mse_m / (mse_p + 1e-6)),
                    }
    except Exception:
        pass
    return diag


def _v13_build_state_baseline(X_raw: np.ndarray, max_sample: int = 200) -> Optional[Dict]:
    """Ver13.21: 現在の状態がモデルの学習経験の中でどれくらい「普通」かを測るための
    ベースラインを構築する。K近傍密度に基づくノベルティ検出の考え方を流用した簡易版:
      1. 特徴量を標準化 (train_area_predictors側のscalerとは別に、生の学習履歴全体で計算)
      2. 学習履歴内の各点について、他の最近傍点までの距離(leave-one-out)を求める
         → これが「学習データが通常どれくらい密集しているか」の分布になる
      3. 新しい状態ベクトルが来たら、学習履歴内の最近傍距離を求め、上記分布の中で
         percentileがどこに位置するかを「経験度」として返す（100%=典型的、0%=最も孤立=未知状態寄り）
    重い計算(自己最近傍探索)はtrain_area_predictors()の再学習タイミングでのみ行う。
    """
    n = len(X_raw)
    if n < 10:
        return None
    mean = X_raw.mean(axis=0)
    std = X_raw.std(axis=0)
    std[std < 1e-9] = 1e-9
    z = (X_raw - mean) / std

    if n > max_sample:
        rng = np.random.default_rng(0)
        sample_idx = rng.choice(n, max_sample, replace=False)
    else:
        sample_idx = np.arange(n)

    self_nn = []
    for i in sample_idx:
        d = np.linalg.norm(z - z[i], axis=1)
        d[i] = np.inf
        self_nn.append(float(np.min(d)))
    self_nn_sorted = np.sort(np.array(self_nn))
    return {'mean': mean, 'std': std, 'z_hist': z, 'self_nn_sorted': self_nn_sorted, 'n': n}


def compute_state_coverage(area: str, fvec: List[float]) -> Optional[float]:
    """Ver13.21: 現在の状態(fvec)がエリアareaの学習履歴の中でどれだけ「経験済み」かを
    0〜100%で返す。100%に近いほど典型的な状態、0%に近いほど過去に類似状態がほぼ無い
    (＝Driftとは別の「未知状態=Novelty」寄り)。ベースライン未構築の場合はNone。
    """
    base = _area_state_baseline.get(area)
    if base is None:
        return None
    try:
        z_cur = (np.array(fvec, dtype=float) - base['mean']) / base['std']
        d = np.linalg.norm(base['z_hist'] - z_cur, axis=1)
        nn_dist = float(np.min(d))
        self_nn_sorted = base['self_nn_sorted']
        rank = float(np.sum(self_nn_sorted >= nn_dist)) / len(self_nn_sorted)
        return float(np.clip(rank * 100.0, 0.0, 100.0))
    except Exception:
        return None


def train_area_predictors() -> Dict[str, float]:
    """全エリアのRidgeモデルを（再）学習する（仕様書3章）。
    戻り値: {area: r2}  ※ 運用中の選択・ゲーティングロジックはVer13.19までと不変。
    """
    if not _SKLEARN_AVAILABLE:
        return {}
    results = {}
    for area, hist in _area_forecast_history.items():
        if len(hist) < 10:
            continue
        X_raw = np.array([h[0] for h in hist])
        y     = np.array([h[1] for h in hist])
        comp_flags = np.array([bool(h[2]) if len(h) > 2 else False for h in hist])
        try:
            scaler = StandardScaler()
            X = scaler.fit_transform(X_raw)
            ridge = Ridge(alpha=RIDGE_ALPHA)
            ridge.fit(X, y)
            y_pred = ridge.predict(X)
            ss_res = float(np.sum((y - y_pred)**2))
            ss_tot = float(np.sum((y - y.mean())**2))
            r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else 0.0
            _area_ridge_models[area] = (ridge, scaler, AREA_FEATURE_NAMES, r2)
            results[area] = r2
            _area_holdout_diag[area] = _v13_walk_forward_skill_score(X_raw, y, comp_flags)
            _base = _v13_build_state_baseline(X_raw)
            if _base is not None:
                _area_state_baseline[area] = _base
        except Exception:
            pass
    return results


def area_pred_confidence(area: str, r2: float) -> str:
    """エリア予測の信頼度ラベル（仕様書10章）。

    Rev14.1修正(#6): 従来はR²とサンプル数のみで判定しており、Holdout Skill
    (_area_holdout_diag、persistenceベースラインとの相対性能)は参考診断として
    別辞書に保存されるのみで採用判定には反映されていなかった。その結果、
    見かけ上のR²が高くても実運用ではpersistence(現状維持予測)にすら勝てない
    モデルがHIGH/MEDIUM/LOW扱いされ得る問題があった。
    ここではholdout_skillが計算済みでかつAREA_PRED_MIN_HOLDOUT_SKILL未満
    (persistenceをほとんど、または全く改善していない)の場合、R²やサンプル数の
    条件を満たしていてもLEARNINGへ据え置く。
    holdout_skillが未計算(None、サンプル不足等)の場合は、従来通りR²/サンプル数
    条件のみで判定する(学習初期にholdout評価が間に合わないケースを塞がないため)。

    Rev14.1再指摘(#1): 当初はholdout_skill<=0のみをゲートしていたが、
    「+0.01でもHIGH判定になり得る」という指摘を受け、実用的な下限
    AREA_PRED_MIN_HOLDOUT_SKILLを新設した。
    """
    if (_area_total_samples < AREA_PRED_LEARN_START_TOTAL or
            _area_es_samples < AREA_PRED_LEARN_START_ES):
        return 'LEARNING'
    holdout_skill = (_area_holdout_diag.get(area, {}) or {}).get('skill_score')
    if holdout_skill is not None and holdout_skill < AREA_PRED_MIN_HOLDOUT_SKILL:
        return 'LEARNING'
    if (_area_total_samples >= AREA_PRED_HIGH_TOTAL and
            _area_es_samples >= AREA_PRED_HIGH_ES and r2 >= AREA_PRED_HIGH_R2):
        return 'HIGH'
    if (_area_total_samples >= AREA_PRED_MED_TOTAL and
            _area_es_samples >= AREA_PRED_MED_ES and r2 >= AREA_PRED_MED_R2):
        return 'MEDIUM'
    if (_area_total_samples >= AREA_PRED_LOW_TOTAL and
            _area_es_samples >= AREA_PRED_LOW_ES and r2 >= AREA_PRED_LOW_R2):
        return 'LOW'
    return 'LEARNING'



def compute_ensemble_model_weight(area: str, fvec: Optional[List[float]]) -> Dict[str, Optional[float]]:
    """Rev14.3(#②): 「Skill(モデル性能)」と「Coverage(適用可能性)」を、
    それぞれ別々に0-1の重みへ変換してから掛け合わせる。

    設計上の判断: skill×coverageを生の値のまま乗算しない(指摘②で名指しされた
    問題点)。理由は、skill_score(persistence比のMSE改善率)とcoverage_pct(0-100の
    百分位)はスケールも意味も異なる量であり、生のまま掛けても数値としての意味が
    保証されないため。ここでは:
      skill_weight    = clip(effective_skill / SKILL_REFERENCE, 0.0, 1.0)
      coverage_weight = clip(state_coverage_pct / 100.0, 0.0, 1.0)
      model_weight     = skill_weight * coverage_weight
    という2段階変換にしてから乗算することで、「モデル性能として十分信頼できるか」
    と「今回の状態がその性能を発揮できる状況か」を別々に評価してから統合する。

    戻り値: {'skill_weight','coverage_weight','model_weight',
             'effective_skill','state_coverage_pct'} (状態不明な項目はNone)
    """
    out = {'skill_weight': None, 'coverage_weight': None, 'model_weight': None,
           'effective_skill': None, 'state_coverage_pct': None}
    diag = _area_holdout_diag.get(area, {}) or {}
    effective_skill = diag.get('effective_skill')
    out['effective_skill'] = effective_skill
    if effective_skill is None:
        skill_weight = 0.0
    else:
        skill_weight = float(np.clip(effective_skill / SKILL_REFERENCE, 0.0, 1.0))
    out['skill_weight'] = skill_weight

    coverage_pct = compute_state_coverage(area, fvec) if fvec is not None else None
    out['state_coverage_pct'] = coverage_pct
    if coverage_pct is None:
        coverage_weight = 0.0
    else:
        coverage_weight = float(np.clip(coverage_pct / 100.0, 0.0, 1.0))
    out['coverage_weight'] = coverage_weight

    out['model_weight'] = skill_weight * coverage_weight
    return out


def geometry_ridge_consistency(eta_min: Optional[float], ridge_pred_rpi: Optional[float],
                                forecast_horizon_min: float = FORECAST_HORIZON_MIN_APPROX,
                                rpi_threshold: float = GEOMETRY_CONSISTENCY_RPI_THRESHOLD) -> str:
    """Rev14.3(#③): 幾何ETA(独立した物理的証拠)とRidge予測(統計モデル)が
    同じ結論を指しているかを判定する。

    設計上の判断: Geometryは「30分後72%」のようなRidgeと同一スケールの数値
    予測を持たない(ETAは到達時刻の見積もりであり、到達時の強度は表さない)ため、
    まず ETA を forecast_horizon_min との大小関係で3値化し(NO_ETA/
    WITHIN_HORIZON/BEYOND_HORIZON)、それとRidgeの「開通相当かどうか」
    (rpi_threshold基準)を組み合わせて判定する。

      WITHIN_HORIZON(ホライズン内に到達見込) + Ridge>=閾値(開通予測)  → CONSISTENT
      WITHIN_HORIZON                        + Ridge<閾値(弱い予測)    → CONFLICT
      BEYOND_HORIZON(ホライズン内は未到達見込) + Ridge<閾値            → CONSISTENT
      BEYOND_HORIZON                        + Ridge>=閾値            → CONFLICT
      NO_ETA(速度ベクトル不定などでETA算出不可) または Ridge予測が無い(LEARNING等)
                                                                      → NEUTRAL
    戻り値: 'CONSISTENT' | 'CONFLICT' | 'NEUTRAL'
    """
    if ridge_pred_rpi is None or (isinstance(ridge_pred_rpi, float) and pd.isna(ridge_pred_rpi)):
        return 'NEUTRAL'
    if eta_min is None or (isinstance(eta_min, float) and pd.isna(eta_min)):
        return 'NEUTRAL'
    ridge_says_open = ridge_pred_rpi >= rpi_threshold
    eta_within_horizon = eta_min <= forecast_horizon_min
    if eta_within_horizon and ridge_says_open:
        return 'CONSISTENT'
    if eta_within_horizon and not ridge_says_open:
        return 'CONFLICT'
    if (not eta_within_horizon) and (not ridge_says_open):
        return 'CONSISTENT'
    return 'CONFLICT'


def climatology_anomaly_pt(ridge_pred_rpi: Optional[float],
                            climatology: Optional[Dict]) -> Optional[float]:
    """Rev14.3(#④): 「今のRidge予測は、その季節・時間帯として普通なのか異常なのか」
    をポイント差(%pt)で表す。 anomaly = ridge_pred_rpi*100 - climatology['prob_pct']

    正: 気候学的な平年値より強い予測(その時間帯としては珍しく良い状態)
    負: 平年値より弱い予測(通常なら有利な時間帯なのにモデルは弱いと判断)
    Ridge予測が無い(LEARNING等)、またはclimatologyが無い場合はNone。
    Ensembleの重みには使わない(独立した異常検知的な参考情報として提示するのみ)。
    """
    if ridge_pred_rpi is None or (isinstance(ridge_pred_rpi, float) and pd.isna(ridge_pred_rpi)):
        return None
    if not climatology:
        return None
    prob_pct = climatology.get('prob_pct')
    if prob_pct is None:
        return None
    return float(ridge_pred_rpi) * 100.0 - float(prob_pct)


def compute_shadow_ensemble(area: str, fvec: Optional[List[float]],
                             ridge_pred_rpi: Optional[float], persistence_rpi: Optional[float],
                             eta_min: Optional[float],
                             geometry_consistency: Optional[str] = None) -> Dict:
    """Rev14.3(#⑤): pred_rpi(本番)は変更せず、その横にpred_rpi_ensembleを算出する。

    model_weight = compute_ensemble_model_weight() の skill_weight×coverage_weight。
    さらにgeometry_ridge_consistency()がCONFLICT(独立した物理的証拠と食い違う)の
    場合、GEOMETRY_CONFLICT_WEIGHT_DAMPINGを掛けてRidge重みを追加減衰させる。

    設計上の判断: レビューで提案された「ridge/persistence/geometryの3方式を
    独立に数値ブレンドする」案(weights={'ridge':0.55,'persistence':0.25,
    'geometry':0.20}のような形)は今回採用しなかった。理由は、Geometryが
    RPIと同一スケールの数値予測を持たず(#③参照)、無理に数値化すると
    その変換方法自体が検証されていない新たな仮定を持ち込むことになるため。
    代わりに、Geometryは「Ridge予測の妥当性チェック」としてRidge重みを
    調整する形で間接的にensembleへ統合する(ridge/persistenceの2方式ブレンド＋
    geometryによる重み減衰)。

    戻り値: {'pred_rpi_ensemble','ensemble_mode','weights':{'ridge','persistence'},
             'geometry_adjustment'}
    """
    out = {'pred_rpi_ensemble': None, 'ensemble_mode': 'SHADOW',
           'weights': {'ridge': None, 'persistence': None}, 'geometry_adjustment': None}
    if ridge_pred_rpi is None or (isinstance(ridge_pred_rpi, float) and pd.isna(ridge_pred_rpi)):
        return out
    if persistence_rpi is None or (isinstance(persistence_rpi, float) and pd.isna(persistence_rpi)):
        return out

    mw = compute_ensemble_model_weight(area, fvec)
    model_weight = mw['model_weight'] if mw['model_weight'] is not None else 0.0

    consistency = geometry_consistency
    if consistency is None:
        consistency = geometry_ridge_consistency(eta_min, ridge_pred_rpi)
    geometry_adjustment = 1.0
    if consistency == 'CONFLICT':
        geometry_adjustment = GEOMETRY_CONFLICT_WEIGHT_DAMPING
    model_weight_final = float(np.clip(model_weight * geometry_adjustment, 0.0, 1.0))

    out['pred_rpi_ensemble'] = model_weight_final * float(ridge_pred_rpi) + \
        (1.0 - model_weight_final) * float(persistence_rpi)
    out['weights'] = {'ridge': model_weight_final, 'persistence': 1.0 - model_weight_final}
    out['geometry_adjustment'] = geometry_adjustment
    return out


def summarize_ensemble_shadow_eval(area: str) -> Dict:
    """Rev14.3(#⑤): resolve_pending_samples()が蓄積したShadow Ensembleの事後誤差
    (_area_ensemble_shadow_eval)を、Evidence出力用に要約する。

    「ensembleを入れる価値があるか」を確認するための唯一の判断材料。
    n(サンプル数)が小さいうちは参考程度にしかならない点に注意。
    戻り値: {'n','ridge_mae','ensemble_mae','ensemble_better_frac'}
    (未蓄積時は全てNone/0)
    """
    dq = _area_ensemble_shadow_eval.get(area)
    out = {'n': 0, 'ridge_mae': None, 'ensemble_mae': None, 'ensemble_better_frac': None,
           'persistence_mae': None, 'ensemble_vs_persistence_better_frac': None,
           'ridge_vs_persistence_better_frac': None,
           'comparison_axes': ['persistence', 'ridge', 'ensemble'],
           'geometric_eta_mae': None,
           'geometric_eta_note': 'Rev14.10時点で幾何ETAはpending登録時にeta_minを保持して'
                                  'いないため4者比較の対象外(未実装、今後の拡張課題)'}
    if not dq:
        return out
    ridge_errs = np.array([d['ridge_err'] for d in dq])
    ens_errs = np.array([d['ensemble_err'] for d in dq])
    out['n'] = int(len(dq))
    out['ridge_mae'] = float(np.mean(np.abs(ridge_errs)))
    out['ensemble_mae'] = float(np.mean(np.abs(ens_errs)))
    out['ensemble_better_frac'] = float(np.mean(np.abs(ens_errs) < np.abs(ridge_errs)))
    persist_dq = [d for d in dq if 'persistence_err' in d]
    if persist_dq:
        persist_errs = np.array([d['persistence_err'] for d in persist_dq])
        ridge_errs_p = np.array([d['ridge_err'] for d in persist_dq])
        ens_errs_p = np.array([d['ensemble_err'] for d in persist_dq])
        out['persistence_mae'] = float(np.mean(np.abs(persist_errs)))
        out['ridge_vs_persistence_better_frac'] = float(
            np.mean(np.abs(ridge_errs_p) < np.abs(persist_errs)))
        out['ensemble_vs_persistence_better_frac'] = float(
            np.mean(np.abs(ens_errs_p) < np.abs(persist_errs)))
    return out


def _v13_pooled_r2(y_true: np.ndarray, y_pred: np.ndarray) -> Optional[float]:
    """予測値と実測値からプールされた決定係数(R²)を計算する。
    Rev14.10(#重大な問題①): Ver11.3のr2_gain本格実装で、baselineとクラスタ別の
    両方に共通して使う。サンプルが少なすぎる、またはy_trueの分散が0(定数)の
    場合は評価不能としてNoneを返す。
    """
    if y_true is None or len(y_true) < 2:
        return None
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    if ss_tot <= 1e-9:
        return None
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    return 1.0 - ss_res / ss_tot


VER113_CLUSTER_MIN_TRAIN = 8
VER113_CLUSTER_MIN_TEST  = 2


def _v13_cluster_r2_gain(scaler_es, km4,
                          n_folds: int = WALK_FORWARD_N_FOLDS) -> Tuple[Optional[float], Dict]:
    """Ver11.3 r2_gain の本格実装(Rev14.10 #重大な問題①対応)。

    【修正前の問題】
      旧実装は `r2_clustered_list.append(min(sil * 0.8, 0.95))` という、
      silhouette係数から線形に作った代理値をクラスタ別R²として扱っており、
      クラスタ別モデルを実際には一度も学習・評価していなかった。これでは
      「クラスタリングによって予測精度が改善した」ことの証明にならない。

    【本実装】
      仕様どおり、以下の手順で実測のR²改善量を算出する:
        全データ   → baselineモデル(1本)      → rolling-origin walk-forward R²
        クラスタ別 → 各クラスタで個別にRidge再学習 → *同じtest期間*で予測
                     → クラスタ別R² → サンプル数加重平均
        r2_gain = weighted(クラスタ別R²) - baseline R²

    【データソースと「同じtest期間」の再現】
      _area_forecast_history(全エリア分)をプールして使う。各サンプルのfvec
      (AREA_FEATURE_NAMES)には既にEs空間状態(center_lat/center_lon/spread/speed、
      正規化済み)が列として含まれているため、正規化を解いて評価元
      (evaluate_ver113_unlock)で既にfitされているscaler_es/km4へ通すだけで、
      「そのサンプルが観測された時点のEsクラスタ」を復元できる。これにより
      _es_spatial_history(サイクル単位・全エリア共通)と_area_forecast_history
      (エリア単位)をタイムスタンプで突合せる必要がない。
      各エリアの学習履歴はリスト内では追加順=時系列順だが、複数エリアを
      跨いだ厳密なタイムスタンプは保持していないため、エリア内インデックスを
      疑似時系列キーとして全エリア分を1本の擬似時系列にプールする
      (resolve_pending_samples()が全エリアを同一サイクルでまとめて昇格させる
      運用のため、この近似は妥当)。fold境界(train_end/test_end)はbaseline・
      クラスタ別の両方でこの同じ擬似時系列上の同じ位置を使うため、
      「同じtest期間で予測」という要件を満たす。

      戻り値: (r2_gain または データ不足/評価不能時はNone, 診断dict)
    """
    diag = {'baseline_r2': None, 'clustered_r2': None, 'n_folds_used': 0,
            'cluster_r2': {}, 'cluster_n_test': {}, 'pooled_n': 0,
            'method': 'per_cluster_walk_forward_r2'}
    if not _area_forecast_history:
        return None, diag

    pooled = []
    for area, hist in _area_forecast_history.items():
        for i, h in enumerate(hist):
            fvec = h[0]
            if len(fvec) != len(AREA_FEATURE_NAMES):
                continue
            pooled.append((i, area, fvec, h[1]))
    if len(pooled) < (WALK_FORWARD_MIN_TRAIN + WALK_FORWARD_MIN_TEST_TOTAL):
        return None, diag
    pooled.sort(key=lambda t: (t[0], t[1]))
    X_all = np.array([p[2] for p in pooled], dtype=float)
    y_all = np.array([p[3] for p in pooled], dtype=float)
    n = len(y_all)
    diag['pooled_n'] = n

    try:
        lat_idx = AREA_FEATURE_NAMES.index('center_lat')
        lon_idx = AREA_FEATURE_NAMES.index('center_lon')
        spread_idx = AREA_FEATURE_NAMES.index('spread')
        speed_idx = AREA_FEATURE_NAMES.index('speed')
        raw_es_state = np.column_stack([
            X_all[:, lat_idx] * 90.0,
            X_all[:, lon_idx] * 180.0,
            X_all[:, spread_idx] * 2000.0,
            X_all[:, speed_idx] * 500.0,
        ])
        Xc_all = scaler_es.transform(raw_es_state)
        labels_all = km4.predict(Xc_all)
    except Exception:
        return None, diag

    test_pool = n - max(WALK_FORWARD_MIN_TRAIN, int(n * 0.5))
    test_size = max(WALK_FORWARD_MIN_TEST_TOTAL, test_pool // max(1, n_folds))
    folds = max(1, min(n_folds, test_pool // test_size)) if test_size > 0 else 1
    train_end0 = n - folds * test_size
    if train_end0 < WALK_FORWARD_MIN_TRAIN:
        train_end0 = WALK_FORWARD_MIN_TRAIN
        folds = 1
        test_size = n - train_end0
        if test_size < WALK_FORWARD_MIN_TEST_TOTAL:
            return None, diag

    baseline_ytrue: List[float] = []
    baseline_ypred: List[float] = []
    cluster_ytrue: Dict[int, List[float]] = {}
    cluster_ypred: Dict[int, List[float]] = {}

    for f in range(folds):
        test_start = train_end0 + f * test_size
        test_end = n if f == folds - 1 else test_start + test_size
        if test_start >= n or test_end <= test_start or test_start < WALK_FORWARD_MIN_TRAIN:
            continue
        X_train_all, y_train_all = X_all[:test_start], y_all[:test_start]
        X_test_all, y_test_all = X_all[test_start:test_end], y_all[test_start:test_end]

        try:
            scaler_b = StandardScaler()
            Xb = scaler_b.fit_transform(X_train_all)
            ridge_b = Ridge(alpha=RIDGE_ALPHA)
            ridge_b.fit(Xb, y_train_all)
            y_pred_b = ridge_b.predict(scaler_b.transform(X_test_all))
            baseline_ytrue.extend(y_test_all.tolist())
            baseline_ypred.extend(y_pred_b.tolist())
        except Exception:
            pass

        labels_train = labels_all[:test_start]
        labels_test = labels_all[test_start:test_end]
        for cl in set(labels_all):
            tr_mask = labels_train == cl
            te_mask = labels_test == cl
            if (int(np.sum(tr_mask)) < VER113_CLUSTER_MIN_TRAIN or
                    int(np.sum(te_mask)) < VER113_CLUSTER_MIN_TEST):
                continue
            try:
                Xc_train = X_train_all[tr_mask]
                yc_train = y_train_all[tr_mask]
                Xc_test = X_test_all[te_mask]
                yc_test = y_test_all[te_mask]
                scaler_c = StandardScaler()
                Xcs = scaler_c.fit_transform(Xc_train)
                ridge_c = Ridge(alpha=RIDGE_ALPHA)
                ridge_c.fit(Xcs, yc_train)
                y_pred_c = ridge_c.predict(scaler_c.transform(Xc_test))
                cluster_ytrue.setdefault(int(cl), []).extend(yc_test.tolist())
                cluster_ypred.setdefault(int(cl), []).extend(y_pred_c.tolist())
            except Exception:
                continue
        diag['n_folds_used'] += 1

    if not baseline_ytrue:
        return None, diag
    r2_baseline = _v13_pooled_r2(np.array(baseline_ytrue), np.array(baseline_ypred))
    diag['baseline_r2'] = r2_baseline

    weighted_num, weighted_den = 0.0, 0.0
    for cl, yt in cluster_ytrue.items():
        yp = cluster_ypred[cl]
        r2_c = _v13_pooled_r2(np.array(yt), np.array(yp))
        diag['cluster_r2'][cl] = r2_c
        diag['cluster_n_test'][cl] = len(yt)
        if r2_c is not None:
            weighted_num += r2_c * len(yt)
            weighted_den += len(yt)
    r2_clustered = (weighted_num / weighted_den) if weighted_den > 0 else None
    diag['clustered_r2'] = r2_clustered

    if r2_baseline is None or r2_clustered is None:
        return None, diag
    return float(r2_clustered - r2_baseline), diag



def evaluate_ver113_unlock() -> Tuple[bool, Dict]:
    """Ver11.3解禁条件を評価する。
    戻り値: (unlock:bool, metrics:dict)
    データ不足時はすべてFalseを返す（unlock=False）。
    """
    global _ver113_silhouette, _ver113_ari, _ver113_r2_gain, _ver113_season_start
    metrics = {
        'sample_count': _area_total_samples,
        'es_samples': _area_es_samples,
        'silhouette': np.nan,
        'r2_gain': np.nan,
        'ari': np.nan,
        'season_days': 0.0,
    }
    if _area_total_samples < VER113_SAMPLE_MIN or _area_es_samples < VER113_ES_MIN:
        return False, metrics
    if not _SKLEARN_AVAILABLE or not _SKLEARN_CLUSTER_OK:
        return False, metrics
    es_rows = [r for r in _es_spatial_history
               if not pd.isna(r.get('es_area_center_lat', np.nan))]
    if len(es_rows) < 50:
        return False, metrics
    feat_mat = np.array([[
        r.get('es_area_center_lat', 0.0),
        r.get('es_area_center_lon', 0.0),
        r.get('es_area_spread', 0.0),
        r.get('es_area_speed', 0.0) if not pd.isna(r.get('es_area_speed', np.nan)) else 0.0,
    ] for r in es_rows])
    mask = np.all(np.isfinite(feat_mat), axis=1)
    feat_clean = feat_mat[mask]
    if len(feat_clean) < 30:
        return False, metrics
    try:
        from sklearn.preprocessing import StandardScaler as _SS
        scaler_es = _SS()
        Xc = scaler_es.fit_transform(feat_clean)
        km4 = KMeans(n_clusters=4, n_init=10, random_state=42).fit(Xc)
        labels4 = km4.labels_
        sil = float(_sk_silhouette(Xc, labels4)) if len(set(labels4)) > 1 else 0.0
        _ver113_silhouette = sil
        metrics['silhouette'] = sil
        if sil < VER113_SILHOUETTE_MIN:
            return False, metrics
        r2_gain, r2_diag = _v13_cluster_r2_gain(scaler_es, km4)
        metrics['r2_gain_diag'] = r2_diag
        if r2_gain is None:
            return False, metrics
        _ver113_r2_gain = r2_gain
        metrics['r2_gain'] = r2_gain
        if r2_gain < VER113_R2GAIN_MIN:
            return False, metrics
        ari = np.nan
        if len(_ver113_label_history) >= 2:
            try:
                ari = float(_sk_ari(_ver113_label_history[-2], _ver113_label_history[-1]))
            except Exception:
                ari = 0.0
        _ver113_ari = ari if not pd.isna(ari) else 0.0
        metrics['ari'] = _ver113_ari
        _ver113_label_history.append(list(labels4[:100]))
        if len(_ver113_label_history) > 10:
            _ver113_label_history.pop(0)
        if pd.isna(ari) or ari < VER113_ARI_MIN:
            return False, metrics
        if _ver113_season_start is None:
            _ver113_season_start = now_jst()
        season_days = (now_jst() - _ver113_season_start).days
        metrics['season_days'] = float(season_days)
        if season_days < VER113_SEASON_DAYS:
            return False, metrics
        return True, metrics
    except Exception:
        return False, metrics



def update_model_version() -> Tuple[str, float, bool]:
    """モデルバージョンとスコアを算出して更新する（毎サイクル呼び出し）。
    戻り値: (model_version, model_score, ver113_active)
    """
    global _model_version, _model_score, _ver113_unlock
    unlock, metrics = evaluate_ver113_unlock()
    _ver113_unlock = unlock
    sp = min(_area_total_samples / VER113_SAMPLE_MIN, 1.0)
    ep = min(_area_es_samples   / VER113_ES_MIN,      1.0)
    si = float(np.nan_to_num(float(_ver113_silhouette)))
    ar = float(np.nan_to_num(float(_ver113_ari)))
    rg = min(float(np.nan_to_num(float(_ver113_r2_gain))) / VER113_R2GAIN_MIN, 1.0)
    score = (MODEL_SCORE_W_SAMPLE * sp +
             MODEL_SCORE_W_ES     * ep +
             MODEL_SCORE_W_SIL    * si +
             MODEL_SCORE_W_ARI    * ar +
             MODEL_SCORE_W_R2GAIN * rg) * 100.0
    _model_score = clamp(score, 0.0, 100.0)
    if unlock:
        mv = "11.3"
    elif (_area_total_samples >= AREA_PRED_LEARN_START_TOTAL and
          _area_es_samples >= AREA_PRED_LEARN_START_ES):
        mv = "11.2"
    elif _area_total_samples >= 30:
        mv = "11.1"
    else:
        mv = "11.0"
    _model_version = mv
    return mv, _model_score, unlock


def area_pred_summary_for_display(area: str, current_ratio: float,
                                   fvec: List[float]) -> Dict:
    """エリア予測の表示用サマリを返す（仕様書8章の各表示モードに対応）。"""
    r2 = _area_ridge_models.get(area, (None,None,None,0.0))[3] if area in _area_ridge_models else 0.0
    conf = area_pred_confidence(area, r2)
    if conf == 'LEARNING':
        return {'conf': 'LEARNING', 'pred': np.nan, 'bias': np.nan,
                'total': _area_total_samples, 'es': _area_es_samples}
    pred_corr, pred_raw, bias = predict_area_reach(area, fvec)
    return {'conf': conf, 'pred': pred_corr, 'bias': bias,
            'total': _area_total_samples, 'es': _area_es_samples}


def ver113_es_type_label() -> str:
    """Ver11.3有効時のEsタイプラベル（停滞型/移動型/拡大型/多層型）。"""
    if not _ver113_unlock or not _es_spatial_history:
        return '不明'
    last = _es_spatial_history[-1]
    speed = last.get('es_area_speed', np.nan)
    spread = last.get('es_area_spread', np.nan)
    count = last.get('es_area_count', 0.0)
    if pd.isna(speed):
        return '不明'
    if speed < 50:
        return '停滞型'
    if not pd.isna(spread) and spread > 800:
        return '拡大型' if count > 6 else '移動型'
    if not pd.isna(count) and count > 8:
        return '多層型'
    return '移動型'


    if _area_total_samples >= VER113_SAMPLE_MIN - 10:
        return f"Es分類学習 特徴量:{_area_total_samples}/{VER113_SAMPLE_MIN} 分類評価中"
    return f"Es分類学習 特徴量:{_area_total_samples}/{VER113_SAMPLE_MIN} Ver11.3準備中"



def _acquire_file_lock(lock_path: Path, timeout: float = EDFS_LOCK_TIMEOUT_SEC,
                        poll: float = EDFS_LOCK_POLL_SEC) -> bool:
    """Ver12.1 機能B: 排他ロックを取得する。

    fcntl は使わず、O_CREAT|O_EXCL によるロックファイル作成のみで実装する
    （Android共有ストレージ／SDカード等、fcntlのアドバイザリロックが正しく
    機能しない場合があるため、より移植性の高い方式を採用）。
    複数端末が同一の累積CSV（マージ後のマスターファイル含む）へ同時に
    書き込んでも行が破損・混在しないようにするための仕組み。
    """
    start = time.time()
    while True:
        try:
            fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            return True
        except FileExistsError:
            if time.time() - start > timeout:
                return False
            time.sleep(poll)
        except OSError:
            return False


def _release_file_lock(lock_path: Path) -> None:
    """Ver12.1 機能B: ロックファイルを解放する。"""
    try:
        os.remove(str(lock_path))
    except OSError:
        pass


def _csv_append(path: Path, df_new: pd.DataFrame) -> None:
    """CSVへの安全追記。存在しない場合はheader付きで新規作成。

    Ver12.1: 複数端末が同一パス（統合先マスターファイル）へ同時に追記しても
    行が壊れないよう、<path>.lock による排他制御を追加（機能B）。
    ロック取得に失敗した場合は今回分の追記をスキップする
    （致命的なデータ損失にはならず、次サイクルで再試行される）。
    このロック規約(<file>.lock)は es_data_merge_accelerator.py 側の
    merge_cumulative_logs() とも共通化しており、マージツールによる
    読み込み・書き戻し中の同時追記も防ぐ。

    Ver13.17修正: 従来 path.exists() の真偽で mode='w'(新規上書き)か
    mode='a'(追記)かを切り替えていたが、Android共有ストレージでは
    直前に書いたファイルへの exists() が一時的に False を返すことがあり、
    その瞬間に mode='w' が選ばれると既存データが丸ごと消える重大な
    データロス経路になっていた。mode は常に 'a' に固定し、ヘッダー要否
    だけをファイルサイズ(0バイト=新規)で判定する。万一サイズ判定を
    誤っても「ヘッダー行が1行余分に混ざる」だけで済み、既存データの
    上書き消失は原理的に起こらない。
    """
    lock_path = path.with_suffix(path.suffix + '.lock')
    if not _acquire_file_lock(lock_path):
        return
    try:
        need_header = True
        try:
            need_header = (not path.exists()) or path.stat().st_size == 0
        except OSError:
            need_header = not path.exists()
        df_new.to_csv(path, mode='a', header=need_header, index=False)
    finally:
        _release_file_lock(lock_path)


_edfs_last_merge_ts: float = 0.0
_edfs_merge_module = None
_edfs_merge_module_load_failed: bool = False


def _find_merge_accelerator_module():
    """本体と同じフォルダ等から es_data_merge_accelerator.py を探し、動的import する。

    Ver12.1対応版の auto_merge_once() が定義されているものだけを有効とみなす。
    見つからない／読み込みに失敗した場合は None を返し、本体は単体運用として
    従来どおり動作を継続する（マージ機能は必須ではない）。
    """
    global _edfs_merge_module, _edfs_merge_module_load_failed
    if _edfs_merge_module is not None:
        return _edfs_merge_module
    if _edfs_merge_module_load_failed:
        return None

    candidates = []
    try:
        here = Path(__file__).resolve().parent
        candidates.append(here / MERGE_ACCELERATOR_FILENAME)
    except Exception:
        pass
    try:
        candidates.append(Path.cwd() / MERGE_ACCELERATOR_FILENAME)
    except Exception:
        pass

    for p in candidates:
        try:
            if not p.exists():
                continue
            import importlib.util
            spec = importlib.util.spec_from_file_location("es_data_merge_accelerator", p)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            if hasattr(mod, "auto_merge_once"):
                _edfs_merge_module = mod
                return mod
        except Exception as e:
            print(f"  警告: マージツール読み込みに失敗しました（{p}）: {e}")

    _edfs_merge_module_load_failed = True
    return None


def maybe_run_scheduled_merge(log_dir: Path, use_color: bool = True, force: bool = False) -> None:
    """Ver12.1 機能A: EDFS_MERGE_INTERVAL_SEC（既定15分）ごとに他端末ログとの
    自動統合を試みる。es_data_merge_accelerator.py が同フォルダに無ければ
    何もしない（単体運用との後方互換性を維持）。

    本体のメインループから毎サイクル呼び出す想定で、内部で経過時間を見て
    間引く（force=True の場合は間隔を無視して即実行＝初回起動時処理用）。
    """
    global _edfs_last_merge_ts
    now = time.time()
    if not force and (now - _edfs_last_merge_ts) < EDFS_MERGE_INTERVAL_SEC:
        return
    _edfs_last_merge_ts = now

    mod = _find_merge_accelerator_module()
    if mod is None:
        return

    try:
        result = mod.auto_merge_once(data_dir=log_dir, edfs_module=sys.modules[__name__], silent=True)
    except Exception as e:
        print(f"  警告: 自動マージの実行に失敗しました ({e})")
        return

    if result and result.get("merged_any"):
        n_cum = sum(1 for r in result.get("cumulative", {}).values() if r.get("sources_merged", 0) > 0)
        print(color_by_state(
            f"[自動マージ] 他端末ログを統合しました "
            f"(セッション+{result.get('session_copied', 0)}件 / 累積ログ{n_cum}種を更新)",
            use_color=use_color))


def run_first_launch_migration_if_needed(log_dir: Path, use_color: bool = True) -> None:
    """Ver12.1 初回起動時処理。

    log_dir 直下に移行マーカー(EDFS_V12_MIGRATION_MARKER)が無い場合のみ、
    起動時に一度だけ以下を行う:
      1) es_data_merge_accelerator.py があれば自動マージを1回実行し、
         他端末ログ・旧バージョンの個別ログを統合先マスターファイルへ集約する。
      2) マージ済みで import_logs/_imported_done 配下に退避された旧ログ
         （統合済みで再度参照する意味がないもの）を削除して整理する。
    2回目以降の起動ではマーカーの存在により何もしない。
    マージツールが見つからない場合は統合をスキップし、マーカーのみ作成する
    （単体運用でも初回起動が失敗しないようにするため）。
    """
    marker = log_dir / EDFS_V12_MIGRATION_MARKER
    if marker.exists():
        return

    print("=" * 40)
    print(" Ver12.1 初回起動処理: 旧ログの統合・整理を確認しています...")

    mod = _find_merge_accelerator_module()
    if mod is not None:
        try:
            result = mod.auto_merge_once(data_dir=log_dir, edfs_module=sys.modules[__name__], silent=False)
            if result.get("merged_any"):
                print(" 初回マージ: 他端末／旧ログの統合が完了しました。")
            else:
                print(" 初回マージ: 統合対象の他端末ログは見つかりませんでした。")
        except Exception as e:
            print(f" 警告: 初回マージに失敗しました ({e})")

        if hasattr(mod, "purge_merged_archive"):
            try:
                n_deleted, size_deleted = mod.purge_merged_archive(data_dir=log_dir)
                if n_deleted:
                    print(f" 統合済み旧ログの削除: {n_deleted}件（約{size_deleted/1024:.1f}KB解放）")
                else:
                    print(" 削除対象の旧ログはありませんでした。")
            except Exception as e:
                print(f" 警告: 旧ログ削除に失敗しました ({e})")
    else:
        print(" es_data_merge_accelerator.py が見つからないため、統合はスキップしました。")
        print(" （単体運用の場合はこのメッセージのみで問題ありません）")

    try:
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(datetime.now().isoformat(), encoding='utf-8')
    except Exception as e:
        print(f" 警告: 移行マーカーの作成に失敗しました ({e})")
    print("=" * 40)


def save_edfs_area_history(area_df: pd.DataFrame, ft8_feats: Dict,
                            obs_dt: datetime, log_dir: Path) -> None:
    """機能①: edfs_area_history.csv に全エリアの現在状態を追記する。
    1サイクル×12エリア = 12行追記。
    """
    if area_df is None or area_df.empty:
        return
    rpi = ft8_feats.get('rpi', np.nan)
    rows = []
    for _, row in area_df.iterrows():
        rows.append({
            'timestamp':            obs_dt.strftime('%Y-%m-%d %H:%M:%S'),
            'area':                 row['area'],
            'rpi':                  rpi,
            'reach_rate':           row['reach_ratio'],
            'spot_count':           int(row['heard_count']),
            'fixed_station_count':  int(row['fixed_count']),
        })
    if rows:
        _csv_append(log_dir / EDFS_AREA_HISTORY_FILE, pd.DataFrame(rows))


CENTROID_MULTISAMPLE_HORIZONS = (0, 15, 30, 45, 60)


def save_edfs_centroid_history(spatial: Dict, obs_dt: datetime, log_dir: Path) -> None:
    """機能④⑤: edfs_centroid_history.csv に重心座標と速度を追記する。

    Ver13.13: マルチサンプリング化。Es を検出したサイクル（clat が有効）で
    現在値＋4ホライズン外挿の計5行を保存する。horizon_min 列を追加し、
    horizon_min==0 が実観測、それ以外は velocity による線形外挿。
    既存 CSV（horizon_min 列なし）とは _csv_append 側で自然にマージされる
    （旧行は読み込み時 horizon_min=NaN=実観測扱いになり後方互換）。
    """
    global _edfs_centroid_history_count, _edfs_prev_centroid
    clat = spatial.get('es_area_center_lat', np.nan)
    clon = spatial.get('es_area_center_lon', np.nan)
    if pd.isna(clat):
        return
    vel_lat = np.nan; vel_lon = np.nan
    if _edfs_prev_centroid:
        plat = _edfs_prev_centroid.get('clat', np.nan)
        plon = _edfs_prev_centroid.get('clon', np.nan)
        if not pd.isna(plat):
            vel_lat = clat - plat
            vel_lon = clon  - plon
    _edfs_prev_centroid = {'clat': clat, 'clon': clon}

    rows = []
    for h_min in CENTROID_MULTISAMPLE_HORIZONS:
        frac = h_min / 15.0
        if pd.isna(vel_lat):
            proj_lat = clat
            proj_lon = clon
        else:
            proj_lat = clat + vel_lat * frac
            proj_lon = clon + vel_lon * frac
        rows.append({
            'timestamp':    (obs_dt + timedelta(minutes=h_min)).strftime('%Y-%m-%d %H:%M:%S'),
            'centroid_lat': proj_lat,
            'centroid_lon': proj_lon,
            'velocity_lat': vel_lat,
            'velocity_lon': vel_lon,
            'horizon_min':  h_min,
        })
    _csv_append(log_dir / EDFS_CENTROID_HISTORY_FILE, pd.DataFrame(rows))
    _edfs_centroid_history_count += len(rows)


EVIDENCE_JSON_SUBDIR = "evidence"
EDFS_EVIDENCE_ANDROID_DIR = Path("/storage/emulated/0/EDFS/evidence")

LOG_CLEANUP_INTERVAL_SEC = 3600.0
EVIDENCE_JSON_MAX_AGE_HOURS = 6.0
EVIDENCE_JSON_MAX_COUNT = 150
LOG_CLEANUP_MAX_ROWS = {
    'ft8_spots_history.csv': 200000,
    'effective_contribution_history.csv': 100000,
    'feature_ranking_history.csv': 5000,
    'prediction_feature.csv': 100000,
    'teacher.csv': 20000,
    'prediction.csv': 20000,
    'prediction_archive.csv': 300000,
    'prediction_calibration_archive.csv': 300000,
    'edfs_ridge_training.csv': 500000,
}


def _cleanup_evidence_json(evidence_dir: Path, max_age_hours: float, max_count: int) -> int:
    if not evidence_dir.exists():
        return 0
    removed = 0
    try:
        files = list(evidence_dir.glob('evidence_*.json'))
    except Exception:
        return 0
    threshold = time.time() - max_age_hours * 3600.0
    remaining = []
    for f in files:
        try:
            if f.stat().st_mtime < threshold:
                f.unlink()
                removed += 1
            else:
                remaining.append(f)
        except Exception:
            continue
    if len(remaining) > max_count:
        try:
            remaining.sort(key=lambda p: p.stat().st_mtime)
        except Exception:
            pass
        for f in remaining[:len(remaining) - max_count]:
            try:
                f.unlink()
                removed += 1
            except Exception:
                continue
    return removed


def _tail_truncate_csv(path: Path, max_rows: int) -> int:
    if not path.exists():
        return 0
    try:
        with open(path, 'r', encoding='utf-8', errors='replace', newline='') as f:
            header = f.readline()
            if not header:
                return 0
            tail = _deque(maxlen=max_rows)
            total = 0
            for line in f:
                tail.append(line)
                total += 1
        if total <= max_rows:
            return 0
        tmp_path = path.with_suffix(path.suffix + '.cleanup_tmp')
        with open(tmp_path, 'w', encoding='utf-8', newline='') as f:
            f.write(header)
            f.writelines(tail)
        os.replace(tmp_path, path)
        return total - max_rows
    except Exception:
        return 0


def cleanup_unbounded_logs(log_dir: Path) -> Dict[str, int]:
    result: Dict[str, int] = {}
    try:
        evidence_dir = EDFS_EVIDENCE_ANDROID_DIR if is_android() else Path(log_dir) / EVIDENCE_JSON_SUBDIR
        n = _cleanup_evidence_json(evidence_dir, EVIDENCE_JSON_MAX_AGE_HOURS, EVIDENCE_JSON_MAX_COUNT)
        if n:
            result['evidence_json'] = n
    except Exception:
        pass
    for fname, cap in LOG_CLEANUP_MAX_ROWS.items():
        try:
            n = _tail_truncate_csv(Path(log_dir) / fname, cap)
            if n:
                result[fname] = n
        except Exception:
            pass
    return result


KNOWN_LEARNING_FILES = {
    EDFS_AREA_HISTORY_FILE, EDFS_CENTROID_HISTORY_FILE, EDFS_RIDGE_TRAINING_FILE,
    EDFS_PREDICTION_ERROR_FILE, TEACHER_CSV, PRED_FEATURE_CSV, PREDICTION_CSV,
    "dpp_calibration.json",
}


def _learning_merge_safe_root() -> Path:
    """フォルダブラウザで閲覧を許可する安全なルートディレクトリを決定する。"""
    if is_android():
        p = Path('/storage/emulated/0')
        if p.exists():
            return p
    try:
        home = Path.home()
        if home.exists():
            return home
    except Exception:
        pass
    return Path(DEFAULT_DATA_DIR).resolve()


def _is_within_safe_root(p: Path, root: Path) -> bool:
    """パストラバーサル対策: 指定パスが安全なルート配下にあることを検証する。"""
    try:
        p_real = str(p.resolve())
        root_real = str(root.resolve())
        return p_real == root_real or p_real.startswith(root_real + os.sep)
    except Exception:
        return False


def list_browse_dir(path_str: Optional[str]) -> Dict[str, Any]:
    """フォルダブラウザ用にディレクトリ内容を一覧化する(Rev16.1新設)。"""
    root = _learning_merge_safe_root()
    target = Path(path_str).expanduser() if path_str else root
    if not _is_within_safe_root(target, root):
        return {'ok': False, 'error': 'アクセスが許可されていないパスです', 'root': str(root)}
    if not target.exists() or not target.is_dir():
        return {'ok': False, 'error': 'フォルダが存在しないか、フォルダではありません', 'root': str(root)}
    entries = []
    try:
        for child in sorted(target.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower())):
            try:
                st = child.stat()
                entries.append({
                    'name': child.name,
                    'path': str(child),
                    'is_dir': child.is_dir(),
                    'size': st.st_size if child.is_file() else None,
                    'mtime': datetime.fromtimestamp(st.st_mtime).strftime('%Y-%m-%d %H:%M'),
                    'is_learning_file': (child.is_file() and child.name in KNOWN_LEARNING_FILES),
                })
            except Exception:
                continue
    except PermissionError:
        return {'ok': False, 'error': 'このフォルダへのアクセス権限がありません', 'root': str(root)}
    except Exception as e:
        return {'ok': False, 'error': str(e), 'root': str(root)}
    parent = None
    if target != root and _is_within_safe_root(target.parent, root):
        parent = str(target.parent)
    return {'ok': True, 'path': str(target), 'root': str(root), 'parent': parent, 'entries': entries}


def _generic_learning_merge_key_columns(columns: List[str]) -> Optional[List[str]]:
    candidates = ['timestamp', 'area', 'horizon_min', 'station']
    key = [c for c in candidates if c in columns]
    return key if 'timestamp' in key else None


def merge_learning_file(current_path: Path, incoming_path: Path) -> Dict[str, Any]:
    """他端末から選択された学習コアCSVを、現在の端末のファイルへ安全に統合する。
    timestamp(+area/horizon_min/station)をキーに重複行を除去し、現在の端末に
    まだ存在しない新規サンプルのみを追加する(=「学習が加速する部分だけ」)。
    """
    import shutil
    try:
        cur_df = pd.read_csv(current_path) if current_path.exists() else pd.DataFrame()
    except Exception as e:
        return {'before': 0, 'after': 0, 'added': 0, 'error': f'既存ファイル読込失敗: {e}'}
    try:
        new_df = pd.read_csv(incoming_path)
    except Exception as e:
        return {'before': len(cur_df), 'after': len(cur_df), 'added': 0,
                'error': f'統合元ファイル読込失敗: {e}'}

    before = len(cur_df)
    if new_df.empty:
        return {'before': before, 'after': before, 'added': 0, 'error': None}

    all_cols = list(dict.fromkeys(list(cur_df.columns) + list(new_df.columns)))
    for c in all_cols:
        if c not in cur_df.columns:
            cur_df[c] = np.nan
        if c not in new_df.columns:
            new_df[c] = np.nan
    cur_df = cur_df[all_cols] if not cur_df.empty else pd.DataFrame(columns=all_cols)
    new_df = new_df[all_cols]

    combined = pd.concat([cur_df, new_df], ignore_index=True)
    key_cols = _generic_learning_merge_key_columns(all_cols)
    if key_cols:
        combined = combined.drop_duplicates(subset=key_cols, keep='first')
        if 'timestamp' in combined.columns:
            try:
                combined = combined.sort_values('timestamp', kind='stable')
            except Exception:
                pass
    else:
        combined = combined.drop_duplicates(keep='first')

    after = len(combined)
    added = after - before

    if added > 0:
        try:
            if current_path.exists():
                backup_path = current_path.with_suffix(current_path.suffix + '.bak')
                shutil.copy2(current_path, backup_path)
        except Exception:
            pass
        try:
            tmp_path = current_path.with_suffix(current_path.suffix + '.merge_tmp')
            combined.to_csv(tmp_path, index=False)
            os.replace(tmp_path, current_path)
        except Exception as e:
            return {'before': before, 'after': before, 'added': 0,
                    'error': f'書き戻し失敗: {e}'}

    return {'before': before, 'after': int(after), 'added': int(added), 'error': None}


def merge_duct_calibration_file(current_path: Path, incoming_path: Path) -> Dict[str, Any]:
    """ESDUCT(整理版): 他端末のDUCT較正データ(dpp_calibration.json)を、現在の
    端末の較正データへ安全に統合する。CSVベースの merge_learning_file() とは
    データ構造が根本的に異なる(行の集合ではなく、セル別Bayesianパラメータの
    辞書)ため、専用のマージロジックを用意した。

    【マージ方式】各セルのalpha/beta(ベータ分布のパラメータ)は、事前分布
    (prior.alpha0/beta0)からの積み増し量として意味を持つ。したがって、
    2台の端末で独立に観測されたセルのalpha/beta増分は単純加算してよい
    (ベータ分布の共役性により、独立な観測をまとめて反映したのと数学的に
    等価になる)。この考え方に基づき、
      増分 = 相手セルのalpha - alpha0 (0未満は切り捨て)
    を自セルへ加算する形で統合する(既存セルの学習内容を上書き・破棄しない)。
    n_obs(観測回数)は単純合算、last_updateは新しい方を採用する。
    """
    import shutil
    try:
        mine = load_calibration(str(current_path))
    except Exception as e:
        return {'before': 0, 'after': 0, 'added': 0, 'error': f'既存ファイル読込失敗: {e}'}
    try:
        with open(incoming_path, 'r', encoding='utf-8') as f:
            other = json.load(f)
    except Exception as e:
        n0 = len(mine.get('cells', {}))
        return {'before': n0, 'after': n0, 'added': 0, 'error': f'統合元ファイル読込失敗: {e}'}
    if not isinstance(other, dict) or 'cells' not in other:
        n0 = len(mine.get('cells', {}))
        return {'before': n0, 'after': n0, 'added': 0,
                'error': '統合元がDUCT較正データ(dpp_calibration.json)の形式ではありません'}

    with _CAL_IO_LOCK:
        prior = mine.get('prior', {'alpha0': 0.6, 'beta0': 6.0})
        a0, b0 = prior.get('alpha0', 0.6), prior.get('beta0', 6.0)
        cells = mine.setdefault('cells', {})
        before_n = len(cells)
        before_obs = sum(c.get('n_obs', 0) for c in cells.values())
        new_cnt = upd_cnt = 0
        for k, oc in (other.get('cells', {}) or {}).items():
            try:
                oa = oc.get('alpha', a0) - a0
                ob = oc.get('beta', b0) - b0
            except Exception:
                continue
            if k not in cells:
                cells[k] = {'alpha': a0 + max(0.0, oa), 'beta': b0 + max(0.0, ob),
                            'last_update': oc.get('last_update'), 'n_obs': oc.get('n_obs', 0)}
                new_cnt += 1
            else:
                mc = cells[k]
                mc['alpha'] = mc.get('alpha', a0) + max(0.0, oa)
                mc['beta'] = mc.get('beta', b0) + max(0.0, ob)
                mc['n_obs'] = mc.get('n_obs', 0) + oc.get('n_obs', 0)
                ot = oc.get('last_update'); mt = mc.get('last_update')
                if ot and (not mt or ot > mt):
                    mc['last_update'] = ot
                upd_cnt += 1
        after_n = len(cells)
        after_obs = sum(c.get('n_obs', 0) for c in cells.values())
        if new_cnt or upd_cnt:
            try:
                if current_path.exists():
                    backup_path = current_path.with_suffix(current_path.suffix + '.bak')
                    shutil.copy2(current_path, backup_path)
            except Exception:
                pass
            try:
                save_calibration(mine, str(current_path))
            except Exception as e:
                return {'before': before_n, 'after': before_n, 'added': 0,
                        'error': f'書き戻し失敗: {e}'}

    return {'before': before_n, 'after': after_n, 'added': int(new_cnt),
            'updated': int(upd_cnt), 'before_obs': round(before_obs, 1),
            'after_obs': round(after_obs, 1), 'error': None}


def run_learning_merge(log_dir: Path, file_paths: List[str]) -> Dict[str, Any]:
    """選択された複数ファイルをまとめて統合し、結果一覧を返す(Rev16.1新設)。
    ESDUCT(整理版): dpp_calibration.json(DUCT較正データ)はCSVではなくJSON
    構造のため、専用のmerge_duct_calibration_file()へ振り分ける。
    """
    root = _learning_merge_safe_root()
    results = []
    total_added = 0
    duct_cal_name = "dpp_calibration.json"
    for fp in file_paths:
        try:
            incoming = Path(fp).expanduser()
        except Exception:
            results.append({'file': str(fp), 'ok': False, 'error': '不正なパスです'})
            continue
        if not _is_within_safe_root(incoming, root):
            results.append({'file': str(fp), 'ok': False, 'error': 'アクセスが許可されていないパスです'})
            continue
        if not incoming.exists() or not incoming.is_file():
            results.append({'file': str(fp), 'ok': False, 'error': 'ファイルが見つかりません'})
            continue
        basename = incoming.name
        if basename not in KNOWN_LEARNING_FILES:
            results.append({'file': basename, 'ok': False,
                             'error': '学習コアデータとして認識できないファイル名のため対象外です'})
            continue
        if basename == duct_cal_name:
            current_path = Path(CALIBRATION_FILE)
            r = merge_duct_calibration_file(current_path, incoming)
        else:
            current_path = Path(log_dir) / basename
            r = merge_learning_file(current_path, incoming)
        r['file'] = basename
        r['ok'] = (r.get('error') is None)
        results.append(r)
        total_added += r.get('added', 0)
    return {'results': results, 'total_added': total_added}


def _sanitize_for_json(obj):
    """Ver13.23: NaN/InftyをJSON仕様準拠のnullへ変換する(再帰的)。
    Python標準jsonはNaN/Infinityを裸トークンとして出力するが、これはRFC 8259上は
    不正なJSONであり、Gemini API等の厳格なJSONパーサーでパースエラーになる。
    Evidence JSONは今後LLMへ渡すことが前提のため、ここで正規化しておく。
    """
    if isinstance(obj, float):
        return None if (pd.isna(obj) or math.isinf(obj)) else obj
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    return obj


def write_evidence_json(snapshot: Dict, obs_dt: datetime, log_dir: Path) -> None:
    """Ver13.23: EDFS/EICA共通研究基盤向けEvidence JSONを1サイクル1ファイルで書き出す。

    設計方針(ionosphere-research-platform):
    - LLM(Gemma/Gemini)に計算させず、EDFS側で「確定した事実」のみをここでJSON化する。
    - 現場(Android/RPi)はこのJSONを蓄積するだけ。USB接続時にWindows側へ集約される。
    - 書き込みはSessionLogger.flush()と同じ .tmp→os.replace() の原子的置換パターンを踏襲し、
      Pydroid3のOSキル等で書き込み中断してもevidenceディレクトリ全体を壊さないようにする。
    - ファイル数最小化のコンセプトに合わせ、evidenceディレクトリ配下に1サイクル1ファイルのみ生成し、
      集約・重複排除・DB化はWindows側(SQLite取り込み層)の責務として持たせない。
    - Android上では出力先をlog_dirに委ねず、GPSブリッジ(/storage/emulated/0/EDFS/)と
      同じ共有ストレージ配下に固定する。ランチャー起動(CLI引数なし)ではDEFAULT_DATA_DIR
      が相対パスのままアプリ内部サンドボックスへ解決される場合があり、その場合Filesアプリ
      からもUSB/MTP経由のWindows側からも一切見えなくなるため。RPi/Windows側ではlog_dir配下
      にそのまま出力する(サンドボックス制約がないため)。
    """
    try:
        evidence_dir = EDFS_EVIDENCE_ANDROID_DIR if is_android() else (log_dir / EVIDENCE_JSON_SUBDIR)
        evidence_dir.mkdir(parents=True, exist_ok=True)
        fname = f"evidence_{obs_dt.strftime('%Y%m%d_%H%M%S')}.json"
        final_path = evidence_dir / fname
        tmp_path = final_path.with_suffix('.tmp')
        with open(tmp_path, 'w', encoding='utf-8') as f:
            json.dump(_sanitize_for_json(snapshot), f, ensure_ascii=False, allow_nan=False, default=str)
        os.replace(tmp_path, final_path)
    except Exception:
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass


def build_area_forecast(area_df: pd.DataFrame, area_fvecs: Dict[str, List[float]],
                         area_geo_eta: Optional[Dict[str, Dict]] = None,
                         area_climatology: Optional[Dict[str, Dict]] = None) -> Dict[str, Dict]:
    """Ver13.25: エリア別将来予測をEvidence Schema向けにまとめる。

    predict_area_reach()自体はUI表示(print_area_reachability_panel)で既に呼ばれている
    学習済みモデルの推論結果であり、ここで新規に学習・計算するものではない
    (Evidence JSONは出力専用。学習ループへの影響はない)。

    Rev14.0: 以下を追加(いずれも参考情報であり、Ridgeの選択・ゲーティング
    ロジックには一切影響しない):
      - eta_min_geometric / eta_approach_cos: 学習器に依存しない幾何学的
        到達目安(geometric_eta_estimate)。conf=='LEARNING'でも算出できる
        ため、学習済みモデルが無いエリアでも「何もない」状態を避けられる。
      - is_compensated / compensation_basis: 表示中の到達率がIDW補間/
        ブレンドによるものかどうかの透明性情報(interpolate_area_reachability_idw)。
      - climatology: Rev14.1の②自己ログ/③NICT統計による季節気候値
        (get_climatology_prior()の結果、{'prob_pct','source','n','detail'})。
        conf=='LEARNING'で数値予測が無い場合の参考値として利用できる。

    Rev14.3: 以下を追加(いずれもSHADOW=本番pred_rpiには一切影響しない診断・実験系統。
    詳細は各関数のdocstring参照):
      - geometry_ridge_consistency: Ridge予測と幾何ETAの結論が一致しているか
        (CONSISTENT/CONFLICT/NEUTRAL)。
      - climatology_anomaly_pt: Ridge予測が季節気候値からどれだけ乖離しているか(%pt)。
      - pred_rpi_ensemble / ensemble_mode / weights / geometry_adjustment:
        Skill×CoverageベースのモデルウェイトでRidgeとpersistenceを混合した
        シャドー予測値。本番pred_rpiとは別フィールドであり、EDFS自身の判断や
        学習ループには一切使われない(実運用データでの追跡評価専用)。
    """
    forecast = {}
    if area_df is None or area_df.empty:
        return forecast
    area_geo_eta = area_geo_eta or {}
    area_climatology = area_climatology or {}
    for _, row in area_df.iterrows():
        area = row.get('area')
        if not area:
            continue
        geo = area_geo_eta.get(area, {})
        common_extra = {
            "eta_min_geometric": geo.get("eta_min"),
            "eta_approach_cos": geo.get("approach_cos"),
            "eta_confidence": geo.get("eta_confidence", "LOW"),
            "eta_geometric_method": "geometric_linear_extrapolation",
            "eta_geometric_is_forecast_probability": False,
            "is_compensated": bool(row.get("is_compensated", False)),
            "compensation_basis": row.get("compensation_basis", "raw"),
            "ratio_source": "idw" if bool(row.get("is_compensated", False)) else "observed",
            "idw_used": bool(row.get("is_compensated", False)),
            "idw_source_count": int(row.get("compensation_n_sources", 0) or 0),
            "climatology": area_climatology.get(area),
            "holdout_method": (_area_holdout_diag.get(area, {}) or {}).get('method'),
            "holdout_n_folds": (_area_holdout_diag.get(area, {}) or {}).get('n_folds'),
        }
        r2 = _area_ridge_models.get(area, (None, None, None, 0.0))[3] if area in _area_ridge_models else 0.0
        conf = area_pred_confidence(area, r2)
        persistence_rpi = row.get("reach_ratio_compensated", row.get("reach_ratio"))
        if isinstance(persistence_rpi, float) and pd.isna(persistence_rpi):
            persistence_rpi = None
        eta_min = geo.get("eta_min")

        def _shadow_fields(ridge_pred_rpi: Optional[float], fvec_: Optional[List[float]]) -> Dict:
            """Rev14.3(#③④⑤): Ridge予測が確定した時点でのみ意味を持つ
            シャドー系フィールドをまとめて計算する(LEARNING等でridge_pred_rpiが
            Noneの場合は各関数側がNEUTRAL/None/空を返す設計)。"""
            consistency = geometry_ridge_consistency(eta_min, ridge_pred_rpi)
            anomaly = climatology_anomaly_pt(ridge_pred_rpi, area_climatology.get(area))
            ensemble = compute_shadow_ensemble(area, fvec_, ridge_pred_rpi, persistence_rpi,
                                                eta_min, consistency)
            return {
                "geometry_ridge_consistency": consistency,
                "climatology_anomaly_pt": anomaly,
                "ensemble_shadow_eval": summarize_ensemble_shadow_eval(area),
                **ensemble,
            }

        if conf == 'LEARNING':
            forecast[area] = {"pred_rpi": None, "confidence": "LEARNING", "bias": None,
                               "is_forecast": False, **common_extra, **_shadow_fields(None, None)}
            continue
        fvec = area_fvecs.get(area)
        if not fvec:
            forecast[area] = {"pred_rpi": None, "confidence": conf, "bias": None,
                               "is_forecast": False, **common_extra, **_shadow_fields(None, None)}
            continue
        pred_c, _, bias = predict_area_reach(area, fvec)
        pred_is_valid = not pd.isna(pred_c)
        forecast[area] = {
            "pred_rpi": pred_c if pred_is_valid else None,
            "confidence": conf,
            "bias": None if pd.isna(bias) else bias,
            "is_forecast": pred_is_valid,
            **common_extra,
            **_shadow_fields(pred_c if pred_is_valid else None, fvec),
        }
    return forecast


def build_evidence_snapshot(obs_dt: datetime, spatial: Dict, ft8_feats: Dict,
                             eqi_repr_val: float, eqi_conf_val: float,
                             sw_score_val: float, trace_quality_val: float,
                             prediction_confidence: float, es_phase: str,
                             cycle_idx: int, forecast: Optional[Dict] = None,
                             directional_field: Optional[Dict] = None) -> Dict:
    """Ver13.23: 1サイクル分のEvidence Schemaを組み立てる。

    フィールド構成はionosphere-research-platform構想のEvidence Schemaに準拠。
    ここに含めるのはEDFSが既に確定させた値のみ(LLMへ計算させる余地を残さない)。

    Rev14.5追加: directional_field(combine_directional_evidence()の結果)を
    受け取り、"directional"キーとして出力する。命名方針: FT8観測が無い場合
    でも常に出力されるのは'directional_opening_expectancy'(NICT4局の1次元
    勾配のみに基づく粗い期待度)であり、これを「伝播方位」とは呼ばない。
    FT8観測が加わって初めて'propagation_directional_confidence'が入る
    (Noneでなくなる)。詳細はcombine_directional_evidence()のdocstring参照。
    """
    return {
        "schema_version": "1.0",
        "timestamp": obs_dt.strftime('%Y-%m-%d %H:%M:%S'),
        "cycle_idx": cycle_idx,
        "domain": "edfs",
        "state": {
            "intensity": eqi_repr_val,
            "rpi": ft8_feats.get('rpi_norm', np.nan),
            "area_reach": ft8_feats.get('es_mode_score', np.nan),
            "centroid": [spatial.get('es_area_center_lat', np.nan),
                         spatial.get('es_area_center_lon', np.nan)],
            "centroid_speed_kmh": spatial.get('es_area_speed', np.nan),
            "centroid_direction": spatial.get('es_area_direction', np.nan),
            "spread": spatial.get('es_area_spread', np.nan),
            "es_count": spatial.get('es_area_count', np.nan),
            "es_type": ver113_es_type_label(),
            "phase": es_phase,
            "space_weather_score": sw_score_val,
            "eqi": eqi_repr_val,
        },
        "prediction": {
            "confidence": prediction_confidence,
        },
        "forecast": forecast or {},
        "directional": directional_field or {
            'directional_opening_expectancy': {name: 0.5 for name, _ in COMPASS_DIRECTIONS_8},
            'propagation_directional_confidence': None,
            'primary_evidence_source': 'none',
            'nict_gradient': None, 'ft8_anisotropy': None,
        },
        "observation_quality": {
            "eqi_confidence": eqi_conf_val,
            "trace_quality": trace_quality_val,
        },
    }


def build_gemma_prompt(snapshot: Dict) -> str:
    """Ver13.23: Evidence SnapshotからGemma(現場層)向けの日本語プロンプトを生成する。
    AI Edge Galleryは現時点で外部アプリからのIntent/API連携を公開していないため、
    自動送信ではなく「コピー&ペーストしてAI Chat/Prompt Labへ貼るだけ」を想定した文面にする。
    将来LiteRT-LM等での自動化に置き換える場合も、プロンプト生成部はそのまま再利用できる。
    """
    st = snapshot.get("state", {})
    pred = snapshot.get("prediction", {})

    def _fmt(v, unit=""):
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return "不明"
        return f"{v}{unit}"

    lat, lon = (st.get("centroid") or [None, None])[:2]
    lat_is_nan = isinstance(lat, float) and pd.isna(lat)
    lon_is_nan = isinstance(lon, float) and pd.isna(lon)
    centroid_str = "不明" if (lat is None or lon is None or lat_is_nan or lon_is_nan) else f"({lat}, {lon})"

    forecast = snapshot.get("forecast") or {}
    if forecast:
        lines = []
        for area, f in forecast.items():
            if f.get("confidence") == "LEARNING" or f.get("pred_rpi") is None:
                lines.append(f"  {area}: 学習中")
            else:
                lines.append(f"  {area}: 到達率予測 {f['pred_rpi']*100:.0f}%(信頼度:{f.get('confidence')})")
        forecast_str = "\n".join(lines) if lines else "不明"
    else:
        forecast_str = "不明"

    directional = snapshot.get("directional") or {}
    _dir_field = directional.get("propagation_directional_confidence") or \
        directional.get("directional_opening_expectancy")
    if _dir_field:
        _best_dir = max(_dir_field, key=lambda k: _dir_field[k])
        _src = directional.get("primary_evidence_source")
        if _src == "ft8_observed":
            directional_str = f"{_best_dir}方向(FT8実測を含む推定)"
        elif _src == "nict_model_only":
            directional_str = f"{_best_dir}方向寄り(NICT観測局データからの粗い推定、実測ではない)"
        else:
            directional_str = "不明(情報不足)"
    else:
        directional_str = "不明"

    return (
        "以下はEs(スポラディックE層)伝播の観測状態です。無線家向けに2〜3文で、"
        "何が起きていて今後どうなりそうかを平易な日本語で解説してください。\n"
        f"- 観測時刻: {snapshot.get('timestamp', '不明')}\n"
        f"- Esタイプ: {_fmt(st.get('es_type'))}\n"
        f"- Phase(状態遷移段階): {_fmt(st.get('phase'))}\n"
        f"- Es強度(EQI): {_fmt(st.get('intensity'))}\n"
        f"- RPI: {_fmt(st.get('rpi'))}\n"
        f"- 伝播モードスコア(area_reach): {_fmt(st.get('area_reach'))}\n"
        f"- 重心位置: {centroid_str}\n"
        f"- 重心移動速度: {_fmt(st.get('centroid_speed_kmh'), 'km/h')}\n"
        f"- 広がり(spread): {_fmt(st.get('spread'))}\n"
        f"- 開けやすい方位の目安: {directional_str}\n"
        f"- 宇宙天気スコア: {_fmt(st.get('space_weather_score'))}\n"
        f"- 予測信頼度: {_fmt(pred.get('confidence'), '%')}\n"
        f"- エリア別将来予測(到達率):\n{forecast_str}\n"
    )


def write_gemma_prompt(snapshot: Dict, log_dir: Path) -> None:
    """Ver13.23: 最新1サイクル分のGemma向けプロンプトを単一ファイルへ上書き保存する。
    ファイル数最小化のコンセプトに合わせ、evidenceのような1サイクル1ファイルにせず
    常に同じファイル名(gemma_prompt_latest.txt)を上書きする。
    """
    try:
        target_dir = EDFS_EVIDENCE_ANDROID_DIR if is_android() else (log_dir / EVIDENCE_JSON_SUBDIR)
        target_dir.mkdir(parents=True, exist_ok=True)
        final_path = target_dir / "gemma_prompt_latest.txt"
        tmp_path = final_path.with_suffix('.tmp')
        with open(tmp_path, 'w', encoding='utf-8') as f:
            f.write(build_gemma_prompt(snapshot))
        os.replace(tmp_path, final_path)
    except Exception:
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass


_LATEST_EVIDENCE_SNAPSHOT: Dict = {}
_MCP_SERVER_STARTED = False
_MCP_HOST = "127.0.0.1"
_MCP_PORT = 8765
_MCP_TOOL_NAME = "get_latest_es_evidence"
_MCP_TOOL_DESCRIPTION = (
    "EDFSが観測した最新のEs(スポラディックE層)伝播状態と、エリア別の将来到達率予測を返す。"
    "強度(EQI)、RPI、伝播モードスコア、重心位置・移動速度・広がり、"
    "Esタイプ(停滞型/移動型/拡大型/多層型)、Phase(状態遷移段階)、"
    "予測信頼度、観測品質、エリア別forecast(到達率予測・信頼度・バイアス)を含む。引数は不要。"
)


def _mcp_get_latest_evidence_text() -> str:
    """メモリ上の最新evidence_snapshotをそのままJSON文字列で返す(ファイルI/O不要)。"""
    if not _LATEST_EVIDENCE_SNAPSHOT:
        return json.dumps(
            {"error": "evidenceがまだ生成されていません。EDFSが1サイクル以上動作しているか確認してください。"},
            ensure_ascii=False,
        )
    snap = dict(_LATEST_EVIDENCE_SNAPSHOT)
    snap["_mcp_fetched_at"] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return json.dumps(snap, ensure_ascii=False)


class _MCPRequestHandler(http.server.BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):
        pass

    def _send_json(self, obj, status=200):
        body = json.dumps(obj, ensure_ascii=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path in ('/', '/health'):
            self._send_json({"status": "ok", "server": "edfs-evidence-mcp"})
        else:
            self._send_json({"error": "not found"}, status=404)

    def do_POST(self):
        try:
            length = int(self.headers.get('Content-Length', 0))
            raw = self.rfile.read(length) if length > 0 else b'{}'
            req = json.loads(raw.decode('utf-8'))
        except Exception:
            self._send_json({"jsonrpc": "2.0", "id": None,
                              "error": {"code": -32700, "message": "Parse error"}}, status=400)
            return

        method = req.get('method')
        req_id = req.get('id')
        if req_id is None:
            self.send_response(202)
            self.send_header('Content-Length', '0')
            self.end_headers()
            return

        if method == 'initialize':
            self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {
                "protocolVersion": "2026-07-28",
                "capabilities": {"tools": {}},
                "serverInfo": {"name": "edfs-evidence-mcp", "version": "1.0"},
            }})
        elif method == 'tools/list':
            self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {"tools": [{
                "name": _MCP_TOOL_NAME, "description": _MCP_TOOL_DESCRIPTION,
                "inputSchema": {"type": "object", "properties": {}, "required": []},
            }]}})
        elif method == 'tools/call':
            name = req.get('params', {}).get('name')
            if name == _MCP_TOOL_NAME:
                text = _mcp_get_latest_evidence_text()
                self._send_json({"jsonrpc": "2.0", "id": req_id,
                                  "result": {"content": [{"type": "text", "text": text}], "isError": False}})
            else:
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {
                    "content": [{"type": "text", "text": f"未知のtool: {name}"}], "isError": True}})
        else:
            self._send_json({"jsonrpc": "2.0", "id": req_id,
                              "error": {"code": -32601, "message": f"Method not found: {method}"}}, status=404)


def start_mcp_server_once() -> None:
    """Ver13.24: MCPサーバーを一度だけdaemonスレッドとして起動する(冪等)。
    起動失敗(ポート使用中等)はEDFS本体の予測動作に影響させず、ログのみ出す。
    """
    global _MCP_SERVER_STARTED
    if _MCP_SERVER_STARTED:
        return
    _MCP_SERVER_STARTED = True
    try:
        httpd = http.server.ThreadingHTTPServer((_MCP_HOST, _MCP_PORT), _MCPRequestHandler)
        t = threading.Thread(target=httpd.serve_forever, daemon=True, name="edfs-mcp-server")
        t.start()
        print(f"[EDFS] MCPサーバー起動: http://{_MCP_HOST}:{_MCP_PORT}/mcp (Gemma連携用)")
    except Exception as e:
        print(f"[EDFS] MCPサーバー起動失敗(evidence機能自体には影響なし): {e}")


def save_edfs_prediction_error(area: str, predicted_rpi: float, actual_rpi: float,
                                obs_dt: datetime, log_dir: Path) -> None:
    """機能⑨: edfs_prediction_error.csv に予測誤差を追記する。"""
    _csv_append(log_dir / EDFS_PREDICTION_ERROR_FILE, pd.DataFrame([{
        'timestamp':     obs_dt.strftime('%Y-%m-%d %H:%M:%S'),
        'area':          area,
        'predicted_rpi': predicted_rpi,
        'actual_rpi':    actual_rpi,
        'error':         actual_rpi - predicted_rpi,
    }]))


def load_edfs_csv_counts(log_dir: Path) -> Tuple[int, int]:
    """起動時に既存CSVの行数を読み込んでカウンタを初期化する。
    戻り値: (area_teachers, centroid_history_count)
    """
    at = 0; ch = 0
    try:
        p = log_dir / EDFS_AREA_HISTORY_FILE
        if p.exists():
            n = sum(1 for _ in open(p, encoding='utf-8')) - 1
            at = max(0, n - 12)
    except Exception:
        pass
    try:
        p = log_dir / EDFS_CENTROID_HISTORY_FILE
        if p.exists():
            try:
                cdf = pd.read_csv(p)
                ch = len(cdf)
            except Exception:
                ch = sum(1 for _ in open(p, encoding='utf-8')) - 1
    except Exception:
        pass
    return at, ch



def _robust_read_csv_dropping_tail(path: Path, parse_dates: Optional[List[str]] = None,
                                    max_drop_lines: int = 5) -> Optional[pd.DataFrame]:
    """Rev14.8(堅牢性): CSV全体のパースに失敗した場合、末尾の壊れた行
    (アプリ強制終了・プロセスkill・ストレージ異常等による書き込み中断で
    発生しうる)を1行ずつ捨てながら再試行する。

    【想定する障害】攻撃ではなく通常運用中に普通に起こりうる事象:
    Pydroid3強制終了、Androidのプロセスkill、Windows終了、ストレージI/O異常、
    CSV書き込み途中の電源断など。追記型ログ(_csv_append())の構造上、
    壊れうるのは基本的に最後に書き込み中だった行のみであるため、
    末尾から少数行を落とすだけで残りの資産(学習履歴)を救済できる。

    通常のパース(1回目)が成功すればそのまま返す。失敗した場合のみ、
    テキストとして読み直し、末尾からmax_drop_lines行までを順に落として
    再パースを試みる。どうしても成功しなければNoneを返す
    (呼び出し側で「復元0件」として扱う、既存動作と同じ)。
    """
    try:
        return pd.read_csv(path, parse_dates=parse_dates)
    except Exception:
        pass
    try:
        raw_lines = path.read_text(encoding='utf-8', errors='replace').splitlines()
    except Exception:
        return None
    if len(raw_lines) < 2:
        return None
    for drop_n in range(1, min(max_drop_lines, len(raw_lines) - 1) + 1):
        trial_text = '\n'.join(raw_lines[:-drop_n])
        try:
            df = pd.read_csv(io.StringIO(trial_text), parse_dates=parse_dates)
            print(f"  [EDFS] {path.name}: 末尾{drop_n}行が破損していたため切り捨てて復元しました"
                  f"(残り{len(df)}行)")
            return df
        except Exception:
            continue
    return None


def bootstrap_ridge_training_from_csv(log_dir: Path,
                                       max_rows_per_area: int = 2000
                                       ) -> Tuple[int, int]:
    """Rev14.4新設: edfs_ridge_training.csv(実際にRidgeへ入力された特徴量
    ベクトルそのものを保存する専用ログ、resolve_pending_samples()参照)から
    _area_forecast_history を復元する。

    _bootstrap_area_history_from_csv_legacy_approx()と異なり、rate_change/
    acceleration等を「推定」する必要がない(そのサイクルで実際に使われた値が
    そのまま入っている)。空間特徴量やcentroid特徴量も0埋めではなく実値。

    スキーマ自己記述形式による安全なマイグレーション:
      各行は書き込み時点のAREA_FEATURE_NAMES(schema_names列、JSON配列)と
      そのときの値(feature_json列)を両方保持する。読み込み時に:
        - schema_names == 現在のAREA_FEATURE_NAMES        → そのまま使用
        - schema_names が現在のAREA_FEATURE_NAMESの前方一致
          (=既存次元の並びは変えず末尾に新次元を追加しただけ、
           実際にVer14.0のcentroid特徴量追加もこのパターン)
                                                            → 末尾を0埋めして次元を揃える
        - それ以外(並び替え・削除など安全に整列できない変更)  → その行は使わない
      「わからないものは0埋めで無理に使う」のではなく「安全に整列できない
      ものは使わない」という判断基準にすることで、異なる意味の特徴量を
      誤って同じ次元として学習させる事故を防ぐ。

    戻り値: (復元教師数, スキップ行数)
    """
    global _area_forecast_history, _area_total_samples, _area_es_samples
    path = log_dir / EDFS_RIDGE_TRAINING_FILE
    if not path.exists():
        return 0, 0
    df = _robust_read_csv_dropping_tail(path, parse_dates=['timestamp'])
    if df is None:
        return 0, 0
    required_cols = {'area', 'schema_names', 'feature_json', 'target'}
    if df.empty or not required_cols.issubset(df.columns):
        return 0, 0
    df = df.sort_values('timestamp').reset_index(drop=True)
    current_schema = AREA_FEATURE_NAMES

    restored_total = 0
    skipped = 0
    es_like = 0
    by_area: Dict[str, List[Tuple[List[float], float, bool]]] = {}
    for _, row in df.iterrows():
        try:
            stored_names = json.loads(row['schema_names'])
            stored_vals = json.loads(row['feature_json'])
            if not isinstance(stored_names, list) or not isinstance(stored_vals, list) or \
                    len(stored_names) != len(stored_vals):
                skipped += 1
                continue
            if stored_names == current_schema:
                fvec = [float(v) for v in stored_vals]
            elif len(stored_names) < len(current_schema) and \
                    stored_names == current_schema[:len(stored_names)]:
                fvec = [float(v) for v in stored_vals] + [0.0] * (len(current_schema) - len(stored_names))
            else:
                skipped += 1
                continue
            if len(fvec) != len(current_schema):
                skipped += 1
                continue
            if not all(math.isfinite(v) for v in fvec):
                skipped += 1
                continue
            target = row.get('target')
            if pd.isna(target):
                skipped += 1
                continue
            target = float(target)
            if not (math.isfinite(target) and 0.0 <= target <= 1.0):
                skipped += 1
                continue
            area = row['area']
            is_comp = bool(row.get('is_compensated', False))
            by_area.setdefault(area, []).append((fvec, target, is_comp))
        except Exception:
            skipped += 1
            continue

    for area, pairs in by_area.items():
        pairs = pairs[-max_rows_per_area:]
        _area_forecast_history.setdefault(area, []).extend(pairs)
        restored_total += len(pairs)
        es_like += sum(1 for _, y, _ in pairs if y >= 0.5)

    _area_total_samples = restored_total
    _area_es_samples = es_like
    return restored_total, skipped


def _bootstrap_area_history_from_csv_legacy_approx(log_dir: Path,
                                                     max_rows_per_area: int = 2000
                                                     ) -> Tuple[int, int]:
    """Ver13.12由来: edfs_area_history.csv(reach_rate等の観測記録のみ、
    特徴量ベクトル自体は保存していない)から _area_forecast_history を近似復元する。

    - CSV は読むだけで書き換えない（既存学習データを一切壊さない）
    - RAM 上の学習履歴だけを再構築する
    - FORECAST_CYCLES=2 サイクル後の実績と突き合わせ、教師データを作る
    - 特徴量は rate_change / acceleration など時系列依存分を 0 埋めした簡易版
      （最初の数サイクル走れば正規特徴量で自然に上書きされる）

    Rev14.4: この関数は「実際にRidgeへ入力された値」ではなく、reach_rateの
    時系列から事後的に近似計算した代理値であることに注意(この設計はVer13.12
    導入時からの既知の制約であり、Rev14.0での特徴量14次元化がこの性質を
    新たに生んだわけではない — 末尾2次元の0埋めが追加されただけ)。
    bootstrap_area_history_from_csv()は、edfs_ridge_training.csv(実値)が
    利用可能な場合はそちらを優先し、本関数は「実値ログがまだ蓄積されていない」
    場合のフォールバックとしてのみ使う。

    戻り値: (復元教師数, 学習成功エリア数)
    """
    global _area_forecast_history, _area_total_samples, _area_es_samples
    path = log_dir / EDFS_AREA_HISTORY_FILE
    if not path.exists():
        return 0, 0
    df = _robust_read_csv_dropping_tail(path, parse_dates=['timestamp'])
    if df is None:
        return 0, 0
    if df.empty or 'area' not in df.columns or 'reach_rate' not in df.columns:
        return 0, 0
    df = df.sort_values('timestamp').reset_index(drop=True)

    restored_total = 0
    es_like = 0
    for area in df['area'].unique():
        adf = df[df['area'] == area].reset_index(drop=True)
        if len(adf) < FORECAST_CYCLES + 1:
            continue
        rr = pd.to_numeric(adf['reach_rate'], errors='coerce').fillna(0.0).values
        rpi_series = pd.to_numeric(adf.get('rpi', pd.Series([np.nan]*len(adf))),
                                    errors='coerce').fillna(0.0).values / 100.0
        pairs = []
        for i in range(len(adf) - FORECAST_CYCLES):
            cur = float(rr[i])
            future = float(rr[i + FORECAST_CYCLES])
            rate_change = float(rr[i] - rr[i - 1]) if i >= 1 else 0.0
            acceleration = (float((rr[i] - rr[i - 1]) - (rr[i - 1] - rr[i - 2]))
                            if i >= 2 else 0.0)
            fvec = [
                cur, rate_change, acceleration,
                float(rpi_series[i]),
                0.0,
                clamp(cur, 0.0, 1.0),
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0,
            ]
            pairs.append((fvec, future, False))
        pairs = pairs[-max_rows_per_area:]
        _area_forecast_history.setdefault(area, []).extend(pairs)
        restored_total += len(pairs)
        es_like += sum(1 for _, y, _ in pairs if y >= 0.5)

    _area_total_samples = restored_total
    _area_es_samples = es_like
    return restored_total, 0


def bootstrap_area_history_from_csv(log_dir: Path,
                                    max_rows_per_area: int = 2000
                                    ) -> Tuple[int, int, str]:
    """起動時に永続化ログから _area_forecast_history を復元して即座に
    train_area_predictors() を走らせる。

    Rev14.4: 内部を2階層化した。
      1. edfs_ridge_training.csv(実際にRidgeへ入力された特徴量ベクトルそのもの、
         bootstrap_ridge_training_from_csv()参照)が存在し1件以上復元できれば、
         そちらを優先して使う(近似ではなく厳密な復元)。
      2. 存在しない、または0件の場合のみ、従来のedfs_area_history.csvベースの
         近似復元(_bootstrap_area_history_from_csv_legacy_approx())にフォールバック
         する(Rev14.3以前から運用中で、まだedfs_ridge_training.csvが
         蓄積されていない環境との後方互換のため)。
      二重計上を避けるため、両方を足し合わせることはしない(どちらか一方のみ採用)。
    Rev14.4で戻り値に3番目の要素(method: 'exact'|'legacy_approx')を追加した
    (呼び出し側で「今回はどちらの経路で復元したか」を表示できるようにするため。
    唯一の呼び出し元であるrun_cycle起動処理側も同時に更新済み)。

    戻り値: (復元教師数, 学習成功エリア数, method)
    """
    restored_total, _skipped = bootstrap_ridge_training_from_csv(log_dir, max_rows_per_area)
    used_exact = restored_total > 0
    if not used_exact:
        restored_total, _ = _bootstrap_area_history_from_csv_legacy_approx(log_dir, max_rows_per_area)

    trained = 0
    if _SKLEARN_AVAILABLE and restored_total >= 10:
        try:
            r2s = train_area_predictors()
            trained = len(r2s)
        except Exception:
            trained = 0
    return restored_total, trained, ('exact' if used_exact else 'legacy_approx')



def load_motion_training_data(log_dir: Path) -> pd.DataFrame:
    """edfs_centroid_history.csv から Motion 予測器の教師データを読み込む。

    ★重要（Ver13.13 設計契約）★
    重心履歴はマルチサンプリングで 1 検出サイクルにつき 5 行
    （horizon_min = 0/15/30/45/60）保存される。このうち実観測は
    horizon_min==0 の行だけで、残り 4 行は velocity による線形外挿の
    「解禁カウント押し上げ用」水増し行である。

    したがって Motion 予測器を学習・評価する際は、必ず本関数を通して
    horizon_min==0 の行だけを教師として使用すること。外挿行を教師に
    混ぜると、モデルが「自分が引いた直線」を再学習する自己参照になり、
    見かけの R² だけ上がって実際の移動予測精度が劣化する。

    後方互換: horizon_min 列を持たない旧 CSV は、全行を実観測(0)とみなす。

    戻り値: horizon_min==0 のみに絞った重心履歴 DataFrame
    """
    path = log_dir / EDFS_CENTROID_HISTORY_FILE
    if not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path)
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return df
    if 'horizon_min' in df.columns:
        hz = pd.to_numeric(df['horizon_min'], errors='coerce').fillna(0)
        df = df[hz == 0].copy()
    return df.reset_index(drop=True)


def count_real_centroid_observations(log_dir: Path) -> int:
    """Motion 解禁カウントの実態確認用: horizon_min==0 の実観測行数を返す。
    表示・デバッグ用途（解禁判定そのものは _edfs_centroid_history_count を使う）。
    """
    df = load_motion_training_data(log_dir)
    return int(len(df))



def generate_future_teachers(log_dir: Path) -> Tuple[pd.DataFrame, int]:
    """機能②③⑧: edfs_area_history.csvから未来RPI教師データを生成する。
    1時刻 = 12教師（全エリア展開）。
    戻り値: (teacher_df, valid_teacher_count)
    """
    path = log_dir / EDFS_AREA_HISTORY_FILE
    if not path.exists():
        return pd.DataFrame(), 0
    df = _robust_read_csv_dropping_tail(path, parse_dates=['timestamp'])
    if df is None:
        return pd.DataFrame(), 0
    df = df.sort_values('timestamp').reset_index(drop=True)
    result_parts = []
    for area in df['area'].unique():
        adf = df[df['area'] == area].copy().reset_index(drop=True)
        adf['area_rpi'] = pd.to_numeric(adf['rpi'], errors='coerce')
        for lag_min, shift in FUTURE_LAGS.items():
            col = f'future_rpi_{lag_min}'
            adf[col] = adf['area_rpi'].shift(-shift)
        for lag_min in FUTURE_LAGS:
            rcol = f'future_rpi_{lag_min}'
            ocol = f'future_open_{lag_min}'
            adf[ocol] = np.where(adf[rcol].isna(), np.nan,
                                 (adf[rcol] >= OPEN_RPI_THRESHOLD).astype(float))
        result_parts.append(adf)
    if not result_parts:
        return pd.DataFrame(), 0
    teacher_df = pd.concat(result_parts, ignore_index=True)
    valid_count = int(teacher_df['future_rpi_15'].notna().sum())
    return teacher_df, valid_count


def compute_area_centroid_features(centroid_lat: float, centroid_lon: float) -> Dict[str, Dict]:
    """機能⑥⑦: 各エリア代表座標からEs重心までの距離と方位を計算する。
    AREA_CENTERSを使用（BUILTIN_FIXED_MONITORのエリア別平均座標）。

    Rev14.0: bearing_from_centroid(重心からエリアを見た方位)を追加。
    Es重心の移動方位(es_area_direction)と比較することで、「重心がこの
    エリアへ向かっているか」の判定(geometric_eta_estimate等)に使う。
    既存キー(distance_to_centroid/bearing_to_centroid)は後方互換のため
    そのまま維持する。

    ※ この関数自体はVer11.3A(機能⑥⑦)で既に実装されていたが、Rev14.0
      以前は算出結果がどこからも参照されない死んだコードだった
      (run_cycle内で計算だけして代入先変数が以降一度も読まれていなかった)。
      Rev14.0でRidge特徴量(_area_feature_vector)・幾何学的ETA
      (geometric_eta_estimate)・表示パネルの3箇所へ配線し、実際に活用する。
    """
    result: Dict[str, Dict] = {}
    if pd.isna(centroid_lat):
        return result
    for area, (alat, alon) in AREA_CENTERS.items():
        dist = haversine_km(alat, alon, centroid_lat, centroid_lon)
        bear_to = bearing_deg(alat, alon, centroid_lat, centroid_lon)
        bear_from = bearing_deg(centroid_lat, centroid_lon, alat, alon)
        result[area] = {
            'distance_to_centroid': dist,
            'bearing_to_centroid': bear_to,
            'bearing_from_centroid': bear_from,
        }
    return result


def geometric_eta_estimate(distance_km: float, bearing_from_centroid_deg: float,
                            es_speed_kmh: float, es_direction_deg: float) -> Optional[float]:
    """Rev14.0: Ridge学習(統計モデル)に一切依存しない、幾何学的な到達時刻
    (分)の推定(提案②(c))。学習サンプルが0件の起動直後(コールドスタート)
    や、EDFS Level0の状態でも計算可能。

    Es重心の移動方位(es_direction_deg)と、重心から見た対象エリアの方位
    (bearing_from_centroid_deg)との角度差から、「重心が正味どれだけの
    速度成分でこのエリアへ近づいているか(射影速度)」を求め、
    距離/射影速度で到達までの概算時間を出す。

    あくまで瞬間速度ベクトルをそのまま延長する直線運動近似であり、
    Es雲の発生・消滅、加減速までは表現しない(粗い外挿である点に注意)。
    """
    if any(v is None or pd.isna(v) for v in
           (distance_km, bearing_from_centroid_deg, es_speed_kmh, es_direction_deg)):
        return None
    if es_speed_kmh <= 0:
        return None
    angle_diff = angular_diff_deg(es_direction_deg, bearing_from_centroid_deg)
    if angle_diff is None or pd.isna(angle_diff):
        return None
    projected_speed = es_speed_kmh * math.cos(math.radians(angle_diff))
    if projected_speed < GEO_ETA_MIN_PROJECTED_SPEED_KMH:
        return None
    eta_min = distance_km / projected_speed * 60.0
    if eta_min > GEO_ETA_MAX_MINUTES:
        return None
    return float(eta_min)


def geometric_eta_confidence(approach_cos: Optional[float]) -> str:
    """Rev14.1再指摘(#2): 幾何ETAの「値」と「信頼度」を分離する。

    速度ベクトル(es_area_speed/es_area_direction)自体は前サイクルとの2点
    差分から毎サイクル算出される瞬間値であり、「十分な履歴から得た速度」なのか
    「たまたま2点だけで偶然出た速度」なのかをETA値自体は区別しない。
    ここでは以下2軸の粗い代理指標からHIGH/MEDIUM/LOWを決める:
      - _edfs_centroid_history_count(重心履歴の蓄積量。EDFS_LV3_CENTROID_HIST
        =Motion Predictor解禁ラインを基準に、速度推定が安定してくる目安とする)
      - approach_cos(重心がそのエリアへどれだけ正面から向かっているか。
        1.0に近いほど射影速度の推定誤差の影響を受けにくい)
    あくまで参考表示用の粗い分類であり、Ridgeの選択・ゲーティングには使わない。
    """
    ch = _edfs_centroid_history_count
    if ch < 30:
        return 'LOW'
    if approach_cos is None or pd.isna(approach_cos):
        return 'LOW'
    ac = abs(float(approach_cos))
    if ch >= EDFS_LV3_CENTROID_HIST and ac >= 0.7:
        return 'HIGH'
    if ch >= 60 and ac >= 0.4:
        return 'MEDIUM'
    return 'LOW'


def compute_geometric_eta_all(area_centroid_feats: Dict[str, Dict],
                               spatial: Dict) -> Dict[str, Dict]:
    """Rev14.0: 全エリア分の幾何学的ETAをまとめて算出する。
    戻り値: {area: {'eta_min': float|None, 'approach_cos': float|None,
                     'eta_confidence': str}}
    approach_cos: 1.0に近いほど重心がそのエリアへ真っ直ぐ向かっている、
    -1.0に近いほど遠ざかっている(Ridge特徴量centroid_approach_cosと同じ定義)。
    eta_confidence: Rev14.1再指摘(#2)。ETA値自体の信頼性の粗い目安
    (HIGH/MEDIUM/LOW)。詳細はgeometric_eta_confidence()参照。
    """
    result: Dict[str, Dict] = {}
    es_speed = spatial.get('es_area_speed', np.nan)
    es_dir = spatial.get('es_area_direction', np.nan)
    for area, feats in area_centroid_feats.items():
        dist = feats.get('distance_to_centroid')
        bearing_from = feats.get('bearing_from_centroid')
        eta = geometric_eta_estimate(dist, bearing_from, es_speed, es_dir)
        approach_cos = None
        if not (bearing_from is None or pd.isna(bearing_from) or pd.isna(es_dir)):
            angle_diff = angular_diff_deg(es_dir, bearing_from)
            if angle_diff is not None and not pd.isna(angle_diff):
                approach_cos = float(math.cos(math.radians(angle_diff)))
        result[area] = {
            'eta_min': eta,
            'approach_cos': approach_cos,
            'eta_confidence': geometric_eta_confidence(approach_cos),
        }
    return result



def update_edfs_maturity() -> int:
    """機能⑩: EDFS成熟レベルを算出して更新する（毎サイクル自動判定）。
    Level0 Learning     → Level5 Dynamics Engine の6段階。
    """
    global _edfs_level, _edfs_level_label
    at  = _edfs_area_teachers
    ch  = _edfs_centroid_history_count
    ls  = _edfs_lifetime_samples
    ts  = _area_total_samples
    es  = _area_es_samples
    if at >= EDFS_LV5_AREA_TEACHERS and ch >= EDFS_LV5_CENTROID and ls >= EDFS_LV5_LIFETIME:
        lv, lb = 5, "Dynamics Engine"
    elif ls >= EDFS_LV4_LIFETIME:
        lv, lb = 4, "Lifetime Predictor"
    elif ch >= EDFS_LV3_CENTROID_HIST:
        lv, lb = 3, "Motion Predictor"
    elif at >= EDFS_LV2_AREA_TEACHERS:
        lv, lb = 2, "Area Predictor"
    elif ts >= EDFS_LV1_TOTAL and es >= EDFS_LV1_ES:
        lv, lb = 1, "Area Learning"
    else:
        lv, lb = 0, "Learning"
    _edfs_level = lv
    _edfs_level_label = lb
    return lv


def update_predictor_activation() -> None:
    """機能⑪: 自動予測器解禁 - しきい値到達で各予測器をONにする。"""
    at  = _edfs_area_teachers
    ch  = _edfs_centroid_history_count
    ls  = _edfs_lifetime_samples
    _edfs_predictor_active['area']     = at >= EDFS_LV2_AREA_TEACHERS
    _edfs_predictor_active['motion']   = ch >= EDFS_LV3_CENTROID_HIST
    _edfs_predictor_active['lifetime'] = ls >= EDFS_LV4_LIFETIME
    _edfs_predictor_active['dynamics'] = (
        at >= FEAT_DYNAMICS_AREA and
        ch >= FEAT_DYNAMICS_CENT and
        ls >= FEAT_DYNAMICS_LIFE
    )


def update_feature_activation() -> str:
    """機能⑫: 自動特徴量解禁 - データ量に応じてFeature Setを昇格させる。"""
    global _edfs_feature_set
    at = _edfs_area_teachers
    ch = _edfs_centroid_history_count
    ls = _edfs_lifetime_samples
    if at >= FEAT_DYNAMICS_AREA and ch >= FEAT_DYNAMICS_CENT and ls >= FEAT_DYNAMICS_LIFE:
        fs = "Dynamics"
    elif ch >= FEAT_ADV_MOTION_MIN:
        fs = "Advanced Motion"
    elif ch >= FEAT_MOTION_MIN:
        fs = "Motion"
    elif at >= FEAT_AREA_MIN:
        fs = "Area"
    else:
        fs = "Basic"
    _edfs_feature_set = fs
    return fs


def update_edfs_all(log_dir: Path, teacher_df: pd.DataFrame) -> None:
    """毎サイクル呼び出し: 教師数更新→成熟度→予測器→特徴量の順に自動更新。"""
    global _edfs_area_teachers, _edfs_lifetime_samples
    if not teacher_df.empty:
        _edfs_area_teachers = int(teacher_df['future_rpi_15'].notna().sum())
    _edfs_lifetime_samples = _area_es_samples // 5
    update_edfs_maturity()
    update_predictor_activation()
    update_feature_activation()



def _edfs_training_state_path(log_dir: Path) -> Path:
    return Path(log_dir) / EDFS_TRAINING_STATE_FILE


def load_training_state(log_dir: Path) -> None:
    """log_dir内の永続化ファイルから成熟モード状態を復元する（起動後1回だけ）。"""
    global _edfs_training_mode, _edfs_training_baseline, _edfs_drift_cooldown_until_cycle
    global _edfs_drift_events, _edfs_training_state_loaded, _edfs_drift_pending, _edfs_sfi_recent
    if _edfs_training_state_loaded:
        return
    _edfs_training_state_loaded = True
    p = _edfs_training_state_path(log_dir)
    try:
        if p.exists():
            data = json.loads(p.read_text(encoding='utf-8'))
            _edfs_training_mode = data.get('training_mode', 'normal')
            _edfs_training_baseline = data.get('baseline')
            _edfs_drift_cooldown_until_cycle = int(data.get('drift_cooldown_until_cycle', 0) or 0)
            _edfs_drift_events = data.get('drift_events', []) or []
            _edfs_drift_pending = data.get('drift_pending', {'location': 0, 'spread': 0, 'sfi': 0})
            _edfs_sfi_recent = data.get('sfi_recent', []) or []
    except Exception:
        pass


def save_training_state(log_dir: Path) -> None:
    """現在の成熟モード状態をlog_dirへ永続化する（再起動後も引き継ぐため）。"""
    state = {
        'training_mode': _edfs_training_mode,
        'baseline': _edfs_training_baseline,
        'drift_cooldown_until_cycle': _edfs_drift_cooldown_until_cycle,
        'drift_events': _edfs_drift_events[-10:],
        'drift_pending': _edfs_drift_pending,
        'sfi_recent': _edfs_sfi_recent[-EDFS_DRIFT_SFI_SMOOTH_N:],
        'saved_jst': now_jst().isoformat(),
    }
    try:
        p = _edfs_training_state_path(log_dir)
        p.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding='utf-8')
    except Exception:
        pass


def _edfs_stage12_done() -> bool:
    """Stage1(自動補正解禁)・Stage2(EDFS Lv2エリア別予測)が両方100%到達済みか。"""
    active_profile = globals().get('_v13_active_profile_data') or {}
    stage1_done = bool(active_profile.get('updated_jst')) or \
        int(active_profile.get('valid_rows_used', 0) or 0) > 0
    stage2_done = _edfs_area_teachers >= EDFS_LV2_AREA_TEACHERS
    return stage1_done and stage2_done


def _edfs_snapshot_stats() -> Optional[Dict]:
    """直近のEs重心分布(緯度/経度/広がり)の統計スナップショットを作る。
    サンプル不足時はNoneを返す（ドリフト判定を保留する）。
    """
    rows = [r for r in _es_spatial_history[-EDFS_DRIFT_SNAPSHOT_WINDOW:]
            if not pd.isna(r.get('es_area_center_lat', np.nan))
            and not pd.isna(r.get('es_area_center_lon', np.nan))]
    if len(rows) < EDFS_DRIFT_MIN_SAMPLES:
        return None
    lats = np.array([r['es_area_center_lat'] for r in rows], dtype=float)
    lons = np.array([r['es_area_center_lon'] for r in rows], dtype=float)
    spreads = np.array([r.get('es_area_spread', np.nan) for r in rows], dtype=float)
    spreads = spreads[np.isfinite(spreads)]
    sfi_val = np.nan
    try:
        sfi_val = float(fetch_realtime_solar_indices().get('sfi', np.nan))
    except Exception:
        pass
    return {
        'lat_mean': float(np.mean(lats)), 'lat_std': float(np.std(lats)) or 0.01,
        'lon_mean': float(np.mean(lons)), 'lon_std': float(np.std(lons)) or 0.01,
        'spread_mean': float(np.mean(spreads)) if len(spreads) else np.nan,
        'sfi': sfi_val,
        'n': len(rows),
        'obs_jst': now_jst().isoformat(),
    }


def _edfs_detect_drift(baseline: Dict, current: Dict) -> Dict[str, Optional[str]]:
    """ベースラインと現在のスナップショットを比較し、カテゴリ別
    ('location'/'spread'/'sfi')にドリフト理由の文字列（該当なしはNone）を返す。
    誤検知対策のため、この結果は単発では確定させず、呼び出し側で
    EDFS_DRIFT_CONFIRM_COUNT回の連続検知をもって初めて発動させること。
    """
    out: Dict[str, Optional[str]] = {'location': None, 'spread': None, 'sfi': None}
    try:
        blat, blon = baseline['lat_mean'], baseline['lon_mean']
        slat, slon = max(baseline['lat_std'], 0.01), max(baseline['lon_std'], 0.01)
        dlat = abs(current['lat_mean'] - blat) / slat
        dlon = abs(current['lon_mean'] - blon) / slon
        if dlat >= EDFS_DRIFT_LOCATION_SIGMA or dlon >= EDFS_DRIFT_LOCATION_SIGMA:
            out['location'] = (
                f"Es重心位置が変化 (緯度{dlat:.1f}σ/経度{dlon:.1f}σ, "
                f"しきい値{EDFS_DRIFT_LOCATION_SIGMA}σ)")
        bsp = baseline.get('spread_mean'); csp = current.get('spread_mean')
        if bsp and csp and not pd.isna(bsp) and not pd.isna(csp) and bsp > 0:
            ratio = csp / bsp
            if ratio >= EDFS_DRIFT_SPREAD_RATIO or ratio <= 1.0 / EDFS_DRIFT_SPREAD_RATIO:
                out['spread'] = (
                    f"Es分布の広がりが変化 (比率{ratio:.2f}倍, "
                    f"しきい値{EDFS_DRIFT_SPREAD_RATIO}倍)")
    except Exception:
        pass
    try:
        bsfi = baseline.get('sfi'); csfi = current.get('sfi')
        if bsfi and csfi and not pd.isna(bsfi) and not pd.isna(csfi) and bsfi > 0:
            diff = abs(csfi - bsfi)
            ratio = csfi / bsfi
            if diff >= EDFS_DRIFT_SFI_ABS_MIN and (
                    ratio >= EDFS_DRIFT_SFI_RATIO or ratio <= 1.0 / EDFS_DRIFT_SFI_RATIO):
                out['sfi'] = (
                    f"太陽活動(SFI移動平均)が変化 ({bsfi:.0f}→{csfi:.0f}sfu, "
                    f"比率{ratio:.2f}倍, しきい値{EDFS_DRIFT_SFI_RATIO}倍)")
    except Exception:
        pass
    return out


def check_maturity_and_drift(log_dir: Path, cur_cycle: int) -> int:
    """毎サイクル呼び出し: 成熟モードへの移行判定・ドリフト監視を行い、
    次回再学習までの間隔（サイクル数）を返す。

    Ver13.19: 誤検知対策のため、
      (1) SFIは点検のたびの単発値ではなく直近EDFS_DRIFT_SFI_SMOOTH_N回の
          移動平均で比較する（太陽自転由来の短期変動を吸収）
      (2) 各カテゴリ(location/spread/sfi)ごとに、EDFS_DRIFT_CONFIRM_COUNT回
          「連続」で検知されて初めて実際に見直し(通常モード復帰+再学習)を
          発動する。単発の外れ値だけでは動かない。
      (3) 連続が途切れたカテゴリはpendingカウントを0にリセットする。
    """
    global _edfs_training_mode, _edfs_training_baseline, _edfs_last_drift_check_cycle
    global _edfs_drift_cooldown_until_cycle, _edfs_drift_events
    global _edfs_drift_pending, _edfs_sfi_recent

    load_training_state(log_dir)

    if _edfs_training_mode == "normal":
        if _edfs_stage12_done() and cur_cycle >= _edfs_drift_cooldown_until_cycle:
            snap = _edfs_snapshot_stats()
            if snap is not None:
                _edfs_training_mode = "reduced"
                _edfs_training_baseline = snap
                _edfs_last_drift_check_cycle = cur_cycle
                _edfs_drift_pending = {'location': 0, 'spread': 0, 'sfi': 0}
                _edfs_sfi_recent = [snap['sfi']] if not pd.isna(snap.get('sfi', np.nan)) else []
                save_training_state(log_dir)
        return EDFS_NORMAL_TRAIN_INTERVAL_CYCLES

    if cur_cycle - _edfs_last_drift_check_cycle >= EDFS_DRIFT_CHECK_INTERVAL_CYCLES:
        _edfs_last_drift_check_cycle = cur_cycle
        snap = _edfs_snapshot_stats()
        if snap is not None and _edfs_training_baseline is not None:
            if not pd.isna(snap.get('sfi', np.nan)):
                _edfs_sfi_recent.append(snap['sfi'])
                _edfs_sfi_recent[:] = _edfs_sfi_recent[-EDFS_DRIFT_SFI_SMOOTH_N:]
            snap_for_compare = dict(snap)
            if _edfs_sfi_recent:
                snap_for_compare['sfi'] = float(np.mean(_edfs_sfi_recent))

            cat_reasons = _edfs_detect_drift(_edfs_training_baseline, snap_for_compare)

            for cat in ('location', 'spread', 'sfi'):
                if cat_reasons.get(cat):
                    _edfs_drift_pending[cat] = _edfs_drift_pending.get(cat, 0) + 1
                else:
                    _edfs_drift_pending[cat] = 0

            confirmed = [cat for cat, cnt in _edfs_drift_pending.items()
                         if cnt >= EDFS_DRIFT_CONFIRM_COUNT]
            if confirmed:
                reasons = [cat_reasons[cat] for cat in confirmed if cat_reasons.get(cat)]
                _edfs_training_mode = "normal"
                _edfs_drift_cooldown_until_cycle = cur_cycle + EDFS_DRIFT_COOLDOWN_CYCLES
                _edfs_drift_events.append({
                    'cycle': cur_cycle, 'obs_jst': now_jst().isoformat(),
                    'reasons': reasons,
                    'confirm_count': EDFS_DRIFT_CONFIRM_COUNT,
                })
                _edfs_drift_events[:] = _edfs_drift_events[-10:]
                _edfs_drift_pending = {'location': 0, 'spread': 0, 'sfi': 0}
                _edfs_sfi_recent = []
                try:
                    train_area_predictors()
                except Exception:
                    pass
                save_training_state(log_dir)
                return EDFS_NORMAL_TRAIN_INTERVAL_CYCLES
        save_training_state(log_dir)
    return EDFS_MATURE_TRAIN_INTERVAL_CYCLES



def print_edfs_status(android: bool = False) -> None:
    """EDFS Status表示（仕様書GUI追加表示: Level/Predictor/Statistics）。"""
    W = DISP_W_ANDROID if android else DISP_W_NORMAL
    lvstr = f"Level{_edfs_level} {_edfs_level_label}"
    if android:
        print(f"[EDFS] {lvstr}")
        active_preds = [k.capitalize() for k, v in _edfs_predictor_active.items() if v]
        if active_preds:
            print("  Active: " + " ".join(active_preds))
        print(f"  Teachers:{_edfs_area_teachers}  Centroid:{_edfs_centroid_history_count}")
        print(f"  {_edfs_feature_set} Feature Set")
        mode_s = "軽量(成熟)" if _edfs_training_mode == "reduced" else "通常"
        print(f"  学習モード:{mode_s}")
        return
    print("-" * W)
    print("EDFS Status")
    print(f"  EDFS Level  : {lvstr}")
    print(f"  Feature Set : {_edfs_feature_set}")
    label_map = {
        'area':     'Area Predictor',
        'motion':   'Motion Predictor',
        'lifetime': 'Lifetime Predictor',
        'dynamics': 'Dynamics Engine',
    }
    print("  Predictor Status:")
    for k, label in label_map.items():
        status = "\033[92mACTIVE\033[0m" if _edfs_predictor_active[k] else "待機中"
        needed = {
            'area':     f"Teachers>={EDFS_LV2_AREA_TEACHERS}",
            'motion':   f"Centroid>={EDFS_LV3_CENTROID_HIST}",
            'lifetime': f"Lifetime>={EDFS_LV4_LIFETIME}",
            'dynamics': f"Teachers>={FEAT_DYNAMICS_AREA}&Centroid>={FEAT_DYNAMICS_CENT}",
        }[k]
        suffix = "" if _edfs_predictor_active[k] else f"  ({needed})"
        print(f"    {label}: {status}{suffix}")
    print("  Statistics:")
    print(f"    Area Teachers    : {_edfs_area_teachers:>6}")
    print(f"    Centroid Samples : {_edfs_centroid_history_count:>6}")
    print(f"    Lifetime Samples : {_edfs_lifetime_samples:>6}")
    print(f"    Total Samples    : {_area_total_samples:>6}")
    print(f"    Es Samples       : {_area_es_samples:>6}")

    if _edfs_training_mode == "reduced":
        interval_s = f"{EDFS_MATURE_TRAIN_INTERVAL_CYCLES}サイクルに1回"
        mode_s = "\033[92m軽量学習(成熟モード)\033[0m"
    else:
        interval_s = f"{EDFS_NORMAL_TRAIN_INTERVAL_CYCLES}サイクルに1回"
        mode_s = "通常学習"
    print(f"  学習モード       : {mode_s}  (再学習: {interval_s})")
    if _edfs_training_mode == "reduced" and any(_edfs_drift_pending.values()):
        pend_s = " ".join(f"{k}:{v}/{EDFS_DRIFT_CONFIRM_COUNT}"
                           for k, v in _edfs_drift_pending.items() if v > 0)
        print(f"    ドリフト検知中(未確定): {pend_s}")
    if _edfs_drift_events:
        last_ev = _edfs_drift_events[-1]
        print(f"    直近ドリフト検知: {last_ev.get('obs_jst','--')}")
        for r in last_ev.get('reasons', []):
            print(f"      - {r}")

    print(f"  {SOURCE_VERSION} Causal Model:")
    sw_s = fmt(_latest_sw_features.get("kp_latest", np.nan), 1) if _latest_sw_features else "--"
    bz_s = fmt(_latest_sw_features.get("bz_gsm_delta_60m", np.nan), 2) if _latest_sw_features else "--"
    print(f"    Space Weather    : Kp={sw_s}  ΔBz(60m)={bz_s}nT")
    for st in DISPLAY_ORDER:
        eqi_d = _latest_eqi_by_station.get(st, {})
        eqi_v = eqi_d.get("eqi", np.nan)
        phase_v = _latest_es_phase_by_station.get(st, "---")
        conf_v = eqi_d.get("eqi_confidence", np.nan)
        print(f"    {STATIONS[st]['jp']:<4}  EQI={fmt(eqi_v,1):>5}  Phase={phase_v:<10}  EQI信頼={fmt(conf_v,2)}")
    print("-" * W)


def ver113_monitor_line() -> str:
    """Ver11.3モニター表示（EDFS Levelに更新）。"""
    if _edfs_level >= 5:
        return f"EDFS Level5 Dynamics Engine ACTIVE"
    return f"EDFS Level{_edfs_level} {_edfs_level_label}  Teachers:{_edfs_area_teachers}"


def reach_trend_symbol(area: str) -> str:
    """到達率履歴からトレンド記号(▲/→/▼)のみを返す。
    数値(30分後○%)は出さない。今計算しているのは到達率の傾きであり、
    「傾きの外挿 = 30分後予測」は説明可能性の観点で不正直（タカユキ指摘）。
    エリア別学習器が育つまでは傾き方向のみを報告する。
    """
    hist = _area_reach_history.get(area, [])
    if len(hist) < 2:
        return '→'
    t0 = hist[0][0]
    times = np.array([(t - t0).total_seconds() / 60.0 for t, _ in hist])
    vals = np.array([v for _, v in hist])
    if len(hist) >= 3 and (times[-1] - times[0]) > 0:
        slope = float(np.polyfit(times, vals, 1)[0])
    else:
        dt_min = max(times[-1] - times[-2], 1.0)
        slope = float((vals[-1] - vals[-2]) / dt_min)
    if slope > AREA_REACH_TREND_EPS:
        return '▲'
    if slope < -AREA_REACH_TREND_EPS:
        return '▼'
    return '→'


def update_psk_freshness(success: bool) -> None:
    global _last_psk_success_dt
    if success:
        _last_psk_success_dt = now_jst()


def psk_freshness_minutes() -> float:
    if _last_psk_success_dt is None:
        return float('inf')
    return (now_jst() - _last_psk_success_dt).total_seconds() / 60.0


def ft8_confidence_label(fresh_min: float, fixed_db_empty: bool, rpi: float, area_df: Optional[pd.DataFrame] = None) -> str:
    """信頼度(LOW/MEDIUM/HIGH)。Ver11.0J: 固定局5局未満(area_confidence=='LOW')の
    エリアがRPI算出に混ざっている場合はHIGHを名乗らずMEDIUMに抑える
    （局数の薄いエリアの到達率は統計的にブレやすいため）。
    """
    if fresh_min > FT8_STALE_MIN_THRESHOLD or pd.isna(rpi):
        return "LOW"
    if fixed_db_empty:
        return "MEDIUM"
    if area_df is not None and not area_df.empty and "area_confidence" in area_df.columns:
        valid = area_df.dropna(subset=["reach_ratio"])
        if (valid["area_confidence"] == "LOW").any():
            return "MEDIUM"
    return "HIGH"


def build_mock_psk_spots(step: int = 0) -> pd.DataFrame:
    """--mock時のPSKReporter応答合成データ。ネットワーク無しでVer11.0J表示の
    全パスを検証できるようにする（既存のbuild_mock_latest_sample()と同じ思想）。
    """
    base_dt = now_jst().replace(second=0, microsecond=0)
    rows = [
        ("JA1ABC", "PM95dq", -8), ("JA1XYZ", "QM05hk", -12),
        ("JA7DEF", "QN02ab", -15), ("8N8TOK", "QN03cd", -18),
        ("JR6GHI", "PL36ab", -6 + step), ("JS6JKL", "PL36cd", -10),
        ("JA6MNO", "PL52ab", -14), ("JA3PQR", "PM74cd", -20),
        ("BV1STU", "PL05ef", -22 - step),
    ]
    recs = []
    for cs, grid, snr in rows:
        recs.append({
            "receiver_callsign": cs, "receiver_locator": grid,
            "sender_locator": "QN03ks",
            "freq_hz": 28074000.0, "snr": float(snr),
            "obs_jst": base_dt.isoformat(),
        })
    return pd.DataFrame(recs)



def build_mstid_nict_url(dt_utc: datetime, kind: str = MSTID_NICT_KIND) -> str:
    """Ver13.14: NICT GEONET リアルタイムGPS全電子数マップの直リンクURLを構築する。
       dt_utc は UT(協定世界時)の datetime（分は10分単位に切り捨て済みであること）。
       kind は MSTID_NICT_KIND_TABLE のキー ('A'/'H60'/'H15'/'R')。"""
    folder, stem_word = MSTID_NICT_KIND_TABLE.get(kind, MSTID_NICT_KIND_TABLE['H15'])
    y = dt_utc.year
    doy = dt_utc.timetuple().tm_yday
    hh = f"{dt_utc.hour:02d}"
    mm = f"{dt_utc.minute:02d}"
    mmdd = f"{dt_utc.month:02d}{dt_utc.day:02d}"
    fname = f"{hh}_{mm}_{mmdd}_{y}_{stem_word}.jpg"
    return f"{MSTID_NICT_BASE_URL}/{folder}/{y}/{doy:03d}/{fname}"


def _mstid_nict_latest_slot_utc(now_utc: Optional[datetime] = None,
                                 publish_lag_min: int = MSTID_NICT_PUBLISH_LAG_MIN) -> datetime:
    """直近で公開済みと見なせる10分枠(UT)を返す。配信ラグ分だけ遡ってから
       10分単位に切り捨てることで、まだ生成されていない枠を掴みにいくのを防ぐ。"""
    now_utc = now_utc or datetime.now(timezone.utc)
    target = now_utc - timedelta(minutes=publish_lag_min)
    floored_minute = (target.minute // 10) * 10
    return target.replace(minute=floored_minute, second=0, microsecond=0, tzinfo=None)


def list_mstid_images(image_dir: Path) -> List[Path]:
    if not image_dir.exists():
        return []
    files = []
    for ext in ['*.png','*.jpg','*.jpeg','*.bmp','*.webp']:
        files.extend(image_dir.glob(ext))
    return sorted(files, key=lambda p: p.stat().st_mtime)


_mstid_fetch_fail_count = 0
_mstid_fetch_last_error = None


def fetch_mstid_image_once(save_dir: Path,
                            keep_latest_n: int = MSTID_KEEP_LATEST_N,
                            timeout: int = HTTP_TIMEOUT_SEC,
                            kind: Optional[str] = None) -> bool:
    """Ver13.14: MSTID画像を1枚取得し save_dir へ保存する。
       優先順位: 1) MSTID_AUTO_FETCH_URL（--mstid-fetch-url で固定URL指定時）
                 2) MSTID_NICT_AUTO_ENABLED が真なら NICT GEONET の直近10分枠を自動生成して取得
                 3) いずれも不可なら何もせず False（--mstid-image-dir への手動配置にフォールバック）
       - kind 省略時は現在の MSTID_NICT_KIND（--mstid-source）をそのつど参照する
         （関数定義時の既定値に固定されないよう、呼び出し時に global を読む）。
       - 同じ10分枠は決定論的なファイル名にすることで重複ダウンロードを避ける
         （既に保存済みならリクエストせず True を返す）。
       - 取得成功時は save_dir 配下に mtime 昇順で keep_latest_n 枚だけ残し、
         古いファイルは自動削除する（無制限に増え続けるのを防止）。
       - Ver13.14修正: 失敗は握りつぶして False を返す点は変わらないが（1枚の取得
         失敗でサイクル自体を止めない）、以前は完全に無音だったため NICT 側が
         URL形式を変えた等の不調に気付く手段が無かった。初回失敗・エラー内容が
         変化した時・一定回数ごと、はコンソールへ警告を出すようにした。"""
    global _mstid_fetch_fail_count, _mstid_fetch_last_error
    kind = kind or MSTID_NICT_KIND
    url = MSTID_AUTO_FETCH_URL
    slot_tag = None
    if not url:
        if not MSTID_NICT_AUTO_ENABLED:
            return False
        try:
            slot = _mstid_nict_latest_slot_utc()
            url = build_mstid_nict_url(slot, kind=kind)
            slot_tag = slot.strftime('%Y%m%d_%H%M')
        except Exception:
            return False
    try:
        save_dir.mkdir(parents=True, exist_ok=True)
        if slot_tag:
            fname_guess = save_dir / f"mstid_nict_{kind}_{slot_tag}.jpg"
            if fname_guess.exists():
                _mstid_fetch_fail_count = 0
                _mstid_fetch_last_error = None
                return True
        headers = {"User-Agent": "EDFS-Ver13.14-Pydroid3-MSTID/1.0"}
        r = requests.get(url, timeout=timeout, headers=headers)
        r.raise_for_status()
        ctype = (r.headers.get('Content-Type') or '').lower()
        ext = '.png'
        if 'jpeg' in ctype or 'jpg' in ctype:
            ext = '.jpg'
        elif 'gif' in ctype:
            ext = '.gif'
        elif 'bmp' in ctype:
            ext = '.bmp'
        elif 'webp' in ctype:
            ext = '.webp'
        if slot_tag:
            fname = f"mstid_nict_{kind}_{slot_tag}{ext}"
        else:
            fname = f"mstid_{time.strftime('%Y%m%d_%H%M%S')}{ext}"
        (save_dir / fname).write_bytes(r.content)
        files = list_mstid_images(save_dir)
        if len(files) > keep_latest_n:
            for old in files[:-keep_latest_n]:
                try:
                    old.unlink()
                except Exception:
                    pass
        _mstid_fetch_fail_count = 0
        _mstid_fetch_last_error = None
        return True
    except Exception as e:
        _mstid_fetch_fail_count += 1
        err_msg = f"{type(e).__name__}: {e}"
        should_log = (
            _mstid_fetch_fail_count == 1
            or err_msg != _mstid_fetch_last_error
            or _mstid_fetch_fail_count in (5, 20)
            or _mstid_fetch_fail_count % 50 == 0
        )
        if should_log:
            print(f"⚠ MSTID画像取得失敗 (連続{_mstid_fetch_fail_count}回目・kind={kind}): {err_msg}")
            print(f"  URL: {url}")
        _mstid_fetch_last_error = err_msg
        return False


def _load_image_gray_small(path: Path, size=(256,256),
                            crop_frac: Optional[Tuple[float,float,float,float]]=None) -> np.ndarray:
    """Ver13.14: crop_frac=(left,top,right,bottom) を指定すると、画像サイズに
       対する比率でその範囲だけを切り出してから縮小する。タイトル文字・軸ラベル・
       カラーバー・凡例枠など、常に同じ位置に出る「地図の飾り」を除外し、
       エッジ検出/PCAが実データ以外に引っ張られるのを防ぐための前処理。"""
    img = Image.open(path).convert('L')
    if crop_frac:
        w, h = img.size
        l, t, r, b = crop_frac
        box = (int(w*l), int(h*t), int(w*r), int(h*b))
        if box[2] > box[0] and box[3] > box[1]:
            img = img.crop(box)
    img = img.resize(size)
    img = img.filter(ImageFilter.GaussianBlur(radius=1.0))
    return np.asarray(img, dtype=np.float32)/255.0


def _extract_edge_points(arr: np.ndarray, q: float=0.92) -> np.ndarray:
    gy, gx = np.gradient(arr)
    mag = np.sqrt(gx*gx + gy*gy)
    thr = np.quantile(mag, q)
    ys, xs = np.where(mag >= thr)
    if len(xs) == 0:
        return np.empty((0,2), dtype=np.float32)
    return np.column_stack([xs.astype(np.float32), ys.astype(np.float32)])


def _pca_orientation_deg(pts: np.ndarray, sx: float = 1.0, sy: float = 1.0) -> float:
    """主成分方向を返す。sx/sy はそれぞれ x(経度)・y(緯度)方向の実距離換算係数
       (km/px 等、相対比のみ意味を持つ)。既定 sx=sy=1.0 なら従来どおり
       ピクセル空間そのままの角度(等方仮定)。Ver13.14: 経度側が緯度側より
       km/px が小さい(=同じ画素数でも実距離が短い)場合、この補正をしないと
       方向が最大で ~atan(sx/sy) 分だけ実際の角度からズレる。"""
    if pts.shape[0] < 5:
        return np.nan
    m = pts.mean(axis=0)
    a = pts - m
    cov = np.cov(a.T)
    vals, vecs = np.linalg.eigh(cov)
    v = vecs[:, np.argmax(vals)]
    vx_real, vy_real = v[0] * sx, v[1] * sy
    ang = math.degrees(math.atan2(vy_real, vx_real))
    if ang < 0:
        ang += 180.0
    return ang


def _centroid_of_strength(arr: np.ndarray) -> Tuple[float,float,float]:
    gy, gx = np.gradient(arr)
    mag = np.sqrt(gx*gx + gy*gy)
    s = mag
    total = float(s.sum())
    if total <= 0:
        return np.nan, np.nan, 0.0
    ys, xs = np.indices(arr.shape)
    cx = float((xs*s).sum()/total)
    cy = float((ys*s).sum()/total)
    strength = float(np.quantile(mag, 0.98))
    return cx, cy, strength


def angular_diff_deg(a: float, b: float) -> float:
    if pd.isna(a) or pd.isna(b):
        return np.nan
    return (a - b + 180.0) % 360.0 - 180.0


def assumed_es_axis_screen_deg() -> float:
    """Rev14.10改名(#重大な問題③): 旧名 estimate_es_screen_vector() は「推定」を
    謳っていたが、実態は画面上のEs軸方向を45°固定で仮定しているだけで、実測でも
    推定でもない。呼び出し側(compute_approach_index)の結果が実測ベースの指標だと
    誤解されないよう、関数名と戻り値の性質を一致させた。
    値そのもの(45.0固定)はRev14.9までと不変(挙動不変の改名のみ)。
    将来、NICT複数局やFT8観測から実際のEs軸方位を求められるようになった時点で、
    この固定値をそちらの実測値に置き換えることを想定している。
    """
    return 45.0


def analyze_mstid_images(image_dir: Path, km_per_px: float=8.0, frame_min: float=15.0,
                         origin_xy: Optional[Tuple[float,float]]=None,
                         km_per_px_lat: Optional[float]=None,
                         km_per_px_lon: Optional[float]=None,
                         crop_frac: Optional[Tuple[float,float,float,float]]=MSTID_CROP_FRAC) -> Dict[str,float]:
    """Ver13.14: 東西(経度)/南北(緯度)の異方性校正に対応。
       km_per_px_lat/km_per_px_lon を指定すると、方向(PCA)・距離・速度の
       いずれもその軸別スケールで計算する（NICT GEONETの等経緯度マップは
       経度1°の実距離が緯度1°より短いため、単一スカラーだと方向が
       最大 ~20% 分ズレていた）。どちらも未指定なら従来どおり km_per_px を
       等方(縦横共通)に使う（後方互換）。
       crop_frac は (left,top,right,bottom) の比率タプルで、タイトル・軸ラベル・
       カラーバー・凡例枠を除外してからエッジ検出/PCAにかける。既定は
       MSTID_CROP_FRAC（NICT画像向けの暫定値）。手動配置した別形式の画像等で
       クロップ自体が不適切な場合は crop_frac=None で無効化できる。"""
    sx = km_per_px_lon if km_per_px_lon is not None else km_per_px
    sy = km_per_px_lat if km_per_px_lat is not None else km_per_px
    files = list_mstid_images(image_dir)
    out = {'mstid_dir_deg': np.nan, 'mstid_speed_kmh': np.nan, 'mstid_dist_km': np.nan,
           'mstid_strength': np.nan, 'mstid_cx_px': np.nan, 'mstid_cy_px': np.nan}
    if len(files) == 0:
        return out
    arr = _load_image_gray_small(files[-1], crop_frac=crop_frac)
    pts = _extract_edge_points(arr)
    out['mstid_dir_deg'] = _pca_orientation_deg(pts, sx=sx, sy=sy)
    cx, cy, strength = _centroid_of_strength(arr)
    out['mstid_cx_px'] = cx; out['mstid_cy_px'] = cy; out['mstid_strength'] = strength
    ox, oy = origin_xy if origin_xy else (arr.shape[1]/2.0, arr.shape[0]/2.0)
    if not pd.isna(cx):
        out['mstid_dist_km'] = math.hypot((cx-ox)*sx, (cy-oy)*sy)
    if len(files) >= 2:
        arr2 = _load_image_gray_small(files[-2], crop_frac=crop_frac); cx2, cy2, _ = _centroid_of_strength(arr2)
        if not (pd.isna(cx) or pd.isna(cx2)):
            dist_km = math.hypot((cx-cx2)*sx, (cy-cy2)*sy)
            out['mstid_speed_kmh'] = dist_km / max(frame_min/60.0, 1e-6)
    return out


def compute_approach_index(mstid: Dict[str,float], es_axis_dir_deg: float=45.0) -> float:
    """MSTIDの画面上の移動方向と「仮定したEs軸方向(es_axis_dir_deg、既定は
    assumed_es_axis_screen_deg()の固定45°)」との整合度を返す。

    Rev14.10注記(#重大な問題③): 内部識別子(mstid_approach_index等)は既存の
    CSVスキーマ・Ridge特徴量名との互換性維持のため変更していないが、これは
    「MSTIDとEsの実測ベースの接近指数」ではない。es_axis_dir_degが実測ではなく
    固定45°の仮定値である以上、本関数の出力は「MSTIDの画面方向と、仮定した
    Es軸との整合度」以上の意味を持たない。UI表示やEvidence説明でこの値を
    見せる際は過大評価した名称・説明を避けること。
    """
    d = mstid.get('mstid_dir_deg', np.nan); v = mstid.get('mstid_speed_kmh', np.nan); dist = mstid.get('mstid_dist_km', np.nan); strength = mstid.get('mstid_strength', np.nan)
    if any(pd.isna(x) for x in [d,v,dist]):
        return np.nan
    diff = angular_diff_deg(d, es_axis_dir_deg)
    proj = math.cos(math.radians(diff))
    s = 1.0 if pd.isna(strength) else max(0.1, float(strength)/0.2)
    return proj * float(v) * s / max(float(dist), 1.0)



def empirical_probability(fx_for_prob, trend, cscore, ns_score, shear_feats, vstruct, gw_feats, ef_feats, mstid_feats, cycles: int, profile: Dict, ft8_feats: Optional[Dict] = None) -> float:
    if pd.isna(fx_for_prob):
        return np.nan
    z = float(profile.get('prob_intercept', -2.0))
    z += float(profile.get('w_fxes', 0.32)) * (float(fx_for_prob)-6.0)
    z += float(profile.get('w_trend', 0.12)) * clamp(trend if not pd.isna(trend) else 0.0, -4, 4)
    z += float(profile.get('w_cluster', 0.75)) * (0.0 if pd.isna(cscore) else cscore)
    z += float(profile.get('w_ns', 0.25)) * (0.0 if pd.isna(ns_score) else ns_score)
    z += float(profile.get('w_shear', 0.45)) * (0.0 if pd.isna(shear_feats.get('convergence_index')) else min(1.5, shear_feats['convergence_index']))
    z += float(profile.get('w_layer', 0.25)) * (0.0 if pd.isna(vstruct.get('blanketing_index')) else vstruct['blanketing_index'])
    z += float(profile.get('w_descent', 0.20)) * (0.0 if pd.isna(vstruct.get('height_descent_rate')) else min(2.0, max(0.0, vstruct['height_descent_rate']/10.0)))
    z -= float(profile.get('w_gw_damp', 0.20)) * (0.0 if pd.isna(gw_feats.get('forecast_damping')) else min(2.0, gw_feats['forecast_damping']))
    z -= float(profile.get('w_ef_penalty', 0.15)) * (0.0 if pd.isna(ef_feats.get('masking_risk_index')) else ef_feats['masking_risk_index'])
    if not pd.isna(mstid_feats.get('mstid_approach_index', np.nan)):
        z += float(profile.get('w_mstid_approach', 0.0)) * float(mstid_feats.get('mstid_approach_index', 0.0))
    if not pd.isna(trend) and trend > 0.8:
        z += float(profile.get('w_rising', 0.20))
    w_fof2 = float(profile.get('w_fof2', 0.0))
    if w_fof2 != 0.0:
        fes = ef_feats.get('foF2_es_ratio', np.nan)
        if pd.isna(fes):
            fes = ef_feats.get('foF2_ratio', np.nan)
        z_fes = std_z(fes, 'foF2_es_ratio', profile)
        if not pd.isna(z_fes):
            sign_fof2 = float(profile.get('sign_fof2', 0.0))
            z += w_fof2 * sign_fof2 * z_fes
    w_dfof2 = float(profile.get('w_dfof2', 0.0))
    if w_dfof2 != 0.0:
        z_dfof2 = std_z(ef_feats.get('dfoF2', np.nan), 'dfoF2', profile)
        if not pd.isna(z_dfof2):
            sign_dfof2 = float(profile.get('sign_dfof2', 0.0))
            z += w_dfof2 * sign_dfof2 * z_dfof2
    w_hmf2 = float(profile.get('w_hmf2', 0.0))
    if w_hmf2 != 0.0:
        z_hmf2 = std_z(ef_feats.get('hmF2', np.nan), 'hmF2', profile)
        if not pd.isna(z_hmf2):
            sign_hmf2 = float(profile.get('sign_hmf2', 0.0))
            z += w_hmf2 * sign_hmf2 * z_hmf2
    w_dhmf2 = float(profile.get('w_dhmf2', 0.0))
    if w_dhmf2 != 0.0:
        z_dhmf2 = std_z(ef_feats.get('dhmF2', np.nan), 'dhmF2', profile)
        if not pd.isna(z_dhmf2):
            sign_dhmf2 = float(profile.get('sign_dhmf2', 0.0))
            z += w_dhmf2 * sign_dhmf2 * z_dhmf2
    w_m3000f2 = float(profile.get('w_m3000f2', 0.0))
    if w_m3000f2 != 0.0:
        z_m3000f2 = std_z(ef_feats.get('m3000f2_norm', np.nan), 'm3000f2_norm', profile)
        if not pd.isna(z_m3000f2):
            sign_m3000f2 = float(profile.get('sign_m3000f2', 0.0))
            z += w_m3000f2 * sign_m3000f2 * z_m3000f2
    ft8_feats = ft8_feats or {}
    w_rpi = float(profile.get('w_rpi', 0.0))
    if w_rpi != 0.0:
        z_rpi = std_z(ft8_feats.get('rpi_norm', np.nan), 'rpi_norm', profile)
        if not pd.isna(z_rpi):
            sign_rpi = float(profile.get('sign_rpi', 0.0))
            z += w_rpi * sign_rpi * z_rpi
    w_es_mode = float(profile.get('w_es_mode', 0.0))
    if w_es_mode != 0.0:
        z_es_mode = std_z(ft8_feats.get('es_mode_score', np.nan), 'es_mode_score', profile)
        if not pd.isna(z_es_mode):
            sign_es_mode = float(profile.get('sign_es_mode', 0.0))
            z += w_es_mode * sign_es_mode * z_es_mode
    w_eqi = float(profile.get('w_eqi', 0.0))
    if w_eqi != 0.0:
        z_eqi = std_z(ft8_feats.get('eqi_norm', np.nan), 'eqi_norm', profile)
        if not pd.isna(z_eqi):
            sign_eqi = float(profile.get('sign_eqi', 0.0))
            z += w_eqi * sign_eqi * z_eqi
    w_sw_score = float(profile.get('w_sw_score', 0.0))
    if w_sw_score != 0.0:
        z_sw = std_z(ft8_feats.get('sw_score', np.nan), 'sw_score', profile)
        if not pd.isna(z_sw):
            sign_sw_score = float(profile.get('sign_sw_score', 0.0))
            z += w_sw_score * sign_sw_score * z_sw
    w_es_phase = float(profile.get('w_es_phase', 0.0))
    if w_es_phase != 0.0:
        z_phase = std_z(ft8_feats.get('es_phase_num', np.nan), 'es_phase_num', profile)
        if not pd.isna(z_phase):
            sign_es_phase = float(profile.get('sign_es_phase', 0.0))
            z += w_es_phase * sign_es_phase * z_phase
    for ia_col, ia_wkey in [
        ("ia_eqi_fxes",  "w_ia_eqi_fxes"),
        ("ia_eqi_mstid", "w_ia_eqi_mstid"),
        ("ia_sw_fxes",   "w_ia_sw_fxes"),
        ("ia_sw_mstid",  "w_ia_sw_mstid"),
        ("ia_eqi_sw",    "w_ia_eqi_sw"),
    ]:
        w_ia = float(profile.get(ia_wkey, 0.0))
        if w_ia != 0.0:
            z_ia = std_z(ft8_feats.get(ia_col, np.nan), ia_col, profile)
            if not pd.isna(z_ia):
                sign_ia = float(profile.get("sign_" + ia_wkey[2:], 0.0))
                z += w_ia * sign_ia * z_ia
    p = sigmoid(z) * 100.0
    if cycles < 3:
        p = min(p, 50.0)
    elif cycles < MIN_ANALYSIS_POINTS:
        p = min(p, 85.0)
    return clamp(p, 1.0, PROB_MAX_PERCENT)


def state_from_prob(p, cycles: int, profile: Dict) -> str:
    if cycles < 3:
        return '暫定'
    if pd.isna(p):
        return '判定不可'
    if p >= float(profile.get('strong_open_thr', 75.0)):
        return '強い開通'
    if p >= float(profile.get('open_thr', 55.0)):
        return '開通'
    if p >= float(profile.get('watch_thr', 35.0)):
        return '注意監視'
    if p >= float(profile.get('weak_thr', 20.0)):
        return '弱い兆候'
    return '低調'


def estimate_network_motion(dyn: pd.DataFrame) -> Dict[str,float]:
    """Ver10.8J: 重心緯度(es_centroid_lat)とその変化量(d_centroid_lat)を追加公開。
    沖縄/山川/国分寺/稚内の4点では隣接局差分による勾配(lat_gradient)はノイズが
    大きいため、fxEs加重重心方式を採用（MSTID画像が欠損していても移動推定可能）。
    """
    nan_result = {'direction':'不明','speed_kmh':np.nan,'confidence':0.0,
                  'es_centroid_lat':np.nan,'d_centroid_lat':np.nan}
    if dyn.empty:
        return nan_result
    df = dyn.dropna(subset=['obs_dt','fxEs']).copy()
    if df.empty:
        return nan_result
    df['lat'] = df['station'].map(lambda k: STATIONS[k]['lat'])
    df['w'] = np.clip(df['fxEs'] - 3.0, 0.0, None)
    centers = []
    for t, g in df.groupby('obs_dt'):
        if len(g) < 3 or g['w'].sum() <= 0:
            continue
        center_lat = float(np.average(g['lat'], weights=g['w']))
        centers.append((t, center_lat, len(g), float(g['w'].sum())))
    if len(centers) < 2:
        return nan_result
    cdf = pd.DataFrame(centers, columns=['obs_dt','center_lat','n','w_sum']).sort_values('obs_dt')
    recent = cdf.tail(4)
    es_centroid_lat = float(recent.iloc[-1]['center_lat'])
    dt_h = (recent.iloc[-1]['obs_dt'] - recent.iloc[0]['obs_dt']).total_seconds()/3600.0
    if dt_h <= 0:
        return {'direction':'不明','speed_kmh':np.nan,'confidence':0.0,
                'es_centroid_lat':es_centroid_lat,'d_centroid_lat':np.nan}
    dlat = float(recent.iloc[-1]['center_lat'] - recent.iloc[0]['center_lat'])
    speed_kmh = abs(dlat)*111.0/dt_h
    if speed_kmh < MIN_SPEED_KMH:
        direction = '停滞'
    elif dlat > 0:
        direction = '北上'
    else:
        direction = '南下'
    if speed_kmh > MAX_SPEED_KMH:
        speed_kmh = np.nan; direction='不安定'
    confidence = clamp(len(recent)/4.0 * min(1.0, float(recent['n'].mean())/4.0), 0.0, 1.0)
    return {'direction':direction,'speed_kmh':speed_kmh,'confidence':confidence,
            'es_centroid_lat':es_centroid_lat,'d_centroid_lat':dlat}



def calculate_analysis_rate(dyn: pd.DataFrame):
    if dyn.empty:
        return 0.0, 0, {}
    counts = dyn.groupby('station')['obs_dt'].nunique().to_dict()
    rates = {st: clamp(counts.get(st,0)/MIN_ANALYSIS_POINTS*100.0, 0.0, 100.0) for st in STATIONS}
    overall = float(np.mean(list(rates.values()))) if rates else 0.0
    cycles = int(max(counts.values())) if counts else 0
    return overall, cycles, rates


def forecast_30m(g: pd.DataFrame, stage: int, shear_feats, vstruct, gw_feats, ef_feats, cycles: int):
    if g.empty or cycles < 3:
        return np.nan, '未収束'
    latest = float(g.sort_values('obs_dt').iloc[-1]['fxEs'])
    trend_raw = g.sort_values('obs_dt').iloc[-1].get('d_eff', np.nan)
    trend = 0.0 if pd.isna(trend_raw) else float(trend_raw)
    base = latest + trend*0.5
    conv = 0.0 if pd.isna(shear_feats.get('convergence_index')) else float(shear_feats['convergence_index'])
    descent = 0.0 if pd.isna(vstruct.get('height_descent_rate')) else float(vstruct['height_descent_rate'])
    damping = 0.0 if pd.isna(gw_feats.get('forecast_damping')) else float(gw_feats['forecast_damping'])
    ef_penalty = 0.0 if pd.isna(ef_feats.get('masking_risk_index')) else float(ef_feats['masking_risk_index'])
    pred = base + 0.4*conv + 0.2*max(0.0, descent/5.0) - 0.5*damping - 0.3*ef_penalty
    if stage >= 3 and len(g) >= 8:
        d = g['d_inst'].dropna().tail(5)
        acc = float(d.diff().mean()/0.25) if len(d) >= 3 else 0.0
        pred += 0.15*acc*0.5*0.5
        model = '準物理+加速度'
    else:
        model = '準物理'
    return clamp(pred, 0.0, 25.0), model


def analyze(history: pd.DataFrame, stage: int, wind_file: Path, profile: Dict, force_no_wind: bool=False, mstid_feats: Optional[Dict[str,float]]=None, ft8_feats: Optional[Dict[str,float]]=None):
    dyn = add_dynamics(history)
    sample = latest_network_sample(dyn)
    ns_score, sign_changes = mstid_sign_change_score(sample)
    ns_phrase = north_south_phrase(sign_changes)
    analysis_rate, cycles, station_rates = calculate_analysis_rate(dyn)
    motion = estimate_network_motion(dyn)
    mstid_feats = dict(mstid_feats or {})
    mstid_feats.setdefault('mstid_approach_index', np.nan)
    ft8_feats = dict(ft8_feats or {})
    ft8_feats.setdefault('rpi_norm', np.nan)
    ft8_feats.setdefault('es_mode_score', np.nan)

    sw_score_net = ft8_feats.get("sw_score", np.nan)
    eqi_vals = [_latest_eqi_by_station.get(st, {}).get("eqi", np.nan) for st in DISPLAY_ORDER]
    eqi_valid = [v for v in eqi_vals if not pd.isna(v)]
    eqi_repr = float(np.max(eqi_valid)) if eqi_valid else np.nan
    eqi_norm_repr = clamp(eqi_repr / 100.0, 0.0, 1.0) if not pd.isna(eqi_repr) else np.nan

    phase_nums = [_latest_es_phase_num_by_station.get(st, 0.0) for st in DISPLAY_ORDER
                  if st in _latest_es_phase_num_by_station]
    es_phase_num_repr = float(max(phase_nums)) if phase_nums else 0.0

    fxes_repr_for_ia = np.nan
    if not history.empty:
        dyn_tmp = add_dynamics(history)
        if not dyn_tmp.empty:
            ltime = dyn_tmp['obs_dt'].max()
            cand = dyn_tmp[dyn_tmp['obs_dt'] == ltime]['fxEs'].dropna()
            if not cand.empty:
                fxes_repr_for_ia = float(cand.max())
    fxes_norm_repr = clamp(fxes_repr_for_ia / 15.0, 0.0, 1.0) if not pd.isna(fxes_repr_for_ia) else np.nan

    ia_inputs = {
        "eqi_norm":              eqi_norm_repr,
        "fxEs_norm":             fxes_norm_repr,
        "mstid_approach_index":  mstid_feats.get("mstid_approach_index", np.nan),
        "sw_score":              sw_score_net,
    }
    interaction_terms = compute_interaction_terms(ia_inputs)

    ft8_feats["eqi"]          = eqi_repr
    ft8_feats["eqi_norm"]     = eqi_norm_repr
    ft8_feats["sw_score"]     = sw_score_net
    ft8_feats["es_phase_num"] = es_phase_num_repr
    ft8_feats.update(interaction_terms)

    _mstid_strength_val = mstid_feats.get('mstid_strength', np.nan)
    mstid_wave_clear = (not pd.isna(_mstid_strength_val)) and float(_mstid_strength_val) >= 0.5

    results = []; predlog_rows = []
    for st in DISPLAY_ORDER:
        g = dyn[dyn['station']==st].copy()
        if g.empty:
            continue
        g = g.sort_values('obs_dt')
        latest = g.iloc[-1]
        fx = float(latest['fxEs'])
        trend = latest.get('d_eff', np.nan)
        accel = latest.get('d2_eff', np.nan)
        fxes_std60 = latest.get('fxes_std_60', np.nan)
        cscore = cluster_score(sample, st)
        cluster_text = cluster_phrase(cscore, st)
        vstruct = compute_vertical_structure_features(g)
        prof = load_wind_profile_nearest(latest['obs_dt'], st, wind_file)
        shear_feats = compute_shear_features(prof, vstruct, g, force_no_wind=force_no_wind)
        gw_feats = compute_gw_residual_features(g)
        ef_feats = compute_f_region_background(pd.DataFrame([latest]), station=st, fxEs=fx)
        fx30_raw, forecast_model = forecast_30m(g, stage, shear_feats, vstruct, gw_feats, ef_feats, cycles)
        fx_for_prob = fx if cycles < 3 or pd.isna(fx30_raw) else fx30_raw
        p27 = empirical_probability(fx_for_prob, trend, cscore, ns_score, shear_feats, vstruct, gw_feats, ef_feats, mstid_feats, cycles, profile, ft8_feats)
        state = state_from_prob(p27, cycles, profile)

        local_scores = compute_local_shape_scores(
            cscore, shear_feats, gw_feats, fxes_std60, mstid_wave_clear,
        )
        local_shape = local_shape_label(local_scores)

        def _eff(w_key: str, raw_val, norm_fn=None) -> float:
            """重み×正規化済み特徴量値を返す。NaN時は0.0。
            Ver10.8J: norm_fnがstd_z()の場合、スケーラー未学習(旧プロファイルからの
            アップグレード直後等)だとNaNを返しうるため、norm_fn適用後の値も
            NaNチェックする（raw_valだけのチェックでは不十分）。
            """
            w = float(profile.get(w_key, 0.0))
            if pd.isna(raw_val):
                return 0.0
            v = float(norm_fn(raw_val)) if norm_fn else float(raw_val)
            if pd.isna(v):
                return 0.0
            return abs(w * v)

        eff_contribs = {
            "fxEs_trend":      _eff("w_trend",          trend,          lambda x: clamp(x, -4, 4)),
            "cluster":         _eff("w_cluster",         cscore,         lambda x: clamp(x, 0, 1)),
            "ns_wave":         _eff("w_ns",              ns_score,       lambda x: clamp(x, 0, 1)),
            "shear":           _eff("w_shear",           shear_feats.get("convergence_index"), lambda x: min(1.5, x)),
            "gw_damp":         _eff("w_gw_damp",         gw_feats.get("forecast_damping"),    lambda x: min(2.0, x)),
            "mstid_approach":  _eff("w_mstid_approach",  mstid_feats.get("mstid_approach_index")),
            "foF2_es_ratio":   _eff("w_fof2",            ef_feats.get("foF2_es_ratio"),       lambda x: std_z(x, "foF2_es_ratio", profile)),
            "dfoF2":           _eff("w_dfof2",           ef_feats.get("dfoF2"),               lambda x: std_z(x, "dfoF2", profile)),
            "hmF2":            _eff("w_hmf2",            ef_feats.get("hmF2"),                lambda x: std_z(x, "hmF2", profile)),
            "dhmF2":           _eff("w_dhmf2",           ef_feats.get("dhmF2"),               lambda x: std_z(x, "dhmF2", profile)),
            "m3000f2":         _eff("w_m3000f2",         ef_feats.get("m3000f2_norm"),        lambda x: std_z(x, "m3000f2_norm", profile)),
            "rpi":             _eff("w_rpi",              ft8_feats.get("rpi_norm"),           lambda x: std_z(x, "rpi_norm", profile)),
            "es_mode":         _eff("w_es_mode",          ft8_feats.get("es_mode_score"),      lambda x: std_z(x, "es_mode_score", profile)),
        }
        eff_total = sum(eff_contribs.values())
        if eff_total > 0:
            eff_pct = {k: round(v / eff_total * 100, 1) for k, v in eff_contribs.items()}
        else:
            eff_pct = {k: 0.0 for k in eff_contribs}

        row = {
            'station': st, 'name': STATIONS[st]['jp'], 'obs': latest['obs_dt'], 'fxEs': fx,
            'r18': fx/TARGET_18_MHZ, 'r27': fx/TARGET_27_MHZ, 'trend': trend,
            'cluster': cscore, 'cluster_phrase': cluster_text, 'north_south': ns_phrase,
            'fx30': fx30_raw, 'r27_30': fx30_raw/TARGET_27_MHZ if not pd.isna(fx30_raw) else np.nan,
            'p27': p27, 'state': state, 'forecast_model': forecast_model,
            'analysis_rate': station_rates.get(st,0.0), 'cycles': int(g['obs_dt'].nunique()),
            'vstruct': vstruct, 'shear': shear_feats, 'gw': gw_feats, 'ef': ef_feats, 'mstid': mstid_feats,
            'eff_pct': eff_pct,
            'local_shape': local_shape, 'local_shape_scores': local_scores,
            'd2_eff': accel, 'fxes_std_60': fxes_std60,
            'es_centroid_lat': motion.get('es_centroid_lat', np.nan),
            'd_centroid_lat': motion.get('d_centroid_lat', np.nan),
        }
        results.append(row)
        predlog_rows.append({
            'datetime_jst': latest['obs_dt'].strftime('%Y-%m-%d %H:%M:%S'), 'station': st,
            'fxEs': fx, 'fx30': fx30_raw, 'r18': row['r18'], 'r27': row['r27'], 'r27_30': row['r27_30'],
            'trend': trend, 'cluster': cscore, 'north_south': ns_phrase,
            'shear_peak_strength': shear_feats.get('shear_peak_strength'), 'shear_peak_height': shear_feats.get('shear_peak_height'),
            'convergence_index': shear_feats.get('convergence_index'), 'layer_height': vstruct.get('layer_height'),
            'layer_thickness': vstruct.get('layer_thickness'), 'blanketing_index': vstruct.get('blanketing_index'),
            'height_descent_rate': vstruct.get('height_descent_rate'), 'gw_residual': gw_feats.get('gw_residual'),
            'gw_modulation_index': gw_feats.get('gw_modulation_index'), 'forecast_damping': gw_feats.get('forecast_damping'),
            'f_region_background_index': ef_feats.get('f_region_background_index'), 'ef_coupling_flag': ef_feats.get('ef_coupling_flag'),
            'masking_risk_index': ef_feats.get('masking_risk_index'),
            'd2_eff': accel, 'fxes_std_60': fxes_std60,
            'es_centroid_lat': motion.get('es_centroid_lat', np.nan),
            'd_centroid_lat': motion.get('d_centroid_lat', np.nan),
            'foF2': ef_feats.get('foF2'), 'foF2_ratio': ef_feats.get('foF2_ratio'),
            'foF2_es_ratio': ef_feats.get('foF2_es_ratio'),
            'dfoF2': ef_feats.get('dfoF2'), 'hmF2': ef_feats.get('hmF2'),
            'dhmF2': ef_feats.get('dhmF2'), 'M3000F2': ef_feats.get('M3000F2'),
            'm3000f2_norm': ef_feats.get('m3000f2_norm'),
            'mstid_dir_deg': mstid_feats.get('mstid_dir_deg'),
            'mstid_speed_kmh': mstid_feats.get('mstid_speed_kmh'), 'mstid_dist_km': mstid_feats.get('mstid_dist_km'),
            'mstid_strength': mstid_feats.get('mstid_strength'), 'mstid_approach_index': mstid_feats.get('mstid_approach_index'),
            'rpi': ft8_feats.get('rpi'), 'rpi_norm': ft8_feats.get('rpi_norm'),
            'es_hop_ratio': ft8_feats.get('es_hop_ratio'),
            'es_mode_score': ft8_feats.get('es_mode_score'),
            'propagation_mode': ft8_feats.get('mode_label'),
            'psk_heard_total': ft8_feats.get('heard_total'),
            'psk_ok': ft8_feats.get('psk_ok'), 'psk_err': ft8_feats.get('psk_err'),
            'distance_std': ft8_feats.get('distance_std'),
            'bearing_std': ft8_feats.get('bearing_std'),
            'dist_peak_count': ft8_feats.get('dist_peak_count'),
            'bearing_axial_r': ft8_feats.get('bearing_axial_r'),
            'bearing_orientation_deg': ft8_feats.get('bearing_orientation_deg'),
            'bearing_ns_alignment': ft8_feats.get('bearing_ns_alignment'),
            'dist_in_window_count': ft8_feats.get('dist_in_window_count'),
            'ft8_confidence_label': ft8_feats.get('confidence_label'),
            'ft8_low_conf_area_count': ft8_feats.get('low_conf_area_count'),
            'guide_percent': p27, 'decision': state, 'analysis_phase': phase_label(cycles),
            **{f"eff_pct_{k}": v for k, v in eff_pct.items()},
        })
    if results:
        shapes_ns = [r['local_shape'] for r in results]
        continuity_penalty = spatial_continuity_penalty(shapes_ns)
        for r in results:
            r['local_shape_confidence'] = local_shape_confidence(
                r['local_shape_scores'],
                analysis_rate=r.get('analysis_rate'),
                eqi=ft8_feats.get('eqi'),
                mstid_wave_clear=mstid_wave_clear,
                continuity_penalty=continuity_penalty,
            )

    eff_keys = list(results[0]['eff_pct'].keys()) if results else []
    eff_avg: Dict[str, float] = {}
    if results and eff_keys:
        for k in eff_keys:
            vals = [r['eff_pct'][k] for r in results if not pd.isna(r['eff_pct'].get(k, np.nan))]
            eff_avg[k] = round(float(np.mean(vals)), 1) if vals else 0.0
    meta = {'analysis_rate': analysis_rate, 'cycles': cycles, 'phase': phase_label(cycles),
            'motion': motion, 'north_south': ns_phrase, 'predlog_rows': predlog_rows,
            'mstid_feats': mstid_feats, 'eff_avg': eff_avg,
            'ft8_feats': ft8_feats,
            'quality_hint': profile.get('_quality_hint', '')}
    return results, meta


def build_mock_latest_sample(step: int=0) -> pd.DataFrame:
    base_dt = now_jst().replace(second=0, microsecond=0) - timedelta(minutes=15*max(0, 3-step))
    vals = {'wakkanai': 4.8+0.2*step, 'kokubunji': 10.2+0.6*step, 'yamagawa': 11.8+0.75*step, 'okinawa': 5.2+0.25*step}
    return pd.DataFrame([{'obs_jst': base_dt.isoformat(), 'station': st, 'fxEs': vals[st]} for st in DISPLAY_ORDER])


def prepare_mock_session_history(n_steps: int=6) -> pd.DataFrame:
    hist = empty_history()
    for i in range(n_steps):
        hist, _ = merge_history_in_memory(hist, build_mock_latest_sample(i))
    return hist


def create_mock_mstid_images(image_dir: Path):
    image_dir.mkdir(parents=True, exist_ok=True)
    for idx, shift in enumerate([0, 8]):
        arr = np.zeros((240,320), dtype=np.uint8) + 16
        for k in range(30):
            x = 40 + k*7 + shift; y = 180 - k*4
            if 0 <= x < 320 and 0 <= y < 240:
                arr[max(0,y-1):min(240,y+2), max(0,x-12):min(320,x+13)] = 200
        Image.fromarray(arr).save(image_dir / f'mock_{idx}.png')



def ui_set_manual_location(lat: Optional[float], lon: Optional[float], label: str = "手動") -> None:
    """Ver11.7: GPS/APIが使えない場合の手動位置設定。"""
    global _ui_manual_location, _ui_last_location, _ui_gps_status, _ui_gps_reason, _ui_location_source
    if lat is None or lon is None or pd.isna(lat) or pd.isna(lon):
        return
    _ui_manual_location = (float(lat), float(lon), label)
    _ui_last_location = _ui_manual_location
    _ui_gps_status = "手動座標"
    _ui_gps_reason = "--my-lat/--my-lon または手入力"
    _ui_location_source = label


def ui_geocode_place_online(place: str) -> Optional[Tuple[float, float, str]]:
    """Ver11.7: オンライン地名API(Nominatim)で地名/施設名から座標取得。"""
    q = str(place or '').strip()
    if not q:
        return None
    try:
        url = "https://nominatim.openstreetmap.org/search"
        params = {"q": q, "format": "json", "limit": 1, "countrycodes": "jp"}
        headers = {"User-Agent": "EDFS-Ver11.7-Pydroid3-Geocode/1.0"}
        r = requests.get(url, params=params, headers=headers, timeout=10)
        r.raise_for_status()
        data = r.json()
        if not data:
            return None
        lat = float(data[0]["lat"]); lon = float(data[0]["lon"])
        label = data[0].get("display_name", q).split(',')[0]
        return (lat, lon, label)
    except Exception:
        return None


def ui_try_gps_once(wait_sec: float = 5.0) -> Optional[Tuple[float, float, str]]:
    """Ver11.7: GPSを一度だけ試行。Pydroid3 TerminalではJNIEnv問題を理由表示する。"""
    global _ui_last_location, _ui_gps_status, _ui_gps_reason, _ui_location_source
    if not _ui_gps_enabled:
        _ui_gps_status = "無効"; _ui_gps_reason = "--no-gps 指定"
        return None
    if 'ANDROID_ARGUMENT' not in os.environ:
        _ui_gps_status = "未取得"
        _ui_gps_reason = "Pydroid3 Terminal実行のためGPS不可(JNIEnv制限)"
        return None
    try:
        from plyer import gps
        holder: Dict[str, float] = {}
        def _on_location(**kwargs):
            lat = kwargs.get('lat') or kwargs.get('latitude')
            lon = kwargs.get('lon') or kwargs.get('longitude')
            if lat is not None and lon is not None:
                holder['lat'] = float(lat); holder['lon'] = float(lon)
        gps.configure(on_location=_on_location)
        gps.start(minTime=1000, minDistance=0)
        loops = max(1, int(wait_sec / 0.1))
        for _ in range(loops):
            if 'lat' in holder and 'lon' in holder:
                break
            time.sleep(0.10)
        try:
            gps.stop()
        except Exception:
            pass
        if 'lat' in holder and 'lon' in holder:
            _ui_last_location = (holder['lat'], holder['lon'], 'GPS')
            _ui_gps_status = "取得成功"; _ui_gps_reason = "GPS Fix OK"; _ui_location_source = "GPS"
            return _ui_last_location
        _ui_gps_status = "未取得"; _ui_gps_reason = "Fix待ち時間超過/権限未許可の可能性"
        return None
    except ModuleNotFoundError as e:
        _ui_gps_status = "未取得"; _ui_gps_reason = f"plyer未導入/非対応: {str(e)[:28]}"
        return None
    except Exception as e:
        msg = str(e).replace('\n', ' ')
        _ui_gps_status = "未取得"
        _ui_gps_reason = "JNIEnv unavailable" if "JNIEnv" in msg else f"GPS例外: {msg[:32]}"
        return None


def ui_try_termux_gps_file() -> Optional[Tuple[float, float, str]]:
    """Ver11.7統合版: Termux GPSブリッジJSONファイルを読み取る。

    対応JSON:
      {"lat":..,"lon":..,"timestamp":"UTC ISO8601","source":"termux-location"}
      {"latitude":..,"longitude":..,"time":...}

    既定では10分以内のファイルのみ有効。timestampが無い場合はファイル更新時刻を使う。

    Ver11.8: --no-gps 指定時（_ui_gps_enabled=False）はTermuxブリッジも
    含めてGPS系取得を一括無効化する（従来はui_try_gps_once()側でしか
    判定しておらず、--no-gps指定後もTermuxブリッジ経由の位置取得が
    動作してしまう不整合があった）。
    """
    global _ui_last_location, _ui_gps_status, _ui_gps_reason, _ui_location_source
    if not _ui_gps_enabled:
        _ui_gps_status = "無効"; _ui_gps_reason = "--no-gps 指定"
        return None
    path = _ui_termux_gps_file
    if not os.path.exists(path):
        _ui_gps_status = "未取得"
        _ui_gps_reason = f"Termux GPSファイルなし({path})"
        return None
    try:
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        lat_raw = data.get('lat', data.get('latitude'))
        lon_raw = data.get('lon', data.get('longitude'))
        lat = float(lat_raw)
        lon = float(lon_raw)

        ts_raw = str(data.get('timestamp', '') or data.get('time', '') or '')
        if ts_raw:
            try:
                ts = datetime.strptime(ts_raw.replace('+00:00', 'Z'), "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
                age_sec = (datetime.now(timezone.utc) - ts).total_seconds()
            except Exception:
                age_sec = time.time() - os.path.getmtime(path)
        else:
            age_sec = time.time() - os.path.getmtime(path)

        if age_sec > _ui_termux_gps_max_age_sec:
            _ui_gps_status = "未取得"
            _ui_gps_reason = f"Termux GPSファイル期限切れ({int(age_sec/60)}分前)"
            return None
        _ui_last_location = (lat, lon, 'Termux-GPS')
        _ui_gps_status = "取得成功"
        _ui_gps_reason = f"Termux GPS({int(age_sec)}秒前)"
        _ui_location_source = "Termux-GPS"
        return _ui_last_location
    except Exception as e:
        _ui_gps_status = "未取得"
        _ui_gps_reason = f"Termux GPSファイル読取失敗: {str(e)[:28]}"
        return None

def ui_resolve_location() -> Optional[Tuple[float, float, str]]:
    """Ver11.7: 位置取得優先順位 GPS → Termux GPSブリッジ → 地名API → 緯度経度手入力。"""
    global _ui_last_location, _ui_location_dialog_done, _ui_location_source, _ui_gps_status, _ui_gps_reason
    if _ui_manual_location is not None:
        _ui_last_location = _ui_manual_location
        return _ui_last_location
    if _ui_last_location is not None:
        return _ui_last_location
    loc = ui_try_gps_once(wait_sec=5.0)
    if loc is not None:
        return loc
    loc = ui_try_termux_gps_file()
    if loc is not None:
        return loc
    if _ui_location_dialog_done:
        return _ui_last_location
    _ui_location_dialog_done = True
    try:
        print("\nGPS未取得: " + (_ui_gps_reason or _ui_gps_status))
        if _active_html_bridge is not None:
            try:
                _active_html_bridge.write_snapshot(
                    sys.stdout.snapshot() + "\n(端末側で位置情報の手入力待ちです。Pydroid3の画面を確認してください)\n"
                )
            except Exception:
                pass
        place = input("地名/施設名を入力(空欄で緯度経度入力): ").strip()
        if place:
            loc = ui_geocode_place_online(place)
            if loc is not None:
                _ui_last_location = (loc[0], loc[1], loc[2])
                _ui_location_source = "地名API"
                _ui_gps_status = "地名API"
                _ui_gps_reason = place
                return _ui_last_location
            print("地名APIで座標を取得できませんでした。")
        lat_s = input("緯度を入力(空欄でGPS未取得のまま継続): ").strip()
        if not lat_s:
            return _ui_last_location
        lon_s = input("経度を入力(空欄でGPS未取得のまま継続): ").strip()
        if not lon_s:
            return _ui_last_location
        lat = float(lat_s); lon = float(lon_s)
        label = input("地点名(任意): ").strip() or "手入力"
        ui_set_manual_location(lat, lon, label=label)
        _ui_location_source = "手入力"
        return _ui_last_location
    except Exception as e:
        _ui_gps_reason = f"位置入力失敗: {str(e)[:24]}"
        return _ui_last_location


def ui_bearing_name_ja(deg) -> str:
    if deg is None or pd.isna(deg): return "方位未定"
    return ["北", "北東", "東", "南東", "南", "南西", "西", "北西"][int((float(deg)+22.5)//45)%8]


def ui_direction_arrow(deg) -> str:
    if deg is None or pd.isna(deg): return "→"
    return ["↑", "↗", "→", "↘", "↓", "↙", "←", "↖"][int((float(deg)+22.5)//45)%8]


def ui_stars(score) -> str:
    if score is None or pd.isna(score): return "★"
    s = float(score); s = s/100.0 if s > 1.0 else s
    return "★" * max(1, min(5, int(round(s*5))))


def ui_future_status_from_fxes(fx_mhz, cycles: int) -> str:
    if cycles < 3 or fx_mhz is None or pd.isna(fx_mhz): return "[--不明]"
    fx = float(fx_mhz)
    if fx >= 12.0: return "[★強開★]"
    if fx >= 9.0:  return "[●開通 ]"
    if fx >= 7.0:  return "[▲注意 ]"
    if fx >= 5.0:  return "[△弱兆 ]"
    return "[×低調 ]"


def ui_open_direction_items(area_df: pd.DataFrame, loc: Optional[Tuple[float, float, str]]) -> List[Dict]:
    if area_df is None or area_df.empty: return []
    rows: List[Dict] = []
    for _, row in area_df.iterrows():
        area = row.get('area'); name = row.get('area_name', format_area_label(area))
        ratio = row.get('reach_ratio', np.nan); heard = int(row.get('heard_count', 0) or 0)
        score = (0.0 if heard == 0 else min(0.25, heard/10.0)) if pd.isna(ratio) else float(ratio)
        if score <= 0.0 and heard <= 0: continue
        bear = np.nan
        if loc is not None and area in AREA_CENTERS:
            bear = bearing_deg(loc[0], loc[1], AREA_CENTERS[area][0], AREA_CENTERS[area][1])
        rows.append({'area': area, 'name': name, 'score': score, 'heard': heard,
                     'bearing': bear, 'bearing_name': ui_bearing_name_ja(bear), 'arrow': ui_direction_arrow(bear)})
    rows.sort(key=lambda x: x['score'], reverse=True)
    return rows


def ui_footer() -> str:
    return "N通常 D詳細 F寄与 E学習 GDebug Q終了"


def ui_current_observation_lines(results: List[Dict], meta: Dict) -> List[str]:
    lines: List[str] = ["現在観測"]
    order = ["稚内", "国分寺", "山川", "沖縄"]
    used = set()
    for name in order:
        r = next((x for x in results if name in str(x.get('name', ''))), None)
        if r is None: continue
        used.add(id(r)); fx = r.get('fxEs', np.nan)
        mark = '★' if (not pd.isna(fx) and float(fx)>=12.0) else ('▲' if (not pd.isna(fx) and float(fx)>=8.0) else ('△' if (not pd.isna(fx) and float(fx)>=5.0) else '×'))
        lines.append(f"{name:<4} {fmt_mhz(fx,1):>7} {mark}")
    for r in results:
        if id(r) in used: continue
        lines.append(f"{str(r.get('name','--'))[:4]:<4} {fmt_mhz(r.get('fxEs'),1):>7}")
    return lines


def ui_learning_progress_lines(meta: Dict) -> List[str]:
    lines: List[str] = ["学習", f"EDFS Lv{_edfs_level} {_edfs_level_label}", f"Feature {_edfs_feature_set}"]
    area_need = EDFS_LV2_AREA_TEACHERS
    area_pct = clamp((_edfs_area_teachers/area_need)*100.0, 0.0, 100.0) if area_need > 0 else np.nan
    area_left = max(0, area_need-int(_edfs_area_teachers))
    lines.append(f"Area教師 {_edfs_area_teachers}/{area_need} ({fmt(area_pct,0)}%)")
    lines.append("Area良好 到達済み" if area_left <= 0 else f"Area良好まで あと{area_left}教師")
    cent_need = EDFS_LV3_CENTROID_HIST
    cent_pct = clamp((_edfs_centroid_history_count/cent_need)*100.0, 0.0, 100.0) if cent_need > 0 else np.nan
    cent_left = max(0, cent_need-int(_edfs_centroid_history_count))
    lines.append(f"重心履歴 {_edfs_centroid_history_count}/{cent_need} ({fmt(cent_pct,0)}%)")
    lines.append("Motion良好 到達済み" if cent_left <= 0 else f"Motion良好まで あと{cent_left}件")
    q = str(meta.get('quality_hint', '') or '').strip()
    if q: lines.append(q)
    return lines


def poll_stdin_nonblocking() -> Optional[str]:
    """Ver11.8修正: WindowsのselectはUNIXソケット専用でありsys.stdin
    （ファイルオブジェクト）を渡すとOSErrorになるため、待機中の
    N/D/F/E/G/Q モード切替が常に無反応になっていた。is_windows()時は
    msvcrt.kbhit()/getwch()による1文字非ブロッキング読み取りに切り替える。
    """
    if is_windows() and _MSVCRT_AVAILABLE:
        try:
            if msvcrt.kbhit():
                ch = msvcrt.getwch()
                return ch.strip() if ch else None
        except Exception:
            return None
        return None
    try:
        ready, _, _ = select.select([sys.stdin], [], [], 0)
        if ready: return sys.stdin.readline().strip()
    except Exception:
        return None
    return None


def handle_runtime_command(cmd: str, use_color: bool=True) -> bool:
    global _ui_display_mode
    if not cmd: return True
    c = cmd.strip().upper()[:1]
    maps = {'N':UI_MODE_NORMAL, 'D':UI_MODE_DETAIL, 'F':UI_MODE_FACTOR, 'E':UI_MODE_EDFS, 'G':UI_MODE_DEBUG}
    if c == 'Q':
        print("\n終了要求を受け付けました。"); return False
    if c in maps:
        _ui_display_mode = maps[c]
        print(f"\n表示モード: {_ui_display_mode}")
        redraw_last_runtime_display(use_color=use_color)
    return True


def redraw_last_runtime_display(use_color: bool=True) -> None:
    if _ui_last_results is None or _ui_last_meta is None: return
    print_unified_display(_ui_last_results, _ui_last_meta,
                          detail=(_ui_display_mode == UI_MODE_DETAIL),
                          use_color=use_color, android=_ui_last_android)


def print_android_normal_v117(results: List[Dict], meta: Dict, use_color: bool=True) -> None:
    W = DISP_W_ANDROID
    cycles = int(meta.get('cycles', 0)); ft8 = meta.get('ft8_feats', {}) or {}; motion = meta.get('motion', {}) or {}
    loc = ui_resolve_location()
    best = max(results, key=lambda x: -1 if pd.isna(x.get('p27', np.nan)) else x.get('p27', -1))
    now_p = best.get('p27', np.nan); now_state = color_or_symbol_state(best.get('state','判定不可'), use_color=use_color, android=True)
    future_fxes = best.get('fx30', np.nan); future_state = ui_future_status_from_fxes(future_fxes, cycles)
    delta = np.nan if (pd.isna(future_fxes) or pd.isna(best.get('fxEs',np.nan))) else float(future_fxes)-float(best.get('fxEs'))
    if pd.isna(delta): trend_mk, trend_txt = "", "未収束"
    elif delta >= 0.5: trend_mk, trend_txt = "↑", "上昇中"
    elif delta <= -0.5: trend_mk, trend_txt = "↓", "低下中"
    else: trend_mk, trend_txt = "→", "横ばい"
    print("現在"); print(f"{now_state} {fmt(now_p,0)}%"); print("")
    print("30分後"); print(f"{future_state} {trend_mk}")
    print(f"予測値 {fmt_mhz(future_fxes,1)}  {trend_txt}" if cycles >= 3 and not pd.isna(future_fxes) else "予測値 --MHz  未収束")
    print("-"*W)
    for line in ui_current_observation_lines(results, meta): print(line[:W])
    print("-"*W)
    if loc is None:
        print("位置 GPS未取得")
        print(f"理由 {_ui_gps_reason or _ui_gps_status}"[:W])
    else:
        print(f"位置 {_ui_location_source}: {loc[0]:.5f}, {loc[1]:.5f}"[:W])
    items = ui_open_direction_items(_ui_last_area_df, loc)
    if items:
        primary_threshold = max(items[0]['score']*0.75, 0.50)
        primary = [x for x in items if x['score'] >= primary_threshold] or [items[0]]
        secondary = [x for x in items if x not in primary and x['score'] >= 0.25]
        print(""); print("主開通")
        top = primary[0]; print(f"{top['arrow']} {top['bearing_name']}")
        for it in primary[:3]: print(f"{ui_stars(it['score'])} {it['name']}")
        if secondary:
            print(""); print("副開通")
            for it in secondary[:2]: print(f"{it['arrow']} {it['name']} {it['score']*100:.0f}%")
    else:
        print(""); print("主開通 判定不可")
    print("-"*W)
    print(f"状態 {ft8.get('mode_label','---')}"); print(f"信頼 {ft8.get('confidence_label','---')}")
    print(f"整合 {'一致' if ft8.get('consistency_match',False) else '不一致'}" if ft8.get('psk_ok', True) else "整合 FT8通信障害")
    print("-"*W); print("移動"); print(f"{motion.get('direction','不明')}  {fmt(motion.get('speed_kmh'),0)}km/h")
    print("-"*W)
    for line in ui_learning_progress_lines(meta): print(line[:W])
    print("="*W); print(ui_footer())


def print_android_detail_v117(results: List[Dict], meta: Dict, use_color: bool=True) -> None:
    W = DISP_W_ANDROID; print("DETAIL 局別詳細"); print("-"*W)
    cycles = int(meta.get('cycles',0))
    for r in results:
        state = color_or_symbol_state(r.get('state','判定不可'), use_color=use_color, android=True)
        print(f"{r.get('name','--'):<4} {fmt_mhz(r.get('fxEs'),1):>7} {state}")
        print(f" 30分後 {fmt_mhz(r.get('fx30'),1) if cycles >= 3 else '--MHz'}  目安 {fmt(r.get('p27'),0)}%")
        print(f" 傾向 {trend_label(r.get('trend'))}  形状 {r.get('local_shape','不明')}")
        print("-"*W)
    print(ui_footer())


def print_android_factor_v117(results: List[Dict], meta: Dict, use_color: bool=True) -> None:
    W = DISP_W_ANDROID; print("FACTOR 寄与率"); print("-"*W)
    eff = meta.get('eff_avg', {}) or {}
    if not eff: print("寄与率データなし")
    for k,v in sorted(eff.items(), key=lambda kv: kv[1], reverse=True)[:8]: print(f"{k:<12} {bar_graph(float(v), width=12)} {float(v):>4.0f}%")
    print(ui_footer())


def print_android_edfs_v117() -> None:
    W = DISP_W_ANDROID; print("EDFS 学習状態"); print("-"*W)
    print(f"Level    {_edfs_level} {_edfs_level_label}"); print(f"Feature  {_edfs_feature_set}")
    print(f"Teachers {_edfs_area_teachers}"); print(f"Centroid {_edfs_centroid_history_count}"); print(f"Lifetime {_edfs_lifetime_samples}")
    print(f"Area良好まで {max(0, EDFS_LV2_AREA_TEACHERS-int(_edfs_area_teachers))}教師")
    print(f"Motion良好まで {max(0, EDFS_LV3_CENTROID_HIST-int(_edfs_centroid_history_count))}件")
    print("Predictor")
    for k,v in _edfs_predictor_active.items(): print(f" {k:<8} {'ACTIVE' if v else 'WAIT'}")
    print(ui_footer())


def print_android_debug_v117(results: List[Dict], meta: Dict) -> None:
    W = DISP_W_ANDROID; print("DEBUG"); print("-"*W)
    print(f"mode      {_ui_display_mode}"); print(f"cycles    {meta.get('cycles')}"); print(f"phase     {meta.get('phase')}")
    print(f"GPS状態   {_ui_gps_status}"); print(f"GPS理由   {_ui_gps_reason}"[:W])
    print(f"位置元    {_ui_location_source}"); print(ui_footer())

def print_header(stage: int, android: bool = False):
    W = DISP_W_ANDROID if android else DISP_W_NORMAL
    sep = "=" * W
    if android:
        print(sep)
        print(f"[{analysis_state_label(stage)}] {now_jst().strftime('%m/%d %H:%M JST')}")
        print(sep)
    else:
        print('\n' + sep)
        print(f" Es Monitor by JL7KHN  {display_version()}  解析状態: {analysis_state_label(stage)}  {now_jst().strftime('%Y-%m-%d %H:%M JST')}")
        print(sep)


def print_debug_info(debug: bool, new_rows: int, parsed_rows: int, latest_rows: int, url: str, logger: SessionLogger, wind_file: Path, force_no_wind: bool, mstid_enabled: bool):
    if not debug:
        return
    print('-'*78)
    print(f"デバッグ: 新規={new_rows}行  抽出={parsed_rows}行  採用={latest_rows}行  URL={url}")
    print(f"セッションログ: {'無効' if not logger.enabled else logger.file_path}")
    print(f"風プロファイル: {'無効化' if force_no_wind else wind_file}")
    print(f"MSTID画像解析: {'有効' if mstid_enabled else '無効'}")


def forecast_text(r: Dict, cycles: int) -> str:
    if cycles < 3 or pd.isna(r.get('fx30')):
        return '30分後=--MHz'
    return f"30分後={fmt_mhz(r['fx30'],1)}{'(参考)' if cycles < MIN_ANALYSIS_POINTS else ''}"


def guide_text(p, cycles: int) -> str:
    if cycles < 3:
        return '目安=参考'
    return f"目安={fmt(p,0)}%{'(参考)' if cycles < MIN_ANALYSIS_POINTS else ''}"


def directional_field_display_text(directional_field: Optional[Dict], compact: bool = False) -> str:
    """Rev14.6(#1): print_unified_display()向けに、directional_field
    (Rev14.5/combine_directional_evidence()の出力)を1行のテキストへ整形する。
    FT8観測の有無で表現を使い分ける(build_gemma_prompt()と同じ命名方針:
    FT8が無い場合は「伝播方位」ではなく「方位」+出典明示にとどめる)。
    データが無い場合は空文字を返す(呼び出し側で行ごと省略できるように)。
    """
    if not directional_field:
        return ""
    field = directional_field.get("propagation_directional_confidence") or \
        directional_field.get("directional_opening_expectancy")
    if not field:
        return ""
    best_dir = max(field, key=lambda k: field[k])
    src = directional_field.get("primary_evidence_source")
    if src == "ft8_observed":
        src_tag = "FT8実測" if compact else "FT8実測を含む推定"
    elif src == "nict_model_only":
        src_tag = "NICT推定" if compact else "NICT局データからの粗い推定(実測ではない)"
    else:
        return ""
    trend_tag = ""
    nict_grad = directional_field.get("nict_gradient") or {}
    trend = (nict_grad.get("trend") or {}).get("trend")
    if trend == "SOUTHWARD":
        trend_tag = " 南下傾向" if not compact else "↓"
    elif trend == "NORTHWARD":
        trend_tag = " 北上傾向" if not compact else "↑"
    if compact:
        return f"方位={best_dir}方向({src_tag}){trend_tag}"
    return f"方位別開通期待: {best_dir}方向が優勢({src_tag}){trend_tag}"


def print_unified_display(results: List[Dict], meta: Dict, detail: bool, use_color: bool=True, android: bool=False):
    global _ui_last_results, _ui_last_meta, _ui_last_android, _ui_last_f2_forecast
    _ui_last_results = results; _ui_last_meta = meta; _ui_last_android = android
    try:
        _ui_last_f2_forecast = compute_f2_layer_forecast()
    except Exception:
        _ui_last_f2_forecast = []

    W = DISP_W_ANDROID if android else DISP_W_NORMAL
    sep_thin = "-" * W

    if not results:
        print('判定可能なデータがありません。')
        return

    if android:
        if _ui_display_mode == UI_MODE_NORMAL:
            print_android_normal_v117(results, meta, use_color=use_color); return
        if _ui_display_mode == UI_MODE_DETAIL:
            print_android_detail_v117(results, meta, use_color=use_color); return
        if _ui_display_mode == UI_MODE_FACTOR:
            print_android_factor_v117(results, meta, use_color=use_color); return
        if _ui_display_mode == UI_MODE_EDFS:
            print_android_edfs_v117(); return
        if _ui_display_mode == UI_MODE_DEBUG:
            print_android_debug_v117(results, meta); return
    motion = meta.get('motion', {})
    conf_pct = np.nan if pd.isna(motion.get('confidence', np.nan)) else motion.get('confidence', 0.0)*100.0
    cycles = int(meta.get('cycles', 0))
    mst = meta.get('mstid_feats', {})
    quality_hint = meta.get('quality_hint', '')

    if android:
        print(f"取得={cycles}回  分析={fmt(meta.get('analysis_rate'),0)}%")
        print(f"段階={meta.get('phase','不明')}  信頼={fmt(conf_pct,0)}%")
        if quality_hint:
            print(quality_hint)
        print(f"Es移動={motion.get('direction','不明')} {fmt(motion.get('speed_kmh'),0)}km/h")
        _dir_text = directional_field_display_text(meta.get('directional_field'), compact=True)
        if _dir_text:
            print(_dir_text)
        if not pd.isna(mst.get('mstid_approach_index', np.nan)):
            print(f"MSTID整合度(仮定Es軸45°)={fmt(mst.get('mstid_approach_index'),3)}")
        print(sep_thin)
        for i, r in enumerate(results):
            sym = color_or_symbol_state(r['state'], use_color=use_color, android=android)
            print(f"【{r['name']}】{sym}")
            print(f" fxEs={fmt_mhz(r['fxEs'],1)}  27MHz={ratio_label(r['r27'])}")
            print(f" 目安={fmt(r['p27'],0)}%  {forecast_text(r, cycles)}")
            print(f" 傾向={trend_label(r['trend'])}  広域={r['cluster_phrase']}")
            print(f" Es形状={r['local_shape']}  Confidence={fmt(r['local_shape_confidence'],0)}%")
            if detail:
                vstruct = r['vstruct']; shear = r['shear']; gw = r['gw']; ef = r['ef']; mst2 = r['mstid']
                print(f" 高度={fmt(vstruct.get('layer_height'),1)}km  層厚={layer_thickness_phrase(vstruct.get('layer_thickness'))}")
                print(f" shear={fmt(shear.get('convergence_index'),2)}  GW={fmt(gw.get('gw_modulation_index'),2)}")
                fof2_v = ef.get('foF2', np.nan); hmf2_v = ef.get('hmF2', np.nan)
                dfof2_v = ef.get('dfoF2', np.nan); dhmf2_v = ef.get('dhmF2', np.nan)
                if not (pd.isna(fof2_v) and pd.isna(hmf2_v)):
                    ef_flag = "↓F下降" if ef.get('ef_coupling_flag', 0.0) == 1.0 else ""
                    print(f" foF2={fmt(fof2_v,1)}MHz  hmF2={fmt(hmf2_v,0)}km {ef_flag}")
            if i != len(results) - 1:
                print(sep_thin)
        best = max(results, key=lambda x: -1 if pd.isna(x['p27']) else x['p27'])
        print(sep_thin)
        print(f"★最有力: {best['name']}")
        print(f"  目安={fmt(best['p27'],0)}%  {color_or_symbol_state(best['state'], use_color, android)}")
    else:
        print(f"分析={fmt(meta.get('analysis_rate'),0)}%  取得={cycles}回  段階={meta.get('phase','不明')}  Es移動={motion.get('direction','不明')}  速度={fmt(motion.get('speed_kmh'),0)}km/h  信頼={fmt(conf_pct,0)}%")
        _dir_text = directional_field_display_text(meta.get('directional_field'), compact=False)
        if _dir_text:
            print(_dir_text)
        if quality_hint:
            print(quality_hint)
        if not pd.isna(mst.get('mstid_approach_index', np.nan)):
            print(f"MSTID: 方向={fmt(mst.get('mstid_dir_deg'),0)}°  速度={fmt(mst.get('mstid_speed_kmh'),0)}km/h  距離={fmt(mst.get('mstid_dist_km'),0)}km  仮定Es軸整合度={fmt(mst.get('mstid_approach_index'),3)}")
        print(sep_thin)
        for i, r in enumerate(results):
            state = color_or_symbol_state(r['state'], use_color=use_color, android=android)
            print(f"{r['name']:<4} fxEs={fmt_mhz(r['fxEs'],1)}  27MHz={ratio_label(r['r27'])}  {forecast_text(r, cycles)}  {guide_text(r['p27'], cycles)}  {state}")
            print(f"      傾向={trend_label(r['trend'])}  広域={r['cluster_phrase']}  Es形状={r['local_shape']}  Confidence={fmt(r['local_shape_confidence'],0)}%")
            print(f"      分析={fmt(r['analysis_rate'],0)}%  18MHz目安={ratio_label(r['r18'])}  27MHz予測={ratio_label(r['r27_30']) if cycles >= 3 else '不明'}")
            if detail:
                vstruct = r['vstruct']; shear = r['shear']; gw = r['gw']; ef = r['ef']; mst2 = r['mstid']
                print(f"      構造: Es高度={fmt(vstruct.get('layer_height'),1)}km  層厚={layer_thickness_phrase(vstruct.get('layer_thickness'))}  ブランケット={blanketing_phrase(vstruct.get('blanketing_index'))}")
                print(f"      [根拠] R18={fmt(r['r18'],2)}  R27={fmt(r['r27'],2)}  R27+30={fmt(r['r27_30'],2)}  d={fmt(r['trend'],1)}MHz/h")
                print(f"             shear={fmt(shear.get('convergence_index'),2)}({shear.get('shear_mode')})  GW={fmt(gw.get('gw_modulation_index'),2)}  F層影響={f_region_phrase(ef.get('masking_risk_index'))}  model={r['forecast_model']}")
                fof2_v = ef.get('foF2', np.nan); hmf2_v = ef.get('hmF2', np.nan)
                dfof2_v = ef.get('dfoF2', np.nan); dhmf2_v = ef.get('dhmF2', np.nan)
                if not (pd.isna(fof2_v) and pd.isna(hmf2_v)):
                    ef_flag = "↓F高度下降" if ef.get('ef_coupling_flag', 0.0) == 1.0 else ""
                    print(f"             F層: foF2={fmt(fof2_v,1)}MHz(Δ{fmt(dfof2_v,1)})  hmF2={fmt(hmf2_v,0)}km(Δ{fmt(dhmf2_v,0)}) {ef_flag}")
                if not pd.isna(mst2.get('mstid_approach_index', np.nan)):
                    print(f"             MSTID dir={fmt(mst2.get('mstid_dir_deg'),0)}° speed={fmt(mst2.get('mstid_speed_kmh'),0)}km/h dist={fmt(mst2.get('mstid_dist_km'),0)}km strength={fmt(mst2.get('mstid_strength'),3)} axis_align(assumed45)={fmt(mst2.get('mstid_approach_index'),3)}")
            if i != len(results)-1:
                print()
        best = max(results, key=lambda x: -1 if pd.isna(x['p27']) else x['p27'])
        print(sep_thin)
        print(f"最有力: {best['name']}周辺  目安={fmt(best['p27'],0)}%  判定={color_or_symbol_state(best['state'], use_color=use_color, android=android)}")
        print('注: 目安はNICT fxEsからの経験則。Es移動/速度と南北分布は4観測点からの簡易推定。')
        print('注: 方位別開通期待度はNICT局データ由来の粗い推定であり、FT8実測が無い場合は伝播方位の直接観測ではない。')
        print('注: 風シア・鉛直構造・GW・E-F coupling・MSTID は観測/モデル代理であり、完全物理解ではない。')



def estimate_network_masking_and_trend(history: pd.DataFrame) -> Tuple[float, float, float, float]:
    """judge_propagation_mode()/compute_consistency()用の軽量ネットワーク代表値。
    4局のうち直近fxEsが最大の局を代表局とし、analyze()内の局別ループと同じ
    compute_f_region_background()/add_dynamics()を再利用して算出する
    （定義の不整合を避けるため、独自の簡易式は作らない）。
    戻り値: (masking_risk_index, trend(d_eff), accel(d2_eff), fxEs代表値)
    """
    if history is None or history.empty:
        return np.nan, np.nan, np.nan, np.nan
    dyn = add_dynamics(history)
    if dyn.empty:
        return np.nan, np.nan, np.nan, np.nan
    latest_time = dyn['obs_dt'].max()
    cand = dyn[dyn['obs_dt'] == latest_time].dropna(subset=['fxEs'])
    if cand.empty:
        return np.nan, np.nan, np.nan, np.nan
    rep = cand.loc[cand['fxEs'].idxmax()]
    st = rep['station']; fx = float(rep['fxEs'])
    ef = compute_f_region_background(pd.DataFrame([rep]), station=st, fxEs=fx)
    return ef.get('masking_risk_index', np.nan), rep.get('d_eff', np.nan), rep.get('d2_eff', np.nan), fx


def print_ft8_header_extra(ft8_feats: Dict, android: bool = False) -> None:
    """NICT更新/FT8補正/信頼度/モデルバージョン（仕様書8章ヘッダー部 + Ver11.2J §16）"""
    W = DISP_W_ANDROID if android else DISP_W_NORMAL
    nict_t = ft8_feats.get('nict_update_str', '--')
    psk_t  = ft8_feats.get('psk_update_str', '--')
    conf   = ft8_feats.get('confidence_label', 'LOW')
    psk_ok = ft8_feats.get('psk_ok', True)
    psk_err = ft8_feats.get('psk_err', '')
    mv     = _model_version
    ms     = f"{_model_score:.0f}%"
    if android:
        model_line = f"Model {mv}  Score {ms}"
        if not psk_ok:
            print(f"NICT={nict_t} FT8=通信障害({psk_err}) 信頼={conf}")
        else:
            print(f"NICT={nict_t} FT8={psk_t} 信頼={conf}")
        print(model_line)
    else:
        print(f"NICT更新 : {nict_t}")
        if not psk_ok:
            print(f"FT8補正  : 通信障害 ({psk_err})")
        else:
            print(f"FT8補正  : {psk_t}")
        print(f"信頼度   : {conf}")
        print(f"Model    : {mv}  Score {ms}")
    print("=" * W)


def print_propagation_mode_panel(ft8_feats: Dict, android: bool = False) -> None:
    """現在状態（仕様書6・8章: Es優勢/F層優勢/混合）"""
    mode_label = ft8_feats.get('mode_label', '不明')
    mode_sub = ft8_feats.get('mode_sub', '')
    tag = {"Es優勢": "[E]", "F層優勢": "[F]", "混合": "[M]"}.get(mode_label, "[-]")
    if android:
        print(f"現在状態 {tag}{mode_label} {mode_sub}".rstrip())
        return
    W = DISP_W_NORMAL
    print("現在状態")
    print(f"{tag} {mode_label}")
    if mode_sub:
        print(mode_sub)
    print("-" * W)


def print_area_reachability_panel(area_df: pd.DataFrame,
                                   area_fvecs: Optional[Dict[str, List[float]]] = None,
                                   area_geo_eta: Optional[Dict[str, Dict]] = None,
                                   area_climatology: Optional[Dict[str, Dict]] = None,
                                   android: bool = False) -> None:
    """エリア別現在到達率＋トレンド記号＋2サイクル後予測（Ver11.2J §8/16）。
    学習中は「総:XX/100 Es:XX/30」を表示。
    学習済みは「2サイクル後予測:XX%  Bias:±X%  信頼度:XX」を表示。

    Rev14.0:
      - reach_ratio(生値)がNaNでも、reach_ratio_compensated(IDW補間/
        ブレンド値)があれば「推定」であることを明示した上でバー表示する
        (提案②(a): 観測点が無い場所の空間推定)。
      - 学習中(LEARNING)で数値予測が出せない場合でも、area_geo_etaに
        値があれば「今の動きなら約N分」という学習器に依存しない目安を
        表示する(提案②(c): 幾何学的ETA)。両者は明確に区別して表示する。
    """
    if area_df is None or area_df.empty:
        return
    area_fvecs = area_fvecs or {}
    area_geo_eta = area_geo_eta or {}
    area_climatology = area_climatology or {}

    def _eta_text(area: str) -> str:
        eta = (area_geo_eta.get(area, {}) or {}).get('eta_min')
        if eta is None:
            return ""
        if eta < 60:
            return f" 幾何ETA:約{eta:.0f}分"
        return f" 幾何ETA:約{eta/60.0:.1f}時間"

    def _clim_text(area: str, android_mode: bool) -> str:
        """Rev14.1: LEARNING時の参考値(②自己ログ/③NICT統計による季節平均)。
        あくまでRidge予測の代替ではなく参考値である点を明示する。

        Android版は行長超過による折り返しを避けるため、幾何ETAと同じ行には
        詰め込まず別行に分離し、ラベルも最短表記(自/NIC/参)にする
        (Pixel8a等の狭幅ターミナルでの可読性を優先)。
        """
        c = area_climatology.get(area)
        if c is None:
            return ""
        if android_mode:
            short = {"self_log": "自", "nict": "NIC", "self_log_thin": "参"}.get(c['source'], "?")
            return f"季節({short}):{c['prob_pct']:.0f}%"
        label = {"self_log": "季節平均(自己ログ)", "nict": "季節平均(NICT統計)",
                 "self_log_thin": "季節平均(参考・サンプル少)"}.get(c['source'], "季節平均")
        return f"       {label}: {c['prob_pct']:.0f}%  [{c['detail']}]"

    if android:
        active = area_df[(area_df["heard_count"] > 0) | (area_df.get("is_compensated", False) == True)]
        if active.empty:
            print("伝播: 受信局なし")
            return
        print(f"エリア状況 ({len(active)}/{len(area_df)}エリア)")
        bar_w = 5
        any_low = False
        any_comp = False
        for _, row in active.iterrows():
            area = row["area"]; name = row["area_name"]
            ratio = row["reach_ratio"]; heard = int(row["heard_count"])
            comp_ratio = row.get("reach_ratio_compensated", ratio)
            is_comp = bool(row.get("is_compensated", False))
            comp_basis = row.get("compensation_basis", "raw")
            sym = reach_trend_symbol(area)
            low_mark = "*" if row.get("area_confidence") == "LOW" else ""
            any_low = any_low or bool(low_mark)
            any_comp = any_comp or is_comp
            if pd.isna(ratio) and not pd.isna(comp_ratio):
                pct = comp_ratio * 100.0
                bar = bar_graph(pct, width=bar_w)
                print(f"{name:<4}{bar}{pct:>4.0f}%(推定) {sym}")
            elif pd.isna(ratio):
                bar = bar_graph(0, width=bar_w)
                print(f"{name:<4}{bar}  {heard}局 {sym}")
            elif comp_basis == "idw_blend" and not pd.isna(comp_ratio):
                pct = comp_ratio * 100.0
                raw_pct = ratio * 100.0
                bar = bar_graph(pct, width=bar_w)
                print(f"{name:<4}{bar}{pct:>4.0f}%{low_mark}(実測{raw_pct:.0f}%) {sym}")
            else:
                pct = ratio * 100.0
                bar = bar_graph(pct, width=bar_w)
                print(f"{name:<4}{bar}{pct:>4.0f}%{low_mark} {sym}")
            fvec = area_fvecs.get(area)
            r2 = _area_ridge_models.get(area,(None,None,None,0.0))[3] if area in _area_ridge_models else 0.0
            conf = area_pred_confidence(area, r2)
            eta_s = _eta_text(area)
            if conf == 'LEARNING':
                clim_s = _clim_text(area, android_mode=True)
                print(f"     2サイクル後: 学習中 "
                      f"総:{_area_total_samples}/{AREA_PRED_LEARN_START_TOTAL} "
                      f"Es:{_area_es_samples}/{AREA_PRED_LEARN_START_ES}")
                if eta_s or clim_s:
                    print(f"     {eta_s.strip()} {clim_s}".rstrip())
            elif fvec:
                pred_c, pred_r, bias = predict_area_reach(area, fvec)
                if not pd.isna(pred_c):
                    bias_s = f"+{bias*100:.0f}%" if bias >= 0 else f"{bias*100:.0f}%"
                    unc = (_area_holdout_diag.get(area, {}) or {}).get('resid_std')
                    unc_s = f" ±{unc*100:.0f}%" if unc is not None else ""
                    print(f"     予測:{pred_c*100:.0f}%{unc_s} Bias:{bias_s} [{conf}]{eta_s}")
        if any_low:
            print(f"*固定局<{MIN_FIXED_STATIONS_FOR_CONFIDENCE}局(低信頼)")
        if any_comp:
            print("(推定)=固定局不足を周辺エリアからIDW補間")
        print("▲上昇 →横ばい ▼下降")
        if _ver113_unlock:
            print(f"Type: {ver113_es_type_label()}")
        return
    W = DISP_W_NORMAL
    bar_w = 10
    print("エリア状況")
    any_low = False
    any_comp = False
    for _, row in area_df.iterrows():
        name  = row["area_name"]; ratio = row["reach_ratio"]
        heard = int(row["heard_count"]); area  = row["area"]
        comp_ratio = row.get("reach_ratio_compensated", ratio)
        is_comp = bool(row.get("is_compensated", False))
        comp_basis = row.get("compensation_basis", "raw")
        n_sources = int(row.get("compensation_n_sources", 0) or 0)
        sym = reach_trend_symbol(area)
        low_mark = row.get("area_confidence") == "LOW"
        any_low = any_low or low_mark
        any_comp = any_comp or is_comp
        conf_suffix = f"  低信頼(N={int(row['fixed_count'])})" if low_mark else ""
        if pd.isna(ratio) and not pd.isna(comp_ratio):
            pct = comp_ratio * 100.0
            comp_tag = f"  [推定:近傍{n_sources}エリアよりIDW補間]"
            print(f"{name:<6} {bar_graph(pct, width=bar_w)} {pct:>3.0f}%(推定) {sym}{comp_tag}")
        elif pd.isna(ratio):
            print(f"{name:<6} {bar_graph(0, width=bar_w)} {heard:>3}局(絶対数) {sym}  [推定不可:近傍データ不足]")
        elif comp_basis == "idw_blend" and not pd.isna(comp_ratio):
            pct = comp_ratio * 100.0
            raw_pct = ratio * 100.0
            blend_tag = f"  [実測{raw_pct:.0f}%を近傍{n_sources}エリアでブレンド補正]"
            print(f"{name:<6} {bar_graph(pct, width=bar_w)} {pct:>3.0f}% {sym}{conf_suffix}{blend_tag}")
        else:
            pct = ratio * 100.0
            print(f"{name:<6} {bar_graph(pct, width=bar_w)} {pct:>3.0f}% {sym}{conf_suffix}")
        fvec = area_fvecs.get(area)
        r2 = _area_ridge_models.get(area,(None,None,None,0.0))[3] if area in _area_ridge_models else 0.0
        conf = area_pred_confidence(area, r2)
        eta_s = _eta_text(area)
        if conf == 'LEARNING':
            print(f"       2サイクル後: 学習中 "
                  f"(総:{_area_total_samples}/{AREA_PRED_LEARN_START_TOTAL} "
                  f"Es:{_area_es_samples}/{AREA_PRED_LEARN_START_ES}){eta_s}")
            clim_line = _clim_text(area, android_mode=False)
            if clim_line:
                print(clim_line)
        elif fvec:
            pred_c, _, bias = predict_area_reach(area, fvec)
            if not pd.isna(pred_c):
                bias_s = f"+{bias*100:.0f}%" if bias >= 0 else f"{bias*100:.0f}%"
                unc = (_area_holdout_diag.get(area, {}) or {}).get('resid_std')
                unc_s = f" ±{unc*100:.0f}%" if unc is not None else ""
                print(f"       2サイクル後予測: {pred_c*100:.0f}%{unc_s}  "
                      f"Bias:{bias_s}  信頼度:{conf}{eta_s}")
    if any_low:
        print(f"低信頼 = 固定局登録数が{MIN_FIXED_STATIONS_FOR_CONFIDENCE}局未満")
    if any_comp:
        print("[推定]表示 = Rev14.0のIDW空間補間(観測点が無い/少ないエリアの補完値)")
    print("幾何ETA = 学習(Ridge)に依らない、現在のEs移動速度・方向からの直線外挿の目安")
    print("季節平均 = 学習中(LEARNING)のエリアに表示される参考値。自己ログ実績を優先し、")
    print("           不足時はNICT電離層月報の統計(直近数年)にフォールバックする")
    print("▲上昇傾向  →横ばい  ▼下降傾向")
    if _ver113_unlock:
        print(f"Esタイプ: {ver113_es_type_label()}")
    print("-" * W)


def print_consistency_panel(ft8_feats: Dict, android: bool = False) -> None:
    """整合性判定（仕様書10章: NICT観測値とFT8実伝播の一致確認）
    Ver11.0J改: 通信障害(psk_ok=False)と0局受信を区別して表示する。
    """
    match = ft8_feats.get('consistency_match', False)
    psk_ok = ft8_feats.get('psk_ok', True)
    psk_err = ft8_feats.get('psk_err', '')
    if android:
        if not psk_ok:
            print(f"整合性 通信障害({psk_err}) →判定不可")
            return
        nict = ft8_feats.get('nict_label', '不明')
        ft8 = ft8_feats.get('ft8_label', '不明')
        verdict = "一致" if match else "不一致"
        print(f"整合性 NICT={nict} FT8={ft8} {verdict}({fmt(ft8_feats.get('consistency'), 2)})")
        if not match and not pd.isna(ft8_feats.get('consistency', np.nan)):
            print(" →未観測支配要因あり")
        return
    if not psk_ok:
        print("整合性判定")
        print(f"FT8通信障害: {psk_err}")
        print("判定不可（NICT観測のみで継続）")
        return
    print("整合性判定")
    print()
    print("NICT")
    print(ft8_feats.get('nict_label', '不明'))
    print()
    print("FT8")
    print(ft8_feats.get('ft8_label', '不明'))
    print()
    print("判定")
    print("一致" if match else "不一致")
    print()
    print(f"Consistency  {fmt(ft8_feats.get('consistency'), 2)}")
    if not match and not pd.isna(ft8_feats.get('consistency', np.nan)):
        print("未観測支配要因あり")

def print_ft8_stopped_notice(android: bool = False) -> None:
    """FT8補正停止（仕様書11章: 最終受信確認>15分）"""
    if android:
        print("FT8補正停止(信頼度LOW・NICT主体)")
        return
    print("FT8補正停止")
    print("信頼度 LOW")
    print("予測はNICT主体")



def run_cycle(session_history: pd.DataFrame, stage: int, detail: bool, logger: SessionLogger,
              profile: Dict, dry: bool=False, debug: bool=False, url: str=NICT_FXES_URL_DEFAULT,
              wind_file: Path=DEFAULT_WIND_PROFILE_CSV, force_no_wind: bool=False,
              mstid_enabled: bool=True, mstid_image_dir: Path=DEFAULT_MSTID_DIR,
              mstid_km_per_px: float=8.0, mstid_frame_min: float=15.0,
              mstid_km_per_px_lat: Optional[float]=None, mstid_km_per_px_lon: Optional[float]=None,
              use_color: bool=True, mock: bool=False,
              psk_enabled: bool=True, psk_callsign: str=DEFAULT_PSK_CALLSIGN,
              fixed_monitor_csv: Path=DEFAULT_FIXED_MONITOR_CSV,
              psk_lookback_sec: Optional[int]=None, cycle_interval_sec: int=DEFAULT_INTERVAL_SEC):
    global _run_cycle_mock_mode
    _run_cycle_mock_mode = bool(mock)
    if mock:
        latest_sample = build_mock_latest_sample(0)
        session_history, new_rows = merge_history_in_memory(session_history, latest_sample)
        parsed_rows = len(latest_sample); latest_rows = len(latest_sample)
    else:
        html = fetch_nict_once(url)
        parsed = parse_fxes_html(html, debug=debug)
        latest_sample = latest_observation_only(parsed)
        session_history, new_rows = merge_history_in_memory(session_history, latest_sample)
        parsed_rows = len(parsed); latest_rows = len(latest_sample)
        try:
            fof2_df = fetch_nict_foF2()
            update_fof2_history(fof2_df)
        except Exception:
            pass

    if mstid_enabled:
        if not mock:
            fetch_mstid_image_once(Path(mstid_image_dir))
        mstid_feats = analyze_mstid_images(Path(mstid_image_dir), km_per_px=mstid_km_per_px,
                                            frame_min=mstid_frame_min,
                                            km_per_px_lat=mstid_km_per_px_lat,
                                            km_per_px_lon=mstid_km_per_px_lon)
        mstid_feats['mstid_approach_index'] = compute_approach_index(mstid_feats, es_axis_dir_deg=assumed_es_axis_screen_deg())
    else:
        mstid_feats = {
            'mstid_dir_deg': np.nan, 'mstid_speed_kmh': np.nan, 'mstid_dist_km': np.nan,
            'mstid_strength': np.nan, 'mstid_cx_px': np.nan, 'mstid_cy_px': np.nan,
            'mstid_approach_index': np.nan,
        }

    fixed_db = load_fixed_monitor_db(Path(fixed_monitor_csv))
    area_df = pd.DataFrame()
    _psk_callsign_set = bool(psk_callsign and str(psk_callsign).strip())
    if psk_enabled and not mock and not _psk_callsign_set:
        psk_spots = pd.DataFrame(); psk_ok = False; psk_err = 'コールサイン未設定'
        update_psk_freshness(False)
        psk_spots_enriched = psk_spots
        area_df = compute_area_reachability(psk_spots_enriched, fixed_db)
        area_df = interpolate_area_reachability_idw(area_df)
        update_area_reach_history(area_df, now_jst())
        rpi = np.nan
        dist_feats = {'es_hop_ratio': np.nan, 'dist_in_window_count': np.nan,
                      'distance_std': np.nan, 'bearing_std': np.nan, 'dist_peak_count': np.nan}
        es_hop_ratio = np.nan
        heard_total = 0
        fresh_min = float('inf'); ft8_active = False
        globals()['_ui_last_psk_points'] = []
    elif psk_enabled:
        if mock:
            psk_spots = build_mock_psk_spots(0)
            psk_ok = True; psk_err = ''
        else:
            if is_ft8_bg_active():
                _bg_spots, _bg_ok, _bg_err, _bg_age = get_ft8_bg_latest()
                if _bg_ok and _bg_age <= FT8_BG_STALE_SEC:
                    psk_spots, psk_ok, psk_err = _bg_spots, True, ''
                else:
                    lookback = psk_lookback_sec if psk_lookback_sec else max(int(cycle_interval_sec) + 120, 300)
                    psk_spots, psk_ok, psk_err = fetch_pskreporter_spots(psk_callsign, lookback_sec=lookback)
            else:
                lookback = psk_lookback_sec if psk_lookback_sec else max(int(cycle_interval_sec) + 120, 300)
                psk_spots, psk_ok, psk_err = fetch_pskreporter_spots(psk_callsign, lookback_sec=lookback)
        update_psk_freshness(psk_ok)
        psk_spots_enriched = enrich_psk_spots_with_geometry(psk_spots) if not psk_spots.empty else psk_spots
        try:
            if not psk_spots_enriched.empty and 'distance_km' in psk_spots_enriched.columns:
                _pts = psk_spots_enriched[['distance_km', 'bearing_deg']].dropna()
                if len(_pts) > 80:
                    _pts = _pts.sample(80, random_state=0)
                globals()['_ui_last_psk_points'] = _pts.values.tolist()
            else:
                globals()['_ui_last_psk_points'] = []
        except Exception:
            globals()['_ui_last_psk_points'] = []
        area_df = compute_area_reachability(psk_spots_enriched, fixed_db)
        area_df = interpolate_area_reachability_idw(area_df)
        update_area_reach_history(area_df, now_jst())
        rpi = compute_rpi(area_df)
        dist_feats = compute_distance_distribution_features(psk_spots_enriched)
        es_hop_ratio = dist_feats['es_hop_ratio']
        heard_total = int(area_df["heard_count"].sum()) if not area_df.empty else 0
        fresh_min = psk_freshness_minutes()
        ft8_active = fresh_min <= FT8_STALE_MIN_THRESHOLD
    else:
        psk_spots = pd.DataFrame(); psk_ok = False; psk_err = 'disabled'
        rpi = np.nan; es_hop_ratio = np.nan; heard_total = 0
        dist_feats = {'es_hop_ratio': np.nan, 'dist_in_window_count': np.nan,
                      'distance_std': np.nan, 'bearing_std': np.nan, 'dist_peak_count': np.nan}
        fresh_min = float('inf'); ft8_active = False

    masking_risk_repr, net_trend, net_accel, fx_repr = estimate_network_masking_and_trend(session_history)
    mode_code, mode_label, mode_sub, es_mode_score = judge_propagation_mode(masking_risk_repr, es_hop_ratio, net_accel)
    nict_label, ft8_label, consistency, consistency_match = compute_consistency(fx_repr, rpi)

    global _global_cycle_idx
    _global_cycle_idx += 1
    cur_cycle = _global_cycle_idx

    global _latest_sw_features
    if not mock:
        try:
            sw_raw_feats = fetch_space_weather_features()
        except Exception:
            sw_raw_feats = {}
    else:
        sw_raw_feats = {
            "bz_gsm_delta_60m": -2.5, "kp_latest": 2.0,
            "speed_delta_30m": 15.0, "bz_southward_streak_min": 5.0,
        }
    sw_score_val = compute_sw_score(sw_raw_feats)
    update_sw_history(sw_raw_feats, sw_score_val)
    _latest_sw_features = sw_raw_feats

    try:
        purge_old_cache()
    except Exception:
        pass

    _fetch_ionogram_this_cycle = (cur_cycle % 5 == 1)
    dyn_for_phase = add_dynamics(session_history) if not session_history.empty else pd.DataFrame()
    latest_sample_for_phase = latest_network_sample(dyn_for_phase) if not dyn_for_phase.empty else pd.DataFrame()

    if not latest_sample_for_phase.empty:
        for _, row in latest_sample_for_phase.iterrows():
            st = row["station"]
            fxEs_st = float(row["fxEs"]) if not pd.isna(row.get("fxEs")) else np.nan
            eqi_v, eqi_c_v, ionogram_f = fetch_and_update_eqi(
                st, fxEs_st, foEs=np.nan,
                fetch_image=(not mock and _fetch_ionogram_this_cycle)
            )
            _latest_eqi_by_station[st] = {
                "eqi": eqi_v, "eqi_confidence": eqi_c_v, **ionogram_f
            }
        update_es_phase_all(latest_sample_for_phase, sw_score_val, dyn_for_phase)

    prev_spatial = _es_spatial_history[-1] if _es_spatial_history else None
    spatial = compute_es_spatial_features(
        psk_spots_enriched if psk_enabled and not psk_spots.empty else pd.DataFrame(),
        prev_spatial=prev_spatial,
        dt_min=(cycle_interval_sec / 60.0) if cycle_interval_sec else 15.0,
    )
    if psk_enabled and not pd.isna(spatial.get('es_area_center_lat', np.nan)):
        _es_spatial_history.append({**spatial, 'cycle_idx': cur_cycle,
                                     'obs_dt': now_jst().isoformat()})
        if len(_es_spatial_history) > 2000:
            _es_spatial_history[:] = _es_spatial_history[-2000:]

    _clat_early = spatial.get('es_area_center_lat', np.nan)
    _clon_early = spatial.get('es_area_center_lon', np.nan)
    area_centroid_feats = compute_area_centroid_features(_clat_early, _clon_early)
    area_geo_eta = compute_geometric_eta_all(area_centroid_feats, spatial)

    if psk_enabled and not dry:
        _log_dir_for_resolve = logger.log_dir if logger and hasattr(logger, 'log_dir') else DEFAULT_DATA_DIR
        resolve_pending_samples(area_df, cur_cycle, log_dir=_log_dir_for_resolve)

    area_fvecs: Dict[str, List[float]] = {}
    if psk_enabled and area_df is not None and not area_df.empty:
        ft8_feats_for_sample = {
            'rpi_norm': rpi/100.0 if not pd.isna(rpi) else 0.0,
            'es_mode_score': es_mode_score,
        }
        for _, row in area_df.iterrows():
            area_code = row['area']
            ratio = row.get('reach_ratio_compensated', row['reach_ratio'])
            if pd.isna(ratio):
                continue
            hist = _area_reach_history.get(area_code, [])
            fvec = _area_feature_vector(area_code, float(ratio), hist,
                                        ft8_feats_for_sample, spatial, fx_repr,
                                        centroid_feats=area_centroid_feats.get(area_code))
            area_fvecs[area_code] = fvec
        if not dry:
            register_pending_sample(area_df, ft8_feats_for_sample, spatial, fx_repr, cur_cycle,
                                     centroid_feats_by_area=area_centroid_feats)

    log_dir_early = logger.log_dir if logger and hasattr(logger, 'log_dir') else DEFAULT_DATA_DIR

    area_climatology: Dict[str, Dict] = {}
    if psk_enabled and area_df is not None and not area_df.empty:
        self_clim_for_cycle = get_self_climatology(log_dir_early)
        nict_clim_df_for_cycle = get_nict_climatology_df()
        _now_for_clim = now_jst()
        for _area_code in area_df['area'].tolist():
            _prior = get_climatology_prior(_area_code, _now_for_clim,
                                            self_clim_for_cycle, nict_clim_df_for_cycle)
            if _prior is not None:
                area_climatology[_area_code] = _prior

    if psk_enabled and not dry:
        train_interval = check_maturity_and_drift(log_dir_early, cur_cycle)
        if cur_cycle % train_interval == 0:
            train_area_predictors()

    model_ver, model_sc, ver113_active = update_model_version()

    log_dir = logger.log_dir if logger and hasattr(logger, 'log_dir') else DEFAULT_DATA_DIR
    teacher_df: pd.DataFrame = pd.DataFrame()
    if psk_enabled and not dry:
        if area_df is not None and not area_df.empty:
            save_edfs_area_history(area_df, {'rpi': rpi}, now_jst(), log_dir)
        save_edfs_centroid_history(spatial, now_jst(), log_dir)
        if cur_cycle % 10 == 0:
            teacher_df, _ = generate_future_teachers(log_dir)
    update_edfs_all(log_dir, teacher_df)

    nict_update_str = '--'
    if session_history is not None and not session_history.empty:
        try:
            nict_update_str = pd.to_datetime(session_history['obs_jst'], errors='coerce').max().strftime('%H:%M')
        except Exception:
            pass
    psk_update_str = _last_psk_success_dt.strftime('%H:%M') if (_last_psk_success_dt and psk_enabled) else '--'
    confidence_label = ('LOW' if (psk_enabled and not psk_ok)
                        else ft8_confidence_label(fresh_min, fixed_db.empty, rpi, area_df) if psk_enabled
                        else 'LOW')
    low_conf_area_count = 0
    if area_df is not None and not area_df.empty and "area_confidence" in area_df.columns:
        valid_areas = area_df.dropna(subset=["reach_ratio"])
        low_conf_area_count = int((valid_areas["area_confidence"] == "LOW").sum())

    ft8_feats = {
        'rpi': rpi, 'rpi_norm': (rpi / 100.0 if not pd.isna(rpi) else np.nan),
        'es_hop_ratio': es_hop_ratio, 'es_mode_score': es_mode_score,
        'distance_std': dist_feats['distance_std'],
        'bearing_std': dist_feats['bearing_std'],
        'dist_peak_count': dist_feats['dist_peak_count'],
        'dist_in_window_count': dist_feats['dist_in_window_count'],
        'bearing_axial_r': dist_feats.get('bearing_axial_r', np.nan),
        'bearing_orientation_deg': dist_feats.get('bearing_orientation_deg', np.nan),
        'bearing_ns_alignment': dist_feats.get('bearing_ns_alignment', np.nan),
        'bearing_anisotropy_label': dist_feats.get('bearing_anisotropy_label', 'データ不足'),
        'bearing_hist16': dist_feats.get('bearing_hist16', [0] * 16),
        'heard_total': heard_total, 'fresh_min': fresh_min, 'active': ft8_active,
        'psk_ok': psk_ok, 'psk_err': psk_err,
        'mode_code': mode_code, 'mode_label': mode_label, 'mode_sub': mode_sub,
        'nict_label': nict_label, 'ft8_label': ft8_label,
        'consistency': consistency, 'consistency_match': consistency_match,
        'nict_update_str': nict_update_str, 'psk_update_str': psk_update_str,
        'confidence_label': confidence_label, 'low_conf_area_count': low_conf_area_count,
        'model_version': model_ver, 'model_score': model_sc,
        'ver113_active': ver113_active,
        'area_total_samples': _area_total_samples,
        'area_es_samples': _area_es_samples,
        'edfs_level': _edfs_level,
        'edfs_level_label': _edfs_level_label,
        'edfs_area_teachers': _edfs_area_teachers,
        'edfs_centroid_count': _edfs_centroid_history_count,
        'edfs_feature_set': _edfs_feature_set,
        'sw_score':          sw_score_val,
        'eqi_repr':          _latest_eqi_by_station.get(
                                 max(DISPLAY_ORDER,
                                     key=lambda s: float(_latest_eqi_by_station.get(s, {}).get("eqi", 0.0) or 0.0)
                                 ), {}).get("eqi", np.nan),
        'eqi_confidence':    float(np.nanmean([v.get("eqi_confidence", 0.3)
                                               for v in _latest_eqi_by_station.values()])
                                  ) if _latest_eqi_by_station else np.nan,
        'es_phase_repr':     max(
                                 (_latest_es_phase_by_station.get(st, "FORMATION") for st in DISPLAY_ORDER),
                                 key=lambda p: _ES_PHASE_NUM.get(p, 0.0), default="FORMATION"
                             ),
        'es_phase_num_repr': float(max(
                                 (_latest_es_phase_num_by_station.get(st, 0.0) for st in DISPLAY_ORDER),
                                 default=0.0
                             )),
    }

    eqi_repr_val  = ft8_feats.get("eqi_repr", np.nan)
    eqi_conf_val  = ft8_feats.get("eqi_confidence", 0.3)
    sw_score_for_conf = sw_score_val
    trace_quality_val = float(np.nanmean(
        [_latest_eqi_by_station.get(st, {}).get("trace_state", 0.3) or 0.3
         for st in DISPLAY_ORDER]
    )) if _latest_eqi_by_station else 0.3

    prediction_confidence = compute_prediction_confidence(
        profile, ft8_feats, eqi_repr_val, eqi_conf_val, sw_score_for_conf, trace_quality_val
    )
    ft8_feats["prediction_confidence"] = prediction_confidence

    directional_field = None

    if not dry:
        pf_row = {
            "timestamp": now_jst().strftime("%Y-%m-%d %H:%M:%S"),
            "space_weather_score": sw_score_val,
            "eqi":                 eqi_repr_val,
            "reflection_quality":  trace_quality_val,
            "es_phase":            ft8_feats.get("es_phase_repr", "FORMATION"),
            "prediction_confidence": prediction_confidence,
        }
        try:
            append_to_prediction_feature_csv(log_dir, pf_row)
        except Exception:
            pass

        try:
            area_forecast = build_area_forecast(area_df, area_fvecs, area_geo_eta=area_geo_eta,
                                                 area_climatology=area_climatology)
            _station_fxes = {}
            if not latest_sample_for_phase.empty:
                for _, _r in latest_sample_for_phase.iterrows():
                    _fv = _r.get('fxEs')
                    if not pd.isna(_fv):
                        _station_fxes[_r['station']] = float(_fv)
            if _station_fxes and not latest_sample_for_phase.empty:
                try:
                    _now_naive = now_jst().replace(tzinfo=None)
                    _obs_age_min = (_now_naive - latest_sample_for_phase['obs_dt'].max()
                                     ).total_seconds() / 60.0
                    if _obs_age_min > NICT_GRADIENT_MAX_AGE_MIN:
                        _station_fxes = {}
                except Exception:
                    pass
            _nict_gradient = compute_nict_es_gradient(_station_fxes) if _station_fxes else None
            if _nict_gradient is not None:
                _nict_gradient['trend'] = update_nict_gradient_trend(
                    _nict_gradient['weighted_position_frac'])
            _ft8_anisotropy = None
            if not pd.isna(spatial.get('es_principal_axis_bearing_deg', np.nan)):
                _ft8_anisotropy = {
                    'centroid_lat': spatial.get('es_area_center_lat'),
                    'centroid_lon': spatial.get('es_area_center_lon'),
                    'major_km': spatial.get('es_major_km'),
                    'minor_km': spatial.get('es_minor_km'),
                    'principal_axis_bearing_deg': spatial.get('es_principal_axis_bearing_deg'),
                    'anisotropy_ratio': spatial.get('es_anisotropy_ratio'),
                    'n_points': int(spatial.get('es_area_count', 0) or 0),
                    'source': 'ft8_observed', 'role': 'observed',
                    'is_directly_observed': True,
                }
            directional_field = combine_directional_evidence(_nict_gradient, _ft8_anisotropy)
            evidence_snapshot = build_evidence_snapshot(
                obs_dt=now_jst(), spatial=spatial, ft8_feats=ft8_feats,
                eqi_repr_val=eqi_repr_val, eqi_conf_val=eqi_conf_val,
                sw_score_val=sw_score_val, trace_quality_val=trace_quality_val,
                prediction_confidence=prediction_confidence,
                es_phase=ft8_feats.get("es_phase_repr", "FORMATION"),
                cycle_idx=cur_cycle, forecast=area_forecast,
                directional_field=directional_field,
            )
            write_evidence_json(evidence_snapshot, now_jst(), log_dir)
            write_gemma_prompt(evidence_snapshot, log_dir)
            global _LATEST_EVIDENCE_SNAPSHOT
            _LATEST_EVIDENCE_SNAPSHOT = evidence_snapshot
            start_mcp_server_once()
        except Exception:
            pass

        if cur_cycle % 10 == 0:
            teacher_row = {
                "timestamp": now_jst().strftime("%Y-%m-%d %H:%M:%S"),
                "eqi":               eqi_repr_val,
                "es_phase":          ft8_feats.get("es_phase_repr", "FORMATION"),
                "space_weather_score": sw_score_val,
                "reflection_quality":  trace_quality_val,
            }
            try:
                append_to_teacher_csv(log_dir, teacher_row)
            except Exception:
                pass

    global _ui_last_area_df, _ui_last_area_fvecs, _ui_last_area_geo_eta, _ui_last_area_climatology
    _ui_last_area_df = area_df.copy() if isinstance(area_df, pd.DataFrame) else pd.DataFrame()
    _ui_last_area_fvecs = dict(area_fvecs or {})
    _ui_last_area_geo_eta = dict(area_geo_eta or {})
    _ui_last_area_climatology = dict(area_climatology or {})

    if not is_android():
        print_header(stage, android=False)
        if psk_enabled:
            print_ft8_header_extra(ft8_feats, android=False)
    results, meta = analyze(session_history, stage, wind_file=wind_file, profile=profile,
                            force_no_wind=force_no_wind, mstid_feats=mstid_feats, ft8_feats=ft8_feats)
    meta['directional_field'] = directional_field
    if not is_android():
        if psk_enabled:
            print_propagation_mode_panel(ft8_feats, android=False)
            print_area_reachability_panel(area_df, area_fvecs=area_fvecs,
                                           area_geo_eta=area_geo_eta,
                                           area_climatology=area_climatology, android=False)
    print_unified_display(results, meta, detail=detail, use_color=use_color, android=is_android())

    eff_avg = meta.get('eff_avg', {})
    if eff_avg and not dry:
        try:
            eff_csv = logger.log_dir / "effective_contribution_history.csv"
            eff_row = {
                "datetime_jst": now_jst().strftime("%Y-%m-%d %H:%M:%S"),
                "cycles": meta.get("cycles", 0),
                "phase": meta.get("phase", ""),
            }
            eff_row.update({f"eff_{k}": v for k, v in eff_avg.items()})
            eff_df_new = pd.DataFrame([eff_row])
            if eff_csv.exists():
                eff_df_old = pd.read_csv(eff_csv)
                pd.concat([eff_df_old, eff_df_new], ignore_index=True).to_csv(eff_csv, index=False)
            else:
                eff_df_new.to_csv(eff_csv, index=False)
        except Exception:
            pass
        W = DISP_W_ANDROID if is_android() else DISP_W_NORMAL
        top4 = sorted(eff_avg.items(), key=lambda x: x[1], reverse=True)[:4]
        print("-" * W)
        best_res = max(results, key=lambda r: -1.0 if pd.isna(r.get("p27")) else float(r.get("p27", 0.0)))
        pred_pct = float(best_res.get("p27", np.nan))
        eqi_xai  = ft8_feats.get("eqi_repr", np.nan)
        phase_xai = ft8_feats.get("es_phase_repr", "---")
        sw_xai   = ft8_feats.get("sw_score", np.nan)
        conf_xai = ft8_feats.get("prediction_confidence", np.nan)
        if not pd.isna(pred_pct):
            if is_android():
                print(f"予測 {fmt(pred_pct,0)}%")
                print(f" EQI={fmt(eqi_xai,1)} Phase={phase_xai}")
                print(f" SW={fmt(sw_xai,2)} PConf={fmt(conf_xai,0)}%")
            else:
                print(f"Prediction {fmt(pred_pct,0)}%  Confidence {fmt(conf_xai,0)}%")
                print("  主要因 (内部状態ベース):")
                internal_label = {
                    "eqi": "EQI",
                    "es_mode": "Es Mode",
                    "fxEs_trend": "fxEs傾向",
                    "cluster": "広域Cluster",
                    "shear": "垂直シア",
                    "mstid_approach": "MSTID接近",
                    "rpi": "FT8 RPI",
                    "dfoF2": "ΔfoF2",
                    "es_phase_num": "Es Phase",
                    "sw_score": "SpaceWX",
                }
                bar_w = 8 if is_android() else 10
                for k, v in top4:
                    label = internal_label.get(k, k)
                    bar = bar_graph(v, width=bar_w)
                    print(f"    {label:<14} {bar} {v:4.0f}%")
                print(f"  EQI上昇    +{fmt(eqi_xai,1) if not pd.isna(eqi_xai) else '--'}  Phase={phase_xai}  SpaceWX={fmt(sw_xai,2) if not pd.isna(sw_xai) else '--'}")
        else:
            print("実効寄与:")
            bar_w = 8 if is_android() else 10
            for k, v in top4:
                bar = bar_graph(v, width=bar_w)
                print(f" {k:<14} {bar} {v:4.0f}%")

    if psk_enabled:
        W = DISP_W_ANDROID if is_android() else DISP_W_NORMAL
        print("-" * W)
        print_consistency_panel(ft8_feats, android=is_android())
        if not ft8_active:
            print("-" * W)
            print_ft8_stopped_notice(android=is_android())

    print_edfs_status(android=is_android())

    if psk_enabled and not dry and not psk_spots.empty:
        try:
            spots_csv = logger.log_dir / "ft8_spots_history.csv"
            log_df = psk_spots_enriched.copy() if not psk_spots_enriched.empty else psk_spots.copy()
            log_df.insert(0, "datetime_jst", now_jst().strftime("%Y-%m-%d %H:%M:%S"))
            keep_cols = [c for c in ["datetime_jst", "receiver_callsign", "area", "distance_km",
                                      "bearing_deg", "snr", "freq_hz"] if c in log_df.columns]
            log_df = log_df[keep_cols]
            if spots_csv.exists():
                old = pd.read_csv(spots_csv)
                pd.concat([old, log_df], ignore_index=True).to_csv(spots_csv, index=False)
            else:
                log_df.to_csv(spots_csv, index=False)
        except Exception:
            pass

    if not dry:
        logger.add_rows(meta.get('predlog_rows', []))
    print_debug_info(debug, new_rows, parsed_rows, latest_rows, url, logger, wind_file, force_no_wind, mstid_enabled)
    return session_history


def wait_with_countdown(wait_sec: int, use_color: bool=True):
    """Ver11.7: 待機中に N/D/F/E/G/Q で表示モードを切替可能。
    Ver11.8修正: android実行時はDISP_W_NORMAL(78)固定ではなく
    DISP_W_ANDROID(48)を用いるようにし、Ver11.7の表示幅拡張（30→48）を
    この待機バナーにも反映する。
    """
    wait_sec = max(1, int(wait_sec))
    W = DISP_W_ANDROID if is_android() else DISP_W_NORMAL
    try:
        remaining = wait_sec
        while remaining > 0:
            if _edfs_stop_event.is_set():
                print('\r' + ' '*W + '\r', end='', flush=True)
                raise KeyboardInterrupt
            if ES_SETTINGS_APPLY_NOW.is_set() and remaining > 30:
                print('\r' + ' '*W + '\r', end='', flush=True)
                print('[settings] 新しい設定が保存されたため、次回サイクルを早めます。')
                remaining = 30
            cmd = poll_stdin_nonblocking()
            if cmd:
                if not handle_runtime_command(cmd, use_color=use_color):
                    raise KeyboardInterrupt
            mm, ss = divmod(remaining, 60)
            msg = f"次回 {mm:02d}:{ss:02d} Mode={_ui_display_mode} N/D/F/E/G/Q"
            prefix = C_GRAY if use_color else ''
            suffix = C_RESET if use_color else ''
            print(prefix + '\r' + msg[:W] + suffix, end='', flush=True)
            time.sleep(1)
            remaining -= 1
        print('\r' + ' '*W + '\r', end='', flush=True)
    except KeyboardInterrupt:
        print('\n終了しました。')
        raise

def resolve_use_color(args) -> bool:
    if args.no_color:
        return False
    if args.color:
        return True
    return supports_ansi_color()




def choose_initial_display_mode(default_mode: str = UI_MODE_NORMAL, skip_menu: bool = False, 
                                 profile: Dict = None, use_color: bool = True) -> None:
    """Ver12.1統合版: 起動時の表示モード選択。

    - Enterのみ: 通常表示
    - N: 通常表示
    - D: 詳細表示
    - F: Factor解析
    - E: EDFS学習状態
    - G: Debug
    - --ui-mode指定時/非対話時はメニューをスキップして指定モードへ設定
    - メニュー表示時に、HFバンド別バーグラフを同時表示（Ver12.1新機能）
    """
    global _ui_display_mode
    mode_map = {
        'N': UI_MODE_NORMAL,
        'NORMAL': UI_MODE_NORMAL,
        'D': UI_MODE_DETAIL,
        'DETAIL': UI_MODE_DETAIL,
        'F': UI_MODE_FACTOR,
        'FACTOR': UI_MODE_FACTOR,
        'E': UI_MODE_EDFS,
        'EDFS': UI_MODE_EDFS,
        'G': UI_MODE_DEBUG,
        'DEBUG': UI_MODE_DEBUG,
    }
    default_mode = mode_map.get(str(default_mode).strip().upper(), UI_MODE_NORMAL)
    if skip_menu:
        _ui_display_mode = default_mode
        return
    try:
        if profile is not None:
            try:
                print("\n")
                print_hf_band_bargraph_monitor_silent(profile=profile, use_color=use_color, android=True)
            except Exception as e_bargraph:
                try:
                    print_hf_band_bargraph_monitor_simple(profile=profile, use_color=use_color, android=True)
                except Exception:
                    print(" (警告: バンド状態表示エラー)")
        
        print("\n表示モード選択")
        print("[N] 通常表示(推奨)")
        print("[D] 詳細表示")
        print("[F] Factor解析")
        print("[E] EDFS学習状態")
        print("[G] Debug")

        if _active_html_bridge is not None:
            try:
                _active_html_bridge.write_snapshot(sys.stdout.snapshot())
            except Exception:
                pass

        sel = input("\n> Enter=通常表示 : ").strip().upper()
        _ui_display_mode = mode_map.get(sel, UI_MODE_NORMAL)
    except Exception:
        _ui_display_mode = default_mode


_ANSI_RE = re.compile(r'\x1b\[([0-9;]*)m')
_ANSI_COLOR_MAP = {
    '30': '#000000', '31': '#e74c3c', '32': '#2ecc71', '33': '#f1c40f',
    '34': '#3498db', '35': '#9b59b6', '36': '#1abc9c', '37': '#dcdcdc',
    '90': '#888888', '91': '#ff6b6b', '92': '#2ecc71', '93': '#f1c40f',
    '94': '#5dade2', '95': '#c39bd3', '96': '#48c9b0', '97': '#ffffff',
}


def _ansi_line_to_html(line: str) -> str:
    """1行分のANSIエスケープ付きテキストをHTML(<span>)に変換する。
    本スクリプトが実際に出力するのは C_RESET/C_RED/C_YELLOW/C_GREEN/C_GRAY
    (\\033[0m,91,93,92,90) のみだが、念のため 30-37/90-97 系も汎用対応する。
    """
    out = []
    pos = 0
    open_span = False
    for m in _ANSI_RE.finditer(line):
        chunk = line[pos:m.start()]
        if chunk:
            out.append(_html_mod.escape(chunk))
        codes = m.group(1).split(';') if m.group(1) else ['0']
        for code in codes:
            if code in ('', '0'):
                if open_span:
                    out.append('</span>')
                    open_span = False
            elif code in _ANSI_COLOR_MAP:
                if open_span:
                    out.append('</span>')
                out.append(f'<span style="color:{_ANSI_COLOR_MAP[code]}">')
                open_span = True
        pos = m.end()
    tail = line[pos:]
    if tail:
        out.append(_html_mod.escape(tail))
    if open_span:
        out.append('</span>')
    return ''.join(out)


def ansi_text_to_html(text: str) -> str:
    """複数行テキストをHTML化する（<pre>内表示前提、改行はそのまま保持）"""
    return '\n'.join(_ansi_line_to_html(ln) for ln in text.split('\n'))


class _TeeStdout:
    """標準出力をコンソール表示と内部バッファへ同時に流すラッパー。
    print()の見た目・タイミングは一切変更せず、ミラーするだけ。

    ITPFS Ver1.0(修正): DUCT予測エンジンはEs/F2本体サイクルと完全に別スレッドで
    動作するため、両者のprint()が同時に発生し得る。EDFS本体のwait_with_countdown()は
    `\\r`で同じ行を毎秒上書きする方式(末尾に改行を付けない)のため、そこへ他スレッド
    (DUCT/PSK/学習等)の通常の改行付きメッセージが割り込むと、コンソール上で
    「次回 14:59 Mode=...[wspr] クエリ失敗:...」のように2つのメッセージが同じ行で
    衝突し、DUCT側のログが読み取れない(=「起動シーケンスでDUCTが出ない」ように見える
    不具合の一因)。対策として、直前の出力が改行で終わっていない状態で、かつ今回の
    書き込みが`\\r`(カウントダウン自身の次回更新)で始まらない場合は、先に改行を1つ
    挿入してから書き込むことで、他スレッドのメッセージが必ず独立した行として
    表示されるようにする。カウントダウン自身の見た目・更新タイミングは変えない。
    """
    def __init__(self, real_stdout):
        self._real = real_stdout
        self._buf = io.StringIO()
        self._lock = threading.Lock()
        self._last_ended_with_newline = True

    def write(self, s):
        if not s:
            return 0
        with self._lock:
            if (not self._last_ended_with_newline) and (not s.startswith('\r')) \
               and (not s.startswith('\n')):
                self._real.write('\n')
                self._buf.write('\n')
            self._real.write(s)
            self._buf.write(s)
            self._last_ended_with_newline = s.endswith('\n')
        return len(s)

    def flush(self):
        try:
            self._real.flush()
        except Exception:
            pass

    def isatty(self):
        return getattr(self._real, 'isatty', lambda: False)()

    def snapshot(self) -> str:
        with self._lock:
            return self._buf.getvalue()

    def snapshot_and_reset(self) -> str:
        with self._lock:
            val = self._buf.getvalue()
            self._buf = io.StringIO()
            return val


_HTML_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>ESDUCT | Es &amp; Duct Propagation Forecast</title>
<meta http-equiv="refresh" content="{refresh_sec}">
<style>
  body {{ background:#111; color:#ddd; font-family:"Consolas","Menlo","MS Gothic",monospace; margin:0; padding:12px; }}
  .updated {{ color:#888; font-size:0.85em; margin-bottom:8px; }}
  pre {{ white-space: pre-wrap; word-break: break-word; font-size:14px; line-height:1.35; }}
  .copyright {{ color:#666; font-size:0.7em; margin-top:16px; white-space:pre-wrap; line-height:1.4; }}
</style>
</head>
<body>
<div class="updated">最終更新: {updated_at}（{refresh_sec}秒毎に自動更新）</div>
<pre>{body}</pre>
<div class="copyright">Copyright (c) 2026 JL7KHN/栃技研
[利用・免責条件] 個人的利用における本ソースコードの改変・利用は自由です。本プログラムの結果および使用に伴ういかなる不利益・損害についても、作成者（JL7KHN/栃技研）は一切の責任を負いません。</div>
</body>
</html>
"""


def _get_lan_ip() -> str:
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.settimeout(0.5)
        s.connect(('8.8.8.8', 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return '127.0.0.1'


def _check_reachable(host: str, port: int, timeout: float = 1.0) -> bool:
    """Rev15.9: 指定host:portへTCP接続できるかを簡易確認する。
    同一プロセス(コンテナ)内からの自己接続確認に過ぎず、ChromeOS側の
    ブラウザから実際に届くかどうかまでは保証しないが、明らかに誤った
    IPアドレス(バインドされていないインターフェース等)を弾くには十分。
    """
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(timeout)
        result = s.connect_ex((host, port))
        s.close()
        return result == 0
    except Exception:
        return False


def _pick_dashboard_url(port: int, filename: str) -> Tuple[str, str, bool]:
    """Rev15.9: ダッシュボードの「実際に開くべきURL」を1箇所で決定する。

    従来の問題点: 表示上はlan_ip(またはChromebook ARCブリッジIP)ベースの
    URLを案内しつつ、自動起動(_try_auto_open_browser)には常にハードコードの
    127.0.0.1版URLを渡していたため、Chromebook(DASHBOARD_BIND_HOST="0.0.0.0")
    環境では「手動で正しいURLを開けば表示できるが、自動起動は常に失敗/誤動作する」
    という状態になっていた。本関数は表示・自動起動の両方で同一のURL選定ロジックを
    共有し、この不整合を解消する。

    戻り値: (推奨URL, 選定理由ラベル, 自己接続で到達確認できたか)
    """
    if DASHBOARD_BIND_HOST == "127.0.0.1":
        url = f"http://127.0.0.1:{port}/{filename}"
        return url, "127.0.0.1(既定)", True
    lan_ip = _get_lan_ip()
    if lan_ip != '127.0.0.1' and _check_reachable(lan_ip, port):
        return f"http://{lan_ip}:{port}/{filename}", f"LAN/ARC実IP({lan_ip})", True
    if lan_ip != '127.0.0.1':
        return f"http://{lan_ip}:{port}/{filename}", f"LAN/ARC実IP({lan_ip}、到達未確認)", False
    return f"http://127.0.0.1:{port}/{filename}", "127.0.0.1(フォールバック)", False


_active_html_bridge = None


def _run_am_start(args_list, timeout=6):
    """`am start ...`を1回実行し、(成功したか, 診断文字列)を返す。

    Rev16.5(バグ修正・重要): 従来はos.system()の戻り値(シェルの終了コード)しか
    見ておらず、なぜ失敗したか(コマンド自体が無い/権限エラー/バックグラウンド
    起動制限等)が一切わからないまま全ての手段が失敗して「ブラウザの自動起動には
    対応していません」とだけ表示されていた。subprocess.run()でstdout/stderrを
    捕捉し、失敗時に具体的な理由を診断情報として残せるようにする。
    """
    try:
        import subprocess as _sp
        proc = _sp.run(args_list, capture_output=True, text=True, timeout=timeout)
        out = (proc.stdout or '').strip()
        err = (proc.stderr or '').strip()
        detail = (err or out or '(出力なし)')[:200]
        ok = (proc.returncode == 0) and ('error' not in detail.lower()) \
            and ('exception' not in detail.lower())
        return ok, f"rc={proc.returncode} {detail}"
    except FileNotFoundError:
        return False, "amコマンドが見つかりません(PATH未設定/非対応環境)"
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"


def _android_open_url_via_am(url: str):
    """Android標準の`am`(activity manager)コマンドでACTION_VIEW Intentを
    発行し、ブラウザでURLを開く。

    【採用理由】Pydroid3のようなAndroidアプリ内で動くPythonからは、
    Python標準のwebbrowserモジュールが機能しない(Android非対応)ことに加え、
    pyjnius経由でIntentを発行するにはActivityコンテキストが必要だが、
    Pydroid3はKivy/python-for-android由来の`org.kivy.android.PythonActivity`を
    持たないため、その定石も使えない。
    一方、`am start -a android.intent.action.VIEW -d <url>`はAndroid標準の
    shellコマンドであり、Pydroid3アプリ内から直接呼び出せる(Activityコンテキスト
    取得が不要)。これはPydroid3で実際に動作確認されている既知の方法
    (pre-3.0のwebbrowser代替として案内されている手法)。

    Rev16.5(バグ修正・重要): Rev16.2で追加した`--user 0`フラグが、一部の端末
    (シングルユーザー構成でも、OSバージョン・ROMによっては非rootアプリからの
    `--user`指定にINTERACT_ACROSS_USERS系権限を要求し、SecurityExceptionで
    コマンド自体が失敗するケースがある)で、かえって「以前は動いていたのに動かなく
    なる」退行を招いていた可能性が高いと判断し、フラグ無し版を先に試すよう変更。
    また、os.system()をsubprocess.run()ベースの_run_am_start()に置き換え、
    各試行の成否・エラー内容を診断リストとして返すようにした
    (呼び出し元_try_auto_open_browser()が全滅時に詳細を表示するため)。
    試行順序: ①パッケージ指定なし(--user無し) ②Chrome指定(--user無し)
    ③パッケージ指定なし(--user 0付き) ④Chrome指定(--user 0付き)
    ①②を先にすることで、`--user`フラグが問題になる環境でも最初の1〜2回で
    成功する可能性を高める(Rev16.2の「WebViewハンドラ誤爆」対策自体は、
    ②/④でChrome指定を明示的に試すことで維持している)。
    """
    safe_url = url.replace('"', '')
    diags = []
    attempts = [
        ("generic", ['am', 'start', '-a', 'android.intent.action.VIEW', '-d', safe_url]),
        ("chrome", ['am', 'start', '-a', 'android.intent.action.VIEW',
                    '-p', 'com.android.chrome', '-d', safe_url]),
        ("generic--user0", ['am', 'start', '--user', '0', '-a',
                             'android.intent.action.VIEW', '-d', safe_url]),
        ("chrome--user0", ['am', 'start', '--user', '0', '-a', 'android.intent.action.VIEW',
                           '-p', 'com.android.chrome', '-d', safe_url]),
    ]
    for label, cmd in attempts:
        ok, detail = _run_am_start(cmd)
        diags.append(f"am({label}): {detail}")
        if ok:
            return True, diags
    return False, diags


def _android_open_url_via_jnius(url: str):
    """`am`コマンドが使えない/失敗する機種向けのpyjnius直接呼び出しフォールバック。
    org.kivy.android.PythonActivity(Pydroid3には存在しない)には依存せず、
    android.app.ActivityThreadから直接ApplicationContextを取得することで、
    Activityクラスを持たないアプリ(Pydroid3等)でもIntent発行を試みる。

    Rev16.2: _android_open_url_via_am()と同じ理由で、まずGoogle Chrome本体
    (com.android.chrome)を明示指定して起動を試み、失敗時のみパッケージ未指定に
    フォールバックする(WebViewベースの簡易ハンドラに渡ってsecure context
    判定が誤ることを防ぐ)。
    Rev16.5: 失敗理由を診断リストとして返すよう変更(呼び出し元で表示するため)。
    """
    diags = []
    try:
        from jnius import autoclass
        Intent = autoclass('android.content.Intent')
        Uri = autoclass('android.net.Uri')
        ActivityThread = autoclass('android.app.ActivityThread')
        context = ActivityThread.currentApplication().getApplicationContext()

        intent = Intent(Intent.ACTION_VIEW, Uri.parse(url))
        intent.addFlags(Intent.FLAG_ACTIVITY_NEW_TASK)
        try:
            intent.setPackage('com.android.chrome')
            context.startActivity(intent)
            return True, diags
        except Exception as e:
            diags.append(f"jnius(chrome): {type(e).__name__}: {e}")
        intent2 = Intent(Intent.ACTION_VIEW, Uri.parse(url))
        intent2.addFlags(Intent.FLAG_ACTIVITY_NEW_TASK)
        context.startActivity(intent2)
        return True, diags
    except Exception as e:
        diags.append(f"jnius(generic): {type(e).__name__}: {e}")
        return False, diags


def _try_auto_open_browser(url: str) -> str:
    """可能ならブラウザを自動で開く。

    Rev15.9: 従来はtermux-open-url(Termux専用)とwebbrowser.open()
    (公式に"Availability: not Android"でPydroid3等では機能しない)しか
    試しておらず、Pydroid3上ではほぼ確実に自動起動できていなかった。
    Android(Pydroid3含む)では`am`コマンドによるACTION_VIEW Intent発行を
    最優先で試すよう変更し、これがAndroid端末・Chromebook(ARC)の双方で
    実際にブラウザを起動できる主要な手段となる。

    Rev16.5(バグ修正・重要): 全ての方式が失敗した場合、従来は理由が一切
    分からない「対応していません」の1行のみだった。各試行の診断情報を集約し、
    失敗時にはコンソール(≒HTMLダッシュボードのログにもミラーされる)へ
    具体的な理由(コマンド不在/権限エラー/例外内容等)を出力するようにした。
    """
    import shutil
    import subprocess
    all_diags = []
    if is_android():
        ok, diags = _android_open_url_via_am(url)
        all_diags.extend(diags)
        if ok:
            return 'am-intent'
        ok, diags = _android_open_url_via_jnius(url)
        all_diags.extend(diags)
        if ok:
            return 'jnius-intent'
    if shutil.which('termux-open-url'):
        try:
            subprocess.Popen(
                ['termux-open-url', url],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            return 'termux-open-url'
        except Exception as e:
            all_diags.append(f"termux-open-url: {type(e).__name__}: {e}")
    if shutil.which('xdg-open'):
        try:
            subprocess.Popen(
                ['xdg-open', url],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            return 'xdg-open'
        except Exception as e:
            all_diags.append(f"xdg-open: {type(e).__name__}: {e}")
    try:
        import webbrowser
        if webbrowser.open(url):
            return 'webbrowser'
        all_diags.append("webbrowser.open(): False(対応ブラウザが見つかりませんでした)")
    except Exception as e:
        all_diags.append(f"webbrowser: {type(e).__name__}: {e}")
    if all_diags:
        print("  ブラウザ自動起動 診断詳細(いずれも失敗しました):")
        for d in all_diags:
            print(f"    - {d}")
    return ''


class _V13ReusableTCPServer(socketserver.ThreadingTCPServer):
    """Rev16.4(バグ修正・重要): allow_reuse_address を既定で True にする。
    socketserver.ThreadingTCPServer は既定で allow_reuse_address=False であり、
    Pydroid3等でスクリプトを「正常終了(Ctrl+C/常駐停止ボタン)」ではなく強制終了・
    アプリ切替・OSによるプロセスkillで止めた場合、直前のTCPソケットがTIME_WAIT状態の
    まま残ることがある。この状態でスクリプトを再起動すると、bind()が
    OSError(Address already in use)を出し続け、既定の20ポート探索(8765-8784)を
    全て使い果たしてHTTPサーバーが一切起動できなくなる。
    この場合ブラウザ自動起動(_try_auto_open_browser)も実行されず、ユーザーが
    手動でHTMLファイルをfile://として開くと、JS側のfetch('/api/status')が
    ネットワークエラーで失敗し続け、ダッシュボードが「接続待ち…」5%のまま
    永久に進まなくなる(この経路で実際に報告された症状と一致)。
    allow_reuse_address=Trueにすることで、TIME_WAIT中のソケットが残っていても
    同じポートへ即座に再バインドできるようにし、この経路での恒久停止を防ぐ。
    """
    daemon_threads = True
    allow_reuse_address = True


class HtmlDashboardBridge:
    """EDFSのコンソール表示(GUI)を自動的にHTMLへミラーする。
    ファイル書き出しと簡易HTTP配信を常時両方行う。
    学習・ログ保存処理には一切関与しない完全な後付けミラーであり、
    HTML/HTTP出力の有無やタイミングに関わらず学習データが失われることはない。
    """
    def __init__(self, html_path: Path, port: Optional[int] = 8765,
                 refresh_sec: int = 60, enable_http: bool = True):
        self.html_path = Path(html_path)
        self.port = port
        self.refresh_sec = max(5, refresh_sec)
        self.enable_http = enable_http and (port is not None)
        self._httpd = None
        self._thread = None

    def write_snapshot(self, text: str):
        body_html = ansi_text_to_html(text) if text.strip() else "(起動中… まだ出力はありません)"
        page = _HTML_PAGE_TEMPLATE.format(
            refresh_sec=self.refresh_sec,
            updated_at=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            body=body_html,
        )
        tmp_path = self.html_path.with_suffix(self.html_path.suffix + '.tmp')
        try:
            self.html_path.parent.mkdir(parents=True, exist_ok=True)
            tmp_path.write_text(page, encoding='utf-8')
            os.replace(tmp_path, self.html_path)
        except Exception as e:
            if not getattr(self, '_warned', False):
                self._warned = True
                print(f"  警告: HTMLダッシュボード書き出しに失敗しました "
                      f"({type(e).__name__}: {e}) / 出力先={self.html_path.resolve()}")

    def start_http_server(self) -> Optional[int]:
        if not self.enable_http or self._httpd is not None:
            return None
        directory = str(self.html_path.parent)
        handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=directory)
        last_err = None
        for attempt_port in range(self.port, self.port + 20):
            try:
                self._httpd = _V13ReusableTCPServer(
                    (DASHBOARD_BIND_HOST, attempt_port), handler)
                self.port = attempt_port
                break
            except OSError as e:
                last_err = e
                continue
        if self._httpd is None:
            print(f"  警告: HTTPサーバーの起動に失敗しました "
                  f"(ポート{self.port}-{self.port+19}を試行, 最終エラー: {last_err})")
            return None
        self._thread = threading.Thread(target=self._httpd.serve_forever, daemon=True)
        self._thread.start()
        return self.port

    def stop(self):
        if self._httpd is not None:
            try:
                self._httpd.shutdown()
            except Exception:
                pass


def main():
    global MSTID_AUTO_FETCH_URL, MSTID_NICT_AUTO_ENABLED, MSTID_NICT_KIND
    global DASHBOARD_BIND_HOST
    if DASHBOARD_BIND_HOST == "127.0.0.1" and is_chromebook_arc():
        DASHBOARD_BIND_HOST = "0.0.0.0"
        print("[EDFS] Chromebook(ARC)環境を検出したため、"
              "ダッシュボードの待受アドレスを0.0.0.0に自動変更しました。")
        print("[EDFS] 注意: この構成ではダッシュボードをhttp(非HTTPS)かつ"
              "localhost以外のURLで開くことになるため、ブラウザのGPS機能"
              "(Geolocation API)がsecure context制約により動作しない場合があります。"
              "その場合はダッシュボード下部「緯度/経度を手動設定」をご利用ください"
              "(環境によらず確実に動作します)。")
    ap = argparse.ArgumentParser(description=TITLE)
    ap.add_argument('--once', action='store_true', help='1回だけ取得して終了')
    ap.add_argument('--loop', action='store_true', help='15分周期で継続（互換用）')
    ap.add_argument('--interval', type=int, default=DEFAULT_INTERVAL_SEC, help='周期秒（標準900）')
    ap.add_argument('--stage', type=int, default=2, choices=[1,2,3], help='解析状態 1/2/3')
    ap.add_argument('--detail', action='store_true', help='数値背景を明示')
    ap.add_argument('--debug', action='store_true', help='デバッグ情報を末尾表示')
    ap.add_argument('--dry', action='store_true', help='ログ等を書かない')
    ap.add_argument('--url', default=NICT_FXES_URL_DEFAULT, help='取得URL')
    ap.add_argument('--color', action='store_true', help='カラー出力を強制有効')
    ap.add_argument('--no-color', action='store_true', help='カラー出力を強制無効')

    ap.add_argument('--no-log', action='store_true', help='セッションログ保存を無効化')
    ap.add_argument('--calibration', action='store_true', help='校正用ログとして保存')
    ap.add_argument('--immediate-log', action='store_true',
                     help='[Ver10.8J以降は常時有効のため指定不要・互換性のため残置]')
    ap.add_argument('--log-interval', type=int, default=1, help='Nサイクルごとに全データを上書き保存')
    ap.add_argument('--tmpfs', action='store_true', help='ログ保存先に /tmp を使用')
    ap.add_argument('--log-dir', default=str(DEFAULT_DATA_DIR), help='ログ保存ディレクトリ')
    ap.add_argument('--no-auto-cal', action='store_true', help='起動時自動補正を無効化')
    ap.add_argument('--profile-json', default=str(DEFAULT_PROFILE_JSON), help='補正プロファイルJSON')

    ap.add_argument('--wind-file', default=str(DEFAULT_WIND_PROFILE_CSV), help='風プロファイルCSV')
    ap.add_argument('--no-wind', action='store_true', help='風プロファイルを使わず pseudo 強制')
    ap.add_argument('--mstid-image-dir', default=str(DEFAULT_MSTID_DIR), help='MSTID画像ディレクトリ（自動/手動取得の保存先）')
    ap.add_argument('--mstid-km-per-px', type=float, default=8.2,
                     help='MSTID画像の km/px 近似（等方・後方互換用のフォールバック値。'
                          '--mstid-km-per-px-lat / -lon を指定した場合はそちらが優先される）')
    ap.add_argument('--mstid-km-per-px-lat', type=float, default=9.0,
                     help='Ver13.14: MSTID画像の南北(緯度)方向 km/px 較正値。'
                          '実画像368x368でy軸ラベル間隔から実測（緯度12.3px/度→約9.0km/px）')
    ap.add_argument('--mstid-km-per-px-lon', type=float, default=7.5,
                     help='Ver13.14: MSTID画像の東西(経度)方向 km/px 較正値（緯度34°N付近の概算、約7.5km/px）。'
                          '経度1°の実距離は緯度1°よりcos(緯度)分短いため、南北用の値と分けて'
                          '方向・速度を計算することで単一スカラー使用時の最大20%%程度のズレを解消する')
    ap.add_argument('--mstid-frame-min', type=float, default=10.0,
                     help='MSTID画像のフレーム間隔[分]（Ver13.14: NICT GEONET自動取得は10分間隔が既定）')
    ap.add_argument('--mstid-source', default=MSTID_NICT_KIND, choices=list(MSTID_NICT_KIND_TABLE.keys()),
                     help="Ver13.14: NICT自動取得の画像種別。H15=TEC変動成分15分以下(既定・MSTID向き) "
                          "H60=TEC変動成分60分以下 R=電子密度擾乱指数(ROTI) A=全電子数(TEC)絶対値")
    ap.add_argument('--no-nict-climatology-fetch', action='store_true',
                     help='NICT電離層月報からの季節統計(foEs気候値)自動取得を無効化(③のバックグラウンド自動更新を無効化)')
    ap.add_argument('--no-mstid-auto-fetch', action='store_true',
                     help='Ver13.14: NICT画像の自動取得を無効化し、--mstid-image-dir への手動配置のみで動作させる')
    ap.add_argument('--mstid-fetch-url', default='',
                     help='MSTID画像の自動取得元を固定URLで上書きする（画像を直接返すURL）。'
                          '指定時は NICT自動生成URLより優先される（Ver13.14）')
    ap.add_argument('--no-mstid', action='store_true', help='MSTID画像解析を無効化')

    ap.add_argument('--psk-callsign', default=DEFAULT_PSK_CALLSIGN, help='PSKReporter問い合わせ対象コールサイン（既定: JL7KHN/P）')
    ap.add_argument('--no-pskreporter', action='store_true', help='FT8/PSKReporter連携を無効化しVer10.8J相当の表示に戻す')
    ap.add_argument('--fixed-monitor-csv', default=str(DEFAULT_FIXED_MONITOR_CSV), help='固定観測局DB CSV（callsign,lat,lon,area）')
    ap.add_argument('--psk-lookback-min', type=int, default=0, help='PSKReporter取得のlookback分（0=自動: 直前取得からの経過秒+バッファ）')
    ap.add_argument('--ft8-bg', action='store_true',
                     help='Ver13.12: PSKReporter を独立スレッドで短周期ポーリング（Area Predictor 学習を約5倍加速）')
    ap.add_argument('--ft8-bg-interval', type=int, default=FT8_BG_DEFAULT_INTERVAL_SEC,
                     help=f'Ver13.12: FT8バックグラウンドポーリング周期[秒]（既定{FT8_BG_DEFAULT_INTERVAL_SEC}、最小120）')

    ap.add_argument('--no-gps', action='store_true', help='GPS自動取得を無効化')
    ap.add_argument('--my-lat', type=float, default=None, help='手動現在地 緯度（GPS/API不可時の代替）')
    ap.add_argument('--my-lon', type=float, default=None, help='手動現在地 経度（GPS/API不可時の代替）')
    ap.add_argument('--gps-file', default=None, help='Termux GPSブリッジJSONファイルパス（既定: /storage/emulated/0/EDFS/gps_location.json）')
    ap.add_argument('--gps-max-age', type=int, default=None, help='Termux GPSブリッジファイルの有効期限（分、既定10）')

    ap.add_argument('--mock', action='store_true', help='モックデータで動作確認')
    ap.add_argument('--check-progress', action='store_true',
                     help='学習進捗（有効データ件数・必要件数・除外セグメント等）のみ表示して終了（通信・監視サイクルは実行しない）')
    ap.add_argument('--ui-mode', default='NORMAL', choices=['NORMAL','DETAIL','FACTOR','EDFS','DEBUG','N','D','F','E','G'], help='起動時表示モード（Android/Pydroid3向け）')
    ap.add_argument('--no-ui-menu', action='store_true', help='起動時の表示モード選択メニューを表示しない')
    ap.add_argument('--show-bargraph', action='store_true', help='Ver12.1: Windows/RPi起動時にHFバンド別バーグラフを表示（Android/Pydroid3では常時表示）')

    ap.add_argument('--no-html', action='store_true', help='HTMLダッシュボード書き出しを無効化')
    ap.add_argument('--html-out', default=None, help='HTML出力先パス（既定: ログディレクトリ/edfs_dashboard.html）')
    ap.add_argument('--no-http-server', action='store_true', help='簡易HTTPサーバーによる配信を無効化（ファイル書き出しのみ）')
    ap.add_argument('--http-port', type=int, default=8765, help='HTMLダッシュボード配信用ポート（既定8765、使用中なら自動で+1ずつ探索）')
    ap.add_argument('--no-auto-open', action='store_true', help='起動時にブラウザを自動で開かない')

    ap.add_argument('--no-duct', action='store_true', help='DUCT予測(対流圏ダクト伝搬)を起動しない')
    args = ap.parse_args()

    if args.loop:
        print("[非推奨] --loop は Ver10.8J 以降デフォルトで常時ループ動作のため指定不要です（互換性のため無害）。")
    if args.immediate_log:
        print("[非推奨] --immediate-log は Ver10.8J 以降サイクル毎保存が標準のため指定不要です（互換性のため無害）。")

    _es_settings_startup = load_es_settings()
    _es_settings_needs_setup = (not os.path.exists(ES_SETTINGS_FILE)) or (not _es_settings_startup.get("setup_done"))
    if os.path.exists(ES_SETTINGS_FILE):
        if _es_settings_startup.get("psk_callsign"):
            args.psk_callsign = _es_settings_startup["psk_callsign"]
        args.no_pskreporter = (not _es_settings_startup.get("psk_enabled", True))
        args.ft8_bg = bool(_es_settings_startup.get("ft8_bg", False))
        args.ft8_bg_interval = int(_es_settings_startup.get("ft8_bg_interval", FT8_BG_DEFAULT_INTERVAL_SEC))
        print(f"[settings] 保存済みのEs/FT8設定を読み込みました: "
              f"コールサイン={args.psk_callsign or '(未設定)'} / "
              f"FT8連携={'有効' if not args.no_pskreporter else '無効'} / "
              f"FT8バックグラウンド={'有効' if args.ft8_bg else '無効'}")
    else:
        if args.psk_callsign == DEFAULT_PSK_CALLSIGN:
            args.psk_callsign = ""
        update_es_settings({
            "psk_callsign": args.psk_callsign,
            "psk_enabled": (not args.no_pskreporter),
            "ft8_bg": bool(args.ft8_bg),
            "ft8_bg_interval": int(args.ft8_bg_interval),
        })
    globals()['_es_settings_needs_setup'] = _es_settings_needs_setup

    globals()['_v13_startup_config'] = {
        'interval_sec': args.interval,
        'mstid_source': args.mstid_source,
        'psk_callsign': args.psk_callsign,
        'ft8_bg': bool(args.ft8_bg),
        'ui_mode': args.ui_mode,
        'html_port': args.http_port,
        'no_html': bool(args.no_html),
        'no_duct': bool(args.no_duct),
    }

    if args.mstid_fetch_url:
        MSTID_AUTO_FETCH_URL = args.mstid_fetch_url
    if args.no_mstid_auto_fetch:
        MSTID_NICT_AUTO_ENABLED = False
    global NICT_CLIMATOLOGY_AUTO_ENABLED
    if args.no_nict_climatology_fetch:
        NICT_CLIMATOLOGY_AUTO_ENABLED = False
    if args.mstid_source:
        MSTID_NICT_KIND = args.mstid_source

    global _ui_gps_enabled, _ui_termux_gps_file, _ui_termux_gps_max_age_sec
    _ui_gps_enabled = not args.no_gps
    if args.gps_file:
        _ui_termux_gps_file = args.gps_file
    if args.gps_max_age is not None:
        _ui_termux_gps_max_age_sec = max(1, args.gps_max_age) * 60
    if args.my_lat is not None and args.my_lon is not None:
        ui_set_manual_location(args.my_lat, args.my_lon, label='手動')

    use_color = resolve_use_color(args)

    html_bridge = None
    _orig_stdout = sys.stdout
    if not args.no_html:
        try:
            html_out_path = Path(args.html_out) if args.html_out else (
                (Path('/tmp') if args.tmpfs else Path(args.log_dir)) / 'edfs_dashboard.html'
            )
            html_bridge = HtmlDashboardBridge(
                html_path=html_out_path,
                port=None if args.no_http_server else args.http_port,
                refresh_sec=max(30, min(args.interval, 300)),
                enable_http=not args.no_http_server,
            )
            global _active_html_bridge
            _active_html_bridge = html_bridge
            sys.stdout = _TeeStdout(_orig_stdout)
            html_bridge.write_snapshot('起動中…\n')
            used_port = html_bridge.start_http_server()
            html_abs = html_bridge.html_path.resolve()
            if used_port:
                url, url_reason, url_reachable = _pick_dashboard_url(
                    used_port, html_bridge.html_path.name)
                print(f"HTMLダッシュボード: {html_abs}")
                if DASHBOARD_BIND_HOST == "127.0.0.1":
                    print(f"  ブラウザで閲覧: {url} （自端末内のブラウザのみ。"
                          f"Chromebook等、ブラウザがAndroidアプリと別ネットワーク"
                          f"名前空間で動く環境では127.0.0.1でも届かないことがあります）")
                else:
                    print(f"  ブラウザで閲覧: {url}  [{url_reason}]")
                    if not url_reachable:
                        print("  ※自己接続確認に失敗しています。上記URLで開けない場合は、"
                              "同一Wi-Fi内の他端末からアクセスするか、"
                              "コンソールに表示されるIPアドレスを再確認してください。")
                if not args.no_auto_open:
                    opened_via = _try_auto_open_browser(url)
                    if opened_via:
                        print(f"  ブラウザを自動起動しました（{opened_via}）")
                    else:
                        print("  ブラウザの自動起動には対応していません。"
                              "上記URLを手動で開いてください。")
            else:
                print(f"HTMLダッシュボード: {html_abs} （ファイル書き出しのみ／HTTP配信は無効）")
        except Exception as e:
            sys.stdout = _orig_stdout
            html_bridge = None
            _active_html_bridge = None
            print(f"  警告: HTMLダッシュボード機能の初期化に失敗したため無効化しました "
                  f"({type(e).__name__}: {e})")

    log_enabled = not args.no_log

    logger = SessionLogger(
        enabled=log_enabled and (not args.dry),
        calibration_mode=args.calibration,
        log_interval=args.log_interval,
        log_dir=Path(args.log_dir),
        use_tmpfs=args.tmpfs,
    )

    if not args.dry:
        try:
            run_first_launch_migration_if_needed(
                Path('/tmp') if args.tmpfs else Path(args.log_dir), use_color=use_color,
            )
        except Exception as e:
            print(f"  警告: Ver12.1初回起動処理でエラーが発生しました ({e})")

        try:
            _archived_n, _reports_n = archive_old_session_logs(
                Path('/tmp') if args.tmpfs else Path(args.log_dir),
            )
            if _archived_n or _reports_n:
                print(f"  セッションログ整理: {_archived_n}件をアーカイブに統合、"
                      f"分析専用レポート{_reports_n}件を削除しました（学習データは維持）。")
        except Exception as e:
            print(f"  警告: セッションログ整理でエラーが発生しました ({e})")

    profile_path = Path(args.profile_json)
    if args.no_auto_cal:
        profile = dict(BASE_PROFILE)
        auto_cal_msg = '自動補正: 無効'
    else:
        profile, auto_cal_msg = auto_calibrate_profile(Path('/tmp') if args.tmpfs else Path(args.log_dir), profile_path)

    if is_android():
        try:
            loc = ui_resolve_location()
            if loc:
                ui_set_manual_location(loc[0], loc[1], label=f"入力: {loc[2]}")
        except Exception as e:
            print(f"  警告: GPS入力処理でエラーが発生しました ({e})")
    
    choose_initial_display_mode(default_mode=args.ui_mode, skip_menu=args.no_ui_menu, 
                                 profile=profile, use_color=use_color)

    if html_bridge is not None:
        html_bridge.write_snapshot(sys.stdout.snapshot())

    if args.mock:
        create_mock_mstid_images(Path(args.mstid_image_dir))
        session_history = prepare_mock_session_history(6)
    else:
        session_history = empty_history()

    if args.check_progress:
        print("=" * 40)
        print(f" 学習進捗チェック  {display_version()}")
        print(f" ログディレクトリ: {args.log_dir}")
        print("=" * 40)
        print(auto_cal_msg)
        print("=" * 40)
        return

    try:
        android = is_android()
        W = DISP_W_ANDROID if android else DISP_W_NORMAL
        env_name = 'Android/Pydroid3' if android else platform.system()
        print("=" * W)
        if android:
            print(f" {SYSTEM_NAME} {VERSION_LABEL}")
            print(f" {now_jst().strftime('%Y-%m-%d %H:%M JST')}")
            print(f" OS判定 : Pydroid3 / カラー出力 {'有効' if use_color else '無効'}")
            print(" 操作   : N通常 D詳細 F寄与 E学習 GDebug Q終了")
        else:
            print(f"  {TITLE}")
            print(f"  起動: {now_jst().strftime('%Y-%m-%d %H:%M JST')}")
            print(f"  環境: {env_name}  カラー: {'有効' if use_color else '無効'}")
        print("=" * W)
        print(auto_cal_msg)
        if not android:
            print(f"対応OS想定: Windows / Raspberry Pi / Pydroid3  （現在判定: {env_name} / カラー出力: {'有効' if use_color else '無効'}）")
        global _edfs_area_teachers, _edfs_centroid_history_count
        _edfs_area_teachers, _edfs_centroid_history_count = load_edfs_csv_counts(
            Path(logger.log_dir) if hasattr(logger, 'log_dir') else DEFAULT_DATA_DIR
        )
        if _edfs_area_teachers > 0 or _edfs_centroid_history_count > 0:
            print(f"EDFSログ復元: Teachers={_edfs_area_teachers}  Centroid={_edfs_centroid_history_count}")
            try:
                _real_cent = count_real_centroid_observations(
                    Path(logger.log_dir) if hasattr(logger, 'log_dir') else DEFAULT_DATA_DIR)
                print(f"  └ 重心実観測(horizon0): {_real_cent}件  "
                      f"Motion解禁閾値: {EDFS_LV3_CENTROID_HIST}行")
            except Exception:
                pass
        try:
            _bs_restored, _bs_trained, _bs_method = bootstrap_area_history_from_csv(
                Path(logger.log_dir) if hasattr(logger, 'log_dir') else DEFAULT_DATA_DIR
            )
            if _bs_restored > 0:
                _method_jp = {'exact': '厳密(ridge_training)',
                              'legacy_approx': '近似(area_history、簡易復元)'}.get(_bs_method, _bs_method)
                print(f"Ridge学習履歴復元: Teachers={_bs_restored}  "
                      f"学習成功エリア={_bs_trained}  "
                      f"Total={_area_total_samples}  Es={_area_es_samples}  "
                      f"方式={_method_jp}")
                _ss_vals = [d['skill_score'] for d in _area_holdout_diag.values()
                            if d.get('skill_score') is not None]
                if _ss_vals:
                    print(f"  └ Skill Score(対persistence, ホールドアウト20%): "
                          f"平均={sum(_ss_vals)/len(_ss_vals):+.2f}  "
                          f"評価エリア数={len(_ss_vals)}  "
                          f"(0以下=persistenceに未勝利、参考値)")
        except Exception as _e_bs:
            print(f"  警告: Ridge学習履歴復元に失敗 ({_e_bs})")

        if not args.no_pskreporter:
            _psk_cs_set_startup = bool(args.psk_callsign and str(args.psk_callsign).strip())
            fixed_db_startup = load_fixed_monitor_db(Path(args.fixed_monitor_csv))
            builtin_used = not Path(args.fixed_monitor_csv).exists()
            db_src = "組み込みDB" if builtin_used else f"CSV({args.fixed_monitor_csv})"
            if _psk_cs_set_startup:
                print(f"FT8連携: 有効（対象={args.psk_callsign}）/ 固定局DB: {len(fixed_db_startup)}局 [{db_src}]")
            else:
                print(f"FT8連携: コールサイン未設定のため無効 / 固定局DB: {len(fixed_db_startup)}局 [{db_src}]")
                print("  → --psk-callsign オプションで自局または監視対象のコールサインを指定してください。")
                print("  → 未設定のままでは、伝播方位異方性推定・学習ステージ3(Motion Predictor)には到達できません。")
            SW_CACHE_DIR.mkdir(parents=True, exist_ok=True)
            purge_old_cache()
            if args.ft8_bg and _psk_cs_set_startup:
                try:
                    start_ft8_background_poller(args.psk_callsign,
                                                 poll_interval_sec=args.ft8_bg_interval)
                    print(f"FT8バックグラウンド: 有効（周期{max(120, args.ft8_bg_interval)}秒 / 対象={args.psk_callsign}）")
                except Exception as _e_bg:
                    print(f"  警告: FT8バックグラウンド起動失敗 ({_e_bg})")
            elif args.ft8_bg and not _psk_cs_set_startup:
                print("FT8バックグラウンド: 未起動（コールサイン未設定のためスキップ）")
            if NICT_CLIMATOLOGY_AUTO_ENABLED:
                try:
                    start_nict_climatology_background(Path(args.log_dir))
                    print(f"NICT季節統計(③): 自動取得有効(バックグラウンド、"
                          f"{NICT_CLIMATOLOGY_REFRESH_DAYS}日ごと更新)")
                except Exception as _e_clim:
                    print(f"  警告: NICT季節統計バックグラウンド起動失敗 ({_e_clim})")
            else:
                print("NICT季節統計(③): 無効(--no-nict-climatology-fetch指定)")
            print(f"{SOURCE_VERSION}: Space Weather連携有効 / イオノグラム解析有効 / 因果モデル統合")
        else:
            print("FT8連携: 無効（--no-pskreporter指定、Ver10.8J相当の表示）")

        if html_bridge is not None:
            html_bridge.write_snapshot(sys.stdout.snapshot())

        if args.once:
            _once_settings = get_es_settings()
            run_cycle(
                session_history=session_history,
                stage=args.stage,
                detail=args.detail,
                logger=logger,
                profile=profile,
                dry=args.dry,
                debug=args.debug,
                url=args.url,
                wind_file=Path(args.wind_file),
                force_no_wind=args.no_wind,
                mstid_enabled=(not args.no_mstid),
                mstid_image_dir=Path(args.mstid_image_dir),
                mstid_km_per_px=args.mstid_km_per_px,
                mstid_frame_min=args.mstid_frame_min,
                mstid_km_per_px_lat=args.mstid_km_per_px_lat,
                mstid_km_per_px_lon=args.mstid_km_per_px_lon,
                use_color=use_color,
                mock=args.mock,
                psk_enabled=bool(_once_settings.get("psk_enabled", not args.no_pskreporter)),
                psk_callsign=(_once_settings.get("psk_callsign") or args.psk_callsign),
                fixed_monitor_csv=Path(args.fixed_monitor_csv),
                psk_lookback_sec=(args.psk_lookback_min * 60 if args.psk_lookback_min > 0 else None),
                cycle_interval_sec=args.interval,
            )
            if html_bridge is not None:
                html_bridge.write_snapshot(sys.stdout.snapshot())
            return
        while True:
            if _edfs_stop_event.is_set():
                break
            if html_bridge is not None:
                sys.stdout.snapshot_and_reset()
            _live_es_settings = get_es_settings()
            live_psk_callsign = _live_es_settings.get("psk_callsign", "") or ""
            live_psk_enabled = bool(_live_es_settings.get("psk_enabled", True))
            live_ft8_bg = bool(_live_es_settings.get("ft8_bg", False))
            live_ft8_bg_interval = int(_live_es_settings.get("ft8_bg_interval", FT8_BG_DEFAULT_INTERVAL_SEC))
            if live_ft8_bg and bool(live_psk_callsign) and not is_ft8_bg_active():
                try:
                    start_ft8_background_poller(live_psk_callsign, poll_interval_sec=live_ft8_bg_interval)
                    print(f"[settings] FT8バックグラウンドを起動しました(対象={live_psk_callsign})")
                except Exception as _e_bg2:
                    print(f"  警告: FT8バックグラウンド起動失敗 ({_e_bg2})")
            ES_SETTINGS_APPLY_NOW.clear()
            cycle_start = time.time()
            try:
                session_history = run_cycle(
                    session_history=session_history,
                    stage=args.stage,
                    detail=args.detail,
                    logger=logger,
                    profile=profile,
                    dry=args.dry,
                    debug=args.debug,
                    url=args.url,
                    wind_file=Path(args.wind_file),
                    force_no_wind=args.no_wind,
                    mstid_enabled=(not args.no_mstid),
                    mstid_image_dir=Path(args.mstid_image_dir),
                    mstid_km_per_px=args.mstid_km_per_px,
                    mstid_frame_min=args.mstid_frame_min,
                    mstid_km_per_px_lat=args.mstid_km_per_px_lat,
                    mstid_km_per_px_lon=args.mstid_km_per_px_lon,
                    use_color=use_color,
                    mock=args.mock,
                    psk_enabled=live_psk_enabled,
                    psk_callsign=live_psk_callsign,
                    fixed_monitor_csv=Path(args.fixed_monitor_csv),
                    psk_lookback_sec=(args.psk_lookback_min * 60 if args.psk_lookback_min > 0 else None),
                    cycle_interval_sec=args.interval,
                )
            except Exception as e:
                print(color_by_state(f'エラー: {e}', use_color=use_color))

            if not args.dry and not args.no_auto_cal:
                try:
                    new_profile, new_auto_cal_msg = auto_calibrate_profile(
                        Path('/tmp') if args.tmpfs else Path(args.log_dir), profile_path)
                    if new_auto_cal_msg != auto_cal_msg:
                        print(color_by_state('[自動補正] 再計算しました:', use_color=use_color))
                        print(new_auto_cal_msg)
                    profile, auto_cal_msg = new_profile, new_auto_cal_msg
                except Exception as e:
                    print(f"  警告: 自動補正の再計算に失敗しました ({e})")

            if not args.dry:
                try:
                    maybe_run_scheduled_merge(
                        Path('/tmp') if args.tmpfs else Path(args.log_dir), use_color=use_color,
                    )
                except Exception as e:
                    print(f"  警告: 自動マージ処理でエラーが発生しました ({e})")

            if html_bridge is not None:
                html_bridge.write_snapshot(sys.stdout.snapshot_and_reset())

            elapsed = time.time() - cycle_start
            wait_with_countdown(max(60, int(args.interval - elapsed)), use_color=use_color)
    except KeyboardInterrupt:
        pass
    finally:
        if not args.dry:
            logger.close()
        if logger.enabled:
            print(f"\nセッションログ保存先: {logger.file_path}")
            if logger.corr_csv.exists():
                print(f"相関レポートCSV: {logger.corr_csv}")
            if logger.summary_txt.exists():
                print(f"相関要約TXT: {logger.summary_txt}")
        if profile_path.exists():
            print(f"補正プロファイル: {profile_path}")
        if html_bridge is not None:
            html_bridge.write_snapshot(sys.stdout.snapshot())
            html_bridge.stop()
            sys.stdout = _orig_stdout






import http.server as _v13_http_server
import socketserver as _v13_socketserver
import functools as _v13_functools
import io as _v13_io
import json as _v13_json
import math as _v13_math
import os as _v13_os
import re as _v13_re
import sys as _v13_sys
import threading as _v13_threading
import time as _v13_time
import concurrent.futures as _v13_futures
from collections import deque as _v13_deque
from datetime import datetime as _v13_datetime
from pathlib import Path as _v13_Path
from typing import Any as _v13_Any, Dict as _v13_Dict, List as _v13_List
from typing import Optional as _v13_Optional, Tuple as _v13_Tuple
from urllib.parse import urlparse as _v13_urlparse
from urllib.parse import parse_qs as _v13_parse_qs

_v13_snapshot_lock = _v13_threading.Lock()
_v13_snapshot: _v13_Dict[str, _v13_Any] = {
    'ready': False,
    'startup_stage': '起動中… (NICT観測データを取得しています)',
    'startup_progress_pct': 5,
}
_v13_term_buf = _v13_io.StringIO()
_v13_term_lock = _v13_threading.Lock()

_v13_cycle_state = {
    'last_cycle_end_epoch': None,
    'interval_sec': int(globals().get('DEFAULT_INTERVAL_SEC', 900)),
    'startup_epoch': _v13_time.time(),
}
_v13_cycle_lock = _v13_threading.Lock()

_v13_active_session_logger = None
_v13_active_profile_path = None
_v13_active_profile_data = None
_v13_active_log_dir = None

_v13_history_lock = _v13_threading.Lock()
_v13_history: _v13_Dict[str, _v13_Any] = {
    'timestamps': _v13_deque(maxlen=96),
    'fxes_obs':   _v13_deque(maxlen=96),
    'fx30_pred':  _v13_deque(maxlen=96),
    'p27_now':    _v13_deque(maxlen=96),
    'pending_pred': _v13_deque(maxlen=8),
}

_v13_season_state = {
    'consecutive_below_threshold': 0,
    'is_season_off': False,
    'threshold_mhz': ES_PHASE_FORMATION_FXES,
    'off_cycles_required': 4,
    'window_hours': 72.0,
    'recent_weight': 0.6,
    'tier': None,
    'recent_pct': None,
    'climatology_pct': None,
    'blended_pct': None,
    'history': _v13_deque(maxlen=2000),
}
ES_SEASON_CLIMATOLOGY_PCT = {
    1: 3.0, 2: 3.0, 3: 5.0, 4: 12.0, 5: 28.0, 6: 48.0,
    7: 58.0, 8: 48.0, 9: 18.0, 10: 8.0, 11: 4.0, 12: 3.0,
}
ES_SEASON_TIER_BOUNDS = {'off_late': 10.0, 'late_in': 30.0}
ES_SEASON_HYSTERESIS_PCT = 5.0
_v13_season_lock = _v13_threading.Lock()

_v13_pending_gps: _v13_Dict[str, _v13_Any] = {}
_v13_gps_lock = _v13_threading.Lock()


def _v13_safe(v, default=None):
    try:
        if v is None:
            return default
        f = float(v)
        if _v13_math.isnan(f) or _v13_math.isinf(f):
            return default
        return f
    except (TypeError, ValueError):
        return default


def _v13_rate_label(pct):
    if pct is None:
        return ("測定中", "state-unknown")
    if pct >= 75:
        return ("強い開通", "state-hot")
    if pct >= 50:
        return ("開通", "state-good")
    if pct >= 30:
        return ("注意", "state-warn")
    if pct >= 15:
        return ("弱兆", "state-weak")
    return ("低調", "state-low")


def _v13_fxes_label(fx):
    if fx is None:
        return ("--", "state-unknown")
    if fx >= 8.0:
        return (f"{fx:.1f} MHz — 27MHz 強開通", "state-hot")
    if fx >= 5.0:
        return (f"{fx:.1f} MHz — 27MHz 開通確実", "state-good")
    if fx >= 3.5:
        return (f"{fx:.1f} MHz — 27MHz 条件次第", "state-warn")
    if fx >= 2.5:
        return (f"{fx:.1f} MHz — まだ弱い", "state-weak")
    return (f"{fx:.1f} MHz — 通らず", "state-low")


def _v13_fxes_to_p27(fx):
    if fx is None:
        return None
    if fx <= 2.5:
        return 0.0
    if fx >= 9.0:
        return 100.0
    v = 100.0 / (1.0 + _v13_math.exp(-(fx - 5.0) * 1.4))
    return round(max(0.0, min(100.0, v)), 1)


_v13_AREA_JP = {code: format_area_label(code) for code in AREA_NAME_JP}


def _v13_project_stations(stations, area_lat, area_lon, key='fx30'):
    valid = []
    for s in stations:
        v = s.get(key)
        if v is None or s.get('lat') is None or s.get('lon') is None:
            continue
        d = _v13_math.hypot(s['lat'] - area_lat, s['lon'] - area_lon) + 0.01
        valid.append((v, 1.0 / (d ** 2)))
    if not valid:
        return None
    num = sum(v * w for v, w in valid)
    den = sum(w for _, w in valid)
    return round(num / den, 2) if den > 0 else None


def _v13_compute_area_forecast(stations, area_df, cycles, horizon_min=30):
    out = []
    area_centers = globals().get('AREA_CENTERS', {})
    cur_map, heard_map = {}, {}
    try:
        import pandas as _pd
        if isinstance(area_df, _pd.DataFrame) and not area_df.empty:
            for _, row in area_df.iterrows():
                a = str(row.get('area', ''))
                cur_map[a] = _v13_safe(row.get('reach_ratio')) or 0.0
                heard_map[a] = int(row.get('heard_count', 0) or 0)
    except Exception:
        pass
    scale = horizon_min / 30.0
    for area, (lat, lon) in area_centers.items():
        fx_now = _v13_project_stations(stations, lat, lon, 'fxEs') if cycles >= 1 else None
        fx_30 = _v13_project_stations(stations, lat, lon, 'fx30') if cycles >= 3 else None
        if fx_now is not None and fx_30 is not None:
            fx_target = fx_now + (fx_30 - fx_now) * scale
        elif fx_30 is not None:
            fx_target = fx_30 * scale + (1.0 - scale) * 4.0
        else:
            fx_target = None
        p_future = _v13_fxes_to_p27(fx_target) if fx_target is not None else None
        cur_ratio = cur_map.get(area, 0.0)
        cur_pct = round(cur_ratio * 100)
        if p_future is not None and cur_pct > 0:
            diff = p_future - cur_pct
            if diff >= 15:
                trend = ("良化", "state-good", "↑")
            elif diff <= -15:
                trend = ("悪化", "state-warn", "↓")
            else:
                trend = ("横ばい", "state-weak", "→")
        elif p_future is not None and cur_pct == 0:
            trend = (("新規開通", "state-good", "★") if p_future >= 30
                     else ("継続低調", "state-low", "…"))
        else:
            trend = ("学習中", "state-unknown", "…")
        display_pct = p_future if p_future is not None else cur_pct
        p_label, p_class = _v13_rate_label(display_pct)
        display_name = _v13_AREA_JP.get(area, f"エリア{area}")
        _diag = _area_holdout_diag.get(area, {}) or {}
        _unc = _diag.get('resid_std')
        out.append({
            'area': area, 'display_name': display_name,
            'lat': lat, 'lon': lon,
            'current_pct': cur_pct, 'future_pct': p_future,
            'display_pct': display_pct, 'fx_target': fx_target,
            'heard_count': heard_map.get(area, 0),
            'p_label': p_label, 'p_class': p_class,
            'trend_label': trend[0], 'trend_class': trend[1],
            'trend_mark': trend[2],
            'uncertainty_pct': round(_unc * 100) if _unc is not None else None,
            'skill_score': _diag.get('skill_score'),
        })
    out.sort(key=lambda x: -(x['display_pct'] or 0))
    return out


def _v13_learning_quality(progress, interval_sec):
    cycles = int(progress.get('cycles', 0) or 0)
    teachers = int(progress.get('edfs_area_teachers', 0) or 0)
    centroids = int(progress.get('edfs_centroid_count', 0) or 0)
    min_cal = int(globals().get('MIN_AUTO_CAL_ROWS', 20))
    lv2 = int(globals().get('EDFS_LV2_AREA_TEACHERS', 500))
    lv3 = int(globals().get('EDFS_LV3_CENTROID_HIST', 500))

    def eta(rem):
        if rem <= 0:
            return "達成済み ✅"
        sec = rem * interval_sec
        d, r = divmod(sec, 86400)
        h, r = divmod(r, 3600)
        m = r // 60
        if d > 0:
            return f"あと約 {int(d)}日{int(h)}時間"
        if h > 0:
            return f"あと約 {int(h)}時間{int(m)}分"
        return f"あと約 {int(m)}分"

    active_profile = globals().get('_v13_active_profile_data') or {}
    prior_calibrated = bool(active_profile.get('updated_jst')) or \
        int(active_profile.get('valid_rows_used', 0) or 0) > 0
    stage1_done = prior_calibrated or (cycles >= min_cal)
    if prior_calibrated:
        stage1_pct = 100
        stage1_remaining = 0
        stage1_desc = ('F層特徴量の Ridge 校正。前回校正済みプロファイルを引き継ぎ中'
                        + (f"（前回校正: {active_profile['updated_jst']}）"
                           if active_profile.get('updated_jst') else ''))
        stage1_eta = eta(0)
    else:
        stage1_pct = min(100, round(100 * cycles / max(1, min_cal)))
        stage1_remaining = max(0, min_cal - cycles)
        stage1_desc = 'F層特徴量の Ridge 校正。20サンプルで開始（今回が初回校正）'
        stage1_eta = eta(stage1_remaining)

    stages = [
        {'name': 'Stage 1: 自動補正解禁',
         'desc': stage1_desc,
         'current': (min_cal if prior_calibrated else cycles), 'target': min_cal,
         'pct': stage1_pct,
         'remaining': stage1_remaining,
         'eta': stage1_eta,
         'done': stage1_done},
        {'name': 'Stage 2: エリア別予測 (EDFS Lv2)',
         'desc': f'FT8実測から教師データ蓄積。{lv2}件でエリア別予測稼働',
         'current': teachers, 'target': lv2,
         'pct': min(100, round(100 * teachers / max(1, lv2))),
         'remaining': max(0, lv2 - teachers),
         'eta': eta(max(0, lv2 - teachers)),
         'eta_note': '(概算・FT8スポット数依存)',
         'done': teachers >= lv2},
        {'name': 'Stage 3: Es 移動予測 (EDFS Lv3)',
         'desc': f'Es 重心履歴。{lv3}件で移動予測稼働',
         'current': centroids, 'target': lv3,
         'pct': min(100, round(100 * centroids / max(1, lv3))),
         'remaining': max(0, lv3 - centroids),
         'eta': eta(max(0, lv3 - centroids)),
         'done': centroids >= lv3},
    ]
    overall_pct = round((stages[0]['pct'] + stages[1]['pct'] + stages[2]['pct']) / 3.0)
    if overall_pct >= 90:
        lbl, cls = "十分に成熟", "state-good"
    elif overall_pct >= 60:
        lbl, cls = "実用レベル", "state-good"
    elif overall_pct >= 30:
        lbl, cls = "学習中期", "state-warn"
    elif overall_pct >= 10:
        lbl, cls = "学習初期", "state-weak"
    else:
        lbl, cls = "学習開始直後", "state-low"
    return {'stages': stages, 'overall_pct': overall_pct,
            'overall_label': lbl, 'overall_class': cls}


def _v13_element_readiness(cycles, station_stats, motion, fxes_max, mstid_feats=None):
    """各表示要素のサイクル境界とステータスをリスト化。
       - 表示可: ✓
       - 待機中: あと N サイクル / 条件詳細
       Ver13.14修正: 従来は motion (estimate_network_motion の戻り値) に対して
       motion.get('mstid_dir_deg') を見ていたが、この dict は方向/速度/信頼度など
       Es重心の南北移動情報のみを持ち、mstid_dir_deg キー自体を含んでいなかった。
       そのためMSTID画像を何枚投入しても常に None 判定＝恒久的に未達成という
       実質的なバグになっていた。実際の値は meta['mstid_feats'] 側にあるため、
       ここでは呼び出し元から別引数として明示的に受け取る。"""
    mstid_feats = mstid_feats or {}
    def st(current_ok, need_cycles_from_now, note=""):
        if current_ok:
            return {'status': '✓', 'ready': True, 'note': note or '表示中'}
        return {'status': f'あと{max(0,need_cycles_from_now)}サイクル', 'ready': False, 'note': note}

    items = []
    n_active = int(station_stats.get('n_fxes_gt3', 0))
    items.append({
        'key': 'centroid_lat', 'label': 'Es 重心 (南北位置)',
        **st(n_active >= 3, 1 if n_active < 3 else 0,
             f'fxEs>3MHzの局 {n_active}/3'),
    })
    items.append({
        'key': 'centroid_lon', 'label': 'Es 重心 (東西位置)',
        **st(n_active >= 3, 1 if n_active < 3 else 0,
             f'fxEs>3MHzの局 {n_active}/3 (Ver13.5 追加)'),
    })
    items.append({
        'key': 'shape_simple', 'label': 'Es 形状 (簡易分類)',
        **st(cycles >= 1, 1 - cycles),
    })
    items.append({
        'key': 'shape_accurate', 'label': 'Es 形状 (高精度: fxes_std_60)',
        **st(cycles >= 4, 4 - cycles, '60分移動窓に4サイクル必要'),
    })
    has_motion = motion.get('direction') not in (None, '不明', '停滞')
    items.append({
        'key': 'motion_direction', 'label': 'Es 移動方向 (南北)',
        **st(has_motion, max(0, 2 - cycles),
             '2サイクル目から、4サイクルで安定化'),
    })
    items.append({
        'key': 'motion_speed', 'label': 'Es 移動速度',
        **st(motion.get('speed_kmh') is not None and cycles >= 2,
             max(0, 2 - cycles)),
    })
    mstid_dir_ready = _v13_safe(mstid_feats.get('mstid_dir_deg')) is not None
    items.append({
        'key': 'mstid_dir', 'label': 'MSTID 移動方向',
        'status': '✓' if mstid_dir_ready else '画像未提供',
        'ready': mstid_dir_ready,
        'note': '表示中' if mstid_dir_ready else '画像1枚以上必要 (自動取得または --mstid-image-dir)',
    })
    mstid_speed_ready = _v13_safe(mstid_feats.get('mstid_speed_kmh')) is not None
    items.append({
        'key': 'mstid_speed', 'label': 'MSTID 移動速度',
        'status': '✓' if mstid_speed_ready else ('画像1枚のみ' if mstid_dir_ready else '画像未提供'),
        'ready': mstid_speed_ready,
        'note': '表示中' if mstid_speed_ready else '画像2枚以上必要 (自動取得または --mstid-image-dir)',
    })
    items.append({
        'key': 'fx30', 'label': '30分後 fxEs 予測 (fx30)',
        **st(cycles >= 3, 3 - cycles),
    })
    items.append({
        'key': 'area_forecast', 'label': 'エリア別 開通率予測',
        **st(cycles >= 3, 3 - cycles),
    })
    with _v13_season_lock:
        tier = _v13_season_state['tier']
        blended_pct = _v13_season_state['blended_pct']
        thr = _v13_season_state['threshold_mhz']
    label, _cls, _badge_cls = _v13_season_tier_display(tier, blended_pct)
    items.append({
        'key': 'season', 'label': 'Es シーズン判定',
        'status': '✓' if tier is not None else '判定中',
        'ready': tier is not None,
        'note': f'{label} (閾値{thr:.1f}MHz、直近72h実測+当月気候学ベースラインをブレンド)',
    })
    return items


def _v13_season_tier_display(tier, blended_pct):
    """tier('in_season'/'late_season'/'off_season'/None)と活動率(%)から、
    ダッシュボード表示用のラベル・配色クラスを組み立てる。"""
    pct_txt = f"{blended_pct:.0f}%" if blended_pct is not None else "--%"
    if tier == 'in_season':
        return f"シーズン中(活動率{pct_txt})", 'state-good', 'good'
    elif tier == 'late_season':
        return f"シーズン終盤(活動率{pct_txt})", 'state-warn', 'warn'
    elif tier == 'off_season':
        return f"オフシーズン(活動率{pct_txt})", 'state-low', 'low'
    else:
        return "判定中", 'state-warn', 'warn'


def _v13_update_season_state(fxes_max):
    """4局最大 fxEs から Es シーズン状態を更新する(ITPFS Ver1.0: 直近72時間の実測
    活動率と月別気候学的ベースラインをブレンドし、ヒステリシス付き3段階で判定)。"""
    with _v13_season_lock:
        thr = _v13_season_state['threshold_mhz']
        now_epoch = _v13_time.time()
        hist = _v13_season_state['history']
        hist.append((now_epoch, fxes_max))
        window_sec = _v13_season_state['window_hours'] * 3600.0
        cutoff = now_epoch - window_sec
        while hist and hist[0][0] < cutoff:
            hist.popleft()
        valid = [v for (_, v) in hist if v is not None]
        if valid:
            hits = sum(1 for v in valid if v >= thr)
            recent_pct = 100.0 * hits / len(valid)
        else:
            recent_pct = None
        month = _v13_datetime.now().month
        clim_pct = ES_SEASON_CLIMATOLOGY_PCT.get(month, 10.0)
        w = _v13_season_state['recent_weight']
        if recent_pct is None:
            blended_pct = clim_pct
        else:
            blended_pct = w * recent_pct + (1.0 - w) * clim_pct

        prev_tier = _v13_season_state['tier']
        lo_off_late = ES_SEASON_TIER_BOUNDS['off_late'] - ES_SEASON_HYSTERESIS_PCT
        hi_off_late = ES_SEASON_TIER_BOUNDS['off_late'] + ES_SEASON_HYSTERESIS_PCT
        lo_late_in = ES_SEASON_TIER_BOUNDS['late_in'] - ES_SEASON_HYSTERESIS_PCT
        hi_late_in = ES_SEASON_TIER_BOUNDS['late_in'] + ES_SEASON_HYSTERESIS_PCT
        if prev_tier is None:
            if blended_pct >= ES_SEASON_TIER_BOUNDS['late_in']:
                tier = 'in_season'
            elif blended_pct >= ES_SEASON_TIER_BOUNDS['off_late']:
                tier = 'late_season'
            else:
                tier = 'off_season'
        else:
            tier = prev_tier
            if tier == 'off_season' and blended_pct >= hi_off_late:
                tier = 'late_season'
            elif tier == 'late_season' and blended_pct < lo_off_late:
                tier = 'off_season'
            if tier == 'late_season' and blended_pct >= hi_late_in:
                tier = 'in_season'
            elif tier == 'in_season' and blended_pct < lo_late_in:
                tier = 'late_season'

        _v13_season_state['recent_pct'] = recent_pct
        _v13_season_state['climatology_pct'] = clim_pct
        _v13_season_state['blended_pct'] = blended_pct
        _v13_season_state['tier'] = tier
        _v13_season_state['is_season_off'] = (tier == 'off_season')
        if fxes_max is not None and fxes_max >= thr:
            _v13_season_state['consecutive_below_threshold'] = 0
        else:
            _v13_season_state['consecutive_below_threshold'] += 1
        return dict(_v13_season_state)


def _v13_explain_mode(s):
    if not s or s == '---':
        return "FT8 スポット待ち"
    u = s.upper() if isinstance(s, str) else ''
    if 'ES' in u and 'HOT' in u:
        return "Es 反射で近距離が抜けている典型パターン"
    if 'ES' in u:
        return "Es 反射が主体 (1000–2500km 程度)"
    if 'F2' in u:
        return "F層反射 (長距離向き)"
    if 'GW' in u or '地上波' in s:
        return "地上波中心 (近距離のみ)"
    return "混合／判定中"


def _v13_explain_sw(v):
    if v is None: return "取得中"
    v = float(v)
    if v >= 0.7: return "静穏 (Es に好条件)"
    if v >= 0.4: return "やや擾乱"
    return "擾乱中"


def _v13_explain_eqi(v):
    if v is None: return "測定中"
    v = float(v)
    if v >= 0.7: return "反射面が滑らか (安定)"
    if v >= 0.4: return "普通"
    return "反射面が荒れている"


def _v13_explain_phase(p):
    p = (p or '').upper()
    return {
        'FORMATION': "生成期: Es が育ち始めている",
        'GROWTH': "成長期: これから強くなる",
        'PEAK': "ピーク: 最も開けている",
        'MATURE': "成熟: 安定して開通中",
        'DECAY': "崩壊期: 徐々に閉じていく",
    }.get(p, "判定中")


def _v13_polar_xy(cx, cy, radius, az_deg):
    rad = _v13_math.radians(az_deg)
    return cx + radius * _v13_math.sin(rad), cy - radius * _v13_math.cos(rad)


def _v13_destination_point(lat1_deg, lon1_deg, distance_km, bearing_deg):
    """出発点(lat1,lon1)から方位bearing_deg・距離distance_kmだけ進んだ地点の緯度経度。
    単一ホップEsは反射点がおおよそ送受信点の中間になるため、地図上に実測伝播局を
    描画する際はこの関数で相手局のおおよその位置を逆算する(観測データそのものは
    PSKReporterの距離・方位のみで、相手局の正確な緯度経度を保持していないため)。
    """
    R = 6371.0
    lat1 = _v13_math.radians(lat1_deg)
    lon1 = _v13_math.radians(lon1_deg)
    brg = _v13_math.radians(bearing_deg)
    d_r = distance_km / R
    lat2 = _v13_math.asin(_v13_math.sin(lat1) * _v13_math.cos(d_r) +
                           _v13_math.cos(lat1) * _v13_math.sin(d_r) * _v13_math.cos(brg))
    lon2 = lon1 + _v13_math.atan2(
        _v13_math.sin(brg) * _v13_math.sin(d_r) * _v13_math.cos(lat1),
        _v13_math.cos(d_r) - _v13_math.sin(lat1) * _v13_math.sin(lat2))
    return _v13_math.degrees(lat2), _v13_math.degrees(lon2)


def _v13_observation_confidence_label(n: int) -> str:
    """Ver13.22: RPI/異方性の元になったPSKReporter観測局数から、粗い観測信頼度を返す。
    「RPI=72」という数値そのものと、それが何局の観測に基づくかを分離して見せるための
    簡易ラベル（厳密な統計的信頼区間ではなく目安）。"""
    if n >= 15:
        return 'HIGH'
    elif n >= 5:
        return 'MED'
    elif n >= 1:
        return 'LOW'
    return 'NONE'


def _v13_render_rose_svg(bearing_hist16, axial_r, orientation_deg, size=220):
    """Es伝播方位のローズダイアグラム(風配図)をSVG文字列として生成する。
    等方/異方性を花びらの形で、帯の向き(特に南北軸との角度差)を赤/橙の軸線で
    一目で見えるようにする。ダッシュボードのダークテーマ配色に合わせている。
    """
    cx = cy = size / 2.0
    max_r = size * 0.36
    hist = bearing_hist16 or [0] * 16
    max_bin = max(hist) if hist else 0
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
             f'viewBox="0 0 {size} {size}">']
    for frac in (1 / 3, 2 / 3, 1.0):
        r = max_r * frac
        parts.append(f'<circle cx="{cx}" cy="{cy}" r="{r:.1f}" fill="none" '
                      f'stroke="#2a3556" stroke-width="1" stroke-dasharray="2,3"/>')
    if max_bin > 0:
        bin_width = 360.0 / len(hist)
        for i, w in enumerate(hist):
            if w <= 0:
                continue
            az0, az1 = i * bin_width, (i + 1) * bin_width
            r = max_r * (w / max_bin)
            x1, y1 = _v13_polar_xy(cx, cy, r, az0)
            x2, y2 = _v13_polar_xy(cx, cy, r, az1)
            parts.append(f'<path d="M{cx:.1f},{cy:.1f} L{x1:.1f},{y1:.1f} '
                          f'A{r:.1f},{r:.1f} 0 0 1 {x2:.1f},{y2:.1f} Z" '
                          f'fill="#5fb3ff" fill-opacity="0.55" stroke="#5fb3ff" stroke-width="0.5"/>')
    ex1, ey1 = _v13_polar_xy(cx, cy, max_r, 90); ex2, ey2 = _v13_polar_xy(cx, cy, max_r, 270)
    parts.append(f'<line x1="{ex1:.1f}" y1="{ey1:.1f}" x2="{ex2:.1f}" y2="{ey2:.1f}" '
                  f'stroke="#6b7899" stroke-width="1" stroke-dasharray="4,3"/>')
    nx1, ny1 = _v13_polar_xy(cx, cy, max_r, 0); nx2, ny2 = _v13_polar_xy(cx, cy, max_r, 180)
    parts.append(f'<line x1="{nx1:.1f}" y1="{ny1:.1f}" x2="{nx2:.1f}" y2="{ny2:.1f}" '
                  f'stroke="#9aa7c7" stroke-width="1.5"/>')
    if orientation_deg is not None and axial_r is not None and not _v13_math.isnan(axial_r) and axial_r >= 0.15:
        strength = min(axial_r / 0.6, 1.0)
        line_r = max_r * (0.5 + 0.5 * strength)
        ax1, ay1 = _v13_polar_xy(cx, cy, line_r, orientation_deg)
        ax2, ay2 = _v13_polar_xy(cx, cy, line_r, orientation_deg + 180)
        color = '#ff5252' if axial_r >= 0.40 else '#ffb020'
        parts.append(f'<line x1="{ax1:.1f}" y1="{ay1:.1f}" x2="{ax2:.1f}" y2="{ay2:.1f}" '
                      f'stroke="{color}" stroke-width="3" stroke-linecap="round"/>')
    label_r = size * 0.46
    for label, az in (('N', 0), ('E', 90), ('S', 180), ('W', 270)):
        lx, ly = _v13_polar_xy(cx, cy, label_r, az)
        parts.append(f'<text x="{lx:.1f}" y="{ly:.1f}" fill="#9aa7c7" font-size="12" '
                      f'text-anchor="middle" dominant-baseline="middle">{label}</text>')
    parts.append(f'<circle cx="{cx}" cy="{cy}" r="2" fill="#e6ecff"/>')
    parts.append('</svg>')
    return ''.join(parts)


def _v13_es_settings_snapshot() -> Dict[str, Any]:
    """ESDUCT(整理版): ダッシュボードへ渡すコールサイン/FT8設定の要約。
    /api/status(定期ポーリング)側に載せることで、初回サイクル完了前
    (起動直後)でも「設定画面を自動的に開くべきか」をフロントエンドが
    判定できるようにする。
    """
    try:
        s = get_es_settings()
        needs_setup = bool(globals().get('_es_settings_needs_setup', False)) and not s.get('setup_done')
        return {'psk_callsign': s.get('psk_callsign', ''),
                'psk_enabled': bool(s.get('psk_enabled', True)),
                'ft8_bg': bool(s.get('ft8_bg', False)),
                'ft8_bg_interval': int(s.get('ft8_bg_interval', FT8_BG_DEFAULT_INTERVAL_SEC)),
                'needs_setup': needs_setup}
    except Exception:
        return {'psk_callsign': '', 'psk_enabled': True, 'ft8_bg': False,
                'ft8_bg_interval': FT8_BG_DEFAULT_INTERVAL_SEC, 'needs_setup': False}


def _v13_build_snapshot():
    import pandas as _pd
    import numpy as _np
    results_raw = globals().get('_ui_last_results', None)
    meta_raw = globals().get('_ui_last_meta', None)
    if results_raw is None or meta_raw is None:
        with _v13_cycle_lock:
            interval_sec = int(_v13_cycle_state.get('interval_sec') or 900)
            startup_epoch = _v13_cycle_state.get('startup_epoch') or _v13_time.time()
        return {
            'ready': False,
            'startup_progress_pct': 30,
            'startup_stage': '初回サイクル実行中…',
            'updated_at': _v13_datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'cycles': 0,
            'cycle_info': {
                'interval_sec': interval_sec,
                'last_cycle_end_epoch': None,
                'startup_epoch': startup_epoch,
                'next_cycle_epoch': startup_epoch + interval_sec,
                'server_epoch': _v13_time.time(),
            },
            'season': {
                'is_season_off': False,
                'tier': None,
                'label': '判定中',
                'recent_pct': None, 'climatology_pct': None, 'blended_pct': None,
                'consecutive_below': 0,
                'required': _v13_season_state['off_cycles_required'],
            },
            'map': _v13_build_provisional_map(),
            'element_readiness': _v13_element_readiness(0, {'n_fxes_gt3': 0},
                                                       {}, None),
            'es_settings': _v13_es_settings_snapshot(),
        }

    results = results_raw or []
    meta = meta_raw or {}
    area_df = globals().get('_ui_last_area_df', _pd.DataFrame())
    if not isinstance(area_df, _pd.DataFrame):
        area_df = _pd.DataFrame()

    cycles = int(meta.get('cycles', 0))
    ft8 = meta.get('ft8_feats', {}) or {}
    motion = meta.get('motion', {}) or {}

    best = None
    try:
        best = max(results, key=lambda x: -1 if _pd.isna(x.get('p27', _np.nan))
                   else float(x.get('p27', -1)))
    except Exception:
        best = results[0] if results else {}

    stations = []
    st_display = globals().get('STATIONS', {})
    st_latlon = globals().get('NICT_STATION_LATLON', {})
    jp_to_key = {v['jp']: k for k, v in st_display.items()}
    fxes_values = []
    for r in results:
        st_key = jp_to_key.get(r.get('name', ''), None)
        latlon = st_latlon.get(st_key) if st_key else None
        p27 = _v13_safe(r.get('p27'))
        p_lbl, p_cls = _v13_rate_label(p27)
        fxes_val = _v13_safe(r.get('fxEs'))
        if fxes_val is not None:
            fxes_values.append(fxes_val)
        shape_scores = r.get('local_shape_scores') or {}
        stations.append({
            'name': r.get('name', '--'),
            '_st_key': st_key,
            'lat': latlon[0] if latlon else None,
            'lon': latlon[1] if latlon else None,
            'fxEs': fxes_val,
            'fx30': _v13_safe(r.get('fx30')) if cycles >= 3 else None,
            'p27': p27, 'p27_label': p_lbl, 'p27_class': p_cls,
            'state_raw': r.get('state', '判定不可'),
            'trend_label': trend_label(r.get('trend'))
                           if 'trend_label' in globals() else '',
            'local_shape': r.get('local_shape', '不明'),
            'local_shape_scores': {k: _v13_safe(v) for k, v in shape_scores.items()},
        })

    fxes_max = max(fxes_values) if fxes_values else None
    _v13_update_season_state(fxes_max)
    with _v13_season_lock:
        season_state = dict(_v13_season_state)

    now_p = _v13_safe(best.get('p27') if best else None)
    fx_now = _v13_safe(best.get('fxEs') if best else None)
    fx_30 = _v13_safe(best.get('fx30') if best else None)
    now_lbl, now_cls = _v13_rate_label(now_p)
    fx_lbl, fx_cls = _v13_fxes_label(fx_now)
    if fx_30 is not None and fx_now is not None:
        d = fx_30 - fx_now
        if d >= 0.5: trend_mk, trend_txt = "↑", "上昇中"
        elif d <= -0.5: trend_mk, trend_txt = "↓", "低下中"
        else: trend_mk, trend_txt = "→", "横ばい"
    else:
        trend_mk, trend_txt = "…", "学習中"

    areas_by_horizon = {}
    for h in [0, 15, 30, 45, 60, 90]:
        areas_by_horizon[str(h)] = _v13_compute_area_forecast(
            stations, area_df, cycles, horizon_min=h)

    _f2_list = globals().get('_ui_last_f2_forecast', []) or []
    _f2_by_area: Dict[str, Dict] = {}
    for _f2item in _f2_list:
        if _f2item.get('no_data'):
            continue
        _f2_by_area[_f2item['area']] = {
            'area_name': _f2item.get('area_name'),
            'distance_km': _f2item.get('distance_km'),
            'by_minutes': {h.get('minutes'): h for h in _f2item.get('horizons', [])},
        }
    for _h_key, _a_list in areas_by_horizon.items():
        _h_int = int(_h_key)
        for _a in _a_list:
            _f2info = _f2_by_area.get(_a['area'])
            _h2 = _f2info['by_minutes'].get(_h_int) if _f2info else None
            if _h2:
                _f2_lbl, _f2_cls = _v13_rate_label(_h2.get('expectancy_pct'))
                _a['f2_pct'] = _h2.get('expectancy_pct')
                _a['f2_label'] = _h2.get('label')
                _a['f2_class'] = _f2_cls
                _a['f2_time_jst'] = _h2.get('forecast_time_jst')
                _a['f2_distance_km'] = _f2info.get('distance_km')
            else:
                _a['f2_pct'] = None
                _a['f2_label'] = None
                _a['f2_class'] = 'state-unknown'
                _a['f2_time_jst'] = None
                _a['f2_distance_km'] = None

    areas = areas_by_horizon['30']

    _fvecs = globals().get('_ui_last_area_fvecs', {}) or {}
    _clim = globals().get('_ui_last_area_climatology', {}) or {}
    _geo_eta_dash = globals().get('_ui_last_area_geo_eta', {}) or {}
    for _h_list in areas_by_horizon.values():
        for _a in _h_list:
            _fv = _fvecs.get(_a['area'])
            _cov = compute_state_coverage(_a['area'], _fv) if _fv is not None else None
            _a['state_coverage_pct'] = round(_cov) if _cov is not None else None
            _a['state_coverage_method'] = 'nearest_neighbor_percentile'
            _a['model_skill'] = (_area_holdout_diag.get(_a['area'], {}) or {}).get('skill_score')
            _diag_h = _area_holdout_diag.get(_a['area'], {}) or {}
            _a['idw_skill_score'] = _diag_h.get('idw_skill_score')
            _a['observed_skill_score'] = _diag_h.get('observed_skill_score')
            _a['holdout_method'] = _diag_h.get('method')
            _a['holdout_n_folds'] = _diag_h.get('n_folds')
            _a['climatology'] = _clim.get(_a['area'])
            _geo = _geo_eta_dash.get(_a['area'], {}) or {}
            _a['eta_min_geometric'] = _geo.get('eta_min')
            _a['eta_confidence'] = _geo.get('eta_confidence', 'LOW')

    es_lat = _v13_safe(motion.get('es_centroid_lat'))
    d_lat = _v13_safe(motion.get('d_centroid_lat'))
    speed = _v13_safe(motion.get('speed_kmh'))
    direction = motion.get('direction', '不明')
    es_lon = None
    if stations:
        try:
            latlon_wt = [(s['lon'], max(0, (s['fxEs'] or 0) - 3.0))
                         for s in stations
                         if s['lon'] is not None and s['fxEs'] is not None]
            wsum = sum(w for _, w in latlon_wt)
            if wsum > 0:
                es_lon = sum(lon * w for lon, w in latlon_wt) / wsum
        except Exception:
            pass

    mstid_dir_deg = _v13_safe(best.get('mstid_dir_deg')) if best else None
    mstid_speed = _v13_safe(best.get('mstid_speed_kmh')) if best else None

    map_provisional = False
    if es_lat is None:
        best_st = None
        best_fx = -1
        for s in stations:
            if s['fxEs'] is not None and s['fxEs'] > best_fx:
                best_fx = s['fxEs']
                best_st = s
        if best_st:
            es_lat = best_st['lat']
            es_lon = best_st['lon']
        else:
            es_lat = 36.0
            es_lon = 138.0
        map_provisional = True

    loc = None
    ul = globals().get('_ui_last_location', None)
    if ul is not None:
        loc = {
            'lat': _v13_safe(ul[0]),
            'lon': _v13_safe(ul[1]),
            'label': str(ul[2]) if len(ul) > 2 else '',
            'source': str(globals().get('_ui_location_source', '')),
        }

    progress = {
        'cycles': cycles,
        'phase': meta.get('phase', ''),
        'edfs_level': globals().get('_edfs_level', 0),
        'edfs_level_label': globals().get('_edfs_level_label', ''),
        'edfs_area_teachers': globals().get('_edfs_area_teachers', 0),
        'edfs_centroid_count': globals().get('_edfs_centroid_history_count', 0),
        'model_version': meta.get('model_version', ''),
        'model_score': _v13_safe(meta.get('model_score')),
    }
    with _v13_cycle_lock:
        interval_sec = int(_v13_cycle_state.get('interval_sec') or 900)
        last_end = _v13_cycle_state.get('last_cycle_end_epoch')
        startup_epoch = _v13_cycle_state.get('startup_epoch') or _v13_time.time()
    learning_quality = _v13_learning_quality(progress, interval_sec)

    ref_epoch = last_end if last_end is not None else startup_epoch
    next_epoch = ref_epoch + interval_sec
    cycle_info = {
        'interval_sec': interval_sec,
        'last_cycle_end_epoch': last_end,
        'startup_epoch': startup_epoch,
        'next_cycle_epoch': next_epoch,
        'server_epoch': _v13_time.time(),
    }

    with _v13_history_lock:
        history = {
            'timestamps': list(_v13_history['timestamps']),
            'fxes_obs': list(_v13_history['fxes_obs']),
            'fx30_pred': list(_v13_history['fx30_pred']),
            'p27_now': list(_v13_history['p27_now']),
        }

    eff_top = []
    eff_avg = meta.get('eff_avg', {}) or {}
    if eff_avg:
        internal = {
            "eqi": "電離層品質(EQI)", "es_mode": "Es モード",
            "fxEs_trend": "fxEs 傾向", "cluster": "広域まとまり",
            "shear": "垂直シア", "mstid_approach": "MSTID 接近",
            "rpi": "FT8 到達性能", "dfoF2": "F層 ΔfoF2",
            "es_phase_num": "Es フェーズ", "sw_score": "宇宙天気",
        }
        for k, v in sorted(eff_avg.items(), key=lambda kv: kv[1], reverse=True)[:6]:
            eff_top.append({'key': k, 'label': internal.get(k, k),
                            'pct': _v13_safe(v) or 0.0})

    ft8_info = {
        'mode_label': ft8.get('mode_label', '---'),
        'mode_hint': _v13_explain_mode(ft8.get('mode_label', '')),
        'confidence_label': ft8.get('confidence_label', '---'),
        'consistency_match': bool(ft8.get('consistency_match', False)),
        'psk_ok': bool(ft8.get('psk_ok', True)),
        'no_callsign': (str(ft8.get('psk_err', '')) == 'コールサイン未設定'),
        'rpi': _v13_safe(ft8.get('rpi')),
        'sw_score': _v13_safe(ft8.get('sw_score')),
        'sw_hint': _v13_explain_sw(ft8.get('sw_score')),
        'eqi': _v13_safe(ft8.get('eqi_repr')),
        'eqi_hint': _v13_explain_eqi(ft8.get('eqi_repr')),
        'es_phase_repr': str(ft8.get('es_phase_repr', '---')),
        'phase_hint': _v13_explain_phase(ft8.get('es_phase_repr', '')),
        'prediction_confidence': _v13_safe(ft8.get('prediction_confidence')),
        'anisotropy_r': _v13_safe(ft8.get('bearing_axial_r')),
        'anisotropy_orientation_deg': _v13_safe(ft8.get('bearing_orientation_deg')),
        'anisotropy_ns_alignment': _v13_safe(ft8.get('bearing_ns_alignment')),
        'anisotropy_label': ft8.get('bearing_anisotropy_label', 'データ不足'),
        'anisotropy_rose_svg': _v13_render_rose_svg(
            ft8.get('bearing_hist16'), ft8.get('bearing_axial_r'), ft8.get('bearing_orientation_deg')),
        'observation_n': int(sum(ft8.get('bearing_hist16') or [0])),
        'observation_label': _v13_observation_confidence_label(
            int(sum(ft8.get('bearing_hist16') or [0]))),
    }

    psk_map_points = []
    _raw_pts = globals().get('_ui_last_psk_points', [])
    if loc and loc.get('lat') is not None and loc.get('lon') is not None and _raw_pts:
        for _d_km, _brg in _raw_pts:
            try:
                _plat, _plon = _v13_destination_point(loc['lat'], loc['lon'], float(_d_km), float(_brg))
                psk_map_points.append({'lat': _plat, 'lon': _plon})
            except Exception:
                continue

    _aniso_r = ft8.get('bearing_axial_r')
    _moving = (speed is not None and speed > 15) or (mstid_speed is not None and mstid_speed > 15)
    _strong_aniso = _aniso_r is not None and not (isinstance(_aniso_r, float) and _v13_math.isnan(_aniso_r)) and _aniso_r >= 0.40
    _weak_aniso = _aniso_r is not None and not (isinstance(_aniso_r, float) and _v13_math.isnan(_aniso_r)) and 0.15 <= _aniso_r < 0.40
    if _aniso_r is None or (isinstance(_aniso_r, float) and _v13_math.isnan(_aniso_r)):
        reciprocity_risk = {'level': '不明', 'class': 'state-unknown', 'note': 'データ不足のため未評価'}
    elif _strong_aniso and _moving:
        reciprocity_risk = {'level': '高', 'class': 'state-hot',
                             'note': '帯状構造が強く、移動もあるため単純な双方向性は仮定しにくい'}
    elif _strong_aniso or _moving:
        reciprocity_risk = {'level': '中', 'class': 'state-warn',
                             'note': '帯状構造あり、または移動を検出（片方のみ検出の可能性）'}
    else:
        reciprocity_risk = {'level': '低', 'class': 'state-good',
                             'note': '等方的でほぼ静的な構造（双方向性を期待しやすい）'}

    map_data = {
        'stations': [
            {'name': s['name'], 'lat': s['lat'], 'lon': s['lon'],
             'fxEs': s['fxEs'], 'fx30': s['fx30'],
             'p27': s['p27'], 'p_class': s['p27_class'],
             'local_shape': s['local_shape'],
             'local_shape_scores': s['local_shape_scores']}
            for s in stations if s['lat'] is not None
        ],
        'areas_by_horizon': areas_by_horizon,
        'es_centroid': {
            'lat': es_lat, 'lon': es_lon,
            'd_lat': d_lat,
            'direction': direction, 'speed_kmh': speed,
            'mstid_dir_deg': mstid_dir_deg,
            'mstid_speed_kmh': mstid_speed,
            'provisional': map_provisional,
            'server_epoch': _v13_time.time(),
            'anisotropy_r': _v13_safe(ft8.get('bearing_axial_r')),
            'anisotropy_orientation_deg': _v13_safe(ft8.get('bearing_orientation_deg')),
        },
        'my_location': loc,
        'psk_points': psk_map_points,
        'reciprocity_risk': reciprocity_risk,
    }

    station_stats = {'n_fxes_gt3': sum(1 for s in stations
                                       if s['fxEs'] is not None and s['fxEs'] > 3.0)}
    mstid_feats_for_readiness = meta.get('mstid_feats', {}) or {}
    readiness = _v13_element_readiness(cycles, station_stats, motion, fxes_max,
                                        mstid_feats=mstid_feats_for_readiness)

    season_label, _season_text_cls, season_class = _v13_season_tier_display(
        season_state.get('tier'), season_state.get('blended_pct'))

    return {
        'ready': True,
        'updated_at': _v13_datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'version': f"{globals().get('TITLE', 'EDFS')} + GUI {SOURCE_VERSION}",
        'startup_config': globals().get('_v13_startup_config', {}),
        'es_settings': _v13_es_settings_snapshot(),
        'cycles': cycles,
        'now': {
            'p27': now_p, 'p27_label': now_lbl, 'p27_class': now_cls,
            'state_raw': best.get('state', '判定不可') if best else '判定不可',
            'fxEs': fx_now, 'fxEs_label': fx_lbl, 'fxEs_class': fx_cls,
            'best_station': best.get('name', '--') if best else '--',
            'fxes_max_all': fxes_max,
        },
        'forecast_30m': {
            'fx30': fx_30 if cycles >= 3 else None,
            'trend_mark': trend_mk, 'trend_label': trend_txt,
        },
        'areas': areas,
        'f2_forecast': globals().get('_ui_last_f2_forecast', []),
        'stations': stations,
        'map': map_data,
        'location': loc,
        'location_status': {
            'status': globals().get('_ui_gps_status', '未試行'),
            'reason': globals().get('_ui_gps_reason', ''),
        },
        'ft8': ft8_info,
        'motion': {
            'direction': direction, 'speed_kmh': speed,
            'centroid_lat': es_lat, 'centroid_lon': es_lon,
            'mstid_dir_deg': mstid_dir_deg,
            'mstid_speed_kmh': mstid_speed,
            'provisional': map_provisional,
        },
        'progress': progress,
        'learning_quality': learning_quality,
        'cycle_info': cycle_info,
        'eff_top': eff_top,
        'history': history,
        'element_readiness': readiness,
        'season': {
            'is_season_off': season_state['is_season_off'],
            'tier': season_state.get('tier'),
            'label': season_label,
            'class': season_class,
            'recent_pct': season_state.get('recent_pct'),
            'climatology_pct': season_state.get('climatology_pct'),
            'blended_pct': season_state.get('blended_pct'),
            'threshold_mhz': season_state.get('threshold_mhz'),
            'consecutive_below': season_state['consecutive_below_threshold'],
            'required': season_state['off_cycles_required'],
            'fxes_max': fxes_max,
        },
    }


def _v13_build_provisional_map():
    """Ver13.5: 起動直後でも表示できる暫定マップデータ。
       Es 位置は日本中央 (北緯 36°, 東経 138°) の暫定値を配置。"""
    return {
        'stations': [],
        'areas_by_horizon': {'0': [], '15': [], '30': [], '45': [], '60': [], '90': []},
        'es_centroid': {
            'lat': 36.0, 'lon': 138.0, 'd_lat': None,
            'direction': '不明', 'speed_kmh': None,
            'mstid_dir_deg': None, 'mstid_speed_kmh': None,
            'provisional': True,
            'server_epoch': _v13_time.time(),
        },
        'my_location': None,
    }


def _v13_update_history(snap):
    if not snap.get('ready'):
        return
    now = _v13_time.time()
    fx_now = (snap.get('now') or {}).get('fxEs')
    fx_30 = (snap.get('forecast_30m') or {}).get('fx30')
    p_now = (snap.get('now') or {}).get('p27')
    with _v13_history_lock:
        _v13_history['timestamps'].append(now)
        _v13_history['fxes_obs'].append(fx_now)
        _v13_history['p27_now'].append(p_now)
        matched = None
        remain = []
        for tgt, val in list(_v13_history['pending_pred']):
            if abs(tgt - now) <= 900:
                matched = val
            elif tgt < now - 1800:
                pass
            else:
                remain.append((tgt, val))
        _v13_history['pending_pred'].clear()
        _v13_history['pending_pred'].extend(remain)
        _v13_history['fx30_pred'].append(matched)
        if fx_30 is not None:
            _v13_history['pending_pred'].append((now + 1800, fx_30))


_DUCT_EMBED_CSS = r"""#duct-embed-root{--bg:#0f172a;--panel:#111c33;--fg:#e6edf7;--sub:#93a4bd;--line:#26344d;}
html[data-theme="light"] #duct-embed-root{--bg:#eef2f8;--panel:#fff;--fg:#12203a;--sub:#4a5b76;--line:#cdd7e6;}
#duct-embed-root, #duct-embed-root *{box-sizing:border-box}
#duct-embed-root{margin:0;font-family:'Segoe UI','Hiragino Kaku Gothic ProN',sans-serif;
 background:var(--bg);color:var(--fg);border-radius:10px;overflow:hidden;border:1px solid var(--line)}
#duct-embed-root #hd{padding:9px 14px;background:var(--panel);border-bottom:1px solid var(--line);
 display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px}
#duct-embed-root #hd h1{margin:0;font-size:15px;font-weight:600}
#duct-embed-root .meta{font-size:11px;color:var(--sub);display:flex;gap:10px;flex-wrap:wrap;align-items:center}
#duct-embed-root .badge{background:var(--bg);border:1px solid var(--line);border-radius:11px;padding:2px 9px}
#duct-embed-root .badge.eta{border-color:#5ec962;color:#7fe0a0}
#duct-embed-root #bar{display:flex;gap:12px;align-items:center;padding:7px 14px;background:var(--panel);
 border-bottom:1px solid var(--line);font-size:12px;flex-wrap:wrap}
#duct-embed-root #map{height:min(75vh,820px);min-height:460px;width:100%;background:var(--bg);position:relative}
#duct-embed-root #duct2-placeholder{position:absolute;inset:0;display:none;align-items:center;
 justify-content:center;background:var(--bg);z-index:500;text-align:center;padding:20px}
#duct-embed-root #duct2-placeholder .duct2-ph-msg{max-width:420px;color:var(--sub);font-size:13px;line-height:1.7}
#duct-embed-root button{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:4px 11px;cursor:pointer;font-size:12px}
#duct-embed-root button:hover{border-color:#5b7bb5}
#duct-embed-root #tslider{width:230px;vertical-align:middle}
#duct-embed-root #tlabel{font-weight:600;min-width:96px;display:inline-block}
#duct-embed-root .legend{display:flex;gap:6px;align-items:center;margin-left:auto}
#duct-embed-root .grad{width:150px;height:11px;border-radius:6px;
 background:linear-gradient(90deg,#3b0f70,#3b528b,#21918c,#5ec962,#fde725)}
#duct-embed-root .leaflet-popup-content{margin:9px 11px;color:#0b1220}
#duct-embed-root #ovl{display:none;position:fixed;inset:0;background:rgba(0,0,0,0.5);z-index:1000}
#duct-embed-root #panel{position:absolute;top:60px;right:20px;width:330px;background:var(--panel);
 border:1px solid var(--line);border-radius:10px;padding:14px 16px;
 box-shadow:0 10px 40px rgba(0,0,0,0.4)}
#duct-embed-root #panel .phd{display:flex;justify-content:space-between;align-items:center;
 font-size:14px;margin-bottom:6px}
#duct-embed-root #panel .x{cursor:pointer;color:var(--sub)}
#duct-embed-root #panel label{display:flex;justify-content:space-between;align-items:center;
 font-size:12px;margin:9px 0;gap:10px}
#duct-embed-root #panel select{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:4px 8px;font-size:12px;min-width:120px}
#duct-embed-root #panel input[type=text]{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:4px 8px;font-size:12px}
#duct-embed-root #panel .note{font-size:11px;color:var(--sub);margin:2px 0 8px}
#duct-embed-root #panel .prow{display:flex;gap:8px;align-items:center;margin-top:10px}
#duct-embed-root #panel .smsg{font-size:11px;color:#5ec962}
#duct-embed-root #panel .sstat{font-size:11px;color:var(--sub);margin-top:10px;line-height:1.7;
 border-top:1px solid var(--line);padding-top:8px}
#duct-embed-root #covl{display:none;position:fixed;inset:0;background:rgba(0,0,0,0.5);z-index:1000}
#duct-embed-root #cpanel{position:absolute;top:60px;right:20px;width:310px;background:var(--panel);
 border:1px solid var(--line);border-radius:10px;padding:14px 16px;
 box-shadow:0 10px 40px rgba(0,0,0,0.4)}
#duct-embed-root #cpanel .phd{display:flex;justify-content:space-between;align-items:center;
 font-size:14px;margin-bottom:6px}
#duct-embed-root #cpanel .x{cursor:pointer;color:var(--sub)}
#duct-embed-root #cpanel label{display:block;font-size:12px;margin:10px 0;color:var(--sub)}
#duct-embed-root #cpanel .cprow{display:flex;gap:6px;align-items:center;margin-top:4px}
#duct-embed-root #cpanel input[type=text]{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:5px 8px;font-size:12px;width:130px;flex:1}
#duct-embed-root #cpanel .gpsbtn{font-size:11px;padding:5px 7px;border-radius:6px;border:1px solid var(--line);
 background:var(--bg);color:var(--fg);cursor:pointer;white-space:nowrap}
#duct-embed-root #cpanel .gpsbtn:active{background:#5ec962;color:#111}
#duct-embed-root #cpanel .prow{display:flex;gap:8px;align-items:center;margin-top:12px}
#duct-embed-root #cpanel .smsg{font-size:11px;color:#5ec962;min-height:14px;display:block;margin-top:6px}
#duct-embed-root .radar-star{font-size:22px;line-height:22px;text-align:center;color:#ffd23f;
 text-shadow:0 0 3px #000,0 0 4px #000,0 0 6px #000}
#duct-embed-root .radar-tip{font-size:11px}
#duct-embed-root .bar2{height:7px;border-radius:4px;background:var(--bg);border:1px solid var(--line);
 overflow:hidden;margin-top:3px}
#duct-embed-root .bar2 i{display:block;height:100%;background:linear-gradient(90deg,#3b528b,#5ec962)}
#duct-embed-root #rpanel{position:fixed;top:0;right:-420px;width:380px;max-width:92vw;height:100%;
 background:var(--panel);border-left:1px solid var(--line);
 box-shadow:-10px 0 40px rgba(0,0,0,0.45);z-index:1100;
 transition:right .22s ease;overflow-y:auto;padding:14px 16px 30px}
#duct-embed-root #rpanel.open{right:0}
#duct-embed-root #rpanel .phd{display:flex;justify-content:space-between;align-items:center;
 font-size:14px;margin-bottom:8px;position:sticky;top:-14px;background:var(--panel);
 padding:6px 0}
#duct-embed-root #rpanel .x{cursor:pointer;color:var(--sub);font-size:16px}
#duct-embed-root #rpanel .rp-dpp{font-size:26px;font-weight:700;color:var(--fg)}
#duct-embed-root #rpanel .rp-sub{font-size:11px;color:var(--sub)}
#duct-embed-root #rpanel .rp-sec{margin-top:10px;font-size:11px;color:var(--fg)}
#duct-embed-root #rpanel .rp-sec b{font-size:12px}
#duct-embed-root #rpanel table{width:100%;font-size:11px;margin-top:4px;border-collapse:collapse}
#duct-embed-root #rpanel th,#duct-embed-root #rpanel td{padding:2px 4px}
#duct-embed-root #rpanel .rp-charts{display:flex;gap:6px;flex-wrap:wrap;margin-top:4px}
#duct-embed-root #rpanel .rp-summary{font-size:11px;margin-top:10px;padding:8px;background:#0b1220;
 color:#cbd5e1;border-radius:6px;line-height:1.6}
html[data-theme="light"] #duct-embed-root #rpanel .rp-summary{background:#e6edf7;color:#12203a}
#duct-embed-root #hpanel{position:fixed;left:12px;bottom:12px;z-index:1100;display:none;
 background:var(--panel);border:1px solid var(--line);border-radius:8px;
 box-shadow:0 6px 24px rgba(0,0,0,0.35);padding:8px 12px;font-size:12px;
 color:var(--fg);max-width:88vw}
#duct-embed-root #hpanel .hp-row{display:flex;align-items:center;gap:6px;white-space:nowrap}
#duct-embed-root #hpanel .hp-sw{display:inline-block;width:12px;height:12px;border-radius:2px;
 border:1px solid var(--line)}
#duct-embed-root #hpanel .hp-sub{color:var(--sub);font-size:10px;margin-left:4px}
"""


_v13_HTML_HEAD = r"""<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<meta name="theme-color" content="#0b1220">
<title>ESDUCT | Es &amp; Duct Propagation Forecast (__EDFS_TITLE_VERSION__)</title>
<!-- ITPFS(整理版・大規模UI改修): DUCT詳細地図を「窓の中(iframe/別タブ)」で表示
     する方式を廃止し、ダッシュボード本体へ直接統合した。地図描画に使う
     Leafletライブラリを本体<head>で1回だけ読み込む(旧実装は別ページ内で
     個別に読み込んでいた)。 -->
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css"/>
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js"></script>
<style>
:root{
  --bg:#0b1220; --panel:#131c2e; --panel2:#1a2540; --line:#2a3556;
  --fg:#e6ecff; --sub:#9aa7c7; --accent:#5fb3ff;
  --hot:#ff5252; --good:#26c281; --warn:#ffb020; --weak:#8b7cff; --low:#6b7899;
  --sea:#0d2140; --land:#1e3a5c; --pref:#2a4a70;
  --map-land:#4a5670; --map-land-stroke:#1e2838; --map-grid:#9aa8c0;
}
/* ITPFS Ver1.0(修正): ダーク/ライト切替。DUCT単体版の地図ページ(duct_forecast_map.html)
   には元々ダークモード切替があったが、EDFS本体側のメインダッシュボードには
   一切無かった(常時ダーク固定)。html[data-theme="light"] でCSS変数を上書きし、
   ボタン+localStorage永続化(#theme-toggle)で切り替えられるようにする。
   ITPFS Ver1.0(修正2): 47都道府県マップの陸地・グリッド線はJS側でSVG属性へ
   直接ハードコード('#4a5670'等)されていたため、ライトモードへ切替えても
   その部分だけダークのまま残っていた。--map-land等のCSS変数を新設し、
   drawMap()側でgetComputedStyle経由の値を使うよう変更して連動させる。*/
html[data-theme="light"]{
  --bg:#f4f6fb; --panel:#ffffff; --panel2:#eef1f8; --line:#d7deec;
  --fg:#1b2537; --sub:#5b6b8c; --accent:#1c6fd6;
  --sea:#dbeafe; --land:#e7edf5; --pref:#cbd8ee;
  --map-land:#c3cee2; --map-land-stroke:#8a9ac0; --map-grid:#7c8bb0;
}
*{box-sizing:border-box}
html,body{margin:0;padding:0;background:var(--bg);color:var(--fg);
  font-family:-apple-system,"Hiragino Kaku Gothic ProN","Yu Gothic UI",
  "Meiryo","Noto Sans JP",system-ui,sans-serif;
  font-size:15px;line-height:1.55;-webkit-text-size-adjust:100%;}
.wrap{max-width:1000px;margin:0 auto;padding:10px 10px 80px 10px;}
h1,h2,h3{margin:0 0 6px 0;font-weight:600;letter-spacing:0.02em}
h1{font-size:1.1rem;color:var(--accent)}
h2{font-size:1.0rem;border-left:4px solid var(--accent);padding-left:8px;margin-top:14px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;
  padding:12px 14px;margin:10px 0;box-shadow:0 1px 0 rgba(255,255,255,0.02) inset}
.sub{color:var(--sub);font-size:0.85rem}
.mono{font-family:ui-monospace,Menlo,Consolas,"MS Gothic",monospace}
.badge{display:inline-block;padding:2px 8px;border-radius:999px;font-size:0.75rem;
  border:1px solid var(--line);color:var(--sub);margin-right:4px}
.badge.hot{background:#3a0d0d;color:#ffb0b0;border-color:#ff5252}
.badge.good{background:#0a3520;color:#7ee6a8;border-color:#26c281}
.badge.warn{background:#3a2708;color:#ffcf7a;border-color:#ffb020}
.badge.low{background:#1a1f30;color:#9aa7c7;border-color:#4b587a}
.badge.weak{background:#241a3a;color:#c7b8ff;border-color:#8b7cff}
#duct-card .duct-band-table{width:100%;border-collapse:collapse;font-size:0.85rem;margin-top:6px}
#duct-card .duct-band-table th,#duct-card .duct-band-table td{padding:3px 6px;border-bottom:1px solid var(--line);text-align:center}
#duct-card .duct-point-form{display:flex;flex-wrap:wrap;gap:6px;margin-top:8px;align-items:center}
#duct-card .duct-point-form input{width:9em;background:var(--panel2);border:1px solid var(--line);color:var(--fg);border-radius:6px;padding:4px 6px}
#duct-card .duct-result{margin-top:10px;padding:8px;border:1px solid var(--line);border-radius:8px;background:var(--panel2)}
.pct-big{font-size:2.2rem;font-weight:700;line-height:1.1}
.state-hot{color:var(--hot)} .state-good{color:var(--good)}
.state-warn{color:var(--warn)} .state-weak{color:var(--weak)}
.state-low{color:var(--low)} .state-unknown{color:var(--sub)}
table{width:100%;border-collapse:collapse;font-size:0.9rem}
th,td{padding:6px 4px;border-bottom:1px solid var(--line);text-align:left}
th{color:var(--sub);font-weight:500;font-size:0.78rem}
.bar{background:var(--panel2);border:1px solid var(--line);border-radius:6px;
  height:14px;overflow:hidden;position:relative}
.bar>span{display:block;height:100%;background:linear-gradient(90deg,#5fb3ff,#26c281)}
.bar.warn>span{background:linear-gradient(90deg,#ffb020,#ff5252)}
.bar.low>span{background:#4b587a}
.btn{background:var(--accent);color:#001426;border:0;border-radius:8px;
  padding:6px 12px;font-weight:600;cursor:pointer;font-size:0.9rem}
.btn.ghost{background:transparent;color:var(--accent);border:1px solid var(--accent)}
.btn:disabled{opacity:0.5;cursor:not-allowed}
details{margin-top:10px}
details>summary{cursor:pointer;color:var(--accent);font-weight:600;padding:6px 0}
pre.terminal{white-space:pre-wrap;word-break:break-word;background:#050912;
  color:#c6d0e8;border-radius:8px;padding:10px;font-size:12px;line-height:1.4;
  max-height:60vh;overflow:auto}
.startup{padding:24px;text-align:center}
.startup .pct-big{color:var(--accent);margin:8px 0}
.progress{background:var(--panel2);border-radius:6px;height:8px;overflow:hidden;
  margin:8px 0;border:1px solid var(--line)}
.progress>span{display:block;height:100%;background:linear-gradient(90deg,#5fb3ff,#26c281);
  transition:width 0.9s linear}
.tag-good{background:#0a3520;color:#7ee6a8;padding:1px 8px;border-radius:6px;font-size:0.75rem}
.tag-warn{background:#3a2708;color:#ffcf7a;padding:1px 8px;border-radius:6px;font-size:0.75rem}
.footer{color:var(--sub);font-size:0.75rem;text-align:center;margin-top:20px}
.copyright-footer{color:var(--sub);font-size:0.65rem;text-align:center;margin-top:6px;
  opacity:0.75;white-space:pre-wrap;line-height:1.4}
.small{font-size:0.8rem}.hint{color:var(--sub);font-size:0.82rem;margin-top:2px}
.kpi{display:flex;gap:10px;flex-wrap:wrap}
.kpi>div{background:var(--panel2);border:1px solid var(--line);border-radius:8px;
  padding:8px 12px;min-width:120px;flex:1}
.kpi .label{color:var(--sub);font-size:0.75rem}
.kpi .value{font-size:1.05rem;font-weight:600}
.hero-forecast{background:linear-gradient(135deg,var(--panel),var(--panel2));
  border:1px solid var(--line);border-radius:14px;padding:14px 16px;margin-bottom:12px}
.hero-forecast h1{font-size:1.15rem;margin-bottom:8px;color:var(--accent)}
.hero-metric{display:flex;align-items:baseline;gap:12px;flex-wrap:wrap}
.hero-metric .pct-big{font-size:2.6rem}
.countdown-box{display:inline-flex;align-items:center;gap:8px;background:var(--panel2);
  padding:6px 12px;border-radius:8px;border:1px solid var(--line);
  font-family:ui-monospace,monospace;margin-left:auto}
.countdown-box .val{font-weight:700;color:var(--accent);font-size:1.05rem}
.pulse-dot{display:inline-block;width:10px;height:10px;border-radius:50%;
  background:#26c281;box-shadow:0 0 8px #26c281;
  animation:pulse-blink 1s ease-in-out infinite}
@keyframes pulse-blink{
  0%,100%{opacity:1;transform:scale(1)}
  50%{opacity:0.3;transform:scale(0.7)}
}
.status-alive{color:#26c281;font-size:0.7rem;font-weight:normal;margin-left:6px}
.hero-header{display:flex;align-items:center;flex-wrap:wrap;gap:8px}
.hero-header .grow{flex:1;min-width:180px}
.map-wrap{position:relative;background:var(--sea);border-radius:10px;overflow:hidden;
  border:1px solid var(--line);margin-top:8px}
.map-wrap svg{display:block;width:100%;height:auto;max-height:80vh}
.map-wrap.paper{background:var(--sea)}
.map-wrap.paper svg{background:var(--sea)}
.map-legend{display:flex;flex-wrap:wrap;gap:8px;margin-top:6px;font-size:0.75rem;color:var(--sub)}
.map-legend span{display:inline-flex;align-items:center;gap:4px}
.map-legend .sw{display:inline-block;width:12px;height:12px;border-radius:50%}
.map-tooltip{position:absolute;background:rgba(11,20,40,0.95);border:1px solid var(--accent);
  color:var(--fg);padding:6px 10px;border-radius:6px;font-size:0.8rem;pointer-events:none;
  display:none;max-width:220px;line-height:1.35;z-index:10}
.lq-row{margin:8px 0;padding:6px 0;border-bottom:1px solid var(--line)}
.lq-row:last-child{border-bottom:none}
.lq-head{display:flex;justify-content:space-between;align-items:baseline;font-size:0.9rem}
.lq-head .name{font-weight:600}
.lq-head .cur{color:var(--sub);font-family:ui-monospace,monospace;font-size:0.85rem}
.lq-desc{font-size:0.78rem;color:var(--sub);margin:2px 0 4px 0}
.lq-eta{font-size:0.8rem;margin-top:2px}
.horizon-slider{display:flex;align-items:center;gap:10px;margin:8px 0;flex-wrap:wrap}
.horizon-slider input[type="range"]{flex:1;accent-color:var(--accent);min-width:80px}
.horizon-slider .value{font-weight:600;color:var(--accent);min-width:70px;text-align:right}
.horizon-slider .horizon-play-btn{padding:5px 10px;font-size:0.8rem;white-space:nowrap}
.chart-wrap{background:var(--panel2);border-radius:8px;padding:8px;margin-top:6px;border:1px solid var(--line)}
.chart-wrap svg{width:100%;height:auto;display:block}
.readiness-list{display:grid;grid-template-columns:1fr 1fr;gap:4px 12px;
  font-size:0.78rem;margin-top:6px}
@media (max-width:600px){.readiness-list{grid-template-columns:1fr}}
.readiness-row{display:flex;justify-content:space-between;align-items:center;
  padding:4px 8px;background:var(--panel2);border-radius:4px;gap:6px}
.readiness-row .lbl{color:var(--fg);font-size:0.78rem;flex:1}
.readiness-row .st{font-family:ui-monospace,monospace;font-size:0.75rem;white-space:nowrap}
.readiness-row.ready{border-left:3px solid var(--good)}
.readiness-row.waiting{border-left:3px solid var(--warn)}
.readiness-row.ready .st{color:var(--good)}
.readiness-row.waiting .st{color:var(--warn)}
.readiness-note{color:var(--sub);font-size:0.7rem;margin-left:6px}
.season-off-warn{background:#3a2708;border:1px solid #ffb020;
  padding:8px 12px;border-radius:8px;color:#ffcf7a;font-size:0.85rem;margin-top:8px}
.provisional-tag{background:#4b3308;color:#ffcf7a;padding:1px 6px;
  border-radius:4px;font-size:0.65rem;margin-left:4px;font-weight:600}
.page-title{font-size:1.05rem;font-weight:700;color:var(--accent);margin:4px 0 8px 0;
  display:flex;align-items:center;gap:8px;flex-wrap:wrap}
.page-title-ver{font-size:0.8rem;font-weight:600;color:var(--fg);background:var(--panel2);
  border:1px solid var(--line);border-radius:6px;padding:2px 8px}
.tier-switch{display:grid;grid-template-columns:repeat(3,1fr);gap:2px;margin:6px 0 12px 0;
  background:var(--panel2);border:1px solid var(--line);border-radius:10px;padding:4px;max-width:560px}
.tier-btn{background:transparent;border:0;color:var(--sub);padding:8px 6px;
  border-radius:8px;font-size:0.82rem;font-weight:600;cursor:pointer;white-space:nowrap}
.tier-btn.active{background:var(--accent);color:#001426}
body.mode-beginner .tier-advanced{display:none !important}
body.mode-advanced .tier-beginner-only{display:none !important}
/* ITPFS Ver1.0(修正): 「ダクト予報」タブへ行くとかんたん/くわしい表示ボタンごと
   消えてしまう(=別ページへ丸ごと遷移してしまう)という報告への対応。
   window.open()で完全に別ページへ飛ぶのをやめ、他の2タブと同じ「同一ページ内の
   表示切替」に変更する。#root(Es/F2側の内容)を隠し、#duct-ui-persistent
   (常設のDUCTカード)を目立たせるだけなので、上部のタブ切替ボタン自体は
   他の2タブと同様に常に表示されたまま残る。*/
body.mode-duct #root{display:none !important}
body:not(.mode-duct) #duct-ui-persistent{opacity:0.0; max-height:0; overflow:hidden;
  pointer-events:none; margin:0; transition:none;}
body.mode-duct #duct-ui-persistent{opacity:1; max-height:none; pointer-events:auto;}
/* ITPFS(整理版・大規模UI改修): 「詳細地図を別タブで開く」を押さないと本来欲しい
   表示にならない不具合、およびiframe内Leafletのサイズ誤初期化バグへの対応として、
   単体版DUCT(V6.2.1相当)をiframe/別タブに窓表示する方式そのものを廃止した。
   #duct-embed-rootとして本体DOMへ直接統合し、データはJSON API
   (/api/duct/celldata)経由でフェッチする。旧duct_forecast_map.html由来のCSSを
   #duct-embed-root配下にスコープして追記する(以下)。 */
__DUCT_EMBED_CSS__
</style>
</head>
<body class="mode-beginner">
<div class="wrap">
  <div class="page-title" id="page-title">📡 ESDUCT | Es &amp; Duct Propagation Forecast <span class="page-title-ver">__EDFS_TITLE_VERSION__</span>
    <button id="settings-toggle" class="btn ghost" type="button" style="margin-left:auto;padding:2px 10px;font-size:0.8rem" title="コールサイン・FT8設定">⚙️ 設定</button>
    <button id="theme-toggle" class="btn ghost" type="button" style="padding:2px 10px;font-size:0.8rem" title="ダーク/ライト表示切替">🌓</button>
  </div>
  <!-- ESDUCT(整理版): コールサイン/FT8設定モーダル。右上「⚙️ 設定」ボタン、および
       初回起動時(コールサイン未設定)に自動的に開く。#root外に常設配置することで、
       タブ切替(#rootの丸ごと再構築)の影響を受けない(#merge-ui-persistent等と
       同じ設計方針)。 -->
  <div id="settings-modal-overlay" style="display:none;position:fixed;inset:0;background:rgba(0,0,0,0.55);z-index:1000;align-items:center;justify-content:center;padding:16px">
    <div id="settings-modal-box" class="card" style="max-width:480px;width:100%;max-height:88vh;overflow-y:auto;margin:0">
      <h2 style="margin-top:0">⚙️ コールサイン・FT8設定</h2>
      <div id="settings-modal-firstrun-note" class="hint" style="display:none;color:var(--warn);margin-bottom:8px">
        初回起動時の設定です。ここでコールサインを設定すると、FT8実測データによる
        方面別到達率(RPI)・伝播方位推定・学習ステージが有効になります(未設定でも
        Es予測の基本機能はそのままご利用いただけます)。
      </div>
      <div style="margin:10px 0">
        <label class="sub" for="settings-callsign-input">コールサイン(自局または監視対象局)</label><br>
        <input id="settings-callsign-input" type="text" placeholder="例: JA1ABC またはJA1ABC/P" maxlength="32"
          style="width:100%;background:var(--panel2);border:1px solid var(--line);color:var(--fg);border-radius:8px;padding:8px 10px;font-size:1rem;margin-top:4px">
        <div class="hint">PSKReporterで実際に送信レポートが登録されているコールサインを指定してください。空欄で保存するとFT8連携は無効になります。</div>
      </div>
      <div style="margin:10px 0">
        <label style="display:flex;align-items:center;gap:8px;cursor:pointer">
          <input id="settings-ft8-enabled" type="checkbox" style="width:auto">
          <span>FT8/PSKReporter連携を有効にする</span>
        </label>
      </div>
      <div style="margin:10px 0">
        <label style="display:flex;align-items:center;gap:8px;cursor:pointer">
          <input id="settings-ft8-bg" type="checkbox" style="width:auto">
          <span>FT8バックグラウンド短周期ポーリングを有効にする(学習を加速)</span>
        </label>
        <div class="hint">通常サイクル(数分〜十数分毎)とは別に、独立スレッドで短周期にPSKReporterへ問い合わせ、エリア別予測の学習を加速します。</div>
        <div style="margin-top:6px">
          <label class="sub" for="settings-ft8-bg-interval">ポーリング周期</label><br>
          <select id="settings-ft8-bg-interval" style="background:var(--panel2);border:1px solid var(--line);color:var(--fg);border-radius:8px;padding:6px 8px;margin-top:4px">
            <option value="120">2分(120秒)</option>
            <option value="180">3分(180秒)</option>
            <option value="300">5分(300秒)</option>
            <option value="600">10分(600秒)</option>
          </select>
        </div>
      </div>
      <div id="settings-modal-msg" class="hint" style="min-height:1.2em;margin:6px 0"></div>
      <div style="display:flex;gap:8px;justify-content:flex-end;margin-top:12px">
        <button id="settings-modal-cancel" class="btn ghost" type="button">閉じる</button>
        <button id="settings-modal-save" class="btn" type="button">保存</button>
      </div>
    </div>
  </div>
  <div class="tier-switch" id="tier-switch">
    <button class="tier-btn active" type="button" data-tier="beginner">🔰 かんたん表示</button>
    <button class="tier-btn" type="button" data-tier="advanced">🔧 くわしい表示</button>
    <button class="tier-btn" type="button" data-tier="duct" title="DUCT予報タブへ切り替えます(タブ切替ボタンは消えません)">🌫️ ダクト予報</button>
  </div>
  <div id="root">
    <div class="card startup" id="startup">
      <div class="sub">起動中… NICT観測データを取得しています</div>
      <div class="pct-big" id="startup-pct">5%</div>
      <div class="progress"><span id="startup-bar" style="width:5%"></span></div>
      <div class="sub" id="startup-msg">初期化中</div>
    </div>
  </div>
  <!-- ITPFS Ver1.0(修正): DUCTタブ表示時に「学習データ統合」カードが地図より
       上に出てしまう指摘のため、DOM順序をduct→mergeへ入れ替える。
       通常タブ(かんたん/くわしい)では#duct-ui-persistentは非表示(高さ0)に
       なるため、見た目上の影響は無い。 -->
  <div id="duct-ui-persistent"></div>
  <!-- Rev16.2(不具合修正): 他端末学習統合UIは#root(20秒毎に丸ごと再構築される)の
       "外側"に永続配置する。理由はJS側のinitMergeUI()コメント参照。 -->
  <div id="merge-ui-persistent"></div>
  <div class="footer" id="footer"></div>
  <div class="copyright-footer">Copyright (c) 2026 JL7KHN/栃技研
[利用・免責条件] 個人的利用における本ソースコードの改変・利用は自由です。本プログラムの結果および使用に伴ういかなる不利益・損害についても、作成者（JL7KHN/栃技研）は一切の責任を負いません。</div>
</div>
"""
_v13_HTML_HEAD = _v13_HTML_HEAD.replace("__EDFS_TITLE_VERSION__", VERSION_LABEL)
_v13_HTML_HEAD = _v13_HTML_HEAD.replace("__DUCT_EMBED_CSS__", _DUCT_EMBED_CSS)

_v13_PREFECTURES_JS = r"""
// Ver13.11: 精緻化ポリゴン v2 (389頂点、14陸地)
// 島根半島の突出と九州西岸の折り返しを修正、頂点順序を全て反時計回りで統一
const JAPAN_COAST = {
  hokkaido: [[141.940,45.520],[141.700,45.400],[141.480,45.250],[141.400,44.950],[141.200,44.700],[141.020,44.400],[140.850,44.150],[140.550,43.900],[140.400,43.700],[140.250,43.500],[140.200,43.350],[140.000,43.300],[140.100,43.100],[140.300,42.950],[140.450,42.800],[140.500,42.550],[140.400,42.250],[140.350,42.000],[140.200,41.750],[140.100,41.550],[140.070,41.420],[140.250,41.400],[140.500,41.550],[140.800,41.700],[141.000,41.850],[141.200,42.000],[141.550,42.150],[141.850,42.300],[142.150,42.250],[142.500,42.150],[142.900,42.150],[143.240,41.930],[143.500,42.150],[143.850,42.400],[144.150,42.700],[144.500,42.900],[144.850,43.000],[145.100,43.100],[145.400,43.200],[145.820,43.390],[145.500,43.550],[145.250,43.700],[145.100,43.850],[145.200,44.100],[145.350,44.300],[145.150,44.400],[144.850,44.400],[144.500,44.250],[144.200,44.150],[143.900,44.100],[143.550,44.100],[143.100,44.250],[142.700,44.500],[142.400,44.750],[142.150,45.050],[141.950,45.250],[141.850,45.420],[141.940,45.520]],
  honshu: [[140.910,41.530],[140.750,41.400],[140.550,41.300],[140.400,41.150],[140.200,41.000],[140.100,40.800],[139.980,40.500],[139.910,40.200],[139.900,39.900],[139.850,39.600],[139.700,39.900],[139.750,39.700],[139.850,39.500],[139.900,39.200],[139.900,38.900],[139.780,38.650],[139.700,38.450],[139.550,38.300],[139.400,38.100],[139.200,37.900],[138.950,37.750],[138.750,37.550],[138.500,37.350],[138.250,37.200],[138.000,37.150],[137.800,37.100],[137.600,37.100],[137.400,37.150],[137.250,37.300],[137.360,37.550],[137.150,37.500],[136.900,37.300],[136.750,37.100],[136.700,36.900],[136.750,36.750],[136.550,36.550],[136.250,36.350],[136.200,36.150],[136.200,35.950],[136.150,35.750],[135.950,35.700],[135.750,35.700],[135.500,35.650],[135.300,35.700],[135.100,35.700],[134.900,35.700],[134.700,35.650],[134.500,35.650],[134.250,35.600],[134.000,35.600],[133.700,35.600],[133.400,35.600],[133.100,35.550],[132.900,35.500],[132.750,35.450],[132.620,35.380],[132.500,35.250],[132.300,35.100],[132.100,34.950],[131.900,34.800],[131.700,34.700],[131.500,34.600],[131.300,34.500],[131.100,34.400],[130.950,34.300],[130.900,34.150],[130.940,33.940],[131.150,33.950],[131.400,33.930],[131.700,33.930],[131.950,33.950],[132.150,34.050],[132.400,34.150],[132.700,34.250],[132.900,34.300],[133.100,34.350],[133.400,34.400],[133.700,34.450],[134.000,34.450],[134.250,34.500],[134.500,34.550],[134.700,34.600],[134.900,34.550],[135.050,34.550],[135.150,34.600],[135.250,34.650],[135.400,34.650],[135.400,34.500],[135.350,34.350],[135.250,34.200],[135.150,34.100],[135.100,33.950],[135.200,33.800],[135.350,33.700],[135.500,33.600],[135.760,33.440],[135.950,33.550],[136.150,33.700],[136.300,33.850],[136.500,34.000],[136.700,34.150],[136.850,34.300],[136.850,34.500],[136.850,34.650],[137.100,34.700],[137.500,34.700],[137.850,34.650],[138.100,34.650],[138.300,34.650],[138.400,34.680],[138.550,34.750],[138.750,34.800],[138.850,34.720],[138.950,34.900],[139.100,35.100],[139.300,35.250],[139.500,35.300],[139.650,35.300],[139.750,35.150],[139.850,34.900],[140.100,35.050],[140.350,35.250],[140.550,35.450],[140.700,35.700],[140.870,35.710],[140.800,35.900],[140.750,36.100],[140.800,36.400],[140.850,36.700],[140.950,37.000],[141.050,37.150],[141.020,37.400],[141.000,37.750],[141.050,38.000],[141.150,38.300],[141.400,38.400],[141.700,38.500],[141.900,38.700],[141.950,38.900],[141.950,39.150],[142.000,39.400],[142.050,39.700],[141.950,40.000],[141.850,40.300],[141.700,40.550],[141.550,40.850],[141.450,41.150],[141.430,41.400],[141.300,41.500],[141.150,41.520],[140.910,41.530]],
  shikoku: [[134.650,34.240],[134.700,34.100],[134.750,33.950],[134.700,33.750],[134.600,33.550],[134.400,33.400],[134.180,33.250],[134.000,33.350],[133.800,33.400],[133.550,33.400],[133.300,33.300],[133.150,33.100],[133.020,32.720],[132.850,32.850],[132.700,33.000],[132.550,33.220],[132.500,33.400],[132.500,33.550],[132.550,33.700],[132.650,33.850],[132.700,34.000],[132.750,34.100],[132.900,34.150],[133.100,34.200],[133.300,34.250],[133.550,34.300],[133.850,34.350],[134.100,34.400],[134.300,34.400],[134.500,34.400],[134.650,34.300],[134.650,34.240]],
  kyushu: [[130.950,33.950],[131.150,33.850],[131.350,33.800],[131.550,33.750],[131.650,33.650],[131.850,33.550],[132.000,33.350],[132.050,33.250],[131.950,33.100],[131.850,32.900],[131.800,32.700],[131.700,32.500],[131.650,32.350],[131.550,32.100],[131.500,31.900],[131.430,31.650],[131.340,31.370],[131.150,31.200],[130.950,31.050],[130.800,31.000],[130.660,30.990],[130.500,31.100],[130.450,31.250],[130.550,31.400],[130.550,31.600],[130.400,31.400],[130.300,31.250],[130.200,31.500],[130.150,31.750],[130.150,32.000],[130.100,32.250],[130.050,32.450],[130.100,32.600],[130.150,32.700],[130.100,32.850],[129.900,32.800],[129.750,32.650],[129.800,32.850],[129.700,33.000],[129.600,33.150],[129.700,33.300],[129.850,33.400],[129.950,33.500],[130.050,33.500],[130.250,33.550],[130.350,33.700],[130.500,33.650],[130.650,33.700],[130.850,33.850],[130.950,33.950]],
  okinawa: [[128.250,26.870],[128.200,26.780],[128.100,26.700],[128.050,26.600],[127.950,26.550],[127.900,26.500],[127.850,26.400],[127.850,26.300],[127.800,26.220],[127.780,26.150],[127.720,26.100],[127.710,26.080],[127.680,26.130],[127.660,26.180],[127.680,26.250],[127.760,26.320],[127.870,26.400],[127.980,26.470],[128.080,26.550],[128.170,26.650],[128.250,26.750],[128.290,26.830],[128.250,26.870]],
  sado: [[138.180,38.330],[138.400,38.330],[138.550,38.250],[138.550,38.050],[138.500,37.900],[138.350,37.780],[138.200,37.800],[138.080,37.950],[138.080,38.150],[138.180,38.330]],
  awaji: [[135.050,34.750],[135.100,34.650],[135.100,34.500],[135.050,34.300],[134.980,34.220],[134.850,34.250],[134.780,34.400],[134.750,34.550],[134.850,34.720],[135.050,34.750]],
  tsushima: [[129.200,34.700],[129.400,34.650],[129.480,34.500],[129.480,34.350],[129.400,34.200],[129.300,34.100],[129.200,34.150],[129.180,34.300],[129.180,34.500],[129.200,34.700]],
  tanega: [[130.980,30.830],[131.050,30.700],[131.020,30.550],[130.950,30.450],[130.850,30.500],[130.820,30.700],[130.980,30.830]],
  yaku: [[130.500,30.450],[130.680,30.450],[130.700,30.320],[130.600,30.200],[130.450,30.230],[130.380,30.350],[130.500,30.450]],
  amami: [[129.550,28.550],[129.700,28.500],[129.750,28.350],[129.650,28.180],[129.450,28.150],[129.300,28.250],[129.320,28.400],[129.450,28.520],[129.550,28.550]],
  miyako: [[125.150,24.850],[125.450,24.850],[125.500,24.760],[125.400,24.680],[125.180,24.680],[125.100,24.760],[125.150,24.850]],
  ishigaki: [[124.100,24.500],[124.350,24.500],[124.400,24.380],[124.320,24.280],[124.150,24.290],[124.050,24.400],[124.100,24.500]],
  jeju: [[126.150,33.550],[126.400,33.550],[126.600,33.400],[126.550,33.250],[126.250,33.200],[126.100,33.300],[126.100,33.450],[126.150,33.550]],
};

const PREFECTURES = [
  // 北海道 (4分割: 道央/道南/道東/道北)
  {name:'北海道(道南)', pts:[[140.05,41.42],[140.28,41.28],[140.55,41.37],[140.85,41.44],[141.15,41.62],[141.28,42.02],[141.05,42.20],[140.75,42.32],[140.42,42.15],[140.18,41.85],[140.05,41.55]]},
  {name:'北海道(道央)', pts:[[140.55,42.20],[141.35,42.52],[142.10,42.85],[142.85,43.28],[142.95,43.85],[142.55,44.02],[141.85,43.85],[141.25,43.42],[140.85,42.95],[140.65,42.55]]},
  {name:'北海道(道北)', pts:[[141.85,43.85],[142.55,44.02],[142.95,44.42],[142.70,44.85],[142.15,45.28],[141.55,45.42],[141.20,44.95],[141.35,44.35],[141.68,43.95]]},
  {name:'北海道(道東)', pts:[[142.95,43.28],[143.65,42.95],[144.35,42.98],[145.28,43.30],[145.85,44.28],[145.42,44.55],[144.55,44.45],[143.85,44.15],[143.35,43.75],[142.95,43.42]]},

  // 東北 6県
  {name:'青森', pts:[[140.35,40.55],[140.85,40.42],[141.25,40.38],[141.55,40.62],[141.68,41.05],[141.42,41.35],[141.05,41.42],[140.55,41.35],[140.25,41.00],[140.28,40.72]]},
  {name:'秋田', pts:[[139.85,39.45],[140.28,39.35],[140.55,39.55],[140.55,40.15],[140.42,40.55],[140.15,40.62],[139.90,40.35],[139.75,39.85]]},
  {name:'岩手', pts:[[141.05,39.05],[141.55,38.95],[141.95,39.15],[142.05,39.85],[141.85,40.35],[141.55,40.42],[141.25,40.05],[141.05,39.55]]},
  {name:'山形', pts:[[139.55,38.35],[140.05,38.28],[140.45,38.55],[140.55,39.15],[140.28,39.45],[139.85,39.42],[139.55,39.05],[139.42,38.65]]},
  {name:'宮城', pts:[[140.45,37.85],[141.05,37.85],[141.55,38.05],[141.75,38.55],[141.55,38.95],[141.15,38.98],[140.75,38.75],[140.42,38.35]]},
  {name:'福島', pts:[[139.45,36.85],[140.35,36.85],[140.95,36.95],[141.05,37.55],[140.95,37.85],[140.45,37.85],[139.85,37.75],[139.42,37.35]]},

  // 関東 7都県
  {name:'茨城', pts:[[140.02,35.85],[140.42,35.75],[140.75,35.98],[140.85,36.35],[140.85,36.85],[140.55,36.95],[140.18,36.85],[139.95,36.42]]},
  {name:'栃木', pts:[[139.35,36.25],[139.75,36.20],[140.05,36.35],[140.05,36.85],[139.85,37.05],[139.55,37.02],[139.32,36.75],[139.28,36.42]]},
  {name:'群馬', pts:[[138.45,36.05],[138.85,35.98],[139.28,36.15],[139.45,36.45],[139.42,36.85],[139.15,37.05],[138.75,37.00],[138.42,36.65]]},
  {name:'埼玉', pts:[[138.75,35.75],[139.15,35.72],[139.55,35.78],[139.85,35.85],[139.85,36.18],[139.42,36.22],[139.02,36.15],[138.72,35.98]]},
  {name:'千葉', pts:[[139.75,34.95],[140.05,34.92],[140.35,35.05],[140.85,35.35],[140.85,35.75],[140.42,35.85],[140.05,35.72],[139.85,35.42],[139.75,35.15]]},
  {name:'東京', pts:[[138.95,35.55],[139.25,35.55],[139.55,35.65],[139.85,35.68],[139.85,35.82],[139.55,35.85],[139.15,35.78],[138.92,35.72]]},
  {name:'神奈川', pts:[[139.05,35.15],[139.35,35.15],[139.65,35.25],[139.75,35.42],[139.72,35.62],[139.45,35.65],[139.15,35.55],[139.02,35.35]]},

  // 中部 (北陸+甲信越+東海)
  {name:'新潟', pts:[[137.85,36.75],[138.55,36.85],[139.15,37.15],[139.55,37.55],[139.75,38.05],[139.55,38.45],[139.15,38.55],[138.65,38.35],[138.25,38.02],[137.95,37.45]]},
  {name:'長野', pts:[[137.55,35.25],[138.05,35.22],[138.55,35.35],[138.75,35.85],[138.68,36.45],[138.42,36.95],[137.95,37.00],[137.55,36.65],[137.42,36.05],[137.48,35.55]]},
  {name:'山梨', pts:[[138.25,35.25],[138.55,35.22],[138.85,35.35],[138.95,35.65],[138.85,35.85],[138.55,35.85],[138.28,35.72]]},
  {name:'静岡', pts:[[137.55,34.65],[138.05,34.62],[138.55,34.72],[138.95,34.95],[138.95,35.25],[138.55,35.28],[138.05,35.15],[137.65,34.95]]},
  {name:'愛知', pts:[[136.75,34.65],[137.05,34.62],[137.35,34.72],[137.65,34.85],[137.65,35.20],[137.35,35.35],[137.05,35.28],[136.75,35.10]]},
  {name:'岐阜', pts:[[136.35,35.15],[136.75,35.15],[137.15,35.35],[137.55,35.65],[137.65,36.15],[137.35,36.45],[136.85,36.42],[136.42,36.15],[136.35,35.65]]},
  {name:'富山', pts:[[136.75,36.45],[137.15,36.42],[137.55,36.55],[137.75,36.85],[137.55,36.92],[137.15,36.85],[136.85,36.75]]},
  {name:'石川', pts:[[136.35,36.15],[136.65,36.20],[136.85,36.55],[137.15,37.05],[137.25,37.55],[136.85,37.42],[136.42,36.85],[136.28,36.42]]},
  {name:'福井', pts:[[135.45,35.45],[135.85,35.42],[136.25,35.55],[136.55,35.85],[136.65,36.25],[136.25,36.28],[135.85,36.05],[135.55,35.75]]},

  // 近畿 2府5県
  {name:'滋賀', pts:[[135.85,34.85],[136.15,34.85],[136.42,35.05],[136.42,35.55],[136.15,35.68],[135.85,35.65],[135.75,35.25]]},
  {name:'三重', pts:[[136.05,33.75],[136.35,33.72],[136.65,33.95],[136.85,34.55],[136.85,35.15],[136.55,35.28],[136.15,35.05],[135.95,34.55]]},
  {name:'京都', pts:[[135.05,34.85],[135.45,34.85],[135.75,34.98],[135.85,35.35],[135.75,35.75],[135.42,35.75],[135.05,35.55],[134.95,35.15]]},
  {name:'大阪', pts:[[135.05,34.35],[135.35,34.32],[135.65,34.45],[135.72,34.75],[135.65,34.92],[135.35,34.92],[135.05,34.75]]},
  {name:'兵庫', pts:[[134.25,34.25],[134.75,34.15],[135.15,34.35],[135.35,34.75],[135.42,35.25],[135.15,35.65],[134.65,35.68],[134.28,35.42],[134.15,34.85]]},
  {name:'奈良', pts:[[135.65,33.95],[135.85,33.92],[136.15,34.05],[136.25,34.55],[136.05,34.85],[135.75,34.85],[135.55,34.55]]},
  {name:'和歌山', pts:[[135.05,33.55],[135.45,33.45],[135.85,33.55],[136.05,33.95],[135.85,34.35],[135.42,34.35],[135.05,34.05]]},

  // 中国 5県
  {name:'鳥取', pts:[[133.35,35.15],[133.75,35.15],[134.15,35.28],[134.55,35.45],[134.55,35.65],[134.15,35.65],[133.75,35.55],[133.35,35.42]]},
  {name:'島根', pts:[[131.85,34.55],[132.25,34.55],[132.75,34.75],[133.15,35.05],[133.35,35.42],[133.35,35.55],[132.95,35.55],[132.35,35.35],[131.85,35.05]]},
  {name:'岡山', pts:[[133.35,34.45],[133.75,34.42],[134.15,34.55],[134.45,34.85],[134.45,35.15],[134.05,35.25],[133.65,35.15],[133.35,34.85]]},
  {name:'広島', pts:[[132.15,34.05],[132.55,33.98],[133.05,34.15],[133.55,34.42],[133.55,34.85],[133.15,34.95],[132.65,34.85],[132.25,34.55]]},
  {name:'山口', pts:[[130.75,33.75],[131.15,33.72],[131.55,33.85],[131.95,34.15],[132.15,34.55],[131.85,34.75],[131.35,34.65],[130.95,34.35],[130.75,34.05]]},

  // 四国 4県
  {name:'徳島', pts:[[133.65,33.65],[134.05,33.55],[134.45,33.65],[134.75,33.95],[134.75,34.15],[134.45,34.25],[134.05,34.15],[133.75,33.95]]},
  {name:'香川', pts:[[133.55,34.10],[133.95,34.05],[134.25,34.15],[134.55,34.25],[134.55,34.45],[134.25,34.45],[133.95,34.35],[133.65,34.28]]},
  {name:'愛媛', pts:[[132.35,33.02],[132.85,32.95],[133.35,33.15],[133.75,33.45],[133.75,33.95],[133.35,34.05],[132.85,33.95],[132.42,33.65],[132.28,33.35]]},
  {name:'高知', pts:[[132.55,32.75],[133.05,32.72],[133.55,32.85],[134.05,33.15],[134.25,33.55],[133.85,33.65],[133.35,33.55],[132.85,33.35],[132.55,33.05]]},

  // 九州 7県 + 沖縄
  {name:'福岡', pts:[[130.05,33.15],[130.45,33.05],[130.85,33.15],[131.15,33.45],[131.15,33.85],[130.85,33.98],[130.45,33.85],[130.15,33.65],[130.02,33.35]]},
  {name:'大分', pts:[[130.95,32.85],[131.35,32.75],[131.75,32.95],[132.05,33.35],[132.05,33.75],[131.75,33.85],[131.35,33.75],[131.05,33.45]]},
  {name:'佐賀', pts:[[129.75,33.05],[130.05,32.95],[130.28,33.15],[130.28,33.55],[130.05,33.65],[129.85,33.55],[129.72,33.28]]},
  {name:'長崎', pts:[[129.45,32.65],[129.75,32.55],[130.05,32.75],[130.15,33.15],[130.15,34.05],[129.75,34.15],[129.45,33.85],[129.35,33.35],[129.42,32.95]]},
  {name:'熊本', pts:[[130.05,32.25],[130.45,32.15],[130.85,32.35],[131.15,32.75],[131.15,33.15],[130.85,33.25],[130.42,33.15],[130.15,32.85],[130.02,32.55]]},
  {name:'宮崎', pts:[[131.05,31.45],[131.35,31.35],[131.75,31.55],[131.95,32.15],[131.85,32.75],[131.55,32.85],[131.15,32.65],[131.02,32.15],[131.05,31.75]]},
  {name:'鹿児島', pts:[[130.15,30.45],[130.55,30.35],[130.95,30.55],[131.15,31.15],[131.15,31.85],[131.02,32.15],[130.65,32.15],[130.25,31.85],[130.05,31.35],[130.05,30.85]]},
  {name:'沖縄本島', pts:[[127.65,26.15],[127.85,26.10],[128.05,26.28],[128.25,26.55],[128.35,26.85],[128.15,27.02],[127.85,26.98],[127.68,26.75],[127.62,26.42]]},
  {name:'宮古/八重山', pts:[[124.10,24.30],[124.35,24.25],[124.55,24.40],[124.65,24.55],[124.55,24.65],[124.25,24.60],[124.10,24.45]]},

  // 佐渡島
  {name:'佐渡', pts:[[138.15,37.85],[138.35,37.75],[138.55,37.85],[138.55,38.15],[138.35,38.35],[138.15,38.25],[138.05,38.05]]},
];
"""

_DUCT_EMBED_BODY_HTML = r"""<div id="hd"><h1>🗾 DUCT予測 詳細地図</h1><div class="meta">
<span class="badge" id="genBadge">生成 --:--:--</span>
<span class="badge" id="nptsBadge">地点 --</span>
<span class="badge eta" id="etaBadge">実用まで --</span>
<span class="badge" id="cd">次回更新 --:--</span>
</div></div>
<div id="bar">
<button id="play">▶ 再生</button>
<label>予報時刻 <span id="tlabel"></span></label>
<input type="range" id="tslider" min="0" max="0" value="0" step="1">
<button id="pauseBtn" title="観測・DPP計算エンジンの一時停止/再開（時刻アニメの再生ボタンとは別機能）">⏯ 観測稼働中</button>
<button id="contourBtn" onclick="duct2_toggleContour()" title="ダクト形成可能性が高い経路をコンター線(等値線)で表示">📈 コンター</button>
<button onclick="duct2_openSettings()" title="設定">⚙ 設定</button>
<button onclick="duct2_askPass()" title="自局と相手局、2地点間のPass(伝搬経路)を診断">🛤️ Pass</button>
<button onclick="duct2_askAntenna()" title="自局1地点のみ。アンテナ高を上げた場合の効果を試算">📶 アンテナ高</button>
<button onclick="duct2_askRadar()" title="自局1地点のみ。そこから全方位への到達距離を表示">🧭 移動運用</button>
<span class="badge" id="readout" title="地図をタップ(またはマウスで指した地点)の推定DPP値を表示します">📍 地点DPP: タップして確認</span>
<div class="legend"><span>0</span><div class="grad"></div><span>100 (DPP指数)</span></div>
<span class="badge" id="contourLegend" style="display:none">等値線: 30 / 50 / 70(実線太) / 85(実線最太)</span>
</div>
<div id="ovl" onclick="if(event.target===this)duct2_closeSettings()">
 <div id="panel">
  <div class="phd"><b>⚙ 設定</b><span class="x" onclick="duct2_closeSettings()">✕</span></div>
  <p class="note">変更は<b>保存で即時反映</b>(サーバへ保存し即再計算)。</p>
  <label>更新周期 <select id="s_interval"></select></label>
  <label>常駐時間 <select id="s_total" title="0=無制限。設定時間が経過するとサーバーごと自動停止します"></select></label>
  <p class="note" style="margin-top:-6px">⚠ 0以外を選ぶと、その時間経過でサーバーが自動的に停止します(0=無制限、ずっと動かす場合はこちら)</p>
  <label>気象キャッシュ <select id="s_wx"></select></label>
  <label>スポット遡り <select id="s_back"></select></label>
  <p class="note">🔔 アラート(Discord Webhook通知)</p>
  <label>有効化 <input type="checkbox" id="s_alert_on"></label>
  <label>Webhook URL <input type="text" id="s_webhook" placeholder="https://discord.com/api/webhooks/..." style="width:200px"></label>
  <label>DPP閾値 <select id="s_alert_dpp"></select></label>
  <label>監視半径(km) <select id="s_alert_radius"></select></label>
  <div class="prow">
   <button onclick="duct2_applySettings()">保存して即時更新</button>
   <button onclick="duct2_recomputeNow()">今すぐ再計算</button>
   <button onclick="duct2_closeSettings()">閉じる</button>
  </div>
  <div class="prow"><span id="s_msg" class="smsg"></span></div>
  <div id="s_status" class="sstat"></div>
 </div>
</div>
<div id="rpanel">
 <div class="phd"><b>📡 地点診断レポート</b><span class="x" onclick="duct2_closeReport()">✕</span></div>
 <div id="rpbody"></div>
</div>
<div id="covl" onclick="if(event.target===this)duct2_closeCoordModal()">
 <div id="cpanel">
  <div class="phd"><b id="cp_title">座標入力</b><span class="x" onclick="duct2_closeCoordModal()">✕</span></div>
  <div id="cp_body"></div>
  <div class="prow">
   <button onclick="duct2_cp_submit()">実行</button>
   <button onclick="duct2_closeCoordModal()">閉じる</button>
  </div>
  <span id="cp_msg" class="smsg"></span>
 </div>
</div>
<div id="hpanel"></div>
<div id="map">
 <div id="duct2-placeholder">
  <div class="duct2-ph-msg">📡 DUCT予測データを取得中…</div>
 </div>
</div>
"""
_DUCT_EMBED_SCRIPT_JS = r"""(function(){
'use strict';
// ITPFS(整理版・大規模UI改修): 単体版DUCT(V6.2.1)を「窓の中(iframe/別タブ)」で
// 表示する方式を全面的に廃止し、ダッシュボード本体のDOMへ直接統合した。
// 従来はgenerate_html()がPython文字列置換で__CELLDATA__等をHTMLへ直接焼き込んで
// いたが、統合版では/api/duct/celldataから取得したJSONで以下の変数を都度
// 更新する(初期値は空。duct2FetchAndRender()参照)。
var CELLS=[], FRAMES_T=[], STEP=2.0, REGION=null;
var ANTENNA_DEFAULT_H=10, BAND_LABELS=['50','144','430','1200'];
function viridis(t){t=Math.max(0,Math.min(1,t/100));
 var s=[[59,15,112],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
 var x=t*(s.length-1),i=Math.floor(x),f=x-i;if(i>=s.length-1){i=s.length-2;f=1;}
 var a=s[i],b=s[i+1];return [Math.round(a[0]+(b[0]-a[0])*f),
  Math.round(a[1]+(b[1]-a[1])*f),Math.round(a[2]+(b[2]-a[2])*f)];}
var map=null, fcv=null, _duct2MapInited=false;
var curFr=0;
function initLeafletMapOnce(){
 // ITPFS(整理版): REGIONが未取得(初回fetch前)の間は地図を初期化できないため、
 // 実データ取得後の1回目のみ呼び出される(duct2FetchAndRender参照)。
 if(_duct2MapInited || !REGION)return;
 _duct2MapInited=true;
 map=L.map('map',{zoomControl:true});
 /* Rev16.7(バグ修正): 従来使用していたCARTO Voyager(basemaps.cartocdn.com)は
    これ以外の場所は自由に使えていた無登録の無料タイルだが、CARTO側の運用方針
    変更により、未登録(APIキー無し)でのアクセス時にタイル画像内へ
    「Get your own API key」等の透かし文字が焼き込まれて配信されるようになった。
    対策として、登録・APIキーが一切不要で透かしも入らない標準の
    OpenStreetMap公式タイル(tile.openstreetmap.org)へ切り替える。 */
 L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
  {attribution:'&copy; OpenStreetMap contributors',subdomains:'abc',maxZoom:19}).addTo(map);
 map.fitBounds([[REGION.lat_min,REGION.lon_min],[REGION.lat_max,REGION.lon_max]]);
 fcv=L.DomUtil.create('canvas','',map.getPanes().overlayPane);
 fcv.id='fieldcv';fcv.style.position='absolute';fcv.style.pointerEvents='none';
 map.on('moveend zoomend resize',renderField);
 map.on('click',function(e){var best=-1,bd=1e9;
  for(var i=0;i<CELLS.length;i++){var dl=CELLS[i].lat-e.latlng.lat,
   dn=CELLS[i].lon-e.latlng.lng,d=dl*dl+dn*dn;if(d<bd){bd=d;best=i;}}
  if(best>=0&&bd<STEP*STEP)openPopup(best);
  fetchHepburnAt(e.latlng.lat,e.latlng.lng);
  updateReadout(idwAt(e.containerPoint.x,e.containerPoint.y));});
 map.on('mousemove',function(e){updateReadout(idwAt(e.containerPoint.x,e.containerPoint.y));});
 // ダッシュボード本体側からタブ切替時に呼ばれる(setUiTier()参照)。
 // #duct-ui-persistentが非表示(max-height:0)の間にLeafletが誤ったサイズで
 // 初期化される既知動作を防ぐため、可視化直後に必ずinvalidateSize()する。
 window._ductMapInvalidateSize=function(){
  try{ if(map) map.invalidateSize(); }catch(e){}
 };
 renderField();
}
function cellPixel(){var a=map.latLngToContainerPoint([REGION.lat_min,REGION.lon_min]);
 var b=map.latLngToContainerPoint([REGION.lat_min,REGION.lon_min+STEP]);
 return Math.max(10,Math.hypot(b.x-a.x,b.y-a.y));}
function renderField(){var size=map.getSize();fcv.width=size.x;fcv.height=size.y;
 L.DomUtil.setPosition(fcv,map.containerPointToLayerPoint([0,0]));
 var ctx=fcv.getContext('2d');ctx.clearRect(0,0,size.x,size.y);
 var pts=[];
 for(var i=0;i<CELLS.length;i++){var f=CELLS[i].frames[curFr];if(!f)continue;
  var p=map.latLngToContainerPoint([CELLS[i].lat,CELLS[i].lon]);
  pts.push({x:p.x,y:p.y,v:f.dpp,c:CELLS[i].conf});}
 if(!pts.length)return;
 var R=cellPixel()*2.6,R2=R*R,block=6;
 for(var y=0;y<size.y;y+=block){for(var x=0;x<size.x;x+=block){
  var sw=0,sv=0,sc=0,dmin=1e18;
  for(var k=0;k<pts.length;k++){var dx=x-pts[k].x,dy=y-pts[k].y,d2=dx*dx+dy*dy;
   if(d2>R2)continue;if(d2<dmin)dmin=d2;var w=1.0/(d2+25.0);
   sw+=w;sv+=w*pts[k].v;sc+=w*pts[k].c;}
  if(sw<=0)continue;var val=sv/sw,conf=sc/sw,dm=Math.sqrt(dmin);
  var mask=dm<R*0.30?1.0:Math.max(0,1.0-(dm-R*0.30)/(R*0.70));
  mask*=(0.45+0.55*conf);if(mask<=0.02)continue;
  var cc=viridis(val);
  ctx.fillStyle='rgba('+cc[0]+','+cc[1]+','+cc[2]+','+(0.74*mask).toFixed(3)+')';
  ctx.fillRect(x,y,block,block);}}
 drawContours(ctx,curFr);}
function idwAt(px,py){var sw=0,sv=0,R=cellPixel()*2.6,R2=R*R;
 for(var i=0;i<CELLS.length;i++){var f=CELLS[i].frames[curFr];if(!f)continue;
  var p=map.latLngToContainerPoint([CELLS[i].lat,CELLS[i].lon]);
  var dx=px-p.x,dy=py-p.y,d2=dx*dx+dy*dy;if(d2>R2)continue;
  var w=1.0/(d2+25.0);sw+=w;sv+=w*f.dpp;}
 return sw>0?sv/sw:null;}
/* ---- コンター(等値線)表示: 経度緯度グリッド上でIDW補間しMarching Squaresで等値線を抽出 ---- */
var CONTOUR_LEVELS=[30,50,70,85];
/* ITPFS Ver1.0(修正4・重要): 「統合前は初回からコンターが出ていたのに、統合後は
   一度も出ない」という報告の実際の原因はこれだった。CONTOUR_ONの既定値は
   localStorageの'dpp_contour'キーが文字列'1'の場合のみtrueになる設計だが、
   localStorageはオリジン(スキーム+ホスト+ポート)ごとに完全に分離される。
   統合前の単体版DUCTは別ポート(例:8787)で動いていたため、そちらで一度
   コンターをONにした形跡があってもITPFS側(別ポート、例:8765)には一切
   引き継がれず、常に「未設定」→既定値falseになっていた。コンターボタンの
   緑色っぽい枠線はCSSの既定スタイルであり、初回表示時にCONTOUR_ONの実値を
   反映していなかったため、見た目は「ON」に見えるのに実際は「OFF」という
   食い違いが起きていた。対策として、明示的に'0'が保存されている場合のみ
   OFFとし、それ以外(未設定=新規オリジン含む)は既定でONにする。 */
var CONTOUR_ON=(function(){try{var v=localStorage.getItem('dpp_contour');
 return v===null?true:(v==='1');}catch(e){return true;}})();
var contourCache={};
function idwAtLatLon(lat,lon,fr,R2deg){
 var sw=0,sv=0;
 for(var i=0;i<CELLS.length;i++){var f=CELLS[i].frames[fr];if(!f)continue;
  var dlat=lat-CELLS[i].lat,dlon=(lon-CELLS[i].lon)*Math.cos(lat*Math.PI/180);
  var d2=dlat*dlat+dlon*dlon;if(d2>R2deg)continue;
  var w=1.0/(d2+1e-7);sw+=w;sv+=w*f.dpp;}
 return sw>0?sv/sw:NaN;}
function buildGrid(fr,res){
 res=res||44;
 var latMin=REGION.lat_min,latMax=REGION.lat_max,lonMin=REGION.lon_min,lonMax=REGION.lon_max;
 var dlat=(latMax-latMin)/res,dlon=(lonMax-lonMin)/res;
 // ITPFS Ver1.0(修正): 初回サイクルは高速化のため粗いグリッド(STEPが大きい)で
 // 計算するため、従来通りR=STEP*2.6のままだと補間半径が地域全体を覆うほど
 // 巨大化し、DPP値がほぼ均一化されて等値線が1本も引けなくなる不具合があった。
 // 半径の上限を設け、粗い初回サイクルでも等値線が出るようにする。
 var R=Math.min(STEP*2.6, 6.0),R2=R*R;
 var vals=new Float32Array((res+1)*(res+1));
 for(var iy=0;iy<=res;iy++){var lat=latMin+iy*dlat;
  for(var ix=0;ix<=res;ix++){var lon=lonMin+ix*dlon;
   vals[iy*(res+1)+ix]=idwAtLatLon(lat,lon,fr,R2);}}
 return {res:res,latMin:latMin,lonMin:lonMin,dlat:dlat,dlon:dlon,vals:vals};}
function msInterp(v1,v2,p1,p2,level){
 if(!isFinite(v1)||!isFinite(v2))return p1;
 var t=(level-v1)/(v2-v1);if(!isFinite(t))t=0.5;t=Math.max(0,Math.min(1,t));
 return [p1[0]+(p2[0]-p1[0])*t,p1[1]+(p2[1]-p1[1])*t];}
function marchingSquares(grid,level){
 var res=grid.res,vals=grid.vals,segs=[];
 function vAt(ix,iy){return vals[iy*(res+1)+ix];}
 for(var iy=0;iy<res;iy++){for(var ix=0;ix<res;ix++){
  var tl=vAt(ix,iy),tr=vAt(ix+1,iy),br=vAt(ix+1,iy+1),bl=vAt(ix,iy+1);
  if(!isFinite(tl)||!isFinite(tr)||!isFinite(br)||!isFinite(bl))continue;
  var c=0;if(tl>level)c|=8;if(tr>level)c|=4;if(br>level)c|=2;if(bl>level)c|=1;
  if(c===0||c===15)continue;
  var latT=grid.latMin+iy*grid.dlat,lonL=grid.lonMin+ix*grid.dlon;
  var latB=latT+grid.dlat,lonR=lonL+grid.dlon;
  var pTL=[latT,lonL],pTR=[latT,lonR],pBR=[latB,lonR],pBL=[latB,lonL];
  var top=msInterp(tl,tr,pTL,pTR,level),right=msInterp(tr,br,pTR,pBR,level),
   bottom=msInterp(bl,br,pBL,pBR,level),left=msInterp(tl,bl,pTL,pBL,level);
  switch(c){
   case 1:segs.push([left,bottom]);break;
   case 2:segs.push([bottom,right]);break;
   case 3:segs.push([left,right]);break;
   case 4:segs.push([top,right]);break;
   case 5:segs.push([left,top]);segs.push([bottom,right]);break;
   case 6:segs.push([top,bottom]);break;
   case 7:segs.push([left,top]);break;
   case 8:segs.push([top,left]);break;
   case 9:segs.push([top,bottom]);break;
   case 10:segs.push([top,right]);segs.push([left,bottom]);break;
   case 11:segs.push([top,right]);break;
   case 12:segs.push([left,right]);break;
   case 13:segs.push([bottom,right]);break;
   case 14:segs.push([left,bottom]);break;}}}
 return segs;}
// ITPFS(整理版・コンター表示の最適化): 固定閾値[30,50,70,85]がその時点のDPP分布
// 範囲と一度も交差しない場合(較正が浅い/弱いダクト状況で最大値が20程度に留まる、
// あるいは学習が進んで全域が60超になる等)、固定閾値だけでは等値線が1本も
// 描画されず「コンターが全く表示されない」状態になり得る。この場合に限り、
// 実際の分布の40/65/85パーセンタイルから相対的な等値線レベルを自動算出し、
// 常に意味のある空間構造が見える状態を保証する(固定閾値が機能する通常時は
// 追加計算・描画そのものを行わないため、負荷・見た目とも影響しない)。
function computeAdaptiveLevels(grid){
 var vals=[];
 for(var i=0;i<grid.vals.length;i++){ if(isFinite(grid.vals[i])) vals.push(grid.vals[i]); }
 if(vals.length<8)return null;
 vals.sort(function(a,b){return a-b;});
 var gmin=vals[0], gmax=vals[vals.length-1];
 if(gmax-gmin<1.0)return null;  // 実質フラットな場では意味のある等値線を引けない
 var crossesFixed=false;
 for(var ci=0;ci<CONTOUR_LEVELS.length;ci++){
  if(CONTOUR_LEVELS[ci]>gmin && CONTOUR_LEVELS[ci]<gmax){crossesFixed=true;break;}}
 if(crossesFixed)return null;  // 固定閾値だけで十分表示できるため不要
 function pct(p){var idx=Math.min(vals.length-1,Math.max(0,Math.round(p*(vals.length-1))));return vals[idx];}
 var levels=[pct(0.4),pct(0.65),pct(0.85)];
 levels=levels.filter(function(v,i,a){return a.indexOf(v)===i;}).sort(function(a,b){return a-b;});
 return levels.length?levels:null;}
function getContours(fr){
 if(contourCache[fr])return contourCache[fr];
 var grid=buildGrid(fr),out={};
 for(var i=0;i<CONTOUR_LEVELS.length;i++){out[CONTOUR_LEVELS[i]]=marchingSquares(grid,CONTOUR_LEVELS[i]);}
 var adaptiveLevels=computeAdaptiveLevels(grid),adaptiveOut={};
 if(adaptiveLevels){
  for(var ai=0;ai<adaptiveLevels.length;ai++){
   adaptiveOut[adaptiveLevels[ai].toFixed(1)]=marchingSquares(grid,adaptiveLevels[ai]);}}
 var result={fixed:out,adaptive:adaptiveOut,adaptiveActive:!!adaptiveLevels};
 contourCache[fr]=result;return result;}
function updateContourLegend(adaptiveActive,adaptKeys){
 var lg=document.getElementById('contourLegend');if(!lg)return;
 if(adaptiveActive&&adaptKeys&&adaptKeys.length){
  lg.innerText='相対等値線(自動調整): '+adaptKeys.join(' / ')+' ※固定閾値の範囲外のため相対値で表示中';
  lg.title='DPP値の分布が固定閾値(30/50/70/85)と交差しなかったため、実際の分布から自動算出した目安の等値線を表示しています。';
 }else{
  lg.innerText='等値線: 30 / 50 / 70(実線太) / 85(実線最太)';
  lg.title='';
 }}
function drawContours(ctx,fr){
 if(!CONTOUR_ON)return;
 var cs=getContours(fr);
 for(var li=0;li<CONTOUR_LEVELS.length;li++){
  var lvl=CONTOUR_LEVELS[li],segs=cs.fixed[lvl];if(!segs||!segs.length)continue;
  var col=viridis(lvl);
  ctx.strokeStyle='rgba('+col[0]+','+col[1]+','+col[2]+',0.95)';
  ctx.lineWidth=lvl>=85?3.0:(lvl>=70?2.2:1.3);
  ctx.setLineDash(lvl>=70?[]:[6,4]);
  ctx.beginPath();
  for(var i=0;i<segs.length;i++){
   var a=map.latLngToContainerPoint([segs[i][0][0],segs[i][0][1]]);
   var b=map.latLngToContainerPoint([segs[i][1][0],segs[i][1][1]]);
   ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);}
  ctx.stroke();
  if(segs.length){ctx.setLineDash([]);
   var mid=segs[Math.floor(segs.length/2)];
   var mp=map.latLngToContainerPoint([(mid[0][0]+mid[1][0])/2,(mid[0][1]+mid[1][1])/2]);
   ctx.font='10px sans-serif';ctx.fillStyle='rgba('+col[0]+','+col[1]+','+col[2]+',0.95)';
   ctx.fillText(String(lvl),mp.x+3,mp.y-3);}}
 var adaptKeys=Object.keys(cs.adaptive||{});
 for(var aj=0;aj<adaptKeys.length;aj++){
  var segs2=cs.adaptive[adaptKeys[aj]];if(!segs2||!segs2.length)continue;
  ctx.strokeStyle='rgba(255,165,64,0.85)';ctx.lineWidth=1.4;ctx.setLineDash([3,3]);
  ctx.beginPath();
  for(var i2=0;i2<segs2.length;i2++){
   var a2=map.latLngToContainerPoint([segs2[i2][0][0],segs2[i2][0][1]]);
   var b2=map.latLngToContainerPoint([segs2[i2][1][0],segs2[i2][1][1]]);
   ctx.moveTo(a2.x,a2.y);ctx.lineTo(b2.x,b2.y);}
  ctx.stroke();
  if(segs2.length){ctx.setLineDash([]);
   var mid2=segs2[Math.floor(segs2.length/2)];
   var mp2=map.latLngToContainerPoint([(mid2[0][0]+mid2[1][0])/2,(mid2[0][1]+mid2[1][1])/2]);
   ctx.font='9px sans-serif';ctx.fillStyle='rgba(255,165,64,0.9)';
   ctx.fillText(adaptKeys[aj],mp2.x+3,mp2.y-3);}}
 ctx.setLineDash([]);
 updateContourLegend(cs.adaptiveActive,adaptKeys);}
function duct2_toggleContour(){CONTOUR_ON=!CONTOUR_ON;
 try{localStorage.setItem('dpp_contour',CONTOUR_ON?'1':'0');}catch(e){}
 var b=document.getElementById('contourBtn');if(b)b.style.borderColor=CONTOUR_ON?'#5ec962':'';
 var lg=document.getElementById('contourLegend');if(lg)lg.style.display=CONTOUR_ON?'':'none';
 renderField();}
function fetchHepburnAt(lat,lon){
 var el=document.getElementById('hpanel');
 if(!el)return;
 el.style.display='block';
 el.innerHTML='<div class="hp-row">Hepburn Map 照会中…</div>';
 fetch('/api/duct/hepburn?lat='+lat.toFixed(3)+'&lon='+lon.toFixed(3))
  .then(function(r){return r.json();})
  .then(function(d){
   if(!d.ok){
    var msg=d.error==='not_calibrated (geo)'?'未較正':
     d.error==='out_of_map_bounds'?'マップ範囲外':d.error;
    el.innerHTML='<div class="hp-row">Hepburn: '+msg+'</div>';return;}
   var rgb='rgb('+d.rgb.join(',')+')';
   el.innerHTML='<div class="hp-row">'
    +'<span class="hp-sw" style="background:'+rgb+'"></span>'
    +'<b>'+(d.label||'—')+'</b>'
    +'<span class="hp-sub"> ('+lat.toFixed(2)+', '+lon.toFixed(2)+') '
    +'valid '+d.valid_utc.replace('T',' ').slice(0,16)+'UTC</span></div>';
  })
  .catch(function(){el.innerHTML='<div class="hp-row">Hepburn: 通信エラー</div>';});}
function updateReadout(v){
 var el=document.getElementById('readout');
 if(el)el.innerText=(v==null?'📍 地点DPP: データ無し':'📍 地点DPP: '+v.toFixed(0));}
function draw(fr){curFr=fr;renderField();
 document.getElementById('tlabel').innerText=(FRAMES_T[fr]!=null?FRAMES_T[fr]+' JST':'--:--');
 try{localStorage.setItem('dpp_fr',fr);}catch(e){}
 if(curReportIdx>=0&&document.getElementById('rpanel').classList.contains('open')){
  openPopup(curReportIdx);}}
var curReportIdx=-1;
function openPopup(idx){curReportIdx=idx;var c=CELLS[idx],fr=+slider.value,f=c.frames[fr];
 if(!f)return;var bands=BAND_LABELS;
 var brow=bands.map(function(b,k){return '<tr><td>'+b+'MHz</td><td style=\"text-align:center\">'
  +(f.trap[k]==='1'?'○':'×')+'</td></tr>';}).join('');
 var ds=c.duct||'ダクトなし';
 var wkbHtml='';
 if(c.wkb){var w=c.wkb;
  wkbHtml='<div class=\"rp-sec\"><b>WKB('+w.freq_mhz+'MHz)</b><br>'
   +'捕捉周波数下限: <b>'+(w.f_min_mhz!=null?w.f_min_mhz.toFixed(0)+'MHz':'—')+'</b> / '
   +'捕捉モード: <b>'+w.n_modes+'</b><br>'
   +'結合効率: <b>'+w.coupling.toFixed(2)+'</b> / '
   +'推定漏洩: <b>'+(w.leak_db_per_100km!=null?w.leak_db_per_100km.toFixed(1)+' dB/100km':'—')+'</b></div>';}
 var corrHtml='';
 if(c.corr){var cr=c.corr;
  corrHtml='<div class=\"rp-sec\"><b>PSKReporter補正</b><br>'
   +'合成前 <b>'+cr.raw_dpp.toFixed(0)+'</b> → PSK補正 <b>×'+cr.factor.toFixed(2)+'</b> → '
   +'最終 <b>'+cr.corr_dpp.toFixed(0)+'</b><br>'
   +'観測数: <b>'+cr.n_obs+'</b> / 信頼度κ: <b>'+cr.kappa.toFixed(2)+'</b></div>';}
 var hepHtml='';
 if(c.corr&&c.corr.hepburn){var hp=c.corr.hepburn;
  if(hp.applied){
   hepHtml='<div class=\"rp-sec\"><b>Hepburn Map補正</b><br>'
    +'補正前 <b>'+hp.pre_dpp.toFixed(0)+'</b> → Hepburn期待値 <b>'+hp.expected_from_hepburn.toFixed(0)+'</b>'
    +'(κ='+hp.kappa.toFixed(2)+') → 最終 <b>'+hp.final_dpp.toFixed(0)+'</b><br>'
    +'回帰サンプル数: <b>'+hp.n+'</b> / R²: <b>'+hp.r2.toFixed(2)+'</b></div>';
  }else{
   hepHtml='<div class=\"rp-sec\"><b>Hepburn Map補正</b><br>'
    +'<span class=\"rp-sub\">サンプル蓄積中('+(hp.n||0)+'/'+(hp.min_required||'?')+') — 未適用</span></div>';
  }}
 var sumHtml=c.summary?('<div class=\"rp-summary\">'+c.summary+'</div>'):'';
 var body='<div class=\"rp-sub\">('+c.lat.toFixed(1)+', '+c.lon.toFixed(1)+')・予報時刻 '+FRAMES_T[fr]+'</div>'
  +'<div class=\"rp-dpp\">DPP '+f.dpp.toFixed(1)+'<span class=\"rp-sub\"> /100 (信頼度 '+c.conf.toFixed(2)+')</span></div>'
  +'<div class=\"rp-sub\">予測到達(430MHz): <b>'+f.r.toFixed(0)+' km</b> ・ 構造: '+ds+'</div>'
  +'<div class=\"rp-charts\">'+(c.svg||'')+(c.struct||'')+'</div>'
  +'<div class=\"rp-sec\"><b>DPP内訳</b>'+(c.comp||'')+'</div>'
  +'<table><tr><th>帯域</th><th>捕捉</th></tr>'+brow+'</table>'
  +'<div class=\"rp-sec\"><b>時系列(24h)</b>'+(c.ts||'')+'</div>'
  +wkbHtml+corrHtml+hepHtml+sumHtml;
 document.getElementById('rpbody').innerHTML=body;
 document.getElementById('rpanel').classList.add('open');}
function duct2_closeReport(){document.getElementById('rpanel').classList.remove('open');curReportIdx=-1;}
var slider=document.getElementById('tslider');
slider.oninput=function(){draw(+this.value);};
var playing=false,timer=null;
document.getElementById('play').onclick=function(){
 playing=!playing;this.innerText=playing?'⏸ 停止':'▶ 再生';
 if(playing){timer=setInterval(function(){
   if(!CELLS.length||!CELLS[0].frames)return;  // データ未到着時の安全策
   var v=(+slider.value+1)%(CELLS[0].frames.length);slider.value=v;draw(v);},900);}
 else{clearInterval(timer);}};
(function(){var b=document.getElementById('contourBtn');if(b)b.style.borderColor=CONTOUR_ON?'#5ec962':'';
 var lg=document.getElementById('contourLegend');if(lg)lg.style.display=CONTOUR_ON?'':'none';})();

// ITPFS(整理版・大規模UI改修): 従来の<meta http-equiv="refresh">による全ページ
// 再読み込み方式(単体版・iframe版共通)を廃止し、/api/duct/celldataを
// ソフトフェッチしてCELLS等を差し替える方式に統一した(ダッシュボード本体の
// 他の要素を巻き込んでページ全体を再読込しないため、体感速度・操作性が向上する)。
// 準備未完了(初回サイクル計算中)の間は5秒毎に軽量確認し、完了した瞬間に
// 自動的に地図を表示する(タブを開いているかどうかに関わらず反映される)。
var _duct2Ready=false, _duct2PollTimer=null, _duct2RefreshTimer=null;
var _duct2IntervalSec=300, _duct2Next=null;
function duct2ShowPlaceholder(msg){
 var ph=document.getElementById('duct2-placeholder');
 if(!ph)return; ph.style.display='flex';
 var m=ph.querySelector('.duct2-ph-msg'); if(m)m.innerText=msg;}
function duct2HidePlaceholder(){
 var ph=document.getElementById('duct2-placeholder');
 if(ph)ph.style.display='none';}
function duct2NotReadyMessage(d){
 if(d&&d.bootstrap_ok===false){
  return '⚠ DUCT予測エンジンの起動に失敗しました: '+(d.bootstrap_error||'不明な原因')
   +'。この状態は自動的には解消しません。アプリ/プログラムを再起動してください。';}
 if(d&&d.first_cycle_started&&!d.first_cycle_done){
  var cd=d.cycle_diag||{};
  return '📡 初回サイクルを計算中です(実行回数='+(cd.cycles_run||0)+'回)。'
   +'数分程度かかる場合があります。';}
 return '📡 DUCT予測データを取得中…(初回サイクル完了までお待ちください)';}
function duct2FetchAndRender(){
 return fetch('/api/duct/celldata?ts='+Date.now(),{cache:'no-store'})
  .then(function(r){return r.json();})
  .then(function(d){
   if(!d||!d.ok||!d.ready){
    _duct2Ready=false; duct2ShowPlaceholder(duct2NotReadyMessage(d));
    return false;}
   CELLS=d.cells; FRAMES_T=d.frame_labels; STEP=d.step_deg; REGION=d.region;
   ANTENNA_DEFAULT_H=d.tx_height_m; BAND_LABELS=d.band_list;
   _duct2IntervalSec=d.interval_sec||300;
   contourCache={};  // 新しいデータのため等値線キャッシュを破棄
   initLeafletMapOnce();
   duct2HidePlaceholder();
   var gb=document.getElementById('genBadge'); if(gb)gb.innerText='生成 '+d.generated_at;
   var nb=document.getElementById('nptsBadge'); if(nb)nb.innerText='地点 '+d.n_pts;
   var startFr=0;
   try{var sv=localStorage.getItem('dpp_fr');
    if(sv!==null)startFr=Math.max(0,Math.min(+sv,FRAMES_T.length-1));}catch(e){}
   if(!(+slider.value>=0 && +slider.value<FRAMES_T.length)) slider.value=startFr;
   slider.max=String(Math.max(0,FRAMES_T.length-1));
   draw(+slider.value);
   _duct2Ready=true;
   _duct2Next=new Date(new Date().getTime()+_duct2IntervalSec*1000);
   return true;})
  .catch(function(e){console.warn('[duct2] celldata取得失敗',e);return false;});}
function duct2StartSteadyRefresh(){
 if(_duct2RefreshTimer)clearInterval(_duct2RefreshTimer);
 _duct2RefreshTimer=setInterval(duct2FetchAndRender,Math.max(15000,_duct2IntervalSec*1000));}
function duct2StartPolling(){
 duct2FetchAndRender().then(function(ok){
  if(ok){duct2StartSteadyRefresh();return;}
  _duct2PollTimer=setInterval(function(){
   duct2FetchAndRender().then(function(ok2){
    if(ok2&&_duct2PollTimer){clearInterval(_duct2PollTimer);_duct2PollTimer=null;
     duct2StartSteadyRefresh();}});},5000);});}
function duct2Tick(){
 if(!_duct2Next)return;
 var s=Math.max(0,(_duct2Next-new Date())/1000);
 var m=String(Math.floor(s/60)).padStart(2,'0'),c=String(Math.floor(s%60)).padStart(2,'0');
 var el=document.getElementById('cd'); if(el)el.innerText='次回更新 '+m+':'+c;}
setInterval(duct2Tick,1000);
// タブ切替時にダッシュボード本体(setUiTier())から即時更新を要求された場合に
// 呼び出されるフック(iframe時代のwindow._ductMapForceReloadを踏襲)。
window._ductMapForceReload=function(){ duct2FetchAndRender(); };
duct2StartPolling();
// ===== 設定パネル =====
var CH_LABEL={interval_sec:function(v){return v+'秒 ('+(v/60).toFixed(v%60?1:0)+'分)';},
 total_minutes:function(v){return v===0?'無制限':v+'分';},
 weather_cache_sec:function(v){return (v/60)+'分';},
 seconds_back:function(v){return v+'秒 ('+(v/60)+'分)';}};
function fillSelect(id,key,choices,cur){var el=document.getElementById(id);
 el.innerHTML='';choices.forEach(function(v){var o=document.createElement('option');
  o.value=v;o.text=CH_LABEL[key]?CH_LABEL[key](v):v;
  if(+v===+cur)o.selected=true;el.appendChild(o);});}
function duct2_openSettings(){fetch('/api/duct/settings').then(function(r){return r.json();})
 .then(function(d){var s=d.settings,c=d.choices;
  fillSelect('s_interval','interval_sec',c.interval_sec,s.interval_sec);
  fillSelect('s_total','total_minutes',c.total_minutes,s.total_minutes);
  fillSelect('s_wx','weather_cache_sec',c.weather_cache_sec,s.weather_cache_sec);
  fillSelect('s_back','seconds_back',c.seconds_back,s.seconds_back);
  fillSelect('s_alert_dpp','alert_dpp_threshold',c.alert_dpp_threshold,s.alert_dpp_threshold);
  fillSelect('s_alert_radius','alert_radius_km',c.alert_radius_km,s.alert_radius_km);
  document.getElementById('s_alert_on').checked=!!s.alert_enabled;
  document.getElementById('s_webhook').value=s.discord_webhook_url||'';
  refreshStatus();document.getElementById('ovl').style.display='block';})
 .catch(function(e){alert('設定取得失敗: '+e);});}
function duct2_closeSettings(){document.getElementById('ovl').style.display='none';}
function fmtEta(h){if(h===null||h===undefined)return '推定中…';
 if(h<=0)return '到達済み';if(h<1)return '約'+Math.round(h*60)+'分';
 return '約'+h.toFixed(1)+'時間';}
function duct2_applySettings(){var body={
  interval_sec:+document.getElementById('s_interval').value,
  total_minutes:+document.getElementById('s_total').value,
  weather_cache_sec:+document.getElementById('s_wx').value,
  seconds_back:+document.getElementById('s_back').value,
  alert_enabled:document.getElementById('s_alert_on').checked,
  discord_webhook_url:document.getElementById('s_webhook').value,
  alert_dpp_threshold:+document.getElementById('s_alert_dpp').value,
  alert_radius_km:+document.getElementById('s_alert_radius').value};
 fetch('/api/duct/settings',{method:'POST',headers:{'Content-Type':'application/json'},
  body:JSON.stringify(body)}).then(function(r){return r.json();})
 .then(function(d){var m=document.getElementById('s_msg');
  m.innerText=d.recompute?'✓ 保存しました。即時再計算を開始':'✓ 保存しました(停止中は保留)';
  setTimeout(function(){m.innerText='';},6000);refreshStatus();})
 .catch(function(e){document.getElementById('s_msg').innerText='保存失敗: '+e;});}
var cpConfig=null;
function openCoordModal(cfg){
 cpConfig=cfg;
 document.getElementById('cp_title').innerText=cfg.title;
 var html='';
 cfg.fields.forEach(function(f){
  html+='<label>'+f.label
   +'<span class="cprow"><input type="text" id="'+f.id+'" placeholder="'+f.placeholder+'">'
   +'<button type="button" class="gpsbtn" onclick="duct2_fillGPS(\''+f.id+'\')" title="現在地(GPS)から自動取得">📍 現在地</button></span></label>';});
 if(cfg.heightField){var hf=cfg.heightField;
  html+='<label>'+hf.label
   +'<span class="cprow"><input type="text" id="'+hf.id+'" placeholder="'+hf.placeholder+'"></span></label>';}
 document.getElementById('cp_body').innerHTML=html;
 document.getElementById('cp_msg').innerText='';
 document.getElementById('covl').style.display='block';}
function duct2_closeCoordModal(){document.getElementById('covl').style.display='none';}
function duct2_fillGPS(inputId){
 var msg=document.getElementById('cp_msg');
 if(!navigator.geolocation){msg.innerText='⚠ このブラウザはGPS(位置情報)に対応していません';return;}
 msg.innerText='📍 現在地を取得中…';
 navigator.geolocation.getCurrentPosition(function(pos){
   var lat=pos.coords.latitude, lon=pos.coords.longitude;
   document.getElementById(inputId).value=lat.toFixed(5)+','+lon.toFixed(5);
   msg.innerText='✓ 現在地を取得しました(精度 約'+Math.round(pos.coords.accuracy)+'m)';},
  function(err){
   var reason=err.code===1?'位置情報の利用が許可されていません(ブラウザ/OS設定を確認)'
    :err.code===2?'現在地を取得できませんでした':'タイムアウトしました';
   msg.innerText='⚠ 取得失敗: '+reason;},
  {enableHighAccuracy:true,timeout:10000,maximumAge:60000});}
function looksLikeLatLon(v){
 var p=v.split(',').map(function(s){return s.trim();});
 if(p.length!==2)return false;
 var a=parseFloat(p[0]),b=parseFloat(p[1]);
 return isFinite(a)&&isFinite(b)&&Math.abs(a)<=90&&Math.abs(b)<=180;}
function geocodeResolve(q){
 return fetch('/api/duct/geocode?q='+encodeURIComponent(q))
  .then(function(r){return r.json();})
  .then(function(d){
   if(!d.ok)throw new Error('「'+q+'」が見つかりませんでした('+(d.error||'not_found')+')');
   return d.lat.toFixed(5)+','+d.lon.toFixed(5);});}
function duct2_cp_submit(){
 if(!cpConfig)return;
 var msg=document.getElementById('cp_msg');
 var raw={},ok=true;
 cpConfig.fields.forEach(function(f){var v=document.getElementById(f.id).value.trim();
  if(!v)ok=false;raw[f.id]=v;});
 if(!ok){msg.innerText='⚠ 座標または市町村名を入力するか📍現在地で取得してください';return;}
 if(cpConfig.heightField)raw[cpConfig.heightField.id]=document.getElementById(cpConfig.heightField.id).value.trim();
 msg.innerText='⌛ 座標を解決中…';
 var jobs=cpConfig.fields.map(function(f){
  var v=raw[f.id];
  if(looksLikeLatLon(v))return Promise.resolve();
  return geocodeResolve(v).then(function(ll){raw[f.id]=ll;});});
 Promise.all(jobs).then(function(){
  duct2_closeCoordModal();cpConfig.onSubmit(raw);
 }).catch(function(e){msg.innerText='⚠ '+e.message;});}
function apiErrMsg(code){
 if(code==='no_data'||code==='no_duct')return '計算データがまだありません(サーバー起動直後で初回計算中の可能性)。1〜2分待ってから再度お試しください。';
 if(code==='not_found')return '地名が見つかりませんでした。より具体的な市町村名を試すか、座標を直接入力してください。';
 return code;}
function duct2_askPass(){
 openCoordModal({title:'🛤️ Pass診断(自局⇔相手局)',
  fields:[{id:'cp_a',label:'自局の緯度,経度 または 市町村名',placeholder:'例: 35.68,139.77 または 東京都千代田区'},
   {id:'cp_b',label:'相手局の緯度,経度 または 市町村名',placeholder:'例: 33.59,130.40 または 福岡市'}],
  onSubmit:function(v){
   var pa=v.cp_a.split(',').map(parseFloat),pb=v.cp_b.split(',').map(parseFloat);
   fetch('/api/duct/pass?lat1='+pa[0]+'&lon1='+pa[1]+'&lat2='+pb[0]+'&lon2='+pb[1])
    .then(function(r){return r.json();}).then(function(d){
     if(!d.ok){alert('診断失敗: '+apiErrMsg(d.error));return;}
     alert('🛤️ Pass診断(自局⇔相手局)\n距離: '+d.distance_km+'km\n平均DPP: '+d.avg_dpp
      +'\n最小DPP: '+d.min_dpp+' / 最大DPP: '+d.max_dpp
      +'\n途切れリスク: '+d.break_risk_pct+'%');
    }).catch(function(e){alert('通信エラー: '+e);});}});}
function duct2_askAntenna(){
 openCoordModal({title:'📶 アンテナ高提案(自局単独)',
  fields:[{id:'cp_a',label:'自局の緯度,経度 または 市町村名',placeholder:'例: 38.27,140.87 または 仙台市'}],
  heightField:{id:'cp_h',label:'現在のアンテナ高[m](空欄で既定値'+ANTENNA_DEFAULT_H+'mを使用)',placeholder:String(ANTENNA_DEFAULT_H)},
  onSubmit:function(v){
   var p=v.cp_a.split(',').map(parseFloat);
   var h=(v.cp_h!=null && v.cp_h!=='') ? parseFloat(v.cp_h) : null;
   var url='/api/duct/antenna?lat='+p[0]+'&lon='+p[1]+(h!=null && !isNaN(h) ? '&height='+h : '');
   fetch(url)
    .then(function(r){return r.json();}).then(function(d){
     if(!d.ok){alert('取得失敗: '+apiErrMsg(d.error||'ダクト無し'));return;}
     var msg='📶 アンテナ高提案(自局単独)\nダクト底面: '+d.duct_base_m+'m (上端 '+d.duct_top_m+'m)\n'
      +'現在('+d.current_height_m+'m)の結合効率: '+(d.current_efficiency*100).toFixed(0)+'%';
     msg+=d.suggested_height_m!=null
      ?('\n効率70%以上に必要な高さ: 約'+d.suggested_height_m+'m'
        +(d.height_gap_m>0?(' (あと約'+d.height_gap_m+'m)'):' (既に到達)'))
      :'\nこのダクト構造では70%到達見込みなし';
     alert(msg);
    }).catch(function(e){alert('通信エラー: '+e);});}});}
var radarLayer=null;
function destPoint(lat,lon,bearingDeg,distKm){
 var R=6371.0088,brng=bearingDeg*Math.PI/180;
 var lat1=lat*Math.PI/180,lon1=lon*Math.PI/180;
 var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
 var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),
  Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
 return [lat2*180/Math.PI,((lon2*180/Math.PI)+540)%360-180];}
function showRadarOnMap(d){
 if(radarLayer){map.removeLayer(radarLayer);radarLayer=null;}
 var group=L.layerGroup();
 var starIcon=L.divIcon({className:'radar-star',html:'★',iconSize:[24,24],iconAnchor:[12,12]});
 L.marker([d.lat,d.lon],{icon:starIcon,zIndexOffset:1000})
  .bindTooltip('🧭 移動運用地点 ('+d.lat.toFixed(4)+', '+d.lon.toFixed(4)+')',{className:'radar-tip'})
  .addTo(group);
 var ranges=d.radar.map(function(x){return x.max_range_km;});
 var maxR=Math.max.apply(null,ranges.concat([1]));
 d.radar.forEach(function(x){
  if(x.max_range_km<=0)return;
  var end=destPoint(d.lat,d.lon,x.azimuth,x.max_range_km);
  var col=viridis(x.end_dpp!=null?x.end_dpp:(x.max_range_km/maxR)*100);
  var colStr='rgb('+col[0]+','+col[1]+','+col[2]+')';
  L.polyline([[d.lat,d.lon],end],{color:colStr,weight:1.4,opacity:0.85})
   .bindTooltip(x.azimuth+'° / '+x.max_range_km+'km'+(x.end_dpp!=null?' (DPP '+x.end_dpp+')':''),
    {className:'radar-tip',sticky:true})
   .addTo(group);
  L.circleMarker(end,{radius:2.5,color:colStr,fillColor:colStr,fillOpacity:1,weight:1})
   .addTo(group);});
 group.addTo(map);radarLayer=group;
 map.panTo([d.lat,d.lon]);}
function duct2_askRadar(){
 openCoordModal({title:'🧭 移動運用レーダー(自局から全方位)',
  fields:[{id:'cp_a',label:'自局(現在地)の緯度,経度 または 市町村名',placeholder:'例: 38.27,140.87 または 仙台市'}],
  onSubmit:function(v){
   var p=v.cp_a.split(',').map(parseFloat);
   fetch('/api/duct/radar?lat='+p[0]+'&lon='+p[1])
    .then(function(r){return r.json();}).then(function(d){
     if(!d.ok){alert('取得失敗: '+apiErrMsg(d.error||'ダクト無し'));return;}
     showRadarOnMap(d);
     var lines=d.radar.filter(function(x){return x.max_range_km>0;})
      .sort(function(a,b){return b.max_range_km-a.max_range_km;}).slice(0,8)
      .map(function(x){return x.azimuth+'° → '+x.max_range_km+'km';});
     var msg='🧭 移動運用レーダー(自局から全方位)\n'+(lines.join('\n')||'有効なダクト方向なし');
     if(d.elevation_m!=null)msg+='\n\n標高: '+d.elevation_m+'m';
     if(d.duct_base_m!=null)msg+=' / ダクト底面: '+d.duct_base_m+'m';
     if(d.height_to_duct_m!=null)msg+=' / 差: '+d.height_to_duct_m+'m';
     alert(msg);
    }).catch(function(e){alert('通信エラー: '+e);});}});}
function duct2_recomputeNow(){fetch('/api/duct/recompute',{method:'POST'})
 .then(function(r){return r.json();}).then(function(d){
  var m=document.getElementById('s_msg');m.innerText='⟳ 再計算を要求しました';
  setTimeout(function(){m.innerText='';},6000);})
 .catch(function(e){document.getElementById('s_msg').innerText='要求失敗: '+e;});}
function refreshStatus(){fetch('/api/duct/status').then(function(r){return r.json();})
 .then(function(d){var L=d.learning||{};var mat=L.maturity_pct||0;
  document.getElementById('s_status').innerHTML=
   '学習成熟度: <b>'+mat+'%</b>（目標 平均信頼度 '+(L.kappa_target||0.7)+'）'
   +'<div class="bar2"><i style="width:'+Math.min(100,mat)+'%"></i></div>'
   +'実用まで: <b>'+fmtEta(L.eta_hours)+'</b><br>'
   +'平均信頼度: '+(L.mean_confidence||0)+' / 被覆 '+(L.covered_cells||0)+'セル'
   +' / cone '+(d.cone_width_deg||'—')+'°<br>'
   +'状態: '+(d.settings.paused?'⏯ 観測停止中':'⏯ 観測稼働中');
  var eb=document.getElementById('etaBadge');
  if(eb)eb.innerText='実用まで '+fmtEta(L.eta_hours);}).catch(function(){});}
var pauseBtn=document.getElementById('pauseBtn');
function syncPause(){fetch('/api/duct/status').then(function(r){return r.json();})
 .then(function(d){pauseBtn.innerText=d.settings.paused?'⏯ 観測停止中':'⏯ 観測稼働中';});}
pauseBtn.onclick=function(){fetch('/api/duct/status').then(function(r){return r.json();})
 .then(function(d){return fetch('/api/duct/settings',{method:'POST',
   headers:{'Content-Type':'application/json'},
   body:JSON.stringify({paused:!d.settings.paused})});})
 .then(function(){syncPause();refreshStatus();});};
syncPause();
// ETAバッジをライブ更新(30秒毎)
setInterval(refreshStatus,30000);

// ITPFS(整理版): 静的HTML内のonclick=""属性はグローバルスコープで実行される
// ため、このIIFE内で定義した関数のうちHTMLから直接呼ばれるものだけを
// window直下へ明示的に公開する(他の内部変数・関数は一切漏らさない)。
window.duct2_toggleContour=duct2_toggleContour;
window.duct2_openSettings=duct2_openSettings;
window.duct2_closeSettings=duct2_closeSettings;
window.duct2_askPass=duct2_askPass;
window.duct2_askAntenna=duct2_askAntenna;
window.duct2_askRadar=duct2_askRadar;
window.duct2_applySettings=duct2_applySettings;
window.duct2_recomputeNow=duct2_recomputeNow;
window.duct2_closeReport=duct2_closeReport;
window.duct2_closeCoordModal=duct2_closeCoordModal;
window.duct2_cp_submit=duct2_cp_submit;
window.duct2_fillGPS=duct2_fillGPS;
})();
"""


_v13_JS_BODY = r"""
<script>
// ==================== EDFS Ver13.5 Frontend ====================
const REFRESH_MS = 20000;
// ITPFS(整理版・大規模UI改修): 単体版DUCT(V6.2.1相当)を「窓の中(iframe/別タブ)」
// で表示する方式を全面廃止し、ダッシュボード本体のDOMへ直接統合するための
// HTML断片/JavaScript本体。initDuctUI()がこの2つの定数からDOMを構築する
// (Python側でPLACEHOLDERをjson.dumps()した文字列へ実行時に置換する)。
const DUCT_EMBED_BODY_HTML = __DUCT_EMBED_BODY_HTML_JSON__;
const DUCT_EMBED_SCRIPT_JS = __DUCT_EMBED_SCRIPT_JS_JSON__;
let lastGpsSent = 0;
let latestSnapshot = null;
let currentHorizon = 0;  // 0/15/30/45/60/90

// Ver13.14: 予測地図の説明文 (現在=0分の場合の文言を分ける)
function mapSubText(h){
  const target = (h === 0) ? '緑破線円と同じ現在位置' : '+' + h + '分先のEs推定範囲';
  return '緑破線円: Es 現在推定位置 (時間同期で 1秒毎更新)。 赤破線: ' + target + '。 '
    + 'スライダで 現在/+15/30/45/60/90分を切替、または自動再生で移動方向を確認。 '
    + '三角: NICT 4観測局。 円: エリア別開通率予測 (半径=強度)。 ★: 自局。';
}

// Ver13.14: 自動再生タイマーの停止 (スライダ手動操作時にも呼ぶ)
function stopHorizonPlayback(){
  if(window._horizonPlayTimer){
    clearInterval(window._horizonPlayTimer);
    window._horizonPlayTimer = null;
    const btn = document.querySelector('.horizon-play-btn');
    if(btn) btn.textContent = '▶ 自動再生';
  }
}
// Ver13.5: Es 移動リアルタイムシミュレーション基準
let esSimBase = null;  // {epoch, lat, lon, speed_kmh, dir_deg}

function el(t, cls, txt){
  const e = document.createElement(t);
  if(cls) e.className = cls;
  if(txt !== undefined && txt !== null) e.textContent = txt;
  return e;
}
function svgEl(t, attrs){
  const e = document.createElementNS('http://www.w3.org/2000/svg', t);
  if(attrs){ for(const k in attrs) e.setAttribute(k, attrs[k]); }
  return e;
}
// ITPFS Ver1.0(修正): ダーク/ライト切替に47都道府県マップの陸地・グリッド線色を
// 連動させるためのヘルパー。SVG要素のfill/stroke属性は素のCSSプロパティでは
// ないため、var(--x)を直接指定できず、都度getComputedStyleで実値を読む必要がある。
function cssVar(name, fallback){
  try{
    const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fallback;
  }catch(_){ return fallback; }
}
// Ver13.11: SVG テキストに白ハロー (縁取り) を付ける汎用ヘルパー。
// paint-order:stroke で先にストロークを描き、その上に fill を重ねる。
// これで海色・陸色に関係なく文字が読み取れるようになる。
function svgTextHalo(attrs, text, opts){
  opts = opts || {};
  const fill = opts.fill || '#111';
  const halo = opts.halo || '#ffffff';
  const haloWidth = opts.haloWidth || 3;
  const t = svgEl('text', Object.assign({}, attrs, {
    fill: fill,
    stroke: halo,
    'stroke-width': haloWidth,
    'stroke-linejoin': 'round',
    'stroke-linecap': 'round',
    'paint-order': 'stroke',
  }));
  t.textContent = text;
  return t;
}
function fmtNum(v, nd){
  if(v === null || v === undefined || Number.isNaN(v)) return '--';
  return Number(v).toFixed(nd ?? 1);
}
// Rev15.8③: タイムスライダーで選択中のホライズン(+0/15/30/45/60/90分)が
// 実際に何時何分の予測かを小さく自動表記するためのJST時刻ラベルを計算する。
// サーバの絶対時刻(cycle_info.server_epoch、UTC epoch秒)を基準に、ブラウザの
// ローカルタイムゾーン設定に依存せず常にJST(UTC+9)で表示する。
function updateHorizonTimeLabel(){
  const targets = document.querySelectorAll('.horizon-time-lbl');
  if(!targets.length || !latestSnapshot || !latestSnapshot.cycle_info) return;
  const serverEpoch = latestSnapshot.cycle_info.server_epoch;
  let text = '';
  if(serverEpoch){
    const targetUtcMs = (serverEpoch + currentHorizon*60) * 1000;
    const jstMs = targetUtcMs + 9*3600*1000;
    const d = new Date(jstMs);
    const hh = String(d.getUTCHours()).padStart(2,'0');
    const mm = String(d.getUTCMinutes()).padStart(2,'0');
    text = '(' + hh + ':' + mm + ' JST頃)';
  }
  targets.forEach(t=>{ t.textContent = text; });
}
function pctBar(pct, warn){
  const w = el('div', 'bar' + (warn ? ' warn' : ''));
  const sp = el('span');
  sp.style.width = Math.max(0, Math.min(100, pct||0)) + '%';
  w.appendChild(sp);
  return w;
}
function pctColor(pct){
  if(pct === null || pct === undefined) return '#6b7899';
  if(pct >= 75) return '#ff5252';
  if(pct >= 50) return '#26c281';
  if(pct >= 30) return '#ffb020';
  if(pct >= 15) return '#8b7cff';
  return '#4b587a';
}
// Rev16.2: 「簡単表示」専用の3段階トラフィックライト信号。ユーザー要望により、
// 詳細な状態ラベル(state-hot/good/warn/weak/low等)をそのまま見せるのではなく、
// 赤(期待できる)/オレンジ(兆候あり)/青(静か)の3色シグナルへ単純化し、
// 数値を読まなくても一目で判断できるようにする(詳しい表示側は引き続き
// 数値・詳細ラベルのまま、pctColor()/state-*クラスを使用する)。
function trafficLight(pct){
  if(pct === null || pct === undefined) return {emoji:'⚪', label:'まだ判定できていません', color:'#6b7899'};
  if(pct >= 50) return {emoji:'🔴', label:'期待できます', color:'#ff5252'};
  if(pct >= 25) return {emoji:'🟠', label:'兆候があります', color:'#ffb020'};
  return {emoji:'🔵', label:'今は静かです', color:'#3498db'};
}
// Rev15.9: 初心者/上級者の表示モード切替。localStorageで選択を記憶する。
// ITPFS Ver1.0(修正): 「duct」を第3のタブ状態として追加(既存のbeginner/advancedと
// 同列)。ページ遷移は発生せず、#root(Es/F2側)を隠して#duct-ui-persistentを
// 表示するだけなので、タブ切替ボタン自体は常に残る。
const TIER_STORAGE_KEY = 'edfs_ui_tier';
function getUiTier(){
  try{
    const v = localStorage.getItem(TIER_STORAGE_KEY);
    return (v === 'advanced' || v === 'duct') ? v : 'beginner';
  }catch(e){ return 'beginner'; }
}
function setUiTier(tier){
  try{ localStorage.setItem(TIER_STORAGE_KEY, tier); }catch(e){}
  document.body.classList.toggle('mode-advanced', tier === 'advanced');
  document.body.classList.toggle('mode-duct', tier === 'duct');
  document.body.classList.toggle('mode-beginner', tier !== 'advanced' && tier !== 'duct');
  document.querySelectorAll('.tier-btn[data-tier="beginner"], .tier-btn[data-tier="advanced"], .tier-btn[data-tier="duct"]')
    .forEach(b=>{ b.classList.toggle('active', b.dataset.tier === tier); });
  // ITPFS(整理版・大規模UI改修): ダクトタブへ切り替えた瞬間に、その時点で
  // 分かっている最新状態へ地図データを再取得する(embedded scriptのduct2FetchAndRender
  // 参照。ページ読み込み直後にはまだ未定義のため、存在チェックしてから呼ぶ)。
  if(tier === 'duct' && window._ductMapForceReload){ window._ductMapForceReload(); }
  // ITPFS(整理版・バグ修正): 祖先要素(#duct-ui-persistent)が非表示の間に
  // Leafletが誤ったサイズで初期化される既知動作への対策として、タブが
  // 可視化された直後にも明示的にinvalidateSize()を呼ぶ(念のための二重対策。
  // 詳細はembedded script内のinitLeafletMapOnce()コメント参照)。
  if(tier === 'duct' && window._ductMapInvalidateSize){
    setTimeout(window._ductMapInvalidateSize, 60);
    setTimeout(window._ductMapInvalidateSize, 400);
  }
}
function initTierSwitch(){
  const tier = getUiTier();
  setUiTier(tier);
  document.querySelectorAll('.tier-btn[data-tier="beginner"], .tier-btn[data-tier="advanced"], .tier-btn[data-tier="duct"]')
    .forEach(b=>{ b.addEventListener('click', ()=> setUiTier(b.dataset.tier)); });
}
// Rev15.9/16.2: 初心者向けの要約カード。Es(方面別最有力)・F2(遠隔方面最有力)・
// 位置情報・Esシーズンの4点を、数値だけでなく平易な一言で提示する。
// Rev16.2の変更点:
//   ①Es/F2ともに「方面別予測」(areas_by_horizon)を共通データソースとして  
//     使うことで、タイムスライダー(+0/15/30/45/60/90分)の操作に連動して  
//     内容が変化するようにした(従来はEsが常に代表局の+0分値に固定、F2も  
//     +0分固定で、地図以外の表示がスライダーに追従していなかった)。
//   ②開通期待の表現を、詳細な5段階状態ラベルではなく赤/オレンジ/青の3色
//     トラフィックライト信号(trafficLight())に統一し、数値を読まなくても
//     一目で判断できるようにした。
function renderBeginnerSummary(data){
  const card = el('div', 'card tier-beginner-only');
  card.appendChild(el('h2', null, '📡 かんたん予報'));
  const body = el('div'); body.id = 'beginner-summary-body';
  card.appendChild(body);
  refreshBeginnerSummaryBody(data);
  card.appendChild(el('div', 'hint', 'もっと詳しい情報は上部の「🔧 くわしい表示」を選んでください。'));
  return card;
}
// Rev16.2: タイムスライダー操作時にも呼び出され、#beginner-summary-body の
// 中身だけを差し替える(「かんたん予報」カード自体は作り直さずちらつきを防止)。
function refreshBeginnerSummaryBody(data){
  const body = document.getElementById('beginner-summary-body');
  if(!body) return;
  body.innerHTML = '';
  const areasNow = (data.map && data.map.areas_by_horizon
    && data.map.areas_by_horizon[String(currentHorizon)]) || data.areas || [];
  // Es: 方面別予測の中から選択中ホライズンで最も期待度が高い方面を探す
  // (F2側と同じ「方面ベース」の切り口に統一し、スライダーとの連動を保証する)
  let bestEs = null;
  areasNow.forEach(a=>{
    if(a.display_pct===null || a.display_pct===undefined) return;
    if(!bestEs || a.display_pct > bestEs.pct) bestEs = {name:a.display_name, pct:a.display_pct};
  });
  const row1 = el('div');
  row1.style.display='flex'; row1.style.alignItems='center'; row1.style.gap='10px';
  const tlEs = trafficLight(bestEs ? bestEs.pct : null);
  row1.appendChild(el('span', null, tlEs.emoji));
  row1.appendChild(el('div', null,
    'Es(スポラディックE層): ' + tlEs.label
    + (bestEs ? '（' + bestEs.name + ' 方面、約' + Math.round(bestEs.pct) + '%）' : '（データ収集中）')));
  body.appendChild(row1);
  // F2: 方面別予測(選択中ホライズン)の中から最有力を探す
  let bestF2 = null;
  areasNow.forEach(a=>{
    if(a.f2_pct===null || a.f2_pct===undefined) return;
    if(!bestF2 || a.f2_pct > bestF2.pct) bestF2 = {name:a.display_name, pct:a.f2_pct};
  });
  const row2 = el('div');
  row2.style.display='flex'; row2.style.alignItems='center'; row2.style.gap='10px';
  row2.style.marginTop='8px';
  const tlF2 = trafficLight(bestF2 ? bestF2.pct : null);
  row2.appendChild(el('span', null, tlF2.emoji));
  row2.appendChild(el('div', null,
    'F2層(遠距離): ' + tlF2.label
    + (bestF2 ? '（' + bestF2.name + ' 方面、約' + Math.round(bestF2.pct) + '%）' : '（データ収集中、または対象方面なし）')));
  body.appendChild(row2);
  // 選択中の予測時刻(スライダー連動。.horizon-time-lblはupdateHorizonTimeLabel()が一括更新)
  const timeRow = el('div', 'sub small');
  timeRow.style.marginTop='4px';
  timeRow.appendChild(document.createTextNode(
    (currentHorizon===0 ? '現在の状態です ' : '+' + currentHorizon + '分後の予測です ')));
  timeRow.appendChild(el('span', 'horizon-time-lbl'));
  timeRow.appendChild(document.createTextNode('（上部地図のスライダーで変更できます）'));
  body.appendChild(timeRow);
  // 位置情報
  const row3 = el('div');
  row3.style.display='flex'; row3.style.alignItems='center'; row3.style.gap='10px';
  row3.style.marginTop='8px';
  if(data.location){
    row3.appendChild(el('span', null, '📍'));
    row3.appendChild(el('div', 'sub',
      '現在地: 緯度' + fmtNum(data.location.lat,3) + ' / 経度' + fmtNum(data.location.lon,3)
      + '（' + (data.location.source||data.location.label||'取得済み') + '）'));
  } else {
    row3.appendChild(el('span', null, '⚠'));
    row3.appendChild(el('div', 'state-warn', '位置情報が未取得です。予測精度が下がる場合があります。'));
    const locBtn = el('button', 'btn ghost', '📍 位置を取得');
    locBtn.style.marginLeft='auto';
    locBtn.onclick = ()=>{ lastGpsSent = 0; reqGps(); };
    row3.appendChild(locBtn);
  }
  body.appendChild(row3);
  // Esシーズン状態
  const season = data.season || {};
  const row4 = el('div', 'sub small');
  row4.style.marginTop='8px';
  row4.textContent = 'Esシーズン状態: ' + (season.label || '判定中');
  body.appendChild(row4);
}
function sendGps(lat, lon, source){
  const now = Date.now();
  if(now - lastGpsSent < 30000) return;
  lastGpsSent = now;
  fetch('/api/gps', {method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({lat, lon, source: source||'browser'})}).catch(()=>{});
}
// Rev15.9: ブラウザのGeolocation APIは「secure context」(HTTPS、または
// http://localhost / http://127.0.0.1 のみ)でしか動作しない仕様がある。
// Chromebook(ARC)ではダッシュボードをLAN/ARCブリッジIP経由の素のhttpで
// 開く必要があり、この制約に該当してGPSボタンが常に無反応(無言で失敗)になる。
// ここで理由を明示的に検知し、UI側で手動入力欄への誘導を出せるようにする。
window._gpsLastError = null;
function reqGps(){
  if(!('geolocation' in navigator)){
    window._gpsLastError = 'このブラウザはGeolocation APIに対応していません';
    return;
  }
  if(window.isSecureContext === false){
    window._gpsLastError = 'insecure-context';
    return;
  }
  navigator.geolocation.getCurrentPosition(
    (p)=>{ window._gpsLastError = null; sendGps(p.coords.latitude, p.coords.longitude, 'browser-geolocation'); },
    (err)=>{ window._gpsLastError = (err && err.message) ? err.message : 'GPS取得に失敗しました'; },
    {enableHighAccuracy:true, timeout:8000, maximumAge:300000});
}
async function flushLearning(btn){
  btn.disabled = true;
  const orig = btn.textContent;
  btn.textContent = '保存中…';
  try{
    const r = await fetch('/api/flush', {method:'POST'});
    const d = await r.json();
    btn.textContent = d.ok ? '✅ 保存しました' : '⚠️ 失敗: ' + (d.error||'');
  } catch(e){ btn.textContent = '⚠️ 通信失敗'; }
  setTimeout(()=>{ btn.textContent = orig; btn.disabled = false; }, 3000);
}
function renderStartup(pct, msg){
  const s = document.getElementById('startup');
  if(!s) return;
  document.getElementById('startup-pct').textContent = (pct||5) + '%';
  document.getElementById('startup-bar').style.width = (pct||5) + '%';
  document.getElementById('startup-msg').textContent = msg || '初期化中';
}

// ==================== 日本地図 (Ver13.11: 広域ビューポート) ====================
// Es 伝搬は日本より広い範囲 (1000-2500km 程度)。 中国大陸〜西太平洋の
// 広域ビューポートで海域を大きく取り、Japan は真ん中に小さく配置する。
// これにより:
//   1) Es 伝搬エリアが視覚化される
//   2) AI シルエットの細かな歪みが目立たなくなる
//   3) 全マーカーが正確な緯度経度で表示される (Okinawa 等の実座標維持)
// Ver13.11: ビューポートを日本+近接海域に絞り込み (地図を大きく表示)
// 中国大陸大部分を切り、Es 単ホップ範囲 (~1500km) は縁まで収まる
const MAP_W = 1000, MAP_H = 1140;
// 地理範囲: 緯度 22°N〜46.5°N (24.5°), 経度 122°E〜148°E (26°)
const MAP_LAT_MIN = 22, MAP_LAT_MAX = 46.5;
const MAP_LON_MIN = 122, MAP_LON_MAX = 148;
// Ver13.11: PNG不使用。 JAPAN_COAST ポリゴン (実座標) で正確に描画。
function lat2y(lat){ return MAP_H * (1 - (lat - MAP_LAT_MIN)/(MAP_LAT_MAX - MAP_LAT_MIN)); }
function lon2x(lon){ return MAP_W * (lon - MAP_LON_MIN)/(MAP_LON_MAX - MAP_LON_MIN); }

// Ver13.5: Es 移動の推定現在位置を計算 (時間同期シミュレーション)
function computeEsCurrentPosition(mapData, elapsedSec){
  const c = mapData.es_centroid || {};
  if(c.lat === null || c.lat === undefined) return null;
  // 速度が無い/停滞なら固定位置
  let dir_deg = null, speed = null;
  if(c.mstid_dir_deg !== null && c.mstid_dir_deg !== undefined
     && c.mstid_speed_kmh !== null && c.mstid_speed_kmh !== undefined){
    dir_deg = c.mstid_dir_deg;
    speed = c.mstid_speed_kmh;
  } else if(c.speed_kmh && c.direction && c.direction !== '不明' && c.direction !== '停滞'){
    speed = c.speed_kmh;
    dir_deg = (c.direction === '北上') ? 0 : 180;
  }
  if(dir_deg === null || speed === null || speed < 1){
    return {lat: c.lat, lon: c.lon, elapsed_min: elapsedSec / 60};
  }
  const dt_h = elapsedSec / 3600;
  const dist_km = speed * dt_h;
  const th = dir_deg * Math.PI / 180;
  const dlat = dist_km * Math.cos(th) / 111.0;
  const kmPerLon = 111 * Math.cos(c.lat * Math.PI / 180);
  const dlon = dist_km * Math.sin(th) / kmPerLon;
  return {lat: c.lat + dlat, lon: c.lon + dlon, elapsed_min: elapsedSec / 60};
}

// Ver13.11: Catmull-Rom スプライン → Cubic Bezier で滑らかな海岸線に変換
// 閉じたポリゴンの頂点列 pts を SVG path 文字列 (M/C/... Z) に変換
function smoothClosedPath(pts, tension){
  if(!pts || pts.length < 3) return '';
  const T = tension === undefined ? 0.5 : tension;
  const N = pts.length;
  // 各頂点を SVG 座標に事前変換
  const P = pts.map(pt => [lon2x(pt[0]), lat2y(pt[1])]);
  let d = 'M ' + P[0][0].toFixed(1) + ',' + P[0][1].toFixed(1) + ' ';
  for(let i=0; i<N; i++){
    const p0 = P[(i - 1 + N) % N];
    const p1 = P[i];
    const p2 = P[(i + 1) % N];
    const p3 = P[(i + 2) % N];
    // Catmull-Rom → Cubic Bezier 制御点
    const c1x = p1[0] + (p2[0] - p0[0]) * T / 3;
    const c1y = p1[1] + (p2[1] - p0[1]) * T / 3;
    const c2x = p2[0] - (p3[0] - p1[0]) * T / 3;
    const c2y = p2[1] - (p3[1] - p1[1]) * T / 3;
    d += 'C ' + c1x.toFixed(1) + ',' + c1y.toFixed(1) + ' '
             + c2x.toFixed(1) + ',' + c2y.toFixed(1) + ' '
             + p2[0].toFixed(1) + ',' + p2[1].toFixed(1) + ' ';
  }
  d += 'Z';
  return d;
}

function drawMap(mapData, horizonMin){
  const svg = svgEl('svg', {viewBox:`0 0 ${MAP_W} ${MAP_H}`,
    xmlns:'http://www.w3.org/2000/svg', preserveAspectRatio:'xMidYMid meet'});
  // 海 (淡い青灰)
  svg.appendChild(svgEl('rect', {x:0, y:0, width:MAP_W, height:MAP_H, fill:cssVar('--sea','#e8ecf3')}));

  // Ver13.11: 精緻化ポリゴン (389頂点、14陸地) を直線で明瞭に描画。
  // Bezier で丸めすぎる違和感を排除し、岬や湾のシャープさを保持する。
  // JAPAN_COAST の全キーを描画するので、離島 (佐渡・淡路・対馬・種子・屋久・
  // 奄美・宮古・石垣) や近隣 (済州島) も自動で表示される。
  const gLand = svgEl('g', {fill:cssVar('--map-land','#4a5670'), stroke:cssVar('--map-land-stroke','#1e2838'),
    'stroke-width':1.2, 'stroke-linejoin':'round'});
  Object.keys(JAPAN_COAST).forEach(key=>{
    const pts = JAPAN_COAST[key];
    if(!pts || pts.length < 3) return;
    // 直線ポリゴン (Bezierなし) で明瞭な海岸線
    let d = 'M ';
    pts.forEach((pt, i) => {
      const x = lon2x(pt[0]).toFixed(1);
      const y = lat2y(pt[1]).toFixed(1);
      d += (i === 0 ? '' : 'L ') + x + ',' + y + ' ';
    });
    d += 'Z';
    const path = svgEl('path', {d});
    path.setAttribute('data-tip', key);
    gLand.appendChild(path);
  });
  svg.appendChild(gLand);

  // グリッド (Ver13.11: 線は薄く、ラベルは白ハロー付き濃色で確実に読める)
  const gGrid = svgEl('g', {stroke:cssVar('--map-grid','#9aa8c0'), 'stroke-width':0.7,
    fill:'none', opacity:0.5});
  for(let lat=25; lat<=45; lat+=5){
    const y = lat2y(lat);
    gGrid.appendChild(svgEl('line', {x1:0, y1:y, x2:MAP_W, y2:y,
      'stroke-dasharray':'5,5'}));
  }
  for(let lon=125; lon<=145; lon+=5){
    const x = lon2x(lon);
    gGrid.appendChild(svgEl('line', {x1:x, y1:0, x2:x, y2:MAP_H,
      'stroke-dasharray':'5,5'}));
  }
  svg.appendChild(gGrid);
  // グリッドラベルは別レイヤ (地図描画後に上書き)
  const gGridLbl = svgEl('g');
  for(let lat=25; lat<=45; lat+=5){
    const y = lat2y(lat);
    gGridLbl.appendChild(svgTextHalo({x:8, y:y-4, 'font-size':14, 'font-weight':'600'},
      lat+'°N', {fill:'#2a3548', halo:'#ffffff', haloWidth:3}));
  }
  for(let lon=125; lon<=145; lon+=5){
    const x = lon2x(lon);
    gGridLbl.appendChild(svgTextHalo({x:x+4, y:MAP_H-8, 'font-size':14, 'font-weight':'600'},
      lon+'°E', {fill:'#2a3548', halo:'#ffffff', haloWidth:3}));
  }
  svg.appendChild(gGridLbl);

  // 矢印マーカー
  const defs = svgEl('defs');
  ['#ffcccc','#ffd700','#8b7cff','#26c281'].forEach((c,i)=>{
    const m = svgEl('marker', {id:'ah'+i, markerWidth:8, markerHeight:8,
      refX:6, refY:4, orient:'auto'});
    m.appendChild(svgEl('polygon', {points:'0,0 8,4 0,8', fill:c}));
    defs.appendChild(m);
  });
  svg.appendChild(defs);

  // --- Es 電離雲 (Ver13.5: 暫定表示対応) ---
  const c = mapData.es_centroid || {};
  const provisional = c.provisional === true;
  if(c.lat !== null && c.lat !== undefined){
    // ホライズン分先の推定位置 (静的スライダ表示用)
    const targetPos = computeEsCurrentPosition(mapData, horizonMin * 60);
    const cx = lon2x(targetPos ? targetPos.lon : c.lon);
    const cy = lat2y(targetPos ? targetPos.lat : c.lat);
    const targetLat = targetPos ? targetPos.lat : c.lat;

    // 形状スコアから歪み係数
    let stretch_ew = 1.0, stretch_ns = 0.5, rotate = 0;
    let dominant = provisional ? '暫定' : '不明';
    const stations = mapData.stations || [];
    const agg = {};
    let wsum = 0;
    stations.forEach(s=>{
      const w = Math.max(0.1, (s.fxEs || 3) - 3);
      wsum += w;
      const sc = s.local_shape_scores || {};
      for(const k in sc){ agg[k] = (agg[k] || 0) + (sc[k] || 0) * w; }
      if(s.local_shape && !provisional) dominant = s.local_shape;
    });
    const arch = (agg['アーチ'] || agg['arch'] || 0) / (wsum || 1);
    const band = (agg['帯状'] || agg['band'] || 0) / (wsum || 1);
    const wave = (agg['波状'] || agg['wave'] || 0) / (wsum || 1);
    const fault = (agg['断層'] || agg['fault'] || 0) / (wsum || 1);
    const blob = (agg['塊'] || agg['blob'] || 0) / (wsum || 1);
    if(band > 0.3){ stretch_ew = 1.8; stretch_ns = 0.4; }
    if(arch > 0.3){ stretch_ew = 1.5; stretch_ns = 0.6; }
    if(fault > 0.3){ stretch_ew = 0.6; stretch_ns = 1.4; rotate = 30; }
    if(blob > 0.3){ stretch_ew = 1.0; stretch_ns = 0.9; }
    // Ver13.21: 実測FT8伝播方位の統計(bearing_axial_r/orientation_deg)が十分な
    // サンプルを持っていれば、局所形状分類より優先して実測角度で傾きを決める。
    // (orientation_degは0=南北,90=東西。楕円の初期長軸はrx側=東西のため変換する)
    if(c.anisotropy_r !== null && c.anisotropy_r !== undefined && c.anisotropy_r >= 0.15
       && c.anisotropy_orientation_deg !== null && c.anisotropy_orientation_deg !== undefined){
      rotate = ((90 - c.anisotropy_orientation_deg) % 180 + 180) % 180;
      if(c.anisotropy_r >= 0.40){
        stretch_ew = Math.max(stretch_ew, 1.6);
        stretch_ns = Math.min(stretch_ns, 0.45);
      }
    }

    const kmPerLon = 111 * Math.cos(targetLat * Math.PI/180);
    const rx = (800 * stretch_ew / 2) / kmPerLon * (MAP_W/(MAP_LON_MAX-MAP_LON_MIN));
    const ry = (400 * stretch_ns / 2) / 111 * (MAP_H/(MAP_LAT_MAX-MAP_LAT_MIN));

    const fillColor = provisional ? 'rgba(155,155,155,0.15)' : 'rgba(255,82,82,0.15)';
    const strokeColor = provisional ? 'rgba(180,180,180,0.7)' : 'rgba(255,82,82,0.7)';

    if(wave > 0.3){
      let d = '';
      const N = 40;
      for(let i=0; i<=N; i++){
        const t = i/N * 2*Math.PI;
        const px = cx + rx * Math.cos(t);
        const py = cy + ry * Math.sin(t) + 15 * Math.sin(3*t);
        d += (i===0 ? 'M ' : 'L ') + px.toFixed(1) + ',' + py.toFixed(1) + ' ';
      }
      d += 'Z';
      svg.appendChild(svgEl('path', {d, fill:fillColor,
        stroke:strokeColor, 'stroke-width':2, 'stroke-dasharray':'6,4'}));
    } else {
      svg.appendChild(svgEl('ellipse', {cx, cy, rx, ry,
        fill:fillColor, stroke:strokeColor,
        'stroke-width':2, 'stroke-dasharray':'6,4',
        transform: `rotate(${rotate} ${cx} ${cy})`}));
    }
    // 中心マーカー (Ver13.5: リアルタイム移動のため id 付き)
    const centerColor = provisional ? '#aaa' : '#ff5252';
    const centerDot = svgEl('circle', {cx, cy, r:6, fill:centerColor,
      stroke:'#fff', 'stroke-width':1.5, id:'es-center-dot'});
    centerDot.setAttribute('data-tip', provisional ? '暫定 Es 位置' :
      `Es 重心 [${dominant}]\\n${horizonMin}分後の推定位置`);
    svg.appendChild(centerDot);

    const lblText = provisional
      ? `Es 位置 (暫定・観測データ収集中)`
      : `Es 重心 [${dominant}] +${horizonMin}分`;
    // Ver13.11: ハロー付きで背景色に依存せず読める
    const esColor = provisional ? '#666' : '#c22';
    const lbl = svgTextHalo({x:cx+10, y:cy-8, 'font-size':14,
      'font-weight':'bold', id:'es-center-label'},
      lblText, {fill:esColor, halo:'#ffffff', haloWidth:3.5});
    svg.appendChild(lbl);

    // Ver13.5: リアルタイム移動マーカー (別マーカー、JS で 1秒毎に位置更新)
    // 現在の推定位置 (前サイクルからの経過時間分だけ進んだ位置)
    if(!provisional && !c.provisional){
      const liveDot = svgEl('circle', {cx, cy, r:8, fill:'none',
        stroke:'#0e8a55', 'stroke-width':2.8, id:'es-live-dot',
        'stroke-dasharray':'4,3'});
      liveDot.style.transformOrigin = 'center';
      svg.appendChild(liveDot);
      const liveLbl = svgTextHalo({x:cx+14, y:cy+20, 'font-size':12,
        'font-weight':'bold', id:'es-live-label'},
        '● 現在推定位置', {fill:'#0e8a55', halo:'#ffffff', haloWidth:3});
      svg.appendChild(liveLbl);
    }

    // 移動ベクトル (30分先まで矢印)
    if(!provisional){
      const dt_h = horizonMin / 60.0;
      let dlat = 0, dlon = 0;
      if(c.mstid_dir_deg !== null && c.mstid_dir_deg !== undefined
         && c.mstid_speed_kmh !== null && c.mstid_speed_kmh !== undefined){
        const th = c.mstid_dir_deg * Math.PI/180;
        const dist_km = c.mstid_speed_kmh * dt_h;
        dlat = dist_km * Math.cos(th) / 111.0;
        dlon = dist_km * Math.sin(th) / kmPerLon;
      } else if(c.speed_kmh && c.direction && c.direction !== '不明' && c.direction !== '停滞'){
        const sign = (c.direction === '北上') ? 1 : -1;
        dlat = sign * c.speed_kmh * dt_h / 111.0;
      }
      if(Math.abs(dlat) > 0.01 || Math.abs(dlon) > 0.01){
        // ベースからターゲットまでの矢印
        const cx0 = lon2x(c.lon), cy0 = lat2y(c.lat);
        const x2 = lon2x(c.lon + dlon), y2 = lat2y(c.lat + dlat);
        const arrow = svgEl('line', {x1:cx0, y1:cy0, x2, y2,
          stroke:'#ffcccc', 'stroke-width':2, opacity:0.6});
        arrow.setAttribute('marker-end', 'url(#ah0)');
        svg.appendChild(arrow);
      }
    }
  }

  // --- Ver13.21/22: PSKReporter実測伝播局(片方向矢印) — 推定Es帯(破線)と明確に区別。
  // 注意: PSKReporterへの問い合わせは自局送信→相手局受信(TXレポート)のみを取得しており
  // (senderCallsign=自局, rronly=1)、逆方向(自局が受信したレポート)は現状取得していない。
  // したがって矢印は常に「自局→相手」の片方向として描画し、双方向を装わない。
  const pskPts = mapData.psk_points || [];
  const myLoc = mapData.my_location;
  if(pskPts.length && myLoc && myLoc.lat && myLoc.lon){
    const mx = lon2x(myLoc.lon), my_ = lat2y(myLoc.lat);
    const gPsk = svgEl('g', {opacity:0.75});
    pskPts.forEach(p=>{
      if(p.lat === null || p.lon === null || p.lat === undefined || p.lon === undefined) return;
      const px = lon2x(p.lon), py = lat2y(p.lat);
      const ln = svgEl('line', {x1:mx, y1:my_, x2:px, y2:py,
        stroke:'#26c281', 'stroke-width':1, opacity:0.45});
      ln.setAttribute('marker-end', 'url(#ah3)');
      gPsk.appendChild(ln);
      gPsk.appendChild(svgEl('circle', {cx:px, cy:py, r:2.5,
        fill:'#26c281', stroke:'#0a3', 'stroke-width':0.6}));
    });
    svg.appendChild(gPsk);
  }

  // --- NICT 4局 (三角、Ver13.11: 白ハロー付き濃色ラベルで背景色に依存せず可読) ---
  (mapData.stations||[]).forEach(s=>{
    if(!s.lat || !s.lon) return;
    const x = lon2x(s.lon), y = lat2y(s.lat);
    const size = 13;
    const color = pctColor(s.p27);
    const tri = svgEl('polygon', {
      points:`${x},${y-size} ${x-size},${y+size} ${x+size},${y+size}`,
      fill:color, stroke:'#000000', 'stroke-width':1.8});
    tri.style.cursor='pointer';
    tri.setAttribute('data-tip',
      `NICT ${s.name}\\nfxEs: ${fmtNum(s.fxEs,1)} MHz\\n30分後: ${fmtNum(s.fx30,1)} MHz\\n開通率: ${s.p27===null?'--':Math.round(s.p27)+'%'}\\n形状: ${s.local_shape||'--'}`);
    svg.appendChild(tri);
    // 局名 (濃黒 + 白ハロー)
    svg.appendChild(svgTextHalo(
      {x:x+16, y:y+5, 'font-size':16, 'font-weight':'bold'},
      s.name, {fill:'#0a0a0a', halo:'#ffffff', haloWidth:4}));
    // fxEs 値 (濃青 + 白ハロー)
    svg.appendChild(svgTextHalo(
      {x:x+16, y:y+25, 'font-size':13, 'font-weight':'600'},
      fmtNum(s.fxEs,1)+' MHz',
      {fill:'#1e3a8a', halo:'#ffffff', haloWidth:3}));
  });

  // --- エリア別予測 (Ver13.11: 数値を白ハロー付き濃紺で表示) ---
  const areas = (mapData.areas_by_horizon || {})[String(horizonMin)] || [];
  areas.forEach(a=>{
    const x = lon2x(a.lon), y = lat2y(a.lat);
    const p = a.display_pct;
    if(p === null || p === undefined) return;
    const r = 8 + (p/100) * 22;
    const color = pctColor(p);
    const circ = svgEl('circle', {cx:x, cy:y, r,
      fill:color, 'fill-opacity':0.55, stroke:color, 'stroke-width':2.2});
    circ.style.cursor='pointer';
    circ.setAttribute('data-tip',
      `${a.display_name}\\n+${horizonMin}分後: ${Math.round(p)}%` +
      (a.uncertainty_pct !== null && a.uncertainty_pct !== undefined ? ` ±${a.uncertainty_pct}%` : '') +
      `\\n現在到達: ${a.current_pct}%\\n傾向: ${a.trend_label}` +
      (a.skill_score !== null && a.skill_score !== undefined
        ? `\\nSkill(参考): ${a.skill_score >= 0 ? '+' : ''}${a.skill_score.toFixed(2)}` : '') +
      (a.state_coverage_pct !== null && a.state_coverage_pct !== undefined
        ? `\\n状態経験度: ${a.state_coverage_pct}%` : ''));
    svg.appendChild(circ);
    // 数値ラベル (白ハロー付き濃黒で可読性向上)
    svg.appendChild(svgTextHalo(
      {x, y:y+4, 'font-size':13, 'text-anchor':'middle', 'font-weight':'bold'},
      Math.round(p)+'%',
      {fill:'#0a0a0a', halo:'#ffffff', haloWidth:3.5}));
  });

  // --- 自局 (★、Ver13.11: 白ハロー付きラベル) ---
  const my = mapData.my_location;
  if(my && my.lat && my.lon){
    const x = lon2x(my.lon), y = lat2y(my.lat);
    const star = svgEl('polygon', {
      points:`${x},${y-14} ${x+4},${y-4} ${x+14},${y-4} ${x+6},${y+3} ${x+9},${y+13} ${x},${y+7} ${x-9},${y+13} ${x-6},${y+3} ${x-14},${y-4} ${x-4},${y-4}`,
      fill:'#ffd700', stroke:'#000', 'stroke-width':1.8});
    star.style.cursor='pointer';
    star.setAttribute('data-tip', `自局\\n${fmtNum(my.lat,3)}, ${fmtNum(my.lon,3)}\\n(${my.source||''})`);
    svg.appendChild(star);
    svg.appendChild(svgTextHalo(
      {x:x+18, y:y+5, 'font-size':15, 'font-weight':'bold'},
      'あなた', {fill:'#7a5a00', halo:'#ffffff', haloWidth:4}));
  }

  return svg;
}

// Ver13.5: リアルタイム位置更新 (1秒周期)
function updateEsLiveMarker(){
  if(!latestSnapshot || !latestSnapshot.map) return;
  const c = latestSnapshot.map.es_centroid;
  if(!c || c.provisional || c.lat === null || c.lat === undefined) return;
  const dot = document.getElementById('es-live-dot');
  const lbl = document.getElementById('es-live-label');
  if(!dot) return;
  // サーバ epoch と現在時刻から経過秒を算出
  const nowClient = Date.now() / 1000;
  const serverEpoch = c.server_epoch || nowClient;
  const drift = serverEpoch - (nowClient - (Date.now() - (window._snapArrivedAt || Date.now()))/1000);
  const elapsed = nowClient - serverEpoch + drift;
  const pos = computeEsCurrentPosition(latestSnapshot.map, Math.max(0, elapsed));
  if(pos){
    const x = lon2x(pos.lon);
    const y = lat2y(pos.lat);
    dot.setAttribute('cx', x);
    dot.setAttribute('cy', y);
    if(lbl){
      lbl.setAttribute('x', x + 14);
      lbl.setAttribute('y', y + 18);
      const mins = Math.floor(elapsed / 60);
      lbl.textContent = `● 現在推定 (前サイクル+${mins}分)`;
    }
  }
}

function bindTips(container){
  const tip = el('div', 'map-tooltip'); tip.id = 'maptip';
  container.appendChild(tip);
  container.addEventListener('mouseover', (e)=>{
    const t = e.target && e.target.getAttribute && e.target.getAttribute('data-tip');
    if(t){ tip.textContent = t.replace(/\\n/g,' / '); tip.style.display='block'; }
  });
  container.addEventListener('mousemove', (e)=>{
    const rect = container.getBoundingClientRect();
    tip.style.left = (e.clientX - rect.left + 12) + 'px';
    tip.style.top = (e.clientY - rect.top + 12) + 'px';
  });
  container.addEventListener('mouseout', ()=>{ tip.style.display='none'; });
  container.addEventListener('touchstart', (e)=>{
    const t = e.target && e.target.getAttribute && e.target.getAttribute('data-tip');
    if(t){
      tip.textContent = t.replace(/\\n/g,' / ');
      tip.style.display='block';
      const rect = container.getBoundingClientRect();
      const touch = e.touches[0];
      tip.style.left = (touch.clientX - rect.left + 12) + 'px';
      tip.style.top = (touch.clientY - rect.top + 12) + 'px';
      setTimeout(()=>{ tip.style.display='none'; }, 3000);
    }
  });
}

// ==================== 履歴グラフ ====================
function drawHistoryChart(history){
  const ts = history.timestamps || [];
  if(ts.length < 2){
    const p = el('div', 'sub'); p.textContent = '履歴データ収集中 (2サンプル以上必要)';
    return p;
  }
  const W = 1000, H = 240, PAD_L = 40, PAD_R = 10, PAD_T = 20, PAD_B = 30;
  const svg = svgEl('svg', {viewBox:`0 0 ${W} ${H}`, xmlns:'http://www.w3.org/2000/svg',
    preserveAspectRatio:'xMidYMid meet'});
  svg.appendChild(svgEl('rect', {x:0, y:0, width:W, height:H, fill:'#0b1428'}));
  let maxObsVal = 0;
  ts.forEach((tv,i)=>{
    const vo = history.fxes_obs[i], vp = history.fx30_pred[i];
    if(vo !== null && vo !== undefined && vo > maxObsVal) maxObsVal = vo;
    if(vp !== null && vp !== undefined && vp > maxObsVal) maxObsVal = vp;
  });
  const yStep = 3;
  const yMax = Math.max(yStep, Math.ceil((maxObsVal * 1.15) / yStep) * yStep);
  const t0 = ts[0], tN = ts[ts.length-1];
  const spanT = Math.max(60, tN - t0);
  function x(tv){ return PAD_L + (W-PAD_L-PAD_R) * (tv-t0)/spanT; }
  function y(v){ return H - PAD_B - (H-PAD_T-PAD_B) * (v/yMax); }
  for(let v=0; v<=yMax; v+=yStep){
    svg.appendChild(svgEl('line', {x1:PAD_L, y1:y(v), x2:W-PAD_R, y2:y(v),
      stroke:'#2a3556', 'stroke-width':0.5}));
    const t = svgEl('text', {x:5, y:y(v)+3, fill:'#6b7899', 'font-size':10});
    t.textContent = v+' MHz'; svg.appendChild(t);
  }
  [0, Math.floor(ts.length/2), ts.length-1].forEach(i=>{
    const d = new Date(ts[i]*1000);
    const label = String(d.getHours()).padStart(2,'0') + ':' + String(d.getMinutes()).padStart(2,'0');
    const t = svgEl('text', {x:x(ts[i]), y:H-5, fill:'#6b7899', 'font-size':10, 'text-anchor':'middle'});
    t.textContent = label; svg.appendChild(t);
  });
  let dObs = '';
  ts.forEach((tv,i)=>{
    const v = history.fxes_obs[i];
    if(v === null || v === undefined) return;
    dObs += (dObs ? 'L ' : 'M ') + x(tv).toFixed(1) + ',' + y(v).toFixed(1) + ' ';
  });
  if(dObs) svg.appendChild(svgEl('path', {d:dObs, fill:'none', stroke:'#5fb3ff', 'stroke-width':2}));
  let dPred = '';
  ts.forEach((tv,i)=>{
    const v = history.fx30_pred[i];
    if(v === null || v === undefined) return;
    dPred += (dPred ? 'L ' : 'M ') + x(tv).toFixed(1) + ',' + y(v).toFixed(1) + ' ';
  });
  if(dPred) svg.appendChild(svgEl('path', {d:dPred, fill:'none', stroke:'#ffb020',
    'stroke-width':2, 'stroke-dasharray':'6,4'}));
  const lg = svgEl('g');
  lg.appendChild(svgEl('line', {x1:W-180, y1:15, x2:W-155, y2:15, stroke:'#5fb3ff', 'stroke-width':2}));
  const t1 = svgEl('text', {x:W-150, y:19, fill:'#c9d5f0', 'font-size':11});
  t1.textContent = '実測 fxEs'; lg.appendChild(t1);
  lg.appendChild(svgEl('line', {x1:W-180, y1:30, x2:W-155, y2:30, stroke:'#ffb020',
    'stroke-width':2, 'stroke-dasharray':'6,4'}));
  const t2 = svgEl('text', {x:W-150, y:34, fill:'#c9d5f0', 'font-size':11});
  t2.textContent = '30分前に立てた予測'; lg.appendChild(t2);
  svg.appendChild(lg);

  let errs = [];
  ts.forEach((tv,i)=>{
    if(history.fxes_obs[i] !== null && history.fx30_pred[i] !== null
       && history.fxes_obs[i] !== undefined && history.fx30_pred[i] !== undefined){
      errs.push(Math.abs(history.fxes_obs[i] - history.fx30_pred[i]));
    }
  });
  const wrap = el('div', 'chart-wrap');
  wrap.appendChild(svg);
  if(errs.length > 0){
    const avgErr = errs.reduce((a,b)=>a+b,0) / errs.length;
    const info = el('div', 'sub small');
    info.style.marginTop = '4px';
    info.textContent = `検証: 30分予測 vs 実測 の平均誤差 = ${avgErr.toFixed(2)} MHz (${errs.length}件)`;
    wrap.appendChild(info);
  } else {
    wrap.appendChild(el('div', 'sub small', '検証データ蓄積中'));
  }
  return wrap;
}

// ==================== 学習品質 ====================
function renderLQ(lq){
  const c = el('div', 'card');
  c.appendChild(el('h2', null, '📊 学習の品質と成熟度'));
  const head = el('div');
  head.style.display='flex'; head.style.alignItems='baseline';
  head.style.gap='10px'; head.style.flexWrap='wrap';
  head.appendChild(el('div', 'pct-big ' + lq.overall_class, lq.overall_pct + '%'));
  head.appendChild(el('div', lq.overall_class, lq.overall_label));
  head.appendChild(el('div', 'sub small', '3段階マイルストーン平均'));
  c.appendChild(head);
  (lq.stages||[]).forEach(s=>{
    const row = el('div', 'lq-row');
    const h = el('div', 'lq-head');
    const n = el('div', 'name');
    n.appendChild(document.createTextNode(s.name + ' '));
    n.appendChild(el('span', s.done ? 'tag-good' : 'tag-warn', s.done ? '達成' : '未達'));
    h.appendChild(n);
    h.appendChild(el('div', 'cur', s.current + ' / ' + s.target + ' (' + s.pct + '%)'));
    row.appendChild(h);
    row.appendChild(el('div', 'lq-desc', s.desc));
    row.appendChild(pctBar(s.pct, !s.done));
    if(!s.done){
      row.appendChild(el('div', 'lq-eta state-warn',
        '⏳ ' + s.eta + (s.eta_note ? ' ' + s.eta_note : '')));
    } else {
      row.appendChild(el('div', 'lq-eta state-good', '✅ この段階の予測モデル稼働中'));
    }
    c.appendChild(row);
  });
  return c;
}

// Ver13.5: 各要素の表示境界リスト
function renderReadiness(items){
  if(!items || !items.length) return null;
  const c = el('div', 'card');
  c.appendChild(el('h2', null, '📋 各要素の表示可否 (サイクル境界)'));
  c.appendChild(el('div', 'sub',
    '各予測要素がいつから表示可能かをリスト表示。 ✓ = 表示中、それ以外は待機中の残サイクル数。'));
  const list = el('div', 'readiness-list');
  items.forEach(it=>{
    const row = el('div', 'readiness-row ' + (it.ready ? 'ready' : 'waiting'));
    const lbl = el('div', 'lbl');
    lbl.textContent = it.label;
    if(it.note){
      const note = el('span', 'readiness-note', ' (' + it.note + ')');
      lbl.appendChild(note);
    }
    row.appendChild(lbl);
    row.appendChild(el('div', 'st', it.status));
    list.appendChild(row);
  });
  c.appendChild(list);
  return c;
}

// ==================== メイン描画 ====================
// Rev16.2(②対応): 従来はこの関数がEs/F2を1つの表に統合表示しており、
// 「くわしい表示」を選んでも常にEs/F2が混在した表のままだった。
// ここではEs(スポラディックE層)専用・数値詳細の表とし(「くわしい表示」限定)、
// F2層は別カード(下部の🌐F2層 伝搬状態推定カード)で完全に独立して数値表示する。
function renderEsAreaTableAdvanced(){
  const container = document.getElementById('area-table-container');
  if(!container || !latestSnapshot) return;
  const areasNow = (latestSnapshot.map && latestSnapshot.map.areas_by_horizon
    && latestSnapshot.map.areas_by_horizon[String(currentHorizon)]) || latestSnapshot.areas || [];
  container.innerHTML = '';
  if(!areasNow.length){
    container.appendChild(el('div', 'sub', 'エリアデータ収集中'));
    return;
  }
  const tbl = el('table');
  const th = el('tr');
  ['方面','Es +'+currentHorizon+'分後','Es現在','傾向','FT8受信'].forEach(h=>{
    th.appendChild(el('th', null, h));
  });
  tbl.appendChild(th);
  areasNow.forEach(a=>{
    const tr = el('tr');
    tr.appendChild(el('td', null, a.display_name));
    const td1 = el('td');
    const wr = el('div');
    wr.style.display='flex'; wr.style.gap='6px'; wr.style.alignItems='center';
    const b = pctBar(a.display_pct||0, (a.display_pct||0)<30);
    b.style.flex='1'; b.style.minWidth='60px';
    wr.appendChild(b);
    wr.appendChild(el('span', a.p_class,
      (a.display_pct === null || a.display_pct === undefined)
        ? '--' : Math.round(a.display_pct)+'%'));
    td1.appendChild(wr);
    tr.appendChild(td1);
    tr.appendChild(el('td', 'mono', a.current_pct + '%'));
    tr.appendChild(el('td', a.trend_class, a.trend_mark + ' ' + a.trend_label));
    tr.appendChild(el('td', 'mono sub', String(a.heard_count)));
    tbl.appendChild(tr);
  });
  container.appendChild(tbl);
  container.appendChild(el('div', 'sub small',
    'Es(スポラディックE層)のみの方面別到達率(数値)です。F2層(遠距離)の数値詳細は'
    + '下部の「🌐 F2層 伝搬状態推定(開通期待)」カードで独立してご覧いただけます。'));
}
// Rev16.2(②③対応): 「かんたん表示」専用。Es/F2それぞれの方面別予測を、
// 数値ではなく赤/オレンジ/青のトラフィックライト信号で一覧表示する。
// タイムスライダー操作に連動して再描画される(renderEsAreaTableAdvanced()と同じ
// areas_by_horizonデータソースを使うため、選択中ホライズンの内容が反映される)。
function renderAreaSignalTable(){
  const container = document.getElementById('area-signal-container');
  if(!container || !latestSnapshot) return;
  const areasNow = (latestSnapshot.map && latestSnapshot.map.areas_by_horizon
    && latestSnapshot.map.areas_by_horizon[String(currentHorizon)]) || latestSnapshot.areas || [];
  container.innerHTML = '';
  if(!areasNow.length){
    container.appendChild(el('div', 'sub', 'エリアデータ収集中'));
    return;
  }
  const tbl = el('table');
  const th = el('tr');
  ['方面', 'Es(スポラディックE層)', 'F2層(遠距離)'].forEach(h=>{
    th.appendChild(el('th', null, h));
  });
  tbl.appendChild(th);
  areasNow.forEach(a=>{
    const tr = el('tr');
    tr.appendChild(el('td', null, a.display_name));
    const tlEs = trafficLight(a.display_pct);
    const tdEs = el('td');
    tdEs.appendChild(document.createTextNode(tlEs.emoji + ' '));
    tdEs.appendChild(el('span', null, tlEs.label));
    tr.appendChild(tdEs);
    const tdF2 = el('td');
    if(a.f2_pct === null || a.f2_pct === undefined){
      tdF2.appendChild(el('span', 'sub', '対象外'));
    } else {
      const tlF2 = trafficLight(a.f2_pct);
      tdF2.appendChild(document.createTextNode(tlF2.emoji + ' '));
      tdF2.appendChild(el('span', null, tlF2.label));
    }
    tr.appendChild(tdF2);
    tbl.appendChild(tr);
  });
  container.appendChild(tbl);
  container.appendChild(el('div', 'sub small',
    '🔴=期待できます　🟠=兆候があります　🔵=今は静かです　⚪=データ不足。 '
    + '「対象外」は現在地から近すぎる、またはF2遠隔判定対象(最大5方面)に選ばれていない方面です。'));
}
function render(data){
  const root = document.getElementById('root');
  root.innerHTML = '';
  window._snapArrivedAt = Date.now();

  // ★ シーズンオフ警告バナー
  const season = data.season || {};
  if(season.is_season_off){
    const banner = el('div', 'season-off-warn');
    const pctTxt = (season.blended_pct!=null) ? season.blended_pct.toFixed(0)+'%' : '--%';
    const recentTxt = (season.recent_pct!=null) ? season.recent_pct.toFixed(0)+'%' : '--%';
    const climTxt = (season.climatology_pct!=null) ? season.climatology_pct.toFixed(0)+'%' : '--%';
    banner.innerHTML = '⚠️ <b>Es オフシーズン判定中</b>: 直近72時間の活動率('+recentTxt
      +')と当月の気候学的ベースライン('+climTxt+')をブレンドした活動率が'+pctTxt
      +'です。Es 予測は参考値としてご覧ください。';
    root.appendChild(banner);
  }

  // Rev15.0: PSKReporterコールサイン未設定の警告バナー
  // ESDUCT(整理版・バグ修正): 従来は「起動オプション --psk-callsign で指定して
  // ください」という、GUI利用者には実行できない案内のみだった。右上の
  // 「⚙️ 設定」ボタンを直接開けるリンクを添え、その場で設定できるようにする。
  const ft8i = data.ft8 || {};
  if(ft8i.no_callsign){
    const cbanner = el('div', 'season-off-warn');
    cbanner.innerHTML = '⚠️ <b>PSKReporterコールサイン未設定</b>: FT8連携は動作していません。'
      + ' 伝播方位の異方性推定・学習ステージ3(Motion Predictor)には到達できません。'
      + ' <a href="#" id="ft8-banner-settings-link" style="color:var(--accent);text-decoration:underline">画面右上の「⚙️ 設定」から指定してください</a>。';
    root.appendChild(cbanner);
    const link = cbanner.querySelector('#ft8-banner-settings-link');
    if(link) link.addEventListener('click', (e)=>{ e.preventDefault(); openSettingsModal(false); });
  }

  // Rev15.9: 初心者向け要約カード(常時表示、モードに関わらず先頭に出す)
  root.appendChild(renderBeginnerSummary(data));

  // ① Hero: 予測ファースト + カウントダウン刷新
  const hero = el('div', 'card hero-forecast');
  const hh = el('div', 'hero-header');
  const hL = el('div', 'grow');
  hL.appendChild(el('h1', null, '🎯 30分後の予測 (自動更新)'));
  const seasonBadge = el('span', 'badge ' + (season.class || 'low'));
  seasonBadge.textContent = 'Es: ' + (season.label || '判定中');
  hL.appendChild(el('div', 'sub', '本ソフトは "Es 発生の予測" が目的。'));
  hL.appendChild(seasonBadge);
  hh.appendChild(hL);

  // カウントダウン (Ver13.5: パルスドット + "動作中" 表示)
  const cd = el('div', 'countdown-box');
  cd.appendChild(el('span', 'pulse-dot'));
  cd.appendChild(el('span', 'sub', '次サイクル'));
  const cdVal = el('span', 'val'); cdVal.id = 'cdv'; cdVal.textContent = '計算中';
  cd.appendChild(cdVal);
  cd.appendChild(el('span', 'status-alive', '⏱ 動作中'));
  hh.appendChild(cd);
  hero.appendChild(hh);
  const cdBar = el('div', 'progress'); cdBar.style.marginTop='4px';
  const cdSpan = el('span'); cdSpan.id = 'cdb'; cdSpan.style.width='0%';
  cdBar.appendChild(cdSpan); hero.appendChild(cdBar);

  const now = data.now || {};
  const fc = data.forecast_30m || {};
  const metric = el('div', 'hero-metric');
  const b1 = el('div');
  b1.appendChild(el('div', 'sub', '代表局 (' + (now.best_station||'--') + ') の 30分後 fxEs'));
  const bV = el('div', 'pct-big ' + (now.fxEs_class||'state-unknown'));
  bV.textContent = (fc.fx30 === null || fc.fx30 === undefined) ? '学習中'
    : fmtNum(fc.fx30,1) + ' MHz';
  b1.appendChild(bV);
  b1.appendChild(el('div', 'sub', '傾向: ' + (fc.trend_mark||'…') + ' ' + (fc.trend_label||'')));
  metric.appendChild(b1);
  const b2 = el('div');
  b2.appendChild(el('div', 'sub', '現在の代表局 開通率'));
  const bP = el('div', 'pct-big ' + (now.p27_class||'state-unknown'));
  bP.textContent = (now.p27 === null) ? '--' : Math.round(now.p27) + '%';
  b2.appendChild(bP);
  b2.appendChild(el('div', 'sub', now.p27_label || ''));
  metric.appendChild(b2);
  hero.appendChild(metric);
  root.appendChild(hero);

  // ② 日本地図 (Ver13.5 精緻化) + タイムスライダ
  const mapCard = el('div', 'card');
  const mapTitle = el('h2', null, '🗾 予測地図 (47都道府県)');
  mapCard.appendChild(mapTitle);
  const mapProv = (data.map && data.map.es_centroid && data.map.es_centroid.provisional);
  if(mapProv){
    const tag = el('span', 'provisional-tag', '暫定表示');
    mapTitle.appendChild(tag);
  }

  const sliderBox = el('div', 'horizon-slider');
  sliderBox.appendChild(el('span', 'sub', '予測時間軸:'));
  const slider = el('input');
  slider.type = 'range';
  const horizons = [0, 15, 30, 45, 60, 90];
  slider.min = 0; slider.max = horizons.length - 1; slider.step = 1;
  slider.value = horizons.indexOf(currentHorizon);
  const val = el('span', 'value', (currentHorizon === 0 ? '現在' : '+' + currentHorizon + '分'));
  // Rev15.8③: 実時刻(JST)の小さな自動表記。 .horizon-time-lbl クラスを持つ
  // 要素は複数箇所(スライダー横・エリア表見出し)に置け、updateHorizonTimeLabel()
  // が一括更新する。
  const timeLbl = el('span', 'sub small horizon-time-lbl', '');
  const updateMapForHorizon = ()=>{
    // Rev16.2(③対応): 従来は地図(drawMap)とエリア表(renderAreaTable)のみが
    // スライダーに連動し、「かんたん予報」カードの内容は+0分に固定されたまま
    // だった。ここで初心者/上級者の両表示をまとめて更新し、選択中ホライズンの
    // 内容が地図以外にも確実に反映されるようにする。
    val.textContent = (currentHorizon === 0 ? '現在' : '+' + currentHorizon + '分');
    const sub = document.getElementById('map-horizon-sub');
    if(sub) sub.textContent = mapSubText(currentHorizon);
    const wrap = document.getElementById('mapwrap');
    if(wrap && latestSnapshot){
      wrap.innerHTML = '';
      wrap.appendChild(drawMap(latestSnapshot.map || {}, currentHorizon));
      bindTips(wrap);
    }
    const areaTitle = document.getElementById('area-table-title');
    if(areaTitle) areaTitle.textContent = '📍 Es 方面別 +' + currentHorizon + '分後 予測 (数値詳細)';
    renderEsAreaTableAdvanced();
    renderAreaSignalTable();
    if(latestSnapshot) refreshBeginnerSummaryBody(latestSnapshot);
    updateHorizonTimeLabel();
  };
  slider.oninput = ()=>{
    stopHorizonPlayback();
    currentHorizon = horizons[parseInt(slider.value)];
    updateMapForHorizon();
  };
  sliderBox.appendChild(slider);
  sliderBox.appendChild(val);
  sliderBox.appendChild(timeLbl);

  // Ver13.14: 自動送り再生ボタン (0→90分を順に表示し、Es の移動方向を視認しやすくする)
  const playBtn = el('button', 'btn ghost horizon-play-btn', '▶ 自動再生');
  playBtn.type = 'button';
  playBtn.onclick = ()=>{
    if(window._horizonPlayTimer){
      stopHorizonPlayback();
      return;
    }
    playBtn.textContent = '⏸ 停止';
    let idx = horizons.indexOf(currentHorizon);
    window._horizonPlayTimer = setInterval(()=>{
      idx = (idx + 1) % horizons.length;
      currentHorizon = horizons[idx];
      slider.value = idx;
      updateMapForHorizon();
      if(idx === horizons.length - 1){
        // 90分まで到達したら1回停止して分かりやすくする
        clearInterval(window._horizonPlayTimer);
        window._horizonPlayTimer = null;
        playBtn.textContent = '▶ 自動再生';
      }
    }, 1200);
  };
  sliderBox.appendChild(playBtn);
  mapCard.appendChild(sliderBox);
  const mapSub = el('div', 'sub'); mapSub.id = 'map-horizon-sub';
  mapSub.textContent = mapSubText(currentHorizon);
  mapCard.appendChild(mapSub);
  const mapWrap = el('div', 'map-wrap paper'); mapWrap.id = 'mapwrap';
  mapWrap.appendChild(drawMap(data.map || {}, currentHorizon));
  bindTips(mapWrap);
  mapCard.appendChild(mapWrap);
  const legend = el('div', 'map-legend');
  [['#ff5252','75%以上 強開通'],['#26c281','50-75% 開通'],
   ['#ffb020','30-50% 注意'],['#8b7cff','15-30% 弱兆'],
   ['#4b587a','15%未満 低調']].forEach(([c,t])=>{
    const s = el('span');
    const sw = el('span', 'sw'); sw.style.background = c;
    s.appendChild(sw); s.appendChild(document.createTextNode(' '+t));
    legend.appendChild(s);
  });
  mapCard.appendChild(legend);

  // Ver13.21: 実測/推定の凡例 + 非相反性リスクバッジ
  const legend2 = el('div', 'map-legend');
  [['line', '#26c281', '実測伝播 (自局→相手、PSKReporter TXレポートのみ・片方向)'],
   ['dash', '#ff5252', '推定Es帯 (実測から推定、雲の実形状ではない)']].forEach(([kind,c,t])=>{
    const s = el('span');
    const sw = el('span', 'sw'); sw.style.background = c;
    if(kind === 'dash'){ sw.style.opacity = '0.5'; }
    s.appendChild(sw); s.appendChild(document.createTextNode(' '+t));
    legend2.appendChild(s);
  });
  mapCard.appendChild(legend2);

  const rr = (data.map && data.map.reciprocity_risk) || null;
  if(rr){
    const rrBadge = el('div', 'sub');
    rrBadge.style.marginTop = '6px';
    rrBadge.style.fontWeight = 'bold';
    const rrIcon = rr.level === '高' ? '⚠ ' : (rr.level === '不明' ? '？ ' : '');
    rrBadge.textContent = rrIcon + '非相反性リスク: ' + rr.level + '　(' + (rr.note || '') + ')';
    rrBadge.className = 'sub ' + (rr.class || '');
    mapCard.appendChild(rrBadge);
  }

  root.appendChild(mapCard);

  // ③ 予測 vs 実測 履歴グラフ (くわしい表示のみ)
  const histCard = el('div', 'card tier-advanced');
  histCard.appendChild(el('h2', null, '📈 予測 vs 実測 履歴 (直近24時間)'));
  histCard.appendChild(el('div', 'sub',
    '青実線: 実測 fxEs。 橙破線: その時刻の 30分前に立てた予測。'));
  histCard.appendChild(drawHistoryChart(data.history || {}));
  root.appendChild(histCard);

  // ④ エリア別詳細
  // Rev16.2(②③対応): 従来はEs/F2を1つの表に統合し、初心者/上級者どちらの
  // 表示でも同じ内容だった。ここでは表示モードに応じて完全に別カードへ分離する:
  //   ・かんたん表示: 赤/オレンジ/青のトラフィックライト信号で方面別の目安を表示
  //     (renderAreaSignalTable())。Es/F2は同じ表内だが視覚的には信号色のみ。
  //   ・くわしい表示: Es専用の数値詳細表(renderEsAreaTableAdvanced())。
  //     F2の数値詳細は下部の別カード(🌐F2層 伝搬状態推定)に完全分離。
  // いずれもタイムスライダー操作(updateMapForHorizon())に連動して再描画される。
  const areaSignalCard = el('div', 'card tier-beginner-only');
  areaSignalCard.appendChild(el('h2', null, '📍 方面別 開通の目安'));
  const areaSignalTimeLbl = el('div', 'sub small');
  areaSignalTimeLbl.appendChild(document.createTextNode(
    (currentHorizon === 0 ? '現在' : '+' + currentHorizon + '分後') + ' '));
  areaSignalTimeLbl.appendChild(el('span', 'horizon-time-lbl'));
  areaSignalCard.appendChild(areaSignalTimeLbl);
  areaSignalCard.appendChild(el('div', 'sub small',
    '色の信号で方面ごとの開通しやすさを一目で確認できます'
    + '（🔴期待できます／🟠兆候あり／🔵静かです）。'
    + '上部地図のスライダーで時間を変えると、この表も連動して変わります。'));
  const areaSignalContainer = el('div'); areaSignalContainer.id = 'area-signal-container';
  areaSignalCard.appendChild(areaSignalContainer);
  root.appendChild(areaSignalCard);
  renderAreaSignalTable();

  const areaCard = el('div', 'card tier-advanced');
  const areaTitle = el('h2', null, '📍 Es 方面別 +' + currentHorizon + '分後 予測 (数値詳細)');
  areaTitle.id = 'area-table-title';
  areaCard.appendChild(areaTitle);
  const areaTimeLbl = el('span', 'sub small horizon-time-lbl', '');
  areaCard.appendChild(areaTimeLbl);
  areaCard.appendChild(el('div', 'sub small',
    'Es(スポラディックE層)のみの方面別到達率(数値)です。'
    + 'F2層(遠距離)の数値詳細は下部の「🌐 F2層 伝搬状態推定(開通期待)」カードで'
    + '独立してご覧いただけます(Es/F2は完全に別カードです)。'));
  const areaTableContainer = el('div'); areaTableContainer.id = 'area-table-container';
  areaCard.appendChild(areaTableContainer);
  root.appendChild(areaCard);
  renderEsAreaTableAdvanced();
  updateHorizonTimeLabel();

  // ④-2 Rev16.0/16.1: F2層 伝搬状態推定(開通期待) — 詳細情報は引き続き
  // 別カードで表示する(上表④はこのカードの数値を方面別に要約したもの)。
  const f2fc = data.f2_forecast || [];
  if(f2fc.length){
    const cf2 = el('div', 'card tier-advanced');
    cf2.appendChild(el('h2', null, '🌐 F2層 伝搬状態推定(開通期待) [現在地からの遠隔方面]'));
    cf2.appendChild(el('div', 'sub',
      'Es予測とは独立した別評価系です。NICT実測foF2/hmF2/M(3000)F2から、'
      + '球面幾何(Earth curvature)を考慮したMUFモデルと実測M(3000)F2ベースの'
      + 'MUFモデルの両方を独立に評価し、整合性を確認した上で統合評価値を算出します。'
      + 'LUF(D層吸収限界)はSFI/Ap実測値と経路上の太陽天頂角による簡易吸収指標'
      + '(VOACAP/ITU-R本体の物理演算ではない)です。'
      + '開通期待度(f2_expectancy)は統計的確率ではなく、MUF/LUFマージンを正規化した'
      + '運用指標であり、モデル信頼度(f2_confidence)とは別軸のスコアです。'
      + '⬆ 上表「方面別予測」のF2列はこのカードの現在(+0分)〜選択中ホライズンの値を要約したものです。'));
    const consistencyJp = {
      'CONSISTENT': '両モデル整合', 'CAUTION': '要注意(やや乖離)',
      'DIVERGENT': '不一致大(保守値採用)', 'GEOMETRY_ONLY': '幾何モデルのみ(M(3000)F2欠測)',
      'UNKNOWN': '評価不能',
    };
    f2fc.forEach(item=>{
      const box = el('div');
      box.style.marginTop='10px'; box.style.paddingTop='10px';
      box.style.borderTop='1px solid var(--line)';
      const head = el('div');
      head.style.display='flex'; head.style.justifyContent='space-between'; head.style.flexWrap='wrap';
      head.appendChild(el('div', null, '📡 ' + item.area_name));
      if(item.no_data){
        head.appendChild(el('span', 'state-warn', item.reason || 'データ不足'));
        box.appendChild(head);
        cf2.appendChild(box);
        return;
      }
      head.appendChild(el('span', 'sub mono',
        '距離約' + Math.round(item.distance_km) + 'km (参照局:' + (item.source_station || '--') + ')'));
      box.appendChild(head);
      let ageStr = '--';
      const stAges = item.station_ages || {};
      const stKeys = Object.keys(stAges);
      if(stKeys.length){
        ageStr = stKeys.map(st=>{
          const age = stAges[st];
          return st + (age===null||age===undefined ? '--' : Math.round(age)+'分前');
        }).join(' / ');
      }
      box.appendChild(el('div', 'sub small',
        'foF2=' + item.foF2 + 'MHz  hmF2=' + (item.hmF2!==null&&item.hmF2!==undefined ? item.hmF2+'km' : '--')
        + '  M(3000)F2=' + (item.M3000F2!==null&&item.M3000F2!==undefined ? item.M3000F2 : '--')
        + '  (観測鮮度: ' + ageStr + ')'));
      box.appendChild(el('div', 'sub small',
        '有効局数=' + (item.valid_station_count!==undefined ? item.valid_station_count : '--')
        + '局 / 陳腐化局数=' + (item.stale_station_count!==undefined ? item.stale_station_count : '--')
        + '局  ホップ数=' + (item.hop_count!==undefined ? item.hop_count : '--')
        + '(単一ホップ最大約' + (item.max_single_hop_km!==undefined&&item.max_single_hop_km!==null ? Math.round(item.max_single_hop_km) : '--') + 'km)'));
      const mufBox = el('div', 'sub small');
      mufBox.style.marginTop = '4px';
      const consLabel = consistencyJp[item.muf_consistency] || item.muf_consistency || '--';
      mufBox.appendChild(document.createTextNode(
        'MUF診断: 幾何=' + (item.muf_geometry!==undefined&&item.muf_geometry!==null?item.muf_geometry+'MHz':'--')
        + ' / M(3000)F2実測=' + (item.muf_3000_reference!==undefined&&item.muf_3000_reference!==null?item.muf_3000_reference+'MHz':'--')
        + ' / 統合評価=' + (item.muf_estimate!==undefined&&item.muf_estimate!==null?item.muf_estimate+'MHz':'--')));
      const consBadgeCls = item.muf_consistency==='DIVERGENT' ? 'state-warn'
        : (item.muf_consistency==='CAUTION' ? 'state-warn' : 'state-good');
      const consSpan = el('span', consBadgeCls, ' [' + consLabel + ']');
      consSpan.style.marginLeft = '4px';
      mufBox.appendChild(consSpan);
      box.appendChild(mufBox);
      if((item.muf_nict_oblique!==undefined && item.muf_nict_oblique!==null)
         || (item.luf_nict_oblique!==undefined && item.luf_nict_oblique!==null)){
        box.appendChild(el('div', 'sub small',
          'NICT実測クロスチェック(斜め伝搬可能周波数): MUF='
          + (item.muf_nict_oblique!==null&&item.muf_nict_oblique!==undefined?item.muf_nict_oblique+'MHz':'--')
          + ' / LUF=' + (item.luf_nict_oblique!==null&&item.luf_nict_oblique!==undefined?item.luf_nict_oblique+'MHz':'--')
          + '（参考表示、主計算には不使用）'));
      }
      const confPct = item.f2_confidence;
      if(confPct !== undefined && confPct !== null){
        const confRow = el('div', 'sub small');
        confRow.style.marginTop = '2px';
        confRow.appendChild(document.createTextNode('モデル信頼度(f2_confidence): '));
        const confBar = pctBar(confPct, confPct<50);
        confBar.style.width = '80px'; confBar.style.display='inline-block'; confBar.style.verticalAlign='middle';
        confRow.appendChild(confBar);
        confRow.appendChild(document.createTextNode(' ' + Math.round(confPct) + '%'));
        box.appendChild(confRow);
      }
      const tbl2 = el('table');
      const th2 = el('tr');
      (item.horizons||[]).forEach(h=>{
        const timeS = h.forecast_time_jst ? '(' + h.forecast_time_jst + ')' : '';
        th2.appendChild(el('th', null, '+' + h.minutes + '分' + timeS));
      });
      tbl2.appendChild(th2);
      const tr2 = el('tr');
      (item.horizons||[]).forEach(h=>{
        const td = el('td');
        const wr = el('div');
        wr.style.display='flex'; wr.style.gap='4px'; wr.style.alignItems='center'; wr.style.flexDirection='column';
        const pct = (h.f2_expectancy!==undefined ? h.f2_expectancy : (h.expectancy_pct!==undefined ? h.expectancy_pct : h.probability_pct));
        const b = pctBar(pct||0, (pct||0)<35);
        b.style.width='60px';
        wr.appendChild(b);
        const margin = h.margin_mhz!==undefined && h.margin_mhz!==null ? ' M'+(h.margin_mhz>=0?'+':'')+Number(h.margin_mhz).toFixed(1) : '';
        const marginLuf = h.margin_luf_mhz!==undefined && h.margin_luf_mhz!==null ? ' L'+(h.margin_luf_mhz>=0?'+':'')+Number(h.margin_luf_mhz).toFixed(1) : '';
        const role = h.horizon_role ? ' '+h.horizon_role : '';
        wr.appendChild(el('span', 'mono small',
          (pct===null||pct===undefined ? '--' : Math.round(pct)+'%') + ' ' + h.label + margin + marginLuf + role));
        if(h.muf!==undefined && h.muf!==null){
          wr.appendChild(el('span', 'sub small',
            'MUF' + h.muf + ' / LUF' + (h.luf!==undefined&&h.luf!==null?h.luf:'--')));
        }
        td.appendChild(wr);
        tr2.appendChild(td);
      });
      tbl2.appendChild(tr2);
      box.appendChild(tbl2);
      cf2.appendChild(box);
    });
    cf2.appendChild(el('div', 'sub small',
      'M=MUFマージン、L=LUFマージン(いずれも27MHz帯からの余裕[MHz])。'
      + 'MUFは幾何モデル/M(3000)F2実測モデルを統合した評価値(muf_estimate)を使用。'
      + 'LUF(D層吸収限界)は簡易吸収指標であり、VOACAP/ITU-R本体の物理演算ではありません。'
      + '+45/+60分は参考(外挿)、+90分はトレンド外挿のため予測として強く扱いません。'));
    root.appendChild(cf2);
  }

  // ⑤ Ver13.5: Es 形状/移動 情報 (くわしい表示のみ)
  const mot = data.motion || {};
  const c3 = el('div', 'card tier-advanced');
  const shapeTitle = el('h2', null, '↔️ Es 電離雲の形状と移動');
  c3.appendChild(shapeTitle);
  if(mot.provisional){
    shapeTitle.appendChild(el('span', 'provisional-tag', '暫定'));
  }
  const kpi = el('div', 'kpi');
  function kpiItem(l, v, h){
    const d = el('div');
    d.appendChild(el('div', 'label', l));
    d.appendChild(el('div', 'value', v));
    if(h) d.appendChild(el('div', 'sub small', h));
    return d;
  }
  kpi.appendChild(kpiItem('重心 移動方向', mot.direction || '不明',
    'NICT 4局のfxEs加重重心緯度変化'));
  kpi.appendChild(kpiItem('重心 速度',
    mot.speed_kmh !== null ? Math.round(mot.speed_kmh)+' km/h' : '--'));
  kpi.appendChild(kpiItem('MSTID 方向',
    (mot.mstid_dir_deg !== null && mot.mstid_dir_deg !== undefined)
      ? Math.round(mot.mstid_dir_deg)+'°' : '--',
    '衛星画像 (北=0° 東=90°)'));
  kpi.appendChild(kpiItem('MSTID 速度',
    (mot.mstid_speed_kmh !== null && mot.mstid_speed_kmh !== undefined)
      ? Math.round(mot.mstid_speed_kmh)+' km/h' : '--'));
  const ft8 = data.ft8 || {};
  kpi.appendChild(kpiItem('Es フェーズ', ft8.es_phase_repr || '---', ft8.phase_hint));
  c3.appendChild(kpi);
  root.appendChild(c3);

  // ⑤b Ver13.20: Es伝播方位の異方性(粗い指標) — ローズダイアグラム
  if(ft8.anisotropy_rose_svg){
    const c3b = el('div', 'card tier-advanced');
    c3b.appendChild(el('h2', null, '🌐 Es伝播方位の異方性 (粗い指標)'));
    c3b.appendChild(el('div', 'sub',
      'FT8伝播レポートの方位分布。花びらの偏りが等方/異方性を、' +
      '赤(強)・橙(弱)の軸線が帯の推定方向と南北基準線からの傾きを示す。' +
      '雲の形そのものを再構成するものではなく、あくまで簡易判定。'));
    const roseWrap = el('div');
    roseWrap.style.display = 'flex'; roseWrap.style.gap = '14px'; roseWrap.style.alignItems = 'center';
    roseWrap.style.flexWrap = 'wrap';
    const roseImg = el('div'); roseImg.innerHTML = ft8.anisotropy_rose_svg;
    roseWrap.appendChild(roseImg);
    const roseKpi = el('div', 'kpi');
    roseKpi.style.flex = '1'; roseKpi.style.minWidth = '160px';
    roseKpi.appendChild(kpiItem('判定', ft8.anisotropy_label || 'データ不足'));
    roseKpi.appendChild(kpiItem('R値 (0=等方 / 1=強い帯状)',
      (ft8.anisotropy_r === null || ft8.anisotropy_r === undefined) ? '--' : ft8.anisotropy_r.toFixed(2)));
    roseKpi.appendChild(kpiItem('帯の向き',
      (ft8.anisotropy_orientation_deg === null || ft8.anisotropy_orientation_deg === undefined)
        ? '--' : Math.round(ft8.anisotropy_orientation_deg) + '° (0=南北,90=東西)'));
    roseKpi.appendChild(kpiItem('南北整合度',
      (ft8.anisotropy_ns_alignment === null || ft8.anisotropy_ns_alignment === undefined)
        ? '--' : (ft8.anisotropy_ns_alignment >= 0 ? '+' : '') + ft8.anisotropy_ns_alignment.toFixed(2)));
    roseKpi.appendChild(kpiItem('観測信頼度 (Observation)',
      ft8.observation_n + '局 [' + (ft8.observation_label || '--') + ']',
      'RPI/異方性判定の元になった観測局数。値そのものと観測量は別物として見る'));
    roseWrap.appendChild(roseKpi);
    c3b.appendChild(roseWrap);
    root.appendChild(c3b);
  }

  // ⑥ Ver13.5: 各要素の表示境界リスト (くわしい表示のみ)
  const readiness = renderReadiness(data.element_readiness);
  if(readiness){ readiness.classList.add('tier-advanced'); root.appendChild(readiness); }

  // ⑦ Ver13.5: 学習品質 (くわしい表示のみ)
  const lqCard = renderLQ(data.learning_quality || {stages:[]});
  lqCard.classList.add('tier-advanced');
  root.appendChild(lqCard);

  // ⑧ NICT 4局詳細 (くわしい表示のみ)
  if(data.stations && data.stations.length){
    const c = el('div', 'card tier-advanced');
    c.appendChild(el('h2', null, '🛰 NICT 4局 数値詳細'));
    const tbl = el('table');
    const th = el('tr');
    ['観測局','fxEs 今','fxEs 30分後','開通%','傾向','形状'].forEach(h=>{
      th.appendChild(el('th', null, h));
    });
    tbl.appendChild(th);
    data.stations.forEach(s=>{
      const tr = el('tr');
      tr.appendChild(el('td', null, s.name));
      tr.appendChild(el('td', 'mono', s.fxEs === null ? '--' : fmtNum(s.fxEs,1)));
      tr.appendChild(el('td', 'mono', s.fx30 === null ? '--' : fmtNum(s.fx30,1)));
      tr.appendChild(el('td', 'mono ' + (s.p27_class||''),
        s.p27 === null ? '--' : Math.round(s.p27)+'%'));
      tr.appendChild(el('td', null, s.trend_label || ''));
      tr.appendChild(el('td', null, s.local_shape || ''));
      tbl.appendChild(tr);
    });
    c.appendChild(tbl);
    root.appendChild(c);
  }

  // ⑨ FT8 / 伝播モード (くわしい表示のみ)
  const c2 = el('div', 'card tier-advanced');
  c2.appendChild(el('h2', null, '📻 FT8 実測と伝播モード'));
  const kpi2 = el('div', 'kpi');
  kpi2.appendChild(kpiItem('伝播モード', ft8.mode_label || '---', ft8.mode_hint));
  kpi2.appendChild(kpiItem('反射面品質 EQI',
    ft8.eqi !== null ? fmtNum(ft8.eqi,2) : '--', ft8.eqi_hint));
  kpi2.appendChild(kpiItem('宇宙天気',
    ft8.sw_score !== null ? fmtNum(ft8.sw_score,2) : '--', ft8.sw_hint));
  kpi2.appendChild(kpiItem('予測信頼',
    ft8.prediction_confidence !== null
      ? Math.round(ft8.prediction_confidence) + '%' : '--'));
  const cons = ft8.consistency_match
    ? el('span', 'tag-good', '一致')
    : el('span', 'tag-warn', 'ズレあり');
  const cb = el('div');
  cb.appendChild(el('div', 'label', 'NICT⇔FT8 整合'));
  const cv = el('div', 'value'); cv.appendChild(cons); cb.appendChild(cv);
  kpi2.appendChild(cb);
  c2.appendChild(kpi2);
  root.appendChild(c2);

  // ⑩ 実効寄与 (くわしい表示のみ)
  if(data.eff_top && data.eff_top.length){
    const c5 = el('div', 'card tier-advanced');
    c5.appendChild(el('h2', null, '🔍 いま予測に効いている要素'));
    data.eff_top.forEach(e=>{
      const row = el('div');
      row.style.display='flex'; row.style.gap='10px';
      row.style.alignItems='center'; row.style.margin='3px 0';
      const lb = el('div', null, e.label);
      lb.style.minWidth='140px'; lb.style.fontSize='0.9rem';
      row.appendChild(lb);
      const b = pctBar(e.pct, false); b.style.flex='1';
      row.appendChild(b);
      const v = el('div', 'mono sub');
      v.textContent = Math.round(e.pct)+'%';
      v.style.minWidth='40px'; v.style.textAlign='right';
      row.appendChild(v);
      c5.appendChild(row);
    });
    root.appendChild(c5);
  }

  // ⑪ Ver13.5: GPS/保存 (一番下)
  const locCard = el('div', 'card');
  locCard.appendChild(el('h2', null, '📍 自局位置 & 保存管理'));
  const gps = el('div');
  gps.style.display='flex'; gps.style.gap='8px';
  gps.style.alignItems='center'; gps.style.flexWrap='wrap';
  if(data.location){
    const loc = data.location;
    gps.appendChild(el('div', null,
      '緯度 ' + fmtNum(loc.lat,4) + ' / 経度 ' + fmtNum(loc.lon,4)));
    gps.appendChild(el('span', 'badge', loc.source || loc.label));
  } else {
    const st = data.location_status || {};
    gps.appendChild(el('div', 'state-warn',
      '位置未取得 — ' + (st.reason || st.status || '')));
  }
  const btn = el('button', 'btn ghost', '位置を再取得');
  btn.onclick = ()=>{ lastGpsSent = 0; reqGps();
    setTimeout(()=>{
      if(window._gpsLastError === 'insecure-context'){
        secWarn.style.display = 'block';
      }
    }, 300);
  };
  gps.appendChild(btn);
  locCard.appendChild(gps);

  // Rev16.2: secure context制約によりGeolocation APIが動作しない場合の警告表示。
  // 従来はこの制約を「Chromebook特有」と決め打ちした文言・対処法(chrome://flags)
  // を案内していたが、実際には以下のように発生要因が複数あり、Chromebookでない
  // 通常のAndroid端末(Pydroid3含む)でも発生しうることが判明した:
  //   ・LAN IP経由でダッシュボードを開いている場合(Chromebook以外でも起こりうる)
  //   ・Android WebViewなど、localhost/127.0.0.1のsecure context例外を
  //     正しく実装していないブラウザ/ハンドラで開かれた場合
  // 「Chromebook等でよく発生します」という決め打ち表現、およびモバイルでは
  // 実行が難しい/紛らわしいchrome://flagsの案内を主対策から外し、
  // 環境を問わず常に確実に使える「緯度/経度を手動設定」を第一の対処法として
  // 案内するよう変更した(chrome://flagsの案内はPC版Chrome限定の補足情報として
  // 残す程度に留める)。
  const secWarn = el('div', 'sub');
  secWarn.style.display = (window.isSecureContext === false) ? 'block' : 'none';
  secWarn.style.marginTop = '6px';
  secWarn.style.color = '#ffb020';
  secWarn.innerHTML = '⚠ このページは非HTTPS(http://)かつlocalhost以外の'
    + 'アドレス、またはGeolocation APIのsecure context判定を正しく行えない'
    + 'ブラウザ/アプリで開かれているため、ブラウザのGPS機能が動作しません。'
    + '<b>下の「緯度/経度を手動設定」をご利用ください（環境によらず確実に動作します）。</b>'
    + '<br><span class="hint">(参考: PC版のGoogle Chromeをお使いの場合に限り、'
    + '<code>chrome://flags/#unsafely-treat-insecure-origin-as-secure</code>'
    + 'にこのURL(' + location.origin + ')を追加して再起動すると、'
    + 'ブラウザGPSが利用できるようになる場合があります)</span>';
  locCard.appendChild(secWarn);

  // Rev15.9: 緯度/経度の手動入力(secure context制約に依存しないGPS代替経路)。
  // /api/gps は既存のsource許可リストに'manual'を含むため、バックエンド変更なしで
  // 利用できる。
  const manBox = el('div');
  manBox.style.marginTop = '10px';
  manBox.style.paddingTop = '10px';
  manBox.style.borderTop = '1px solid var(--line)';
  manBox.appendChild(el('div', null, '📝 緯度/経度を手動設定'));
  manBox.appendChild(el('div', 'hint',
    'GPSが使えない環境(Chromebook等)向けの代替入力です。'
    + '地図アプリ等で調べた緯度・経度を入力してください。'));
  const manRow = el('div');
  manRow.style.display = 'flex'; manRow.style.gap = '8px';
  manRow.style.flexWrap = 'wrap'; manRow.style.marginTop = '6px';
  const latInput = el('input'); latInput.type = 'text'; latInput.placeholder = '緯度 例: 35.681';
  latInput.style.cssText = 'background:var(--panel2);border:1px solid var(--line);color:var(--fg);'
    + 'border-radius:8px;padding:6px 10px;font-size:0.9rem;width:140px';
  const lonInput = el('input'); lonInput.type = 'text'; lonInput.placeholder = '経度 例: 139.767';
  lonInput.style.cssText = latInput.style.cssText;
  if(data.location){
    latInput.value = fmtNum(data.location.lat, 4);
    lonInput.value = fmtNum(data.location.lon, 4);
  }
  manRow.appendChild(latInput); manRow.appendChild(lonInput);
  const manBtn = el('button', 'btn', '設定');
  const manMsg = el('span', 'sub small');
  manMsg.style.marginLeft = '8px';
  manBtn.onclick = async ()=>{
    const lat = parseFloat(latInput.value);
    const lon = parseFloat(lonInput.value);
    if(!isFinite(lat) || !isFinite(lon) || lat < -90 || lat > 90 || lon < -180 || lon > 180){
      manMsg.textContent = '緯度/経度の形式が正しくありません';
      manMsg.style.color = '#ff5252';
      return;
    }
    manBtn.disabled = true;
    manMsg.textContent = '設定中…'; manMsg.style.color = '';
    try{
      const r = await fetch('/api/gps', {method:'POST', headers:{'Content-Type':'application/json'},
        body: JSON.stringify({lat, lon, source: 'manual'})});
      const d = await r.json();
      manMsg.textContent = d.ok ? '✅ 設定しました(次回更新で反映)' : '⚠ 失敗: ' + (d.error||'');
      manMsg.style.color = d.ok ? '#26c281' : '#ff5252';
    }catch(e){
      manMsg.textContent = '⚠ 通信エラー';
      manMsg.style.color = '#ff5252';
    }
    manBtn.disabled = false;
  };
  manRow.appendChild(manBtn); manRow.appendChild(manMsg);
  manBox.appendChild(manRow);
  locCard.appendChild(manBox);

  const fBox = el('div');
  fBox.style.marginTop='10px'; fBox.style.paddingTop='10px';
  fBox.style.borderTop='1px solid var(--line)';
  fBox.appendChild(el('div', null, '💾 学習値の手動保存'));
  fBox.appendChild(el('div', 'hint',
    '通常は毎サイクル (' + Math.round((data.cycle_info?.interval_sec||900)/60)
    + '分毎) に自動保存されます。 ブラウザ閉じる直前や、Pydroid3 切替直前に'
    + '押すと、現在サイクル途中データを即座に CSV/JSON に書き出します。'));
  const fBtn = el('button', 'btn', '💾 今すぐ学習値を保存');
  fBtn.style.marginTop='6px';
  fBtn.onclick = ()=>flushLearning(fBtn);
  fBox.appendChild(fBtn);
  locCard.appendChild(fBox);

  // Rev16.2(不具合修正): 他端末学習統合UIはinitMergeUI()として#root外へ
  // 独立させたため、render()内では何も生成しない。
  // 【原因】以前はrender()の呼び出しのたびに(=tick()による20秒毎の定期全体
  // 再描画のたびに)mergeBox一式をここで作り直しており、root.innerHTML=''で
  // #root配下が丸ごと消去される仕様と衝突していた。フォルダ閲覧中や
  // ファイル選択中に20秒経過すると表示が強制的に閉じ状態へリセットされ、
  // 選択内容(mergeSelected)も消えるため、「統合ボタンを押しても動作しない」
  // ように見える不具合の直接原因になっていた。詳細はinitMergeUI()参照。


  // Rev15.9(常駐停止): EDFSのバックグラウンド動作(NICT/FT8取得サイクル)を
  // 安全に停止する。実行中のサイクルは完了させてから停止するため、押した
  // 直後は反映まで数秒〜数十秒かかることがある(サイクル途中でのデータ
  // 不整合を避けるための仕様)。かんたん表示でも常に見えるようにする
  // (万一の停止操作をモード切替なしで行えるようにするため)。
  const stopBox = el('div');
  stopBox.style.marginTop='10px'; stopBox.style.paddingTop='10px';
  stopBox.style.borderTop='1px solid var(--line)';
  stopBox.appendChild(el('div', null, '⏹ 常駐停止'));
  stopBox.appendChild(el('div', 'hint',
    'EDFSのバックグラウンド動作(NICT/FT8取得サイクル)を停止します。'
    + '停止前に現在までの学習データ(セッションログ・校正プロファイル)を'
    + '自動的に保存してから停止します。'
    + 'APK版ではアプリ内の「常駐停止」ボタンと同じ機能です。'
    + '停止後にダッシュボードを再開するには、アプリ/プログラムの再起動が必要です。'));
  const stopBtn = el('button', 'btn ghost', '⏹ 常駐停止');
  stopBtn.style.marginTop='6px'; stopBtn.style.borderColor='var(--hot)'; stopBtn.style.color='var(--hot)';
  stopBtn.onclick = async ()=>{
    if(!confirm('EDFSの常駐動作を停止します。よろしいですか？\n(停止前に学習データは自動保存されます)')) return;
    stopBtn.disabled = true;
    stopBtn.textContent = '停止処理中…';
    try{
      const r = await fetch('/api/stop', {method:'POST'});
      const d = await r.json();
      stopBtn.textContent = d.ok ? '⏹ 停止しました(学習データ保存済み)' : '⚠️ 失敗: ' + (d.error||'');
    }catch(e){
      stopBtn.textContent = '⚠️ 通信失敗';
    }
  };
  stopBox.appendChild(stopBtn);
  locCard.appendChild(stopBox);
  root.appendChild(locCard);

  // Ver13.20: 起動設定(読み取り専用) — pythonランチャー等CLIが見えない環境向け
  const cfg = data.startup_config || {};
  if(Object.keys(cfg).length){
    const cfgDet = el('details', 'tier-advanced');
    cfgDet.appendChild(el('summary', null, '⚙️ 起動設定 (読み取り専用)'));
    const cfgBody = el('div', 'sub');
    cfgBody.style.padding = '6px 2px';
    const cfgLines = [
      ['周期 (--interval)', cfg.interval_sec !== undefined ? cfg.interval_sec + '秒' : '--'],
      ['MSTID画像種別 (--mstid-source)', cfg.mstid_source || '--'],
      ['PSKReporter対象局 (--psk-callsign)', cfg.psk_callsign || '--'],
      ['FT8バックグラウンドポーリング (--ft8-bg)', cfg.ft8_bg ? '有効' : '無効'],
      ['起動時表示モード (--ui-mode)', cfg.ui_mode || '--'],
      ['HTML配信ポート (--http-port)', cfg.html_port !== undefined ? cfg.html_port : '--'],
    ];
    cfgLines.forEach(([k, v])=>{
      const row = el('div');
      row.style.display = 'flex'; row.style.justifyContent = 'space-between'; row.style.gap = '10px';
      row.appendChild(el('span', null, k));
      row.appendChild(el('span', 'mono', String(v)));
      cfgBody.appendChild(row);
    });
    cfgBody.appendChild(el('div', 'sub small',
      '※これは起動コマンド引数(CLI)の記録のみを表示する読み取り専用パネルです。'
      +'PSKReporter対象局・FT8連携の変更は、画面右上の「⚙️ 設定」ボタンから行ってください'
      +'(再起動不要で次回サイクルから反映されます)。'));
    cfgDet.appendChild(cfgBody);
    root.appendChild(cfgDet);
  }

  // 詳細 (くわしい表示のみ)
  const det = el('details', 'tier-advanced');
  det.appendChild(el('summary', null, '⚙️ ターミナル生ログ'));
  const pre = el('pre', 'terminal'); pre.id = 'termlog';
  det.appendChild(pre);
  root.appendChild(det);

  document.getElementById('footer').textContent =
    (data.version || 'EDFS') + ' — 最終更新 ' + (data.updated_at || '');
}

async function loadLog(){
  try{
    const r = await fetch('/api/log?tail=200');
    if(!r.ok) return;
    const t = await r.text();
    const e = document.getElementById('termlog');
    if(e){ e.textContent = t; e.scrollTop = e.scrollHeight; }
  } catch(e){}
}

// Ver13.5: カウントダウン (即カウント開始 + パルスドット表示)
let cdTimer = null;
function startCd(){
  if(cdTimer) clearInterval(cdTimer);
  cdTimer = setInterval(()=>{
    const v = document.getElementById('cdv');
    const b = document.getElementById('cdb');
    if(!v || !b) return;
    if(!latestSnapshot || !latestSnapshot.cycle_info){
      // データ未受信中でも "初期化中" を秒単位で表示
      v.textContent = '初期化中…';
      return;
    }
    const ci = latestSnapshot.cycle_info;
    if(!ci.next_cycle_epoch){
      v.textContent = '待機中';
      return;
    }
    const now = Date.now()/1000;
    const drift = (ci.server_epoch || now) - now;
    const rem = ci.next_cycle_epoch - now - drift;
    if(rem <= 0){
      v.textContent = '次サイクル取得中…';
      b.style.width = '95%';
    } else {
      const m = Math.floor(rem/60), s = Math.floor(rem%60);
      v.textContent = `${m}分${String(s).padStart(2,'0')}秒`;
      b.style.width = (100*(1-rem/ci.interval_sec)) + '%';
    }
    // Es リアルタイム位置更新
    updateEsLiveMarker();
  }, 1000);
}

async function tick(){
  try{
    const r = await fetch('/api/status?ts=' + Date.now());
    if(!r.ok) throw new Error('status ' + r.status);
    const d = await r.json();
    // ESDUCT(整理版): 初回起動(コールサイン未設定)を検知したら、起動完了を
    // 待たずに(d.ready === false の段階でも)自動的に設定モーダルを開く。
    checkFirstRunSettings(d);
    // Ver13.5: ready 前でも cycle_info と暫定 map を受け取れるようにする
    if(!d.ready){
      renderStartup(d.startup_progress_pct || 10, d.startup_stage || '起動中…');
      // カウントダウンだけは動かす
      latestSnapshot = d;
    } else {
      latestSnapshot = d;
      render(d);
      loadLog();
    }
    startCd();
  } catch(e){
    console.log('fetch err', e);
    renderStartup(5, '接続待ち…');
  }
}
// Rev16.2(不具合修正): 他端末学習統合UIを#root(tick()の度に丸ごと再構築される
// 領域)の外側にある永続コンテナ(#merge-ui-persistent、静的HTMLに用意済み)へ
// 一度だけ構築する。これにより、フォルダ閲覧中や複数ファイル選択中に定期
// 更新(20秒毎)が来ても、開いている状態・選択内容が消えずに保持される。
// 【根本原因だった旧実装】従来はrender()関数内(#root配下)でmergeBox一式を
// 毎回作り直しており、tick()がrender()を呼ぶたびに(=20秒毎に)フォルダ
// 一覧・チェック状態がリセットされていた。ユーザーがフォルダ移動や複数
// ファイルの選択に20秒以上かけると、その間に発生した定期更新で表示が
// 強制的に閉じ状態へ巻き戻り、選択内容(mergeSelected)も消えるため、
// 「統合ボタンを押しても動作しない」ように見えていた(実際にはボタン自体は
// 動作していたが、UI状態が定期的に消去されていた)。Pydroid3ではフォルダ
// 階層が深い/タップ操作に時間がかかることが多く、特に発現しやすかったと
// 考えられる。
function initMergeUI(){
  const container = document.getElementById('merge-ui-persistent');
  if(!container || container.childElementCount > 0) return;  // 二重初期化防止

  const card = el('div', 'card');
  card.appendChild(el('h2', null, '🔄 他端末の学習データを統合'));
  card.appendChild(el('div', 'hint',
    '他端末(スマートフォン・タブレット等)からコピーしてきた学習データファイル'
    + '(edfs_area_history.csv等)を選択すると、今の端末にまだ無い新規サンプルのみを'
    + '安全に統合し、学習を加速できます。統合前に自動でバックアップ(.bak)を作成します。'));

  const mergeOpenBtn = el('button', 'btn ghost', '📂 フォルダを開いて選択');
  mergeOpenBtn.style.marginTop='6px';
  const mergeBrowserWrap = el('div');
  mergeBrowserWrap.style.display='none';
  mergeBrowserWrap.style.marginTop='8px';
  card.appendChild(mergeOpenBtn);
  card.appendChild(mergeBrowserWrap);
  container.appendChild(card);

  let mergeSelected = new Set();
  let mergeCurrentPath = null;

  async function mergeBrowseTo(path){
    mergeBrowserWrap.innerHTML = '';
    mergeBrowserWrap.appendChild(el('div', 'sub small', '読み込み中…'));
    let d;
    try{
      const url = path ? ('/api/browse?path=' + encodeURIComponent(path)) : '/api/browse';
      const r = await fetch(url);
      d = await r.json();
    }catch(e){
      mergeBrowserWrap.innerHTML = '';
      mergeBrowserWrap.appendChild(el('div', 'state-warn', '⚠ フォルダ一覧の取得に失敗しました'));
      return;
    }
    mergeBrowserWrap.innerHTML = '';
    if(!d.ok){
      mergeBrowserWrap.appendChild(el('div', 'state-warn', '⚠ ' + (d.error||'エラー')));
      return;
    }
    mergeCurrentPath = d.path;
    const pathRow = el('div', 'sub mono small');
    pathRow.textContent = '📍 ' + d.path;
    mergeBrowserWrap.appendChild(pathRow);
    const list = el('div');
    list.style.maxHeight='260px'; list.style.overflowY='auto';
    list.style.border='1px solid var(--line)'; list.style.borderRadius='6px';
    list.style.marginTop='6px';
    if(d.parent){
      const upRow = el('div', 'browse-row');
      upRow.style.padding='6px 8px'; upRow.style.cursor='pointer'; upRow.style.borderBottom='1px solid var(--line)';
      upRow.textContent = '⬆ ..(上の階層)';
      upRow.onclick = ()=>mergeBrowseTo(d.parent);
      list.appendChild(upRow);
    }
    (d.entries||[]).forEach(ent=>{
      const row = el('div', 'browse-row');
      row.style.padding='6px 8px'; row.style.display='flex'; row.style.alignItems='center';
      row.style.gap='8px'; row.style.borderBottom='1px solid var(--line)';
      if(ent.is_dir){
        row.style.cursor='pointer';
        row.appendChild(el('span', null, '📁'));
        row.appendChild(el('span', null, ent.name));
        row.onclick = ()=>mergeBrowseTo(ent.path);
      } else {
        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.disabled = !ent.is_learning_file;
        cb.checked = mergeSelected.has(ent.path);
        cb.onchange = ()=>{
          if(cb.checked) mergeSelected.add(ent.path); else mergeSelected.delete(ent.path);
        };
        row.appendChild(cb);
        row.appendChild(el('span', null, ent.is_learning_file ? '📄✅' : '📄'));
        const nameSpan = el('span', ent.is_learning_file ? null : 'sub', ent.name);
        row.appendChild(nameSpan);
        const metaSpan = el('span', 'sub small');
        metaSpan.style.marginLeft='auto';
        metaSpan.textContent = (ent.size!==null&&ent.size!==undefined ? Math.round(ent.size/1024)+'KB ' : '') + (ent.mtime||'');
        row.appendChild(metaSpan);
      }
      list.appendChild(row);
    });
    mergeBrowserWrap.appendChild(list);
    mergeBrowserWrap.appendChild(el('div', 'hint',
      '✅マーク付きのファイルのみ学習データとして統合できます(その他のファイルは選択できません)。'));
    const actionRow = el('div');
    actionRow.style.marginTop='8px'; actionRow.style.display='flex'; actionRow.style.gap='8px';
    const mergeRunBtn = el('button', 'btn', '🔄 選択したファイルを統合');
    const mergeCloseBtn = el('button', 'btn ghost', '閉じる');
    mergeCloseBtn.onclick = ()=>{ mergeBrowserWrap.style.display='none'; };
    mergeRunBtn.onclick = async ()=>{
      if(mergeSelected.size===0){
        alert('統合するファイルを1つ以上選択してください。');
        return;
      }
      mergeRunBtn.disabled = true;
      const origText = mergeRunBtn.textContent;
      mergeRunBtn.textContent = '統合中…';
      try{
        const r = await fetch('/api/learning-merge', {
          method:'POST', headers:{'Content-Type':'application/json'},
          body: JSON.stringify({files: Array.from(mergeSelected)})
        });
        const d2 = await r.json();
        mergeBrowserWrap.style.display = 'none';
        // 統合結果(学習がどれだけ加速したか)を1回だけ表示する専用領域。
        // このカードも#merge-ui-persistent配下(root外)に挿入するため、
        // 定期全体再描画で消えることはない。
        const resultBox = el('div', 'card');
        resultBox.style.marginTop='10px';
        resultBox.appendChild(el('h2', null, '🔄 学習データ統合 結果'));
        if(!d2.ok){
          resultBox.appendChild(el('div', 'state-warn', '⚠ 統合処理に失敗しました: ' + (d2.error||'')));
        } else {
          resultBox.appendChild(el('div', null,
            '合計 ' + (d2.total_added||0) + ' 件の新規学習サンプルを追加し、学習を加速しました。'));
          const tbl = el('table');
          const th = el('tr');
          ['ファイル','状態','追加前','追加後','加速件数'].forEach(h=>th.appendChild(el('th', null, h)));
          tbl.appendChild(th);
          (d2.results||[]).forEach(r2=>{
            const tr = el('tr');
            tr.appendChild(el('td', null, r2.file||'--'));
            tr.appendChild(el('td', r2.ok?'state-good':'state-warn', r2.ok ? 'OK' : (r2.error||'失敗')));
            tr.appendChild(el('td', 'mono', r2.before!==undefined?String(r2.before):'--'));
            tr.appendChild(el('td', 'mono', r2.after!==undefined?String(r2.after):'--'));
            tr.appendChild(el('td', 'mono', r2.added!==undefined?('+'+r2.added):'--'));
            tbl.appendChild(tr);
          });
          resultBox.appendChild(tbl);
        }
        card.parentNode.insertBefore(resultBox, card.nextSibling);
        mergeSelected = new Set();
      }catch(e){
        alert('統合処理の通信に失敗しました。');
      }
      mergeRunBtn.disabled = false;
      mergeRunBtn.textContent = origText;
    };
    actionRow.appendChild(mergeRunBtn);
    actionRow.appendChild(mergeCloseBtn);
    mergeBrowserWrap.appendChild(actionRow);
  }

  mergeOpenBtn.onclick = ()=>{
    const willShow = (mergeBrowserWrap.style.display === 'none');
    mergeBrowserWrap.style.display = willShow ? 'block' : 'none';
    if(willShow) mergeBrowseTo(mergeCurrentPath);
  };
}
// ==============================================================================
// ITPFS Ver1.0: DUCT予測カード(対流圏ダクト伝搬)
// EDFS本体の#root(20秒毎に丸ごと再構築)とは独立して、#duct-ui-persistent配下に
// 一度だけ構築し、自前のsetIntervalで更新する(initMergeUI()と同じ設計理由)。
// 「現在地から面が開通する」というEDFS本体の考え方とは異なり、DUCTは
// 任意に指定した地点そのものの見込みを調べる、という設計を明確に分けている。
// ==============================================================================
function ductLevelBadge(level){
  if(!level) return el('span','badge low','データ無');
  const pct = (level.dpp!=null) ? (' '+level.dpp+'%') : '';
  return el('span','badge '+(level.cls||'low'), (level.label||'—')+pct);
}
function ductRenderBandTable(container, timeline){
  container.innerHTML = '';
  const first = (timeline||[]).find(f=>f && f.bands);
  if(!first){ container.appendChild(el('div','hint','時間帯別データがありません')); return; }
  const bands = Object.keys(first.bands);
  const table = document.createElement('table'); table.className='duct-band-table';
  const thead = document.createElement('tr');
  thead.appendChild(el('th', null, '時刻'));
  thead.appendChild(el('th', null, 'DPP'));
  bands.forEach(b=>thead.appendChild(el('th', null, b)));
  table.appendChild(thead);
  (timeline||[]).forEach(f=>{
    if(!f) return;
    const tr = document.createElement('tr');
    tr.appendChild(el('td', null, f.t||''));
    const dppTd = document.createElement('td'); dppTd.appendChild(ductLevelBadge(f.level)); tr.appendChild(dppTd);
    bands.forEach(b=>{
      const bd = (f.bands||{})[b] || {};
      const txt = bd.trapped ? ('◎'+(bd.f_min_mhz!=null?(' '+bd.f_min_mhz+'MHz~'):'')) : '—';
      tr.appendChild(el('td', null, txt));
    });
    table.appendChild(tr);
  });
  container.appendChild(table);
}
function initDuctUI(){
  const container = document.getElementById('duct-ui-persistent');
  if(!container || container.childElementCount > 0) return;  // 二重初期化防止

  const card = el('div', 'card'); card.id = 'duct-card';
  card.appendChild(el('h2', null, '🌫️ DUCT予測(対流圏ダクト伝搬)'));
  card.appendChild(el('div', 'hint',
    '大気の温度・湿度境界による電波の異常伝搬(ダクト)の見込みです。'
    +'50/144/351/430/1200MHz帯を対象に、EDFS(Es/F2)とは完全に並行してデータ収集・学習しています。'));

  const homeWrap = el('div'); homeWrap.style.margin = '8px 0';
  card.appendChild(homeWrap);

  // ITPFS(整理版・大規模UI改修): 単体版DUCT(V6.2.1相当)を「窓の中(iframe/別
  // タブ)」で表示する方式を全面的に廃止し、ダッシュボード本体のDOMへ直接統合
  // した。DUCT_EMBED_BODY_HTML(HTML断片)をinnerHTMLで挿入した後、
  // DUCT_EMBED_SCRIPT_JSは<script>要素をcreateElement+appendChildで動的生成
  // して実行する(innerHTML経由で挿入された<script>タグは仕様上実行されない
  // ため、この方式を用いる)。
  const embedRoot = document.createElement('div');
  embedRoot.id = 'duct-embed-root';
  embedRoot.innerHTML = DUCT_EMBED_BODY_HTML;
  card.appendChild(embedRoot);
  const embedScriptTag = document.createElement('script');
  embedScriptTag.textContent = DUCT_EMBED_SCRIPT_JS;
  embedRoot.appendChild(embedScriptTag);

  // ITPFS(整理版): 学習度スケール(成熟度の視覚的ゲージ)を、折りたたみ式
  // <details>の中(=クリックして展開しないと見えない場所)から解放し、
  // 地図のすぐ下に常時表示するカードへ独立させた。
  // 【旧実装の問題点】学習成熟度・平均信頼度κ・実用までの残り時間・直近
  // サイクルの気象取得/解析成功率といった情報一式が、既定で閉じている
  // <details>(「🔧 詳しい表示 / 任意の座標を調べる」)の内側に置かれており、
  // 初回起動時にダクトタブを開いても、ユーザーが自分でクリックして展開
  // しない限り一切見えなかった。他端末データ統合(#merge-ui-persistent)・
  // 常駐停止ボタンと同様、「開かなくても常に見える」情報として扱うべき
  // ものであるため、この位置へ移動する。
  const maturityCard = el('div');
  maturityCard.style.marginTop = '10px';
  maturityCard.style.paddingTop = '10px';
  maturityCard.style.borderTop = '1px solid var(--line)';
  maturityCard.appendChild(el('div', null, '📊 DUCT学習度スケール'));
  const learnWrap = el('div'); learnWrap.style.margin = '8px 0';
  maturityCard.appendChild(learnWrap);
  card.appendChild(maturityCard);

  // ITPFS(整理版): 常駐停止ボタン。Es/F2タブ(#root配下)にのみ存在しており、
  // ダクトタブ表示中(body.mode-duct時は#rootがdisplay:noneで丸ごと隠れる)は
  // 一切アクセスできなかった不具合を修正し、ダクトタブ単独でも停止できる
  // ようにする(処理内容はEs/F2タブ側の実装と同一、/api/stopを共通利用)。
  const ductStopBox = el('div');
  ductStopBox.style.marginTop = '10px'; ductStopBox.style.paddingTop = '10px';
  ductStopBox.style.borderTop = '1px solid var(--line)';
  ductStopBox.appendChild(el('div', null, '⏹ 常駐停止'));
  ductStopBox.appendChild(el('div', 'hint',
    'Es/F2予測・DUCT予測の両方を含む、バックグラウンド動作全体を停止します。'
    +'停止前に現在までの学習データ(Es側セッションログ・校正プロファイル)を'
    +'自動的に保存してから停止します。再開にはアプリ/プログラムの再起動が必要です。'));
  const ductStopBtn = el('button', 'btn ghost', '⏹ 常駐停止');
  ductStopBtn.style.marginTop = '6px';
  ductStopBtn.style.borderColor = 'var(--hot)'; ductStopBtn.style.color = 'var(--hot)';
  ductStopBtn.onclick = async ()=>{
    if(!confirm('EDFS/DUCTの常駐動作を停止します。よろしいですか？\n(停止前に学習データは自動保存されます)')) return;
    ductStopBtn.disabled = true;
    ductStopBtn.textContent = '停止処理中…';
    try{
      const r = await fetch('/api/stop', {method:'POST'});
      const d = await r.json();
      ductStopBtn.textContent = d.ok ? '⏹ 停止しました(学習データ保存済み)' : '⚠️ 失敗: ' + (d.error||'');
    }catch(e){
      ductStopBtn.textContent = '⚠️ 通信失敗';
    }
  };
  ductStopBox.appendChild(ductStopBtn);
  card.appendChild(ductStopBox);

  // ITPFS(整理版): 「他端末の学習データを統合」機能は、このカードの直下に
  // 常設表示される#merge-ui-persistentカード(initMergeUI()で構築、Es/DUCT
  // 両対応)からご利用いただけます(DUCT較正データ dpp_calibration.json も
  // 選択可能)。タブ切替に関係なく常時同じ場所に表示されるため、ここでは
  // 案内のみ表示する。
  const ductMergeHint = el('div', 'hint');
  ductMergeHint.style.marginTop = '10px'; ductMergeHint.style.paddingTop = '10px';
  ductMergeHint.style.borderTop = '1px solid var(--line)';
  ductMergeHint.innerHTML = '🔄 他端末のDUCT学習データ(<code>dpp_calibration.json</code>)を統合するには、'
    +'このすぐ下に表示されている「他端末の学習データを統合」カードをご利用ください。';
  card.appendChild(ductMergeHint);

  const detailsEl = document.createElement('details');
  detailsEl.style.marginTop = '8px';
  const summaryEl = document.createElement('summary');
  summaryEl.textContent = '🔧 任意の座標を調べる(詳しい表示)';
  summaryEl.style.cursor = 'pointer';
  detailsEl.appendChild(summaryEl);

  detailsEl.appendChild(el('div', 'hint',
    '任意の緯度経度を入力して調べられます。EDFS本体の「現在地から面が開通していく」'
    +'予測とは別の考え方で、指定した地点そのものの見込みだけを表示します。'));
  const formWrap = el('div', 'duct-point-form');
  const latInput = document.createElement('input'); latInput.placeholder = '緯度 例:35.68';
  const lonInput = document.createElement('input'); lonInput.placeholder = '経度 例:139.77';
  const gpsBtn = el('button', 'btn ghost', '📍 現在地');
  const searchBtn = el('button', 'btn', 'この地点を調べる');
  formWrap.appendChild(latInput); formWrap.appendChild(lonInput);
  formWrap.appendChild(gpsBtn); formWrap.appendChild(searchBtn);
  detailsEl.appendChild(formWrap);

  // ITPFS Ver1.0(修正): 単体版DUCTにあった地名検索(/api/geocode)が統合版で
  // 未収録だったため、緯度経度の直接入力に加えて地名からの検索も追加する。
  const geoWrap = el('div', 'duct-point-form');
  const placeInput = document.createElement('input'); placeInput.placeholder = '地名 例:那須塩原市';
  placeInput.style.width = '14em';
  const geoBtn = el('button', 'btn ghost', '🔎 地名から検索');
  geoWrap.appendChild(placeInput); geoWrap.appendChild(geoBtn);
  detailsEl.appendChild(geoWrap);

  const resultWrap = el('div', 'duct-result'); resultWrap.style.display = 'none';
  detailsEl.appendChild(resultWrap);

  card.appendChild(detailsEl);
  container.appendChild(card);

  geoBtn.onclick = async ()=>{
    const q = placeInput.value.trim();
    if(!q){ alert('地名を入力してください'); return; }
    geoBtn.disabled = true; geoBtn.textContent = '検索中…';
    try{
      const res = await fetch('/api/duct/geocode?q='+encodeURIComponent(q));
      const data = await res.json();
      if(!data || !data.ok){ alert('見つかりませんでした: '+(data&&data.error||q)); return; }
      latInput.value = data.lat.toFixed(4); lonInput.value = data.lon.toFixed(4);
      searchBtn.click();
    }catch(e){
      alert('地名検索でエラーが発生しました: '+e);
    }finally{
      geoBtn.disabled = false; geoBtn.textContent = '🔎 地名から検索';
    }
  };

  gpsBtn.onclick = ()=>{
    if(!navigator.geolocation){ alert('このブラウザは位置情報取得に対応していません'); return; }
    navigator.geolocation.getCurrentPosition(
      pos=>{ latInput.value = pos.coords.latitude.toFixed(4); lonInput.value = pos.coords.longitude.toFixed(4); },
      err=>{ alert('現在地の取得に失敗しました: '+err.message); }
    );
  };

  searchBtn.onclick = async ()=>{
    const lat = parseFloat(latInput.value), lon = parseFloat(lonInput.value);
    if(!isFinite(lat) || !isFinite(lon)){ alert('緯度・経度を数値で入力してください'); return; }
    resultWrap.style.display = ''; resultWrap.innerHTML = '調べています…';
    try{
      const res = await fetch('/api/duct/point?lat='+encodeURIComponent(lat)+'&lon='+encodeURIComponent(lon));
      const data = await res.json();
      if(!data.ok){ resultWrap.innerHTML = ''; resultWrap.appendChild(el('div','hint', data.error||'取得に失敗しました')); return; }
      resultWrap.innerHTML = '';
      const head = el('div'); head.style.display = 'flex'; head.style.gap = '8px'; head.style.alignItems = 'center';
      head.appendChild(el('b', null, lat.toFixed(3)+', '+lon.toFixed(3)));
      head.appendChild(ductLevelBadge(data.now));
      resultWrap.appendChild(head);
      if(data.summary) resultWrap.appendChild(el('div', 'hint', data.summary));
      resultWrap.appendChild(el('div', 'hint', data.reliability_note||''));
      if(data.components){
        const c = data.components;
        resultWrap.appendChild(el('div', 'hint',
          '構造:'+c.structure+' / 安定性:'+c.stability+' / 結合:'+c.coupling+' / 乱流:'+c.turbulence
          +' (信頼度κ='+(data.confidence_kappa!=null?data.confidence_kappa:'—')+')'));
      }
      const tblWrap = document.createElement('div');
      resultWrap.appendChild(tblWrap);
      ductRenderBandTable(tblWrap, data.timeline);
    }catch(e){
      resultWrap.innerHTML = ''; resultWrap.appendChild(el('div', 'hint', '通信エラー: '+e));
    }
  };

  async function refreshHome(){
    try{
      const res = await fetch('/api/duct/status');
      const data = await res.json();
      homeWrap.innerHTML = '';
      if(!data.ok){ homeWrap.appendChild(el('div', 'hint', 'DUCT予測: 起動中…')); return; }
      // Rev16.9: 起動処理(duct_start_background)自体が失敗した場合、
      // bootstrap_ok=falseが返る。この場合は「データ収集中…」ではなく、
      // 具体的な失敗理由と「自動では解消しない」旨を明示する
      // (従来はcycles_computed=0のままいつまでも「データ収集中…」としか
      // 表示されず、利用者が「起動失敗で永久に進まない」ことに気付けなかった)。
      if(data.bootstrap_ok === false){
        const errBox = el('div');
        errBox.style.background = '#3a0d0d'; errBox.style.border = '1px solid #ff5252';
        errBox.style.borderRadius = '8px'; errBox.style.padding = '8px 10px';
        errBox.appendChild(el('div', null, '⚠ DUCT予測エンジンの起動に失敗しました'));
        errBox.appendChild(el('div', 'hint', '原因: ' + (data.bootstrap_error || '不明')));
        errBox.appendChild(el('div', 'hint', 'この状態は自動的には解消しません。アプリ/プログラムの再起動が必要です。'));
        homeWrap.appendChild(errBox);
        return;
      }
      const head = el('div'); head.style.display = 'flex'; head.style.gap = '8px'; head.style.alignItems = 'center';
      head.appendChild(el('span', null, '自局(現在地)付近:'));
      head.appendChild(data.home ? ductLevelBadge(data.home.level) : el('span', 'hint', 'データ収集中…'));
      homeWrap.appendChild(head);
      if(data.home && data.home.summary) homeWrap.appendChild(el('div', 'hint', data.home.summary));
      const lq = data.learning || {};
      learnWrap.innerHTML = '';
      // ITPFS(整理版): 学習成熟度を、Es/F2タブの「学習の品質と成熟度」カードと
      // 同じ視覚言語(.progress/<span>のゲージバー)で表示する(旧実装は文字列
      // のみでスケール感が伝わりにくかった)。
      const maturityPct = (lq.maturity_pct!=null) ? Math.max(0, Math.min(100, lq.maturity_pct)) : 0;
      const gaugeRow = el('div');
      gaugeRow.style.display = 'flex'; gaugeRow.style.alignItems = 'baseline'; gaugeRow.style.gap = '8px';
      gaugeRow.appendChild(el('div', 'pct-big ' + (maturityPct>=70?'state-good':maturityPct>=30?'state-warn':'state-low'),
        (lq.maturity_pct!=null?lq.maturity_pct+'%':'—')));
      gaugeRow.appendChild(el('div', 'sub', '学習成熟度(κ=0.7到達で実用レベル)'));
      learnWrap.appendChild(gaugeRow);
      const gaugeBar = el('div', 'progress');
      const gaugeSpan = el('span'); gaugeSpan.style.width = maturityPct + '%';
      gaugeBar.appendChild(gaugeSpan);
      learnWrap.appendChild(gaugeBar);
      learnWrap.appendChild(el('div', 'hint',
        '平均信頼度κ='+(lq.mean_confidence!=null?lq.mean_confidence:'—')
        +' / 実用(κ0.7)まで:'+(lq.eta_hours!=null?(lq.eta_hours<=0?'到達済み ✅':'あと約'+lq.eta_hours+'時間'):'推定中')
        +' / 被覆セル数:'+(lq.covered_cells!=null?lq.covered_cells:'—')));
      // ITPFS Ver1.0(修正): 「データ取得・表示が遅れる」という報告に対し、
      // 直近サイクルの気象取得成功率・解析成功率をダッシュボード上でも
      // 直接確認できるようにする(/api/duct/mapを開かなくても分かる)。
      const cd = data.cycle_diag || {};
      if(cd.cycles_run){
        learnWrap.appendChild(el('div', 'hint',
          '直近サイクル: 実行'+cd.cycles_run+'回 / 解像度='
          +(cd.last_step_deg!=null?cd.last_step_deg+'度':'—')
          +(cd.last_is_coarse_pass?'(初回・粗い)':'(通常)')
          +' / 気象取得='
          +cd.last_weather_ok+'/'+cd.last_points_total+'地点 / 解析成功='
          +cd.last_analysis_valid+'/'+cd.last_points_total+'地点 / 所要'
          +cd.last_duration_sec+'秒'));
        if(cd.last_weather_ok === 0){
          const warnBox = el('div');
          warnBox.style.color = '#ffb020';
          warnBox.textContent = '⚠ 気象データが全地点で取得できていません(通信環境またはOpen-Meteo側の制限の可能性)。この状態が続くと地図は更新されません。';
          learnWrap.appendChild(warnBox);
        }
      }else if(data.first_cycle_started){
        learnWrap.appendChild(el('div', 'hint', '初回サイクルを実行中です…'));
      }
    }catch(e){
      homeWrap.innerHTML = ''; homeWrap.appendChild(el('div', 'hint', 'DUCT予測: 通信エラー'));
    }
  }
  refreshHome();
  setInterval(refreshHome, REFRESH_MS);
}

// ITPFS Ver1.0(修正): ダーク/ライト切替ボタンの初期化。localStorageに保存し、
// 次回起動時も選択した表示モードを維持する(既定はダーク=従来通り)。
function initThemeToggle(){
  const KEY = 'itpfs-theme';
  const root = document.documentElement;
  const btn = document.getElementById('theme-toggle');
  function apply(theme){
    if(theme === 'light'){ root.setAttribute('data-theme', 'light'); }
    else{ root.removeAttribute('data-theme'); }
  }
  let saved = null;
  try{ saved = localStorage.getItem(KEY); }catch(_){}
  apply(saved === 'light' ? 'light' : 'dark');
  if(btn){
    btn.addEventListener('click', ()=>{
      const now = root.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
      apply(now);
      try{ localStorage.setItem(KEY, now); }catch(_){}
      // ITPFS Ver1.0(修正): 47都道府県マップはJS描画のSVGでCSS変数を
      // getComputedStyle経由で読むため、切替直後すぐには反映されない
      // (次の20秒ごとの自動再構築まで旧配色のまま)。tick()を即時に
      // 呼び直して見た目を即座に切替に追従させる。
      try{ tick(); }catch(_){}
    });
  }
}
// ==============================================================================
// ESDUCT(整理版): 右上「⚙️ 設定」ボタン — コールサイン/FT8連携設定モーダル
// ------------------------------------------------------------------------------
// 【解決した不具合】従来、コールサイン/FT8連携を変更する手段がGUI上に一切
// 存在せず(読み取り専用の「起動設定」パネルのみ)、指示されていた「右上の
// 設定ボタン」自体が未実装だった。本モーダルはバックエンドの
// /api/es/settings (GET/POST) と連動し、保存内容はes_settings.json経由で
// 永続化され、次回サイクル(数十秒以内、main()側の即時反映ロジック)から
// 確実に反映される。
// ==============================================================================
let _settingsModalFirstOpenDone = false;
function openSettingsModal(isFirstRun){
  const overlay = document.getElementById('settings-modal-overlay');
  const note = document.getElementById('settings-modal-firstrun-note');
  const msg = document.getElementById('settings-modal-msg');
  if(!overlay) return;
  msg.textContent = '';
  note.style.display = isFirstRun ? 'block' : 'none';
  fetch('/api/es/settings').then(r=>r.json()).then(d=>{
    const s = (d && d.settings) || {};
    document.getElementById('settings-callsign-input').value = s.psk_callsign || '';
    document.getElementById('settings-ft8-enabled').checked = (s.psk_enabled !== false);
    document.getElementById('settings-ft8-bg').checked = !!s.ft8_bg;
    const sel = document.getElementById('settings-ft8-bg-interval');
    sel.value = String(s.ft8_bg_interval || 180);
  }).catch(()=>{});
  overlay.style.display = 'flex';
}
function closeSettingsModal(){
  const overlay = document.getElementById('settings-modal-overlay');
  if(overlay) overlay.style.display = 'none';
}
function initSettingsModal(){
  const openBtn = document.getElementById('settings-toggle');
  const cancelBtn = document.getElementById('settings-modal-cancel');
  const saveBtn = document.getElementById('settings-modal-save');
  const overlay = document.getElementById('settings-modal-overlay');
  if(openBtn) openBtn.addEventListener('click', ()=>openSettingsModal(false));
  if(cancelBtn) cancelBtn.addEventListener('click', closeSettingsModal);
  if(overlay) overlay.addEventListener('click', (e)=>{ if(e.target === overlay) closeSettingsModal(); });
  if(saveBtn){
    saveBtn.addEventListener('click', async ()=>{
      const msg = document.getElementById('settings-modal-msg');
      const callsign = document.getElementById('settings-callsign-input').value.trim().toUpperCase();
      if(callsign && !/^[A-Z0-9/-]{1,32}$/.test(callsign)){
        msg.textContent = '⚠ コールサインの形式が不正です(英数字と/-のみ)';
        msg.style.color = '#ff5252';
        return;
      }
      const payload = {
        psk_callsign: callsign,
        psk_enabled: document.getElementById('settings-ft8-enabled').checked,
        ft8_bg: document.getElementById('settings-ft8-bg').checked,
        ft8_bg_interval: parseInt(document.getElementById('settings-ft8-bg-interval').value, 10),
      };
      saveBtn.disabled = true;
      msg.style.color = '';
      msg.textContent = '保存中…';
      try{
        const r = await fetch('/api/es/settings', {method:'POST', headers:{'Content-Type':'application/json'},
          body: JSON.stringify(payload)});
        const d = await r.json();
        if(d.ok){
          msg.style.color = '#26c281';
          msg.textContent = '✅ 保存しました。次回サイクルから反映されます(通常30秒〜数分以内)。';
          setTimeout(closeSettingsModal, 1500);
        } else {
          msg.style.color = '#ff5252';
          msg.textContent = '⚠ 保存失敗: ' + (d.error || '');
        }
      }catch(e){
        msg.style.color = '#ff5252';
        msg.textContent = '⚠ 通信エラー';
      }
      saveBtn.disabled = false;
    });
  }
}
function checkFirstRunSettings(data){
  // tick()の定期取得結果から、初回起動(needs_setup)かどうかを判定し、
  // セッション中1回だけ自動的に設定モーダルを開く。
  if(_settingsModalFirstOpenDone) return;
  const es = data && data.es_settings;
  if(es && es.needs_setup){
    _settingsModalFirstOpenDone = true;
    openSettingsModal(true);
  }
}
initThemeToggle();
initTierSwitch();
initSettingsModal();
initMergeUI();
initDuctUI();
reqGps();
setTimeout(reqGps, 5000);
tick();
setInterval(tick, REFRESH_MS);
window.addEventListener('beforeunload', ()=>{
  try{ fetch('/api/flush', {method:'POST', keepalive:true}); } catch(_){}
});
</script>
</body>
</html>
"""

_v13_JS_BODY = _v13_JS_BODY.replace(
    "__DUCT_EMBED_BODY_HTML_JSON__",
    json.dumps(_DUCT_EMBED_BODY_HTML, ensure_ascii=True).replace("</script", "<\\/script")
).replace(
    "__DUCT_EMBED_SCRIPT_JS_JSON__",
    json.dumps(_DUCT_EMBED_SCRIPT_JS, ensure_ascii=True).replace("</script", "<\\/script")
)

_v13_HTML = (_v13_HTML_HEAD + "<script>" + _v13_PREFECTURES_JS
             + "\n</script>\n" + _v13_JS_BODY)


class _V13HttpHandler(_v13_http_server.SimpleHTTPRequestHandler):
    server_version = "EDFS-GUI/13.5"

    _MAX_BODY = 64 * 1024
    _MAX_TAIL = 200
    _ALLOWED_GPS_SOURCES = {'browser', 'browser-geolocation', 'manual', 'termux'}
    _ALLOWED_GET_PATHS = {'/', '/index.html', '/edfs_dashboard.html', '/api/status', '/api/log', '/api/browse',
                         '/api/duct/status', '/api/duct/point', '/api/duct/map', '/api/duct/settings',
                         '/api/duct/hepburn', '/api/duct/geocode', '/api/duct/pass',
                         '/api/duct/antenna', '/api/duct/radar',
                         '/api/duct/celldata',
                         '/api/es/settings'}

    def log_message(self, fmt, *args):
        pass

    def _send_json(self, obj, code=200):
        body = _v13_json.dumps(obj, ensure_ascii=False).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        self.wfile.write(body)

    def _send_text(self, text, code=200, ctype='text/plain; charset=utf-8'):
        body = text.encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        try:
            path = _v13_urlparse(self.path).path
            if path not in self._ALLOWED_GET_PATHS:
                self._send_json({'ok': False, 'error': 'not found'}, code=404)
                return
            if path == '/api/status':
                try:
                    with _v13_snapshot_lock:
                        payload = _sanitize_for_json(dict(_v13_snapshot))
                    self._send_json(payload)
                except Exception as e:
                    self._send_json({
                        'ready': False,
                        'startup_stage': f'一時的な内部エラー(次回更新で回復見込み): '
                                         f'{type(e).__name__}: {e}',
                        'startup_progress_pct': 10,
                        'startup_progress_error': True,
                    })
                return
            if path == '/api/log':
                with _v13_term_lock:
                    txt = _v13_term_buf.getvalue()
                txt = _v13_re.sub(r'\x1b\[[0-9;]*m', '', txt)
                q = _v13_urlparse(self.path).query
                m = _v13_re.search(r'tail=(\d+)', q or '')
                if m:
                    n = min(int(m.group(1)), self._MAX_TAIL)
                    lines = txt.splitlines()
                    txt = '\n'.join(lines[-n:]) if len(lines) > n else txt
                self._send_text(txt)
                return
            if path in ('/', '/index.html', '/edfs_dashboard.html'):
                self._send_text(_v13_HTML, ctype='text/html; charset=utf-8')
                return
            if path == '/api/browse':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                target = qs.get('path', [None])[0]
                result = list_browse_dir(target)
                self._send_json(result)
                return
            if path == '/api/duct/status':
                try:
                    self._send_json(_sanitize_for_json(duct_status_summary()))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/point':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                try:
                    lat = float(qs.get('lat', ['nan'])[0])
                    lon = float(qs.get('lon', ['nan'])[0])
                    if not (_v13_math.isfinite(lat) and _v13_math.isfinite(lon)):
                        raise ValueError()
                    if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
                        raise ValueError()
                except (ValueError, IndexError):
                    self._send_json({'ok': False, 'error': 'lat/lonが不正です'}, code=400)
                    return
                try:
                    self._send_json(_sanitize_for_json(duct_point_query(lat, lon)))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/celldata':
                try:
                    with _RESULTS_LOCK:
                        results = list(_LATEST_RESULTS)
                    cal = load_calibration(CALIBRATION_FILE)
                    cdiag = dict(globals().get('_duct_cycle_diag', {}) or {})
                    region = dict(DEFAULT_REGION)
                    if cdiag.get('last_step_deg') is not None:
                        region['step_deg'] = cdiag['last_step_deg']
                    cells, frame_labels, lq, eta_s = build_duct_celldata(results, region, cal)
                    if cells is None:
                        bootstrap = dict(globals().get('_duct_bootstrap_diag', {}) or {})
                        self._send_json({
                            'ok': False, 'ready': False,
                            'reason': 'no_valid_points',
                            'bootstrap_ok': bootstrap.get('ok'),
                            'bootstrap_error': bootstrap.get('error'),
                            'first_cycle_started': bootstrap.get('first_cycle_started'),
                            'first_cycle_done': bootstrap.get('first_cycle_done'),
                            'cycle_diag': cdiag,
                        })
                        return
                    self._send_json(_sanitize_for_json({
                        'ok': True, 'ready': True,
                        'cells': cells, 'frame_labels': frame_labels,
                        'region': region, 'step_deg': region['step_deg'],
                        'tx_height_m': TX_HEIGHT_M,
                        'band_list': [b.replace('MHz', '') for b in BANDS],
                        'maturity_pct': lq['maturity_pct'], 'eta_s': eta_s,
                        'covered_cells': lq['covered_cells'],
                        'mean_confidence': lq['mean_confidence'],
                        'n_pts': len(cells),
                        'generated_at': datetime.now().strftime('%m/%d %H:%M:%S'),
                        'interval_sec': get_settings().get('interval_sec', 300),
                    }))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/map':
                _DUCT_BUILD_MARKER = 'DUCT_BUILD_MARKER=OSM_TILES_REV16_7'
                try:
                    html_text = None
                    if os.path.exists(OUTPUT_HTML):
                        with open(OUTPUT_HTML, 'rb') as f:
                            candidate = f.read().decode('utf-8', errors='replace')
                        if _DUCT_BUILD_MARKER in candidate:
                            html_text = candidate
                        else:
                            global _duct_stale_map_last_logged_mtime
                            _mtime = os.path.getmtime(OUTPUT_HTML)
                            if _duct_stale_map_last_logged_mtime != _mtime:
                                _duct_stale_map_last_logged_mtime = _mtime
                                print(f"  [DUCT] {os.path.basename(OUTPUT_HTML)}は古い形式"
                                      f"(ビルドマーカー無し)のため配信をスキップします。"
                                      f"次回の予測サイクルで新しい地図に上書きされるまで"
                                      f"お待ちください。")
                    if html_text is not None:
                        self._send_text(html_text, ctype='text/html; charset=utf-8')
                    else:
                        diag = dict(globals().get('_duct_bootstrap_diag', {}) or {})
                        if diag.get('ok') is False:
                            self._send_text(
                                '<html><body style="background:#0b1220;color:#e6ecff;'
                                'font-family:sans-serif;padding:2em">'
                                '<h3 style="color:#ff5252">⚠ DUCT予測エンジンの起動に失敗しました</h3>'
                                f'<p>原因: {_html_mod.escape(str(diag.get("error") or "不明"))}</p>'
                                '<p>この状態は自動的には解消しません。アプリ/プログラムを'
                                '再起動してください。改善しない場合は、'
                                'ログディレクトリへの書き込み権限やストレージ空き容量を'
                                'ご確認ください。</p>'
                                '</body></html>', ctype='text/html; charset=utf-8')
                        elif diag.get('first_cycle_started') and not diag.get('first_cycle_done'):
                            cdiag = dict(globals().get('_duct_cycle_diag', {}) or {})
                            detail = (
                                f"<p>直近のサイクル状況: 実行回数={cdiag.get('cycles_run', 0)}回 / "
                                f"解像度={cdiag.get('last_step_deg', '—')}度"
                                f"{'(初回・粗い)' if cdiag.get('last_is_coarse_pass') else '(通常)'} / "
                                f"気象取得成功={cdiag.get('last_weather_ok', '—')}"
                                f"/{cdiag.get('last_points_total', '—')}地点 / "
                                f"解析成功={cdiag.get('last_analysis_valid', '—')}"
                                f"/{cdiag.get('last_points_total', '—')}地点 / "
                                f"前回所要時間={cdiag.get('last_duration_sec', '—')}秒</p>"
                                if cdiag.get('cycles_run') else
                                "<p>まだ1回もサイクルが完了していません(初回実行中)。</p>"
                            )
                            weather_hint = ""
                            if cdiag.get('cycles_run') and cdiag.get('last_weather_ok') == 0:
                                weather_hint = (
                                    '<p style="color:#ffb020">⚠ 前回サイクルは気象データが'
                                    '全地点で取得できていません。通信環境、または'
                                    'Open-Meteo側の一時的な制限の可能性があります。'
                                    'コンソールの「[気象] バッチ失敗」行に理由が出ています。</p>')
                            self._send_text(
                                '<html><body style="background:#0b1220;color:#e6ecff;'
                                'font-family:sans-serif;padding:2em">'
                                '<h3>DUCT地図を生成中です。数分後に再読込してください</h3>'
                                '<p>予測サイクルは起動済みで、現在初回計算を実行中です。'
                                '気象データ取得・初回較正には数分かかる場合があります。'
                                'この画面が1時間以上続く場合は、コンソールに'
                                '「[DUCT] 予測サイクルで例外発生」または'
                                '「[予測] 地点解析で例外発生」という行が繰り返し出ていないか'
                                'ご確認ください(出ていれば、その内容が真の原因です)。</p>'
                                + detail + weather_hint +
                                '</body></html>', ctype='text/html; charset=utf-8')
                        else:
                            cdiag = dict(globals().get('_duct_cycle_diag', {}) or {})
                            detail = ""
                            weather_hint = ""
                            if cdiag.get('cycles_run'):
                                detail = (
                                    f"<p>直近のサイクル状況: 実行回数={cdiag.get('cycles_run', 0)}回 / "
                                    f"解像度={cdiag.get('last_step_deg', '—')}度"
                                    f"{'(初回・粗い)' if cdiag.get('last_is_coarse_pass') else '(通常)'} / "
                                    f"気象取得成功={cdiag.get('last_weather_ok', '—')}"
                                    f"/{cdiag.get('last_points_total', '—')}地点 / "
                                    f"解析成功={cdiag.get('last_analysis_valid', '—')}"
                                    f"/{cdiag.get('last_points_total', '—')}地点 / "
                                    f"前回所要時間={cdiag.get('last_duration_sec', '—')}秒</p>")
                                if cdiag.get('last_weather_ok') == 0:
                                    weather_hint = (
                                        '<p style="color:#ffb020">⚠ 前回サイクルは気象データが'
                                        '全地点で取得できていません。通信環境、または'
                                        'Open-Meteo側の一時的な制限の可能性があります。'
                                        'コンソールの「[気象] バッチ失敗」行に理由が出ています。'
                                        'この状態が続く場合、地図は更新されません。</p>')
                            self._send_text(
                                '<html><body style="background:#0b1220;color:#e6ecff;'
                                'font-family:sans-serif;padding:2em">'
                                '<h3>DUCT地図を生成中です。数分後に再読込してください</h3>'
                                '<p>起動直後は気象データ取得・初回較正に少し時間がかかります。'
                                'また、以前のバージョンで生成された古い地図ファイルが残って'
                                'いた場合は、最初の予測サイクルが完了するまでこの画面が'
                                '表示されます(コンソールに[予測]の進捗ログが出ます)。</p>'
                                + detail + weather_hint +
                                '</body></html>', ctype='text/html; charset=utf-8')
                except Exception as e:
                    self._send_text(f'error: {e}', code=500)
                return
            if path == '/api/duct/settings':
                try:
                    self._send_json({'ok': True, 'settings': get_settings(),
                                     'choices': SETTINGS_CHOICES})
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/es/settings':
                try:
                    s = get_es_settings()
                    self._send_json({
                        'ok': True,
                        'settings': s,
                        'needs_setup': bool(globals().get('_es_settings_needs_setup', False)) and not s.get('setup_done'),
                        'ft8_bg_active': is_ft8_bg_active(),
                    })
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/hepburn':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                try:
                    lat = float(qs.get('lat', ['nan'])[0])
                    lon = float(qs.get('lon', ['nan'])[0])
                    if _v13_math.isnan(lat) or _v13_math.isnan(lon):
                        raise ValueError()
                except (ValueError, IndexError):
                    self._send_json({'ok': False, 'error': 'lat/lonが必要です'}, code=400)
                    return
                try:
                    self._send_json(_sanitize_for_json(hepburn_lookup(lat, lon)))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/geocode':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                query_str = (qs.get('q', [''])[0] or '').strip()
                if not query_str:
                    self._send_json({'ok': False, 'error': 'qが必要です'}, code=400)
                    return
                try:
                    self._send_json(_sanitize_for_json(geocode_place_jp(query_str)))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/pass':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                try:
                    lat1 = float(qs.get('lat1', ['nan'])[0]); lon1 = float(qs.get('lon1', ['nan'])[0])
                    lat2 = float(qs.get('lat2', ['nan'])[0]); lon2 = float(qs.get('lon2', ['nan'])[0])
                    if any(_v13_math.isnan(v) for v in (lat1, lon1, lat2, lon2)):
                        raise ValueError()
                except (ValueError, IndexError):
                    self._send_json({'ok': False, 'error': 'lat1/lon1/lat2/lon2が必要です'}, code=400)
                    return
                try:
                    with _RESULTS_LOCK:
                        res = list(_LATEST_RESULTS)
                    self._send_json(_sanitize_for_json(analyze_pass(res, lat1, lon1, lat2, lon2)))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/antenna':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                try:
                    lat = float(qs.get('lat', ['nan'])[0]); lon = float(qs.get('lon', ['nan'])[0])
                    if _v13_math.isnan(lat) or _v13_math.isnan(lon):
                        raise ValueError()
                    h = qs.get('height')
                    height_m = float(h[0]) if h else None
                except (ValueError, IndexError):
                    self._send_json({'ok': False, 'error': 'lat/lon/heightが不正です'}, code=400)
                    return
                try:
                    with _RESULTS_LOCK:
                        r = _nearest_result(_LATEST_RESULTS, lat, lon)
                    self._send_json(_sanitize_for_json(antenna_height_curve(r, height_m)))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            if path == '/api/duct/radar':
                q = _v13_urlparse(self.path).query
                qs = _v13_parse_qs(q or '')
                try:
                    lat = float(qs.get('lat', ['nan'])[0]); lon = float(qs.get('lon', ['nan'])[0])
                    if _v13_math.isnan(lat) or _v13_math.isnan(lon):
                        raise ValueError()
                except (ValueError, IndexError):
                    self._send_json({'ok': False, 'error': 'lat/lonが必要です'}, code=400)
                    return
                try:
                    with _RESULTS_LOCK:
                        res = list(_LATEST_RESULTS)
                    self._send_json(_sanitize_for_json(field_radar(res, lat, lon)))
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return
            self._send_json({'ok': False, 'error': 'not found'}, code=404)
        except Exception as e:
            self._send_text(f'error: {e}', code=500)

    def do_POST(self):
        try:
            path = _v13_urlparse(self.path).path
            length = int(self.headers.get('Content-Length', '0') or '0')
            if length > self._MAX_BODY:
                self._send_json({'ok': False, 'error': 'request too large'}, code=413)
                return
            raw = self.rfile.read(length) if length > 0 else b''

            if path == '/api/gps':
                try:
                    data = _v13_json.loads(raw.decode('utf-8'))
                    lat = float(data['lat'])
                    lon = float(data['lon'])
                except Exception:
                    self._send_json({'ok': False, 'error': 'invalid request body'}, code=400)
                    return
                if not (_v13_math.isfinite(lat) and _v13_math.isfinite(lon)):
                    self._send_json({'ok': False, 'error': 'non-finite GPS value'}, code=400)
                    return
                if not (-90.0 <= lat <= 90.0):
                    self._send_json({'ok': False, 'error': 'latitude out of range'}, code=400)
                    return
                if not (-180.0 <= lon <= 180.0):
                    self._send_json({'ok': False, 'error': 'longitude out of range'}, code=400)
                    return
                source = str(data.get('source', 'browser'))[:32]
                if source not in self._ALLOWED_GPS_SOURCES:
                    source = 'browser'
                with _v13_gps_lock:
                    _v13_pending_gps.update({'lat': lat, 'lon': lon,
                                             'source': source, 'ts': _v13_time.time()})
                try:
                    label = "手動入力(ブラウザ)" if source == 'manual' else "ブラウザGPS"
                    ui_set_manual_location(lat, lon, label=label)
                    globals()['_ui_location_source'] = label
                    globals()['_ui_gps_status'] = "取得成功(browser)" if source != 'manual' else "手動設定済み"
                    globals()['_ui_gps_reason'] = source
                    globals()['_ui_location_dialog_done'] = True
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                    return
                self._send_json({'ok': True, 'lat': lat, 'lon': lon})
                return

            if path == '/api/flush':
                result = {'ok': True, 'items': []}
                try:
                    if _v13_active_session_logger is not None:
                        try:
                            _v13_active_session_logger.flush()
                            result['items'].append('session_log flushed')
                        except Exception as e:
                            result['items'].append(f'session_log FAILED: {e}')
                    try:
                        if _v13_active_profile_path is not None and _v13_active_profile_data:
                            _v13_active_profile_path.write_text(
                                _v13_json.dumps(_v13_active_profile_data,
                                                ensure_ascii=False, indent=2),
                                encoding='utf-8')
                            result['items'].append('calibration_profile.json saved')
                    except Exception as e:
                        result['items'].append(f'profile FAILED: {e}')
                    try:
                        if 'update_edfs_maturity' in globals():
                            update_edfs_maturity()
                            result['items'].append('edfs_maturity updated')
                    except Exception as e:
                        result['items'].append(f'edfs FAILED: {e}')
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                    return
                self._send_json(result)
                return

            if path == '/api/stop':
                flush_items = []
                try:
                    if _v13_active_session_logger is not None:
                        _v13_active_session_logger.flush()
                        flush_items.append('session_log flushed')
                except Exception as e:
                    flush_items.append(f'session_log FAILED: {e}')
                try:
                    if _v13_active_profile_path is not None and _v13_active_profile_data:
                        _v13_active_profile_path.write_text(
                            _v13_json.dumps(_v13_active_profile_data,
                                            ensure_ascii=False, indent=2),
                            encoding='utf-8')
                        flush_items.append('calibration_profile.json saved')
                except Exception as e:
                    flush_items.append(f'profile FAILED: {e}')
                try:
                    self._send_json({'ok': True, 'message': '常駐停止を受け付けました',
                                      'flush_items': flush_items})
                except Exception:
                    pass

                def _delayed_stop():
                    _v13_time.sleep(0.3)
                    try:
                        request_edfs_stop()
                    except Exception:
                        pass
                _v13_threading.Thread(target=_delayed_stop, daemon=True).start()
                return

            if path == '/api/learning-merge':
                try:
                    data = _v13_json.loads(raw.decode('utf-8'))
                    file_paths = data.get('files', [])
                    if not isinstance(file_paths, list) or not file_paths:
                        self._send_json({'ok': False, 'error': 'files未指定です'}, code=400)
                        return
                except Exception:
                    self._send_json({'ok': False, 'error': 'invalid request body'}, code=400)
                    return
                log_dir = _v13_active_log_dir if _v13_active_log_dir is not None else _v13_Path(DEFAULT_DATA_DIR)
                result = run_learning_merge(log_dir, file_paths)
                self._send_json({'ok': True, **result})
                return

            if path == '/api/duct/settings':
                try:
                    patch = _v13_json.loads(raw.decode('utf-8')) if raw else {}
                except Exception:
                    self._send_json({'ok': False, 'error': 'invalid request body'}, code=400)
                    return
                try:
                    applied = update_settings(patch)
                    self._send_json({'ok': True, 'applied': applied, 'settings': get_settings()})
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return

            if path == '/api/es/settings':
                try:
                    patch = _v13_json.loads(raw.decode('utf-8')) if raw else {}
                except Exception:
                    self._send_json({'ok': False, 'error': 'invalid request body'}, code=400)
                    return
                try:
                    if 'psk_callsign' in patch:
                        cs = str(patch.get('psk_callsign') or '')
                        if len(cs) > 32 or not re.match(r'^[A-Za-z0-9/\-]*$', cs):
                            self._send_json({'ok': False, 'error':
                                             'コールサインの形式が不正です(英数字と/-のみ、32文字以内)'}, code=400)
                            return
                    patch['mark_setup_done'] = True
                    applied = update_es_settings(patch)
                    self._send_json({'ok': True, 'applied': applied, 'settings': get_es_settings()})
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return

            if path == '/api/duct/recompute':
                try:
                    _RECOMPUTE_NOW.set()
                    self._send_json({'ok': True})
                except Exception as e:
                    self._send_json({'ok': False, 'error': str(e)}, code=500)
                return

            self._send_json({'ok': False, 'error': 'unknown path'}, code=404)
        except Exception as e:
            self._send_json({'ok': False, 'error': str(e)}, code=500)


class _V13EnhancedBridge(HtmlDashboardBridge):
    def write_snapshot(self, text):
        try:
            self.html_path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.html_path.with_suffix(self.html_path.suffix + '.tmp')
            tmp.write_text(_v13_HTML, encoding='utf-8')
            _v13_os.replace(tmp, self.html_path)
        except Exception as e:
            if not getattr(self, '_warned', False):
                self._warned = True
                print(f"  警告: HTMLダッシュボード書き出しに失敗 ({e})")

        try:
            snap = _v13_build_snapshot()
            with _v13_snapshot_lock:
                _v13_snapshot.clear()
                _v13_snapshot.update(snap)
            if snap.get('ready'):
                with _v13_cycle_lock:
                    _v13_cycle_state['last_cycle_end_epoch'] = _v13_time.time()
                _v13_update_history(snap)
        except Exception as e:
            with _v13_snapshot_lock:
                _v13_snapshot['ready'] = False
                _v13_snapshot['startup_stage'] = f'初期化中... ({e})'
                _v13_snapshot['startup_progress_pct'] = 20

        if text:
            with _v13_term_lock:
                _v13_term_buf.write(text[-4000:] if len(text) > 4000 else text)
                v = _v13_term_buf.getvalue()
                if len(v) > 100000:
                    _v13_term_buf.seek(0); _v13_term_buf.truncate(0)
                    _v13_term_buf.write(v[-80000:])

    def start_http_server(self):
        if not self.enable_http or self._httpd is not None:
            return None
        directory = str(self.html_path.parent)
        handler = _v13_functools.partial(_V13HttpHandler, directory=directory)
        last_err = None
        for attempt_port in range(self.port, self.port + 20):
            try:
                self._httpd = _V13ReusableTCPServer(
                    (DASHBOARD_BIND_HOST, attempt_port), handler)
                self.port = attempt_port
                break
            except OSError as e:
                last_err = e
                continue
        if self._httpd is None:
            print(f"  警告: HTTPサーバーの起動に失敗しました "
                  f"(ポート{self.port}-{self.port+19}を試行, 最終エラー: {last_err})")
            return None
        self._thread = _v13_threading.Thread(target=self._httpd.serve_forever, daemon=True)
        self._thread.start()
        return self.port


def _v13_patched_ui_resolve_location():
    if globals().get('_ui_manual_location', None) is not None:
        globals()['_ui_last_location'] = globals()['_ui_manual_location']
        return globals()['_ui_last_location']
    if globals().get('_ui_last_location', None) is not None:
        return globals()['_ui_last_location']

    deadline = _v13_time.time() + 30
    while _v13_time.time() < deadline:
        with _v13_gps_lock:
            if _v13_pending_gps:
                lat = _v13_pending_gps['lat']
                lon = _v13_pending_gps['lon']
                src = _v13_pending_gps.get('source', 'browser')
                ui_set_manual_location(lat, lon, label='ブラウザGPS')
                globals()['_ui_location_source'] = "ブラウザGPS"
                globals()['_ui_gps_status'] = "取得成功(browser)"
                globals()['_ui_gps_reason'] = src
                return globals()['_ui_last_location']
        _v13_time.sleep(0.5)

    globals()['_ui_gps_status'] = "未取得"
    globals()['_ui_gps_reason'] = "ブラウザから位置未受信 — ブラウザで位置を許可してください"
    globals()['_ui_location_dialog_done'] = True
    return None


def _v13_patched_choose_mode(default_mode='NORMAL', skip_menu=False,
                              profile=None, use_color=True):
    try:
        globals()['_ui_display_mode'] = globals().get('UI_MODE_NORMAL', 'NORMAL')
    except Exception:
        pass
    return


_v13_original_session_init = SessionLogger.__init__
_v13_original_auto_cal = auto_calibrate_profile


def _v13_patched_session_init(self, *args, **kwargs):
    _v13_original_session_init(self, *args, **kwargs)
    global _v13_active_session_logger, _v13_active_log_dir
    _v13_active_session_logger = self
    _v13_active_log_dir = self.log_dir


def _v13_patched_auto_cal(log_dir, profile_path, min_rows=None):
    if min_rows is None:
        result = _v13_original_auto_cal(log_dir, profile_path)
    else:
        result = _v13_original_auto_cal(log_dir, profile_path, min_rows)
    global _v13_active_profile_path, _v13_active_profile_data
    _v13_active_profile_path = _v13_Path(profile_path)
    try:
        _v13_active_profile_data = result[0] if isinstance(result, tuple) else None
    except Exception:
        _v13_active_profile_data = None
    return result


def _v13_apply_patches():
    globals()['HtmlDashboardBridge'] = _V13EnhancedBridge
    globals()['_HTML_PAGE_TEMPLATE'] = _v13_HTML
    globals()['ui_resolve_location'] = _v13_patched_ui_resolve_location
    globals()['choose_initial_display_mode'] = _v13_patched_choose_mode
    SessionLogger.__init__ = _v13_patched_session_init
    globals()['auto_calibrate_profile'] = _v13_patched_auto_cal


def _v13_startup_ticker():
    stages = [
        (10, 'NICT観測データを待機中…'),
        (25, 'foF2/hmF2 履歴を取得中…'),
        (40, '宇宙天気データを取得中…'),
        (55, 'FT8 実測レポートを取得中…'),
        (70, 'MSTID 画像を解析中…'),
        (85, '初回サイクル完了直前…'),
    ]
    start = _v13_time.time()
    while True:
        if _edfs_stop_event.is_set():
            return
        if globals().get('_ui_last_meta', None) is not None:
            try:
                snap = _v13_build_snapshot()
                with _v13_snapshot_lock:
                    _v13_snapshot.clear()
                    _v13_snapshot.update(snap)
            except Exception:
                pass
        else:
            try:
                elapsed = _v13_time.time() - start
                i = min(int(elapsed // 15), len(stages) - 1)
                pct, msg = stages[i]
                with _v13_cycle_lock:
                    interval_sec = int(_v13_cycle_state.get('interval_sec') or 900)
                    startup_epoch = _v13_cycle_state.get('startup_epoch') or _v13_time.time()
                with _v13_snapshot_lock:
                    _v13_snapshot.update({
                        'ready': False,
                        'startup_progress_pct': pct,
                        'startup_stage': msg,
                        'updated_at': _v13_datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                        'cycle_info': {
                            'interval_sec': interval_sec,
                            'last_cycle_end_epoch': None,
                            'startup_epoch': startup_epoch,
                            'next_cycle_epoch': startup_epoch + interval_sec,
                            'server_epoch': _v13_time.time(),
                        },
                        'map': _v13_build_provisional_map(),
                        'season': {
                            'is_season_off': False,
                            'tier': None,
                            'label': '判定中',
                            'class': 'warn',
                            'recent_pct': None, 'climatology_pct': None, 'blended_pct': None,
                            'consecutive_below': 0,
                            'required': _v13_season_state['off_cycles_required'],
                        },
                    })
            except Exception as e:
                with _v13_snapshot_lock:
                    _v13_snapshot['ready'] = False
                    _v13_snapshot['startup_stage'] = f'起動処理で例外発生(再試行します): {e}'
                    _v13_snapshot['startup_progress_pct'] = max(
                        5, int(_v13_snapshot.get('startup_progress_pct') or 5))
        _v13_time.sleep(3)


def _v13_log_cleanup_ticker(log_dir_str: str):
    """ハイブリッド型ログ管理(cleanup_unbounded_logs())を定期実行するバックグラウンド
    スレッド。起動直後に1回、以降はLOG_CLEANUP_INTERVAL_SEC(既定1時間)ごとに実行する。
    """
    log_dir = _v13_Path(log_dir_str)
    try:
        result = cleanup_unbounded_logs(log_dir)
        if result:
            print(f'[EDFS] ログ自動整理(起動時): {result}')
    except Exception as e:
        print(f'[EDFS] 警告: 起動時ログ整理に失敗 ({e})')
    while True:
        for _ in range(int(LOG_CLEANUP_INTERVAL_SEC)):
            if _edfs_stop_event.is_set():
                return
            _v13_time.sleep(1)
        try:
            result = cleanup_unbounded_logs(log_dir)
            if result:
                print(f'[EDFS] ログ自動整理: {result}')
        except Exception as e:
            print(f'[EDFS] 警告: 定期ログ整理に失敗 ({e})')


def edfs_v13_run_main():
    """Ver13.5 エントリ。 パッチ適用後、原本 main() を呼ぶ。"""
    _v13_apply_patches()

    argv = list(_v13_sys.argv[1:])
    log_dir_arg = str(DEFAULT_DATA_DIR)
    for i, a in enumerate(argv):
        if a == '--interval' and i + 1 < len(argv):
            try:
                with _v13_cycle_lock:
                    _v13_cycle_state['interval_sec'] = int(argv[i + 1])
            except ValueError:
                pass
        elif a.startswith('--interval='):
            try:
                with _v13_cycle_lock:
                    _v13_cycle_state['interval_sec'] = int(a.split('=', 1)[1])
            except ValueError:
                pass
        elif a == '--log-dir' and i + 1 < len(argv):
            log_dir_arg = argv[i + 1]
        elif a.startswith('--log-dir='):
            log_dir_arg = a.split('=', 1)[1]
        elif a == '--tmpfs':
            log_dir_arg = '/tmp'

    with _v13_cycle_lock:
        _v13_cycle_state['startup_epoch'] = _v13_time.time()

    t = _v13_threading.Thread(target=_v13_startup_ticker, daemon=True)
    t.start()

    t_cleanup = _v13_threading.Thread(
        target=_v13_log_cleanup_ticker, args=(log_dir_arg,), daemon=True)
    t_cleanup.start()

    print("=" * 50)
    print(f" {SYSTEM_NAME} {VERSION_LABEL}  (Es & Duct Propagation Forecast Program)")
    print(" All-in-One 単一ファイル版")
    print("=" * 50)
    print("[起動] 以下のサブシステムを、それぞれ完全に独立したバックグラウンド")
    print("       スレッドとして起動します(互いに待ち合わせません):")
    print("[起動]  1/4 Es層予測(電離圏Es層)      … このあとのメインループで開始")
    print("[起動]  2/4 F2層予測(MUF)             … 廃止済み(NICTデータ取得不可のため無効化)")
    print(f"[起動]  3/4 DUCT予測(対流圏ダクト)     … "
          + ("起動します" if '--no-duct' not in argv else "起動しません(--no-duct指定)"))
    print("[起動]  4/4 PSKReporter/WSPR観測・学習 … Es側/DUCT側それぞれ独立スレッド")

    if '--no-duct' not in argv:
        try:
            _v13_threading.Thread(
                target=duct_start_background, args=(_edfs_stop_event,),
                name='duct-bootstrap', daemon=True).start()
            print("[起動]  3/4 DUCT予測: スレッドを起動しました"
                  " (詳細は以降に出る[DUCT]ログ、または/api/duct/statusで確認可能)")
        except Exception as e:
            print(f"[起動]  3/4 DUCT予測: 起動に失敗しました(EDFS本体は継続します): {e}")
    else:
        print("[起動]  3/4 DUCT予測: --no-duct指定のため起動しません")


    if '--no-ui-menu' not in argv:
        argv.append('--no-ui-menu')
    _v13_sys.argv = [_v13_sys.argv[0]] + argv

    main()






"""
ダクト伝播予測プログラム  Ver 6.2.1
==========================================================================
・予測が主目的。自局ウォッチ廃止、全局PSKReporterスポットを面センサー化。
・気象(1時間キャッシュ)とスポット較正(数分)を並行、予測は毎周期再計算。
・HTMLを自動起動しローカルサーバ(127.0.0.1)で配信。設定はUIから変更、
  保存時に即時再計算。地図はIDW連続場+時間スライダー+DPP定量表示。
Ver6.2.1: (実運用フィードバック反映)
  ・[Errno 2]較正ファイル競合クラッシュ&ロストアップデートを修正:
    ライブ観測ループ/WSPRバックフィル/PSK起動時バックフィルの3スレッドが
    dpp_calibration.jsonへ排他制御なしに同時読み書きしていたのが原因
    (コンター消失の実害も確認)。_CAL_IO_LOCKで読込→更新→保存を直列化し、
    WSPRバックフィルはチャンクごとに毎回較正を読み直す方式に変更。
  ・PSKReporterの503連続時に自動バックオフ(15→30→60→…上限300秒)を追加。
    成功すれば即座にウォームアップの通常間隔へ復帰。
  ・WSPR Live統合の可視性向上: 0件でも毎サイクル1行ログ、バックフィルも
    チャンクごとに進捗ログを出すよう変更(「動いているが0件」と「動いて
    いない」を区別できるように)。
  ・「実用まで」のETA推定が「—」に張り付く不具合を修正: 直近24点固定
    (5分間隔想定=約2時間)ではなく直近2時間という実時間窓で回帰するよう
    変更。ウォームアップ中(60秒間隔)は24点=24分しかカバーできず、
    セル減衰の短期ノイズで傾き判定が頻繁に失敗していたため。
  ・SOURCE_VERSION定数を新設し、コンソール起動バナー/生成HTMLのタイトル
    ・バッジ/`/api/status`の3箇所を単一ソースに統一(タイトルだけ旧
    バージョン表記のまま取り残される不具合がここにもあったため修正)。
Ver6.2:
  ・学習速度向上(観測数を実際に増やす方向のみ採用。obs_strength等の
    「見かけ上の学習率」を追加でいじる対応は既にVer6.1で上限到達済みの
    ため見送り):
    ①WSPR Live統合(wspr.live ClickHouse HTTP API, db1.wspr.live): 6m/2m/
      70cm/23cm帯のWSPRスポットを毎サイクルPSKReporterスポットへ追加。
      PSKReporterも一部WSPR局を含むため、(tx_call, rx_call, band,
      5分バケット)キーで重複排除してから ingest_observation() に投入。
      設定 wspr_enabled(既定ON)でOFF可。バンド範囲・地域(region±2°)・
      時間窓でSQL側フィルタしDL量を抑制、20req/分のレート制限を踏まえ
      MIN間隔3.5秒のスロットリングを実装。
    ②WSPR historical backfill: 起動時にバックグラウンドスレッドで
      wspr_backfill_days(既定14日)分の過去WSPRデータを6時間チャンクで
      取得しBayesian較正へ一括投入。進捗はdpp_calibration.jsonの
      wspr_backfillに保存し中断・再起動しても続きから再開。1回完了後は
      done_untilが「現在-数時間」以内なら自動的にスキップ(重複実行防止)。
      設定 wspr_backfill_enabled(既定ON)、wspr_backfill_days で調整可。
    ③PSKReporter起動時バックフィル: 通常サイクルのseconds_back(最大
      3600秒)とは別に、サービス起動直後の1回だけ
      pskreporter_startup_backfill_sec(既定21600秒=6h)で広めに問い合わせ、
      取得できた分だけ較正に反映(API側が拒否/空応答でも通常運転は継続)。
      ※PSKReporter側の実際の遡及可能上限は環境により異なる可能性があり、
      本サンドボックスはpskreporter.info/wspr.liveへのネットワーク到達性
      がないため実機での動作確認が必要。
    ・上記どちらも「観測数を実際に増やす」施策のみ採用し、confidence_n0
      やobs_strength_boostのような閾値操作型の高速化は統計的信頼性を
      損なうため見送り(ユーザーとの相談を踏まえた設計方針)。
Ver6.1:
  ・学習性能加重κ_fm: feature_modelの信頼度を「学習件数」だけでなく
    5-fold CVのAUC/Brier Skill Scoreから算出したperf_factorで重み付け
    (kappa_fm = n_factor × perf_factor)。性能の悪いモデルは件数が
    十分でも自動的に重みが下がる。
  ・feature_logをFIFO上限からReservoir Sampling(Algorithm R)へ変更。
    直近データに偏らず、季節等の長期分布を保ったままサンプルを保持
    (feature_log_seenで全期間の総観測数を管理)。ログに観測時刻を追加。
  ・ウォームアップモード(warmup_mode: auto/on/off)を追加。平均κが
    WARMUP_KAPPA_THRESHOLD(0.35)未満の間はseconds_back/interval_secを
    自動的に最大加速(3600秒窓・60秒周期)し、実用(κ=0.7)到達までの
    時間を短縮。κ_targetに達すれば自動でユーザー設定周期へ復帰する。
    較正はサイクル毎にディスク保存されるため、数時間の連続実行でも
    複数回に分けたセッションでも学習は失われず積算される。
Ver6.0:
  ・学習収束の高速化5案+忘却係数の2段階化を実装。
    ①Surprise Learning: 更新量を|予測-結果|のsurprise_gamma乗で極端化。
    ②可変学習率: 全体平均κに応じてobs_strengthを動的算出(未成熟時は
      最大5倍、成熟でbaseへ収束)。
    ③Adaptive Mesh: 観測密度が閾値を超えたセルのみ動的に0.5°へ細分化
      (split_threshold_n、fine_divisor)。
    ④時間帯別学習: night/morning/afternoon/eveningの4区分でセルを分岐、
      予報の各フレームがその時刻のbucketを参照。
    ⑤特徴量ベース学習: ΔM/厚さ/高度/乱流/結合効率/WKBモード数を観測ごとに
      ログし、numpyのみのロジスティック回帰で「条件→成功率」を学習、
      セル別補正とのブレンドに反映(未観測セルにも汎化可能)。
    ・忘却係数2段階化: 高速忘却(セル単位、half_life_hours=168h)に加え、
      季節・機材ドリフトを吸収する低速忘却(global_drift、約6ヶ月)を追加。
Ver5.2:
  ・🛤️ Pass診断/📶 アンテナ高/🧭 移動運用の座標入力で、緯度経度に加えて
    市町村名などの地名でも入力可能に(Nominatim/OpenStreetMapでジオコー
    ディング、/api/geocode)。「35.68,139.77」のような数値2つが入力欄に
    無ければ自動的に地名として解決してから各APIへ渡す。同一クエリは
    プロセス内メモリでキャッシュし過度な照会を避ける。
Ver5.1:
  ・地図上コンター(等値線)表示: ダクト形成可能性が高いエリア・経路を
    Marching Squares法で等値線化(DPP=30/50/70/85)、UIの📈コンター
    ボタンでON/OFF切替(閾値はCONTOUR_LEVELSで変更可)。
  ・座標入力にGPS自動取得を追加: 🛤️ Pass/📶 アンテナ高/🧭 移動運用の
    各モーダルに「📍 現在地」ボタンを設置、ブラウザのGeolocation APIで
    現在地(緯度,経度)を自動入力(HTTPS/localhost限定の点に注意)。
  ・移動運用レーダーの地図可視化: 🧭 移動運用で取得した地点に★マーカーを
    地図上へ表示し、そこから各方位の到達距離をベクトル(線)として重畳
    表示(色は到達点の実DPP値=既存の凡例スケールに準拠、線をタップで
    方位/距離/DPP値を表示)。
  ・【不具合修正】idw_dpp_at()に最近傍点からの距離上限(max_dist_km)が無く、
    実データの無い遠方まで平坦に外挿されていたため、移動運用レーダーが
    信頼度の低い状況で全方位600km近くまで到達判定してしまう問題を修正。
    field_radar()に外挿距離の上限と、resultsの実データ範囲外へ出たら
    打ち切る境界チェックを追加。
Ver5.0:
  ・帯域に50MHzを追加(BANDS: 50/144/430/1200MHz)。
  ・アラート&自動通知: エリア(半径km)×帯域×DPP閾値の条件成立時、および
    「成長中(+Xhでピーク予測)」検出時にDiscord Webhookへ通知。設定は
    既存のdpp_settings.jsonに統合(新規ファイルなし)、通知クールダウンは
    プロセス内メモリで管理。
  ・特定Pass(2地点間)監視: 大円経路をIDW補間でスコア化し、平均/最小DPP・
    途切れリスクを算出(analyze_pass, /api/pass, UIの🛤️ Passボタン)。
    常時監視したいパスは WATCHED_PASSES に登録するとサイクルごとに自動
    チェック&アラート対象になる。
  ・アンテナ高提案: 既存のcoupling_efficiency()物理式を高さ方向にスイープし、
    結合効率70%到達に必要な高さを提案(antenna_height_curve, /api/antenna,
    UIの📶 アンテナ高ボタン)。
  ・移動運用モード: 24方位への到達距離をレーダー状に算出し、Open-Meteo
    Elevation APIで自局標高とダクト層基底高度の差を算出
    (field_radar, /api/radar, UIの🧭 移動運用ボタン)。
Ver4.0:
  ・Hepburn Tropo Map(DX Info Centre)によるDPP補正を追加(単一ファイルへ統合)。
    1回の画像取得でDPP格子全点を一括サンプリングしログ蓄積、直近14日分で
    OLS回帰(Hepburn強度→DPP)をフィット。サンプル数閾値(30件)未満は
    無補正、到達後はκ=n/(n+40)のベイズ的シュリンケージでブレンド
    (PSKReporter補正と同一思想)。R²を常時算出しレポートパネルに表示。
  ・地図クリックで任意地点のHepburn値をオンデマンド照会するAPI
    (/api/hepburn)を追加。
Ver3.2:
  ・学習成熟度を「信頼度(平均κ)主体」に再定義(被覆定数/200の不整合を廃止)。
    maturity = min(1, 平均κ / KAPPA_TARGET)  (KAPPA_TARGET=0.7)
  ・「実用まで あと何時間」ETA を κ履歴の線形回帰から外挿し UI 表示。
依存: requests, numpy, Pillow(Hepburn補正機能用)。表示は Leaflet(CDN)。
==========================================================================
"""
import os, sys, math, json, time, signal, threading, webbrowser, random
from typing import Optional, Union, Any
from urllib.parse import urlencode
from urllib.request import urlopen, Request
import xml.etree.ElementTree as ET
import numpy as np

try:
    import requests
except ImportError:
    requests = None

try:
    from PIL import Image as _HepImage
except ImportError:
    _HepImage = None

HOME_LAT, HOME_LON = 38.27, 140.87
HEPBURN_COMPARE_ENABLED = True
_HERE = str(DEFAULT_DATA_DIR)
CALIBRATION_FILE = os.path.join(_HERE, "dpp_calibration.json")
OUTPUT_HTML = os.path.join(_HERE, "duct_forecast_map.html")

DUCT_SOURCE_VERSION = "6.2.1"

PRESSURE_LEVELS = [1000, 975, 950, 925, 900, 850, 800, 700, 600, 500]
DEFAULT_REGION = {"lat_min": 26.0, "lat_max": 45.0, "lon_min": 128.0,
                  "lon_max": 145.0, "step_deg": 2.0}
OPEN_METEO_URL = "https://api.open-meteo.com/v1/forecast"
BATCH_SIZE = 40
REQUEST_TIMEOUT = 30
BANDS = {"50MHz": 50.0, "144MHz": 144.0, "351MHz": 351.0, "430MHz": 430.0, "1200MHz": 1296.0}
TX_HEIGHT_M = RX_HEIGHT_M = 10.0

ALERT_BANDS = ["50MHz", "144MHz", "351MHz", "430MHz", "1200MHz"]
ALERT_CENTER_LAT, ALERT_CENTER_LON = HOME_LAT, HOME_LON
ALERT_COOLDOWN_MIN = 60
ALERT_PEAK_LOOKAHEAD = True
ALERT_PASS_BREAK_DPP = 30.0

WATCHED_PASSES = [
]

_ALERT_STATE = {}
_LATEST_RESULTS = []
_RESULTS_LOCK = threading.Lock()

FRAME_HOURS = list(range(0, 25, 3))
REFERENCE_BAND_MHZ = 430.0
TURBULENCE_SHEAR_CAP = 0.05

QUERY_URL = "https://retrieve.pskreporter.info/query"
BAND_INFO = {"6m": (50.0, 500.0), "4m": (70.0, 450.0), "2m": (144.0, 350.0),
             "70cm": (432.0, 250.0), "23cm": (1296.0, 180.0)}
DEFAULT_CONE_DEG = 35.0
EARTH_KM = 6371.0

WSPR_DB_URL = "https://db1.wspr.live/"
WSPR_BAND_CODES = {"6m": 50, "2m": 144, "70cm": 432, "23cm": 1296}
WSPR_QUERY_TIMEOUT = 20
WSPR_MIN_INTERVAL_SEC = 3.5
WSPR_REGION_MARGIN_DEG = 2.0
WSPR_DEDUP_BUCKET_MIN = 5
WSPR_BACKFILL_CHUNK_HOURS = 6
WSPR_BACKFILL_RETRY_SLEEP_SEC = 30
_WSPR_LAST_REQ = {"t": 0.0}
_WSPR_THROTTLE_LOCK = threading.Lock()
_WSPR_BACKFILL_STARTED = threading.Event()

KAPPA_TARGET = 0.7
WARMUP_KAPPA_THRESHOLD = 0.35

_LATEST_DPP_MAP = {}
_LATEST_FEATURE_MAP = {}
_DPP_LOCK = threading.Lock()
_CAL_IO_LOCK = threading.Lock()
WEATHER_CACHE_SEC = 3600
_WX_CACHE = {"ts": 0.0, "key": None, "points": None, "weather": None}
_WX_LOCK = threading.Lock()

SETTINGS_FILE = os.path.join(_HERE, "dpp_settings.json")
SERVER_PORT_RANGE = (8787, 8807)
_SETTINGS = {"interval_sec": 300, "total_minutes": 0, "weather_cache_sec": 3600,
             "seconds_back": 900, "paused": False,
             "alert_enabled": False, "discord_webhook_url": "",
             "alert_dpp_threshold": 70, "alert_radius_km": 300,
             "warmup_mode": "auto",
             "wspr_enabled": True, "wspr_backfill_enabled": True,
             "wspr_backfill_days": 30,
             "pskreporter_startup_backfill_sec": 43200}
_SET_LOCK = threading.Lock()
_RECOMPUTE_NOW = threading.Event()
SETTINGS_CHOICES = {"interval_sec": [60, 120, 180, 300, 600],
                    "total_minutes": [15, 30, 60, 120, 0],
                    "weather_cache_sec": [1800, 3600, 7200],
                    "seconds_back": [600, 900, 1800, 3600],
                    "alert_dpp_threshold": [50, 60, 70, 80, 90],
                    "alert_radius_km": [100, 200, 300, 500, 800],
                    "warmup_mode": ["auto", "on", "off"],
                    "wspr_backfill_days": [3, 7, 14, 30, 60],
                    "pskreporter_startup_backfill_sec":
                        [3600, 10800, 21600, 43200, 86400]}
_BOOL_SETTINGS = {"paused", "alert_enabled", "wspr_enabled", "wspr_backfill_enabled"}
_STRING_SETTINGS = {"discord_webhook_url", "warmup_mode"}


def load_settings():
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
                data = json.load(f)
            with _SET_LOCK:
                for k in _SETTINGS:
                    if k in data:
                        _SETTINGS[k] = data[k]
        except Exception:
            pass
    return get_settings()


def save_settings():
    with _SET_LOCK:
        snap = dict(_SETTINGS)
    try:
        tmp = SETTINGS_FILE + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(snap, f, ensure_ascii=False, indent=1)
        os.replace(tmp, SETTINGS_FILE)
    except Exception as e:
        print(f"[settings] 保存失敗: {e}")


def get_settings():
    with _SET_LOCK:
        return dict(_SETTINGS)


def update_settings(patch, trigger_recompute=True):
    applied = {}
    with _SET_LOCK:
        for k, v in (patch or {}).items():
            if k not in _SETTINGS:
                continue
            if k in _BOOL_SETTINGS:
                _SETTINGS[k] = bool(v); applied[k] = _SETTINGS[k]
            elif k in _STRING_SETTINGS:
                _SETTINGS[k] = str(v)[:300].strip(); applied[k] = _SETTINGS[k]
            else:
                try:
                    _SETTINGS[k] = max(0, int(v)); applied[k] = _SETTINGS[k]
                except (TypeError, ValueError):
                    pass
    save_settings()
    if trigger_recompute and applied and not _SETTINGS.get("paused"):
        request_recompute()
    return applied


def request_recompute():
    _RECOMPUTE_NOW.set()


def _region_key(region):
    return (round(region["lat_min"], 3), round(region["lat_max"], 3),
            round(region["lon_min"], 3), round(region["lon_max"], 3),
            round(region["step_deg"], 3))


def get_weather_cached(region, max_age_sec=WEATHER_CACHE_SEC, force=False):
    key = _region_key(region); now = time.time()
    with _WX_LOCK:
        fresh = (_WX_CACHE["weather"] is not None and _WX_CACHE["key"] == key
                 and (now - _WX_CACHE["ts"]) < max_age_sec)
        if fresh and not force:
            return _WX_CACHE["points"], _WX_CACHE["weather"], True
    points = build_grid(region)
    print(f"[気象] {len(points)}点 取得中(キャッシュ更新)...")
    weather = fetch_weather_all(points)
    with _WX_LOCK:
        _WX_CACHE.update({"ts": now, "key": key, "points": points,
                          "weather": weather})
    return points, weather, False


def saturation_vapor_pressure_hpa(t: float) -> float:
    return 6.1121 * math.exp((18.678 - t / 234.5) * (t / (257.14 + t)))


def refractivity_N(t: float, rh: float, p: float) -> float:
    tk = t + 273.15
    e = max(0.0, rh) / 100.0 * saturation_vapor_pressure_hpa(t)
    return 77.6 * (p / tk) + 3.73e5 * (e / (tk ** 2))


def modified_refractivity_M(n: float, h: float) -> float:
    return n + 0.157 * h


def detect_all_ducts(heights: "np.ndarray | list | None",
                      M_vals: "np.ndarray | list | None",
                      min_thickness_m: float = 10.0) -> list[dict]:
    if heights is None or len(heights) < 3:
        return []
    heights = np.asarray(heights, float); M_vals = np.asarray(M_vals, float)
    dMdz = np.gradient(M_vals, heights); d2 = np.gradient(dMdz, heights)
    mask = dMdz < 0.0
    if not mask.any():
        return []
    idx = np.where(mask)[0]
    groups = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    _trapz = getattr(np, "trapezoid", None) or np.trapz
    ducts = []
    for g in groups:
        if len(g) < 2:
            continue
        b, t = int(g[0]), int(g[-1])
        base_m, top_m = float(heights[b]), float(heights[t])
        thick = top_m - base_m
        if thick < min_thickness_m:
            continue
        seg_M = M_vals[b:t + 1]; seg_h = heights[b:t + 1]
        delta_M = float(M_vals[b] - seg_M.min())
        area = float(_trapz(np.maximum(M_vals[b] - seg_M, 0.0), seg_h))
        deep = b + int(np.argmin(seg_M))
        ducts.append({"base_m": base_m, "top_m": top_m, "thickness_m": thick,
                      "delta_M": delta_M, "min_dMdz": float(dMdz[g].min()),
                      "curvature": float(d2[deep]), "trapped_area": area,
                      "kind": "surface" if base_m <= 50.0 else "elevated"})
    ducts.sort(key=lambda d: d["delta_M"] * d["thickness_m"], reverse=True)
    for i, d in enumerate(ducts):
        d["rank"] = i
    return ducts


def max_trapping_angle_rad(delta_M: float) -> float:
    return math.sqrt(max(0.0, 2e-6 * delta_M))


def coupling_efficiency(antenna_height_m: float, duct: Optional[dict],
                         elevation_angle_deg: float = 0.5) -> float:
    if duct is None:
        return 0.0
    theta_max = max_trapping_angle_rad(duct["delta_M"])
    theta = math.radians(max(0.0, elevation_angle_deg))
    ang = 1.0 if theta <= theta_max else \
        math.exp(-((theta - theta_max) / max(theta_max, 1e-4)) ** 2)
    base, top = duct["base_m"], duct["top_m"]
    if base <= antenna_height_m <= top:
        ht = 1.0
    else:
        dist = (base - antenna_height_m) if antenna_height_m < base \
            else (antenna_height_m - top)
        ht = math.exp(-(dist / max(0.25 * duct["thickness_m"], 20.0)))
    return float(max(0.0, min(1.0, ang * ht)))


def wkb_trapped_modes(heights: np.ndarray, M_vals: np.ndarray,
                       duct: Optional[dict], freq_mhz: float) -> dict:
    if duct is None:
        return {"n_modes": 0, "f_min_mhz": float("inf"), "trapped": False}
    heights = np.asarray(heights, float); M_vals = np.asarray(M_vals, float)
    b = int(np.argmin(np.abs(heights - duct["base_m"])))
    M_c = float(M_vals[b]); seg_h = [heights[b]]; seg_M = [M_vals[b]]
    top_rec = heights[-1]
    for i in range(b + 1, len(heights)):
        seg_h.append(heights[i]); seg_M.append(M_vals[i])
        if M_vals[i] >= M_c:
            h0, h1 = heights[i - 1], heights[i]; m0, m1 = M_vals[i - 1], M_vals[i]
            if m1 != m0:
                top_rec = h0 + (M_c - m0) / (m1 - m0) * (h1 - h0)
            seg_h[-1] = top_rec; seg_M[-1] = M_c
            break
    seg_h = np.asarray(seg_h, float); seg_M = np.asarray(seg_M, float)
    _trapz = getattr(np, "trapezoid", None) or np.trapz
    I = float(_trapz(np.sqrt(np.maximum(2e-6 * (M_c - seg_M), 0.0)), seg_h))
    lam_max = 4.0 * I
    f_min = 300.0 / lam_max if lam_max > 0 else float("inf")
    lam = 300.0 / freq_mhz
    n = max(0, int(math.floor(2.0 * I / lam + 0.5)))
    return {"n_modes": n, "f_min_mhz": f_min, "trapped": freq_mhz >= f_min}


def leakage_db_per_km(freq_mhz: float, wkb: dict) -> Optional[float]:
    if not wkb["trapped"] or wkb["f_min_mhz"] <= 0:
        return None
    ratio = freq_mhz / wkb["f_min_mhz"]
    if ratio <= 1.0:
        return None
    return float(min(0.05 + 2.0 / (ratio - 1.0 + 0.2) ** 2, 20.0))


def refraction_state(dNdz: float) -> dict:
    denom = 1.0 + 6.371e6 * dNdz * 1e-6
    k = 1.0 / denom if abs(denom) > 1e-9 else float("inf")
    if dNdz > 0:
        state = "sub"
    elif dNdz > -0.079:
        state = "normal"
    elif dNdz > -0.157:
        state = "super"
    else:
        state = "trapping"
    k_std = 1.0 / (1.0 + 6.371e6 * (-0.039e-6))
    hg = math.sqrt(max(k, 0.1) / k_std) if k > 0 else 3.0
    return {"k_factor": float(k), "state": state,
            "horizon_gain": float(min(hg, 4.0))}


def radio_horizon_km(h: float, k: float = 1.33) -> float:
    return 3.57 * math.sqrt(max(k, 0.1) * max(h, 0.1))


def duct_haversine_km(lat1: float, lon1: float,
                  lat2: "float | np.ndarray", lon2: "float | np.ndarray"
                  ) -> "float | np.ndarray":
    """2点間の大円距離(km)。lat2/lon2はnumpy配列も可(その場合は要素ごとの距離配列を返す)。"""
    R = EARTH_KM
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(np.subtract(lat2, lat1)); dl = np.radians(np.subtract(lon2, lon1))
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    d = 2 * R * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    return float(d) if np.isscalar(a) or np.ndim(d) == 0 else d


def great_circle_points(lat1: float, lon1: float, lat2: float, lon2: float,
                         n: int = 21) -> list[tuple[float, float]]:
    """2点間の大円上をslerpで等間隔に n 点サンプリングする。"""
    p1 = math.radians(lat1), math.radians(lon1)
    p2 = math.radians(lat2), math.radians(lon2)

    def to_xyz(la, lo):
        return (math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la))
    x1 = to_xyz(*p1); x2 = to_xyz(*p2)
    dot = max(-1.0, min(1.0, sum(a * b for a, b in zip(x1, x2))))
    omega = math.acos(dot)
    if omega < 1e-9:
        return [(lat1, lon1)] * n
    pts = []
    for i in range(n):
        t = i / (n - 1)
        a = math.sin((1 - t) * omega) / math.sin(omega)
        b = math.sin(t * omega) / math.sin(omega)
        v = tuple(a * x1[k] + b * x2[k] for k in range(3))
        lat = math.degrees(math.asin(max(-1.0, min(1.0, v[2]))))
        lon = math.degrees(math.atan2(v[1], v[0]))
        pts.append((lat, lon))
    return pts


def destination_point(lat: float, lon: float, bearing_deg: float,
                       dist_km: float) -> tuple[float, float]:
    """出発点・方位角・距離から到達点(lat,lon)を求める(順測地線問題)。"""
    R = EARTH_KM
    br = math.radians(bearing_deg)
    lat1 = math.radians(lat); lon1 = math.radians(lon)
    dR = dist_km / R
    lat2 = math.asin(math.sin(lat1) * math.cos(dR) + math.cos(lat1) * math.sin(dR) * math.cos(br))
    lon2 = lon1 + math.atan2(math.sin(br) * math.sin(dR) * math.cos(lat1),
                              math.cos(dR) - math.sin(lat1) * math.sin(lat2))
    return math.degrees(lat2), math.degrees(lon2)


def idw_dpp_at(results: list, lat: float, lon: float,
               power: float = 2.0, k: int = 6, min_dist_km: float = 1.0,
               max_dist_km: "float | None" = None
               ) -> "float | None":
    """resultsグリッド(各点frames[0].dpp)から任意地点のDPPを逆距離加重補間する。
    距離計算はnumpyでベクトル化(格子点数が多い場合もPythonループなしで高速)。
    max_dist_km指定時、最近傍点がそれより遠ければNoneを返す(実データの無い
    遠方への平坦な外挿を防ぐガード。無指定時は従来通り制限なし)。"""
    valid = [r for r in results if r and r.get("frames") and r["frames"][0]]
    if not valid:
        return None
    lat_arr = np.array([r["lat"] for r in valid], dtype=float)
    lon_arr = np.array([r["lon"] for r in valid], dtype=float)
    dpp_arr = np.array([r["frames"][0]["dpp"] for r in valid], dtype=float)
    d_arr = duct_haversine_km(lat, lon, lat_arr, lon_arr)
    order = np.argsort(d_arr)[:k]
    d_k = d_arr[order]; v_k = dpp_arr[order]
    if max_dist_km is not None and d_k[0] > max_dist_km:
        return None
    if d_k[0] < min_dist_km:
        return float(v_k[0])
    weights = 1.0 / (d_k ** power)
    return float(np.sum(weights * v_k) / np.sum(weights))


def band_reachability(heights: Optional[np.ndarray], M_vals: Optional[np.ndarray],
                       ducts: list[dict]) -> dict:
    out = {}
    dNdz = 0.0
    if heights is not None and len(heights) >= 2:
        dMdz0 = (M_vals[1] - M_vals[0]) / max(heights[1] - heights[0], 1.0)
        dNdz = dMdz0 - 0.157
    refr = refraction_state(dNdz)
    los = radio_horizon_km(TX_HEIGHT_M, refr["k_factor"]) + \
        radio_horizon_km(RX_HEIGHT_M, refr["k_factor"])
    for band, f_mhz in BANDS.items():
        best = None
        for duct in ducts:
            wkb = wkb_trapped_modes(heights, M_vals, duct, f_mhz)
            coup = coupling_efficiency(TX_HEIGHT_M, duct) * \
                coupling_efficiency(RX_HEIGHT_M, duct)
            leak = leakage_db_per_km(f_mhz, wkb)
            score = (1.0 if wkb["trapped"] else 0.0) * coup * \
                (1.0 / (1.0 + (leak or 20.0)))
            cand = {"trapped": wkb["trapped"], "n_modes": wkb["n_modes"],
                    "f_min_mhz": round(wkb["f_min_mhz"], 1)
                    if math.isfinite(wkb["f_min_mhz"]) else None,
                    "coupling": round(coup, 3), "_score": score}
            if best is None or score > best["_score"]:
                best = cand
        if best is None:
            best = {"trapped": False, "n_modes": 0, "f_min_mhz": None,
                    "coupling": 0.0, "_score": 0.0}
        best["horizon_km"] = round(los * refr["horizon_gain"], 1)
        out[band] = best
    return out, refr


def default_calibration(step_deg=2.0, region=None):
    region = region or DEFAULT_REGION
    return {"schema_version": 3, "model": "bayesian_beta_path_weighted",
            "grid": {"step_deg": step_deg, "lat_min": region["lat_min"],
                     "lat_max": region["lat_max"], "lon_min": region["lon_min"],
                     "lon_max": region["lon_max"]},
            "prior": {"alpha0": 0.6, "beta0": 6.0},
            "half_life_hours": 168.0,
            "half_life_hours_slow": 4380.0,
            "obs_strength_base": 4.0,
            "obs_strength_boost": 5.0,
            "surprise_gamma": 1.8,
            "confidence_n0": 48.0, "bias_clip": [0.5, 1.8],
            "cells": {}, "cone_width_deg": DEFAULT_CONE_DEG, "cone_history": [],
            "kappa_history": [], "history": [],
            "global_drift": {"alpha": 1.0, "beta": 1.0, "last_update": None},
            "time_buckets_enabled": True,
            "adaptive_mesh": {"enabled": True, "split_threshold_n": 60,
                               "fine_divisor": 4},
            "cell_density": {}, "split_coarse_cells": [],
            "feature_log": [], "feature_log_max": 4000, "feature_log_seen": 0,
            "feature_model": None, "feature_model_min_n": 200,
            "wspr_backfill": {"done_until": None, "target_days": None,
                              "completed": False, "last_run": None}}


def load_calibration(path):
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                c = json.load(f)
            if c.get("schema_version") == 3:
                c.setdefault("kappa_history", [])
                c.setdefault("half_life_hours_slow", 4380.0)
                c.setdefault("obs_strength_base", c.get("obs_strength", 4.0))
                c.setdefault("obs_strength_boost", 5.0)
                c.setdefault("surprise_gamma", 1.8)
                c.setdefault("global_drift",
                             {"alpha": 1.0, "beta": 1.0, "last_update": None})
                c.setdefault("time_buckets_enabled", True)
                c.setdefault("adaptive_mesh",
                              {"enabled": True, "split_threshold_n": 60,
                               "fine_divisor": 4})
                c.setdefault("cell_density", {})
                c.setdefault("split_coarse_cells", [])
                c.setdefault("feature_log", [])
                c.setdefault("feature_log_max", 4000)
                c.setdefault("feature_log_seen", len(c.get("feature_log", [])))
                c.setdefault("feature_model", None)
                c.setdefault("feature_model_min_n", 200)
                c.setdefault("wspr_backfill", {"done_until": None, "target_days": None,
                                               "completed": False, "last_run": None})
                if c.get("half_life_hours") in (None, 72.0):
                    c["half_life_hours"] = 168.0
                return c
        except Exception:
            pass
    return default_calibration()


def save_calibration(cal, path):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(cal, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)


def _cell_key(lat, lon, step, bucket=None):
    clat = round(math.floor(lat / step) * step + step / 2.0, 4)
    clon = round(math.floor(lon / step) * step + step / 2.0, 4)
    base = f"{clat},{clon}"
    return f"{base}|{bucket}" if bucket else base


TIME_BUCKETS = [("night", 0, 6), ("morning", 6, 12),
                ("afternoon", 12, 18), ("evening", 18, 24)]


def _time_bucket(dt):
    h = dt.hour
    for name, lo, hi in TIME_BUCKETS:
        if lo <= h < hi:
            return name
    return "night"


def _fine_step(cal):
    am = cal.get("adaptive_mesh", {})
    return cal["grid"]["step_deg"] / max(1, am.get("fine_divisor", 4))


def _mark_split(cal, coarse_key):
    s = cal.setdefault("split_coarse_cells", [])
    if coarse_key not in s:
        s.append(coarse_key)


def _bump_density(cal, coarse_key):
    am = cal.get("adaptive_mesh", {})
    if not am.get("enabled", True):
        return
    d = cal.setdefault("cell_density", {})
    d[coarse_key] = d.get(coarse_key, 0) + 1
    if d[coarse_key] >= am.get("split_threshold_n", 60):
        _mark_split(cal, coarse_key)


def great_circle_km(a, b):
    p1, p2 = math.radians(a[0]), math.radians(b[0])
    dphi = math.radians(b[0] - a[0]); dl = math.radians(b[1] - a[1])
    h = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * EARTH_KM * math.asin(min(1.0, math.sqrt(h)))


def path_cell_weights(tx, rx, cal_or_step, n=200, bucket=None):
    """cal_or_stepにcal(dict)を渡すと適応メッシュ(分割済み地域は自動でfine step)対応。
    従来通りstep(float)を渡した場合は単純な固定gridとして動作(後方互換)。"""
    if isinstance(cal_or_step, dict):
        cal = cal_or_step
        step = cal["grid"]["step_deg"]
        am_enabled = cal.get("adaptive_mesh", {}).get("enabled", True)
        fine_step = _fine_step(cal) if am_enabled else None
        split_set = set(cal.get("split_coarse_cells", [])) if am_enabled else set()
    else:
        step = cal_or_step
        am_enabled = False; fine_step = None; split_set = set()
    per, total = {}, 0.0
    prev = (tx[0], tx[1])
    for i in range(1, n + 1):
        f = i / n
        cur = (tx[0] + (rx[0] - tx[0]) * f, tx[1] + (rx[1] - tx[1]) * f)
        seg = great_circle_km(prev, cur)
        mid = ((prev[0] + cur[0]) / 2, (prev[1] + cur[1]) / 2)
        coarse_base = _cell_key(mid[0], mid[1], step)
        if am_enabled and coarse_base in split_set:
            k = _cell_key(mid[0], mid[1], fine_step, bucket)
        else:
            k = f"{coarse_base}|{bucket}" if bucket else coarse_base
        per[k] = per.get(k, 0.0) + seg; total += seg; prev = cur
    if total <= 0:
        return {}, 0.0
    return {k: v / total for k, v in per.items()}, total


def _coarse_density_wmap(tx, rx, step, n=200):
    """密度カウント専用(bucket/adaptive無視の純粋な粗gridパス通過セル集計)。"""
    per, total = {}, 0.0
    prev = (tx[0], tx[1])
    for i in range(1, n + 1):
        f = i / n
        cur = (tx[0] + (rx[0] - tx[0]) * f, tx[1] + (rx[1] - tx[1]) * f)
        seg = great_circle_km(prev, cur)
        mid = ((prev[0] + cur[0]) / 2, (prev[1] + cur[1]) / 2)
        k = _cell_key(mid[0], mid[1], step)
        per[k] = per.get(k, 0.0) + seg; total += seg; prev = cur
    return per


def dpp_to_success_prob(dpp):
    return 1.0 / (1.0 + math.exp(-0.09 * (dpp - 45.0)))


def _decay_cell(cell, prior, now, half_life_h):
    last = cell.get("last_update")
    if not last:
        return
    try:
        dt_h = (now - datetime.fromisoformat(last)).total_seconds() / 3600.0
    except Exception:
        return
    if dt_h <= 0:
        return
    f = 0.5 ** (dt_h / max(half_life_h, 1.0))
    cell["alpha"] = prior["alpha0"] + (cell["alpha"] - prior["alpha0"]) * f
    cell["beta"] = prior["beta0"] + (cell["beta"] - prior["beta0"]) * f


def dynamic_obs_strength(cal):
    """②学習率可変: 全体成熟度(平均κ)が低いほど大きく、成熟でbaseへ収束。"""
    base = cal.get("obs_strength_base", 4.0)
    boost = cal.get("obs_strength_boost", 5.0)
    k = mean_kappa(cal)
    return base * (boost - (boost - 1.0) * k)


def _update_global_drift(cal, p_pred, y, now):
    """忘却係数(低速側): 季節・機材変更等による全体ドリフトをゆっくり吸収。"""
    half_life_h = cal.get("half_life_hours_slow", 4380.0)
    gd = cal.setdefault("global_drift",
                        {"alpha": 1.0, "beta": 1.0, "last_update": None})
    last = gd.get("last_update")
    if last:
        try:
            dt_h = (now - datetime.fromisoformat(last)).total_seconds() / 3600.0
        except Exception:
            dt_h = 0.0
        if dt_h > 0:
            f = 0.5 ** (dt_h / max(half_life_h, 1.0))
            gd["alpha"] = 1.0 + (gd["alpha"] - 1.0) * f
            gd["beta"] = 1.0 + (gd["beta"] - 1.0) * f
    if y >= 0.5:
        gd["alpha"] += max(0.0, y - p_pred)
    else:
        gd["beta"] += max(0.0, p_pred - y)
    gd["last_update"] = now.isoformat(timespec="seconds")
    cal["global_drift"] = gd


def global_drift_factor(cal):
    gd = cal.get("global_drift")
    if not gd:
        return 1.0
    a, b = gd.get("alpha", 1.0), gd.get("beta", 1.0)
    lo, hi = cal.get("bias_clip", [0.5, 1.8])
    return float(min(hi, max(lo, a / b)))


FEATURE_NAMES = ["delta_M", "thickness_m", "base_m", "turbulence", "coupling", "n_modes"]


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))


def _fit_logreg(Xz, y, epochs=200, lr=0.3, l2=0.01):
    w = np.zeros(Xz.shape[1]); b = 0.0
    n = len(y)
    if n == 0:
        return w, b
    for _ in range(epochs):
        p = _sigmoid(Xz.dot(w) + b)
        grad_w = Xz.T.dot(p - y) / n + l2 * w
        grad_b = float(np.mean(p - y))
        w -= lr * grad_w; b -= lr * grad_b
    return w, b


def _brier(y, p):
    return float(np.mean((np.asarray(p) - np.asarray(y)) ** 2))


def _auc(y, p):
    y = np.asarray(y); p = np.asarray(p)
    pos = p[y == 1]; neg = p[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return 0.5
    order = np.argsort(p, kind="mergesort")
    ranks = np.empty(len(p)); ranks[order] = np.arange(1, len(p) + 1)
    n_pos, n_neg = len(pos), len(neg)
    sum_ranks_pos = ranks[y == 1].sum()
    return float((sum_ranks_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def train_feature_model(cal, epochs=200, lr=0.3, k_folds=5, seed=42):
    """⑤特徴量ベース学習: ロジスティック回帰(numpyのみ)。
    k-fold CVでAUC/Brier Skill Scoreを評価し、性能係数(perf_factor)として
    feature_modelに保存する。件数だけでなく「予測精度」もκ_fmに反映させるため。"""
    log = cal.get("feature_log", [])
    n_min = cal.get("feature_model_min_n", 200)
    if len(log) < n_min:
        return False
    X = np.array([[e["f"].get(k, 0.0) for k in FEATURE_NAMES] for e in log], dtype=float)
    y = np.array([e["y"] for e in log], dtype=float)
    n = len(y)
    if len(np.unique(y)) < 2:
        return False
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    folds = np.array_split(idx, max(2, min(k_folds, n // 30)))
    briers, aucs = [], []
    for i in range(len(folds)):
        val_idx = folds[i]
        train_idx = np.concatenate([folds[j] for j in range(len(folds)) if j != i])
        if len(val_idx) == 0 or len(np.unique(y[train_idx])) < 2:
            continue
        mu_f = X[train_idx].mean(axis=0); sd_f = X[train_idx].std(axis=0)
        sd_f = np.where(sd_f < 1e-6, 1.0, sd_f)
        Xz_tr = (X[train_idx] - mu_f) / sd_f
        Xz_va = (X[val_idx] - mu_f) / sd_f
        w_f, b_f = _fit_logreg(Xz_tr, y[train_idx], epochs=epochs, lr=lr)
        p_va = _sigmoid(Xz_va.dot(w_f) + b_f)
        briers.append(_brier(y[val_idx], p_va)); aucs.append(_auc(y[val_idx], p_va))
    if not briers:
        return False
    cv_brier = float(np.mean(briers)); cv_auc = float(np.mean(aucs))
    p_base = min(max(float(y.mean()), 1e-6), 1.0 - 1e-6)
    brier_baseline = float(np.mean((p_base - y) ** 2))
    bss = 1.0 - cv_brier / brier_baseline if brier_baseline > 1e-9 else 0.0
    auc_score = max(0.0, (cv_auc - 0.5) * 2.0)
    bss_score = max(0.0, min(1.0, bss))
    perf_factor = float(max(0.0, min(1.0, 0.5 * auc_score + 0.5 * bss_score)))
    mu = X.mean(axis=0); sd = X.std(axis=0)
    sd = np.where(sd < 1e-6, 1.0, sd)
    Xz = (X - mu) / sd
    w, b = _fit_logreg(Xz, y, epochs=epochs, lr=lr)
    cal["feature_model"] = {
        "feature_names": FEATURE_NAMES, "mu": mu.tolist(), "sd": sd.tolist(),
        "weights": w.tolist(), "bias": float(b), "p0_baseline": float(y.mean()),
        "n_trained": n,
        "perf": {"cv_auc": round(cv_auc, 4), "cv_brier": round(cv_brier, 4),
                 "brier_skill": round(bss, 4), "perf_factor": round(perf_factor, 4),
                 "n_folds": len(briers)},
        "trained_at": datetime.now().isoformat(timespec="seconds")}
    return True


def feature_model_correction(cal, features):
    """学習済み特徴量モデルによる補正係数とその信頼度κ_fmを返す。未学習なら中立(1.0,0.0)。
    κ_fm = 件数係数 × 性能係数(CV AUC/Brier Skill Scoreベース)。
    性能が悪いモデル(perf_factor≈0)はいくら件数があっても実質重み0になる。"""
    fm = cal.get("feature_model")
    if not fm or not features:
        return 1.0, 0.0
    names = fm["feature_names"]
    x = np.array([features.get(k, 0.0) for k in names], dtype=float)
    mu = np.array(fm["mu"]); sd = np.array(fm["sd"])
    sd = np.where(sd < 1e-6, 1.0, sd)
    xz = (x - mu) / sd
    z = float(np.dot(xz, np.array(fm["weights"])) + fm["bias"])
    p_hat = 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z))))
    p0 = max(fm.get("p0_baseline", 0.5), 1e-3)
    lo, hi = cal.get("bias_clip", [0.5, 1.8])
    factor = min(hi, max(lo, p_hat / p0))
    n0 = cal.get("confidence_n0", 48.0) * 2.0
    n_trained = fm.get("n_trained", 0)
    n_factor = n_trained / (n_trained + n0)
    perf_factor = fm.get("perf", {}).get("perf_factor", 0.0)
    kappa_fm = n_factor * perf_factor
    return factor, float(kappa_fm)


def _reservoir_add(cal, entry):
    """Reservoir Sampling(Algorithm R): FIFOだと直近偏り(例: 夏の観測で埋まり冬を忘れる)が
    起きるため、全期間から一様ランダムにサンプルを保持し季節分布を維持する。"""
    log = cal.setdefault("feature_log", [])
    maxlen = cal.get("feature_log_max", 4000)
    seen = cal.get("feature_log_seen", 0) + 1
    cal["feature_log_seen"] = seen
    if len(log) < maxlen:
        log.append(entry)
    else:
        j = random.randint(1, seen)
        if j <= maxlen:
            log[j - 1] = entry


def ingest_observation(cal, tx, rx, received, freq_mhz, snr_db, dpp_lookup,
                       timestamp=None, feature_lookup=None):
    prior = cal["prior"]; step = cal["grid"]["step_deg"]
    now = timestamp or datetime.now()
    bucket = _time_bucket(now) if cal.get("time_buckets_enabled", True) else None
    for cbase in _coarse_density_wmap(tx, rx, step):
        _bump_density(cal, cbase)
    wmap, dist = path_cell_weights(tx, rx, cal, bucket=bucket)
    if not wmap:
        return 0
    for c, w in wmap.items():
        if w < 0.05:
            continue
        cell = cal["cells"].get(c, {"alpha": prior["alpha0"],
                                    "beta": prior["beta0"], "last_update": None,
                                    "n_obs": 0})
        cell["n_obs"] = cell.get("n_obs", 0) + 1
        cal["cells"][c] = cell
    dpps = [(dpp_lookup(c.split("|")[0]), w) for c, w in wmap.items()]
    dpps = [(d, w) for d, w in dpps if d is not None]
    dpp_bar = sum(d * w for d, w in dpps) / sum(w for _, w in dpps) if dpps else 0.0
    p_pred = dpp_to_success_prob(dpp_bar)
    y = 1.0 if received else 0.0
    snr_term = min(1.0, max(0.0, (snr_db + 5.0) / 25.0)) if received else 1.0
    surprise = abs(y - p_pred) ** cal.get("surprise_gamma", 1.0)
    k_obs = dynamic_obs_strength(cal) * surprise * (0.5 + 0.5 * snr_term)
    _update_global_drift(cal, p_pred, y, now)
    if feature_lookup is not None and wmap:
        primary_c = max(wmap, key=wmap.get).split("|")[0]
        feat = feature_lookup(primary_c)
        if feat:
            _reservoir_add(cal, {"f": feat, "y": y,
                                 "t": now.isoformat(timespec="seconds")})
            if cal.get("feature_log_seen", 0) % 25 == 0:
                train_feature_model(cal)
    if k_obs <= 1e-3:
        return 0
    updated = 0
    for c, w in wmap.items():
        cell = cal["cells"].get(c, {"alpha": prior["alpha0"],
                                    "beta": prior["beta0"], "last_update": None,
                                    "n_obs": 0})
        _decay_cell(cell, prior, now, cal["half_life_hours"])
        if received:
            cell["alpha"] += k_obs * w
        else:
            cell["beta"] += k_obs * w
        cell["last_update"] = now.isoformat(timespec="seconds")
        cal["cells"][c] = cell; updated += 1
    return updated


def _posterior(cell, prior):
    a = cell.get("alpha", prior["alpha0"]); b = cell.get("beta", prior["beta0"])
    return a / (a + b), (a + b) - (prior["alpha0"] + prior["beta0"])


def _lookup_cell(cal, lat, lon, bucket=None):
    """③適応メッシュ対応セル参照: 分割済み地域はfineキーを優先し、無ければ粗キーへfallback。"""
    step = cal["grid"]["step_deg"]
    coarse_base = _cell_key(lat, lon, step)
    am = cal.get("adaptive_mesh", {})
    if am.get("enabled", True) and coarse_base in set(cal.get("split_coarse_cells", [])):
        fine_key = _cell_key(lat, lon, _fine_step(cal), bucket)
        fcell = cal["cells"].get(fine_key)
        if fcell:
            return fcell
    coarse_key = f"{coarse_base}|{bucket}" if bucket else coarse_base
    return cal["cells"].get(coarse_key)


def apply_correction(cal, lat, lon, bucket=None, features=None):
    prior = cal["prior"]; p0 = prior["alpha0"] / (prior["alpha0"] + prior["beta0"])
    lo, hi = cal["bias_clip"]
    cell = _lookup_cell(cal, lat, lon, bucket)
    local, kappa = 1.0, 0.0
    if cell:
        _decay_cell(cell, prior, datetime.now(), cal["half_life_hours"])
        mean, n_eff = _posterior(cell, prior)
        kappa = n_eff / (n_eff + cal.get("confidence_n0", 48.0))
        bias = min(hi, max(lo, mean / p0))
        local = (1.0 - kappa) * 1.0 + kappa * bias
    fm_factor, fm_kappa = feature_model_correction(cal, features)
    w_local, w_fm = kappa, fm_kappa
    w_sum = w_local + w_fm
    if w_sum > 1e-6:
        blended = (w_local * local + w_fm * fm_factor + max(0.0, 1.0 - w_sum) * 1.0)
    else:
        blended = 1.0
    g = global_drift_factor(cal)
    final = min(hi, max(lo, blended * g))
    return float(final), float(max(kappa, fm_kappa))


def mean_kappa(cal):
    cells = cal.get("cells", {})
    if not cells:
        return 0.0
    prior = cal["prior"]; n0 = cal.get("confidence_n0", 48.0)
    ks = []
    for cell in cells.values():
        _, n_eff = _posterior(cell, prior)
        ks.append(n_eff / (n_eff + n0))
    return sum(ks) / len(ks)


ETA_WINDOW_HOURS = 2.0
ETA_MIN_SPAN_HOURS = 0.25


def record_kappa(cal, max_hist=288):
    """平均κの時系列を記録(ETA外挿用)。~24h @5分。ウォームアップ中は60秒周期で
    呼ばれるため実際にはより短い時間しかカバーしないが、estimate_eta_hours側で
    実時間窓による抽出に切り替えているため問題ない。"""
    hist = cal.setdefault("kappa_history", [])
    hist.append({"t": datetime.now().isoformat(timespec="seconds"),
                 "k": round(mean_kappa(cal), 4)})
    cal["kappa_history"] = hist[-max_hist:]


def estimate_eta_hours(cal, cur_kappa, target=KAPPA_TARGET):
    """
    κ履歴の線形回帰から κ=target 到達までの時間[h]を外挿。
    戻り値: 0.0(到達済) / 数値[h] / None(データ不足or増加傾向なし)
    Ver6.2修正: 直近"24点"ではなく直近"ETA_WINDOW_HOURS時間分"を対象にする。
    ウォームアップ中(60秒周期)は24点だと実時間わずか~24分しかカバーできず、
    セル減衰による短周期ノイズで傾きがすぐ0以下になり「—」に張り付いていたため。"""
    if cur_kappa >= target:
        return 0.0
    hist = cal.get("kappa_history", [])
    if len(hist) < 3:
        return None
    now = datetime.now()
    pts = hist[-24:]
    try:
        recent = [p for p in hist if
                  (now - datetime.fromisoformat(p["t"])).total_seconds()
                  <= ETA_WINDOW_HOURS * 3600]
        if len(recent) >= 3:
            pts = recent
    except Exception:
        pass
    try:
        t0 = datetime.fromisoformat(pts[0]["t"])
        xs = [(datetime.fromisoformat(p["t"]) - t0).total_seconds() / 3600.0
              for p in pts]
        ys = [p["k"] for p in pts]
    except Exception:
        return None
    span_hours = xs[-1] - xs[0] if len(xs) >= 2 else 0.0
    if span_hours < ETA_MIN_SPAN_HOURS:
        return None
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-9:
        return None
    slope = (n * sxy - sx * sy) / denom
    if slope <= 1e-3:
        return None
    return round((target - cur_kappa) / slope, 1)


def _warmup_active(cal):
    """ウォームアップモード判定: 平均κが閾値未満の間はデータ取得を加速する。"""
    mode = get_settings().get("warmup_mode", "auto")
    if mode == "off":
        return False
    if mode == "on":
        return True
    return mean_kappa(cal) < WARMUP_KAPPA_THRESHOLD


def _effective_fetch_params(cal, s):
    """ウォームアップ中はseconds_back/interval_secをUI選択肢の最大加速値まで自動的に強める。
    実用(κ=KAPPA_TARGET)までの時間を、数時間〜複数セッションでの積算により短縮する狙い。
    平常時(κ十分)はユーザー設定どおりに戻る。較正はサイクル毎にディスク保存されるため、
    セッションをまたいでも学習は失われず積算される。"""
    warm = _warmup_active(cal)
    if warm:
        eff_seconds_back = max(s["seconds_back"], max(SETTINGS_CHOICES["seconds_back"]))
        eff_interval = min(s["interval_sec"], min(SETTINGS_CHOICES["interval_sec"]))
    else:
        eff_seconds_back = s["seconds_back"]; eff_interval = s["interval_sec"]
    return eff_seconds_back, eff_interval, warm


def learning_quality(cal):
    """信頼度(平均κ)主体の成熟度 + ETA。被覆は補助表示。"""
    cells = cal.get("cells", {})
    warm = _warmup_active(cal)
    if not cells:
        return {"covered_cells": 0, "mean_confidence": 0.0,
                "maturity_pct": 0.0, "eta_hours": None,
                "kappa_target": KAPPA_TARGET, "warmup_active": warm}
    mk = mean_kappa(cal)
    maturity = min(1.0, mk / KAPPA_TARGET)
    eta = estimate_eta_hours(cal, mk, KAPPA_TARGET)
    return {"covered_cells": len(cells), "mean_confidence": round(mk, 3),
            "maturity_pct": round(maturity * 100, 1), "eta_hours": eta,
            "kappa_target": KAPPA_TARGET, "warmup_active": warm}


def duct_maidenhead_to_latlon(loc):
    if not loc or len(loc.strip()) < 4:
        return None
    loc = loc.strip().upper()
    try:
        lon = (ord(loc[0]) - 65) * 20.0 - 180.0
        lat = (ord(loc[1]) - 65) * 10.0 - 90.0
        lon += int(loc[2]) * 2.0; lat += int(loc[3]) * 1.0
        if len(loc) >= 6:
            lon += (ord(loc[4]) - 65) * (2.0 / 24.0) + (1.0 / 24.0)
            lat += (ord(loc[5]) - 65) * (1.0 / 24.0) + (0.5 / 24.0)
        else:
            lon += 1.0; lat += 0.5
        if -90 <= lat <= 90 and -180 <= lon <= 180:
            return (round(lat, 4), round(lon, 4))
    except (ValueError, IndexError):
        return None
    return None


def _band_from_freq(f):
    for name, (fc, _) in BAND_INFO.items():
        if abs(f - fc) / fc < 0.15:
            return name
    return "HF" if f < 30 else None


def normalize_spot(raw):
    def g(*ks):
        for k in ks:
            if k in raw and raw[k] not in (None, ""):
                return raw[k]
        return None
    tx = duct_maidenhead_to_latlon(g("sl", "senderLocator"))
    rx = duct_maidenhead_to_latlon(g("rl", "receiverLocator"))
    if tx is None or rx is None:
        return None
    try:
        freq = float(g("f", "frequency"))
        f_mhz = freq / 1e6 if freq > 1e5 else freq
    except (TypeError, ValueError):
        return None
    try:
        snr = float(g("rp", "sNR", "snr") or 0.0)
    except (TypeError, ValueError):
        snr = 0.0
    t = g("t", "flowStartSeconds")
    try:
        ts = (datetime.fromtimestamp(int(t), timezone.utc)
              .replace(tzinfo=None)) if t else datetime.now()
    except (TypeError, ValueError, OverflowError):
        ts = datetime.now()
    return {"tx_call": g("sc", "senderCallsign"),
            "rx_call": g("rc", "receiverCallsign"), "tx": tx, "rx": rx,
            "snr_db": snr, "freq_mhz": f_mhz,
            "band": g("b") or _band_from_freq(f_mhz), "time": ts}


def _bearing_deg(a, b):
    lat1, lat2 = math.radians(a[0]), math.radians(b[0])
    dlon = math.radians(b[1] - a[1])
    y = math.sin(dlon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    return (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def _ang_diff(a, b):
    d = abs(a - b) % 360.0
    return d if d <= 180 else 360 - d


def in_region(pt, region, margin=1.0):
    return (region["lat_min"] - margin <= pt[0] <= region["lat_max"] + margin and
            region["lon_min"] - margin <= pt[1] <= region["lon_max"] + margin)


def is_anomalous(spot):
    band = spot.get("band")
    if band not in BAND_INFO:
        return False
    return great_circle_km(spot["tx"], spot["rx"]) > BAND_INFO[band][1]


def filter_spots(spots, region):
    return [s for s in spots
            if (in_region(s["tx"], region) or in_region(s["rx"], region))
            and is_anomalous(s)]


def _wspr_utcnow():
    """naive UTC の現在時刻(normalize_spot()のPSKReporter側と表現を揃える)。
    datetime.utcnow()はPython 3.12+で非推奨のためここに集約。"""
    return datetime.now(timezone.utc).replace(tzinfo=None)


def _wspr_throttle():
    """wspr.liveのレート制限(20req/分)対策。呼び出し間隔を最低限確保する。"""
    with _WSPR_THROTTLE_LOCK:
        now = time.monotonic()
        wait = WSPR_MIN_INTERVAL_SEC - (now - _WSPR_LAST_REQ["t"])
        if wait > 0:
            time.sleep(wait)
        _WSPR_LAST_REQ["t"] = time.monotonic()


def wspr_live_query(sql, timeout=WSPR_QUERY_TIMEOUT):
    """wspr.liveのClickHouse HTTPインターフェースへSELECTを投げ、行のlist[dict]を返す。
    失敗時は空リスト(呼び出し側は通常運転を継続する)。"""
    if requests is None:
        return []
    _wspr_throttle()
    try:
        r = requests.get(WSPR_DB_URL, params={"query": sql + " FORMAT JSON"},
                         timeout=timeout)
        r.raise_for_status()
        data = r.json()
        return data.get("data", [])
    except Exception as e:
        print(f"[wspr] クエリ失敗: {e}")
        return []


def _wspr_row_to_spot(row, band_name):
    try:
        tx = (float(row["tx_lat"]), float(row["tx_lon"]))
        rx = (float(row["rx_lat"]), float(row["rx_lon"]))
        freq_hz = float(row["frequency"])
        ts = datetime.strptime(row["time"][:19], "%Y-%m-%d %H:%M:%S")
    except (KeyError, ValueError, TypeError):
        return None
    if tx == (0.0, 0.0) or rx == (0.0, 0.0):
        return None
    return {"tx_call": row.get("tx_sign"), "rx_call": row.get("rx_sign"),
            "tx": tx, "rx": rx, "snr_db": float(row.get("snr") or 0.0),
            "freq_mhz": freq_hz / 1e6, "band": band_name, "time": ts,
            "source": "wspr"}


def wspr_fetch_window(region, start_dt, end_dt, band_codes=None, max_rows=100000):
    """指定期間・地域(±WSPR_REGION_MARGIN_DEG)のWSPRスポットを取得。
    start_dt/end_dtはnaive UTC(normalize_spot()のPSKReporter側と揃える)。"""
    band_codes = band_codes or WSPR_BAND_CODES
    lat_min = region["lat_min"] - WSPR_REGION_MARGIN_DEG
    lat_max = region["lat_max"] + WSPR_REGION_MARGIN_DEG
    lon_min = region["lon_min"] - WSPR_REGION_MARGIN_DEG
    lon_max = region["lon_max"] + WSPR_REGION_MARGIN_DEG
    out = []
    for band_name, code in band_codes.items():
        sql = (
            "SELECT time, tx_sign, tx_lat, tx_lon, rx_sign, rx_lat, rx_lon, "
            "frequency, snr FROM wspr.rx WHERE band = {code} "
            "AND time >= toDateTime('{s}') AND time < toDateTime('{e}') "
            "AND ((tx_lat BETWEEN {a} AND {b} AND tx_lon BETWEEN {c} AND {d}) "
            "OR (rx_lat BETWEEN {a} AND {b} AND rx_lon BETWEEN {c} AND {d})) "
            "LIMIT {lim}"
        ).format(code=code, s=start_dt.strftime("%Y-%m-%d %H:%M:%S"),
                 e=end_dt.strftime("%Y-%m-%d %H:%M:%S"),
                 a=lat_min, b=lat_max, c=lon_min, d=lon_max, lim=max_rows)
        for row in wspr_live_query(sql):
            sp = _wspr_row_to_spot(row, band_name)
            if sp:
                out.append(sp)
    return out


def _spot_dedup_key(s, bucket_min=WSPR_DEDUP_BUCKET_MIN):
    """(送信局, 受信局, バンド, 時刻バケット)でスポットを同一視するためのキー。
    PSKReporterはWSPR局からの報告も含み得るため、WSPR Live側との重複投入を防ぐ。"""
    t = s.get("time")
    tb = int(t.timestamp() // (bucket_min * 60)) if isinstance(t, datetime) else 0
    return (s.get("tx_call"), s.get("rx_call"), s.get("band"), tb)


def merge_wspr_live(spots, region, seconds_back):
    """毎サイクル、直近seconds_back分のWSPR Liveスポットをspotsへ重複排除して追加する。
    件数0の場合も含めて必ず1行ログを出す(「動いているが0件」と「動いていない」を
    区別できるようにするため)。"""
    if not get_settings().get("wspr_enabled", True):
        return spots
    try:
        now = _wspr_utcnow()
        start = now - timedelta(seconds=int(seconds_back))
        wspots = wspr_fetch_window(region, start, now)
    except Exception as e:
        print(f"[wspr] ライブ取得失敗(通常運転は継続): {e}")
        return spots
    raw_n = len(wspots)
    wspots = filter_spots(wspots, region) if wspots else []
    existing_keys = {_spot_dedup_key(s) for s in spots}
    added = [s for s in wspots if _spot_dedup_key(s) not in existing_keys]
    print(f"[wspr] 確認: 生取得{raw_n}件 → 地域内異常{len(wspots)}件 → "
          f"新規追加{len(added)}件(重複{len(wspots) - len(added)}件)")
    return spots + added


def wspr_backfill_worker(cal_path, region, stop_event, days=None):
    """起動時バックグラウンドスレッド: 過去days日分のWSPRデータを
    WSPR_BACKFILL_CHUNK_HOURS時間刻みで取得し較正へ一括投入する。
    中断されても dpp_calibration.json の wspr_backfill.done_until から再開できる。
    Ver6.2追記: 較正ファイルはライブ観測ループ(_spots_daemon)も並行して読み書きするため、
    このワーカーはチャンクごとに cal を都度再読込→更新→保存する(_CAL_IO_LOCKで直列化)。
    起動時に一度だけ cal を読んでチャンクの間ずっと使い回すと、その間にライブループが
    書き込んだ内容をワーカーの最終保存が丸ごと上書きして消してしまう(ロストアップデート)。"""
    if _WSPR_BACKFILL_STARTED.is_set():
        return
    _WSPR_BACKFILL_STARTED.set()
    try:
        days = days if days is not None else get_settings().get("wspr_backfill_days", 14)
        now = _wspr_utcnow()
        target_start = now - timedelta(days=int(days))
        with _CAL_IO_LOCK:
            cal0 = load_calibration(cal_path)
        wb0 = cal0.get("wspr_backfill") or {}
        if wb0.get("completed") and wb0.get("target_days") == days and wb0.get("done_until"):
            try:
                done_until = datetime.fromisoformat(wb0["done_until"])
                if (now - done_until).total_seconds() < 6 * 3600:
                    print("[wspr] バックフィルは直近に完了済みのためスキップ")
                    return
            except (ValueError, TypeError):
                pass
        cursor = target_start
        if wb0.get("done_until") and wb0.get("target_days") == days:
            try:
                saved = datetime.fromisoformat(wb0["done_until"])
                if saved > cursor:
                    cursor = saved
            except (ValueError, TypeError):
                pass
        total_hours = max(1, int((now - cursor).total_seconds() // 3600))
        print(f"[wspr] バックフィル開始: 過去{days}日分, 残り約{total_hours}時間分を"
              f"{WSPR_BACKFILL_CHUNK_HOURS}時間刻みで取得します(バックグラウンド)")
        n_total = 0
        n_chunks_total = max(1, -(-total_hours // WSPR_BACKFILL_CHUNK_HOURS))
        chunk_i = 0
        while cursor < now and not stop_event.is_set():
            if not get_settings().get("wspr_backfill_enabled", True):
                print("[wspr] バックフィルが設定でOFFにされたため中断"); return
            chunk_end = min(now, cursor + timedelta(hours=WSPR_BACKFILL_CHUNK_HOURS))
            chunk_i += 1
            try:
                wspots = wspr_fetch_window(region, cursor, chunk_end)
                wspots = filter_spots(wspots, region)
                with _DPP_LOCK:
                    dpp_map = dict(_LATEST_DPP_MAP)
                    feat_map = dict(_LATEST_FEATURE_MAP)
                look = (lambda c: dpp_map.get(c)) if dpp_map else (lambda c: None)
                flook = (lambda c: feat_map.get(c)) if feat_map else (lambda c: None)
                with _CAL_IO_LOCK:
                    cal = load_calibration(cal_path)
                    n_ingested = 0
                    for s in wspots:
                        if ingest_observation(cal, s["tx"], s["rx"], True, s["freq_mhz"],
                                              s["snr_db"], look, timestamp=s.get("time"),
                                              feature_lookup=flook):
                            n_ingested += 1
                    wb = cal.setdefault(
                        "wspr_backfill",
                        {"done_until": None, "target_days": None,
                         "completed": False, "last_run": None})
                    wb["target_days"] = days
                    wb["completed"] = False
                    wb["done_until"] = chunk_end.isoformat(timespec="seconds")
                    wb["last_run"] = now.isoformat(timespec="seconds")
                    save_calibration(cal, cal_path)
                n_total += n_ingested
                print(f"[wspr] backfill進捗 {chunk_i}/{n_chunks_total}: "
                      f"{cursor.strftime('%m/%d %H:%M')}〜{chunk_end.strftime('%m/%d %H:%M')} "
                      f"取得{len(wspots)}件 投入{n_ingested}件 (累計{n_total}件)")
            except Exception as e:
                print(f"[wspr] バックフィルchunk失敗({cursor}〜{chunk_end}): {e}"
                      f" → {WSPR_BACKFILL_RETRY_SLEEP_SEC}秒後リトライ")
                if stop_event.wait(WSPR_BACKFILL_RETRY_SLEEP_SEC):
                    return
                continue
            cursor = chunk_end
        if not stop_event.is_set():
            with _CAL_IO_LOCK:
                cal = load_calibration(cal_path)
                wb = cal.setdefault(
                    "wspr_backfill",
                    {"done_until": None, "target_days": None,
                     "completed": False, "last_run": None})
                wb["completed"] = True
                save_calibration(cal, cal_path)
            print(f"[wspr] バックフィル完了: 累計{n_total}件を投入しました")
    finally:
        _WSPR_BACKFILL_STARTED.clear()


def pskreporter_startup_backfill(cal_path, region):
    """起動直後の1回だけ、通常サイクルより広いseconds_backでPSKReporterへ問い合わせ、
    取得できた分だけ較正へ反映する。API側が広い窓を拒否/空応答でも通常運転は継続する。
    ※PSKReporter側の実際の遡及可能上限は未確認(本サンドボックスは接続不可のため)。"""
    sb = get_settings().get("pskreporter_startup_backfill_sec", 21600)
    if sb <= 0:
        return
    print(f"[psk] 起動時バックフィル試行: seconds_back={sb}秒")
    try:
        with _CAL_IO_LOCK:
            cal = load_calibration(cal_path)
            st = spots_cycle(cal, region, sb)
            save_calibration(cal, cal_path)
        print(f"[psk] 起動時バックフィル完了: +{st['pos']}/-{st['neg']} "
              f"(取得{st['spots']}件)")
    except Exception as e:
        print(f"[psk] 起動時バックフィル失敗(通常運転は継続): {e}")


def estimate_cone_geometric(spots, region, min_samples=8):
    devs = []; groups = {}
    for s in spots:
        if s["band"] not in BAND_INFO:
            continue
        if not (in_region(s["tx"], region) or in_region(s["rx"], region)):
            continue
        groups.setdefault((s["rx_call"], s["band"]), []).append(
            _bearing_deg(s["rx"], s["tx"]))
    for _, brgs in groups.items():
        if len(brgs) < 3:
            continue
        ang = np.radians(brgs)
        axis = math.degrees(math.atan2(np.sin(ang).mean(), np.cos(ang).mean()))
        devs += [_ang_diff(b, axis) for b in brgs]
    if len(devs) < min_samples:
        return None
    return float(np.percentile(np.asarray(devs, float), 90.0))


def estimate_cone_stratified(spots, region, min_bins=3):
    bins = {}
    for s in spots:
        h = s["time"].hour if isinstance(s.get("time"), datetime) else 0
        bins.setdefault(h, []).append(s)
    est = [e for e in (estimate_cone_geometric(g, region, 6)
                       for g in bins.values()) if e is not None]
    return float(np.median(est)) if len(est) >= min_bins else None


def update_cone_history(cal, new_est, max_hist=60):
    if new_est is None:
        return cal.get("cone_width_deg", DEFAULT_CONE_DEG)
    hist = cal.setdefault("cone_history", [])
    hist.append({"t": datetime.now().isoformat(timespec="seconds"),
                 "est": round(new_est, 2)})
    cal["cone_history"] = hist[-max_hist:]
    cal["cone_width_deg"] = round(float(np.median(
        [h["est"] for h in cal["cone_history"]])), 2)
    return cal["cone_width_deg"]


def build_negatives(spots, region, cone_deg, max_per_rx=3, max_negatives=1500):
    heard_by_rx, active_tx = {}, {}
    for s in spots:
        if s["band"] not in BAND_INFO:
            continue
        rec = heard_by_rx.setdefault((s["rx_call"], s["band"]),
                                     {"pt": s["rx"], "heard": {}})
        rec["heard"][s["tx_call"]] = (great_circle_km(s["tx"], s["rx"]),
                                      _bearing_deg(s["rx"], s["tx"]))
        active_tx[(s["tx_call"], s["band"])] = s["tx"]
    negatives = []
    for (rc, band), rec in heard_by_rx.items():
        if not in_region(rec["pt"], region) or not rec["heard"]:
            continue
        d_far = max(v[0] for v in rec["heard"].values())
        _, los = BAND_INFO[band]; n_added = 0
        for (tc, tb), tx_pt in active_tx.items():
            if tb != band or tc == rc or tc in rec["heard"]:
                continue
            dist = great_circle_km(tx_pt, rec["pt"])
            if dist <= los or dist >= d_far:
                continue
            brg = _bearing_deg(rec["pt"], tx_pt)
            if not any(_ang_diff(brg, v[1]) <= cone_deg and v[0] >= dist
                       for v in rec["heard"].values()):
                continue
            negatives.append({"tx": tx_pt, "rx": rec["pt"],
                              "freq_mhz": BAND_INFO[band][0]})
            n_added += 1
            if n_added >= max_per_rx:
                break
        if len(negatives) >= max_negatives:
            break
    return negatives


_PSK_BACKOFF = {"fail_streak": 0}


def fetch_query_http(seconds_back=900, timeout=REQUEST_TIMEOUT):
    if requests is None:
        return [], []
    params = {"flowStartSeconds": -abs(int(seconds_back)), "rronly": 1,
              "nolocator": 0, "appcontact": "ROA-DPP",
              "frange": "50000000-1300000000"}
    try:
        r = requests.get(QUERY_URL, params=params, timeout=timeout)
        r.raise_for_status()
        _PSK_BACKOFF["fail_streak"] = 0
        return parse_query_xml(r.text)
    except Exception as e:
        _PSK_BACKOFF["fail_streak"] = min(_PSK_BACKOFF["fail_streak"] + 1, 8)
        print(f"[spots] 取得失敗: {e}")
        return [], []


def psk_backoff_extra_sec():
    """PSKReporterが503等を連発している間、ウォームアップの60秒固定間隔に
    追加の待機を上乗せしてサーバー負荷を下げる(指数バックオフ、上限300秒)。
    1回成功すればfail_streakが0に戻り、通常のウォームアップ間隔へ即復帰する。"""
    n = _PSK_BACKOFF["fail_streak"]
    return 0 if n <= 0 else min(300, 15 * (2 ** (n - 1)))


def parse_query_xml(xml_text):
    spots, active = [], []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"[spots] XMLパース失敗: {e}")
        return [], []
    for el in root.iter():
        tag = el.tag.split("}")[-1]; a = el.attrib
        if tag == "receptionReport":
            s = normalize_spot({
                "senderCallsign": a.get("senderCallsign"),
                "senderLocator": a.get("senderLocator"),
                "receiverCallsign": a.get("receiverCallsign"),
                "receiverLocator": a.get("receiverLocator"),
                "sNR": a.get("sNR"), "frequency": a.get("frequency"),
                "flowStartSeconds": a.get("flowStartSeconds")})
            if s:
                spots.append(s)
        elif tag == "activeReceiver":
            pt = duct_maidenhead_to_latlon(a.get("locator"))
            if pt:
                active.append((a.get("callsign"), pt))
    return spots, active


def spots_cycle(cal, region, seconds_back, spots_fetch=None):
    if spots_fetch is not None:
        raw, _active = spots_fetch(seconds_back)
        spots = [r if (isinstance(r, dict) and "tx" in r) else normalize_spot(r)
                 for r in raw]
        spots = [s for s in spots if s]
    else:
        spots, _active = fetch_query_http(seconds_back)
    spots = filter_spots(spots, region)
    if spots_fetch is None:
        spots = merge_wspr_live(spots, region, seconds_back)
    est = estimate_cone_stratified(spots, region) or \
        estimate_cone_geometric(spots, region)
    cone = update_cone_history(cal, est)
    negatives = build_negatives(spots, region, cone)
    with _DPP_LOCK:
        dpp_map = dict(_LATEST_DPP_MAP)
        feat_map = dict(_LATEST_FEATURE_MAP)
    look = (lambda c: dpp_map.get(c)) if dpp_map else (lambda c: None)
    flook = (lambda c: feat_map.get(c)) if feat_map else (lambda c: None)
    npos = nneg = 0
    for s in spots:
        if ingest_observation(cal, s["tx"], s["rx"], True, s["freq_mhz"],
                              s["snr_db"], look, timestamp=s.get("time"),
                              feature_lookup=flook):
            npos += 1
    for g in negatives:
        if dpp_map:
            wmap, _ = path_cell_weights(g["tx"], g["rx"], cal["grid"]["step_deg"])
            dbar = tw = 0.0
            for c, w in wmap.items():
                d = dpp_map.get(c)
                if d is not None:
                    dbar += d * w; tw += w
            if tw > 0 and dbar / tw < 40.0:
                continue
        if ingest_observation(cal, g["tx"], g["rx"], False, g["freq_mhz"],
                              0.0, look, feature_lookup=flook):
            nneg += 1
    record_kappa(cal)
    return {"spots": len(spots), "pos": npos, "neg": nneg, "cone": cone}


def build_grid(region: dict) -> list[tuple[float, float]]:
    lats = np.arange(region["lat_min"], region["lat_max"] + 1e-6, region["step_deg"])
    lons = np.arange(region["lon_min"], region["lon_max"] + 1e-6, region["step_deg"])
    return [(round(float(la), 3), round(float(lo), 3)) for la in lats for lo in lons]


def _hourly_vars() -> list[str]:
    v = ["temperature_2m", "relative_humidity_2m", "surface_pressure",
         "wind_speed_10m", "wind_direction_10m"]
    for p in PRESSURE_LEVELS:
        v += [f"temperature_{p}hPa", f"relative_humidity_{p}hPa",
              f"geopotential_height_{p}hPa", f"wind_speed_{p}hPa",
              f"wind_direction_{p}hPa"]
    return v


def fetch_weather_batch(points: list[tuple[float, float]]) -> list[Optional[dict]]:
    params = {"latitude": ",".join(str(p[0]) for p in points),
              "longitude": ",".join(str(p[1]) for p in points),
              "hourly": ",".join(_hourly_vars()), "forecast_days": 2,
              "timezone": "auto", "wind_speed_unit": "ms"}
    r = requests.get(f"{OPEN_METEO_URL}?{urlencode(params)}", timeout=REQUEST_TIMEOUT)
    r.raise_for_status(); data = r.json()
    return [data] if isinstance(data, dict) else data


def fetch_weather_all(points: list[tuple[float, float]]) -> list[Optional[dict]]:
    res = []
    for i in range(0, len(points), BATCH_SIZE):
        chunk = points[i:i + BATCH_SIZE]
        try:
            res.extend(fetch_weather_batch(chunk))
        except Exception as e:
            print(f"  [気象] バッチ失敗 ({i}): {e}")
            res.extend([None] * len(chunk))
    return res


def build_M_profile(hourly: dict, hi: int) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    heights, M = [], []
    try:
        t0 = hourly["temperature_2m"][hi]; rh0 = hourly["relative_humidity_2m"][hi]
        p0 = hourly["surface_pressure"][hi]
        if None not in (t0, rh0, p0):
            heights.append(2.0)
            M.append(modified_refractivity_M(refractivity_N(t0, rh0, p0), 2.0))
    except (KeyError, IndexError, TypeError):
        pass
    for p in PRESSURE_LEVELS:
        try:
            t = hourly[f"temperature_{p}hPa"][hi]
            rh = hourly[f"relative_humidity_{p}hPa"][hi]
            gh = hourly[f"geopotential_height_{p}hPa"][hi]
            if None in (t, rh, gh):
                continue
            heights.append(gh)
            M.append(modified_refractivity_M(refractivity_N(t, rh, p), gh))
        except (KeyError, IndexError, TypeError):
            continue
    if len(heights) < 3:
        return None, None
    order = np.argsort(heights)
    heights = np.array(heights)[order]; M = np.array(M)[order]
    _, u = np.unique(heights, return_index=True)
    return heights[np.sort(u)], M[np.sort(u)]


def build_wind_profile(hourly: dict, hi: int) -> tuple:
    """高度別の風速/風向を返す(乱流指数=風シアdV/dzの算出用)。"""
    heights, spd, drc = [], [], []
    try:
        s0 = hourly["wind_speed_10m"][hi]; d0 = hourly["wind_direction_10m"][hi]
        if None not in (s0, d0):
            heights.append(10.0); spd.append(s0); drc.append(d0)
    except (KeyError, IndexError, TypeError):
        pass
    for p in PRESSURE_LEVELS:
        try:
            s = hourly[f"wind_speed_{p}hPa"][hi]; d = hourly[f"wind_direction_{p}hPa"][hi]
            gh = hourly[f"geopotential_height_{p}hPa"][hi]
            if None in (s, d, gh):
                continue
            heights.append(gh); spd.append(s); drc.append(d)
        except (KeyError, IndexError, TypeError):
            continue
    if len(heights) < 2:
        return None, None, None
    order = np.argsort(heights)
    heights = np.array(heights)[order]; spd = np.array(spd)[order]; drc = np.array(drc)[order]
    _, u = np.unique(heights, return_index=True)
    return heights[np.sort(u)], spd[np.sort(u)], drc[np.sort(u)]


def _primary(ducts):
    return ducts[0] if ducts else None


def _structure_idx(d):
    if d is None:
        return 0.0
    ts = min(100.0, d["thickness_m"] / 400.0 * 100.0)
    ds = min(100.0, d["delta_M"] / 25.0 * 100.0)
    ar = min(100.0, d["trapped_area"] / 3000.0 * 100.0)
    cv = min(100.0, abs(d["curvature"]) * 5e5)
    return float(max(0.0, min(100.0, 0.35 * ds + 0.30 * ar + 0.20 * ts + 0.15 * cv)))


def _trend(d0, d1):
    p0, p1 = _primary(d0), _primary(d1)
    dm0 = p0["delta_M"] if p0 else 0.0; dm1 = p1["delta_M"] if p1 else 0.0
    dd = dm1 - dm0
    return "growing" if dd > 0.5 else ("collapsing" if dd < -0.5 else "steady")


def _coupling_idx(ducts):
    d = _primary(ducts)
    if d is None:
        return 0.0
    return float(100.0 * coupling_efficiency(TX_HEIGHT_M, d) *
                 coupling_efficiency(RX_HEIGHT_M, d))


def _stability_idx(ducts_series):
    """複数時刻(現在/+3h/+6h等)のΔMのばらつきから持続性を評価。
    ばらつきが小さいほど安定=高スコア。"""
    dms = [(_primary(d)["delta_M"] if _primary(d) else 0.0) for d in ducts_series]
    if not dms or max(dms) <= 0:
        return 0.0
    mean = sum(dms) / len(dms)
    var = sum((x - mean) ** 2 for x in dms) / len(dms)
    cv = math.sqrt(var) / max(mean, 1e-6)
    return float(max(0.0, min(100.0, 100.0 * (1.0 - min(cv, 1.0)))))


def _turbulence_idx(hourly, hi, duct):
    """ダクト層をまたぐ風シア |dV/dz| を乱流の代理指標として0-100で返す。
    値が大きいほど乱流(=ダクトを壊す要因)が強いことを示す(逆評価)。"""
    if duct is None:
        return 0.0
    h, spd, drc = build_wind_profile(hourly, hi)
    if h is None:
        return 0.0
    u = spd * np.sin(np.radians(drc)); v = spd * np.cos(np.radians(drc))
    lo, hi_m = duct["base_m"], duct["top_m"]
    mask = (h >= lo - 50.0) & (h <= hi_m + 50.0)
    idxs = np.where(mask)[0]
    if len(idxs) < 2:
        near = int(np.argmin(np.abs(h - (lo + hi_m) / 2.0)))
        i0, i1 = max(0, near - 1), min(len(h) - 1, near + 1)
    else:
        i0, i1 = int(idxs[0]), int(idxs[-1])
    if i1 == i0:
        return 0.0
    dz = max(float(h[i1] - h[i0]), 1.0)
    shear = math.hypot(float(u[i1] - u[i0]), float(v[i1] - v[i0])) / dz
    return float(max(0.0, min(100.0, shear / TURBULENCE_SHEAR_CAP * 100.0)))


def _wkb_detail(heights, M, duct, freq_mhz=REFERENCE_BAND_MHZ):
    """代表帯域でのWKB捕捉結果(技術者向け表示用)。"""
    if duct is None or heights is None:
        return None
    wkb = wkb_trapped_modes(heights, M, duct, freq_mhz)
    leak = leakage_db_per_km(freq_mhz, wkb)
    coup = coupling_efficiency(TX_HEIGHT_M, duct) * coupling_efficiency(RX_HEIGHT_M, duct)
    return {"freq_mhz": freq_mhz,
            "f_min_mhz": round(wkb["f_min_mhz"], 1) if math.isfinite(wkb["f_min_mhz"]) else None,
            "n_modes": wkb["n_modes"], "trapped": wkb["trapped"],
            "coupling": round(coup, 3),
            "leak_db_per_100km": round(leak * 100.0, 2) if leak is not None else None}


def _generate_summary(ducts_now, trend, trap_str, kappa, raw_dpp, corr_dpp):
    """規則ベースの一文サマリ(要点を自然文で提示)。"""
    d = _primary(ducts_now)
    if d is None:
        return "この地点では有意なダクト層は検出されていません。伝搬は標準大気に近い状態です。"
    kind = "接地ダクト" if d["kind"] == "surface" else "上空ダクト"
    trend_s = {"growing": "さらに成長する", "collapsing": "弱まっていく",
               "steady": "現状を維持する"}[trend]
    band_order = list(BANDS.keys())
    reachable = [b for b, t in zip(band_order, trap_str) if t == "1"]
    band_s = (reachable[-1] + "以上は十分捕捉可能" if reachable else "いずれの帯域も捕捉は困難")
    if corr_dpp > raw_dpp + 3.0:
        obs_s = "現在のPSKReporter実績はこれを上回っており、予測は上方修正されています"
    elif corr_dpp < raw_dpp - 3.0:
        obs_s = "現在のPSKReporter実績はやや下回っており、予測は下方修正されています"
    else:
        obs_s = "現在のPSKReporter実績とも一致しています"
    conf_s = "高め" if kappa >= 0.6 else ("中程度" if kappa >= 0.3 else "低め")
    return f"この地点では{kind}が形成されており、{trend_s}見込みです。{band_s}で、{obs_s}。信頼度は{conf_s}です。"


def _frame_bucket(times, hi, fallback_dt=None):
    try:
        return _time_bucket(datetime.fromisoformat(times[hi]))
    except Exception:
        return _time_bucket(fallback_dt or datetime.now())


def _current_hour_index(hourly):
    """hourly['time']配列(forecast_days=2で当日0:00始まり)の中から、
    現在時刻に最も近い(直近未来含む)インデックスを返す。
    これをフレーム0(=「現在」表示)の基準にすることで、常に0:00表示になる
    問題を回避する。"""
    times = hourly.get("time", [])
    if not times:
        return 0
    target = datetime.now().strftime("%Y-%m-%dT%H:00")
    for i, t in enumerate(times):
        if t >= target:
            return i
    return max(0, len(times) - 1)


def analyze_point_frames(lat, lon, wj, cal):
    if wj is None or "hourly" not in wj:
        return None
    hourly = wj["hourly"]; nt = len(hourly.get("time", []))
    if nt < 2:
        return None
    time_buckets_on = cal.get("time_buckets_enabled", True)
    base_hi = _current_hour_index(hourly)
    times = hourly.get("time", [])
    cell0 = _lookup_cell(cal, lat, lon,
                         _frame_bucket(times, base_hi) if time_buckets_on else None)
    n_obs = int(round(cell0.get("n_obs", 0))) if cell0 else 0
    frames = []; ducts0 = None; mprof0 = None
    early_ducts = []; early_mprofs = []
    raw_dpp0 = None; wkb0 = None; turb0 = None
    corr0 = None; kappa0 = None
    feat0 = None
    for off in FRAME_HOURS:
        hi = base_hi + off
        if hi >= nt:
            frames.append(None); continue
        h, M = build_M_profile(hourly, hi)
        if h is None:
            frames.append(None); continue
        ducts = detect_all_ducts(h, M)
        hn = hi + 3 if hi + 3 < nt else min(hi + 1, nt - 1)
        hh, MM = build_M_profile(hourly, hn)
        ducts_n = detect_all_ducts(hh, MM) if hh is not None else []
        bands, refr = band_reachability(h, M, ducts)
        d_now = _primary(ducts)
        s = _structure_idx(d_now)
        st = {"growing": 75.0, "collapsing": 35.0,
              "steady": 55.0}[_trend(ducts, ducts_n)]
        c = _coupling_idx(ducts)
        turb = _turbulence_idx(hourly, hi, d_now)
        bucket = _frame_bucket(times, hi) if time_buckets_on else None
        n_modes = bands.get("144MHz", {}).get("n_modes", 0) if d_now else 0
        feat = None
        if d_now is not None:
            feat = {"delta_M": float(d_now.get("delta_M", 0.0)),
                    "thickness_m": float(d_now.get("thickness_m", 0.0)),
                    "base_m": float(d_now.get("base_m", 0.0)),
                    "turbulence": float(turb), "coupling": float(c),
                    "n_modes": float(n_modes)}
        corr, kappa = apply_correction(cal, lat, lon, bucket=bucket, features=feat)
        raw = 0.35 * s + 0.20 * st + 0.20 * (100.0 - turb) + 0.25 * c
        dpp = float(max(0.0, min(100.0, raw * corr)))
        r430 = bands.get("430MHz", {}).get("horizon_km", 0.0)
        trapped = "".join("1" if bands[b]["trapped"] else "0" for b in BANDS)
        fmin = {b: bands[b].get("f_min_mhz") for b in BANDS}
        frames.append({"dpp": round(dpp, 1), "r": round(r430, 0),
                       "trap": trapped, "fmin": fmin,
                       "t": times[hi][-5:] if hi < len(times) else str(hi)})
        if ducts0 is None:
            ducts0 = ducts; mprof0 = [(float(z), float(m)) for z, m in zip(h, M)]
        if len(early_ducts) < 3:
            early_ducts.append(ducts)
            early_mprofs.append({"t": times[hi][-5:] if hi < len(times) else str(hi),
                                 "pts": [(float(z), float(m)) for z, m in zip(h, M)]})
        if raw_dpp0 is None:
            raw_dpp0 = float(max(0.0, min(100.0, raw)))
            wkb0 = _wkb_detail(h, M, d_now)
            turb0 = turb
            corr0 = corr; kappa0 = kappa; feat0 = feat
    if all(f is None for f in frames):
        return None
    d0 = _primary(ducts0)
    components = {"structure": round(_structure_idx(d0), 1),
                  "stability": round(_stability_idx(early_ducts or [ducts0]), 1),
                  "coupling": round(_coupling_idx(ducts0), 1),
                  "turbulence": round(turb0 if turb0 is not None else 0.0, 1)}
    kappa = kappa0 if kappa0 is not None else 0.0
    corr = corr0 if corr0 is not None else 1.0
    correction = {"raw_dpp": round(raw_dpp0 or 0.0, 1),
                  "corr_dpp": round((frames[0] or {}).get("dpp", raw_dpp0 or 0.0), 1),
                  "factor": round(corr, 3),
                  "kappa": round(kappa, 3), "n_obs": n_obs}
    trend0 = _trend(ducts0, early_ducts[1] if len(early_ducts) > 1 else ducts0)
    trap0 = (frames[0] or {}).get("trap", "0" * len(BANDS))
    summary = _generate_summary(ducts0, trend0, trap0, kappa,
                                correction["raw_dpp"], correction["corr_dpp"])
    return {"lat": lat, "lon": lon, "confidence": round(kappa, 3),
            "ducts0": ducts0 or [], "mprof": mprof0 or [], "frames": frames,
            "mprofs3": early_mprofs, "components": components,
            "wkb": wkb0, "correction": correction, "summary": summary,
            "_features": feat0}


def _publish_dpp_map(results, step):
    m = {}; fm = {}
    for r in results:
        if not r or not r["frames"]:
            continue
        f0 = next((f for f in r["frames"] if f), None)
        if not f0:
            continue
        key = _cell_key(r["lat"], r["lon"], step)
        if f0["dpp"] >= m.get(key, -1.0):
            m[key] = f0["dpp"]
            if r.get("_features"):
                fm[key] = r["_features"]
    with _DPP_LOCK:
        _LATEST_DPP_MAP.clear(); _LATEST_DPP_MAP.update(m)
        _LATEST_FEATURE_MAP.clear(); _LATEST_FEATURE_MAP.update(fm)


_HTML_TEMPLATE = r"""<!DOCTYPE html><!-- DUCT_BUILD_MARKER=OSM_TILES_REV16_7 -->
<html lang="ja"><head><meta charset="utf-8">
<title>__TITLE__</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="refresh" content="__REFRESH__">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.css"/>
<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.min.js"></script>
<style>
:root{--bg:#0f172a;--panel:#111c33;--fg:#e6edf7;--sub:#93a4bd;--line:#26344d;}
body.light{--bg:#eef2f8;--panel:#fff;--fg:#12203a;--sub:#4a5b76;--line:#cdd7e6;}
*{box-sizing:border-box}
body{margin:0;font-family:'Segoe UI','Hiragino Kaku Gothic ProN',sans-serif;
 background:var(--bg);color:var(--fg);}
#hd{padding:9px 14px;background:var(--panel);border-bottom:1px solid var(--line);
 display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px}
#hd h1{margin:0;font-size:15px;font-weight:600}
.meta{font-size:11px;color:var(--sub);display:flex;gap:10px;flex-wrap:wrap;align-items:center}
.badge{background:var(--bg);border:1px solid var(--line);border-radius:11px;padding:2px 9px}
.badge.eta{border-color:#5ec962;color:#7fe0a0}
#bar{display:flex;gap:12px;align-items:center;padding:7px 14px;background:var(--panel);
 border-bottom:1px solid var(--line);font-size:12px;flex-wrap:wrap}
#map{height:calc(100vh - 150px);width:100%;background:var(--bg)}
button{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:4px 11px;cursor:pointer;font-size:12px}
button:hover{border-color:#5b7bb5}
#tslider{width:230px;vertical-align:middle}
#tlabel{font-weight:600;min-width:96px;display:inline-block}
.legend{display:flex;gap:6px;align-items:center;margin-left:auto}
.grad{width:150px;height:11px;border-radius:6px;
 background:linear-gradient(90deg,#3b0f70,#3b528b,#21918c,#5ec962,#fde725)}
.leaflet-popup-content{margin:9px 11px;color:#0b1220}
#ovl{display:none;position:fixed;inset:0;background:rgba(0,0,0,0.5);z-index:1000}
#panel{position:absolute;top:60px;right:20px;width:330px;background:var(--panel);
 border:1px solid var(--line);border-radius:10px;padding:14px 16px;
 box-shadow:0 10px 40px rgba(0,0,0,0.4)}
#panel .phd{display:flex;justify-content:space-between;align-items:center;
 font-size:14px;margin-bottom:6px}
#panel .x{cursor:pointer;color:var(--sub)}
#panel label{display:flex;justify-content:space-between;align-items:center;
 font-size:12px;margin:9px 0;gap:10px}
#panel select{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:4px 8px;font-size:12px;min-width:120px}
#panel input[type=text]{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:4px 8px;font-size:12px}
#panel .note{font-size:11px;color:var(--sub);margin:2px 0 8px}
#panel .prow{display:flex;gap:8px;align-items:center;margin-top:10px}
#panel .smsg{font-size:11px;color:#5ec962}
#panel .sstat{font-size:11px;color:var(--sub);margin-top:10px;line-height:1.7;
 border-top:1px solid var(--line);padding-top:8px}
#covl{display:none;position:fixed;inset:0;background:rgba(0,0,0,0.5);z-index:1000}
#cpanel{position:absolute;top:60px;right:20px;width:310px;background:var(--panel);
 border:1px solid var(--line);border-radius:10px;padding:14px 16px;
 box-shadow:0 10px 40px rgba(0,0,0,0.4)}
#cpanel .phd{display:flex;justify-content:space-between;align-items:center;
 font-size:14px;margin-bottom:6px}
#cpanel .x{cursor:pointer;color:var(--sub)}
#cpanel label{display:block;font-size:12px;margin:10px 0;color:var(--sub)}
#cpanel .cprow{display:flex;gap:6px;align-items:center;margin-top:4px}
#cpanel input[type=text]{background:var(--bg);color:var(--fg);border:1px solid var(--line);
 border-radius:6px;padding:5px 8px;font-size:12px;width:130px;flex:1}
#cpanel .gpsbtn{font-size:11px;padding:5px 7px;border-radius:6px;border:1px solid var(--line);
 background:var(--bg);color:var(--fg);cursor:pointer;white-space:nowrap}
#cpanel .gpsbtn:active{background:#5ec962;color:#111}
#cpanel .prow{display:flex;gap:8px;align-items:center;margin-top:12px}
#cpanel .smsg{font-size:11px;color:#5ec962;min-height:14px;display:block;margin-top:6px}
.radar-star{font-size:22px;line-height:22px;text-align:center;color:#ffd23f;
 text-shadow:0 0 3px #000,0 0 4px #000,0 0 6px #000}
.radar-tip{font-size:11px}
.bar2{height:7px;border-radius:4px;background:var(--bg);border:1px solid var(--line);
 overflow:hidden;margin-top:3px}
.bar2 i{display:block;height:100%;background:linear-gradient(90deg,#3b528b,#5ec962)}
#rpanel{position:fixed;top:0;right:-420px;width:380px;max-width:92vw;height:100%;
 background:var(--panel);border-left:1px solid var(--line);
 box-shadow:-10px 0 40px rgba(0,0,0,0.45);z-index:1100;
 transition:right .22s ease;overflow-y:auto;padding:14px 16px 30px}
#rpanel.open{right:0}
#rpanel .phd{display:flex;justify-content:space-between;align-items:center;
 font-size:14px;margin-bottom:8px;position:sticky;top:-14px;background:var(--panel);
 padding:6px 0}
#rpanel .x{cursor:pointer;color:var(--sub);font-size:16px}
#rpanel .rp-dpp{font-size:26px;font-weight:700;color:var(--fg)}
#rpanel .rp-sub{font-size:11px;color:var(--sub)}
#rpanel .rp-sec{margin-top:10px;font-size:11px;color:var(--fg)}
#rpanel .rp-sec b{font-size:12px}
#rpanel table{width:100%;font-size:11px;margin-top:4px;border-collapse:collapse}
#rpanel th,#rpanel td{padding:2px 4px}
#rpanel .rp-charts{display:flex;gap:6px;flex-wrap:wrap;margin-top:4px}
#rpanel .rp-summary{font-size:11px;margin-top:10px;padding:8px;background:#0b1220;
 color:#cbd5e1;border-radius:6px;line-height:1.6}
body.light #rpanel .rp-summary{background:#e6edf7;color:#12203a}
#hpanel{position:fixed;left:12px;bottom:12px;z-index:1100;display:none;
 background:var(--panel);border:1px solid var(--line);border-radius:8px;
 box-shadow:0 6px 24px rgba(0,0,0,0.35);padding:8px 12px;font-size:12px;
 color:var(--fg);max-width:88vw}
#hpanel .hp-row{display:flex;align-items:center;gap:6px;white-space:nowrap}
#hpanel .hp-sw{display:inline-block;width:12px;height:12px;border-radius:2px;
 border:1px solid var(--line)}
#hpanel .hp-sub{color:var(--sub);font-size:10px;margin-left:4px}
</style></head><body class="dark">
<div id="hd"><h1>🗾 __TITLE__</h1><div class="meta">
<span class="badge">生成 __GEN__</span>
<span class="badge" id="verBadge" title="このHTMLを生成したPythonのバージョン">v__VER__</span>
<span class="badge">地点 __NPTS__</span>
<span class="badge">学習成熟度 __MATURITY__%</span>
<span class="badge eta" id="etaBadge">実用まで __ETA__</span>
<span class="badge">平均信頼度 __MCONF__</span>
<span class="badge">被覆 __COVERED__セル</span>
<span class="badge" id="cd">次回更新 --:--</span>
</div></div>
<div id="bar">
<button id="play">▶ 再生</button>
<label>予報時刻 <span id="tlabel"></span></label>
<input type="range" id="tslider" min="0" max="__MAXFR__" value="0" step="1">
<button id="pauseBtn" title="観測・DPP計算エンジンの一時停止/再開（時刻アニメの再生ボタンとは別機能）">⏯ 観測稼働中</button>
<button onclick="toggleTheme()">🌓</button>
<button id="contourBtn" onclick="toggleContour()" title="ダクト形成可能性が高い経路をコンター線(等値線)で表示">📈 コンター</button>
<button onclick="openSettings()" title="設定">⚙ 設定</button>
<button onclick="askPass()" title="自局と相手局、2地点間のPass(伝搬経路)を診断">🛤️ Pass</button>
<button onclick="askAntenna()" title="自局1地点のみ。アンテナ高を上げた場合の効果を試算">📶 アンテナ高</button>
<button onclick="askRadar()" title="自局1地点のみ。そこから全方位への到達距離を表示">🧭 移動運用</button>
<span class="badge" id="readout" title="地図をタップ(またはマウスで指した地点)の推定DPP値を表示します">📍 地点DPP: タップして確認</span>
<div class="legend"><span>0</span><div class="grad"></div><span>100 (DPP指数)</span></div>
<span class="badge" id="contourLegend" style="display:none">等値線: 30 / 50 / 70(実線太) / 85(実線最太)</span>
</div>
<div id="ovl" onclick="if(event.target===this)closeSettings()">
 <div id="panel">
  <div class="phd"><b>⚙ 設定</b><span class="x" onclick="closeSettings()">✕</span></div>
  <p class="note">変更は<b>保存で即時反映</b>(サーバへ保存し即再計算)。</p>
  <label>更新周期 <select id="s_interval"></select></label>
  <label>常駐時間 <select id="s_total" title="0=無制限。設定時間が経過するとサーバーごと自動停止します"></select></label>
  <p class="note" style="margin-top:-6px">⚠ 0以外を選ぶと、その時間経過でサーバーが自動的に停止します(0=無制限、ずっと動かす場合はこちら)</p>
  <label>気象キャッシュ <select id="s_wx"></select></label>
  <label>スポット遡り <select id="s_back"></select></label>
  <p class="note">🔔 アラート(Discord Webhook通知)</p>
  <label>有効化 <input type="checkbox" id="s_alert_on"></label>
  <label>Webhook URL <input type="text" id="s_webhook" placeholder="https://discord.com/api/webhooks/..." style="width:200px"></label>
  <label>DPP閾値 <select id="s_alert_dpp"></select></label>
  <label>監視半径(km) <select id="s_alert_radius"></select></label>
  <div class="prow">
   <button onclick="applySettings()">保存して即時更新</button>
   <button onclick="recomputeNow()">今すぐ再計算</button>
   <button onclick="closeSettings()">閉じる</button>
  </div>
  <div class="prow"><span id="s_msg" class="smsg"></span></div>
  <div id="s_status" class="sstat"></div>
 </div>
</div>
<div id="rpanel">
 <div class="phd"><b>📡 地点診断レポート</b><span class="x" onclick="closeReport()">✕</span></div>
 <div id="rpbody"></div>
</div>
<div id="covl" onclick="if(event.target===this)closeCoordModal()">
 <div id="cpanel">
  <div class="phd"><b id="cp_title">座標入力</b><span class="x" onclick="closeCoordModal()">✕</span></div>
  <div id="cp_body"></div>
  <div class="prow">
   <button onclick="cp_submit()">実行</button>
   <button onclick="closeCoordModal()">閉じる</button>
  </div>
  <span id="cp_msg" class="smsg"></span>
 </div>
</div>
<div id="hpanel"></div>
<div id="map"></div>
<script>
var CELLS=__CELLDATA__, FRAMES_T=__FRAMET__, STEP=__STEP__, REGION=__REGION__;
var ANTENNA_DEFAULT_H=__TXHEIGHT__;
var BAND_LABELS=__BANDLIST__;
function viridis(t){t=Math.max(0,Math.min(1,t/100));
 var s=[[59,15,112],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
 var x=t*(s.length-1),i=Math.floor(x),f=x-i;if(i>=s.length-1){i=s.length-2;f=1;}
 var a=s[i],b=s[i+1];return [Math.round(a[0]+(b[0]-a[0])*f),
  Math.round(a[1]+(b[1]-a[1])*f),Math.round(a[2]+(b[2]-a[2])*f)];}
var map=L.map('map',{zoomControl:true});
/* Rev16.7(バグ修正): 従来使用していたCARTO Voyager(basemaps.cartocdn.com)は
   これ以外の場所は自由に使えていた無登録の無料タイルだが、CARTO側の運用方針
   変更により、未登録(APIキー無し)でのアクセス時にタイル画像内へ
   「Get your own API key」等の透かし文字が焼き込まれて配信されるようになった。
   このためDUCT地図の背景にAPIキーを促す文言が常時表示されてしまっていた。
   対策として、登録・APIキーが一切不要で透かしも入らない標準の
   OpenStreetMap公式タイル(tile.openstreetmap.org)へ切り替える。
   個人の趣味利用・低頻度アクセス(このダッシュボードのリフレッシュ周期程度)は
   OSMの利用規約(Tile Usage Policy)の範囲内として問題ない。 */
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
 {attribution:'&copy; OpenStreetMap contributors',subdomains:'abc',maxZoom:19}).addTo(map);
map.fitBounds([[REGION.lat_min,REGION.lon_min],[REGION.lat_max,REGION.lon_max]]);
var fcv=L.DomUtil.create('canvas','',map.getPanes().overlayPane);
fcv.id='fieldcv';fcv.style.position='absolute';fcv.style.pointerEvents='none';
var curFr=0;
function cellPixel(){var a=map.latLngToContainerPoint([REGION.lat_min,REGION.lon_min]);
 var b=map.latLngToContainerPoint([REGION.lat_min,REGION.lon_min+STEP]);
 return Math.max(10,Math.hypot(b.x-a.x,b.y-a.y));}
function renderField(){var size=map.getSize();fcv.width=size.x;fcv.height=size.y;
 L.DomUtil.setPosition(fcv,map.containerPointToLayerPoint([0,0]));
 var ctx=fcv.getContext('2d');ctx.clearRect(0,0,size.x,size.y);
 var pts=[];
 for(var i=0;i<CELLS.length;i++){var f=CELLS[i].frames[curFr];if(!f)continue;
  var p=map.latLngToContainerPoint([CELLS[i].lat,CELLS[i].lon]);
  pts.push({x:p.x,y:p.y,v:f.dpp,c:CELLS[i].conf});}
 if(!pts.length)return;
 var R=cellPixel()*2.6,R2=R*R,block=6;
 for(var y=0;y<size.y;y+=block){for(var x=0;x<size.x;x+=block){
  var sw=0,sv=0,sc=0,dmin=1e18;
  for(var k=0;k<pts.length;k++){var dx=x-pts[k].x,dy=y-pts[k].y,d2=dx*dx+dy*dy;
   if(d2>R2)continue;if(d2<dmin)dmin=d2;var w=1.0/(d2+25.0);
   sw+=w;sv+=w*pts[k].v;sc+=w*pts[k].c;}
  if(sw<=0)continue;var val=sv/sw,conf=sc/sw,dm=Math.sqrt(dmin);
  var mask=dm<R*0.30?1.0:Math.max(0,1.0-(dm-R*0.30)/(R*0.70));
  mask*=(0.45+0.55*conf);if(mask<=0.02)continue;
  var cc=viridis(val);
  ctx.fillStyle='rgba('+cc[0]+','+cc[1]+','+cc[2]+','+(0.74*mask).toFixed(3)+')';
  ctx.fillRect(x,y,block,block);}}
 drawContours(ctx,curFr);}
function idwAt(px,py){var sw=0,sv=0,R=cellPixel()*2.6,R2=R*R;
 for(var i=0;i<CELLS.length;i++){var f=CELLS[i].frames[curFr];if(!f)continue;
  var p=map.latLngToContainerPoint([CELLS[i].lat,CELLS[i].lon]);
  var dx=px-p.x,dy=py-p.y,d2=dx*dx+dy*dy;if(d2>R2)continue;
  var w=1.0/(d2+25.0);sw+=w;sv+=w*f.dpp;}
 return sw>0?sv/sw:null;}
/* ---- コンター(等値線)表示: 経度緯度グリッド上でIDW補間しMarching Squaresで等値線を抽出 ---- */
var CONTOUR_LEVELS=[30,50,70,85];
/* ITPFS Ver1.0(修正4・重要): 「統合前は初回からコンターが出ていたのに、統合後は
   一度も出ない」という報告の実際の原因はこれだった。CONTOUR_ONの既定値は
   localStorageの'dpp_contour'キーが文字列'1'の場合のみtrueになる設計だが、
   localStorageはオリジン(スキーム+ホスト+ポート)ごとに完全に分離される。
   統合前の単体版DUCTは別ポート(例:8787)で動いていたため、そちらで一度
   コンターをONにした形跡があってもITPFS側(別ポート、例:8765)には一切
   引き継がれず、常に「未設定」→既定値falseになっていた。コンターボタンの
   緑色っぽい枠線はCSSの既定スタイルであり、初回表示時にCONTOUR_ONの実値を
   反映していなかったため、見た目は「ON」に見えるのに実際は「OFF」という
   食い違いが起きていた。対策として、明示的に'0'が保存されている場合のみ
   OFFとし、それ以外(未設定=新規オリジン含む)は既定でONにする。 */
var CONTOUR_ON=(function(){try{var v=localStorage.getItem('dpp_contour');
 return v===null?true:(v==='1');}catch(e){return true;}})();
var contourCache={};
function idwAtLatLon(lat,lon,fr,R2deg){
 var sw=0,sv=0;
 for(var i=0;i<CELLS.length;i++){var f=CELLS[i].frames[fr];if(!f)continue;
  var dlat=lat-CELLS[i].lat,dlon=(lon-CELLS[i].lon)*Math.cos(lat*Math.PI/180);
  var d2=dlat*dlat+dlon*dlon;if(d2>R2deg)continue;
  var w=1.0/(d2+1e-7);sw+=w;sv+=w*f.dpp;}
 return sw>0?sv/sw:NaN;}
function buildGrid(fr,res){
 res=res||44;
 var latMin=REGION.lat_min,latMax=REGION.lat_max,lonMin=REGION.lon_min,lonMax=REGION.lon_max;
 var dlat=(latMax-latMin)/res,dlon=(lonMax-lonMin)/res;
 // ITPFS Ver1.0(修正): 初回サイクルは高速化のため粗いグリッド(STEPが大きい)で
 // 計算するため、従来通りR=STEP*2.6のままだと補間半径が地域全体を覆うほど
 // 巨大化し、DPP値がほぼ均一化されて等値線が1本も引けなくなる不具合があった。
 // 半径の上限を設け、粗い初回サイクルでも等値線が出るようにする。
 var R=Math.min(STEP*2.6, 6.0),R2=R*R;
 var vals=new Float32Array((res+1)*(res+1));
 for(var iy=0;iy<=res;iy++){var lat=latMin+iy*dlat;
  for(var ix=0;ix<=res;ix++){var lon=lonMin+ix*dlon;
   vals[iy*(res+1)+ix]=idwAtLatLon(lat,lon,fr,R2);}}
 return {res:res,latMin:latMin,lonMin:lonMin,dlat:dlat,dlon:dlon,vals:vals};}
function msInterp(v1,v2,p1,p2,level){
 if(!isFinite(v1)||!isFinite(v2))return p1;
 var t=(level-v1)/(v2-v1);if(!isFinite(t))t=0.5;t=Math.max(0,Math.min(1,t));
 return [p1[0]+(p2[0]-p1[0])*t,p1[1]+(p2[1]-p1[1])*t];}
function marchingSquares(grid,level){
 var res=grid.res,vals=grid.vals,segs=[];
 function vAt(ix,iy){return vals[iy*(res+1)+ix];}
 for(var iy=0;iy<res;iy++){for(var ix=0;ix<res;ix++){
  var tl=vAt(ix,iy),tr=vAt(ix+1,iy),br=vAt(ix+1,iy+1),bl=vAt(ix,iy+1);
  if(!isFinite(tl)||!isFinite(tr)||!isFinite(br)||!isFinite(bl))continue;
  var c=0;if(tl>level)c|=8;if(tr>level)c|=4;if(br>level)c|=2;if(bl>level)c|=1;
  if(c===0||c===15)continue;
  var latT=grid.latMin+iy*grid.dlat,lonL=grid.lonMin+ix*grid.dlon;
  var latB=latT+grid.dlat,lonR=lonL+grid.dlon;
  var pTL=[latT,lonL],pTR=[latT,lonR],pBR=[latB,lonR],pBL=[latB,lonL];
  var top=msInterp(tl,tr,pTL,pTR,level),right=msInterp(tr,br,pTR,pBR,level),
   bottom=msInterp(bl,br,pBL,pBR,level),left=msInterp(tl,bl,pTL,pBL,level);
  switch(c){
   case 1:segs.push([left,bottom]);break;
   case 2:segs.push([bottom,right]);break;
   case 3:segs.push([left,right]);break;
   case 4:segs.push([top,right]);break;
   case 5:segs.push([left,top]);segs.push([bottom,right]);break;
   case 6:segs.push([top,bottom]);break;
   case 7:segs.push([left,top]);break;
   case 8:segs.push([top,left]);break;
   case 9:segs.push([top,bottom]);break;
   case 10:segs.push([top,right]);segs.push([left,bottom]);break;
   case 11:segs.push([top,right]);break;
   case 12:segs.push([left,right]);break;
   case 13:segs.push([bottom,right]);break;
   case 14:segs.push([left,bottom]);break;}}}
 return segs;}
function getContours(fr){
 if(contourCache[fr])return contourCache[fr];
 var grid=buildGrid(fr),out={};
 for(var i=0;i<CONTOUR_LEVELS.length;i++){out[CONTOUR_LEVELS[i]]=marchingSquares(grid,CONTOUR_LEVELS[i]);}
 contourCache[fr]=out;return out;}
function drawContours(ctx,fr){
 if(!CONTOUR_ON)return;
 var cs=getContours(fr);
 for(var li=0;li<CONTOUR_LEVELS.length;li++){
  var lvl=CONTOUR_LEVELS[li],segs=cs[lvl];if(!segs||!segs.length)continue;
  var col=viridis(lvl);
  ctx.strokeStyle='rgba('+col[0]+','+col[1]+','+col[2]+',0.95)';
  ctx.lineWidth=lvl>=85?3.0:(lvl>=70?2.2:1.3);
  ctx.setLineDash(lvl>=70?[]:[6,4]);
  ctx.beginPath();
  for(var i=0;i<segs.length;i++){
   var a=map.latLngToContainerPoint([segs[i][0][0],segs[i][0][1]]);
   var b=map.latLngToContainerPoint([segs[i][1][0],segs[i][1][1]]);
   ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);}
  ctx.stroke();
  if(segs.length){ctx.setLineDash([]);
   var mid=segs[Math.floor(segs.length/2)];
   var mp=map.latLngToContainerPoint([(mid[0][0]+mid[1][0])/2,(mid[0][1]+mid[1][1])/2]);
   ctx.font='10px sans-serif';ctx.fillStyle='rgba('+col[0]+','+col[1]+','+col[2]+',0.95)';
   ctx.fillText(String(lvl),mp.x+3,mp.y-3);}}
 ctx.setLineDash([]);}
function toggleContour(){CONTOUR_ON=!CONTOUR_ON;
 try{localStorage.setItem('dpp_contour',CONTOUR_ON?'1':'0');}catch(e){}
 var b=document.getElementById('contourBtn');if(b)b.style.borderColor=CONTOUR_ON?'#5ec962':'';
 var lg=document.getElementById('contourLegend');if(lg)lg.style.display=CONTOUR_ON?'':'none';
 renderField();}
map.on('moveend zoomend resize',renderField);
map.on('click',function(e){var best=-1,bd=1e9;
 for(var i=0;i<CELLS.length;i++){var dl=CELLS[i].lat-e.latlng.lat,
  dn=CELLS[i].lon-e.latlng.lng,d=dl*dl+dn*dn;if(d<bd){bd=d;best=i;}}
 if(best>=0&&bd<STEP*STEP)openPopup(best);
 fetchHepburnAt(e.latlng.lat,e.latlng.lng);
 updateReadout(idwAt(e.containerPoint.x,e.containerPoint.y));});
function fetchHepburnAt(lat,lon){
 var el=document.getElementById('hpanel');
 if(!el)return;
 el.style.display='block';
 el.innerHTML='<div class="hp-row">Hepburn Map 照会中…</div>';
 fetch('/api/duct/hepburn?lat='+lat.toFixed(3)+'&lon='+lon.toFixed(3))
  .then(function(r){return r.json();})
  .then(function(d){
   if(!d.ok){
    var msg=d.error==='not_calibrated (geo)'?'未較正':
     d.error==='out_of_map_bounds'?'マップ範囲外':d.error;
    el.innerHTML='<div class="hp-row">Hepburn: '+msg+'</div>';return;}
   var rgb='rgb('+d.rgb.join(',')+')';
   el.innerHTML='<div class="hp-row">'
    +'<span class="hp-sw" style="background:'+rgb+'"></span>'
    +'<b>'+(d.label||'—')+'</b>'
    +'<span class="hp-sub"> ('+lat.toFixed(2)+', '+lon.toFixed(2)+') '
    +'valid '+d.valid_utc.replace('T',' ').slice(0,16)+'UTC</span></div>';
  })
  .catch(function(){el.innerHTML='<div class="hp-row">Hepburn: 通信エラー</div>';});}
function updateReadout(v){
 var el=document.getElementById('readout');
 if(el)el.innerText=(v==null?'📍 地点DPP: データ無し':'📍 地点DPP: '+v.toFixed(0));}
map.on('mousemove',function(e){updateReadout(idwAt(e.containerPoint.x,e.containerPoint.y));});
function draw(fr){curFr=fr;renderField();
 document.getElementById('tlabel').innerText=FRAMES_T[fr]+' JST';
 try{localStorage.setItem('dpp_fr',fr);}catch(e){}
 if(curReportIdx>=0&&document.getElementById('rpanel').classList.contains('open')){
  openPopup(curReportIdx);}}
var curReportIdx=-1;
function openPopup(idx){curReportIdx=idx;var c=CELLS[idx],fr=+slider.value,f=c.frames[fr];
 if(!f)return;var bands=BAND_LABELS;
 var brow=bands.map(function(b,k){return '<tr><td>'+b+'MHz</td><td style=\"text-align:center\">'
  +(f.trap[k]==='1'?'○':'×')+'</td></tr>';}).join('');
 var ds=c.duct||'ダクトなし';
 var wkbHtml='';
 if(c.wkb){var w=c.wkb;
  wkbHtml='<div class=\"rp-sec\"><b>WKB('+w.freq_mhz+'MHz)</b><br>'
   +'捕捉周波数下限: <b>'+(w.f_min_mhz!=null?w.f_min_mhz.toFixed(0)+'MHz':'—')+'</b> / '
   +'捕捉モード: <b>'+w.n_modes+'</b><br>'
   +'結合効率: <b>'+w.coupling.toFixed(2)+'</b> / '
   +'推定漏洩: <b>'+(w.leak_db_per_100km!=null?w.leak_db_per_100km.toFixed(1)+' dB/100km':'—')+'</b></div>';}
 var corrHtml='';
 if(c.corr){var cr=c.corr;
  corrHtml='<div class=\"rp-sec\"><b>PSKReporter補正</b><br>'
   +'合成前 <b>'+cr.raw_dpp.toFixed(0)+'</b> → PSK補正 <b>×'+cr.factor.toFixed(2)+'</b> → '
   +'最終 <b>'+cr.corr_dpp.toFixed(0)+'</b><br>'
   +'観測数: <b>'+cr.n_obs+'</b> / 信頼度κ: <b>'+cr.kappa.toFixed(2)+'</b></div>';}
 var hepHtml='';
 if(c.corr&&c.corr.hepburn){var hp=c.corr.hepburn;
  if(hp.applied){
   hepHtml='<div class=\"rp-sec\"><b>Hepburn Map補正</b><br>'
    +'補正前 <b>'+hp.pre_dpp.toFixed(0)+'</b> → Hepburn期待値 <b>'+hp.expected_from_hepburn.toFixed(0)+'</b>'
    +'(κ='+hp.kappa.toFixed(2)+') → 最終 <b>'+hp.final_dpp.toFixed(0)+'</b><br>'
    +'回帰サンプル数: <b>'+hp.n+'</b> / R²: <b>'+hp.r2.toFixed(2)+'</b></div>';
  }else{
   hepHtml='<div class=\"rp-sec\"><b>Hepburn Map補正</b><br>'
    +'<span class=\"rp-sub\">サンプル蓄積中('+(hp.n||0)+'/'+(hp.min_required||'?')+') — 未適用</span></div>';
  }}
 var sumHtml=c.summary?('<div class=\"rp-summary\">'+c.summary+'</div>'):'';
 var body='<div class=\"rp-sub\">('+c.lat.toFixed(1)+', '+c.lon.toFixed(1)+')・予報時刻 '+FRAMES_T[fr]+'</div>'
  +'<div class=\"rp-dpp\">DPP '+f.dpp.toFixed(1)+'<span class=\"rp-sub\"> /100 (信頼度 '+c.conf.toFixed(2)+')</span></div>'
  +'<div class=\"rp-sub\">予測到達(430MHz): <b>'+f.r.toFixed(0)+' km</b> ・ 構造: '+ds+'</div>'
  +'<div class=\"rp-charts\">'+(c.svg||'')+(c.struct||'')+'</div>'
  +'<div class=\"rp-sec\"><b>DPP内訳</b>'+(c.comp||'')+'</div>'
  +'<table><tr><th>帯域</th><th>捕捉</th></tr>'+brow+'</table>'
  +'<div class=\"rp-sec\"><b>時系列(24h)</b>'+(c.ts||'')+'</div>'
  +wkbHtml+corrHtml+hepHtml+sumHtml;
 document.getElementById('rpbody').innerHTML=body;
 document.getElementById('rpanel').classList.add('open');}
function closeReport(){document.getElementById('rpanel').classList.remove('open');curReportIdx=-1;}
var slider=document.getElementById('tslider');
slider.oninput=function(){draw(+this.value);};
var playing=false,timer=null;
document.getElementById('play').onclick=function(){
 playing=!playing;this.innerText=playing?'⏸ 停止':'▶ 再生';
 if(playing){timer=setInterval(function(){
   var v=(+slider.value+1)%(CELLS[0].frames.length);slider.value=v;draw(v);},900);}
 else{clearInterval(timer);}};
function toggleTheme(){document.body.classList.toggle('light');
 document.body.classList.toggle('dark');
 try{localStorage.setItem('dpp_theme',document.body.className);}catch(e){}}
try{var th=localStorage.getItem('dpp_theme');if(th)document.body.className=th;}catch(e){}
(function(){var b=document.getElementById('contourBtn');if(b)b.style.borderColor=CONTOUR_ON?'#5ec962':'';
 var lg=document.getElementById('contourLegend');if(lg)lg.style.display=CONTOUR_ON?'':'none';})();
var startFr=0;try{var sv=localStorage.getItem('dpp_fr');if(sv!==null)startFr=Math.min(+sv,CELLS[0].frames.length-1);}catch(e){}
slider.value=startFr;draw(startFr);
var NEXT=new Date(new Date().getTime()+__REFRESH__*1000);
function tick(){var s=Math.max(0,(NEXT-new Date())/1000);
 var m=String(Math.floor(s/60)).padStart(2,'0'),c=String(Math.floor(s%60)).padStart(2,'0');
 document.getElementById('cd').innerText='次回更新 '+m+':'+c;}
setInterval(tick,1000);tick();
// ===== 設定パネル =====
var CH_LABEL={interval_sec:function(v){return v+'秒 ('+(v/60).toFixed(v%60?1:0)+'分)';},
 total_minutes:function(v){return v===0?'無制限':v+'分';},
 weather_cache_sec:function(v){return (v/60)+'分';},
 seconds_back:function(v){return v+'秒 ('+(v/60)+'分)';}};
function fillSelect(id,key,choices,cur){var el=document.getElementById(id);
 el.innerHTML='';choices.forEach(function(v){var o=document.createElement('option');
  o.value=v;o.text=CH_LABEL[key]?CH_LABEL[key](v):v;
  if(+v===+cur)o.selected=true;el.appendChild(o);});}
function openSettings(){fetch('/api/duct/settings').then(function(r){return r.json();})
 .then(function(d){var s=d.settings,c=d.choices;
  fillSelect('s_interval','interval_sec',c.interval_sec,s.interval_sec);
  fillSelect('s_total','total_minutes',c.total_minutes,s.total_minutes);
  fillSelect('s_wx','weather_cache_sec',c.weather_cache_sec,s.weather_cache_sec);
  fillSelect('s_back','seconds_back',c.seconds_back,s.seconds_back);
  fillSelect('s_alert_dpp','alert_dpp_threshold',c.alert_dpp_threshold,s.alert_dpp_threshold);
  fillSelect('s_alert_radius','alert_radius_km',c.alert_radius_km,s.alert_radius_km);
  document.getElementById('s_alert_on').checked=!!s.alert_enabled;
  document.getElementById('s_webhook').value=s.discord_webhook_url||'';
  refreshStatus();document.getElementById('ovl').style.display='block';})
 .catch(function(e){alert('設定取得失敗: '+e);});}
function closeSettings(){document.getElementById('ovl').style.display='none';}
function fmtEta(h){if(h===null||h===undefined)return '推定中…';
 if(h<=0)return '到達済み';if(h<1)return '約'+Math.round(h*60)+'分';
 return '約'+h.toFixed(1)+'時間';}
function applySettings(){var body={
  interval_sec:+document.getElementById('s_interval').value,
  total_minutes:+document.getElementById('s_total').value,
  weather_cache_sec:+document.getElementById('s_wx').value,
  seconds_back:+document.getElementById('s_back').value,
  alert_enabled:document.getElementById('s_alert_on').checked,
  discord_webhook_url:document.getElementById('s_webhook').value,
  alert_dpp_threshold:+document.getElementById('s_alert_dpp').value,
  alert_radius_km:+document.getElementById('s_alert_radius').value};
 fetch('/api/duct/settings',{method:'POST',headers:{'Content-Type':'application/json'},
  body:JSON.stringify(body)}).then(function(r){return r.json();})
 .then(function(d){var m=document.getElementById('s_msg');
  m.innerText=d.recompute?'✓ 保存しました。即時再計算を開始':'✓ 保存しました(停止中は保留)';
  setTimeout(function(){m.innerText='';},6000);refreshStatus();})
 .catch(function(e){document.getElementById('s_msg').innerText='保存失敗: '+e;});}
var cpConfig=null;
function openCoordModal(cfg){
 cpConfig=cfg;
 document.getElementById('cp_title').innerText=cfg.title;
 var html='';
 cfg.fields.forEach(function(f){
  html+='<label>'+f.label
   +'<span class="cprow"><input type="text" id="'+f.id+'" placeholder="'+f.placeholder+'">'
   +'<button type="button" class="gpsbtn" onclick="fillGPS(\''+f.id+'\')" title="現在地(GPS)から自動取得">📍 現在地</button></span></label>';});
 if(cfg.heightField){var hf=cfg.heightField;
  html+='<label>'+hf.label
   +'<span class="cprow"><input type="text" id="'+hf.id+'" placeholder="'+hf.placeholder+'"></span></label>';}
 document.getElementById('cp_body').innerHTML=html;
 document.getElementById('cp_msg').innerText='';
 document.getElementById('covl').style.display='block';}
function closeCoordModal(){document.getElementById('covl').style.display='none';}
function fillGPS(inputId){
 var msg=document.getElementById('cp_msg');
 if(!navigator.geolocation){msg.innerText='⚠ このブラウザはGPS(位置情報)に対応していません';return;}
 msg.innerText='📍 現在地を取得中…';
 navigator.geolocation.getCurrentPosition(function(pos){
   var lat=pos.coords.latitude, lon=pos.coords.longitude;
   document.getElementById(inputId).value=lat.toFixed(5)+','+lon.toFixed(5);
   msg.innerText='✓ 現在地を取得しました(精度 約'+Math.round(pos.coords.accuracy)+'m)';},
  function(err){
   var reason=err.code===1?'位置情報の利用が許可されていません(ブラウザ/OS設定を確認)'
    :err.code===2?'現在地を取得できませんでした':'タイムアウトしました';
   msg.innerText='⚠ 取得失敗: '+reason;},
  {enableHighAccuracy:true,timeout:10000,maximumAge:60000});}
function looksLikeLatLon(v){
 var p=v.split(',').map(function(s){return s.trim();});
 if(p.length!==2)return false;
 var a=parseFloat(p[0]),b=parseFloat(p[1]);
 return isFinite(a)&&isFinite(b)&&Math.abs(a)<=90&&Math.abs(b)<=180;}
function geocodeResolve(q){
 return fetch('/api/duct/geocode?q='+encodeURIComponent(q))
  .then(function(r){return r.json();})
  .then(function(d){
   if(!d.ok)throw new Error('「'+q+'」が見つかりませんでした('+(d.error||'not_found')+')');
   return d.lat.toFixed(5)+','+d.lon.toFixed(5);});}
function cp_submit(){
 if(!cpConfig)return;
 var msg=document.getElementById('cp_msg');
 var raw={},ok=true;
 cpConfig.fields.forEach(function(f){var v=document.getElementById(f.id).value.trim();
  if(!v)ok=false;raw[f.id]=v;});
 if(!ok){msg.innerText='⚠ 座標または市町村名を入力するか📍現在地で取得してください';return;}
 if(cpConfig.heightField)raw[cpConfig.heightField.id]=document.getElementById(cpConfig.heightField.id).value.trim();
 msg.innerText='⌛ 座標を解決中…';
 var jobs=cpConfig.fields.map(function(f){
  var v=raw[f.id];
  if(looksLikeLatLon(v))return Promise.resolve();
  return geocodeResolve(v).then(function(ll){raw[f.id]=ll;});});
 Promise.all(jobs).then(function(){
  closeCoordModal();cpConfig.onSubmit(raw);
 }).catch(function(e){msg.innerText='⚠ '+e.message;});}
function apiErrMsg(code){
 if(code==='no_data'||code==='no_duct')return '計算データがまだありません(サーバー起動直後で初回計算中の可能性)。1〜2分待ってから再度お試しください。';
 if(code==='not_found')return '地名が見つかりませんでした。より具体的な市町村名を試すか、座標を直接入力してください。';
 return code;}
function askPass(){
 openCoordModal({title:'🛤️ Pass診断(自局⇔相手局)',
  fields:[{id:'cp_a',label:'自局の緯度,経度 または 市町村名',placeholder:'例: 35.68,139.77 または 東京都千代田区'},
   {id:'cp_b',label:'相手局の緯度,経度 または 市町村名',placeholder:'例: 33.59,130.40 または 福岡市'}],
  onSubmit:function(v){
   var pa=v.cp_a.split(',').map(parseFloat),pb=v.cp_b.split(',').map(parseFloat);
   fetch('/api/duct/pass?lat1='+pa[0]+'&lon1='+pa[1]+'&lat2='+pb[0]+'&lon2='+pb[1])
    .then(function(r){return r.json();}).then(function(d){
     if(!d.ok){alert('診断失敗: '+apiErrMsg(d.error));return;}
     alert('🛤️ Pass診断(自局⇔相手局)\n距離: '+d.distance_km+'km\n平均DPP: '+d.avg_dpp
      +'\n最小DPP: '+d.min_dpp+' / 最大DPP: '+d.max_dpp
      +'\n途切れリスク: '+d.break_risk_pct+'%');
    }).catch(function(e){alert('通信エラー: '+e);});}});}
function askAntenna(){
 openCoordModal({title:'📶 アンテナ高提案(自局単独)',
  fields:[{id:'cp_a',label:'自局の緯度,経度 または 市町村名',placeholder:'例: 38.27,140.87 または 仙台市'}],
  heightField:{id:'cp_h',label:'現在のアンテナ高[m](空欄で既定値'+ANTENNA_DEFAULT_H+'mを使用)',placeholder:String(ANTENNA_DEFAULT_H)},
  onSubmit:function(v){
   var p=v.cp_a.split(',').map(parseFloat);
   var h=(v.cp_h!=null && v.cp_h!=='') ? parseFloat(v.cp_h) : null;
   var url='/api/duct/antenna?lat='+p[0]+'&lon='+p[1]+(h!=null && !isNaN(h) ? '&height='+h : '');
   fetch(url)
    .then(function(r){return r.json();}).then(function(d){
     if(!d.ok){alert('取得失敗: '+apiErrMsg(d.error||'ダクト無し'));return;}
     var msg='📶 アンテナ高提案(自局単独)\nダクト底面: '+d.duct_base_m+'m (上端 '+d.duct_top_m+'m)\n'
      +'現在('+d.current_height_m+'m)の結合効率: '+(d.current_efficiency*100).toFixed(0)+'%';
     msg+=d.suggested_height_m!=null
      ?('\n効率70%以上に必要な高さ: 約'+d.suggested_height_m+'m'
        +(d.height_gap_m>0?(' (あと約'+d.height_gap_m+'m)'):' (既に到達)'))
      :'\nこのダクト構造では70%到達見込みなし';
     alert(msg);
    }).catch(function(e){alert('通信エラー: '+e);});}});}
var radarLayer=null;
function destPoint(lat,lon,bearingDeg,distKm){
 var R=6371.0088,brng=bearingDeg*Math.PI/180;
 var lat1=lat*Math.PI/180,lon1=lon*Math.PI/180;
 var lat2=Math.asin(Math.sin(lat1)*Math.cos(distKm/R)+Math.cos(lat1)*Math.sin(distKm/R)*Math.cos(brng));
 var lon2=lon1+Math.atan2(Math.sin(brng)*Math.sin(distKm/R)*Math.cos(lat1),
  Math.cos(distKm/R)-Math.sin(lat1)*Math.sin(lat2));
 return [lat2*180/Math.PI,((lon2*180/Math.PI)+540)%360-180];}
function showRadarOnMap(d){
 if(radarLayer){map.removeLayer(radarLayer);radarLayer=null;}
 var group=L.layerGroup();
 var starIcon=L.divIcon({className:'radar-star',html:'★',iconSize:[24,24],iconAnchor:[12,12]});
 L.marker([d.lat,d.lon],{icon:starIcon,zIndexOffset:1000})
  .bindTooltip('🧭 移動運用地点 ('+d.lat.toFixed(4)+', '+d.lon.toFixed(4)+')',{className:'radar-tip'})
  .addTo(group);
 var ranges=d.radar.map(function(x){return x.max_range_km;});
 var maxR=Math.max.apply(null,ranges.concat([1]));
 d.radar.forEach(function(x){
  if(x.max_range_km<=0)return;
  var end=destPoint(d.lat,d.lon,x.azimuth,x.max_range_km);
  var col=viridis(x.end_dpp!=null?x.end_dpp:(x.max_range_km/maxR)*100);
  var colStr='rgb('+col[0]+','+col[1]+','+col[2]+')';
  L.polyline([[d.lat,d.lon],end],{color:colStr,weight:1.4,opacity:0.85})
   .bindTooltip(x.azimuth+'° / '+x.max_range_km+'km'+(x.end_dpp!=null?' (DPP '+x.end_dpp+')':''),
    {className:'radar-tip',sticky:true})
   .addTo(group);
  L.circleMarker(end,{radius:2.5,color:colStr,fillColor:colStr,fillOpacity:1,weight:1})
   .addTo(group);});
 group.addTo(map);radarLayer=group;
 map.panTo([d.lat,d.lon]);}
function askRadar(){
 openCoordModal({title:'🧭 移動運用レーダー(自局から全方位)',
  fields:[{id:'cp_a',label:'自局(現在地)の緯度,経度 または 市町村名',placeholder:'例: 38.27,140.87 または 仙台市'}],
  onSubmit:function(v){
   var p=v.cp_a.split(',').map(parseFloat);
   fetch('/api/duct/radar?lat='+p[0]+'&lon='+p[1])
    .then(function(r){return r.json();}).then(function(d){
     if(!d.ok){alert('取得失敗: '+apiErrMsg(d.error||'ダクト無し'));return;}
     showRadarOnMap(d);
     var lines=d.radar.filter(function(x){return x.max_range_km>0;})
      .sort(function(a,b){return b.max_range_km-a.max_range_km;}).slice(0,8)
      .map(function(x){return x.azimuth+'° → '+x.max_range_km+'km';});
     var msg='🧭 移動運用レーダー(自局から全方位)\n'+(lines.join('\n')||'有効なダクト方向なし');
     if(d.elevation_m!=null)msg+='\n\n標高: '+d.elevation_m+'m';
     if(d.duct_base_m!=null)msg+=' / ダクト底面: '+d.duct_base_m+'m';
     if(d.height_to_duct_m!=null)msg+=' / 差: '+d.height_to_duct_m+'m';
     alert(msg);
    }).catch(function(e){alert('通信エラー: '+e);});}});}
function recomputeNow(){fetch('/api/duct/recompute',{method:'POST'})
 .then(function(r){return r.json();}).then(function(d){
  var m=document.getElementById('s_msg');m.innerText='⟳ 再計算を要求しました';
  setTimeout(function(){m.innerText='';},6000);})
 .catch(function(e){document.getElementById('s_msg').innerText='要求失敗: '+e;});}
function refreshStatus(){fetch('/api/duct/status').then(function(r){return r.json();})
 .then(function(d){var L=d.learning||{};var mat=L.maturity_pct||0;
  document.getElementById('s_status').innerHTML=
   '学習成熟度: <b>'+mat+'%</b>（目標 平均信頼度 '+(L.kappa_target||0.7)+'）'
   +'<div class="bar2"><i style="width:'+Math.min(100,mat)+'%"></i></div>'
   +'実用まで: <b>'+fmtEta(L.eta_hours)+'</b><br>'
   +'平均信頼度: '+(L.mean_confidence||0)+' / 被覆 '+(L.covered_cells||0)+'セル'
   +' / cone '+(d.cone_width_deg||'—')+'°<br>'
   +'状態: '+(d.settings.paused?'⏯ 観測停止中':'⏯ 観測稼働中');
  var eb=document.getElementById('etaBadge');
  if(eb)eb.innerText='実用まで '+fmtEta(L.eta_hours);}).catch(function(){});}
var pauseBtn=document.getElementById('pauseBtn');
function syncPause(){fetch('/api/duct/status').then(function(r){return r.json();})
 .then(function(d){pauseBtn.innerText=d.settings.paused?'⏯ 観測停止中':'⏯ 観測稼働中';});}
pauseBtn.onclick=function(){fetch('/api/duct/status').then(function(r){return r.json();})
 .then(function(d){return fetch('/api/duct/settings',{method:'POST',
   headers:{'Content-Type':'application/json'},
   body:JSON.stringify({paused:!d.settings.paused})});})
 .then(function(){syncPause();refreshStatus();});};
syncPause();
// ETAバッジをライブ更新(30秒毎)
setInterval(refreshStatus,30000);
</script></body></html>"""


def _mprofile_svg(mprof, ducts, w=210, h=120):
    """後方互換: 単一時刻のM(z)ミニプロファイル(旧ポップアップ用)。"""
    if not mprof or len(mprof) < 2:
        return ""
    return _mprofile_svg_multi([{"t": "now", "pts": mprof}], ducts, w=w, h=h)


_MPROF_COLORS = ["#38bdf8", "#5ec962", "#fde725"]


def _mprofile_svg_multi(mprofs3, ducts_now, w=260, h=170):
    """現在/+3h/+6h のM(z)を重ね描画し、dM/dz<0(ダクト層)を陰影表示する。"""
    series = [s for s in (mprofs3 or []) if s.get("pts") and len(s["pts"]) >= 2]
    if not series:
        return ""
    all_h = [z for s in series for z, m in s["pts"]]
    all_m = [m for s in series for z, m in s["pts"]]
    hmin, hmax = min(all_h), max(all_h); mmin, mmax = min(all_m), max(all_m)
    pad_l, pad_r, pad_t, pad_b = 26, 6, 6, 16

    def X(m):
        return pad_l + (m - mmin) / max(mmax - mmin, 1e-6) * (w - pad_l - pad_r)

    def Y(z):
        return h - pad_b - (z - hmin) / max(hmax - hmin, 1e-6) * (h - pad_t - pad_b)

    bands = ""
    for d in ducts_now or []:
        y1, y2 = Y(d["top_m"]), Y(d["base_m"])
        c = "#5ec962" if d.get("kind") == "surface" else "#f4a261"
        bands += ('<rect x="%d" y="%.1f" width="%d" height="%.1f" fill="%s" opacity="0.18"/>'
                  % (pad_l, min(y1, y2), w - pad_l - pad_r, abs(y2 - y1), c))
    lines = ""; legend = ""
    for i, s in enumerate(series[:3]):
        col = _MPROF_COLORS[i % len(_MPROF_COLORS)]
        pts = " ".join("%.1f,%.1f" % (X(m), Y(z)) for z, m in s["pts"])
        lines += '<polyline points="%s" fill="none" stroke="%s" stroke-width="1.6"/>' % (pts, col)
        legend += ('<circle cx="%d" cy="%d" r="3" fill="%s"/>'
                   '<text x="%d" y="%d" fill="#cbd5e1" font-size="8">%s</text>'
                   % (pad_l + 4 + i * 70, h - 5, col, pad_l + 10 + i * 70, h - 2, s["t"]))
    axis = ('<text x="2" y="%d" fill="#64748b" font-size="7">%.0fm</text>'
            '<text x="2" y="%d" fill="#64748b" font-size="7">%.0fm</text>'
            % (Y(hmin), hmin, Y(hmax) + 6, hmax))
    return ('<svg width="%d" height="%d" style="background:#0b1220;border-radius:5px;margin-top:4px">'
            '%s%s%s%s<text x="%d" y="11" fill="#94a3b8" font-size="8">M(z) 現在/+3h/+6h(陰影=ダクト層)</text></svg>'
            % (w, h, bands, lines, legend, axis, pad_l))


def _duct_structure_svg(ducts, hmax_hint=None, w=110, h=170):
    """高度別ダクト構造を縦棒(柱)で直感表示。専門知識がなくても
    「接地/上空」「厚さ」「どの高さか」が一目でわかることを狙う。"""
    pad_t, pad_b, pad_l = 8, 18, 34
    col_w = 26
    hmax = hmax_hint or (max((d["top_m"] for d in ducts), default=500.0) * 1.3 if ducts else 500.0)
    hmax = max(hmax, 100.0)

    def Y(z):
        return h - pad_b - max(0.0, min(1.0, z / hmax)) * (h - pad_t - pad_b)

    ticks = ""
    for i in range(6):
        zt = hmax * i / 5.0
        yt = Y(zt)
        ticks += ('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#26344d" stroke-width="1"/>'
                  '<text x="2" y="%.1f" fill="#64748b" font-size="7">%.0fm</text>'
                  % (pad_l - 4, yt, pad_l + col_w + 4, yt, yt + 3, zt))
    col_x = pad_l
    body = ('<rect x="%d" y="%d" width="%d" height="%d" fill="#111c33" '
            'stroke="#26344d" stroke-width="1" rx="3"/>'
            % (col_x, pad_t, col_w, h - pad_t - pad_b))
    blocks = ""; labels = ""
    for i, d in enumerate(ducts or []):
        y1, y2 = Y(d["top_m"]), Y(d["base_m"])
        col = "#5ec962" if d.get("kind") == "surface" else "#f4a261"
        blocks += ('<rect x="%d" y="%.1f" width="%d" height="%.1f" fill="%s" opacity="0.85" rx="2"/>'
                   % (col_x, min(y1, y2), col_w, max(abs(y2 - y1), 2.0), col))
        if i < 2:
            kind_s = "接地" if d.get("kind") == "surface" else "上空"
            labels += ('<text x="%d" y="%.1f" fill="#cbd5e1" font-size="8">%s %.0f-%.0fm(厚%.0fm)</text>'
                       % (pad_l + col_w + 6, (min(y1, y2) + max(y1, y2)) / 2.0 + 3, kind_s,
                          d["base_m"], d["top_m"], d["thickness_m"]))
    ground = '<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#94a3b8" stroke-width="1.5"/>' \
             % (col_x - 4, Y(0), col_x + col_w + 4, Y(0))
    return ('<svg width="%d" height="%d" style="background:#0b1220;border-radius:5px;margin-top:4px">'
            '<text x="4" y="10" fill="#94a3b8" font-size="8">高度別ダクト構造</text>'
            '<g transform="translate(0,10)">%s%s%s%s%s</g></svg>'
            % (w, h + 14, ticks, body, blocks, ground, labels))


def _components_bar_svg(components, w=170, h=90):
    if not components:
        return ""
    order = [("structure", "構造", "#38bdf8"), ("stability", "安定性", "#5ec962"),
             ("coupling", "結合", "#facc15"), ("turbulence", "乱流*", "#f87171")]
    n = len(order); bw = (w - 10) / n; bars = ""
    for i, (key, label, col) in enumerate(order):
        v = max(0.0, min(100.0, components.get(key, 0.0)))
        bh = v / 100.0 * (h - 22)
        x = 5 + i * bw
        bars += ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" rx="2"/>'
                 '<text x="%.1f" y="%d" fill="#e6edf7" font-size="8" text-anchor="middle">%.0f</text>'
                 '<text x="%.1f" y="%d" fill="#94a3b8" font-size="7" text-anchor="middle">%s</text>'
                 % (x + bw * 0.15, h - 14 - bh, bw * 0.7, bh, col,
                    x + bw * 0.5, h - 16 - bh, v, x + bw * 0.5, h - 4, label))
    return ('<svg width="%d" height="%d" style="background:#0b1220;border-radius:5px">'
            '%s</svg><div style="font-size:9px;color:#64748b">*乱流は逆評価(小さいほど良好)</div>'
            % (w, h, bars))


def _timeseries_svg(frames, w=260, h=70):
    vals = [(i, f["dpp"]) for i, f in enumerate(frames) if f]
    if len(vals) < 2:
        return ""
    pad = 5
    xmax = len(frames) - 1
    def X(i): return pad + i / max(xmax, 1) * (w - 2 * pad)
    def Y(v): return h - pad - v / 100.0 * (h - 2 * pad)
    pts = " ".join("%.1f,%.1f" % (X(i), Y(v)) for i, v in vals)
    dots = "".join('<circle cx="%.1f" cy="%.1f" r="1.6" fill="#38bdf8"/>' % (X(i), Y(v))
                   for i, v in vals)
    return ('<svg width="%d" height="%d" style="background:#0b1220;border-radius:5px">'
            '<polyline points="%s" fill="none" stroke="#38bdf8" stroke-width="1.4"/>%s'
            '<text x="4" y="10" fill="#94a3b8" font-size="8">DPP推移(24h)</text></svg>'
            % (w, h, pts, dots))


def build_duct_celldata(results, region, cal):
    """ITPFS(整理版): generate_html()内にあった「セル一覧+学習成熟度」構築ロジックを
    共通関数として独立させた。旧実装ではこの処理はgenerate_html()の内部にしか
    無く、standaloneのHTMLファイル(OUTPUT_HTML)を経由しないと利用できなかった。
    ダッシュボード統合版(inline地図)向けのJSONエンドポイント(/api/duct/celldata)
    からも同じロジックを再利用できるようにするため、HTML生成部分と分離した。
    戻り値: (cells, frame_labels, lq, eta_s) または 有効地点が無い場合 (None, None, None, None)
    """
    valid = [r for r in results if r]
    if not valid:
        return None, None, None, None
    frame_labels = []
    for r in valid:
        frame_labels = [f["t"] if f else "--" for f in r["frames"]]
        break
    cells = []
    for r in valid:
        d0 = r["ducts0"]
        if d0:
            ds = " / ".join(
                "%s%.0f-%.0fm(ΔM%.0f)" % ("接地" if d["kind"] == "surface" else "上空",
                                          d["base_m"], d["top_m"], d["delta_M"])
                for d in d0[:2])
        else:
            ds = "ダクトなし"
        hmax_hint = max((z for z, m in r["mprof"]), default=500.0) if r["mprof"] else 500.0
        cells.append({"lat": r["lat"], "lon": r["lon"], "conf": r["confidence"],
                      "duct": ds, "svg": _mprofile_svg_multi(r.get("mprofs3"), d0, w=190, h=150),
                      "struct": _duct_structure_svg(d0, hmax_hint, w=100, h=150),
                      "comp": _components_bar_svg(r.get("components")),
                      "ts": _timeseries_svg(r["frames"]),
                      "wkb": r.get("wkb"), "corr": r.get("correction"),
                      "summary": r.get("summary", ""),
                      "frames": [({"dpp": f["dpp"], "r": f["r"], "trap": f["trap"]}
                                  if f else None) for f in r["frames"]]})
    lq = learning_quality(cal)
    eta = lq.get("eta_hours")
    if eta is None:
        eta_s = "推定中…"
    elif eta <= 0:
        eta_s = "到達済み"
    elif eta < 1:
        eta_s = f"約{round(eta*60)}分"
    else:
        eta_s = f"約{eta:.1f}時間"
    return cells, frame_labels, lq, eta_s


def generate_html(results, region, cal, refresh_sec, out_path):
    cells, frame_labels, lq, eta_s = build_duct_celldata(results, region, cal)
    if cells is None:
        print("[html] 有効地点なし")
        return None
    html = (_HTML_TEMPLATE
            .replace("__TITLE__", f"{SYSTEM_NAME} {VERSION_LABEL} | Duct Propagation Forecast (DPP engine build {DUCT_SOURCE_VERSION})")
            .replace("__VER__", DUCT_SOURCE_VERSION)
            .replace("__REFRESH__", str(int(refresh_sec)))
            .replace("__GEN__", datetime.now().strftime("%m/%d %H:%M:%S"))
            .replace("__NPTS__", str(len(cells)))
            .replace("__MATURITY__", str(lq["maturity_pct"]))
            .replace("__ETA__", eta_s)
            .replace("__COVERED__", str(lq["covered_cells"]))
            .replace("__MCONF__", str(lq["mean_confidence"]))
            .replace("__MAXFR__", str(len(frame_labels) - 1))
            .replace("__STEP__", str(region["step_deg"]))
            .replace("__TXHEIGHT__", str(TX_HEIGHT_M))
            .replace("__BANDLIST__", json.dumps([b.replace("MHz", "") for b in BANDS]))
            .replace("__REGION__", json.dumps(region))
            .replace("__FRAMET__", json.dumps(frame_labels, ensure_ascii=False)))
    html = html.replace("__CELLDATA__", json.dumps(cells, ensure_ascii=False))
    tmp = out_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(html)
    os.replace(tmp, out_path)
    return out_path


HEPBURN_REGION_CODE = "eas"
HEPBURN_BASE_URL = "https://www.dxinfocentre.com/tr_map/fcst/{code}{hh}.png"
HEPBURN_FRAME_HOURS_EAS = [6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36,
                           42, 48, 54, 60, 66, 72, 78, 84, 90, 96,
                           102, 108, 114, 120, 126, 132, 138, 144]
HEPBURN_LEGEND_LABELS = ["NIL", "MARGINAL", "FAIR", "MODERATE", "HIGH", "STRONG",
                         "VERY_STRONG", "INTENSE", "VERY_INTENSE", "EXTREME1",
                         "EXTREME2"]
HEPBURN_REGION_AXIS_LABELS = {
    "eas": {"lon": list(range(105, 156, 5)), "lat": list(range(51, 20, -3))},
}
HEPBURN_CACHE_DIR = os.path.join(_HERE, "hepburn_cache")
HEPBURN_CALIB_FILE = os.path.join(_HERE, "hepburn_calib.json")
HEPBURN_LOG_FILE = os.path.join(_HERE, "hepburn_dpp_compare.jsonl")
HEPBURN_FIT_FILE = os.path.join(_HERE, "hepburn_fit.json")
HEPBURN_CACHE_MAX_AGE_SEC = 3 * 3600
HEPBURN_MIN_SAMPLES_FOR_CORRECTION = 30
HEPBURN_CORRECTION_KAPPA_K = 40
HEPBURN_CORRECTION_LOOKBACK_DAYS = 14

os.makedirs(HEPBURN_CACHE_DIR, exist_ok=True)


def hepburn_run_base_utc(now=None):
    """地図は毎日1800UTC頃更新。直近のrun基準時刻(UTC)を返す。"""
    now = now or datetime.now(timezone.utc)
    run = now.replace(hour=18, minute=0, second=0, microsecond=0)
    if now < run:
        run -= timedelta(days=1)
    return run


def hepburn_nearest_frame_hour(target_utc=None):
    """target_utc に最も近いフレームのオフセット時間(h)とrun基準時刻を返す。"""
    run = hepburn_run_base_utc(target_utc)
    target_utc = target_utc or datetime.now(timezone.utc)
    delta_h = (target_utc - run).total_seconds() / 3600.0
    best = min(HEPBURN_FRAME_HOURS_EAS, key=lambda h: abs(h - delta_h))
    return best, run


def hepburn_fetch_map(hh, force=False, timeout=20):
    if _HepImage is None:
        raise RuntimeError("Pillow が必要です: pip install pillow --break-system-packages")
    fname = os.path.join(HEPBURN_CACHE_DIR, f"{HEPBURN_REGION_CODE}{hh:03d}.png")
    stale = (not os.path.exists(fname)) or \
            (time.time() - os.path.getmtime(fname) > HEPBURN_CACHE_MAX_AGE_SEC)
    if force or stale:
        url = HEPBURN_BASE_URL.format(code=HEPBURN_REGION_CODE, hh=f"{hh:03d}") + f"?v{int(time.time())}"
        req = Request(url, headers={"User-Agent": "Mozilla/5.0 (personal hobby use; DPP compare tool)"})
        with urlopen(req, timeout=timeout) as r:
            data = r.read()
        with open(fname, "wb") as f:
            f.write(data)
    return _HepImage.open(fname).convert("RGB")


def hepburn_load_calib():
    if os.path.exists(HEPBURN_CALIB_FILE):
        with open(HEPBURN_CALIB_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    return None


def hepburn_save_calib(calib):
    with open(HEPBURN_CALIB_FILE, "w", encoding="utf-8") as f:
        json.dump(calib, f, ensure_ascii=False, indent=1)


def _hepburn_detect_border_box(arr):
    """地図プロット枠(白い矩形)のピクセル座標 (x0,y0,x1,y1) を検出する。"""
    h, w, _ = arr.shape
    white = (arr[:, :, 0] > 230) & (arr[:, :, 1] > 230) & (arr[:, :, 2] > 230)
    row_counts = white.sum(axis=1)
    col_counts = white.sum(axis=0)
    row_thresh = 0.8 * w
    col_thresh = 0.8 * h
    rows = np.where(row_counts > row_thresh)[0]
    cols = np.where(col_counts > col_thresh)[0]
    if len(rows) < 2 or len(cols) < 2:
        raise RuntimeError("プロット枠(白い矩形の境界線)を検出できませんでした。"
                            "画像が想定と違う形式かもしれません。")
    y0, y1 = int(rows.min()), int(rows.max())
    x0, x1 = int(cols.min()), int(cols.max())
    return x0, y0, x1, y1


def _hepburn_cluster_1d(xs, n_expected, gap_thresh=5):
    """1次元座標列を、まず隙間で粗くグループ化し、期待クラスタ数になるまで
    最も近い隣接グループ同士を併合していく(文字幅のばらつき対策)。"""
    xs = sorted(xs)
    groups = [[xs[0]]]
    for x in xs[1:]:
        if x - groups[-1][-1] > gap_thresh:
            groups.append([x])
        else:
            groups[-1].append(x)
    while len(groups) > n_expected:
        gaps = [groups[i + 1][0] - groups[i][-1] for i in range(len(groups) - 1)]
        i = int(np.argmin(gaps))
        groups[i] = groups[i] + groups[i + 1]
        del groups[i + 1]
    if len(groups) != n_expected:
        raise RuntimeError(f"軸ラベルの検出数が想定と合いません: 検出{len(groups)} "
                            f"/ 期待{n_expected}。画像を確認してください。")
    return [sum(g) / len(g) for g in groups]


def hepburn_calibrate_geo_from_image(image_path, region_code="eas", lon_labels=None,
                                      lat_labels=None, save=True):
    """
    画像に焼き込まれた軸ラベル(105E...155E / 51N...21N など)を自動検出し、
    手動でピクセル座標を読み取ることなく緯度経度較正を行う。
    """
    if np is None or _HepImage is None:
        raise RuntimeError("numpy と Pillow が必要です")
    cfg = HEPBURN_REGION_AXIS_LABELS.get(region_code, {})
    lon_labels = lon_labels or cfg.get("lon")
    lat_labels = lat_labels or cfg.get("lat")
    if not lon_labels or not lat_labels:
        raise ValueError(f"region_code='{region_code}' の軸ラベル定義がありません。"
                          f"lon_labels/lat_labels を明示的に渡してください。")

    img = _HepImage.open(image_path).convert("RGB")
    arr = np.array(img)
    x0, y0, x1, y1 = _hepburn_detect_border_box(arr)

    strip = arr[y1 + 1:min(y1 + 30, arr.shape[0]), :, :]
    white = (strip[:, :, 0] > 200) & (strip[:, :, 1] > 200) & (strip[:, :, 2] > 200)
    xs = np.where(white.any(axis=0))[0]
    if len(xs) == 0:
        raise RuntimeError("経度ラベル行を検出できませんでした")
    x_centers = _hepburn_cluster_1d(list(xs), len(lon_labels))

    strip2 = arr[:, max(0, x0 - 40):x0, :]
    white2 = (strip2[:, :, 0] > 200) & (strip2[:, :, 1] > 200) & (strip2[:, :, 2] > 200)
    ys = np.where(white2.any(axis=1))[0]
    if len(ys) == 0:
        raise RuntimeError("緯度ラベル列を検出できませんでした")
    y_centers = _hepburn_cluster_1d(list(ys), len(lat_labels))

    A_lon = np.vstack([lon_labels, np.ones(len(lon_labels))]).T
    mx, bx = np.linalg.lstsq(A_lon, x_centers, rcond=None)[0]
    A_lat = np.vstack([lat_labels, np.ones(len(lat_labels))]).T
    my, by = np.linalg.lstsq(A_lat, y_centers, rcond=None)[0]

    resid_x = [round(abs((mx * lo + bx) - xc), 2) for lo, xc in zip(lon_labels, x_centers)]
    resid_y = [round(abs((my * la + by) - yc), 2) for la, yc in zip(lat_labels, y_centers)]

    calib = hepburn_load_calib() or {}
    calib.update({
        "region_code": region_code,
        "coefx": [float(mx), 0.0, float(bx)],
        "coefy": [0.0, float(my), float(by)],
        "border_box": [x0, y0, x1, y1],
        "geo_auto_calibrated_at": datetime.now().isoformat(timespec="seconds"),
        "geo_residual_px": {"lon": resid_x, "lat": resid_y},
        "image_size": list(img.size),
    })
    if save:
        hepburn_save_calib(calib)
    return calib


def hepburn_calibrate_geo(points, existing=None):
    """
    points: [{"lat":.., "lon":.., "px":.., "py":..}, ...] 最低3点(4〜5点推奨)
    手動でピクセル座標を与える版(calibrate_geo_from_imageが使えない場合用)。
    """
    if np is None:
        raise RuntimeError("numpy が必要です")
    if len(points) < 3:
        raise ValueError("基準点は最低3点必要です(4〜5点推奨)")
    A = np.array([[p["lon"], p["lat"], 1.0] for p in points])
    bx = np.array([p["px"] for p in points])
    by = np.array([p["py"] for p in points])
    coefx, *_ = np.linalg.lstsq(A, bx, rcond=None)
    coefy, *_ = np.linalg.lstsq(A, by, rcond=None)
    resid = []
    for p in points:
        px = coefx[0] * p["lon"] + coefx[1] * p["lat"] + coefx[2]
        py = coefy[0] * p["lon"] + coefy[1] * p["lat"] + coefy[2]
        resid.append(round(math.hypot(px - p["px"], py - p["py"]), 2))
    calib = existing or {}
    calib.update({
        "coefx": coefx.tolist(), "coefy": coefy.tolist(),
        "geo_points": points, "geo_residual_px": resid,
        "geo_calibrated_at": datetime.now().isoformat(timespec="seconds"),
    })
    hepburn_save_calib(calib)
    return calib


def hepburn_latlon_to_pixel(calib, lat, lon):
    cx, cy = calib["coefx"], calib["coefy"]
    px = cx[0] * lon + cx[1] * lat + cx[2]
    py = cy[0] * lon + cy[1] * lat + cy[2]
    return px, py


def hepburn_record_legend(swatches, existing_calib=None):
    """swatches: [{"label":"NIL","px":..,"py":..}, ...] 11点"""
    hh, _ = hepburn_nearest_frame_hour()
    img = hepburn_fetch_map(hh)
    legend_colors = {}
    for sw in swatches:
        rgb = hepburn_sample_rgb(img, sw["px"], sw["py"])
        legend_colors[sw["label"]] = rgb
    calib = existing_calib or hepburn_load_calib() or {}
    calib["legend_colors"] = legend_colors
    calib["legend_calibrated_at"] = datetime.now().isoformat(timespec="seconds")
    hepburn_save_calib(calib)
    return legend_colors


def hepburn_record_legend_manual(rgb_list, existing_calib=None):
    """rgb_list: [[r,g,b], ...] HEPBURN_LEGEND_LABELS と同じ順で11個。"""
    if len(rgb_list) != len(HEPBURN_LEGEND_LABELS):
        raise ValueError(f"{len(HEPBURN_LEGEND_LABELS)}色分(0〜10+)を渡してください")
    calib = existing_calib or hepburn_load_calib() or {}
    calib["legend_colors"] = dict(zip(HEPBURN_LEGEND_LABELS, [list(map(int, c)) for c in rgb_list]))
    calib["legend_calibrated_at"] = datetime.now().isoformat(timespec="seconds")
    hepburn_save_calib(calib)
    return calib["legend_colors"]


def hepburn_record_legend_linear(top_pt, bottom_pt, existing_calib=None, labels=None):
    """凡例が等間隔の直線状カラーバー(0=NILの中心〜10+=EXTREMEの中心)である場合、
    上端・下端2点の中心ピクセル座標だけから11段階ぶんの座標を等分補間して
    自動サンプリングする(11点を1つずつ手打ちしなくて済む簡易版)。
    top_pt/bottom_pt: (px, py) のタプル。
    """
    labels = labels or HEPBURN_LEGEND_LABELS
    n = len(labels)
    swatches = []
    for i, label in enumerate(labels):
        t = i / (n - 1)
        px = top_pt[0] + (bottom_pt[0] - top_pt[0]) * t
        py = top_pt[1] + (bottom_pt[1] - top_pt[1]) * t
        swatches.append({"label": label, "px": px, "py": py})
    return hepburn_record_legend(swatches, existing_calib)


def hepburn_save_debug_grid(image_path=None, hh=None, out_path=None, step=25):
    """凡例スウォッチや軸ラベルのピクセル座標を目視で読み取れるよう、
    グリッド線と座標数値を書き込んだ確認用画像を保存する(較正作業の補助)。
    image_path未指定時はhh(既定: 直近フレーム)の地図を取得して使う。
    """
    if _HepImage is None:
        raise RuntimeError("Pillow が必要です: pip install pillow --break-system-packages")
    from PIL import ImageDraw
    if image_path:
        img = _HepImage.open(image_path).convert("RGB")
    else:
        hh = hh if hh is not None else hepburn_nearest_frame_hour()[0]
        img = hepburn_fetch_map(hh)
    img2 = img.copy()
    draw = ImageDraw.Draw(img2)
    w, h = img2.size
    for x in range(0, w, step):
        col = (255, 0, 0) if x % (step * 4) == 0 else (255, 120, 120)
        draw.line([(x, 0), (x, h)], fill=col, width=1)
        if x % (step * 4) == 0:
            draw.text((x + 2, 2), str(x), fill=(255, 255, 0))
    for y in range(0, h, step):
        col = (255, 0, 0) if y % (step * 4) == 0 else (255, 120, 120)
        draw.line([(0, y), (w, y)], fill=col, width=1)
        if y % (step * 4) == 0:
            draw.text((2, y + 2), str(y), fill=(255, 255, 0))
    out_path = out_path or os.path.join(_HERE, "hepburn_grid_debug.png")
    img2.save(out_path)
    return out_path


def hepburn_sample_rgb(img, px, py, radius=1):
    """(px,py)周辺 radius ピクセルの平均RGBを返す(アンチエイリアス対策)。"""
    w, h = img.size
    xs0, xs1 = max(0, int(px) - radius), min(w, int(px) + radius + 1)
    ys0, ys1 = max(0, int(py) - radius), min(h, int(py) + radius + 1)
    pix = img.load()
    rs, gs, bs, n = 0, 0, 0, 0
    for y in range(ys0, ys1):
        for x in range(xs0, xs1):
            r, g, b = pix[x, y]
            rs += r; gs += g; bs += b; n += 1
    if n == 0:
        return [0, 0, 0]
    return [round(rs / n), round(gs / n), round(bs / n)]


def hepburn_classify_intensity(rgb, legend_colors):
    """凡例色との最近傍(ユークリッド距離)でラベル/0-10のindexを返す。"""
    if not legend_colors:
        return None, None
    best_label, best_d = None, float("inf")
    for label, ref in legend_colors.items():
        d = sum((a - b) ** 2 for a, b in zip(rgb, ref)) ** 0.5
        if d < best_d:
            best_d, best_label = d, label
    idx = HEPBURN_LEGEND_LABELS.index(best_label) if best_label in HEPBURN_LEGEND_LABELS else None
    return best_label, idx


def hepburn_sample_grid(calib, points, hh=None, target_utc=None):
    """
    points: [(lat,lon), ...]
    画像取得を1回だけ行い、複数地点を一括サンプリングする
    (Hepburn Mapという「1枚の面情報」を最大限に活かすための関数)。
    """
    if hh is None:
        hh, run = hepburn_nearest_frame_hour(target_utc)
    else:
        run = hepburn_run_base_utc(target_utc)
    img = hepburn_fetch_map(hh)
    legend = calib.get("legend_colors")
    out = []
    for lat, lon in points:
        px, py = hepburn_latlon_to_pixel(calib, lat, lon)
        if not (0 <= px < img.size[0] and 0 <= py < img.size[1]):
            continue
        rgb = hepburn_sample_rgb(img, px, py)
        label, idx = hepburn_classify_intensity(rgb, legend)
        out.append({"lat": lat, "lon": lon, "rgb": rgb, "label": label,
                     "intensity": idx, "hepburn_frame_offset_h": hh,
                     "hepburn_valid_utc": (run + timedelta(hours=hh)).isoformat(timespec="seconds")})
    return out


def hepburn_append_batch(records):
    """複数のログ行を1回のファイルI/Oでまとめて追記する。"""
    if not records:
        return
    with open(HEPBURN_LOG_FILE, "a", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def hepburn_load_log(limit=None):
    if not os.path.exists(HEPBURN_LOG_FILE):
        return []
    with open(HEPBURN_LOG_FILE, "r", encoding="utf-8") as f:
        lines = [json.loads(l) for l in f if l.strip()]
    return lines[-limit:] if limit else lines


def hepburn_fit_correction(log=None, min_n=HEPBURN_MIN_SAMPLES_FOR_CORRECTION,
                            lookback_days=HEPBURN_CORRECTION_LOOKBACK_DAYS):
    """
    蓄積ログから (Hepburn強度 -> DPP) のOLS回帰を行い、補正パラメータを返す。
    直近lookback_days以内のサンプルのみを使う。標本数が min_n 未満なら無効。
    """
    if np is None:
        return {"applied": False, "n": 0, "reason": "numpy not available"}
    log = hepburn_load_log() if log is None else log
    cutoff = datetime.now(timezone.utc) - timedelta(days=lookback_days)
    xs, ys = [], []
    for rec in log:
        h = rec.get("hepburn_intensity_0_10")
        dpp_field = rec.get("dpp")
        d = dpp_field.get("dpp") if isinstance(dpp_field, dict) else dpp_field
        if h is None or d is None:
            continue
        try:
            t = datetime.fromisoformat(rec.get("logged_at_local", ""))
        except ValueError:
            t = None
        if t is not None and t.tzinfo is None:
            t = t.replace(tzinfo=timezone.utc)
        if t is not None and t < cutoff:
            continue
        xs.append(float(h)); ys.append(float(d))
    n = len(xs)
    if n < min_n:
        return {"applied": False, "n": n, "min_required": min_n}
    X = np.vstack([xs, np.ones(n)]).T
    y = np.array(ys)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    slope, intercept = float(coef[0]), float(coef[1])
    pred = X @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2)) or 1e-9
    r2 = round(1.0 - ss_res / ss_tot, 3)
    kappa = round(n / (n + HEPBURN_CORRECTION_KAPPA_K), 3)
    return {"applied": True, "n": n, "slope": round(slope, 4), "intercept": round(intercept, 2),
            "r2": r2, "kappa": kappa,
            "fitted_at": datetime.now().isoformat(timespec="seconds")}


def hepburn_save_correction_fit(fit):
    with open(HEPBURN_FIT_FILE, "w", encoding="utf-8") as f:
        json.dump(fit, f, ensure_ascii=False, indent=1)
    return fit


def hepburn_load_correction_fit():
    if os.path.exists(HEPBURN_FIT_FILE):
        with open(HEPBURN_FIT_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    return {"applied": False, "n": 0}


def hepburn_blend_correction(dpp_value, hepburn_idx, fit=None):
    """
    fit(hepburn_fit_correction()の結果)を使い、DPP生値をHepburnの期待値
    方向に kappa比率でブレンドする。fit未適用 or hepburn_idx欠損時は無補正。
    """
    fit = fit or hepburn_load_correction_fit()
    if not fit.get("applied") or hepburn_idx is None:
        return {"applied": False, "pre_dpp": dpp_value, "final_dpp": dpp_value,
                "n": fit.get("n", 0), "min_required": fit.get("min_required")}
    expected = fit["slope"] * hepburn_idx + fit["intercept"]
    kappa = fit["kappa"]
    final = dpp_value * (1 - kappa) + expected * kappa
    final = max(0.0, min(100.0, final))
    return {"applied": True, "pre_dpp": round(dpp_value, 1),
            "expected_from_hepburn": round(expected, 1), "kappa": kappa,
            "n": fit["n"], "r2": fit.get("r2"), "final_dpp": round(final, 1)}


def hepburn_compare_now(dpp_value, lat, lon, extra=None, force_fetch=False, target_utc=None):
    """
    直近フレームのHepburn Mapを取得し、(lat,lon)地点の強度を分類して
    DPP側の値と一緒にJSONLへ追記する。戻り値は記録したレコード。
    """
    calib = hepburn_load_calib()
    if not calib or "coefx" not in calib:
        raise RuntimeError("先に hepburn_calibrate_geo_from_image() で較正してください")
    hh, run = hepburn_nearest_frame_hour(target_utc)
    img = hepburn_fetch_map(hh, force=force_fetch)
    px, py = hepburn_latlon_to_pixel(calib, lat, lon)
    if not (0 <= px < img.size[0] and 0 <= py < img.size[1]):
        raise ValueError(f"較正結果が画像範囲外です px={px:.1f},py={py:.1f} "
                          f"(image={img.size}) -- 較正点を見直してください")
    rgb = hepburn_sample_rgb(img, px, py)
    label, idx = hepburn_classify_intensity(rgb, calib.get("legend_colors"))
    rec = {
        "logged_at_local": datetime.now().isoformat(timespec="seconds"),
        "hepburn_run_utc": run.isoformat(timespec="seconds"),
        "hepburn_frame_offset_h": hh,
        "hepburn_valid_utc": (run + timedelta(hours=hh)).isoformat(timespec="seconds"),
        "lat": lat, "lon": lon, "pixel": [round(px, 1), round(py, 1)],
        "hepburn_rgb": rgb, "hepburn_label": label, "hepburn_intensity_0_10": idx,
        "dpp": dpp_value,
    }
    if extra:
        rec["extra"] = extra
    with open(HEPBURN_LOG_FILE, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return rec


def hepburn_lookup(lat, lon):
    """任意の緯度経度でHepburn Mapを即時照会する(ログには書かない)。
    地図クリックによるオンデマンド照会用。較正済みであれば場所を選ばず
    日本全域(Far Eastマップの範囲内)どこでも呼べる。"""
    if _HepImage is None:
        return {"ok": False, "error": "Pillow not installed"}
    calib = hepburn_load_calib()
    if not calib or "coefx" not in calib:
        return {"ok": False, "error": "not_calibrated (geo)"}
    try:
        hh, run = hepburn_nearest_frame_hour()
        img = hepburn_fetch_map(hh)
        px, py = hepburn_latlon_to_pixel(calib, lat, lon)
        if not (0 <= px < img.size[0] and 0 <= py < img.size[1]):
            return {"ok": False, "error": "out_of_map_bounds"}
        rgb = hepburn_sample_rgb(img, px, py)
        label, idx = hepburn_classify_intensity(rgb, calib.get("legend_colors"))
        return {"ok": True, "lat": lat, "lon": lon, "rgb": rgb,
                "label": label, "intensity_0_10": idx,
                "valid_utc": (run + timedelta(hours=hh)).isoformat(timespec="seconds"),
                "frame_offset_h": hh}
    except Exception as e:
        return {"ok": False, "error": str(e)}


def _nearest_result(results, lat, lon):
    best, best_d = None, float("inf")
    for r in results:
        if not r:
            continue
        d = (r["lat"] - lat) ** 2 + (r["lon"] - lon) ** 2
        if d < best_d:
            best_d, best = d, r
    return best


def _hepburn_apply_correction(results):
    """
    Hepburn Tropo Mapを使ったDPP補正(Ver4.0)。

    設計方針:
      1) Web上の情報(1枚のHepburn Map画像)を最大限活用するため、HOME地点
         1点だけでなく「今回の全DPPグリッド点」を1回の画像取得で一括
         サンプリングする(sample_grid)。
      2) サンプルを蓄積ログに追記し、直近14日分の蓄積から
         (Hepburn強度→DPP)のOLS回帰を毎サイクル再フィットする
         (fit_correction: R^2を常時算出し説明力を可視化)。
      3) サンプル数が閾値(MIN_SAMPLES_FOR_CORRECTION)未満のうちは
         補正を無効化し、生DPP(PSKReporter補正後)をそのまま使う
         (過学習・早すぎる外挿を避ける理論的な安全策)。
      4) 有効化後も kappa=n/(n+K) によるベイズ的シュリンケージで
         「DPP」と「Hepburnから期待される値」を段階的にブレンドする
         (既存のPSKReporter補正と同一思想でアーキテクチャを統一)。
      5) 結果は r["correction"]["hepburn"] に透明性を保って記録し、
         最終的な表示値 r["frames"][0]["dpp"] にのみ反映する
         (raw_dpp/PSK補正後の値はcorrection内に保持し、検証可能にする)。

    未較正・モジュール未導入・ネットワーク不可のいずれでも本体の動作は
    止めない(失敗時は無補正のまま継続)。
    """
    if not (HEPBURN_COMPARE_ENABLED and _HepImage is not None):
        return
    calib = hepburn_load_calib()
    if not calib or "coefx" not in calib or "legend_colors" not in calib:
        return

    valid = [r for r in results if r and r["frames"] and r["frames"][0]]
    if not valid:
        return

    try:
        points = [(r["lat"], r["lon"]) for r in valid]
        samples = hepburn_sample_grid(calib, points)
        sample_by_latlon = {(s["lat"], s["lon"]): s for s in samples}

        now_iso = datetime.now().isoformat(timespec="seconds")
        log_recs = []
        for r in valid:
            s = sample_by_latlon.get((r["lat"], r["lon"]))
            if not s or s.get("intensity") is None:
                continue
            pre_dpp = (r.get("correction") or {}).get("corr_dpp", r["frames"][0]["dpp"])
            log_recs.append({
                "logged_at_local": now_iso, "lat": r["lat"], "lon": r["lon"],
                "hepburn_intensity_0_10": s["intensity"], "hepburn_label": s["label"],
                "hepburn_rgb": s["rgb"], "hepburn_valid_utc": s["hepburn_valid_utc"],
                "dpp": {"dpp": pre_dpp},
            })
        hepburn_append_batch(log_recs)

        fit = hepburn_fit_correction()
        hepburn_save_correction_fit(fit)
        n_req = fit.get("min_required", HEPBURN_MIN_SAMPLES_FOR_CORRECTION)
        if fit.get("applied"):
            print(f"[hepburn] 補正有効 n={fit['n']} r2={fit['r2']} kappa={fit['kappa']}")
        else:
            print(f"[hepburn] 補正未達(n={fit.get('n', 0)}/{n_req}) — 生DPPのまま継続")

        for r in valid:
            s = sample_by_latlon.get((r["lat"], r["lon"]))
            pre_dpp = r["frames"][0]["dpp"]
            hcorr = hepburn_blend_correction(pre_dpp, s.get("intensity") if s else None, fit=fit)
            r.setdefault("correction", {})["hepburn"] = hcorr
            if hcorr.get("applied"):
                r["frames"][0]["dpp"] = hcorr["final_dpp"]
    except Exception:
        import traceback
        print("[hepburn] 補正処理で例外発生、今回は無補正のまま継続します:")
        traceback.print_exc()
        return


def _send_discord(lines, webhook_url=None):
    """Discord Webhookへ通知を送る。requests未導入/URL未設定なら何もしない。"""
    webhook_url = webhook_url or get_settings().get("discord_webhook_url")
    if not webhook_url or requests is None or not lines:
        return
    content = "🌀 **Duct Alert**\n" + "\n".join(lines[:20])
    try:
        requests.post(webhook_url, json={"content": content[:1900]}, timeout=10)
    except Exception as e:
        print(f"[alert] Discord送信失敗: {e}")


def _detect_growing_peak(dpp_series, min_gain=2.0):
    """フレーム時系列(0,3,...24h)から「成長中でこの先ピークを迎える」パターンを検出する。"""
    vals = [v for v in dpp_series if v is not None]
    if len(vals) < 3:
        return None
    imax = max(range(len(vals)), key=lambda i: vals[i])
    if imax == 0 or vals[0] >= vals[imax] - min_gain:
        return None
    hours = FRAME_HOURS[imax] if imax < len(FRAME_HOURS) else imax * 3
    return {"hours": hours, "peak_dpp": vals[imax]}


def _check_and_send_alerts(results):
    """
    エリア×帯域のDPP/f_min条件、ピーク予告、Pass監視をチェックしDiscordへ通知する。
    HEPBURN補正適用後のresultsを使う(最終表示値と一致させるため)。
    同一条件は ALERT_COOLDOWN_MIN 分は再通知しない(プロセス内メモリでクールダウン管理)。
    """
    s = get_settings()
    if not s.get("alert_enabled") or not s.get("discord_webhook_url"):
        return
    threshold = float(s.get("alert_dpp_threshold", 70))
    radius_km = float(s.get("alert_radius_km", 300))
    now = time.time()
    events = []
    band_names = list(BANDS.keys())

    for r in results:
        if not r or not r.get("frames") or not r["frames"][0]:
            continue
        d = duct_haversine_km(ALERT_CENTER_LAT, ALERT_CENTER_LON, r["lat"], r["lon"])
        if d > radius_km:
            continue
        f0 = r["frames"][0]
        trap0 = f0.get("trap", "")
        for band in ALERT_BANDS:
            if band not in band_names:
                continue
            bidx = band_names.index(band)
            trapped0 = bidx < len(trap0) and trap0[bidx] == "1"
            if trapped0 and f0["dpp"] >= threshold:
                key = f"band:{round(r['lat'],1)}:{round(r['lon'],1)}:{band}"
                if now - _ALERT_STATE.get(key, 0) >= ALERT_COOLDOWN_MIN * 60:
                    _ALERT_STATE[key] = now
                    fmin = (f0.get("fmin") or {}).get(band)
                    fmin_txt = f" f_min≈{fmin}MHz" if fmin else ""
                    events.append(f"📡 **{band}** DPP={f0['dpp']:.0f} "
                                  f"({r['lat']:.1f},{r['lon']:.1f}){fmin_txt}")
        if ALERT_PEAK_LOOKAHEAD:
            dpp_series = [f["dpp"] if f else None for f in r["frames"]]
            peak = _detect_growing_peak(dpp_series)
            if peak:
                key = f"peak:{round(r['lat'],1)}:{round(r['lon'],1)}"
                if now - _ALERT_STATE.get(key, 0) >= ALERT_COOLDOWN_MIN * 60:
                    _ALERT_STATE[key] = now
                    events.append(f"📈 成長中 ({r['lat']:.1f},{r['lon']:.1f}) "
                                  f"+{peak['hours']}hでピーク予測(DPP {peak['peak_dpp']:.0f})")

    for wp in WATCHED_PASSES:
        info = analyze_pass(results, wp["a"][0], wp["a"][1], wp["b"][0], wp["b"][1])
        if not info.get("ok"):
            continue
        key = f"pass:{wp['name']}"
        broken = info["min_dpp"] < ALERT_PASS_BREAK_DPP
        strong = info["avg_dpp"] >= threshold
        if (broken or strong) and now - _ALERT_STATE.get(key, 0) >= ALERT_COOLDOWN_MIN * 60:
            _ALERT_STATE[key] = now
            tag = "⚠️ 途切れリスク" if broken else "✅ 良好"
            events.append(f"🛤️ Pass「{wp['name']}」{tag} 平均DPP={info['avg_dpp']:.0f} "
                          f"最小DPP={info['min_dpp']:.0f} 途切れ率={info['break_risk_pct']:.0f}%")

    if events:
        _send_discord(events)


_GEOCODE_CACHE = {}


def geocode_place_jp(query, timeout=10):
    """
    市町村名など地名文字列を緯度経度に変換する(Nominatim/OpenStreetMap使用)。
    日本国内を優先(countrycodes=jp)。個人の趣味利用を想定し高頻度照会は
    避ける(同一クエリはプロセス内メモリでキャッシュ)。
    """
    q = (query or "").strip()
    if not q:
        return {"ok": False, "error": "empty_query"}
    if q in _GEOCODE_CACHE:
        return _GEOCODE_CACHE[q]
    if requests is None:
        return {"ok": False, "error": "requests not installed"}
    try:
        resp = requests.get(
            "https://nominatim.openstreetmap.org/search",
            params={"q": q, "format": "json", "limit": 1, "countrycodes": "jp",
                    "accept-language": "ja"},
            headers={"User-Agent": "Mozilla/5.0 (personal hobby use; DPP forecast tool)"},
            timeout=timeout)
        resp.raise_for_status()
        arr = resp.json()
        if not arr:
            result = {"ok": False, "error": "not_found"}
        else:
            hit = arr[0]
            result = {"ok": True, "lat": round(float(hit["lat"]), 5),
                       "lon": round(float(hit["lon"]), 5),
                       "display_name": hit.get("display_name", q)}
    except Exception as e:
        result = {"ok": False, "error": str(e)}
    _GEOCODE_CACHE[q] = result
    return result


def analyze_pass(results, lat1, lon1, lat2, lon2, n=21):
    """2地点間の大円パス上をIDW補間でスコア化する(特定Pass監視)。"""
    pts = great_circle_points(lat1, lon1, lat2, lon2, n)
    vals = [idw_dpp_at(results, la, lo) for la, lo in pts]
    valid = [v for v in vals if v is not None]
    if not valid:
        return {"ok": False, "error": "no_data"}
    weak_ratio = sum(1 for v in valid if v < ALERT_PASS_BREAK_DPP) / len(valid)
    return {"ok": True, "distance_km": round(duct_haversine_km(lat1, lon1, lat2, lon2), 1),
            "avg_dpp": round(sum(valid) / len(valid), 1), "min_dpp": round(min(valid), 1),
            "max_dpp": round(max(valid), 1), "break_risk_pct": round(weak_ratio * 100, 1),
            "profile": [round(v, 1) if v is not None else None for v in vals],
            "points": [[round(la, 3), round(lo, 3)] for la, lo in pts]}


def antenna_height_curve(r, target_height_m=None, h_max=None, step=2.0):
    """
    primaryダクトのcoupling_efficiency()を高さ方向にスイープし、
    どこまでアンテナ高を上げれば結合効率が実用域(0.7以上)に届くかを返す。
    """
    if not r or not r.get("ducts0"):
        return {"ok": False, "error": "no_duct"}
    d0 = r["ducts0"][0]
    cur = target_height_m if target_height_m is not None else TX_HEIGHT_M
    h_max = h_max or (d0["top_m"] + 30.0)
    heights = list(np.arange(0.0, h_max + step, step))
    eff = [round(coupling_efficiency(h, d0), 3) for h in heights]
    cur_eff = round(coupling_efficiency(cur, d0), 3)
    target_h = next((h for h, e in zip(heights, eff) if e >= 0.7), None)
    return {"ok": True, "duct_base_m": round(d0["base_m"], 1), "duct_top_m": round(d0["top_m"], 1),
            "current_height_m": cur, "current_efficiency": cur_eff,
            "heights": [round(h, 1) for h in heights], "efficiency": eff,
            "suggested_height_m": (round(target_h, 1) if target_h is not None else None),
            "height_gap_m": (round(target_h - cur, 1)
                              if target_h is not None and target_h > cur else 0.0)}


def get_elevation_m(lat, lon, timeout=10):
    """Open-Meteo Elevation APIで標高(m)を取得する。"""
    if requests is None:
        return None
    try:
        resp = requests.get("https://api.open-meteo.com/v1/elevation",
                             params={"latitude": lat, "longitude": lon}, timeout=timeout)
        resp.raise_for_status()
        elevs = resp.json().get("elevation")
        return float(elevs[0]) if elevs else None
    except Exception as e:
        print(f"[elevation] 取得失敗: {e}")
        return None


def field_radar(results, lat, lon, n_az=24, max_range_km=600, step_km=25,
                 dpp_floor=15.0, elevation=True, max_extrapolation_km=250.0):
    """
    移動運用モード: 現在地から放射状に方位角ごとの最大到達距離を算出し、
    レーダーチャート用データを返す。elevation=Trueなら標高差も付加する。
    ・実データ格子点(nearest)からmax_extrapolation_kmより遠い地点はNone扱いとし、
      平坦な外挿による過大な到達距離を防ぐ。
    ・resultsの緯度経度の実データ範囲(+マージン)を出た地点も打ち切る
      (未計算領域/海外への飛び出し防止)。
    """
    valid = [r for r in results if r and r.get("frames") and r["frames"][0]]
    if valid:
        margin = 1.0
        lat_min = min(r["lat"] for r in valid) - margin
        lat_max = max(r["lat"] for r in valid) + margin
        lon_min = min(r["lon"] for r in valid) - margin
        lon_max = max(r["lon"] for r in valid) + margin
    else:
        lat_min = lat_max = lon_min = lon_max = None
    out = []
    for i in range(n_az):
        az = i * (360.0 / n_az)
        best_range = 0.0
        best_dpp = 0.0
        d = step_km
        while d <= max_range_km:
            plat, plon = destination_point(lat, lon, az, d)
            if lat_min is not None and not (lat_min <= plat <= lat_max and lon_min <= plon <= lon_max):
                break
            v = idw_dpp_at(results, plat, plon, max_dist_km=max_extrapolation_km)
            if v is None:
                break
            if v >= dpp_floor:
                best_range = d
                best_dpp = v
            elif v < dpp_floor * 0.5:
                break
            d += step_km
        out.append({"azimuth": round(az, 1), "max_range_km": round(best_range, 0),
                     "end_dpp": round(best_dpp, 1)})
    result = {"ok": True, "lat": lat, "lon": lon, "radar": out}
    if elevation:
        home_elev = get_elevation_m(lat, lon)
        duct_r = _nearest_result(results, lat, lon)
        duct_base = None
        if duct_r and duct_r.get("ducts0"):
            duct_base = round(duct_r["ducts0"][0]["base_m"], 1)
        result["elevation_m"] = round(home_elev, 1) if home_elev is not None else None
        result["duct_base_m"] = duct_base
        if home_elev is not None and duct_base is not None:
            result["height_to_duct_m"] = round(duct_base - home_elev, 1)
    return result


def _spots_daemon(cal_path, region, stop_event, spots_fetch=None):
    print("[spots] 全局スポット観測 開始 (周期はUI設定を参照)")
    c = 0; warm_prev = None
    while not stop_event.is_set():
        s = get_settings()
        eff_iv = max(10, s["interval_sec"])
        if not s["paused"]:
            try:
                with _CAL_IO_LOCK:
                    cal = load_calibration(cal_path)
                    eff_sb, eff_iv, warm = _effective_fetch_params(cal, s)
                    eff_iv = max(10, eff_iv)
                    if warm != warm_prev:
                        print(f"[spots] ウォームアップモード: "
                              f"{'ON(取得を加速し集中学習中)' if warm else 'OFF(通常周期に復帰)'} "
                              f"κ={mean_kappa(cal):.3f}(閾値{WARMUP_KAPPA_THRESHOLD})")
                        warm_prev = warm
                    st = spots_cycle(cal, region, eff_sb, spots_fetch)
                    save_calibration(cal, cal_path)
                backoff = psk_backoff_extra_sec()
                if backoff:
                    eff_iv += backoff
                lq = learning_quality(cal)
                eta = lq["eta_hours"]
                eta_s = ("到達" if (eta is not None and eta <= 0)
                         else (f"{eta}h" if eta is not None else "—"))
                tag = " 🔥warmup" if warm else ""
                if backoff:
                    tag += f" ⏳backoff+{backoff}s(PSK{_PSK_BACKOFF['fail_streak']}連続失敗)"
                print(f"[spots] c{c}: +{st['pos']}/-{st['neg']} 異常{st['spots']} "
                      f"cone={st['cone']}° 成熟{lq['maturity_pct']}% "
                      f"κ={lq['mean_confidence']} 実用まで{eta_s}{tag}")
            except Exception as e:
                print(f"[spots] エラー: {e}")
            c += 1
        if stop_event.wait(eff_iv):
            break
    print("[spots] 観測終了")


def forecast_cycle(region, cal_path, weather_override=None, open_browser=False,
                   refresh_sec=300, weather_cache_sec=WEATHER_CACHE_SEC):
    global _duct_cycle_diag
    _cycle_t0 = time.time()
    cal = load_calibration(cal_path)
    if weather_override is not None:
        points, weather = weather_override; wx_cached = True
    else:
        if requests is None:
            print("requests 未導入。pip install requests"); return None
        points, weather, wx_cached = get_weather_cached(region, weather_cache_sec)
    weather_ok = sum(1 for w in weather if w)
    tag = "気象=キャッシュ流用" if wx_cached else "気象=新規取得"
    print(f"[予測] 多層/WKB/較正/フレーム計算中... ({tag}, 較正=最新, "
          f"気象取得成功={weather_ok}/{len(points)}地点)")
    if (not wx_cached) and weather_ok == 0 and len(points) > 0:
        print("[予測] 警告: 気象データが全地点で取得できませんでした"
              "(通信環境/Open-Meteo側の制限の可能性)。この場合、解析は"
              "0/全地点で失敗し、地図は前回サイクルのまま更新されません。")
    results = []
    _first_err_shown = False
    for (lat, lon), wj in zip(points, weather):
        try:
            results.append(analyze_point_frames(lat, lon, wj, cal))
        except Exception as e:
            results.append(None)
            if not _first_err_shown:
                _first_err_shown = True
                import traceback
                print(f"[予測] 地点解析で例外発生 (lat={lat}, lon={lon}): "
                      f"{type(e).__name__}: {e}")
                traceback.print_exc()
    valid = sum(1 for r in results if r)
    _publish_dpp_map(results, region["step_deg"])
    _hepburn_apply_correction(results)
    with _RESULTS_LOCK:
        _LATEST_RESULTS.clear(); _LATEST_RESULTS.extend(results)
    _check_and_send_alerts(results)
    out = generate_html(results, region, cal, refresh_sec, OUTPUT_HTML)
    print(f"[予測] 有効{valid}/{len(points)} → {os.path.basename(out) if out else '生成失敗'}"
          f" (所要{time.time()-_cycle_t0:.1f}秒)")
    if not out:
        print("[予測] 警告: 今回のサイクルではHTML地図を書き出せませんでした"
              "(有効な解析地点が0件、または内部エラー)。次回サイクルで再試行します。")
    _duct_cycle_diag.update({
        'cycles_run': _duct_cycle_diag['cycles_run'] + 1,
        'last_started_at': datetime.fromtimestamp(_cycle_t0).isoformat(timespec='seconds'),
        'last_finished_at': datetime.now().isoformat(timespec='seconds'),
        'last_duration_sec': round(time.time() - _cycle_t0, 1),
        'last_points_total': len(points),
        'last_weather_ok': weather_ok,
        'last_analysis_valid': valid,
        'last_wrote_map': bool(out),
        'last_step_deg': region.get('step_deg'),
        'last_is_coarse_pass': region.get('step_deg', 0) > DEFAULT_REGION.get('step_deg', 2.0),
    })
    if out and open_browser:
        webbrowser.open("file://" + os.path.abspath(out))
    return results


def run_service(region=None, weather_override=None, spots_fetch=None,
                stop_event=None, base_url=""):
    region = region or DEFAULT_REGION
    if not os.path.exists(CALIBRATION_FILE):
        save_calibration(default_calibration(region["step_deg"], region),
                         CALIBRATION_FILE)
    stop_event = stop_event or threading.Event()
    th = threading.Thread(target=_spots_daemon,
                          args=(CALIBRATION_FILE, region, stop_event, spots_fetch),
                          name="dpp-spots", daemon=True)
    th.start()
    s0 = get_settings()
    if s0.get("wspr_backfill_enabled", True):
        threading.Thread(target=wspr_backfill_worker,
                         args=(CALIBRATION_FILE, region, stop_event),
                         name="dpp-wspr-backfill", daemon=True).start()
    if spots_fetch is None:
        threading.Thread(target=pskreporter_startup_backfill,
                         args=(CALIBRATION_FILE, region),
                         name="dpp-psk-backfill", daemon=True).start()
    start = time.time(); s = get_settings()
    print(f"[service] 常駐開始: {'無制限' if s['total_minutes']==0 else str(s['total_minutes'])+'分'}"
          f" / {s['interval_sec']}秒周期 (設定はUIから変更可)")
    try:
        while not stop_event.is_set():
            _RECOMPUTE_NOW.clear()
            s = get_settings(); total_min = s["total_minutes"]
            if total_min and (time.time() - start) >= total_min * 60:
                print("=" * 60)
                print(f"[service] ⏹ 常駐時間({total_min}分)に到達したため終了します"
                      f"(サーバーも停止し、ブラウザは接続エラーになります)")
                print("           ずっと動かし続けたい場合は⚙設定の「常駐時間」を"
                      "0(無制限)にしてください")
                print("=" * 60)
                break
            if not s["paused"]:
                try:
                    forecast_cycle(region, CALIBRATION_FILE,
                                   weather_override=weather_override,
                                   open_browser=False, refresh_sec=s["interval_sec"],
                                   weather_cache_sec=s["weather_cache_sec"])
                except Exception:
                    import traceback
                    print("[service] 今回のサイクルで例外発生、サーバは継続します:")
                    traceback.print_exc()
            interval = max(10, get_settings()["interval_sec"])
            if _RECOMPUTE_NOW.wait(timeout=interval):
                print("[service] 即時再計算要求を受信")
    finally:
        stop_event.set(); th.join(timeout=8)
        print("[service] 終了")







DUCT_LEVEL_THRESHOLDS = [
    (85, "hot",  "絶好調"),
    (70, "good", "良好"),
    (50, "warn", "普通"),
    (30, "weak", "弱い"),
    (0,  "low",  "低い"),
]


def duct_dpp_level(dpp):
    """DPP値(0-100)をEDFS本体と共通の5段階(hot/good/warn/weak/low)へ変換する。
    しきい値はDUCT側の地図等値線(CONTOUR_LEVELS=[30,50,70,85])と同一。"""
    if dpp is None:
        return {"cls": "low", "label": "データ不足", "dpp": None}
    for th, cls, label in DUCT_LEVEL_THRESHOLDS:
        if dpp >= th:
            return {"cls": cls, "label": label, "dpp": round(float(dpp), 1)}
    return {"cls": "low", "label": "低い", "dpp": round(float(dpp), 1)}


def duct_apply_home_location(lat, lon, label=None):
    """EDFS本体が解決した現在地(ui_resolve_location等)をDUCT側の基準点へ反映する。
    起動時に一度、以後は位置が大きく変わった場合のみ呼び出す想定。
    Hepburn比較基準点・デフォルトアラート中心のみを更新し、
    DEFAULT_REGION(観測グリッド範囲=日本全域)は変更しない
    (現在地が範囲内であれば十分であり、グリッドを動かすと較正学習データの
    セルキーが変わってしまい、蓄積済みの学習データが無駄になるため)。"""
    global HOME_LAT, HOME_LON, ALERT_CENTER_LAT, ALERT_CENTER_LON
    try:
        if lat is not None and lon is not None and \
           -90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0:
            HOME_LAT, HOME_LON = float(lat), float(lon)
            ALERT_CENTER_LAT, ALERT_CENTER_LON = HOME_LAT, HOME_LON
            print(f"[DUCT] 基準点をEDFS現在地へ反映: {HOME_LAT:.3f},{HOME_LON:.3f}"
                  + (f" ({label})" if label else ""))
    except Exception as e:
        print(f"[DUCT] 基準点反映に失敗(既定値のまま続行): {e}")


def duct_point_query(lat, lon, max_dist_km=250.0):
    """任意座標(lat, lon)についてダクト伝搬が期待できるかを返す。
    EDFS本体の「現在地から面が開通していく」という考え方とは異なり、
    入力された地点そのものの見込みだけを answer する点に注意。
    直近サイクルの計算済みグリッド(_LATEST_RESULTS)が無い場合はエラーを返す。"""
    with _RESULTS_LOCK:
        results = list(_LATEST_RESULTS)
    if not results:
        return {"ok": False, "error": "まだ1回もDUCT予測サイクルが完了していません"
                                       "(起動直後は数分お待ちください)"}
    nearest = _nearest_result(results, lat, lon)
    if nearest is None:
        return {"ok": False, "error": "有効な計算結果がありません"}
    dist_km = duct_haversine_km(lat, lon, nearest["lat"], nearest["lon"])
    idw_now = idw_dpp_at(results, lat, lon, max_dist_km=max_dist_km)
    reliable = (idw_now is not None) and (dist_km <= max_dist_km)
    frames_out = []
    for f in (nearest.get("frames") or []):
        if f is None:
            frames_out.append(None)
            continue
        frames_out.append({
            "t": f.get("t"),
            "level": duct_dpp_level(f.get("dpp")),
            "bands": {
                b: {"trapped": (f.get("trap", "")[i] == "1")
                    if i < len(f.get("trap", "")) else None,
                    "f_min_mhz": (f.get("fmin") or {}).get(b)}
                for i, b in enumerate(BANDS)
            },
        })
    cal = load_calibration(CALIBRATION_FILE)
    lq = learning_quality(cal)
    return {
        "ok": True,
        "query": {"lat": round(float(lat), 4), "lon": round(float(lon), 4)},
        "nearest_grid": {"lat": nearest["lat"], "lon": nearest["lon"],
                         "distance_km": round(dist_km, 1)},
        "reliable": reliable,
        "reliability_note": (
            "グリッド解像度(約2度)内の補間値です" if reliable else
            f"最近傍の計算点まで{round(dist_km)}km離れており、参考値として扱ってください"),
        "now": duct_dpp_level(idw_now if idw_now is not None
                              else (nearest.get("frames") or [{}])[0].get("dpp")),
        "confidence_kappa": nearest.get("confidence"),
        "components": nearest.get("components"),
        "summary": nearest.get("summary"),
        "timeline": frames_out,
        "learning": lq,
        "duct_source_version": DUCT_SOURCE_VERSION,
    }


def duct_status_summary():
    """EDFSダッシュボードのカード表示用。自局(EDFS現在地)の現況+学習成熟度。"""
    cal = load_calibration(CALIBRATION_FILE)
    lq = learning_quality(cal)
    with _RESULTS_LOCK:
        results = list(_LATEST_RESULTS)
    home = None
    if results:
        r = _nearest_result(results, HOME_LAT, HOME_LON)
        if r is not None:
            f0 = (r.get("frames") or [None])[0]
            home = {"level": duct_dpp_level(f0.get("dpp") if f0 else None),
                    "summary": r.get("summary"),
                    "confidence_kappa": r.get("confidence")}
    bootstrap = dict(globals().get('_duct_bootstrap_diag', {}) or {})
    return {
        "ok": True,
        "home_lat": HOME_LAT, "home_lon": HOME_LON,
        "home": home,
        "learning": lq,
        "settings": get_settings(),
        "duct_source_version": DUCT_SOURCE_VERSION,
        "cycles_computed": len(results),
        "bootstrap_ok": bootstrap.get('ok'),
        "bootstrap_error": bootstrap.get('error'),
        "bootstrap_stage": bootstrap.get('stage'),
        "bootstrap_ts": bootstrap.get('ts'),
        "first_cycle_started": bootstrap.get('first_cycle_started'),
        "first_cycle_done": bootstrap.get('first_cycle_done'),
        "cycle_diag": dict(globals().get('_duct_cycle_diag', {}) or {}),
    }


_duct_threads_started = threading.Event()


_duct_bootstrap_diag = {'ok': None, 'error': None, 'ts': None,
                        'stage': None, 'first_cycle_started': False,
                        'first_cycle_done': False}

_duct_cycle_diag = {
    'cycles_run': 0, 'last_started_at': None, 'last_finished_at': None,
    'last_duration_sec': None, 'last_points_total': 0,
    'last_weather_ok': 0, 'last_analysis_valid': 0, 'last_wrote_map': False,
    'last_step_deg': None, 'last_is_coarse_pass': None, 'last_hard_timeout': False,
}
_duct_stale_map_last_logged_mtime = None

_duct_stale_run = {'future': None, 'executor': None, 'started_at': None}


def duct_start_background(stop_event, region=None):
    """DUCTのデータ収集+学習を、EDFS本体サイクルと完全に並行して開始する。
    EDFS本体のstop_event(_edfs_stop_event)を共有し、EDFS停止時に道連れで
    停止する。EDFS本体のHTTPサーバ・ダッシュボードとはポートを分けず、
    EDFS側のdo_GET/do_POSTへ追加したハンドラ(duct_point_query等)経由で
    同一プロセス内から直接関数呼び出しする(独自HTTPサーバ・独自ポートは
    使用しない)。

    Rev16.9(バグ修正・最重要): 従来はこの関数のほぼ全体が一切try/exceptで
    保護されておらず、起動直後のsave_calibration()等が例外を出すと
    duct-bootstrapスレッド全体が無言でクラッシュし、実際に地図を生成する
    _forecast_loopスレッドが二度と起動されない(かつ_duct_threads_started
    フラグにより再試行もされない)という不具合があった。この場合、コンソールに
    一切のエラーも出ず、/api/duct/mapは「生成中です」のまま永久に進まなく
    見える(「1時間以上続く」報告の直接原因)。
    対策: ①最重要スレッド(_forecast_loop)を最優先・最小限の前提で起動する
    (較正ファイルはload_calibration()側で自動フォールバックするため必須では
    ない)。②各起動ステップを個別にtry/exceptで保護する。③診断情報
    (_duct_bootstrap_diag)を随時更新し、/api/duct/status・/api/duct/map・
    ダッシュボードのDUCTカードが「起動失敗(自動回復しない)」と「起動直後で
    まだ計算中」を明確に区別して利用者に伝えられるようにする。
    """
    global _duct_bootstrap_diag
    if _duct_threads_started.is_set():
        return
    _duct_threads_started.set()
    _duct_bootstrap_diag['stage'] = 'start'
    region = region or DEFAULT_REGION

    try:
        load_settings()
    except Exception as e:
        print(f"[DUCT] 警告: 設定読込に失敗しました(既定値で続行): {type(e).__name__}: {e}")

    def _forecast_loop():
        global _duct_bootstrap_diag
        print("[DUCT] 予測サイクル 開始(EDFS本体と並行動作)")
        _duct_bootstrap_diag['first_cycle_started'] = True
        first_pass = True
        while not stop_event.is_set():
            s = get_settings()
            if not s["paused"]:
                cycle_region = region
                if first_pass:
                    cycle_region = dict(region)
                    cycle_region["step_deg"] = max(region.get("step_deg", 2.0) * 1.2, 2.5)
                    print(f"[DUCT] 初回サイクルは粗いグリッド(間隔{cycle_region['step_deg']}度)"
                          f"で高速に実行します。2回目以降は通常解像度({region['step_deg']}度)に戻ります。"
                          f"(F2層機能廃止により従来より高解像度化)")
                _stale_fut = _duct_stale_run.get('future')
                if _stale_fut is not None and not _stale_fut.done():
                    _elapsed = time.time() - (_duct_stale_run.get('started_at') or time.time())
                    print(f"[DUCT] ⚠ 前回タイムアウトしたサイクルが依然として実行中"
                          f"(経過約{_elapsed:.0f}秒)のため、今回のサイクル開始を"
                          f"見送ります(多重実行によるファイル競合を防止)。")
                else:
                    if _stale_fut is not None:
                        try:
                            _duct_stale_run['executor'].shutdown(wait=False)
                        except Exception:
                            pass
                        _duct_stale_run['future'] = None
                        _duct_stale_run['executor'] = None
                    try:
                        HARD_CYCLE_TIMEOUT_SEC = 180
                        _duct_cycle_diag['last_hard_timeout'] = False
                        _ex = _v13_futures.ThreadPoolExecutor(max_workers=1)
                        _cycle_started_at = time.time()
                        fut = _ex.submit(forecast_cycle, cycle_region, CALIBRATION_FILE,
                                         None, False, s["interval_sec"], s["weather_cache_sec"])
                        try:
                            fut.result(timeout=HARD_CYCLE_TIMEOUT_SEC)
                            _duct_bootstrap_diag['first_cycle_done'] = True
                            first_pass = False
                            _ex.shutdown(wait=False)
                        except _v13_futures.TimeoutError:
                            print(f"[DUCT] ⚠ 予測サイクルが{HARD_CYCLE_TIMEOUT_SEC}秒を超えても"
                                  f"完了しませんでした(気象取得または画像取得がネットワーク側で"
                                  f"ハングしている可能性)。このサイクルは諦めて次回に進みますが、"
                                  f"実行中のスレッド自体は停止できないため、完了するまで新規サイクルの"
                                  f"開始を見送ります。")
                            _duct_cycle_diag['last_hard_timeout'] = True
                            _duct_stale_run['future'] = fut
                            _duct_stale_run['executor'] = _ex
                            _duct_stale_run['started_at'] = _cycle_started_at
                    except Exception:
                        import traceback
                        print("[DUCT] 予測サイクルで例外発生(継続します):")
                        traceback.print_exc()
            interval = max(10, get_settings()["interval_sec"])
            if _RECOMPUTE_NOW.wait(timeout=interval):
                _RECOMPUTE_NOW.clear()
        print("[DUCT] 予測サイクル 終了")

    try:
        threading.Thread(target=_forecast_loop, name="duct-forecast", daemon=True).start()
        _duct_bootstrap_diag.update({'ok': True, 'error': None, 'stage': 'forecast_thread_started',
                                     'ts': datetime.now().isoformat(timespec='seconds')})
    except Exception as e:
        import traceback
        print(f"[DUCT] 致命的エラー: 予測サイクルスレッドの起動に失敗しました: {type(e).__name__}: {e}")
        traceback.print_exc()
        _duct_bootstrap_diag.update({'ok': False, 'error': f"{type(e).__name__}: {e}",
                                     'stage': 'forecast_thread_failed',
                                     'ts': datetime.now().isoformat(timespec='seconds')})
        return

    try:
        if not os.path.exists(CALIBRATION_FILE):
            save_calibration(default_calibration(region["step_deg"], region),
                             CALIBRATION_FILE)
    except Exception as e:
        print(f"[DUCT] 警告: 較正ファイルの事前作成に失敗しました"
              f"(予測サイクル自体はメモリ上の既定値で継続します): {type(e).__name__}: {e}")

    try:
        loc = ui_resolve_location()
        if loc is not None:
            duct_apply_home_location(loc[0], loc[1], label=loc[2] if len(loc) > 2 else None)
    except Exception as e:
        print(f"[DUCT] 起動時の現在地取得に失敗(既定値で続行): {e}")

    try:
        th_spots = threading.Thread(target=_spots_daemon,
                                    args=(CALIBRATION_FILE, region, stop_event, None),
                                    name="duct-spots", daemon=True)
        th_spots.start()
    except Exception as e:
        print(f"[DUCT] 警告: スポット収集スレッドの起動に失敗しました: {type(e).__name__}: {e}")

    try:
        s0 = get_settings()
        if s0.get("wspr_backfill_enabled", True):
            threading.Thread(target=wspr_backfill_worker,
                             args=(CALIBRATION_FILE, region, stop_event),
                             name="duct-wspr-backfill", daemon=True).start()
    except Exception as e:
        print(f"[DUCT] 警告: WSPRバックフィルスレッドの起動に失敗しました: {type(e).__name__}: {e}")

    try:
        threading.Thread(target=pskreporter_startup_backfill,
                         args=(CALIBRATION_FILE, region),
                         name="duct-psk-backfill", daemon=True).start()
    except Exception as e:
        print(f"[DUCT] 警告: PSKReporter起動時バックフィルスレッドの起動に失敗しました: {type(e).__name__}: {e}")


if __name__ == '__main__':
    warnings.filterwarnings('ignore', category=FutureWarning)
    edfs_v13_run_main()




