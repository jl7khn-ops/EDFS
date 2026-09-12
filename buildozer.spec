[app]

# ------------------------------------------------------------------------
# 基本情報
# ------------------------------------------------------------------------
title = ESDUCT
package.name = esduct
# reverse-DNS形式。小文字・数字・ドットのみ使用可能。
package.domain = jp.jl7khn

source.dir = .
# ここに列挙した拡張子のファイルのみAPKに同梱される。
# HTMLダッシュボードを使う前提なら html,css,js も含めておく。
source.include_exts = py,png,jpg,kv,atlas,json,txt,html,css,js

# ビルド作業に不要なフォルダ・ファイルをAPKから除外する。
# (.github/ワークフローや.git履歴を誤って同梱しないための安全策)
source.exclude_dirs = .github,.git,.buildozer,bin,__pycache__
source.exclude_patterns = README.md,*.md,.gitignore

version = 1.0

# ------------------------------------------------------------------------
# Python依存パッケージ
# ------------------------------------------------------------------------
# ・scikit-learn は含めない（p4a recipe不安定のため）。
#   ESDUCT_discription.py 側で ImportError 時に単純相関へフォールバック済み。
# ・lxml も含めない（BeautifulSoup は "html.parser" のみ使用）。
# ・numpy / pandas はバージョン固定しない（p4a recipe側に委ねる）。
# ・kivy は 2.3.1 に固定（旧bootstrap contractでの安定動作を優先）。
requirements = python3,kivy==2.3.1,requests,beautifulsoup4,pillow,numpy,pandas,plyer,pyjnius,certifi,charset-normalizer,idna,urllib3,six,python-dateutil,pytz

# GUIはHTMLダッシュボード(ブラウザ)側に一任しているため、
# 画面表示は縦向き固定で十分。
orientation = portrait
fullscreen = 0

# ------------------------------------------------------------------------
# Android固有設定
# ------------------------------------------------------------------------
# INTERNET            : NICT/PSKReporter/NOAA等への通信、HTMLダッシュボード配信
# ACCESS_FINE/COARSE_LOCATION : GPSによる現在地取得(plyer経由)
# WAKE_LOCK           : 定期観測サイクル実行中のスリープ防止
# READ/WRITE_EXTERNAL_STORAGE : 学習データ・Evidence JSON等の共有ストレージ保存
android.permissions = INTERNET,ACCESS_FINE_LOCATION,ACCESS_COARSE_LOCATION,WAKE_LOCK,READ_EXTERNAL_STORAGE,WRITE_EXTERNAL_STORAGE

# 対象アーキテクチャはarm64-v8aのみに限定する。
android.archs = arm64-v8a

# minapiはAndroid 7.0(API24)以降に限定。
android.minapi = 24

# GitHub Actions（kivy/buildozer:stable）のデフォルトに合わせて明示指定。
android.api = 33
android.ndk = 25b

android.accept_sdk_license = True
android.enable_androidx = True
android.allow_backup = True

# Kivyベースの通常アプリなので sdl2 bootstrap を使用。
p4a.bootstrap = sdl2

# p4aのブランチを develop に固定し、レシピ更新との整合性を取りやすくする。
p4a.branch = develop

# CI上でビルド失敗時の原因調査をしやすくするためのログ設定。
log_level = 2
warn_on_root = 1

[buildozer]
log_level = 2
warn_on_root = 1
