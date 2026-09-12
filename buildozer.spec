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
source.include_exts = py,png,jpg,kv,atlas,json,txt

# ビルド作業に不要なフォルダ・ファイルをAPKから除外する。
# (.github/ワークフローや.git履歴を誤って同梱しないための安全策)
source.exclude_dirs = .github,.git,.buildozer,bin,__pycache__
source.exclude_patterns = README.md,*.md,.gitignore

version = 1.0

# ------------------------------------------------------------------------
# Python依存パッケージ
# ------------------------------------------------------------------------
# 【重要・ビルド安定性のための方針】
#   ・scikit-learn は含めない。ESDUCT_discription.py 側は
#     "try: from sklearn... except ImportError: _SKLEARN_AVAILABLE = False"
#     という設計になっており、sklearn が無い環境では自動的に単純相関
#     ベースの学習にフォールバックする(致命的な機能欠落にはならない)。
#     一方 scikit-learn は python-for-android に安定した recipe が無く、
#     ビルド失敗の主要因になりやすいため、意図的に除外している。
#   ・lxml も含めない。BeautifulSoup の呼び出しは全て
#     BeautifulSoup(html, "html.parser") であり、lxml パーサーには
#     依存していない(lxml は Android 向けビルドが特に不安定なため、
#     不要な依存を増やさない)。
#   ・numpy / pandas はバージョンを明示的に固定していない。
#     python-for-android のrecipeは対応バージョンが決まっており、
#     ここでバージョンを固定すると、ビルド時に取得したrecipeのバージョンと
#     食い違ってパッチ適用に失敗するリスクがあるため、あえて無指定にして
#     recipeの既定バージョンに委ねる。
#   ・kivy はバージョンを 2.3.1 に固定する。python-for-android が
#     2026年時点で「Kivy 3系(新しいbootstrap contract)」と
#     「Kivy 2.3.1(旧bootstrap contract)」の両方をサポートしている
#     ことが確認できているため、実績のある 2.3.1 を明示指定し、
#     将来的なKivyの破壊的変更の影響を受けないようにする。
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
#                        (Android 10+では制限されるが、失敗時はアプリ専用
#                        領域へ自動フォールバックする設計のため必須ではない)
android.permissions = INTERNET,ACCESS_FINE_LOCATION,ACCESS_COARSE_LOCATION,WAKE_LOCK,READ_EXTERNAL_STORAGE,WRITE_EXTERNAL_STORAGE

# 対象アーキテクチャはarm64-v8aのみに限定する。
# 【理由】armeabi-v7aも同時ビルドすると、numpy/pandasのような
# ネイティブコードを含むrecipeのビルド時間・失敗確率が単純に倍増する。
# 現行の主要Android端末(Pixel 8a等)は全てarm64-v8a対応のため、
# 実用上の支障はない。
android.archs = arm64-v8a

# minapiは低くしすぎるとp4aの対応が薄い古いAPIレベル向けの分岐に
# 入りやすくなるため、Android 7.0(API24)以降に限定する。
android.minapi = 24

# android.api / android.ndk / android.sdk はあえて指定しない(コメントアウト)。
# 【理由】ここでバージョンを固定すると、buildozer/python-for-androidが
# 将来バージョンアップした際に「指定した値が既にサポート対象外」と
# なってビルドが失敗するリスクがある。無指定であれば、buildozerが
# インストールされているバージョンに応じた「現時点で推奨される
# SDK/NDKバージョン」を自動的に選択するため、経年劣化によるビルド
# 破綻を避けやすい。
# android.api = 33
# android.ndk = 25b

android.accept_sdk_license = True
android.enable_androidx = True
android.allow_backup = True

p4a.bootstrap = sdl2

# ログを詳細に出す(CI上でビルド失敗時の原因調査をしやすくするため)。
log_level = 2
warn_on_root = 1

[buildozer]
log_level = 2
warn_on_root = 1
