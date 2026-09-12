# ESDUCT APK ビルドプロジェクト

`ESDUCT_discription.py`(Es & Duct Propagation Forecast Program)を
Android APK化するための一式です。ビルドは GitHub Actions 上で行います。

## フォルダ構成

```
esduct_apk/
├── .github/
│   └── workflows/
│       └── build-apk.yml      … GitHub ActionsのビルドWorkflow
├── main.py                    … APKのエントリポイント(Kivy製の最小シェル)
├── ESDUCT_discription.py       … 本体(コメント削除版、ロジックは無改変)
├── buildozer.spec              … Android APK化の設定ファイル
├── .gitignore
└── README.md                   … このファイル
```

## 使い方(GitHub Actionsでビルド)

1. このフォルダの中身をそのまま新規(または既存)のGitHubリポジトリの
   ルート直下へ配置し、`main`ブランチへpush する。
2. GitHubリポジトリの「Actions」タブを開くと、`Build ESDUCT APK` という
   Workflowが自動的に実行される(pushで自動実行、または
   「Run workflow」ボタンから手動実行も可能)。
3. ビルドが成功すると、その実行結果(Run)のページ下部
   「Artifacts」欄に `ESDUCT-debug-apk` という名前でAPKファイルが
   アップロードされる。ダウンロードして端末にインストールする。

**初回ビルドは40〜60分程度かかります**(numpy/pandasのコンパイルに
時間がかかるため)。2回目以降は `.buildozer` ディレクトリのキャッシュが
効くため、大幅に短縮されます。

## なぜこの構成にしたか(ビルド安定性のための設計判断)

今回は「GitHub Actionsでのビルドはエラーになりやすい」という前提を
踏まえ、以下の点を意識して構成しています。

### 1. Docker系のビルドAction(`buildozer-action`等)は使わない

`ArtemSBulgakov/buildozer-action` のようなDocker経由でbuildozerを実行する
Actionは手軽ですが、内部で参照している `kivy/buildozer:latest` という
Dockerイメージが破損しており `docker pull` に失敗するという報告が
複数確認できました。この種の外部要因による失敗は自分では直しようが
ないため、依存を避け、Ubuntu ランナー上で buildozer を直接
`pip install` して実行する構成にしています。ビルド過程が完全に
ログへ出るため、失敗した場合の原因調査もしやすくなります。

### 2. `runs-on: ubuntu-22.04` を明示指定(`ubuntu-latest`は使わない)

`ubuntu-latest` はGitHub側の都合で参照先イメージが将来切り替わり得ます。
切り替わった瞬間に、今回インストールしている `libncurses5-dev` のような
特定バージョンのaptパッケージが存在しなくなり、ある日突然ビルドが
壊れる、という事故が典型的に起こります。バージョンを固定した
`ubuntu-22.04` を明示することで、この種の不意打ちを防いでいます。

### 3. Cythonを`0.29.36`に固定する

numpy/pandasのrecipeビルド中に、Cython 3系だと失敗する事例が複数
報告されています。実績のある0.29系バージョンに明示的に固定することで、
この失敗パターンを避けています。

### 4. `scikit-learn` はrequirementsに含めない

`ESDUCT_discription.py` 側は

```python
try:
    from sklearn.linear_model import Ridge
    ...
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False
```

という設計になっており、sklearnが無い環境では自動的に単純相関ベースの
学習にフォールバックします(機能が完全に失われるわけではありません)。
一方 scikit-learn は python-for-android 向けの安定したrecipeが無く、
含めるとビルド失敗の主要因になりやすいため、意図的に除外しています。

### 5. `lxml` も含めない

ソース内の `BeautifulSoup(html, "html.parser")` は全て標準ライブラリの
`html.parser` を使っており、lxmlパーサーには依存していません。
lxmlはAndroid向けビルドが特に不安定なため、そもそも不要な依存として
除外しています。

### 6. numpy/pandas/kivy 以外のバージョンは固定していない

python-for-androidのrecipeは対応バージョンが決まっており、ここで
バージョンを固定すると、ビルド時に取得したrecipeのバージョンと
食い違ってパッチ適用に失敗するリスクがあります。numpy/pandasは
あえて無指定にし、recipeの既定バージョンに委ねています。
Kivyのみ、python-for-androidが現時点でサポートを明言している
`2.3.1`(旧bootstrap contract)に固定し、将来のKivy 3系の破壊的変更の
影響を受けないようにしています。

### 7. `android.api` / `android.ndk` / `android.sdk` は指定しない

ここでバージョンを固定すると、buildozer/python-for-androidが将来
バージョンアップした際に「指定した値が既にサポート対象外」となって
ビルドが失敗するリスクがあります。無指定であれば、buildozerが
インストールされているバージョンに応じた「現時点で推奨される
SDK/NDKバージョン」を自動選択するため、経年劣化によるビルド破綻を
避けやすくなります。

### 8. `android.archs = arm64-v8a` のみに限定

armeabi-v7aも同時ビルドすると、numpy/pandasのようなネイティブコードを
含むrecipeのビルド時間・失敗確率が単純に倍増します。現行の主要
Android端末は全てarm64-v8a対応のため、実用上の支障はありません。

## main.pyについて

`ESDUCT_discription.py` 本体はKivyに依存しない「HTTPサーバー常駐型」の
アプリケーションで、GUIは全てブラウザ側のHTMLダッシュボード
(本体内蔵のHTTPサーバー、既定 `http://127.0.0.1:8765/`)で完結します。

`main.py` はAndroid上でプロセスを常駐させるための最小限のKivyアプリ
(Activity)を提供しているだけであり、`ESDUCT_discription.py` 側の
ロジック・学習データ・HTTPサーバー実装には一切手を加えていません。
アプリを起動すると:

1. 位置情報等のランタイム権限をリクエスト
2. バックグラウンドスレッドで `edfs_v13_run_main()` を起動
3. 自動的にブラウザが開き、ダッシュボードが表示される
   (自動起動しない場合は `http://127.0.0.1:8765/` を手動で開く)
4. 画面下の「常駐停止」ボタンで安全に停止できる

## ビルドが失敗した場合の確認ポイント

1. **GitHub Actionsの「Artifacts」に `buildozer-log` が出ていないか**
   ビルド失敗時は自動的に詳細ログがアーティファクトとして
   アップロードされるよう設定しています。まずこれを確認してください。

2. **numpy/pandasのビルドで失敗している場合**
   `.github/workflows/build-apk.yml` 内の Cython バージョン
   (`0.29.36`)を、その時点での python-for-android 側の推奨値に
   合わせて調整してください。

3. **apt パッケージが見つからない、というエラーが出た場合**
   `ubuntu-22.04` の既定リポジトリ構成が変わった可能性があります。
   `sudo apt-get install` の対象パッケージ名を最新の状況に合わせて
   見直してください。

4. **キャッシュが古くて逆に失敗する場合**
   GitHubリポジトリの「Actions」→「Caches」からキャッシュを手動で
   削除し、再度ビルドを実行してください(`buildozer.spec` を変更
   すればキャッシュキーが変わり自動的に新規ビルドになりますが、
   変更なしで強制的にクリーンビルドしたい場合はキャッシュの手動
   削除が必要です)。

## 注意事項

- 生成されるAPKは**署名なしのdebug版**です。配布・Google Playへの
  登録には別途キーストアによる署名が必要です(本一式のスコープ外)。
- Android 10以降はスコープドストレージの制約により、共有ストレージ
  (`/storage/emulated/0/...`)への書き込みが制限される場合があります。
  `ESDUCT_discription.py` 側はその場合アプリ専用領域へ自動的に
  フォールバックする設計のため、致命的な問題にはなりません。
- バックグラウンドでの長時間動作(定期観測サイクル)を安定させたい
  場合は、端末側の設定で本アプリを「バッテリー最適化の対象外」に
  指定することを推奨します。
