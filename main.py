#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py  -  ESDUCT Android APK エントリポイント

ESDUCT_discription.py 本体は Kivy に依存しない「HTTPサーバー常駐型」の
コンソール/ダッシュボードアプリケーションであり、GUIは全てブラウザ側の
HTMLダッシュボード(edfs_v13_run_main内蔵のHTTPサーバー、既定ポート8765〜)
で完結する。

このファイルの役割は下記の3点のみ:
  1) Android上でプロセスを常駐させるための最小限のKivyアプリ(Activity)を
     提供する(python-for-androidはActivityコンテキストを必要とするため、
     素のPythonスクリプトを単独APK化することはできない)。
  2) 起動直後にランタイム権限(位置情報等)をリクエストする。
  3) ESDUCT_discription.edfs_v13_run_main() をバックグラウンドスレッドで
     起動し、画面には最小限のステータス表示のみ行う。

ESDUCT_discription.py 側のロジック・学習データ・HTTPサーバー実装には
一切手を加えていない。
"""

import os
import sys
import threading
import traceback

# buildozerでパッケージングされた際、main.pyと同じディレクトリに
# ESDUCT_discription.py が配置されるため、念のため明示的にパスへ追加する。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from kivy.app import App
from kivy.clock import Clock
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.button import Button
from kivy.uix.scrollview import ScrollView

# --------------------------------------------------------------------------
# Android ランタイム権限リクエスト
# --------------------------------------------------------------------------
# 'android.permissions' は python-for-android が Android ビルド時にのみ
# 提供するモジュールであり、PC上での動作確認時(python main.py)には存在
# しないため、ImportError を捕捉して無視する(デスクトップでは権限要求
# 自体が不要なため問題ない)。
def _request_android_permissions():
    try:
        from android.permissions import request_permissions, Permission  # noqa
        request_permissions([
            Permission.INTERNET,
            Permission.ACCESS_FINE_LOCATION,
            Permission.ACCESS_COARSE_LOCATION,
            Permission.WRITE_EXTERNAL_STORAGE,
            Permission.READ_EXTERNAL_STORAGE,
        ])
    except Exception:
        # Android以外の環境、または権限モジュール未提供時は無視する。
        pass


_engine_module = None
_engine_thread = None
_engine_started = False
_engine_error = None


def _run_engine():
    """ESDUCT_discription.edfs_v13_run_main() をバックグラウンドで実行する。
    例外はここで捕捉し、画面表示用にグローバル変数へ格納する
    (バックグラウンドスレッドの例外はそのままではUIに伝わらないため)。
    """
    global _engine_module, _engine_error
    try:
        import ESDUCT_discription as edfs
        _engine_module = edfs
        edfs.edfs_v13_run_main()
    except SystemExit:
        pass
    except Exception as e:
        _engine_error = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
        print("[main.py] ESDUCT engine terminated with exception:")
        print(_engine_error)


class ESDUCTApp(App):
    title = "ESDUCT"

    def build(self):
        _request_android_permissions()

        root = BoxLayout(orientation="vertical", padding=16, spacing=10)

        header = Label(
            text="ESDUCT\n(Es & Duct Propagation Forecast)",
            size_hint=(1, None),
            height=80,
            halign="center",
            valign="middle",
        )
        header.bind(size=lambda *_: setattr(header, "text_size", header.size))
        root.add_widget(header)

        scroll = ScrollView(size_hint=(1, 1))
        self.status_label = Label(
            text=self._initial_message(),
            size_hint_y=None,
            halign="left",
            valign="top",
        )
        self.status_label.bind(
            width=lambda *_: setattr(
                self.status_label, "text_size", (self.status_label.width, None)
            )
        )
        self.status_label.bind(
            texture_size=lambda *_: setattr(
                self.status_label, "height", self.status_label.texture_size[1]
            )
        )
        scroll.add_widget(self.status_label)
        root.add_widget(scroll)

        btn_row = BoxLayout(orientation="horizontal", size_hint=(1, None), height=56, spacing=10)
        stop_btn = Button(text="常駐停止")
        stop_btn.bind(on_release=self._on_stop_pressed)
        btn_row.add_widget(stop_btn)
        root.add_widget(btn_row)

        Clock.schedule_once(self._start_engine, 0.5)
        Clock.schedule_interval(self._refresh_status, 3.0)
        return root

    @staticmethod
    def _initial_message():
        return (
            "起動準備中です。\n\n"
            "しばらくすると自動的にブラウザが開き、\n"
            "ダッシュボードが表示されます。\n\n"
            "自動で開かない場合は、ブラウザで\n"
            "http://127.0.0.1:8765/\n"
            "を開いてください。\n\n"
            "本アプリはバックグラウンドで動作し続けます。\n"
            "終了する場合は下の「常駐停止」ボタンを押してください。"
        )

    def _start_engine(self, *_):
        global _engine_thread, _engine_started
        if _engine_started:
            return
        _engine_started = True
        _engine_thread = threading.Thread(target=_run_engine, daemon=True, name="ESDUCT-engine")
        _engine_thread.start()

    def _refresh_status(self, *_):
        if _engine_error:
            self.status_label.text = (
                "ESDUCTの起動中にエラーが発生しました:\n\n" + _engine_error
            )
            return
        if _engine_module is not None:
            self.status_label.text = (
                "ESDUCTは正常に動作しています。\n\n"
                "ブラウザで以下を開いてダッシュボードを確認してください:\n"
                "http://127.0.0.1:8765/\n\n"
                "(自動的にブラウザが開かない場合のみ、手動で開いてください)\n\n"
                "終了する場合は下の「常駐停止」ボタンを押してください。"
            )

    def _on_stop_pressed(self, *_):
        try:
            if _engine_module is not None:
                _engine_module.request_edfs_stop()
                self.status_label.text = (
                    "常駐停止を要求しました。\n"
                    "実行中のサイクルが完了し次第、安全に停止します。"
                )
        except Exception as e:
            self.status_label.text = f"停止要求でエラーが発生しました: {e}"

    def on_pause(self):
        # バックグラウンド動作(NICT/FT8取得・DUCT予測サイクル)を継続させるため
        # ポーズを許可する(Androidのバッテリー最適化により、画面OFF/長時間の
        # バックグラウンド化ではOSがプロセスを一時停止/終了させる場合がある。
        # 常時稼働させたい場合は、端末側の「バッテリー最適化の対象外」設定を
        # 本アプリに対して行うことを推奨する)。
        return True

    def on_stop(self):
        try:
            if _engine_module is not None:
                _engine_module.request_edfs_stop()
        except Exception:
            pass


if __name__ == "__main__":
    ESDUCTApp().run()
