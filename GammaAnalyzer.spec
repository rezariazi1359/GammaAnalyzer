# -*- mode: python ; coding: utf-8 -*-
import sys ; sys.setrecursionlimit(sys.getrecursionlimit() * 5)

from PyInstaller.utils.hooks import collect_submodules

# Collect hidden imports
hiddenimports = []
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('pymedphys')

# Explicitly add missing files
datas = [
    ('C:/Users/rreiazi01/AppData/Roaming/JupyterLab-Portable-3.1.0-3.9/apps/lib/site-packages/pymedphys/_imports/imports.py', 'pymedphys/_imports'),
    ('C:/Users/rreiazi01/AppData/Roaming/JupyterLab-Portable-3.1.0-3.9/apps/lib/site-packages/pymedphys/_trf/decode/config.json', 'pymedphys/_trf/decode')
]

block_cipher = None

a = Analysis(
    ['GammaAnalyzer.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,   # 👉 use the list you built
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='GammaAnalyzer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,   # GUI app, no console
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
