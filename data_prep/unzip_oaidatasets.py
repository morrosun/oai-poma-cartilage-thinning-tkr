# -*- coding: utf-8 -*-
r"""
批量解包 OAI_datasets 下的数据集 zip（可复用）。

用途
----
本机 bash 环境**没有 `unzip` 命令**，统一改用 Python `zipfile`；
同时可顺带剔除 macOS 打包冗余（__MACOSX/、.DS_Store、._*、.~lock*）。

用法
----
    python unzip_oaidatasets.py                      # 解包默认两个数据集
    python unzip_oaidatasets.py --root D:\OAI_datasets
    python unzip_oaidatasets.py --only OAIZIB-CM
    python unzip_oaidatasets.py --keep-zip           # 默认即保留 zip（本脚本从不删 zip）

约定
----
- 目录即目标：zip 内部的顶层目录名会被保留（nnUNet 风格 imagesTr/ labelsTr/ 等）。
- 逐 5000 文件打点写日志，便于长任务观测进度。
- 解包不校验 md5；md5 校验请用 download 脚本（download_oaizibcm_hfmirror.py / zenodo_get_resume.py）。
"""
import os
import sys
import json
import time
import zipfile
import argparse

DEFAULT_ROOT = r'D:\OAI_datasets'


def is_junk(name: str) -> bool:
    if name.startswith('__MACOSX/') or '/__MACOSX/' in name:
        return True
    base = os.path.basename(name.rstrip('/'))
    return base == '.DS_Store' or base.startswith('._') or base.startswith('.~lock')


def extract(zpath, dest, log, step=5000):
    z = zipfile.ZipFile(zpath)
    n = nbytes = 0
    skipped = []
    t0 = time.time()
    for info in z.infolist():
        if is_junk(info.filename):
            skipped.append(info.filename)
            continue
        z.extract(info, dest)
        if not info.is_dir():
            n += 1
            nbytes += info.file_size
        if n and n % step == 0:
            log('   ... %s: %d files, %.2f GB, %.0fs'
                % (os.path.basename(zpath), n, nbytes / 1e9, time.time() - t0))
    z.close()
    log('DONE %s -> %s | files=%d bytes=%.2f GB junk_skipped=%d elapsed=%.0fs'
        % (os.path.basename(zpath), dest, n, nbytes / 1e9, len(skipped), time.time() - t0))
    return {'files': n, 'bytes': nbytes, 'junk_skipped': skipped}


def datasets(root, only=None):
    """返回 {标签: (dest, [zip 文件名...])}。顺序：小文件先，大文件后。"""
    plan = {
        'CTh-Maps': (os.path.join(root, 'CTh-Maps'),
                     ['OAI_CTh-Maps_CTh-Score_Dataset_v1.zip']),
        'OAIZIB-CM': (os.path.join(root, 'OAIZIB-CM'),
                      ['info.zip', 'template_atlas.zip', 'labelsTs.zip',
                       'labelsTr.zip', 'imagesTs.zip', 'imagesTr.zip']),
    }
    if only:
        plan = {k: v for k, v in plan.items() if k in only}
    return plan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default=DEFAULT_ROOT)
    ap.add_argument('--only', nargs='*', default=None,
                    help='只解包指定数据集（CTh-Maps / OAIZIB-CM）')
    args = ap.parse_args()

    logpath = os.path.join(args.root, '_unzip.log')
    with open(logpath, 'w', encoding='utf-8'):
        pass

    def log(msg):
        line = '[%s] %s' % (time.strftime('%H:%M:%S'), msg)
        with open(logpath, 'a', encoding='utf-8') as f:
            f.write(line + '\n')
        print(line, flush=True)

    log('=== unzip start (root=%s) ===' % args.root)
    summary = {}
    for label, (dest, zips) in datasets(args.root, args.only).items():
        log('[%s]' % label)
        os.makedirs(dest, exist_ok=True)
        summary[label] = {}
        for f in zips:
            zp = os.path.join(dest, f)
            if not os.path.exists(zp):
                log('   MISSING %s' % f)
                continue
            summary[label][f] = extract(zp, dest, log)
            summary[label][f].pop('junk_skipped', None)
    log('=== unzip end ===')

    with open(os.path.join(args.root, '_unzip_summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
