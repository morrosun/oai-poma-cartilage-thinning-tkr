# -*- coding: utf-8 -*-
"""OAIZIB-CM 镜像下载器（HuggingFace 镜像 + Zenodo md5 校验）。

用途：当 zenodo.org 对本机出口 IP 返回 403（风控限流）时，
改从 HuggingFace 镜像 https://hf-mirror.com/datasets/YongchengYAO/OAIZIB-CM 下载同字节文件，
并用 Zenodo 页面公布的 md5 逐文件校验（已实测：镜像 md5 == Zenodo md5）。

用法:
    python download_oaizibcm_hfmirror.py --dst "D:\\OAI_datasets\\OAIZIB-CM"

产出: 6 个 zip（info / template_atlas / labelsTs / labelsTr / imagesTs / imagesTr，共 12.65 GB）
      + _download.log（含每文件 md5 与耗时）

来源:
    Zenodo  10.5281/zenodo.14934086   (Yao Yongcheng, 2025-02-27, CC BY-NC 4.0)
    HF      datasets/YongchengYAO/OAIZIB-CM (main 分支与 Zenodo 同文件)
"""
import argparse
import hashlib
import os
import subprocess
import time

CURL = r"C:\Windows\System32\curl.exe"
BASE_DEFAULT = "https://hf-mirror.com/datasets/YongchengYAO/OAIZIB-CM/resolve/main"

# (文件名, Zenodo 公布 md5, 字节数) —— 小文件优先
FILES = [
    ("info.zip",           "6e6f6827e7766b8115aaace579bf8e79", 51680),
    ("template_atlas.zip", "9565b5d85a819bcbf03e3fbd3d3e9d6b", 3625715),
    ("labelsTs.zip",       "f68973a29e0d57fdcab0f847556bf957", 14122969),
    ("labelsTr.zip",       "30950573ede8ddad2f82f4a330186b11", 57538129),
    ("imagesTs.zip",       "0c57b204ca07ca6be656a188b2971827", 2544663978),
    ("imagesTr.zip",       "5d83ce6dcd10a71a330a1cd792710dfa", 10032464145),
]


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dst", required=True)
    ap.add_argument("--base", default=BASE_DEFAULT)
    args = ap.parse_args()
    os.makedirs(args.dst, exist_ok=True)
    logf = open(os.path.join(args.dst, "_download.log"), "a", encoding="utf-8", buffering=1)

    def log(msg):
        line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
        print(line)
        logf.write(line + "\n")

    log("=== START OAIZIB-CM (mirror %s) -> %s ===" % (args.base, args.dst))
    allok = True
    for name, exp_md5, exp_bytes in FILES:
        path = os.path.join(args.dst, name)
        log("---- %s (expect %d bytes) ----" % (name, exp_bytes))
        if os.path.exists(path) and os.path.getsize(path) == exp_bytes and md5(path) == exp_md5:
            log("SKIP  %s (already complete, md5 ok)" % name)
            continue
        attempt = 0
        while True:
            attempt += 1
            cmd = [CURL, "-L", "-C", "-", "--retry", "12", "--retry-delay", "5",
                   "--retry-all-errors", "--speed-limit", "8192", "--speed-time", "120",
                   "--connect-timeout", "30", "-sS",
                   "-o", path, "%s/%s" % (args.base, name)]
            t0 = time.time()
            r = subprocess.run(cmd, capture_output=True, text=True)
            dt = time.time() - t0
            size = os.path.getsize(path) if os.path.exists(path) else 0
            if r.returncode == 0 and size == exp_bytes and md5(path) == exp_md5:
                log("OK    %s  %d bytes  %.1f MB/s  md5=%s" % (name, size, size / 1e6 / max(dt, .001), exp_md5))
                break
            log("RETRY %s  rc=%s size=%d  %.0fs  %s" % (name, r.returncode, size, dt,
                                                        (r.stderr or "").strip()[:160]))
            if attempt >= 40:
                log("GIVEUP %s" % name)
                allok = False
                break
            time.sleep(5)
    log("=== DONE all_ok=%s ===" % allok)
    logf.close()


if __name__ == "__main__":
    main()
