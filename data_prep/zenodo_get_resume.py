# -*- coding: utf-8 -*-
"""Generic resumable Zenodo record downloader with md5 verification.

用法:
    python zenodo_get_resume.py 18745638 --out "D:\\OAI_datasets\\CTh-Maps"
    python zenodo_get_resume.py 10.5281/zenodo.18745638 --out DIR

行为:
  1. 调 https://zenodo.org/api/records/<id> 取文件清单(key/size/checksum/links.self)
  2. 逐个 curl -C - 断点续传下载到 <out>\\<key>（保持 Zenodo 的目录层级）
  3. 每个文件按 Zenodo 公布的 md5 校验
  4. 已存在且 md5 正确的文件自动跳过

前置条件: 网络能访问 zenodo.org（当前环境下该域名对本机出口 IP 返回 403，
需在可访问 Zenodo 的网络/代理下运行）。
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request

CURL = r"C:\Windows\System32\curl.exe"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")


def norm_record_id(s: str) -> str:
    s = s.strip()
    if "/" in s:
        s = s.rsplit("/", 1)[-1]
    if s.lower().startswith("zenodo."):
        s = s[7:]
    return s


def fetch_json(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": UA,
                                              "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


def md5(path: str, chunk: int = 1 << 22) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download(url: str, path: str) -> int:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cmd = [CURL, "-L", "-C", "-", "--retry", "12", "--retry-delay", "5",
           "--retry-all-errors", "--speed-limit", "8192", "--speed-time", "120",
           "--connect-timeout", "30", "-A", UA, "-sS", "-o", path, url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("   curl rc=%s %s" % (r.returncode, (r.stderr or "").strip()[:200]))
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("record", help="Zenodo record id 或 DOI，例如 18745638 / 10.5281/zenodo.18745638")
    ap.add_argument("--out", required=True, help="下载目录")
    ap.add_argument("--api", default="https://zenodo.org/api/records",
                    help="Zenodo API 前缀（默认 zenodo.org）")
    args = ap.parse_args()

    rid = norm_record_id(args.record)
    api = "%s/%s" % (args.api.rstrip("/"), rid)
    print("== record %s ==" % rid)
    try:
        rec = fetch_json(api)
    except Exception as e:
        print("!! 无法获取记录元数据: %s: %s" % (type(e).__name__, e))
        print("   若为 403，说明当前出口 IP 被 Zenodo 风控拦截；请在可访问 Zenodo 的网络下重试。")
        sys.exit(2)

    files = rec.get("files", [])
    print("title : %s" % (rec.get("metadata", {}).get("title", "?")))
    print("files : %d  (total %.2f GB)" % (
        len(files), sum(f.get("size", 0) for f in files) / 1e9))
    ok = True
    for f in files:
        key = f["key"]
        size = f.get("size")
        checksum = f.get("checksum", "")            # "md5:xxxx"
        exp = checksum.split(":", 1)[1] if checksum.startswith("md5:") else None
        url = f.get("links", {}).get("self")
        path = os.path.join(args.out, key)
        print("\n-- %s  (%s bytes)" % (key, size))
        if os.path.exists(path) and size and os.path.getsize(path) == size:
            got = md5(path)
            if not exp or got == exp:
                print("   SKIP (complete, md5 ok)")
                continue
        rc = download(url, path)
        if not os.path.exists(path):
            print("   MISSING after download"); ok = False; continue
        got = md5(path)
        if exp and got != exp:
            print("   md5 MISMATCH got=%s want=%s" % (got, exp)); ok = False
        elif size and os.path.getsize(path) != size:
            print("   SIZE MISMATCH %d != %d" % (os.path.getsize(path), size)); ok = False
        else:
            print("   OK  %d bytes  md5=%s" % (os.path.getsize(path), got))

    print("\n== inventory ==")
    for f in files:
        p = os.path.join(args.out, f["key"])
        print("  %-60s %s" % (f["key"],
              ("%d" % os.path.getsize(p)) if os.path.exists(p) else "MISSING"))
    print("\nALL_OK=%s" % ok)


if __name__ == "__main__":
    main()
