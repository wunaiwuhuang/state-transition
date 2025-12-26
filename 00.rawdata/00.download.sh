#!/bin/bash

# 切换到脚本所在目录（日志会写在同一目录）
cd "$(dirname "$0")"

nohup ~/software/meiji/cli-downloader-linux \
  -k 64bcdb7ce7d420901147c124269970ee \
  -o /disk2/cai113/data/stateTrans/00.rawdata \
  > download.log 2>&1 &

echo "Download started in background. Log: $(pwd)/download.log"

