import numpy as np
import cv2

# https://daizyu.com/posts/2020-05-17-001/

# FPS（Frame Per Second：１秒間に表示するFrame数）
CLIP_FPS = 20.0
filepath = 'test.mp4'

# 画像を読み込む
img = cv2.imread('test-0000.png')

# 読込んだ画像のサイズより、動画サイズを決定する
w = img.shape[1]
h = img.shape[0]
codec = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(filepath, codec, CLIP_FPS, (w, h))

for idx in range(500):
    # 画像ファイルを読み込む
    img = cv2.imread('test-%04d.png'%(idx))

    # カメラから読込んだ映像をファイルに書き込む
    video.write(img)

# Videoを作成時には、開放処理が必要
video.release()
