import cv2
import numpy as np

cv2.namedWindow('res')
img = cv2.imread('./brown-eyes.jpg')
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
img_blur = cv2.medianBlur(gray, 5)
# edges = cv2.Canny(img_blur, 50, 150)
p1 = 270
p2 = 15
min_r = 12
max_r = 30

# test good_result (p1 270 p2 15 min_r 12 max_r 30)

def test():
    global p1, p2, min_r, max_r, img_blur, img
    img_copy = np.copy(img)
    if p1 == 0:
        p1 = 1
    if p2 == 0:
        p2 = 1
    if min_r == 0:
        min_r = 1
    if max_r == 0:
        max_r = 1
    circles = cv2.HoughCircles(img_blur, cv2.HOUGH_GRADIENT, 1, img.shape[0] / 64,
                               param1=p1, param2=p2, minRadius=min_r, maxRadius=max_r)
    if circles is not None:
        circles = np.uint16(np.round(circles))
        for i in circles[0, :]:
            cv2.circle(img_copy, (i[0], i[1]), i[2], (0, 255, 0), 2)

    cv2.imshow('res', img_copy)


def update_param1(value):
    global p1
    p1 = value
    test()


def update_param2(value):
    global p2
    p2 = value
    test()


def update_min_r(value):
    global min_r
    min_r = value
    test()


def update_max_r(value):
    global max_r
    max_r = value
    test()


cv2.createTrackbar('param1', 'res', 270, 300, update_param1)
cv2.createTrackbar('param2', 'res', 15, 100, update_param2)
cv2.createTrackbar('min_r', 'res', 12, 30, update_min_r)
cv2.createTrackbar('max_r', 'res', 30, 30, update_max_r)
test()
cv2.waitKey()
