import { Bitmap } from "./bitmap";
import { Config, DefaultConfig } from "./config";
import { bmToPathArr, processPath } from "./path";





export function process(image: HTMLImageElement, config: Partial<Config> = {}) {
    const bitmap: Bitmap = getBitmap(image)
    const realConfig = { ...DefaultConfig, ...config };
    const pathArr = bmToPathArr(bitmap, realConfig);
    processPath(pathArr);
}

/**
 * 获取灰度图数据
 */
function getBitmap(image: HTMLImageElement) {
    const { width, height } = image;
    const imgCanvas = new OffscreenCanvas(width, height);
    const ctx = imgCanvas.getContext('2d')!;
    const bitmap = new Bitmap(imgCanvas.width, imgCanvas.height);
    const imgData = ctx.getImageData(0, 0, bitmap.w, bitmap.h);
    const { data } = imgData;
    for (let i = 0, j = 0; i < data.length; i += 4, j++) {
        const color = 0.2126 * data[i] + 0.7153 * data[i + 1] + 0.0721 * data[i + 2];
        bitmap.data[j] = Number(color < 128)
    }
    return bitmap;
}

export function loadImage(src: string) {
    return new Promise<HTMLImageElement>((resolve, reject) => {
        const img = new Image();
        img.onload = () => {
            resolve(img);
        }
        img.onerror = reject;
        img.src = src;
    })
}