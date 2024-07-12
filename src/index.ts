import { Bitmap } from "./bitmap";
import { Config, DefaultConfig } from "./config";
import { adjustVertices, bestPolygon, straightPath, calcSums, optimizeCurve, getCurveData } from "./curve";
import { bitmapToPathArr, reversePath } from "./path";
import { pathArrToSVG } from "./svg";

/**
 * Converts an image to an SVG string.
 * @param image 
 * @param config 
 * @returns 
 */
export function imageToSvg(image: HTMLImageElement, config: Partial<Config> = {}) {
    const bitmap: Bitmap = getBitmap(image)
    const realConfig = { ...DefaultConfig, ...config };
    const pathArr = processPath(bitmap, realConfig);
    const w = bitmap.w * realConfig.scale;
    const h = bitmap.h * realConfig.scale;
    return pathArrToSVG(pathArr, { w, h }, config.isStroke);
}

export { bitmapToPathArr }

/**
 * ## 获取灰度图数据, 并二值化
 */
export function getBitmap(image: HTMLImageElement) {
    const { width, height } = image;
    const canvas = new OffscreenCanvas(width, height);
    const ctx = canvas.getContext('2d')!;
    const bitmap = new Bitmap(canvas.width, canvas.height);
    ctx.drawImage(image, 0, 0);
    const imgData = ctx.getImageData(0, 0, bitmap.w, bitmap.h);
    const { data } = imgData;
    for (let i = 0; 4 * i < data.length; i++) {
        const j = i * 4;
        const [r, g, b, a] = data.slice(j, j + 4);
        const color = 0.2126 * r * a + 0.7152 * g * a + 0.0722 * b * a;
        // assume that the background color of the image is white, and the foreground color is black.
        bitmap.data[i] = color < 128 ? 1 : 0;
    }
    return bitmap;
}


/**
 * Processes an array of paths using the provided configuration.
 * The Potrace algorithm transforms a bitmap into a vector outline in several steps. 
 * - In the first step, the bitmap is decomposed into a number of paths, which form the boundaries between black and white areas. 
 * - In the second step, each path is approximated by an optimal polygon.
 * - In the third step, each polygon is transformed into a smooth outline.
 * - In an optional fourth step, the resulting curve is optimized by joining consecutive Bezier curve segments together where this is possible. 
 * @param pathArr - An array of paths to be processed.
 * @param config - The configuration object.
 */
export function processPath(bitmap: Bitmap, config: Config) {
    // first step: decompose the bitmap into paths
    const pathArr = bitmapToPathArr(bitmap, config);
    for (const path of pathArr) {
        calcSums(path);
        straightPath(path);
        bestPolygon(path);


        adjustVertices(path);
        if (path.sign === "-") {
            reversePath(path);
        }

        getCurveData(path, config);

        if (config.optCurve) {
            // fourth step: optimize the curve
            optimizeCurve(path, config);
        }
    }
    return pathArr;
}

export function loadImage(src: string | File) {
    return new Promise<HTMLImageElement>((resolve, reject) => {
        const img = new Image();
        img.onload = () => {
            resolve(img);
        }
        img.onerror = reject;
        img.src = src instanceof File ? URL.createObjectURL(src) : src;
    })
}


