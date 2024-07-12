import { getBitmap, loadImage, bitmapToPathArr } from '../src';
import { Path } from '../src/base';
import { Bitmap } from '../src/bitmap';
import { getPath } from '../src/path';
import imgUrl from "/quarksb.png";

console.time('loadImage');
const image = await loadImage(imgUrl);
const canvas = document.getElementById('outputCanvas') as HTMLCanvasElement;
const ctx = canvas.getContext('2d')!;
canvas.width = image.width;
canvas.height = image.height;
console.timeEnd('loadImage');

console.time('getBitmap');
const bitmap = getBitmap(image);
// showBitMap(bitmap);
console.log('bitmap:', bitmap);

console.timeEnd('getBitmap');

console.time('bitmapToPathArr');
const pathArr = bitmapToPathArr(bitmap, { turnPolicy: 'right', turdSize: 10 });
console.timeEnd('bitmapToPathArr');

// console.table(pathArr);
// pathArr.forEach((path, i) => {
//     console.log(`path ${i}:`, path.points.length, path.area);
// })
pathArr.forEach((path) => {
    ctx.fillStyle = getRandomColor();
    showPath(path);
})
function showPath(path: Path) {
    const { points, minX, minY } = path;
    // ctx.translate(minX, minY);
    points.forEach((point, i) => {
        const size = 2;
        ctx.fillRect(point.x - size / 2, point.y - size / 2, size, size);
    })
    // ctx.translate(-minX, -minY);
}

function getRandomColor() {
    return '#' + Math.random().toString(16).slice(2, 8);
}

function showBitMap(bitmap: Bitmap) {
    const { w, h, data } = bitmap;
    const imgData = new ImageData(w, h);
    for (let i = 0; i < data.length; i++) {
        const x = i % w;
        const y = Math.floor(i / w);
        const idx = (w * y + x) * 4;
        const color = data[i] === 1 ? 0 : 125;
        imgData.data[idx] = color;
        imgData.data[idx + 1] = color;
        imgData.data[idx + 2] = color;
        imgData.data[idx + 3] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
}




document.getElementById('upload')!.addEventListener('change', generateSDF);

async function generateSDF(event) {
    const file = event.target.files[0];
    if (!file) return;

    const img = await loadImage(file);

}



const render = (img: HTMLImageElement) => {
    const inputCanvas = new OffscreenCanvas(img.width, img.height);
    const inputCtx = inputCanvas.getContext('2d')!;

    inputCtx.drawImage(img, 0, 0);

    console.time("get image data")
    const imageData = inputCtx.getImageData(0, 0, img.width, img.height);
    console.timeEnd("get image data")
    const sdfImageData = createSDF(imageData);
    showData(sdfImageData);
}

const showData = (imageData: ImageData) => {

    ctx.putImageData(imageData, 0, 0);
}

// render(image);


function createSDF(imageData: ImageData) {
    const { width, height, data } = imageData;
    const sdf = new Float32Array(width * height).fill(Infinity);
    const boundary = new Uint8Array(width * height).fill(0);

    for (let y = 1; y < height - 1; y++) {
        for (let x = 1; x < width - 1; x++) {
            const idx = (width * y + x) * 4;
            const alpha = data[idx + 3];

            if (alpha > 128) {
                const alphaLeft = data[(width * y + (x - 1)) * 4 + 3];
                const alphaRight = data[(width * y + (x + 1)) * 4 + 3];
                const alphaUp = data[(width * (y - 1) + x) * 4 + 3];
                const alphaDown = data[(width * (y + 1) + x) * 4 + 3];

                if (alphaLeft <= 128 || alphaRight <= 128 || alphaUp <= 128 || alphaDown <= 128) {
                    boundary[width * y + x] = 1;
                    sdf[width * y + x] = 0.0;
                }
            }
        }
    }

    for (let y = 1; y < height; y++) {
        for (let x = 1; x < width; x++) {
            const idx = width * y + x;
            if (sdf[idx] > 0) {
                sdf[idx] = Math.min(
                    sdf[idx],
                    sdf[idx - width] + 1,
                    sdf[idx - 1] + 1,
                    sdf[idx - width - 1] + Math.SQRT2,
                    sdf[idx - width + 1] + Math.SQRT2
                );
            }
        }
    }

    for (let y = height - 2; y >= 0; y--) {
        for (let x = width - 2; x >= 0; x--) {
            const idx = width * y + x;
            if (sdf[idx] > 0) {
                sdf[idx] = Math.min(
                    sdf[idx],
                    sdf[idx + width] + 1,
                    sdf[idx + 1] + 1,
                    sdf[idx + width + 1] + Math.SQRT2,
                    sdf[idx + width - 1] + Math.SQRT2
                );
            }
        }
    }

    let maxDist = 0;
    for (let i = 0; i < sdf.length; i++) {
        if (sdf[i] === Infinity) sdf[i] = 0;
        maxDist = Math.max(maxDist, sdf[i]);
        // console.log('sdf:', sdf[i]);

    }
    console.log('maxDist:', maxDist);

    const sdfImageData = new ImageData(width, height);

    for (let y = 0; y < height; y++) {
        for (let x = 0; x < width; x++) {
            const idx = width * y + x;
            const dist = sdf[idx];
            const value = Math.floor((dist / maxDist) * 255);

            const pngIdx = idx * 4;
            sdfImageData.data[pngIdx] = value;
            sdfImageData.data[pngIdx + 1] = value;
            sdfImageData.data[pngIdx + 2] = value;
            sdfImageData.data[pngIdx + 3] = 255;
        }
    }

    return sdfImageData;
}

