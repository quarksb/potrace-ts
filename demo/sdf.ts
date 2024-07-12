import { Path } from "../src/base";
import { getBitmap, loadImage, bitmapToPathArr, imageToSvg } from '../src';

document.getElementById('upload')!.addEventListener('change', generateSDF);

async function generateSDF(event) {
    const file = event.target.files[0];
    if (!file) return;

    const img = await loadImage(file);
    render(img);
}



export function render(img: HTMLImageElement) {
    const { width, height } = img;
    const bitmap = getBitmap(img);
    const turdSize = width * height * 1E-4;
    const pathArr = bitmapToPathArr(bitmap, { turnPolicy: 'right', turdSize });
    const sdfImageData = createSDF(pathArr, img.width, img.height);
    showData(sdfImageData);
}

const showData = (imageData: ImageData) => {
    const canvas = document.getElementById('outputCanvas') as HTMLCanvasElement;
    const ctx = canvas.getContext('2d')!;
    ctx.putImageData(imageData, 0, 0);
}

// render(image);

function createSDF(paths: Path[], width: number, height: number) {
    const sdf = new Float32Array(width * height).fill(Infinity);
    const boundary = new Uint8Array(width * height).fill(0);

    paths.forEach((path) => {
        path.points.forEach((point) => {
            const { x, y } = point;
            boundary[y * width + x] = 1;
            sdf[y * width + x] = 0;
        })
    });

    // 展示 boundary
    const imageDataBoundary = new ImageData(width, height);
    for (let y = 0; y < height; y++) {
        for (let x = 0; x < width; x++) {
            const idx = width * y + x;
            const value = boundary[idx] * 255;

            const pngIdx = idx * 4;
            imageDataBoundary.data[pngIdx] = 0;
            imageDataBoundary.data[pngIdx + 1] = value;
            imageDataBoundary.data[pngIdx + 2] = value;
            imageDataBoundary.data[pngIdx + 3] = 255;
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


