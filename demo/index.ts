import { getBitmap, loadImage, bitmapToPathArr, imageToSvg } from '../src';
import { Path } from '../src/base';
import { Bitmap } from '../src/bitmap';
import { render } from './sdf';
import imgUrl from "/quarksb.png";

console.time('loadImage');
const image = await loadImage(imgUrl);
const { width, height } = image;

const canvas = document.getElementById('outputCanvas') as HTMLCanvasElement;
const ctx = canvas.getContext('2d')!;
canvas.width = image.width;
canvas.height = image.height;
console.timeEnd('loadImage');

render(image);

// console.time('getBitmap');
// const bitmap = getBitmap(image);
// // showBitMap(bitmap);
// console.log('bitmap:', bitmap);

// // console.timeEnd('getBitmap');

// console.time('bitmapToPathArr');
// const turdSize = width * height / 10000;
// const pathArr = bitmapToPathArr(bitmap, { turnPolicy: 'right', turdSize });
// console.timeEnd('bitmapToPathArr');

// const svg = imageToSvg(image);
// document.getElementById('svg')!.innerHTML = svg;

// console.table(pathArr);
// pathArr.forEach((path, i) => {
//     console.log(`path ${i}:`, path.points.length, path.area);
// })
// pathArr.forEach((path) => {
//     ctx.fillStyle = getRandomColor();
//     showPath(path);
// })
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





