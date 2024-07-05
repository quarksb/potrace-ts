import { Curve, Path } from "./base";
import { Bitmap } from "./bitmap";

export function getSVG(bitmap: Bitmap, pathArr: Path[], size: number, opt_type: string) {
    const w = bitmap.w * size
    const h = bitmap.h * size
    let len = pathArr.length, c;

    let svg = '<svg id="svg" version="1.1" width="' + w + '" height="' + h +
        '" xmlns="http://www.w3.org/2000/svg">';
    svg += '<path d="';
    for (let i = 0; i < len; i++) {
        const { curve } = pathArr[i];
        svg += path(curve, size);
    }
    let strokec, fillc, fillrule;
    if (opt_type === "curve") {
        strokec = "black";
        fillc = "none";
        fillrule = '';
    } else {
        strokec = "none";
        fillc = "black";
        fillrule = ' fill-rule="evenodd"';
    }
    svg += `" stroke="${strokec}" fill="${fillc}"${fillrule}/></svg>`;
    return svg;
}

function path(curve: Curve, size: number) {
    let n = curve.n;
    let p = 'M' + (curve.points[(n - 1) * 3 + 2].x * size).toFixed(3) +
        ' ' + (curve.points[(n - 1) * 3 + 2].y * size).toFixed(3) + ' ';
    for (let i = 0; i < n; i++) {
        if (curve.tag[i] === "CURVE") {
            p += bezier(curve, size, i);
        } else if (curve.tag[i] === "CORNER") {
            p += segment(curve, size, i);
        }
    }
    //p += 
    return p;
}

function bezier(curve: Curve, size: number, i: number) {
    const simple = (x: number) => (x * size).toFixed(3);
    let b = `C ${simple(curve.points[i * 3 + 0].x)} ${(curve.points[i * 3 + 0].y * size).toFixed(3)},`;
    b += (curve.points[i * 3 + 1].x * size).toFixed(3) + ' ' +
        (curve.points[i * 3 + 1].y * size).toFixed(3) + ',';
    b += (curve.points[i * 3 + 2].x * size).toFixed(3) + ' ' +
        (curve.points[i * 3 + 2].y * size).toFixed(3) + ' ';
    return b;
}

function segment(curve: Curve, size: number, i: number) {
    let s = 'L ' + (curve.points[i * 3 + 1].x * size).toFixed(3) + ' ' +
        (curve.points[i * 3 + 1].y * size).toFixed(3) + ' ';
    s += (curve.points[i * 3 + 2].x * size).toFixed(3) + ' ' +
        (curve.points[i * 3 + 2].y * size).toFixed(3) + ' ';
    return s;
}