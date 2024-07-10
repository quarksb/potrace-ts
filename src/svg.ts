import { Curve, Path, Point } from "./base";

export function getSVG(pathArr: Path[], size: { w: number, h: number }, isStroke = false) {
    const { w, h } = size
    let len = pathArr.length;

    let svg = `<svg id="svg" viewBox="0 0 ${w} ${h}" version="1.1" width="${w}" height="${h}" xmlns="http://www.w3.org/2000/svg">`;
    svg += '<path d="';
    for (let i = 0; i < len; i++) {
        const { curve } = pathArr[i];
        svg += getPathStr(curve);
    }
    let stroke, fill, fillRule;
    if (isStroke) {
        stroke = "black";
        fill = "none";
        fillRule = '';
    } else {
        stroke = "none";
        fill = "black";
        fillRule = ' fill-rule="evenodd"';
    }
    svg += `" stroke="${stroke}" fill="${fill}"${fillRule}/></svg>`;
    return svg;
}

const simple = (num: number) => num.toFixed(3);

function getPathStr(curve: Curve) {
    const { n, controlPoints: points, tag } = curve;
    let pathStr = `M${(points[(n - 1) * 3 + 2].x).toFixed(3)} ${(points[(n - 1) * 3 + 2].y).toFixed(3)} `;
    for (let i = 0; i < n; i++) {
        if (tag[i] === "CURVE") {
            pathStr += bezier(points, 3 * i);
        } else if (tag[i] === "CORNER") {
            pathStr += segment(points, 3 * i);
        }
    }
    return pathStr;
}

function bezier(points: Point[], i: number) {
    let b = `C ${simple(points[i + 0].x)} ${(points[i + 0].y).toFixed(3)},`;
    b += `${simple(points[i + 1].x)} ${(points[i + 1].y).toFixed(3)},`;
    b += `${simple(points[i + 2].x)} ${(points[i + 2].y).toFixed(3)} `;
    return b;
}

/** */
function segment(points: Point[], i: number) {
    return `L ${simple(points[i + 2].x)} ${simple(points[i + 2].y)}`;
}