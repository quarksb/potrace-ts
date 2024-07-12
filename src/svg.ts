import { Curve, Path, Point } from "./base";

export function pathArrToSVG(pathArr: Path[], size: { w: number, h: number }, isStroke = false) {
    const { w, h } = size
    let len = pathArr.length;

    let svg = `<svg id="svg" viewBox="0 0 ${w} ${h}" version="1.1" width="${w}" height="${h}" xmlns="http://www.w3.org/2000/svg">`;
    svg += '<path d="';
    for (let i = 0; i < len; i++) {
        const { curve } = pathArr[i];
        svg += getPathStr(curve);
    }
    let stroke = "none";
    let fill = "black";
    let fillRule = ' fill-rule="evenodd"';
    if (isStroke) {
        stroke = "black";
        fill = "none";
        fillRule = '';
    }
    svg += `" stroke="${stroke}" fill="${fill}"${fillRule}/></svg>`;
    return svg;
}


function getPathStr(curve: Curve) {
    const { n, controlPoints: points, tag } = curve;
    let pathStr = `M${points[(n - 1) * 3 + 2].toString()} `;
    for (let i = 0; i < n; i++) {
        if (tag[i] === "CURVE") {
            pathStr += getBezierStr(points, 3 * i);
        } else if (tag[i] === "CORNER") {
            pathStr += getLineStr(points, 3 * i);
        }
    }
    return pathStr;
}

function getBezierStr(points: Point[], i: number) {
    return `C ${points[i + 0].toString()}, ${points[i + 1].toString()}, ${points[i + 2].toString()}`;
}

/** */
function getLineStr(points: Point[], i: number) {
    return `L ${points[i + 2].toString()}`;
}