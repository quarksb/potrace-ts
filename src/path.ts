import { Path, Point } from "./base";
import { Bitmap } from "./bitmap";
import { Config } from "./config";

export function bitmapToPathArr(bitmap: Bitmap, config: Config) {
    const bmCopy = bitmap!.copy();
    const pathArr: Path[] = [];

    let currentPoint: false | Point = new Point(0, 0);
    while ((currentPoint = getNextPoint(bmCopy, currentPoint))) {
        const path = getPath(bmCopy, currentPoint, config);
        xorPath(bmCopy, path);
        // 面积超过杂点判定阈值，认定为有效路径
        if (path.area > config.turdSize) {
            pathArr.push(path);
        }
    }
    return pathArr;
}

/**
 * ## 寻找下个有效点
 * @param bitmap
 * @param point
 * @returns
 */
function getNextPoint(bitmap: Bitmap, point: Point) {
    let i = bitmap.w * point.y + point.x;
    while (i < bitmap.size && bitmap.data[i] !== 1) {
        i++;
    }
    return i < bitmap.size && bitmap.index(i);
}



/**
 * 通过上下左右获取路径的基础数据，即边缘点 
 * 
 */
export function getPath(bitmap: Bitmap, point: Point, { turnPolicy: p }: Config) {
    const path = new Path();

    let [dirX, dirY] = [0, 1];
    let { x, y } = point;

    path.sign = bitmap!.checkPixel(point.x, point.y) ? "+" : "-";

    while (true) {
        path.points.push(new Point(x, y));

        path.maxX = Math.max(x, path.maxX);
        path.minX = Math.min(x, path.minX);
        path.maxY = Math.max(y, path.maxY);
        path.minY = Math.min(y, path.minY);

        x += dirX;
        y += dirY;
        path.area -= x * dirY;

        // 回到起始点，构成一个完整的路径，跳出循环
        if (x === point.x && y === point.y) break;

        const isLeftOccupied = bitmap.checkPixel(x + (dirX + dirY - 1) / 2, y + (dirY - dirX - 1) / 2);
        const isRightOccupied = bitmap.checkPixel(x + (dirX - dirY - 1) / 2, y + (dirY + dirX - 1) / 2);

        if (isRightOccupied && !isLeftOccupied) {
            const isTurnRight =
                p === "right" ||
                (p === "black" && path.sign === "+") ||
                (p === "white" && path.sign === "-") ||
                (p === "majority" && isMajorityPoint(bitmap, x, y)) ||
                (p === "minority" && !isMajorityPoint(bitmap, x, y));
            if (isTurnRight) {
                [dirX, dirY] = [-dirY, dirX];
            } else {
                [dirX, dirY] = [dirY, -dirX];
            }
        } else if (isRightOccupied) {
            [dirX, dirY] = [-dirY, dirX];
        } else if (!isLeftOccupied) {
            [dirX, dirY] = [dirY, -dirX];
        }
    }
    return path;
}

/**
 * Determines if a given point is a majority point in the bitmap.
 * A majority point is a point where the majority of its surrounding pixels are set.
 *
 * @param bitmap - The bitmap to check.
 * @param x - The x-coordinate of the point.
 * @param y - The y-coordinate of the point.
 * @returns 1 if the point is a majority point, 0 otherwise.
 */
function isMajorityPoint(bitmap: Bitmap, x: number, y: number) {
    for (let i = 2; i < 5; i++) {
        let ct = 0;
        for (let j = -i + 1; j <= i - 1; j++) {
            ct += bitmap.checkPixel(x + j, y + i - 1) ? 1 : -1;
            ct += bitmap.checkPixel(x + i - 1, y + j - 1) ? 1 : -1;
            ct += bitmap.checkPixel(x + j - 1, y - i) ? 1 : -1;
            ct += bitmap.checkPixel(x - i, y + j) ? 1 : -1;
        }
        if (ct > 0) {
            return 1;
        } else if (ct < 0) {
            return 0;
        }
    }
    return 0;
}

/**
 * ## 翻转 path 内的 bitmap 数据 
 * 高耗操作，或许只需在后续添加个标记符号就可以省略这一步操作
 */
export function xorPath(bitmap: Bitmap, path: Path) {
    let y1 = path.points[0].y;
    let maxX = -Infinity;
    let minY = Infinity;
    for (const point of path.points) {
        const { x, y } = point;

        if (y !== y1) {
            minY = Math.min(y, y1);
            maxX = path.maxX;
            for (let i = x; i < maxX; i++) {
                bitmap.flip(i, minY);
            }
            y1 = y;
        }
    }
}



/**
 * ## Reverses the order of points in a given path.
 * 这样可以调整路径方向，方便判定路径的内外
 * @param path - The path to reverse.
 */
export function reversePath(path: Path) {
    const { curve } = path;
    const { basePoints: p, n } = curve;

    for (let i = 0; 2 * i < n - 1; i++) {
        const j = n - 1 - i;
        [p[i], p[j]] = [p[j], p[i]];
    }
}
