import { Curve, Opti, Path, Point, Quad, Sum } from "./base";
import { Bitmap } from "./bitmap";
import { Config } from "./config";
import { bezier, cross, cyclic, ddenom, dpara, getDist, interval, iprod, iprod1, mod, quadform, tangent, xprod } from "./math";


export function bmToPathArr(bitmap: Bitmap, config: Config) {
    const bmCopy = bitmap!.copy()
    const pathArr: Path[] = []

    let currentPoint: false | Point = new Point(0, 0)
    while (currentPoint = getNext(bmCopy, currentPoint)) {
        const path = getPath(bmCopy, currentPoint, config);
        xOrPath(bmCopy, path);
        // 面积超过杂点判定阈值，认定为有效路径
        if (path.area > config.turdSize) {
            pathArr.push(path);
        }
    }
    return pathArr;
}

/**
 * ### 寻找下个有效点
 * @param bitmap 
 * @param point 
 * @returns 
 */
function getNext(bitmap: Bitmap, point: Point) {
    let i = bitmap.w * point.y + point.x;
    while (i < bitmap.size && bitmap.data[i] !== 1) {
        i++;
    }
    return i < bitmap.size && bitmap.index(i);
}

/**
 * Determines the majority value of a bitmap at a given position.
 * The majority value is determined by counting the number of set bits (1s) and unset bits (0s) in the surrounding area.
 * If the count is positive, it returns 1. If the count is negative, it returns 0.
 * If the count is zero, it continues to check the surrounding areas with larger radius until a non-zero count is found.
 * @param bitmap - The bitmap to check.
 * @param x - The x-coordinate of the position.
 * @param y - The y-coordinate of the position.
 * @returns The majority value (1 or 0).
 */
function majority(bitmap: Bitmap, x: number, y: number) {
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

function getPath(bitmap: Bitmap, point: Point, { turnPolicy: p }: Config) {
    const path = new Path();

    let dirX = 0, dirY = 1;
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

        if (x === point.x && y === point.y)
            break;

        const l = bitmap.checkPixel(x + (dirX + dirY - 1) / 2, y + (dirY - dirX - 1) / 2);
        const r = bitmap.checkPixel(x + (dirX - dirY - 1) / 2, y + (dirY + dirX - 1) / 2);

        if (r && !l) {
            const isTurnRight = p === "right" || (p === "black" && path.sign === "+") || (p === "white" && path.sign === "-") || (p === "majority" && majority(bitmap, x, y)) || (p === "minority" && !majority(bitmap, x, y));
            if (isTurnRight) {
                [dirX, dirY] = [-dirY, dirX]
            } else {
                [dirX, dirY] = [dirY, -dirX]
            }
        } else if (r) {
            [dirX, dirY] = [-dirY, dirX]
        } else if (!l) {
            [dirX, dirY] = [dirY, -dirX]
        }
    }
    return path;
}

function xOrPath(bitmap: Bitmap, path: Path) {
    let y1 = path.points[0].y
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

/**惩罚函数 */
export function penalty3(path: Path, i: number, j: number) {
    const { points, sums } = path;
    const { length: n } = points

    let x, y, xy, x2, y2,
        k, a, b, c, s,
        px, py, ex, ey,
        r = 0;
    if (j >= n) {
        j -= n;
        r = 1;
    }

    if (r === 0) {
        x = sums[j + 1].x - sums[i].x;
        y = sums[j + 1].y - sums[i].y;
        x2 = sums[j + 1].x2 - sums[i].x2;
        xy = sums[j + 1].xy - sums[i].xy;
        y2 = sums[j + 1].y2 - sums[i].y2;
        k = j + 1 - i;
    } else {
        x = sums[j + 1].x - sums[i].x + sums[n].x;
        y = sums[j + 1].y - sums[i].y + sums[n].y;
        x2 = sums[j + 1].x2 - sums[i].x2 + sums[n].x2;
        xy = sums[j + 1].xy - sums[i].xy + sums[n].xy;
        y2 = sums[j + 1].y2 - sums[i].y2 + sums[n].y2;
        k = j + 1 - i + n;
    }

    px = (points[i].x + points[j].x) / 2.0 - points[0].x;
    py = (points[i].y + points[j].y) / 2.0 - points[0].y;
    ey = (points[j].x - points[i].x);
    ex = -(points[j].y - points[i].y);

    a = ((x2 - 2 * x * px) / k + px * px);
    b = ((xy - x * py - y * px) / k + px * py);
    c = ((y2 - 2 * y * py) / k + py * py);

    s = ex * ex * a + 2 * ex * ey * b + ey * ey * c;

    return Math.sqrt(s);
}

function calcLon(path: Path) {
    const { points } = path;
    const { length: n } = points
    const pivk = new Array(n)
    const nc = new Array(n)
    const ct = new Array(4);
    path.lon = new Array(n);
    let dir: number;

    const constraint = [new Point(), new Point()],
        cur = new Point(),
        off = new Point(),
        dk = new Point()
    let foundk;

    let j, k1, a, b, c, d, k = 0;
    for (let i = n - 1; i >= 0; i--) {
        if (points[i].x != points[k].x && points[i].y != points[k].y) {
            k = i + 1;
        }
        nc[i] = k;
    }

    for (let i = n - 1; i >= 0; i--) {
        ct[0] = ct[1] = ct[2] = ct[3] = 0;
        dir = (3 + 3 * (points[mod(i + 1, n)].x - points[i].x) +
            (points[mod(i + 1, n)].y - points[i].y)) / 2;
        ct[dir]++;

        constraint[0].x = 0;
        constraint[0].y = 0;
        constraint[1].x = 0;
        constraint[1].y = 0;

        k = nc[i];
        k1 = i;
        while (true) {
            foundk = 0;
            dir = (3 + 3 * Math.sign(points[k].x - points[k1].x) +
                Math.sign(points[k].y - points[k1].y)) / 2;
            ct[dir]++;

            if (ct[0] && ct[1] && ct[2] && ct[3]) {
                pivk[i] = k1;
                foundk = 1;
                break;
            }

            cur.x = points[k].x - points[i].x;
            cur.y = points[k].y - points[i].y;

            if (xprod(constraint[0], cur) < 0 || xprod(constraint[1], cur) > 0) {
                break;
            }

            if (Math.abs(cur.x) <= 1 && Math.abs(cur.y) <= 1) {

            } else {
                off.x = cur.x + ((cur.y >= 0 && (cur.y > 0 || cur.x < 0)) ? 1 : -1);
                off.y = cur.y + ((cur.x <= 0 && (cur.x < 0 || cur.y < 0)) ? 1 : -1);
                if (xprod(constraint[0], off) >= 0) {
                    constraint[0].x = off.x;
                    constraint[0].y = off.y;
                }
                off.x = cur.x + ((cur.y <= 0 && (cur.y < 0 || cur.x < 0)) ? 1 : -1);
                off.y = cur.y + ((cur.x >= 0 && (cur.x > 0 || cur.y < 0)) ? 1 : -1);
                if (xprod(constraint[1], off) <= 0) {
                    constraint[1].x = off.x;
                    constraint[1].y = off.y;
                }
            }
            k1 = k;
            k = nc[k1];
            if (!cyclic(k, i, k1)) {
                break;
            }
        }
        if (foundk === 0) {
            dk.x = Math.sign(points[k].x - points[k1].x);
            dk.y = Math.sign(points[k].y - points[k1].y);
            cur.x = points[k1].x - points[i].x;
            cur.y = points[k1].y - points[i].y;

            a = xprod(constraint[0], cur);
            b = xprod(constraint[0], dk);
            c = xprod(constraint[1], cur);
            d = xprod(constraint[1], dk);

            j = 10000000;
            if (b < 0) {
                j = Math.floor(a / -b);
            }
            if (d > 0) {
                j = Math.min(j, Math.floor(-c / d));
            }
            pivk[i] = mod(k1 + j, n);
        }
    }

    j = pivk[n - 1];
    path.lon[n - 1] = j;
    for (let i = n - 2; i >= 0; i--) {
        if (cyclic(i + 1, pivk[i], j)) {
            j = pivk[i];
        }
        path.lon[i] = j;
    }

    for (let i = n - 1; cyclic(mod(i + 1, n), j, path.lon[i]); i--) {
        path.lon[i] = j;
    }
}

function bestPolygon(path: Path) {
    let i, j, m
    const { points } = path;
    const { length: n } = points
    const pen = new Array(n + 1),
        prev = new Array(n + 1),
        clip0 = new Array(n),
        clip1 = new Array(n + 1),
        seg0 = new Array(n + 1),
        seg1 = new Array(n + 1)

    let thispen, best;

    for (let i = 0; i < n; i++) {
        let c = mod(path.lon[mod(i - 1, n)] - 1, n);
        if (c == i) {
            c = mod(i + 1, n);
        }
        if (c < i) {
            clip0[i] = n;
        } else {
            clip0[i] = c;
        }
    }

    j = 1;
    for (let i = 0; i < n; i++) {
        while (j <= clip0[i]) {
            clip1[j] = i;
            j++;
        }
    }

    i = 0;
    for (let j = 0; i < n; j++) {
        seg0[j] = i;
        i = clip0[i];
    }
    seg0[j] = n;
    m = j;

    i = n;
    for (let j = m; j > 0; j--) {
        seg1[j] = i;
        i = clip1[i];
    }
    seg1[0] = 0;

    pen[0] = 0;
    for (let j = 1; j <= m; j++) {
        for (let i = seg1[j]; i <= seg0[j]; i++) {
            best = -1;
            for (let k = seg0[j - 1]; k >= clip1[i]; k--) {
                thispen = penalty3(path, k, i) + pen[k];
                if (best < 0 || thispen < best) {
                    prev[i] = k;
                    best = thispen;
                }
            }
            pen[i] = best;
        }
    }
    path.m = m;
    path.po = new Array(m);

    for (let i = n, j = m - 1; i > 0; j--) {
        i = prev[i];
        path.po[j] = i;
    }
}

function pointSlope(path: Path, i: number, j: number) {
    const { points, sums } = path;
    const { length: n } = points

    let r = 0;
    while (j >= n) {
        j -= n;
        r += 1;
    }
    while (i >= n) {
        i -= n;
        r -= 1;
    }
    while (j < 0) {
        j += n;
        r -= 1;
    }
    while (i < 0) {
        i += n;
        r += 1;
    }


    const x = sums[j + 1].x - sums[i].x + r * sums[n].x;
    const y = sums[j + 1].y - sums[i].y + r * sums[n].y;
    const x2 = sums[j + 1].x2 - sums[i].x2 + r * sums[n].x2;
    const xy = sums[j + 1].xy - sums[i].xy + r * sums[n].xy;
    const y2 = sums[j + 1].y2 - sums[i].y2 + r * sums[n].y2;
    const k = j + 1 - i + r * n;

    const ctr = new Point(x / k, y / k);


    let a = (x2 - x * x / k) / k;
    let b = (xy - x * y / k) / k;
    let c = (y2 - y * y / k) / k;

    const lambda2 = (a + c + Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2;

    a -= lambda2;
    c -= lambda2;

    let l;
    const dir = new Point();
    if (Math.abs(a) >= Math.abs(c)) {
        l = Math.sqrt(a * a + b * b);
        if (l !== 0) {
            dir.x = -b / l;
            dir.y = a / l;
        }
    } else {
        l = Math.sqrt(c * c + b * b);
        if (l !== 0) {
            dir.x = -c / l;
            dir.y = b / l;
        }
    }
    if (l === 0) {
        dir.x = dir.y = 0;
    }

    return [ctr, dir];
}

function adjustVertices(path: Path) {
    let m = path.m, po = path.po;
    const { points } = path;
    const { length: n } = points;

    let { x: x0, y: y0 } = points[0];
    const controlPoints = new Array(m)
    const directions = new Array(m);

    let q = new Array(m),
        v = new Array(3), d, k, l,
        s = new Point();

    path.curve = new Curve(m);

    for (let i = 0; i < m; i++) {
        let j = po[mod(i + 1, m)];
        j = mod(j - po[i], n) + po[i];

        const [controlPoint, direction] = pointSlope(path, po[i], j);
        controlPoints[i] = controlPoint;
        directions[i] = direction;
    }

    for (let i = 0; i < m; i++) {
        q[i] = new Quad();
        d = directions[i].x * directions[i].x + directions[i].y * directions[i].y;
        if (d === 0.0) {
            for (let j = 0; j < 3; j++) {
                for (let k = 0; k < 3; k++) {
                    q[i].data[j * 3 + k] = 0;
                }
            }
        } else {
            v[0] = directions[i].y;
            v[1] = -directions[i].x;
            v[2] = - v[1] * controlPoints[i].y - v[0] * controlPoints[i].x;
            for (let l = 0; l < 3; l++) {
                for (let k = 0; k < 3; k++) {
                    q[i].data[l * 3 + k] = v[l] * v[k] / d;
                }
            }
        }
    }


    for (let i = 0; i < m; i++) {
        const Q = new Quad();
        const w = new Point();

        s.x = points[po[i]].x - x0;
        s.y = points[po[i]].y - y0;

        const j = mod(i - 1, m);

        for (let l = 0; l < 3; l++) {
            for (let k = 0; k < 3; k++) {
                Q.data[l * 3 + k] = q[j].at(l, k) + q[i].at(l, k);
            }
        }

        while (true) {
            const det = Q.at(0, 0) * Q.at(1, 1) - Q.at(0, 1) * Q.at(1, 0);
            if (det !== 0.0) {
                w.x = (-Q.at(0, 2) * Q.at(1, 1) + Q.at(1, 2) * Q.at(0, 1)) / det;
                w.y = (Q.at(0, 2) * Q.at(1, 0) - Q.at(1, 2) * Q.at(0, 0)) / det;
                break;
            }

            if (Q.at(0, 0) > Q.at(1, 1)) {
                v[0] = -Q.at(0, 1);
                v[1] = Q.at(0, 0);
            } else if (Q.at(1, 1)) {
                v[0] = -Q.at(1, 1);
                v[1] = Q.at(1, 0);
            } else {
                v[0] = 1;
                v[1] = 0;
            }
            const d = v[0] * v[0] + v[1] * v[1];
            v[2] = - v[1] * s.y - v[0] * s.x;
            for (let l = 0; l < 3; l++) {
                for (let k = 0; k < 3; k++) {
                    Q.data[l * 3 + k] += v[l] * v[k] / d;
                }
            }
        }
        let dx = Math.abs(w.x - s.x);
        let dy = Math.abs(w.y - s.y);
        if (dx <= 0.5 && dy <= 0.5) {
            path.curve.vertex[i] = new Point(w.x + x0, w.y + y0);
            continue;
        }

        let min = quadform(Q, s);
        let xMin = s.x;
        let yMin = s.y;

        if (Q.at(0, 0) !== 0.0) {
            for (let z = 0; z < 2; z++) {
                w.y = s.y - 0.5 + z;
                w.x = - (Q.at(0, 1) * w.y + Q.at(0, 2)) / Q.at(0, 0);
                dx = Math.abs(w.x - s.x);
                const cand = quadform(Q, w);
                if (dx <= 0.5 && cand < min) {
                    min = cand;
                    xMin = w.x;
                    yMin = w.y;
                }
            }
        }

        if (Q.at(1, 1) !== 0.0) {
            for (let z = 0; z < 2; z++) {
                w.x = s.x - 0.5 + z;
                w.y = - (Q.at(1, 0) * w.x + Q.at(1, 2)) / Q.at(1, 1);
                dy = Math.abs(w.y - s.y);
                const cand = quadform(Q, w);
                if (dy <= 0.5 && cand < min) {
                    min = cand;
                    xMin = w.x;
                    yMin = w.y;
                }
            }
        }

        for (let l = 0; l < 2; l++) {
            for (let k = 0; k < 2; k++) {
                w.x = s.x - 0.5 + l;
                w.y = s.y - 0.5 + k;
                const cand = quadform(Q, w);
                if (cand < min) {
                    min = cand;
                    xMin = w.x;
                    yMin = w.y;
                }
            }
        }

        path.curve.vertex[i] = new Point(xMin + x0, yMin + y0);
    }
}

function reverse(path: Path) {
    let curve = path.curve, m = curve.n;
    const { vertex: v } = curve;

    for (let i = 0, j = m - 1; i < j; i++, j--) {
        const tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
    }
}

function smooth(path: Path, config: Config) {
    let m = path.curve.n;
    const curve = path.curve;

    let j, k, dd, denom, alpha,
        p2, p3, p4;

    for (let i = 0; i < m; i++) {
        j = mod(i + 1, m);
        k = mod(i + 2, m);
        p4 = interval(1 / 2.0, curve.vertex[k], curve.vertex[j]);

        denom = ddenom(curve.vertex[i], curve.vertex[k]);
        if (denom !== 0.0) {
            dd = dpara(curve.vertex[i], curve.vertex[j], curve.vertex[k]) / denom;
            dd = Math.abs(dd);
            alpha = dd > 1 ? (1 - 1.0 / dd) : 0;
            alpha = alpha / 0.75;
        } else {
            alpha = 4 / 3.0;
        }
        curve.alpha0[j] = alpha;

        if (alpha >= config.alphaMax) {
            curve.tag[j] = "CORNER";
            curve.points[3 * j + 1] = curve.vertex[j];
            curve.points[3 * j + 2] = p4;
        } else {
            if (alpha < 0.55) {
                alpha = 0.55;
            } else if (alpha > 1) {
                alpha = 1;
            }
            p2 = interval(0.5 + 0.5 * alpha, curve.vertex[i], curve.vertex[j]);
            p3 = interval(0.5 + 0.5 * alpha, curve.vertex[k], curve.vertex[j]);
            curve.tag[j] = "CURVE";
            curve.points[3 * j + 0] = p2;
            curve.points[3 * j + 1] = p3;
            curve.points[3 * j + 2] = p4;
        }
        curve.alpha[j] = alpha;
        curve.beta[j] = 0.5;
    }
}

function opti_penalty(path: Path, i: number, j: number, res: Opti, opttolerance: number, convc: any[], areac: any[]) {
    let m = path.curve.n, curve = path.curve, vertex = curve.vertex,
        k2, conv,
        area, alpha, d, d1,
        p0, p1, p2, p3, pt,
        A, R, A1, A2, A3, A4,
        s, t;

    if (i == j) {
        return 1;
    }

    const k = i;
    const i1 = mod(i + 1, m);
    let k1 = mod(k + 1, m);
    conv = convc[k1];
    if (conv === 0) {
        return 1;
    }
    d = getDist(vertex[i], vertex[i1]);
    for (let k = k1; k != j; k = k1) {
        k1 = mod(k + 1, m);
        k2 = mod(k + 2, m);
        if (convc[k1] != conv) {
            return 1;
        }
        if (Math.sign(cross(vertex[i], vertex[i1], vertex[k1], vertex[k2])) !=
            conv) {
            return 1;
        }
        if (iprod1(vertex[i], vertex[i1], vertex[k1], vertex[k2]) <
            d * getDist(vertex[k1], vertex[k2]) * -0.999847695156) {
            return 1;
        }
    }

    p0 = curve.points[mod(i, m) * 3 + 2].copy();
    p1 = vertex[mod(i + 1, m)].copy();
    p2 = vertex[mod(j, m)].copy();
    p3 = curve.points[mod(j, m) * 3 + 2].copy();

    area = areac[j] - areac[i];
    area -= dpara(vertex[0], curve.points[i * 3 + 2], curve.points[j * 3 + 2]) / 2;
    if (i >= j) {
        area += areac[m];
    }

    A1 = dpara(p0, p1, p2);
    A2 = dpara(p0, p1, p3);
    A3 = dpara(p0, p2, p3);

    A4 = A1 + A3 - A2;

    if (A2 == A1) {
        return 1;
    }

    t = A3 / (A3 - A4);
    s = A2 / (A2 - A1);
    A = A2 * t / 2.0;

    if (A === 0.0) {
        return 1;
    }

    R = area / A;
    alpha = 2 - Math.sqrt(4 - R / 0.3);

    res.c[0] = interval(t * alpha, p0, p1);
    res.c[1] = interval(s * alpha, p3, p2);
    res.alpha = alpha;
    res.t = t;
    res.s = s;

    p1 = res.c[0].copy();
    p2 = res.c[1].copy();

    res.pen = 0;

    for (let k = mod(i + 1, m); k != j; k = k1) {
        k1 = mod(k + 1, m);
        t = tangent(p0, p1, p2, p3, vertex[k], vertex[k1]);
        if (t < -0.5) {
            return 1;
        }
        pt = bezier(t, p0, p1, p2, p3);
        d = getDist(vertex[k], vertex[k1]);
        if (d === 0.0) {
            return 1;
        }
        d1 = dpara(vertex[k], vertex[k1], pt) / d;
        if (Math.abs(d1) > opttolerance) {
            return 1;
        }
        if (iprod(vertex[k], vertex[k1], pt) < 0 ||
            iprod(vertex[k1], vertex[k], pt) < 0) {
            return 1;
        }
        res.pen += d1 * d1;
    }

    for (let k = i; k != j; k = k1) {
        const k1 = mod(k + 1, m);
        const t = tangent(p0, p1, p2, p3, curve.points[k * 3 + 2], curve.points[k1 * 3 + 2]);
        if (t < -0.5) {
            return 1;
        }
        const pt = bezier(t, p0, p1, p2, p3);
        const d = getDist(curve.points[k * 3 + 2], curve.points[k1 * 3 + 2]);
        if (d === 0.0) {
            return 1;
        }
        let d1 = dpara(curve.points[k * 3 + 2], curve.points[k1 * 3 + 2], pt) / d;
        let d2 = dpara(curve.points[k * 3 + 2], curve.points[k1 * 3 + 2], vertex[k1]) / d;
        d2 *= 0.75 * curve.alpha[k1];
        if (d2 < 0) {
            d1 = -d1;
            d2 = -d2;
        }
        if (d1 < d2 - opttolerance) {
            return 1;
        }
        if (d1 < d2) {
            res.pen += (d1 - d2) * (d1 - d2);
        }
    }

    return 0;
}

function optiCurve(path: Path, config: Config) {

    let curve = path.curve,
        m = curve.n, vert = curve.vertex;
    const pt = new Array(m + 1)
    const pen = new Array(m + 1)
    const len = new Array(m + 1)
    const opt = new Array(m + 1);
    let om, r;
    let o = new Opti()
    let p0, i1, area, alpha;

    const convc = new Array(m);
    const areac = new Array(m + 1);

    for (let i = 0; i < m; i++) {
        if (curve.tag[i] == "CURVE") {
            convc[i] = Math.sign(dpara(vert[mod(i - 1, m)], vert[i], vert[mod(i + 1, m)]));
        } else {
            convc[i] = 0;
        }
    }

    area = 0.0;
    areac[0] = 0.0;
    p0 = curve.vertex[0];
    for (let i = 0; i < m; i++) {
        i1 = mod(i + 1, m);
        if (curve.tag[i1] == "CURVE") {
            alpha = curve.alpha[i1];
            area += 0.3 * alpha * (4 - alpha) *
                dpara(curve.points[i * 3 + 2], vert[i1], curve.points[i1 * 3 + 2]) / 2;
            area += dpara(p0, curve.points[i * 3 + 2], curve.points[i1 * 3 + 2]) / 2;
        }
        areac[i + 1] = area;
    }

    pt[0] = -1;
    pen[0] = 0;
    len[0] = 0;


    for (let j = 1; j <= m; j++) {
        pt[j] = j - 1;
        pen[j] = pen[j - 1];
        len[j] = len[j - 1] + 1;

        for (let i = j - 2; i >= 0; i--) {
            r = opti_penalty(path, i, mod(j, m), o, config.optTolerance, convc, areac);
            if (r) {
                break;
            }
            if (len[j] > len[i] + 1 ||
                (len[j] == len[i] + 1 && pen[j] > pen[i] + o.pen)) {
                pt[j] = i;
                pen[j] = pen[i] + o.pen;
                len[j] = len[i] + 1;
                opt[j] = o;
                o = new Opti();
            }
        }
    }
    om = len[m];
    const newCurve = new Curve(om);
    const s = new Array(om);
    const t = new Array(om);

    let j = m;
    for (let i = om - 1; i >= 0; i--) {
        if (pt[j] == j - 1) {
            newCurve.tag[i] = curve.tag[mod(j, m)];
            newCurve.points[i * 3 + 0] = curve.points[mod(j, m) * 3 + 0];
            newCurve.points[i * 3 + 1] = curve.points[mod(j, m) * 3 + 1];
            newCurve.points[i * 3 + 2] = curve.points[mod(j, m) * 3 + 2];
            newCurve.vertex[i] = curve.vertex[mod(j, m)];
            newCurve.alpha[i] = curve.alpha[mod(j, m)];
            newCurve.alpha0[i] = curve.alpha0[mod(j, m)];
            newCurve.beta[i] = curve.beta[mod(j, m)];
            s[i] = 1.0;
            t[i] = 1.0;
        } else {
            newCurve.tag[i] = "CURVE";
            newCurve.points[i * 3 + 0] = opt[j].c[0];
            newCurve.points[i * 3 + 1] = opt[j].c[1];
            newCurve.points[i * 3 + 2] = curve.points[mod(j, m) * 3 + 2];
            newCurve.vertex[i] = interval(opt[j].s, curve.points[mod(j, m) * 3 + 2],
                vert[mod(j, m)]);
            newCurve.alpha[i] = opt[j].alpha;
            newCurve.alpha0[i] = opt[j].alpha;
            s[i] = opt[j].s;
            t[i] = opt[j].t;
        }
        j = pt[j];
    }

    for (let i = 0; i < om; i++) {
        const i1 = mod(i + 1, om);
        newCurve.beta[i] = s[i] / (s[i] + t[i1]);
    }
    path.curve = newCurve;
}

export function processPath(pathArr: Path[], config: Config) {
    for (const path of pathArr) {
        calcSums(path);
        calcLon(path);
        bestPolygon(path);
        adjustVertices(path);

        if (path.sign === "-") {
            reverse(path);
        }

        smooth(path, config);

        if (config.optCurve) {
            optiCurve(path, config);
        }
    }

}


function calcSums(path: Path) {
    const { points } = path;
    const { x: x0, y: y0 } = points[0];

    path.sums = [];
    const { sums } = path;
    sums.push(new Sum(0, 0, 0, 0, 0));
    for (let i = 0; i < points.length; i++) {
        const x = points[i].x - x0;
        const y = points[i].y - y0;
        const sum = new Sum(sums[i].x + x, sums[i].y + y, sums[i].xy + x * y,
            sums[i].x2 + x * x, sums[i].y2 + y * y)
        sums.push(sum);
    }
}