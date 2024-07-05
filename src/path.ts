import { Curve, Opti, Path, Point, Quad, Sum } from "./base";
import { Bitmap } from "./bitmap";
import { Config, DefaultConfig } from "./config";
import { bezier, cross, cyclic, ddenom, dpara, getDist, interval, iprod, iprod1, mod, quadform, xprod } from "./math";


export function bmToPathArr(bitmap: Bitmap, info: Config) {
    const bmCopy = bitmap!.copy()
    const pathArr: Path[] = []

    let currentPoint: false | Point = new Point(0, 0)
    while (currentPoint = findNext(bmCopy, currentPoint)) {
        const path = findPath(bmCopy, currentPoint);
        xOrPath(bmCopy, path);
        // 面积超过杂点判定阈值，认定为有效路径
        if (path.area > info.turdSize) {
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
function findNext(bitmap: Bitmap, point: Point) {
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
            ct += bitmap.at(x + j, y + i - 1) ? 1 : -1;
            ct += bitmap.at(x + i - 1, y + j - 1) ? 1 : -1;
            ct += bitmap.at(x + j - 1, y - i) ? 1 : -1;
            ct += bitmap.at(x - i, y + j) ? 1 : -1;
        }
        if (ct > 0) {
            return 1;
        } else if (ct < 0) {
            return 0;
        }
    }
    return 0;
}

function findPath(bitmap: Bitmap, point: Point) {
    const path = new Path();

    let dirX = 0, dirY = 1;
    let tmp;
    let { x, y } = point;

    path.sign = bitmap!.at(point.x, point.y) ? "+" : "-";

    while (true) {
        path.points.push(new Point(x, y));

        path.maxX = Math.max(x, path.maxX);
        path.minX = Math.min(x, path.minX);
        path.maxY = Math.max(y, path.maxY);
        path.minY = Math.min(y, path.minY);
        path.len++;

        x += dirX;
        y += dirY;
        path.area -= x * dirY;

        if (x === point.x && y === point.y)
            break;

        const l = bitmap.at(x + (dirX + dirY - 1) / 2, y + (dirY - dirX - 1) / 2);
        const r = bitmap.at(x + (dirX - dirY - 1) / 2, y + (dirY + dirX - 1) / 2);

        if (r && !l) {
            if (DefaultConfig.turnPolicy === "right" ||
                (DefaultConfig.turnPolicy === "black" && path.sign === '+') ||
                (DefaultConfig.turnPolicy === "white" && path.sign === '-') ||
                (DefaultConfig.turnPolicy === "majority" && majority(bitmap, x, y)) ||
                (DefaultConfig.turnPolicy === "minority" && !majority(bitmap, x, y))) {
                tmp = dirX;
                dirX = -dirY;
                dirY = tmp;
            } else {
                tmp = dirX;
                dirX = dirY;
                dirY = -tmp;
            }
        } else if (r) {
            tmp = dirX;
            dirX = -dirY;
            dirY = tmp;
        } else if (!l) {
            tmp = dirX;
            dirX = dirY;
            dirY = -tmp;
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
    const n = path.len, points = path.points, sums = path.sums;
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
    const n = path.len, pt = path.points,
        pivk = new Array(n),
        nc = new Array(n),
        ct = new Array(4);
    path.lon = new Array(n);
    let dir: number;

    const constraint = [new Point(), new Point()],
        cur = new Point(),
        off = new Point(),
        dk = new Point()
    let foundk;

    let i, j, k1, a, b, c, d, k = 0;
    for (let i = n - 1; i >= 0; i--) {
        if (pt[i].x != pt[k].x && pt[i].y != pt[k].y) {
            k = i + 1;
        }
        nc[i] = k;
    }

    for (let i = n - 1; i >= 0; i--) {
        ct[0] = ct[1] = ct[2] = ct[3] = 0;
        dir = (3 + 3 * (pt[mod(i + 1, n)].x - pt[i].x) +
            (pt[mod(i + 1, n)].y - pt[i].y)) / 2;
        ct[dir]++;

        constraint[0].x = 0;
        constraint[0].y = 0;
        constraint[1].x = 0;
        constraint[1].y = 0;

        k = nc[i];
        k1 = i;
        while (1) {
            foundk = 0;
            dir = (3 + 3 * Math.sign(pt[k].x - pt[k1].x) +
                Math.sign(pt[k].y - pt[k1].y)) / 2;
            ct[dir]++;

            if (ct[0] && ct[1] && ct[2] && ct[3]) {
                pivk[i] = k1;
                foundk = 1;
                break;
            }

            cur.x = pt[k].x - pt[i].x;
            cur.y = pt[k].y - pt[i].y;

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
            dk.x = Math.sign(pt[k].x - pt[k1].x);
            dk.y = Math.sign(pt[k].y - pt[k1].y);
            cur.x = pt[k1].x - pt[i].x;
            cur.y = pt[k1].y - pt[i].y;

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
    let n = path.len
    let pen = new Array(n + 1),
        prev = new Array(n + 1),
        clip0 = new Array(n),
        clip1 = new Array(n + 1),
        seg0 = new Array(n + 1),
        seg1 = new Array(n + 1),
        thispen, best, c;

    for (i = 0; i < n; i++) {
        c = mod(path.lon[mod(i - 1, n)] - 1, n);
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
    for (i = 0; i < n; i++) {
        while (j <= clip0[i]) {
            clip1[j] = i;
            j++;
        }
    }

    i = 0;
    for (j = 0; i < n; j++) {
        seg0[j] = i;
        i = clip0[i];
    }
    seg0[j] = n;
    m = j;

    i = n;
    for (j = m; j > 0; j--) {
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

    for (i = n, j = m - 1; i > 0; j--) {
        i = prev[i];
        path.po[j] = i;
    }
}

function pointslope(path: Path, i: number, j: number, ctr: Point, dir: Point) {

    let n = path.len, sums = path.sums,
        x, y, x2, xy, y2,
        k, a, b, c, l, r = 0;

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

    x = sums[j + 1].x - sums[i].x + r * sums[n].x;
    y = sums[j + 1].y - sums[i].y + r * sums[n].y;
    x2 = sums[j + 1].x2 - sums[i].x2 + r * sums[n].x2;
    xy = sums[j + 1].xy - sums[i].xy + r * sums[n].xy;
    y2 = sums[j + 1].y2 - sums[i].y2 + r * sums[n].y2;
    k = j + 1 - i + r * n;

    ctr.x = x / k;
    ctr.y = y / k;

    a = (x2 - x * x / k) / k;
    b = (xy - x * y / k) / k;
    c = (y2 - y * y / k) / k;

    const lambda2 = (a + c + Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2;

    a -= lambda2;
    c -= lambda2;

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
}

function adjustVertices(path: Path) {
    let m = path.m, po = path.po, n = path.len, pt = path.points;

    let x0 = path.x0, y0 = path.y0,
        ctr = new Array(m), dir = new Array(m),
        q = new Array(m),
        v = new Array(3), d, i, j, k, l,
        s = new Point();

    path.curve = new Curve(m);

    for (i = 0; i < m; i++) {
        j = po[mod(i + 1, m)];
        j = mod(j - po[i], n) + po[i];
        ctr[i] = new Point();
        dir[i] = new Point();
        pointslope(path, po[i], j, ctr[i], dir[i]);
    }

    for (let i = 0; i < m; i++) {
        q[i] = new Quad();
        d = dir[i].x * dir[i].x + dir[i].y * dir[i].y;
        if (d === 0.0) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    q[i].data[j * 3 + k] = 0;
                }
            }
        } else {
            v[0] = dir[i].y;
            v[1] = -dir[i].x;
            v[2] = - v[1] * ctr[i].y - v[0] * ctr[i].x;
            for (l = 0; l < 3; l++) {
                for (k = 0; k < 3; k++) {
                    q[i].data[l * 3 + k] = v[l] * v[k] / d;
                }
            }
        }
    }

    let Q, w, dx, dy, det, min, cand, xmin, ymin, z;
    for (let i = 0; i < m; i++) {
        Q = new Quad();
        w = new Point();

        s.x = pt[po[i]].x - x0;
        s.y = pt[po[i]].y - y0;

        j = mod(i - 1, m);

        for (l = 0; l < 3; l++) {
            for (k = 0; k < 3; k++) {
                Q.data[l * 3 + k] = q[j].at(l, k) + q[i].at(l, k);
            }
        }

        while (1) {

            det = Q.at(0, 0) * Q.at(1, 1) - Q.at(0, 1) * Q.at(1, 0);
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
            d = v[0] * v[0] + v[1] * v[1];
            v[2] = - v[1] * s.y - v[0] * s.x;
            for (l = 0; l < 3; l++) {
                for (k = 0; k < 3; k++) {
                    Q.data[l * 3 + k] += v[l] * v[k] / d;
                }
            }
        }
        dx = Math.abs(w.x - s.x);
        dy = Math.abs(w.y - s.y);
        if (dx <= 0.5 && dy <= 0.5) {
            path.curve.vertex[i] = new Point(w.x + x0, w.y + y0);
            continue;
        }

        min = quadform(Q, s);
        xmin = s.x;
        ymin = s.y;

        if (Q.at(0, 0) !== 0.0) {
            for (z = 0; z < 2; z++) {
                w.y = s.y - 0.5 + z;
                w.x = - (Q.at(0, 1) * w.y + Q.at(0, 2)) / Q.at(0, 0);
                dx = Math.abs(w.x - s.x);
                cand = quadform(Q, w);
                if (dx <= 0.5 && cand < min) {
                    min = cand;
                    xmin = w.x;
                    ymin = w.y;
                }
            }
        }

        if (Q.at(1, 1) !== 0.0) {
            for (z = 0; z < 2; z++) {
                w.x = s.x - 0.5 + z;
                w.y = - (Q.at(1, 0) * w.x + Q.at(1, 2)) / Q.at(1, 1);
                dy = Math.abs(w.y - s.y);
                cand = quadform(Q, w);
                if (dy <= 0.5 && cand < min) {
                    min = cand;
                    xmin = w.x;
                    ymin = w.y;
                }
            }
        }

        for (l = 0; l < 2; l++) {
            for (k = 0; k < 2; k++) {
                w.x = s.x - 0.5 + l;
                w.y = s.y - 0.5 + k;
                cand = quadform(Q, w);
                if (cand < min) {
                    min = cand;
                    xmin = w.x;
                    ymin = w.y;
                }
            }
        }

        path.curve.vertex[i] = new Point(xmin + x0, ymin + y0);
    }
}

function reverse(path: Path) {
    let curve = path.curve, m = curve.n, v = curve.vertex, i, j, tmp;

    for (i = 0, j = m - 1; i < j; i++, j--) {
        tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
    }
}

function smooth(path: Path) {
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

        if (alpha >= DefaultConfig.alphaMax) {
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
        k, k1, k2, conv, i1,
        area, alpha, d, d1, d2,
        p0, p1, p2, p3, pt,
        A, R, A1, A2, A3, A4,
        s, t;

    if (i == j) {
        return 1;
    }

    k = i;
    i1 = mod(i + 1, m);
    k1 = mod(k + 1, m);
    conv = convc[k1];
    if (conv === 0) {
        return 1;
    }
    d = getDist(vertex[i], vertex[i1]);
    for (k = k1; k != j; k = k1) {
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

    for (k = mod(i + 1, m); k != j; k = k1) {
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

    for (k = i; k != j; k = k1) {
        k1 = mod(k + 1, m);
        t = tangent(p0, p1, p2, p3, curve.points[k * 3 + 2], curve.points[k1 * 3 + 2]);
        if (t < -0.5) {
            return 1;
        }
        pt = bezier(t, p0, p1, p2, p3);
        d = getDist(curve.points[k * 3 + 2], curve.points[k1 * 3 + 2]);
        if (d === 0.0) {
            return 1;
        }
        d1 = dpara(curve.points[k * 3 + 2], curve.points[k1 * 3 + 2], pt) / d;
        d2 = dpara(curve.points[k * 3 + 2], curve.points[k1 * 3 + 2], vertex[k1]) / d;
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

function optiCurve(path: Path) {

    let curve = path.curve,
        m = curve.n, vert = curve.vertex;
    const pt = new Array(m + 1)
    const pen = new Array(m + 1)
    const len = new Array(m + 1)
    const opt = new Array(m + 1);
    let om, r;
    let o = new Opti()
    let p0, i1, area, alpha, ocurve,
        s, t;

    const convc = new Array(m), areac = new Array(m + 1);

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
            r = opti_penalty(path, i, mod(j, m), o, DefaultConfig.optTolerance, convc,
                areac);
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
    ocurve = new Curve(om);
    s = new Array(om);
    t = new Array(om);

    let j = m;
    for (let i = om - 1; i >= 0; i--) {
        if (pt[j] == j - 1) {
            ocurve.tag[i] = curve.tag[mod(j, m)];
            ocurve.points[i * 3 + 0] = curve.points[mod(j, m) * 3 + 0];
            ocurve.points[i * 3 + 1] = curve.points[mod(j, m) * 3 + 1];
            ocurve.points[i * 3 + 2] = curve.points[mod(j, m) * 3 + 2];
            ocurve.vertex[i] = curve.vertex[mod(j, m)];
            ocurve.alpha[i] = curve.alpha[mod(j, m)];
            ocurve.alpha0[i] = curve.alpha0[mod(j, m)];
            ocurve.beta[i] = curve.beta[mod(j, m)];
            s[i] = t[i] = 1.0;
        } else {
            ocurve.tag[i] = "CURVE";
            ocurve.points[i * 3 + 0] = opt[j].c[0];
            ocurve.points[i * 3 + 1] = opt[j].c[1];
            ocurve.points[i * 3 + 2] = curve.points[mod(j, m) * 3 + 2];
            ocurve.vertex[i] = interval(opt[j].s, curve.points[mod(j, m) * 3 + 2],
                vert[mod(j, m)]);
            ocurve.alpha[i] = opt[j].alpha;
            ocurve.alpha0[i] = opt[j].alpha;
            s[i] = opt[j].s;
            t[i] = opt[j].t;
        }
        j = pt[j];
    }

    for (let i = 0; i < om; i++) {
        i1 = mod(i + 1, om);
        ocurve.beta[i] = s[i] / (s[i] + t[i1]);
    }
    path.curve = ocurve;
}

export function processPath(pathArr: Path[]) {
    for (const path of pathArr) {
        calcSums(path);
        calcLon(path);
        bestPolygon(path);
        adjustVertices(path);

        if (path.sign === "-") {
            reverse(path);
        }

        smooth(path);

        if (DefaultConfig.optCurve) {
            optiCurve(path);
        }
    }

}



function tangent(p0: Point, p1: Point, p2: Point, p3: Point, q0: Point, q1: Point) {
    const A = cross(p0, p1, q0, q1);
    const B = cross(p1, p2, q0, q1);
    const C = cross(p2, p3, q0, q1);

    const a = A - 2 * B + C;
    const b = -2 * A + 2 * B;
    const c = A;

    const d = b * b - 4 * a * c;

    if (a === 0 || d < 0) {
        return -1.0;
    }

    const s = Math.sqrt(d);

    const r1 = (-b + s) / (2 * a);
    const r2 = (-b - s) / (2 * a);

    if (r1 >= 0 && r1 <= 1) {
        return r1;
    } else if (r2 >= 0 && r2 <= 1) {
        return r2;
    } else {
        return -1.0;
    }
}

function calcSums(path: Path) {
    path.x0 = path.points[0].x;
    path.y0 = path.points[0].y;

    path.sums = [];
    const s = path.sums;
    s.push(new Sum(0, 0, 0, 0, 0));
    for (let i = 0; i < path.len; i++) {
        const x = path.points[i].x - path.x0;
        const y = path.points[i].y - path.y0;
        const sum = new Sum(s[i].x + x, s[i].y + y, s[i].xy + x * y,
            s[i].x2 + x * x, s[i].y2 + y * y)
        s.push(sum);
    }
}