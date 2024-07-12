import { Curve, OptData as OptData, Path, Point, Matrix3x3, Sum } from "./base";
import { Config } from "./config";
import { bezier, cross1, cross2, cross0, isCyclic, dot, dot2, getDenom, getDist, interval, mod, quadForm, tangent } from "./math";

/**### 惩罚函数 */
function getPenalty(path: Path, i: number, j: number) {
    const { points, sums } = path;
    const { length: n } = points;

    let r = 0;
    if (j >= n) {
        j -= n;
        r = 1;
    }

    const x = sums[j + 1].x - sums[i].x + r * sums[n].x;
    const y = sums[j + 1].y - sums[i].y + r * sums[n].y;
    const x2 = sums[j + 1].x2 - sums[i].x2 + r * sums[n].x2;
    const xy = sums[j + 1].xy - sums[i].xy + r * sums[n].xy;
    const y2 = sums[j + 1].y2 - sums[i].y2 + r * sums[n].y2;
    const k = j + 1 - i + r * n;

    const px = (points[i].x + points[j].x) / 2.0 - points[0].x;
    const py = (points[i].y + points[j].y) / 2.0 - points[0].y;
    const ey = points[j].x - points[i].x;
    const ex = -(points[j].y - points[i].y);

    const a = (x2 - 2 * x * px) / k + px * px;
    const b = (xy - x * py - y * px) / k + px * py;
    const c = (y2 - 2 * y * py) / k + py * py;

    const s = ex * ex * a + 2 * ex * ey * b + ey * ey * c;

    return Math.sqrt(s);
}

/**
 * ## straight path
 * 通过合并相邻的点，生成直线路径，拉直路径
 * 
 * It is clear from the definition that if a path is straight, then so are all its subpaths.
 * In order to compute whether a given path is straight or not, we use the stronger fact
 * that straightness is a triplewise property, in the following sense. Suppose that a given 
 * path p = {v0,...,vn} does not use all four directions. Then p is straight if and only if 
 * for all triples (i, j,k) of indices such that 0 6 i < j < k 6 n, there exists a point w on 
 * the straight line through vi and vk such that d(vj,w) 6 1. This observation gives rise to 
 * a naive straightness testing algorithm that is of cubic complexity in the worst case; it 
 * proceeds simply by testing the above property for all triples (i, j,k). 
 * 
 * In the Potrace implementation, we use an optimization that allows us to find all 
 * straight subpaths of a given closed path of length n in time O(n2) in the worst case. 
 * Briefly, the trick is to compute, for every pair (i, j), a constraint on the position of all
 * future vk’s. If i is fixed and j is increasing, it suffices to check the constraint once for 
 * each j. Moreover, a constraint consist of at most two inequalities and can be updated 
 * and checked in constant time.
 * @param path - The path object containing the points.
 */
export function straightPath(path: Path) {
    const { points: p } = path;
    const { length: n } = p;

    const nc = new Array(n);

    path.lon = new Array(n);

    let k = 0;
    for (let i = n - 1; i >= 0; i--) {
        if (p[i].x != p[k].x && p[i].y != p[k].y) {
            k = i + 1;
        }
        nc[i] = k;
    }

    /**
     * ## getDirIndex
     * [0, 1] -> 2 ,[0, -1] -> 1 , [1, 0] -> 3, [-1, 0] -> 0
     * 编码方式有很多种，这里是一种， (3 + dx + 3 * dy) / 2 也行
     */
    const getDirIndex = (dx: number, dy: number) => (3 + 3 * dx + dy) / 2;

    const pivk = new Array(n);
    /**方向约束向量，用来记录和 dir 的极值*/
    const cVec = [new Point(), new Point()];
    for (let i = n - 1; i >= 0; i--) {
        const dirCountArr = [0, 0, 0, 0];
        const vec = Point.fromSub(p[mod(i + 1, n)], p[i]);
        const dirIndex = getDirIndex(vec.x, vec.y);
        dirCountArr[dirIndex]++;

        cVec[0].set(0, 0);
        cVec[1].set(0, 0);

        let k = nc[i];
        let k1 = i;
        let isEnd = false;
        while (true) {
            const vec = Point.fromSub(p[k], p[k1]);
            const dirIndex = getDirIndex(Math.sign(vec.x), Math.sign(vec.y));
            dirCountArr[dirIndex]++;

            // 如果四种方向都有，认为不是直线，退出
            if (dirCountArr[0] && dirCountArr[1] && dirCountArr[2] && dirCountArr[3]) {
                pivk[i] = k1;
                isEnd = true;
                break;
            }

            const curV = Point.fromSub(p[k], p[i]);
            if (cross0(cVec[0], curV) < 0 || cross0(cVec[1], curV) > 0) {
                break;
            }

            if (Math.abs(curV.x) <= 1 && Math.abs(curV.y) <= 1) {
            } else {
                const off = new Point();
                off.x = curV.x + (curV.y >= 0 && (curV.y > 0 || curV.x < 0) ? 1 : -1);
                off.y = curV.y + (curV.x <= 0 && (curV.x < 0 || curV.y < 0) ? 1 : -1);
                if (cross0(cVec[0], off) >= 0) {
                    cVec[0].x = off.x;
                    cVec[0].y = off.y;
                }
                off.x = curV.x + (curV.y <= 0 && (curV.y < 0 || curV.x < 0) ? 1 : -1);
                off.y = curV.y + (curV.x >= 0 && (curV.x > 0 || curV.y < 0) ? 1 : -1);
                if (cross0(cVec[1], off) <= 0) {
                    cVec[1].x = off.x;
                    cVec[1].y = off.y;
                }
            }
            k1 = k;
            k = nc[k1];
            if (!isCyclic(k, i, k1)) {
                break;
            }
        }
        if (isEnd) {
            const dk = new Point();
            dk.x = Math.sign(p[k].x - p[k1].x);
            dk.y = Math.sign(p[k].y - p[k1].y);
            const cur = Point.fromSub(p[k1], p[i]);

            const a = cross0(cVec[0], cur);
            const b = cross0(cVec[0], dk);

            let j = 1E7;
            if (b < 0) {
                j = Math.floor(a / -b);
            }

            const c = cross0(cVec[1], cur);
            const d = cross0(cVec[1], dk);
            if (d > 0) {
                j = Math.min(j, Math.floor(-c / d));
            }
            pivk[i] = mod(k1 + j, n);
        }
    }

    let j = pivk[n - 1];
    path.lon[n - 1] = j;
    for (let i = n - 2; i >= 0; i--) {
        if (isCyclic(i + 1, pivk[i], j)) {
            j = pivk[i];
        }
        path.lon[i] = j;
    }

    for (let i = n - 1; isCyclic(mod(i + 1, n), j, path.lon[i]); i--) {
        path.lon[i] = j;
    }
}

export function bestPolygon(path: Path) {

    const { points } = path;
    const { length: n } = points;
    const penalties = new Array(n + 1);
    /** prevArr[i] 记录了 point[i] 所在直线的起点的 index */
    const prevArr = new Array(n + 1);
    const clip0 = new Array(n);
    const clip1 = new Array(n + 1);
    const segments0 = new Array(n + 1);
    const segments1 = new Array(n + 1);

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

    let j = 1;
    for (let i = 0; i < n; i++) {
        while (j <= clip0[i]) {
            clip1[j] = i;
            j++;
        }
    }

    let i = 0;
    for (let j = 0; i < n; j++) {
        segments0[j] = i;
        i = clip0[i];
    }
    segments0[j] = n;

    let m = j;
    i = n;
    for (let j = m; j > 0; j--) {
        segments1[j] = i;
        i = clip1[i];
    }
    segments1[0] = 0;
    penalties[0] = 0;

    for (let j = 1; j <= m; j++) {
        for (let i = segments1[j]; i <= segments0[j]; i++) {
            let best = -1;
            for (let k = segments0[j - 1]; k >= clip1[i]; k--) {
                const thisPen = getPenalty(path, k, i) + penalties[k];
                if (best < 0 || thisPen < best) {
                    prevArr[i] = k;
                    best = thisPen;
                }
            }
            penalties[i] = best;
        }
    }
    path.polygon = new Array(m);

    for (let i = n, j = m - 1; i > 0; j--) {
        i = prevArr[i];
        path.polygon[j] = i;
    }
}

/**
 * this line passes through the 
 * center of gravity (E(xk),E(yk)), where k = ik,...,ik+1, and its slope is given by the 
 * eigenvector of the larger eigenvalue of the matrix  
 * a b 
 * b c 
 * , where 
 * a = E(xj^2)−E(xj)^2,
 * b = E(xj * yj)−E(xj)*E(yj),
 * c = E(yj^2)−E(yj)^2.
 * @param path 
 * @param i 
 * @param j 
 * @returns 
 */
function getLineInfo(path: Path, i: number, j: number) {
    const { points, sums } = path;
    const { length: n } = points;

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

    const point = new Point(x / k, y / k);

    let a = (x2 - (x * x) / k) / k;
    const b = (xy - (x * y) / k) / k;
    let c = (y2 - (y * y) / k) / k;

    /**
     * 矩阵
     * 
     *   [ a, b ]
     * 
     *   [ b, c ]
     * 
     * 较大的特征值, another eigenvalue is
     * 
     * ```lambda1 = (a + c - Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2;```
     */
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

    return { point, dir };
}

/**
 * In particular, if the intersection of Lk−1,k and Lk,k+1 lies in this unit square, 
 * then we place ak at the intersection; else, we place it at a point near vik 
 * that is “close” to the intersection
 * @param path 
 */
export function adjustVertices(path: Path) {
    let { polygon, points } = path;
    const { length: n } = points;
    const { length: m } = polygon;

    let { x: x0, y: y0 } = points[0];

    /**记录多边形的每条边经过的点 */
    const basePoints: Point[] = new Array(m);
    /**记录多边形的每条边的方向向量 */
    const directions: Point[] = new Array(m);


    path.curve = new Curve(m);

    for (let i = 0; i < m; i++) {
        let j = polygon[mod(i + 1, m)];
        j = mod(j - polygon[i], n) + polygon[i];

        const { point, dir } = getLineInfo(path, polygon[i], j);
        basePoints[i] = point;
        directions[i] = dir;
    }

    const matrixArr: Matrix3x3[] = new Array(m);
    const v = new Array(3);
    for (let i = 0; i < m; i++) {
        matrixArr[i] = new Matrix3x3();
        const { x, y } = directions[i];
        const d = x * x + y * y;
        if (d === 0.0) {
            for (let j = 0; j < 3; j++) {
                for (let k = 0; k < 3; k++) {
                    matrixArr[i].data[j * 3 + k] = 0;
                }
            }
        } else {
            v[0] = y;
            v[1] = -x;
            v[2] = x * basePoints[i].y - y * basePoints[i].x;
            for (let j = 0; j < 3; j++) {
                for (let k = 0; k < 3; k++) {
                    matrixArr[i].data[j * 3 + k] = (v[j] * v[k]) / d;
                }
            }
        }
    }

    for (let i = 0; i < m; i++) {
        const vec = Point.fromSub(points[polygon[mod(i + 1, m)]], points[0]);
        const j = mod(i - 1, m);

        const mat = new Matrix3x3();
        for (let l = 0; l < 3; l++) {
            for (let k = 0; k < 3; k++) {
                mat.data[l * 3 + k] = matrixArr[j].at(l, k) + matrixArr[i].at(l, k);
            }
        }

        const p = new Point();
        while (true) {
            const det = mat.at(0, 0) * mat.at(1, 1) - mat.at(0, 1) * mat.at(1, 0);
            if (det !== 0.0) {
                p.x = (-mat.at(0, 2) * mat.at(1, 1) + mat.at(1, 2) * mat.at(0, 1)) / det;
                p.y = (mat.at(0, 2) * mat.at(1, 0) - mat.at(1, 2) * mat.at(0, 0)) / det;
                break;
            }

            if (mat.at(0, 0) > mat.at(1, 1)) {
                v[0] = -mat.at(0, 1);
                v[1] = mat.at(0, 0);
            } else if (mat.at(1, 1)) {
                v[0] = -mat.at(1, 1);
                v[1] = mat.at(1, 0);
            } else {
                v[0] = 1;
                v[1] = 0;
            }
            const d = v[0] * v[0] + v[1] * v[1];
            v[2] = -v[1] * vec.y - v[0] * vec.x;
            for (let l = 0; l < 3; l++) {
                for (let k = 0; k < 3; k++) {
                    mat.data[l * 3 + k] += (v[l] * v[k]) / d;
                }
            }
        }
        let dx = Math.abs(p.x - vec.x);
        let dy = Math.abs(p.y - vec.y);
        if (dx <= 0.5 && dy <= 0.5) {
            path.curve.basePoints[i] = new Point(p.x + x0, p.y + y0);
            continue;
        }

        let min = quadForm(mat, vec);
        let xMin = vec.x;
        let yMin = vec.y;

        if (mat.at(0, 0) !== 0.0) {
            for (let z = 0; z < 2; z++) {
                p.y = vec.y - 0.5 + z;
                p.x = -(mat.at(0, 1) * p.y + mat.at(0, 2)) / mat.at(0, 0);
                const dx = Math.abs(p.x - vec.x);
                const cand = quadForm(mat, p);
                if (dx <= 0.5 && cand < min) {
                    min = cand;
                    xMin = p.x;
                    yMin = p.y;
                }
            }
        }

        if (mat.at(1, 1) !== 0.0) {
            for (let z = 0; z < 2; z++) {
                p.x = vec.x - 0.5 + z;
                p.y = -(mat.at(1, 0) * p.x + mat.at(1, 2)) / mat.at(1, 1);
                const dy = Math.abs(p.y - vec.y);
                const cand = quadForm(mat, p);
                if (dy <= 0.5 && cand < min) {
                    min = cand;
                    xMin = p.x;
                    yMin = p.y;
                }
            }
        }

        for (let l = 0; l < 2; l++) {
            for (let k = 0; k < 2; k++) {
                p.x = vec.x - 0.5 + l;
                p.y = vec.y - 0.5 + k;
                const cand = quadForm(mat, p);
                if (cand < min) {
                    min = cand;
                    xMin = p.x;
                    yMin = p.y;
                }
            }
        }

        path.curve.basePoints[i] = new Point(xMin + x0, yMin + y0);
    }
}

/**
 * ##  cal the bezier control points
 * Smooths the given path using the specified configuration.
 * 
 * @param path - The path to be smoothed.
 * @param config - The configuration for smoothing.
 */
export function getCurveData(path: Path, config: Config) {
    const { curve } = path;
    let { n, basePoints: points, controlPoints } = curve;

    for (let i = 0; i < n; i++) {
        /**下一个点 */
        const j = mod(i + 1, n);
        /**下下个点 */
        const k = mod(i + 2, n);
        /**在下一点和下下点之间插值出一个点 */
        controlPoints[3 * j + 2] = interval(0.5, points[k], points[j]);;

        const denom = getDenom(points[i], points[k]);
        let alpha = 1;
        if (denom !== 0.0) {
            let dd = Math.abs(cross2(points[i], points[j], points[k]) / denom);
            alpha = dd > 1 ? 1 - 1.0 / dd : 0;
        }
        alpha *= 4 / 3.0;
        curve.alpha0[j] = alpha;

        // alpha 大于设定的阈值，认为是拐角，否则是曲线
        if (alpha >= config.alphaMax) {
            curve.tag[j] = "CORNER";
            controlPoints[3 * j + 1] = points[j];

        } else {
            curve.tag[j] = "CURVE";
            alpha = Math.max(0.55, alpha);
            alpha = Math.min(1, alpha);
            const lambda = 0.5 + 0.5 * alpha;
            const p2 = interval(lambda, points[i], points[j]);
            const p3 = interval(lambda, points[k], points[j]);

            controlPoints[3 * j + 0] = p2;
            controlPoints[3 * j + 1] = p3;
        }
        curve.alpha[j] = alpha;
    }
}

/**
 * Calculates the penalty for optimizing a path segment.
 * 
 * @param path - The path object.
 * @param i - The starting index of the path segment.
 * @param j - The ending index of the path segment.
 * @param res - The result object to store the calculated values.
 * @param optTolerance - The optimization tolerance value.
 * @param convexArr - The array of convexity values. 0 表示直线，不凹凸
 * @param areaArr - The array of area values.
 * @returns {boolean} is Optimize failure.
 */
function optimizeWithPenalty(path: Path, i: number, j: number, res: OptData, optTolerance: number, convexArr: number[], areaArr: number[]): boolean {
    const { curve } = path;
    const { controlPoints, n, basePoints: vertex } = curve;

    if (i == j) {
        return true;
    }

    const k = i;
    const i1 = mod(i + 1, n);
    /**下一个点的 index */
    let k1 = mod(k + 1, n);
    const convexNum = convexArr[k1];
    if (convexNum === 0) {
        return true;
    }
    const dist = getDist(vertex[i], vertex[i1]);
    for (let k = k1; k != j; k = k1) {
        const k1 = mod(k + 1, n);
        const k2 = mod(k + 2, n);
        // 前后两条曲线凹凸性不一致，返回 1
        if (convexArr[k1] != convexNum) {
            return true;
        }
        if (Math.sign(cross1(vertex[i], vertex[i1], vertex[k1], vertex[k2])) != convexNum) {
            return true;
        }
        if (dot(vertex[i], vertex[i1], vertex[k1], vertex[k2]) < dist * getDist(vertex[k1], vertex[k2]) * -0.999847695156) {
            return true;
        }
    }


    let area = areaArr[j] - areaArr[i];
    area -= cross2(vertex[0], controlPoints[i * 3 + 2], controlPoints[j * 3 + 2]) / 2;
    if (i >= j) {
        area += areaArr[n];
    }

    const p0 = controlPoints[mod(i, n) * 3 + 2].copy();
    const p3 = controlPoints[mod(j, n) * 3 + 2].copy();
    let p1 = vertex[mod(i + 1, n)].copy();
    let p2 = vertex[mod(j, n)].copy();
    const A1 = cross2(p0, p1, p2);
    const A2 = cross2(p0, p1, p3);
    const A3 = cross2(p0, p2, p3);

    const A4 = A1 + A3 - A2;

    if (A2 == A1) {
        return true;
    }

    const t = A3 / (A3 - A4);
    const s = A2 / (A2 - A1);
    const A = (A2 * t) / 2.0;

    if (A === 0.0) {
        return true;
    }

    const R = area / A;
    const alpha = 2 - Math.sqrt(4 - R / 0.3);

    res.curve[0] = interval(t * alpha, p0, p1);
    res.curve[1] = interval(s * alpha, p3, p2);
    res.alpha = alpha;

    p1 = res.curve[0].copy();
    p2 = res.curve[1].copy();

    res.penalty = 0;

    for (let k = mod(i + 1, n); k != j; k = k1) {
        const k1 = mod(k + 1, n);
        const t = tangent(p0, p1, p2, p3, vertex[k], vertex[k1]);
        if (t < -0.5) {
            return true;
        }
        const pt = bezier(t, p0, p1, p2, p3);
        const d = getDist(vertex[k], vertex[k1]);
        if (d === 0.0) {
            return true;
        }
        const d1 = cross2(vertex[k], vertex[k1], pt) / d;
        if (Math.abs(d1) > optTolerance) {
            return true;
        }
        if (dot2(vertex[k], vertex[k1], pt) < 0 || dot2(vertex[k1], vertex[k], pt) < 0) {
            return true;
        }
        res.penalty += d1 * d1;
    }

    for (let k = i; k != j; k = k1) {
        const k1 = mod(k + 1, n);
        const t = tangent(p0, p1, p2, p3, controlPoints[k * 3 + 2], controlPoints[k1 * 3 + 2]);
        if (t < -0.5) {
            return true;
        }
        const pt = bezier(t, p0, p1, p2, p3);
        const d = getDist(controlPoints[k * 3 + 2], controlPoints[k1 * 3 + 2]);
        if (d === 0.0) {
            return true;
        }
        let d1 = cross2(controlPoints[k * 3 + 2], controlPoints[k1 * 3 + 2], pt) / d;
        let d2 = cross2(controlPoints[k * 3 + 2], controlPoints[k1 * 3 + 2], vertex[k1]) / d;
        d2 *= 0.75 * curve.alpha[k1];
        if (d2 < 0) {
            d1 = -d1;
            d2 = -d2;
        }
        if (d1 < d2 - optTolerance) {
            return true;
        }
        if (d1 < d2) {
            res.penalty += (d1 - d2) * (d1 - d2);
        }
    }

    return false;
}

/** 尝试合并曲线 */
export function optimizeCurve(path: Path, config: Config) {
    let { curve } = path;
    const { n, controlPoints, basePoints: vert } = curve;
    const pt = new Array(n + 1); // 记录最优路径点
    const penalties = new Array(n + 1); // 记录各路径的惩罚值
    const len = new Array(n + 1); // 记录路径长度
    const opt = new Array(n + 1); // 记录优化参数


    /**
     * - 1:  表示凹， 
     * - 0:  表示直线， 
     * - -1: 表示凸
     */
    const convexNumArr = new Array(n); // 记录曲线的凸性信息
    const areaArr: number[] = new Array(n + 1); // 记录面积信息
    areaArr[0] = 0.0; // 初始面积为0

    // 计算每段曲线的凸性
    for (let i = 0; i < n; i++) {
        if (curve.tag[i] == "CURVE") {
            convexNumArr[i] = Math.sign(cross2(vert[mod(i - 1, n)], vert[i], vert[mod(i + 1, n)]));
        } else {
            convexNumArr[i] = 0;
        }
    }

    const p0 = curve.basePoints[0]; // 初始点
    // 计算每段曲线的面积
    for (let i = 0; i < n; i++) {
        const i1 = mod(i + 1, n);
        let area = 0;
        if (curve.tag[i1] == "CURVE") {
            const alpha = curve.alpha[i1];
            area += (0.3 * alpha * (4 - alpha) * cross2(controlPoints[i * 3 + 2], vert[i1], controlPoints[i1 * 3 + 2])) / 2;
            area += cross2(p0, controlPoints[i * 3 + 2], controlPoints[i1 * 3 + 2]) / 2;
        }
        areaArr[i + 1] = area;
    }

    pt[0] = -1;
    penalties[0] = 0;
    len[0] = 0;

    // 计算最优路径和惩罚值
    for (let j = 1; j <= n; j++) {
        pt[j] = j - 1;
        penalties[j] = penalties[j - 1];
        len[j] = len[j - 1] + 1;

        // 初始化优化对象
        const optData = new OptData();
        // 通过惩罚函数寻找最优路径
        for (let i = j - 2; i >= 0; i--) {
            const failed = optimizeWithPenalty(path, i, mod(j, n), optData, config.optTolerance, convexNumArr, areaArr);
            if (failed) break;

            if (len[j] > len[i] + 1 || (len[j] == len[i] + 1 && penalties[j] > penalties[i] + optData.penalty)) {
                pt[j] = i;
                penalties[j] = penalties[i] + optData.penalty;
                len[j] = len[i] + 1;
                opt[j] = optData;
            }
        }
    }

    // 生成新曲线
    const optimizedMount = len[n];
    const newCurve = new Curve(optimizedMount);

    let j = n;
    for (let i = optimizedMount - 1; i >= 0; i--) {
        if (pt[j] == j - 1) {
            newCurve.tag[i] = curve.tag[mod(j, n)];
            newCurve.controlPoints[i * 3 + 0] = controlPoints[mod(j, n) * 3 + 0];
            newCurve.controlPoints[i * 3 + 1] = controlPoints[mod(j, n) * 3 + 1];
            newCurve.controlPoints[i * 3 + 2] = controlPoints[mod(j, n) * 3 + 2];
            newCurve.basePoints[i] = curve.basePoints[mod(j, n)];
            newCurve.alpha[i] = curve.alpha[mod(j, n)];
            newCurve.alpha0[i] = curve.alpha0[mod(j, n)];
        } else {
            newCurve.tag[i] = "CURVE";
            newCurve.controlPoints[i * 3 + 0] = opt[j].c[0];
            newCurve.controlPoints[i * 3 + 1] = opt[j].c[1];
            newCurve.controlPoints[i * 3 + 2] = controlPoints[mod(j, n) * 3 + 2];
            newCurve.basePoints[i] = interval(opt[j].s, controlPoints[mod(j, n) * 3 + 2], vert[mod(j, n)]);
            newCurve.alpha[i] = opt[j].alpha;
            newCurve.alpha0[i] = opt[j].alpha;
        }
        j = pt[j];
    }

    path.curve = newCurve; // 更新路径曲线
}




/**### 预处理数据，方便后续计算 */
export function calcSums(path: Path) {
    const { points } = path;
    const { x: x0, y: y0 } = points[0];

    const sums: Sum[] = new Array(points.length + 1)
    sums[0] = new Sum(0, 0, 0, 0, 0);
    for (let i = 1; i < points.length; i++) {
        const x = points[i].x - x0;
        const y = points[i].y - y0;
        const sum = new Sum(sums[i].x + x, sums[i].y + y, sums[i].xy + x * y, sums[i].x2 + x * x, sums[i].y2 + y * y);
        sums[i] = (sum);
    }

    path.sums = sums;
}