import { Point, Quad } from "./base";

export function mod(a: number, n: number) {
    return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n;
}

/**
 * Calculates the cross product of two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @returns The cross product of the two points.
 */
export function xprod(p1: Point, p2: Point) {
    return p1.x * p2.y - p1.y * p2.x;
}

/**
 * Checks if a number `b` is cyclically between numbers `a` and `c`.
 * 
 * @param a - The first number.
 * @param b - The number to check if it is cyclically between `a` and `c`.
 * @param c - The third number.
 * @returns `true` if `b` is cyclically between `a` and `c`, `false` otherwise.
 */
export function cyclic(a: number, b: number, c: number) {
    if (a <= c) {
        return (a <= b && b < c);
    } else {
        return (a <= b || b < c);
    }
}

export function quadform(Q: Quad, point: Point) {
    const v = [point.x, point.y, 1];
    let sum = 0.0;

    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            sum += v[i] * Q.at(i, j) * v[j];
        }
    }
    return sum;
}

export function interval(lambda: number, a: Point, b: Point) {
    const x = a.x + lambda * (b.x - a.x);
    const y = a.y + lambda * (b.y - a.y);
    return new Point(x, y);
}

export function dorth_infty(p1: Point, p2: Point) {
    const x = -Math.sign(p2.y - p1.y);
    const y = Math.sign(p2.x - p1.x);
    return new Point(x, y);
}

export function ddenom(p0: Point, p2: Point) {
    const r = dorth_infty(p0, p2);

    return r.y * (p2.x - p0.x) - r.x * (p2.y - p0.y);
}

/** ```p0p1 X p1p2``` */
export function dpara(p0: Point, p1: Point, p2: Point) {
    const x1 = p1.x - p0.x;
    const y1 = p1.y - p0.y;
    const x2 = p2.x - p0.x;
    const y2 = p2.y - p0.y;

    return x1 * y2 - x2 * y1;
}

/**
 * Calculates the cross product of two line segments defined by four points.
 *
 * @param p0 - The starting point of the first line segment.
 * @param p1 - The ending point of the first line segment.
 * @param p2 - The starting point of the second line segment.
 * @param p3 - The ending point of the second line segment.
 * @returns The cross product of the two line segments.
 */
export function cross(p0: Point, p1: Point, p2: Point, p3: Point) {
    const x1 = p1.x - p0.x;
    const y1 = p1.y - p0.y;
    const x2 = p3.x - p2.x;
    const y2 = p3.y - p2.y;

    return x1 * y2 - x2 * y1;
}

/**```p0p1 * p0p2``` */
export function iprod(p0: Point, p1: Point, p2: Point) {
    const x1 = p1.x - p0.x;
    const y1 = p1.y - p0.y;
    const x2 = p2.x - p0.x;
    const y2 = p2.y - p0.y;

    return x1 * x2 + y1 * y2;
}

/**```p0p1 * p2p3``` */
export function iprod1(p0: Point, p1: Point, p2: Point, p3: Point) {
    const x1 = p1.x - p0.x;
    const y1 = p1.y - p0.y;
    const x2 = p3.x - p2.x;
    const y2 = p3.y - p2.y;

    return x1 * x2 + y1 * y2;
}

export function getDist(p: Point, q: Point) {
    const dx = p.x - q.x;
    const dy = p.y - q.y;
    return Math.sqrt(dx * dx + dy * dy);
}

export function bezier(t: number, p0: Point, p1: Point, p2: Point, p3: Point) {
    const s = 1 - t;
    const x = s * s * s * p0.x + 3 * (s * s * t) * p1.x + 3 * (t * t * s) * p2.x + t * t * t * p3.x;
    const y = s * s * s * p0.y + 3 * (s * s * t) * p1.y + 3 * (t * t * s) * p2.y + t * t * t * p3.y;
    return new Point(x, y);
}