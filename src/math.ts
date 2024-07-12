import { Point, Matrix3x3 } from "./base";

export function mod(a: number, n: number) {
    return (a + n) % n;
}

/**
 * Checks if number `b` is cyclically between numbers `a` and `c`.
 * 
 * @param a - The first number.
 * @param b - The number to check if it is cyclically between `a` and `c`.
 * @param c - The third number.
 * @returns `true` if `b` is cyclically between `a` and `c`, `false` otherwise.
 */
export function isCyclic(a: number, b: number, c: number) {
    if (a <= c) {
        return (a <= b && b < c);
    } else {
        return (a <= b || b < c);
    }
}

export function quadForm(Q: Matrix3x3, point: Point) {
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

/**
 * ### 返回 |dx| + |dy|
 * 其中 dx = p2.x - p1.x, dy = p2.y - p1.y
 * @param p1 - The starting point of the line segment.
 * @param p2 - The ending point of the line segment.
 * @returns The calculated denominator.
 */
export function getDenom(p1: Point, p2: Point) {
    const dx = p2.x - p1.x;
    const dy = p2.y - p1.y;
    return Math.sign(dx) * dx + Math.sign(dy) * dy;
}

/**
 * Calculates the cross product of two points.
 * @param x1 - The x-coordinate of the first point.
 * @param y1 - The y-coordinate of the first point.
 * @param x2 - The x-coordinate of the second point.
 * @param y2 - The y-coordinate of the second point.
 */
export function cross(x1: number, y1: number, x2: number, y2: number) {
    return x1 * y2 - y1 * x2;
}

/**
 * Calculates the cross product of two points.
 * @param p1 The first point.
 * @param p2 The second point.
 * @returns The cross product of the two points.
 */
export function cross0(p1: Point, p2: Point) {
    return p1.x * p2.y - p1.y * p2.x;
}

/**
 * ## *p0p1 x p2p3*
 * Calculates the cross product of two line segments defined by four points.
 * @param p0 - The starting point of the first line segment.
 * @param p1 - The ending point of the first line segment.
 * @param p2 - The starting point of the second line segment.
 * @param p3 - The ending point of the second line segment.
 * @returns The cross product of the two line segments.
 */
export function cross1(p0: Point, p1: Point, p2: Point, p3: Point) {
    const x1 = p1.x - p0.x;
    const y1 = p1.y - p0.y;
    const x2 = p3.x - p2.x;
    const y2 = p3.y - p2.y;

    return x1 * y2 - x2 * y1;
}

/** ## *p0p1 x p0p2* */
export function cross2(p0: Point, p1: Point, p2: Point) {
    return cross1(p0, p1, p0, p2);
}



/** 
 * ## $p0p1 * p0p2$ 
 */
export function dot2(p0: Point, p1: Point, p2: Point) {
    return dot(p0, p1, p0, p2);
}

/**## ```p0p1 * p2p3``` */
export function dot(p0: Point, p1: Point, p2: Point, p3: Point) {
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

export function getPointOnBezier(t: number, p0: Point, p1: Point, p2: Point, p3: Point) {
    const s = 1 - t;
    const x = s * s * s * p0.x + 3 * (s * s * t) * p1.x + 3 * (t * t * s) * p2.x + t * t * t * p3.x;
    const y = s * s * s * p0.y + 3 * (s * s * t) * p1.y + 3 * (t * t * s) * p2.y + t * t * t * p3.y;
    return new Point(x, y);
}


export function tangent(p0: Point, p1: Point, p2: Point, p3: Point, q0: Point, q1: Point) {
    const A = cross1(p0, p1, q0, q1);
    const B = cross1(p1, p2, q0, q1);
    const C = cross1(p2, p3, q0, q1);

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