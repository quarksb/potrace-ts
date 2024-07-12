export class Point {
    x: number;
    y: number;
    constructor (x = 0, y = 0) {
        this.x = x;
        this.y = y;
    }
    set(x: number, y: number) {
        this.x = x;
        this.y = y;
        return this;
    }
    sub(p: Point) {
        this.x -= p.x;
        this.y -= p.y;
        return this;
    }
    static fromSub(p1: Point, p2: Point) {
        return new Point(p1.x - p2.x, p1.y - p2.y);
    }
    copy() {
        return new Point(this.x, this.y);
    };
    /**点的字符串表示 */
    toStr(fractionDigits = 3) {
        return `(${this.x.toFixed(fractionDigits)}, ${this.y.toFixed(fractionDigits)})`;
    }
}

export class Curve {
    /**点的数量 */
    n: number;
    /**type 直线还是曲线
     * - CURVE 曲线
     * - CORNER 直线
    */
    tag: string[];
    /**曲线控制点数据 */
    controlPoints: Point[];
    alphaCurve: number;
    /**基础点数据, 即多边形的每条边必然经过的点 */
    basePoints: Point[];
    alpha: number[];
    constructor (n: number) {
        this.n = n;
        this.tag = new Array(n);
        this.controlPoints = new Array(n * 3);
        this.alphaCurve = 0;
        this.basePoints = new Array(n);
        this.alpha = new Array(n);
    }
}

/**
 * Represents a path in the potrace library.
 */
export class Path {
    /**
     * The area enclosed by the path.
     */
    area: number = 0;
    /**
     * The curve associated with the path.
     */
    curve: Curve = new Curve(0);
    /**
     * The points that make up the path.
     */
    points: Point[] = [];
    /**
     * The minimum x-coordinate of the path.
     */
    minX: number = 20000;
    /**
     * The minimum y-coordinate of the path.
     */
    minY: number = 20000;
    /**
     * The maximum x-coordinate of the path.
     */
    maxX: number = -1;
    /**
     * The maximum y-coordinate of the path.
     */
    maxY: number = -1;
    /**
     * The sign of the path.
     * - "+" : positive
     * - "-" : negative
     */
    sign: string = "+";
    /**
     * The sums associated with the path.
     */
    sums: Sum[] = [];
    /**
     * The lon values of the path.
     */
    lon: number[] = [];
    /**
     * The polygon data of the path.
     * 记录单条线的终点
     */
    polygon: number[] = [];
}

/**Matrix 2D 数据？ */
export class Matrix3x3 {
    data: number[];
    constructor () {
        this.data = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
    at(x: number, y: number) {
        return this.data[x * 3 + y];
    };
    apply(p: Point) {
        return new Point(
            this.data[0] * p.x + this.data[1] * p.y + this.data[2],
            this.data[3] * p.x + this.data[4] * p.y + this.data[5]
        );
    }
}

export class Sum {
    x: number;
    y: number;
    /**x * y */
    xy: number;
    /**x^2 */
    x2: number;
    /**y^2 */
    y2: number;
    constructor (x: number, y: number, xy: number, x2: number, y2: number) {
        this.x = x;
        this.y = y;
        this.xy = xy;
        this.x2 = x2;
        this.y2 = y2;
    }
}

/**优化数据 */
export class OptData {
    /**优化惩罚, 或者说误差 */
    penalty: number;
    curve: Point[];
    alpha: number;
    constructor () {
        this.penalty = 0;
        this.curve = [new Point(), new Point()];
        this.alpha = 0;
    }
}