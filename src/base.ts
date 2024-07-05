export class Point {
    x: number;
    y: number;
    constructor (x = 0, y = 0) {
        this.x = x;
        this.y = y;
    }
    copy() {
        return new Point(this.x, this.y);
    };
}

export class Curve {
    /**点的数量 */
    n: number;
    /**type 直线还是曲线*/
    tag: string[];
    /**曲线控制点数据 */
    points: Point[];
    alphaCurve: number;
    vertex: Point[];
    alpha: number[];
    alpha0: number[];
    beta: number[];
    constructor (n: number) {
        this.n = n;
        this.tag = new Array(n);
        this.points = new Array(n * 3);
        this.alphaCurve = 0;
        this.vertex = new Array(n);
        this.alpha = new Array(n);
        this.alpha0 = new Array(n);
        this.beta = new Array(n);
    }
}

export class Path {
    area: number;
    curve: Curve;
    points: Point[];
    minX: number;
    minY: number;
    maxX: number;
    maxY: number;
    sign: string = "";

    sums: Sum[] = [];
    lon: number[] = [];
    m: number = 0;
    po: number[] = [];
    constructor () {
        this.area = 0;
        this.curve = new Curve(0);
        this.points = [];
        this.minX = 100000;
        this.minY = 100000;
        this.maxX = -1;
        this.maxY = -1;
    }
}

/**Matrix 2D 数据？ */
export class Quad {
    data: number[];
    constructor () {
        this.data = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
    at(x: number, y: number) {
        return this.data[x * 3 + y];
    };
}

export class Sum {
    x: number;
    y: number;
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

export class Opti {
    pen: number;
    c: Point[];
    t: number;
    s: number;
    alpha: number;
    constructor () {
        this.pen = 0;
        this.c = [new Point(), new Point()];
        this.t = 0;
        this.s = 0;
        this.alpha = 0;
    }
}